import os
import scanpy as sc
import pandas as pd
import warnings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
import numpy as np
from flask import Flask, render_template, request, redirect, url_for
from langchain_openai import ChatOpenAI
import glob
from dotenv import load_dotenv
from helpers.rag_helper import create_document_loader, initialize_vector_store, create_template, create_graph, run_graph

load_dotenv()

global cell_subtypes
global groups
global experimental_design
global annotation_data
global expression_table


def initialize():
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')


def read_annotation_data():
    adata = sc.read_csv(glob.glob(os.path.abspath('data/*counts.csv'))[0]).T
    meta = pd.read_csv(glob.glob(os.path.abspath('data/*metadata.csv'))[0])

    if any(adata.var_names.duplicated()):
        print("Warning: Duplicate gene names found. These will be handled by combining expression.")

    # Create binary matrix
    binary_matrix = (adata.X > 0).astype(int)

    # Convert binary matrix to DataFrame with gene names
    binary_df = pd.DataFrame(binary_matrix, columns=adata.var_names, index=adata.obs_names)

    meta.index = adata.obs_names
    groups = meta['group'].unique()

    # Initialize results dictionary
    group_percentages = {}

    # Calculate percentages for each group
    for group in groups:
        # Get cells belonging to this group
        group_cells = meta['group'] == group

        # Calculate percentage of cells expressing each gene in this group
        group_percentages[group] = (binary_df[group_cells].mean() * 100)

    # Convert results to DataFrame
    result_df = pd.DataFrame(group_percentages)
    result_df["gene"] = binary_df.columns

    # Handle duplicate gene names by taking the maximum percentage
    if any(result_df.index.duplicated()):
        result_df = result_df.groupby(level=0).max()

    global annotation_data
    annotation_data = adata

    global expression_table
    expression_table = result_df

    print(result_df)

    return adata, meta


def meta_filter_and_adata_append(adata, meta):
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM2', 'Tam1')
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM1', 'Tam2')
    adata.obs = meta


def cell_filter(adata, subtype):
    cell_subtype = adata[adata.obs['Cell_subtype'] == subtype]
    print(cell_subtype)
    return cell_subtype


def control_filter(cell_subtype, control_group, knockout_group):
    Ko = knockout_group
    control = control_group

    cell_subtype.obs.index = cell_subtype.obs.index.astype(str)
    enhanced_cell_subtype = cell_subtype[(cell_subtype.obs['group'] == Ko) | (cell_subtype.obs['group'] == control), :]
    return enhanced_cell_subtype


def run_diffexp(enhanced_cell_subtype, control_group, knockout_group):
    # print("Enhanced cell subtype", enhanced_cell_subtype)
    counts = pd.DataFrame(enhanced_cell_subtype.X, columns=enhanced_cell_subtype.var_names)
    # print(counts)
    # print(enhanced_cell_subtype.obs)

    dds = DeseqDataSet(
        counts=counts,
        metadata=enhanced_cell_subtype.obs,
        design_factors="group"
    )

    sc.pp.filter_genes(dds, min_cells=1)
    dds.deseq2()

    Ko = knockout_group
    control = control_group
    stat_res = DeseqStats(dds, contrast=('group', Ko, control))

    stat_res.summary()
    results = stat_res.results_df
    return results


def process_diffexp(results):
    results = results.sort_values('stat', ascending=False)
    index_list = results.index.tolist()
    results['Gene'] = index_list
    return results


def run_gseapy(results):
    results['Rank'] = -np.log10(results.padj) * results.log2FoldChange
    results = results.sort_values('Rank', ascending=False)

    top_50 = results.head(50)

    # Take the bottom 10 rows
    bottom_50 = results.tail(50)

    # Combine the two subsets into one DataFrame
    selected_results = pd.concat([top_50, bottom_50])

    selected_results.to_csv(os.path.abspath('outputs/degs.csv'))
    results_html = results.to_html(classes='table table-bordered table-striped', index=False)

    ranking = results[['Gene', 'Rank']]
    ranking = ranking.reset_index(drop=True)
    ranking['Gene'] = ranking['Gene'].str.upper()
    pre_res = gp.prerank(rnk=ranking, gene_sets='BioPlanet_2019', seed=6)
    out = []

    for term in list(pre_res.results):
        out.append([term,
                    pre_res.results[term]['fdr'],
                    pre_res.results[term]['es'],
                    pre_res.results[term]['nes'],
                    pre_res.results[term]['lead_genes']]),

    out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes', 'lead_genes']).sort_values('nes').reset_index(
        drop=True)
    out_df = out_df[out_df['fdr'] < 0.1]
    out_df = out_df.sort_values('nes', ascending=False)
    out_df['nes'] = out_df['nes'].apply(
        lambda x: 'upregulated' if x > 1 else ('downregulated' if x < -1 else 'neutral'))
    out_df_filtered = out_df.loc[out_df['nes'] != 'neutral']
    out_df_filtered = out_df_filtered[['Term', 'nes', 'lead_genes']]

    first_15 = out_df_filtered.head(15)
    last_15 = out_df_filtered.tail(15)
    out_df_filtered = pd.concat([first_15, last_15])

    out_df_filtered.to_csv(os.path.abspath("outputs/output.csv"), index=False)
    out_df_filtered.to_csv(os.path.abspath("outputs/output.txt"), sep=',', index=False)

    global expression_table
    expression_html = expression_table.to_html(classes='table table-bordered table-striped', index=False)
    return out_df_filtered, results_html, expression_html


def ai_analysis(data, cell_subtypes, experimental_description):
    csv_file = os.path.abspath('outputs/output.csv')  # Replace with your actual CSV file path
    df = pd.read_csv(csv_file)
    # Convert DataFrame to string (if you want to display the entire CSV content as a string)
    csv_string = df.to_string(index=False)
    llm = ChatOpenAI(temperature=0, model="gpt-4o", api_key=os.environ.get("KEY"))
    docs = create_document_loader(cell_subtypes, experimental_description)
    vector_store = initialize_vector_store()
    prompt = create_template()
    graph = create_graph(vector_store, prompt, llm)
    response = run_graph(graph, csv_string, cell_subtypes, experimental_description)
    output_file = os.path.abspath('outputs/llmoutput.txt')
    with open(output_file, 'w') as file:
        file.write(response)


def filter_control(query):
    # Split the genes string into a list
    genes = query.split(",")

    # Initialize an empty list to store filtered DataFrames
    filtered = []

    # Read the DEGs CSV into a DataFrame
    degs = pd.read_csv("outputs/degs.csv")

    # Iterate through the list of genes
    for gene in genes:
        # Filter the DEGs DataFrame for each gene and append to the list
        gene_filtered = degs.loc[degs['Gene'] == gene]
        filtered.append(gene_filtered)

    # Concatenate all filtered DataFrames into one DataFrame
    result_df = pd.concat(filtered, ignore_index=True)
    filtered_df_html = result_df.to_html(classes='table table-bordered table-striped', index=False)

    return filtered_df_html
