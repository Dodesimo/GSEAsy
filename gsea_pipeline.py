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

load_dotenv()

global cell_subtypes
global groups


def initialize():
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')


def read_annotation_data():
    adata = sc.read_csv(glob.glob(os.path.abspath('data/*_counts.csv'))[0]).T
    meta = pd.read_csv(glob.glob(os.path.abspath('data/*_metadata.csv'))[0])
    return adata, meta


def meta_filter_and_adata_append(adata, meta):
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM2', 'Tam1')
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM1', 'Tam2')
    adata.obs = meta


def cell_filter(adata, subtype):
    cell_subtype = adata[adata.obs['Cell_subtype'] == subtype]
    print(cell_subtype)
    return cell_subtype


def control_filter(cell_subtype, gene):

    Ko = 'KO.' + gene
    control = 'control.' + gene

    cell_subtype.obs.index = cell_subtype.obs.index.astype(str)
    enhanced_cell_subtype = cell_subtype[(cell_subtype.obs['group'] == Ko) | (cell_subtype.obs['group'] == control), :]
    return enhanced_cell_subtype


def run_diffexp(enhanced_cell_subtype, gene):
    #print("Enhanced cell subtype", enhanced_cell_subtype)
    counts = pd.DataFrame(enhanced_cell_subtype.X, columns=enhanced_cell_subtype.var_names)
    #print(counts)
    #print(enhanced_cell_subtype.obs)

    dds = DeseqDataSet(
        counts=counts,
        metadata=enhanced_cell_subtype.obs,
        design_factors="group"
    )

    sc.pp.filter_genes(dds, min_cells=1)
    dds.deseq2()

    Ko = 'KO.' + gene
    control = 'control.' + gene
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
    ranking = results[['Gene', 'Rank']]
    ranking = ranking.reset_index(drop=True)
    ranking['Gene'] = ranking['Gene'].str.upper()
    pre_res = gp.prerank(rnk=ranking, gene_sets='BioPlanet_2019', seed=6)

    out = []

    for term in list(pre_res.results):
        out.append([term,
                    pre_res.results[term]['fdr'],
                    pre_res.results[term]['es'],
                    pre_res.results[term]['nes']])
    out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes']).sort_values('nes').reset_index(drop=True)
    out_df = out_df[out_df['fdr'] < 0.1]
    out_df = out_df.sort_values('nes', ascending=False)
    out_df['nes'] = out_df['nes'].apply(
        lambda x: 'upregulated' if x > 1 else ('downregulated' if x < -1 else 'neutral'))
    out_df_filtered = out_df.loc[out_df['nes'] != 'neutral']
    out_df_filtered = out_df_filtered[['Term', 'nes']]
    out_df_filtered.to_csv(os.path.abspath("outputs/output.csv"), index=False)
    out_df_filtered.to_csv(os.path.abspath("outputs/output.txt"), sep=',', index=False)
    return out_df_filtered

def ai_analysis(data, cell_subtypes):
    csv_file = os.path.abspath('outputs/output.csv')  # Replace with your actual CSV file path
    df = pd.read_csv(csv_file)

    # Convert DataFrame to string (if you want to display the entire CSV content as a string)
    csv_string = df.to_string(index=False)
    llm = ChatOpenAI(temperature=0, model="gpt-4o", api_key=os.getenv("KEY"))
    response = llm.invoke(
        f"Given the following CSV of gene labels and whether they are upregulated or downregulated: \n\n{csv_string}\n\n, as well as that the cell type in question is \n\n{cell_subtypes}\n\n, and that TIM3 is knocked out, give a breakdown of biological significance by gene sets enriched and pose a potential mechanism. Write a coherent Markdown paragraph with no title.")
    print(response)
    response = str(response.content)
    output_file = os.path.abspath('outputs/llmoutput.txt')
    with open(output_file, 'w') as file:
        file.write(response)
