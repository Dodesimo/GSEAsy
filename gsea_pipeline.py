import scanpy as sc
import pandas as pd
import warnings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
import numpy as np
from flask import Flask, render_template, request, redirect, url_for
from langchain_openai import ChatOpenAI
from langchain_community.tools import AIPluginTool

def initialize():
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def read_annotation_data():
    adata = sc.read_csv('/Users/atharvainamdar/Desktop/GSEAwebsite /data/TIM3_counts.csv').T
    meta = pd.read_csv('/Users/atharvainamdar/Desktop/GSEAwebsite /data/TIM3_metadata.csv')
    return adata, meta

def meta_filter_and_adata_append(adata, meta):
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM2', 'Tam1')
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM1', 'Tam2')
    adata.obs = meta

def cell_filter(adata):
    with open('/Users/atharvainamdar/Desktop/GSEAwebsite /GSEAwebsite/inputs/cellsubtype.txt', 'r') as file:
        subtype = file.read()
    print(subtype)
    cell_subtype = adata[adata.obs['Cell_subtype'] == subtype]

    print(cell_subtype.shape)
    print(cell_subtype.obs)
    return cell_subtype


def control_filter(cell_subtype):
    with open('/Users/atharvainamdar/Desktop/GSEAwebsite /GSEAwebsite/inputs/group.txt', 'r') as file:
        gene = file.read()
    print(gene)
    Ko = 'KO.' + gene
    control = 'control.' + gene
    print(Ko)
    print(control)
    print(cell_subtype.obs['group'].dtype)  # Should be categorical
    print(cell_subtype.obs['group'].unique())  # Check unique values
    cell_subtype.obs.index = cell_subtype.obs.index.astype(str)
    enhanced_cell_subtype = cell_subtype[(cell_subtype.obs['group'] == Ko) | (cell_subtype.obs['group'] == control), :]
    print(enhanced_cell_subtype.shape)
    return enhanced_cell_subtype

def run_diffexp(enhanced_cell_subtype):
    counts = pd.DataFrame(enhanced_cell_subtype.X, columns=enhanced_cell_subtype.var_names)

    dds = DeseqDataSet(
        counts=counts,
        metadata=enhanced_cell_subtype.obs,
        design_factors="group"
    )

    sc.pp.filter_genes(dds, min_cells=1)
    dds.deseq2()
    with open('/Users/atharvainamdar/Desktop/GSEAwebsite /GSEAwebsite/inputs/group.txt', 'r') as file:
        gene = file.read()
    print(gene)
    Ko = 'KO.' + gene
    control = 'control.' + gene
    stat_res = DeseqStats(dds, contrast = ('group', Ko, control))

    stat_res.summary()
    results = stat_res.results_df
    return results

def process_diffexp(results):
    results = results.sort_values('stat', ascending = False)
    index_list = results.index.tolist()
    results['Gene'] = index_list
    return results

def run_gseapy(results):
    results['Rank'] = -np.log10(results.padj)*results.log2FoldChange
    results = results.sort_values('Rank', ascending = False)
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
    out_df = out_df.sort_values('nes', ascending = False)
    out_df['nes'] = out_df['nes'].apply(lambda x: 'upregulated' if x > 1 else ('downregulated' if x < -1 else 'neutral'))
    out_df_filtered = out_df.loc[out_df['nes'] != 'neutral']
    out_df_filtered = out_df_filtered[['Term', 'nes']]
    out_df_filtered.to_csv('/Users/atharvainamdar/Desktop/GSEAwebsite /outputs/output.csv', index=False)
    out_df_filtered.to_csv('/Users/atharvainamdar/Desktop/GSEAwebsite /outputs/output.txt', sep=',', index=False)
    return out_df_filtered




def getdata(): 
    cell_subtypes = request.form['cell-subtypes']
    group = request.form['groups']
    # Save the text data to a file
    with open('GSEAwebsite/inputs/cellsubtype.txt', 'a') as f:
        f.write(cell_subtypes)  # Add newline after each submission


    with open('GSEAwebsite/inputs/group.txt', 'a') as f:
        f.write(group)  # Add newline after each submission
    # Redirect to home or confirmation page after submission


def ai_analysis(data): 
    csv_file = '/Users/atharvainamdar/Desktop/GSEAwebsite /outputs/output.csv'  # Replace with your actual CSV file path
    df = pd.read_csv(csv_file)
        
    # Convert DataFrame to string (if you want to display the entire CSV content as a string)
    csv_string = df.to_string(index=False)
    with open('/Users/atharvainamdar/Desktop/GSEAwebsite /GSEAwebsite/inputs/cellsubtype.txt', 'r') as file:
        subtype = file.read()
    llm = ChatOpenAI(temperature=0, model="gpt-4o", api_key = "")
    response = llm.invoke(f"Given the following CSV of gene labels and whether they are upregulated or downregulated: \n\n{csv_string}\n\n, as well as that the cell type in question is \n\n{subtype}\n\n, and that TIM3 is knocked out, give a breakdown of biological signifiance by gene sets enriched and pose a potential mechansim ")
    print(response)
    response = str(response)
    output_file = '/Users/atharvainamdar/Desktop/GSEAwebsite /GSEAwebsite/static/llmoutput.txt'
    with open(output_file, 'w') as file:
        file.write(response)
        

    