import scanpy as sc
import pandas as pd
import warnings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
import numpy as np

def initialize():
    sc.settings.verbosity = 3
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

def read_annotation_data(ann_path, meta_path):
    adata = sc.read_csv(ann_path).T
    meta = pd.read_csv(meta_path)
    return adata, meta

def meta_filter_and_adata_append(adata, meta):
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM2', 'Tam1')
    meta['Cell_subtype'] = meta['Cell_subtype'].replace('TAM1', 'Tam2')
    adata.obs = meta

def cell_filter(adata, subtype):
    cell_subtype = adata[adata.obs['Cell_subtype'] == subtype]
    return cell_subtype

def control_filter(cell_subtype, adata, gene):
    enhanced_cell_subtype = cell_subtype[
                      (cell_subtype.obs['group'] == f'KO.{gene}') | (cell_subtype.obs['group'] == f'control.{gene}'), :]
    return enhanced_cell_subtype

def run_diffexp(enhanced_cell_subtype, gene):
    counts = pd.DataFrame(enhanced_cell_subtype.X, columns=enhanced_cell_subtype.var_names)

    dds = DeseqDataSet(
        counts=counts,
        metadata=enhanced_cell_subtype.obs,
        design_factors="group"
    )

    sc.pp.filter_genes(dds, min_cells=1)
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast = ('group', f'KO.{gene}', f'control.{gene}'))

    stat_res.summary()
    results = stat_res.results_df
    return results

def process_diffexp(results):
    results = results.sort_values('stat', ascending = False)
    index_list = results.index.tolist()
    results['Gene'] = index_list

def run_gseapy(results):
    results['Rank'] = -np.log10(results.padj)*results.log2FoldChange
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
    out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop=True)
    out_df
    return out_df








