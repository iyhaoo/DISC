import numpy as np
import pandas as pd
import h5py
import rpy2.robjects
from rpy2.robjects.pandas2ri import ri2py
readRDS = rpy2.robjects.r["readRDS"]
import os
from multiprocessing import Pool
import collections
import argparse
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
sns.set(context='paper', style='white')


def violin_plot(gene_bc_pd, cell_type_series, output_path, cluster_order=False, sort_group=False):
    plot_gene_number = gene_bc_pd.shape[0]
    cluster_result_np = cell_type_series.values.squeeze()
    cell_number_pd = pd.Series(1).repeat(cell_type_series.size).groupby(cluster_result_np, sort=sort_group).sum()
    if cluster_order:
        cell_number_pd = cell_number_pd[np.argsort([int(x.split(" ", 1)[0]) for x in cell_number_pd.index.values.tolist()])]
    sns_pd = pd.DataFrame()
    sns_pd.loc[:, "Gene Symbol"] = gene_bc_pd.index.str.capitalize().repeat(cell_type_series.size)
    sns_pd.loc[:, "Known Cell Types"] = cluster_result_np.tolist() * plot_gene_number
    sns_pd.loc[:, "log expression"] = np.log1p(np.ravel(gene_bc_pd.loc[:, cell_type_series.index].values))
    sns_pd.loc[:, "expression"] = np.ravel(gene_bc_pd.loc[:, cell_type_series.index].values)
    height = 3 + np.unique(cluster_result_np).size * 0.25
    aspect = np.maximum(0.35 - plot_gene_number * 0.025, 1.5 / np.unique(cluster_result_np).size)
    g = sns.FacetGrid(sns_pd, col='Gene Symbol', height=height, aspect=aspect, gridspec_kws=dict(hspace=0, wspace=0), sharex=False)
    g.map(sns.violinplot, 'log expression', 'Known Cell Types', orient='h', scale='width', order=cell_number_pd.index, palette='husl', inner=None, cut=False)
    g.set(xlabel="")
    g.set_xticklabels("")
    g.set_yticklabels(["{}{}{}".format(x, "\n", y) for x, y in zip(cell_number_pd.index, cell_number_pd)], ha='right')
    g.set_titles("")
    g.set_titles("{col_name}", rotation=40, loc='left', va='bottom')
    g.savefig(output_path)


def dataset_scanning(loom_path, gene_index, ix, batch_size, with_outliers=True):
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        input_format = loom_path.rsplit(".", 1)[1]
        if input_format == "loom":
            cell_number = f["col_attrs/CellID"].shape[0]
            this_batch = f["matrix"][:, np.arange(ix, np.minimum(ix + batch_size, cell_number), dtype=np.int32)].astype(np.float32)
        else:
            cell_number = f["cell_id"].shape[0]
            this_batch = f["imputation"][np.arange(ix, np.minimum(ix + batch_size, cell_number), dtype=np.int32), :].astype(np.float32).transpose()
            if with_outliers:
                this_batch = this_batch + f["outliers"][np.arange(ix, np.minimum(ix + batch_size, cell_number), dtype=np.int32), :].astype(np.float32).transpose()
        this_library_size = np.sum(this_batch, 0)
        this_norm = np.divide(this_batch, this_library_size, out=np.zeros_like(this_batch), where=this_library_size != 0)[gene_index, :]
    return this_norm


def read_loom(loom_path, gene_index, scanning_thread_number=None, scanning_batch_size=2048):
    assert os.path.exists(loom_path)
    f = h5py.File(loom_path, "r", libver='latest', swmr=True)
    input_format = loom_path.rsplit(".", 1)[1]
    if input_format == "loom":
        cell_number = f["col_attrs/CellID"].shape[0]
    else:
        cell_number = f["cell_id"].shape[0]
    f.close()
    result_list = []
    scan_pool = Pool(processes=scanning_thread_number)
    for ii in np.arange(np.ceil(cell_number / scanning_batch_size)).tolist():
        result_list.append(
            scan_pool.apply_async(dataset_scanning, args=(loom_path, gene_index, ii * scanning_batch_size, scanning_batch_size, True)))
    scan_pool.close()
    scan_pool.join()
    this_norm_list = []
    for ii in result_list:
        this_norm = ii.get()
        this_norm_list.append(this_norm)
    gene_bc_norm = np.hstack(this_norm_list)
    return gene_bc_norm


if __name__ == "__main__":
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5', required=True, type=str, help="resume hdf5 path")
    parser.add_argument('--cell-type', required=True, type=str, help="RDS formatted file")
    parser.add_argument('--gene-name', required=True, type=str, help="gene name (split by \",\") or set name (PBMC, RETINA, BRAIN_SPLiT_ALL, BRAIN_SPLiT_N, BRAIN_SPLiT_NN)")
    parser.add_argument('--output-path', required=False, type=str, help="output path")
    FLAGS = vars(parser.parse_args())
    """
    """
    FLAGS = {"h5": "E:/DISC/reproducibility/data/PBMC/raw.loom",
             "cell_type": "E:/DISC/reproducibility/results/PBMC/cell_type.rds",
             "gene_name": "PBMC",
             "output_path": "E:/DISC/reproducibility/results/PBMC/violin_plot.pdf"}
    """
    FLAGS = {"h5": "E:/DISC/reproducibility/data/RETINA/raw.loom",
             "cell_type": "E:/DISC/reproducibility/data/RETINA/cell_type.rds",
             "gene_name": "RETINA",
             "output_path": "E:/DISC/reproducibility/results/RETINA/violin_plot.pdf"}
    cell_type_rds = readRDS(FLAGS["cell_type"])
    cell_type_series = pd.Series(ri2py(cell_type_rds), index=cell_type_rds.names)
    h5_path = FLAGS["h5"]
    gene_name = FLAGS["gene_name"]
    cluster_marker_n = ["DEPTOR", "RARB", "SATB2", "TFAP2D", "FIGN", "ARAP1", "PAX3", "NTN1", "PAX2", "SLC6A3", "FN1", "TSPAN18", "PDE11A", "DLX6OS1"]
    cluster_marker_nn = ["MBP", "PDGFRA", "DOCK2", "RGS5", "COL1A2", "ALDH1L1", "DNAH11", "MYBPC1"]
    if gene_name == "PBMC":
        cluster_marker = ["IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP"]
    elif gene_name == "RETINA":
        cluster_marker = ["LHX1", "SLC17A6", "PAX6", "GAD1", "SLC6A9", "OPN1MW", "VSX2", "RLBP1", "GFAP", "PECAM1", "KCNJ8", "CX3CR1"]
    elif gene_name == "BRAIN_SPLiT_ALL":
        cluster_marker = cluster_marker_n + cluster_marker_nn
    elif gene_name == "BRAIN_SPLiT_N":
        cluster_marker = cluster_marker_n
    elif gene_name == "BRAIN_SPLiT_NN":
        cluster_marker = cluster_marker_nn
    else:
        cluster_marker = gene_name.split(",")

    if FLAGS["output_path"] is None:
        output_path = "{}/violin_plot_{}.pdf".format(h5_path.rsplit("/", 1)[0], gene_name)
    else:
        output_path = FLAGS["output_path"]

    input_format = h5_path.rsplit(".", 1)[1]
    if input_format == "loom":
        with h5py.File(h5_path, "r") as f:
            cell_id = f["col_attrs/CellID"][...].astype(np.str)
            gene_name = f["row_attrs/Gene"][...].astype(np.str)
    else:
        with h5py.File(h5_path, "r") as f:
            cell_id = f["cell_id"][...].astype(np.str)
            gene_name = f["gene_name"][...].astype(np.str)

    assert np.max(list(collections.Counter(gene_name).values())) == 1, "{}".format(list(filter(lambda x: x[1] > 1, collections.Counter(gene_name).items())))
    gene_index, = np.where(np.isin(gene_name, cluster_marker))
    gene_bc_norm = read_loom(h5_path, gene_index)
    use_cell = np.intersect1d(cell_id, cell_type_series.index.values)
    cell_type_filt_pd = cell_type_series.loc[use_cell]
    filt_norm_pd = pd.DataFrame(gene_bc_norm * 10000, index=gene_name[gene_index], columns=cell_id)
    violin_plot(filt_norm_pd.loc[cluster_marker, use_cell], cell_type_filt_pd, output_path)
