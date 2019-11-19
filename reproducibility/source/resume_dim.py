import numpy as np
import pandas as pd
import os
import argparse
import shutil
import h5py

parser = argparse.ArgumentParser()
parser.add_argument('--raw-loom', required=True, type=str, help="loom")
parser.add_argument('--impute-h5', required=True, type=str, help="h5")
FLAGS = vars(parser.parse_args())

loom_path = FLAGS["raw_loom"]
h5_path = FLAGS["impute_h5"]

output_h5 = "{}_resume_dim.loom".format(h5_path.rsplit(".", 1)[0])
output_dir = output_h5.rsplit("/", 1)[0]
os.makedirs(output_dir, exist_ok=True)

input_format = h5_path.rsplit(".", 1)[1]
if input_format == "loom":
    with h5py.File(h5_path, "r", libver='latest', swmr=True) as f:
        impute_pd = pd.DataFrame(f["matrix"][...], index=f["row_attrs/Gene"][...].astype(np.str), columns=f["col_attrs/CellID"][...].astype(np.str))
else:
    with h5py.File(h5_path, "r", libver='latest', swmr=True) as f:
        impute_pd = pd.DataFrame(f["imputation"][...].transpose(), index=f["gene_name"][...].astype(np.str), columns=f["cell_id"][...].astype(np.str))

with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
    gene_bc_mat = f["matrix"][...].astype(np.float32)
    cell_id = f["col_attrs/CellID"][...].astype(np.str)
    gene_name = f["row_attrs/Gene"][...].astype(np.str)
    raw_cell_mask = np.isin(cell_id, impute_pd.columns.values)
    raw_gene_mask = np.isin(gene_name, impute_pd.index.values)
    object_cell = cell_id[raw_cell_mask]
    object_gene = gene_name[raw_gene_mask]
    gene_bc_mat_resume_cell = gene_bc_mat[raw_gene_mask, :]
    gene_bc_mat_resume_cell[:, raw_cell_mask] = impute_pd.reindex(object_cell, axis="columns").reindex(object_gene).values
    gene_bc_mat[raw_gene_mask, :] = gene_bc_mat_resume_cell
    #print(gene_bc_mat)
    print(np.isnan(gene_bc_mat).sum())
    gene_bc_mat[np.isnan(gene_bc_mat)] = 0
shutil.copy(loom_path, output_h5)
with h5py.File(output_h5) as f:
    del f["matrix"]
    f_matrix = f.create_dataset("matrix", shape=gene_bc_mat.shape, chunks=(gene_bc_mat.shape[0], 1), dtype=np.float32)
    f_matrix[...] = gene_bc_mat

