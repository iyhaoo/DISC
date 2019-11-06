import scscope as DeepImpute
import scanpy.api as sc
import numpy as np
import h5py
import pandas as pd
import glob
import argparse
import os
import time
from collections import Counter


def read_loom(loom_path):
    assert os.path.exists(loom_path)
    try:
        with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
            gene_name = f["row_attrs/Gene"][...].astype(np.str)
            assert np.max(list(Counter(gene_name).values())) == 1, "{}".format(list(filter(lambda x: x[1] > 1, Counter(gene_name).items())))
            gene_bc_mat = f["matrix"][...]
            cell_id = f["col_attrs/CellID"][...].astype(np.str) if "col_attrs/CellID" in f.keys() else np.arange(gene_bc_mat.shape[1]).astype(np.str)
        return gene_bc_mat, cell_id, gene_name
    except OSError:
        print("{}\nis connecting by other process, waiting...".format(loom_path))
        time.sleep(5)
        return read_loom(loom_path)


parser = argparse.ArgumentParser()
parser.add_argument('--filt-loom', required=True, type=str, help="Filtered data")
parser.add_argument('--num-gpus', required=False, type=int, default=1, help="how many gpus to use")
parser.add_argument('--epoch', required=False, type=int, default=100, help="how many epochs to run")
parser.add_argument('--min-expressed-cell', required=False, type=int, default=10, help="min-expressed-cell")
parser.add_argument("--min-expressed-cell-average-expression", required=False, type=float, default=1, help="min-expressed-cell-average-expression")
FLAGS = vars(parser.parse_args())
output_dir = "{}/imputation".format(FLAGS["filt_loom"].rsplit("/", 1)[0])
os.makedirs(output_dir, exist_ok=True)

starttime = time.time()
gene_bc_mat, cell_id, gene_name = read_loom(FLAGS["filt_loom"])
min_expressed_cell = FLAGS["min_expressed_cell"]
min_expressed_cell_average_expression = FLAGS["min_expressed_cell_average_expression"]

input_loom_name = FLAGS["filt_loom"].rsplit("/", 1)[1]
if FLAGS["epoch"] == 100:
    output_raw_h5 = input_loom_name.replace(".loom", "_scScope_mc_{}_mce_{}_raw.hdf5".format(min_expressed_cell, min_expressed_cell_average_expression))
    output_feature_h5 = input_loom_name.replace(".loom", "_scScope_mc_{}_mce_{}_feature.hdf5".format(min_expressed_cell, min_expressed_cell_average_expression))
    output_filtered_h5 = input_loom_name.replace(".loom", "_scScope_mc_{}_mce_{}.hdf5".format(min_expressed_cell, min_expressed_cell_average_expression))
else:
    output_raw_h5 = input_loom_name.replace(".loom", "_scScope_{}epochs_mc_{}_mce_{}_raw.hdf5".format(FLAGS["epoch"], min_expressed_cell, min_expressed_cell_average_expression))
    output_feature_h5 = input_loom_name.replace(".loom", "_scScope_{}epochs_mc_{}_mce_{}_feature.hdf5".format(FLAGS["epoch"], min_expressed_cell, min_expressed_cell_average_expression))
    output_filtered_h5 = input_loom_name.replace(".loom", "_scScope_{}epochs_mc_{}_mce_{}.hdf5".format(FLAGS["epoch"], min_expressed_cell, min_expressed_cell_average_expression))
expressed_cell = (gene_bc_mat > 0).sum(1)
gene_expression = gene_bc_mat.sum(1)
gene_filter = np.bitwise_and(expressed_cell >= min_expressed_cell, gene_expression > expressed_cell * min_expressed_cell_average_expression)
input_bc_gene_mat = gene_bc_mat[gene_filter, :].transpose()
print(input_bc_gene_mat.shape)
raw_library_size = input_bc_gene_mat.sum(0)
filt_adata = sc.AnnData(input_bc_gene_mat)
# 2. Normalize gene expression based on scanpy (normalize each cell to have same library size)
# normalize each cell to have same count number
sc.pp.normalize_per_cell(filt_adata)
# update datastructure to use normalized data
gene_expression_norm = filt_adata.X
latent_dim = 50
# 3. scScope learning
DI_model = DeepImpute.train(gene_expression_norm, latent_dim, num_gpus=FLAGS["num_gpus"], max_epoch=FLAGS["epoch"])

# 4. latent representations and imputed expressions
latent_code, imputed_val, predicted_batch_effect = DeepImpute.predict(gene_expression_norm, DI_model)
with h5py.File("{}/{}".format(output_dir, output_feature_h5), "w") as f:
    f["cell_id"] = cell_id.astype(h5py.special_dtype(vlen=str))
    ff_dset_feature = f.create_dataset("feature", shape=(cell_id.size, latent_code.shape[1]), chunks=(1, latent_code.shape[1]), dtype=np.float32)
    ff_dset_feature[...] = latent_code

imputed_val_resume = imputed_val * raw_library_size / filt_adata.X.sum(1, keepdims=True)
with h5py.File("{}/{}".format(output_dir, output_raw_h5), "w") as f:
    f["cell_id"] = cell_id.astype(h5py.special_dtype(vlen=str))
    f["gene_name"] = gene_name[gene_filter].astype(h5py.special_dtype(vlen=str))
    if_dset_imputation = f.create_dataset("imputation", shape=(cell_id.size, gene_filter.sum()), chunks=(1, gene_filter.sum()), dtype=np.float32)
    if_dset_imputation[...] = imputed_val_resume

zero_entries = np.bitwise_not((input_bc_gene_mat != 0))
with h5py.File("{}/{}".format(output_dir, output_filtered_h5), "w") as f:
    f["cell_id"] = cell_id.astype(h5py.special_dtype(vlen=str))
    f["gene_name"] = gene_name[gene_filter].astype(h5py.special_dtype(vlen=str))
    if_dset_imputation = f.create_dataset("imputation", shape=(cell_id.size, gene_filter.sum()), chunks=(1, gene_filter.sum()), dtype=np.float32)
    if_dset_imputation[...] = imputed_val_resume * zero_entries + input_bc_gene_mat

print(pd.to_timedelta(time.time() - starttime, unit="s"))


