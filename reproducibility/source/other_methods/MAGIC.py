import magic
import numpy as np
import h5py
import pandas as pd
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
parser.add_argument('--loom', required=True, type=str, help="Loom")
parser.add_argument('--min-expressed-cell', required=False, type=int, default=10, help="min-expressed-cell")
parser.add_argument("--min-expressed-cell-average-expression", required=False, type=float, default=1, help="min-expressed-cell-average-expression")
FLAGS = vars(parser.parse_args())
output_dir = "{}/imputation".format(FLAGS["loom"].rsplit("/", 1)[0])
os.makedirs(output_dir, exist_ok=True)

starttime = time.time()
gene_bc_mat, cell_id, gene_name = read_loom(FLAGS["loom"])
min_expressed_cell = FLAGS["min_expressed_cell"]
min_expressed_cell_average_expression = FLAGS["min_expressed_cell_average_expression"]
expressed_cell = (gene_bc_mat > 0).sum(1)
gene_expression = gene_bc_mat.sum(1)
gene_filter = np.bitwise_and(expressed_cell >= min_expressed_cell, gene_expression > expressed_cell * min_expressed_cell_average_expression)
input_gene_bc_mat = gene_bc_mat[gene_filter, :]
print(input_gene_bc_mat.shape)
bc_gene_filtered = input_gene_bc_mat.transpose()
raw_library_size = np.asarray(bc_gene_filtered.sum(0)).astype(np.float32)
bc_gene_norm = bc_gene_filtered / raw_library_size * np.median(raw_library_size.squeeze())
bc_gene_norm_sqrt = np.sqrt(bc_gene_norm)
pd_input_norm_sqrt = pd.DataFrame(bc_gene_norm_sqrt, columns=gene_name[gene_filter])
magic_operator = magic.MAGIC()
impute_norm_sqrt = magic_operator.fit_transform(bc_gene_norm_sqrt)
impute_norm = np.square(impute_norm_sqrt)
input_loom_name = FLAGS["loom"].rsplit("/", 1)[1]
output_h5 = input_loom_name.replace(".loom", "_MAGIC_mc_{}_mce_{}.hdf5".format(min_expressed_cell, min_expressed_cell_average_expression))
with h5py.File("{}/{}".format(output_dir, output_h5), "w") as f:
    f["cell_id"] = cell_id.astype(h5py.special_dtype(vlen=str))
    f["gene_name"] = gene_name[gene_filter].astype(h5py.special_dtype(vlen=str))
    if_dset_imputation = f.create_dataset("imputation", shape=(cell_id.size, gene_filter.sum()), chunks=(1, gene_filter.sum()), dtype=np.float32)
    if_dset_imputation[...] = impute_norm * raw_library_size / np.median(raw_library_size.squeeze())
print(pd.to_timedelta(time.time() - starttime, unit="s"))
