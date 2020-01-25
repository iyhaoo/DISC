import loompy
from deepimpute.multinet import MultiNet
import numpy as np
import h5py
import pandas as pd
import glob
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--filt-loom', required=True, type=str, help="Filtered data")
parser.add_argument('--n-cores', required=False, type=int, default=1, help="Cores use")
FLAGS = vars(parser.parse_args())
output_dir = "{}/imputation".format(FLAGS["filt_loom"].rsplit("/", 1)[0])
os.makedirs(output_dir, exist_ok=True)

starttime = time.time()
with loompy.connect(FLAGS["filt_loom"]) as ds:
    gene_bc_sparse = ds.sparse().astype(np.float32).tocsr()
    gene_name = ds.ra["Gene"]
    cell_id = ds.ca["CellID"]

min_expressed_cell = 10
min_expressed_cell_average_expression = 1
expressed_cell = np.asarray((gene_bc_sparse > 0).sum(1)).squeeze()
gene_expression = np.asarray(gene_bc_sparse.sum(1)).squeeze()
gene_filter = np.bitwise_and(expressed_cell >= min_expressed_cell, gene_expression > expressed_cell * min_expressed_cell_average_expression)
input_gene_bc_mat = gene_bc_sparse[gene_filter, :]
print(input_gene_bc_mat.shape)
data = pd.DataFrame(input_gene_bc_mat.toarray()).T
multinet = MultiNet(n_cores=FLAGS["n_cores"])
multinet.fit(data, cell_subset=1)
imputedData = multinet.predict(data)
input_loom_name = FLAGS["filt_loom"].rsplit("/", 1)[1]
output_h5 = input_loom_name.replace(".loom", "_deepImpute_mc_{}_mce_{}.hdf5".format(min_expressed_cell, min_expressed_cell_average_expression))
with h5py.File("{}/{}".format(output_dir, output_h5), "w") as f:
    f["cell_id"] = cell_id.astype(h5py.special_dtype(vlen=str))
    f["gene_name"] = gene_name[gene_filter].astype(h5py.special_dtype(vlen=str))
    if_dset_imputation = f.create_dataset("imputation", shape=(cell_id.size, gene_filter.sum()), chunks=(1, gene_filter.sum()), dtype=np.float32)
    if_dset_imputation[...] = imputedData.values
print(pd.to_timedelta(time.time() - starttime, unit="s"))





























