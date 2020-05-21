import loompy
import numpy as np
import seaborn as sns
import pandas as pd


loom_path = "/home/yuanhao/data/fn/neuron1.3m/1M_neurons.loom"
with loompy.connect(loom_path) as ds:
    gene_bc_sparse = ds.sparse().tocsc()
    gene_name = ds.ra["Gene"]
    cell_id = ds.ca["CellID"]

g = sns.violinplot(np.log1p(np.asarray(gene_bc_sparse.sum(0)).squeeze()))
g.figure.savefig("/home/yuanhao/data/fn/neuron1.3m/library_size.pdf")

expressed_cell = np.asarray((gene_bc_sparse > 0).sum(1)).squeeze()
gene_expression = np.asarray(gene_bc_sparse.sum(1)).squeeze()
gene_expression_number = np.asarray((gene_bc_sparse > 0).sum(0)).squeeze()
cell_filter = np.bitwise_and(gene_expression_number >= 500, gene_expression_number <= 5000)
gene_bc_filt = gene_bc_sparse[:, cell_filter]

row_attrs = {"Gene": gene_name}
col_attrs = {"CellID": cell_id[cell_filter], "SampleID": pd.Series(cell_id[cell_filter]).str.split("-").str.get(1).values}
loompy.create("/home/yuanhao/data/fn/neuron1.3m/1M_neurons_filt_gene_500_5000.loom", gene_bc_filt, row_attrs, col_attrs)



