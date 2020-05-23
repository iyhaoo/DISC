import scipy.sparse as sp_sparse
import tables
import numpy as np
import loompy
import pandas as pd



def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            matrix = sp_sparse.csc_matrix((dsets['data'], dsets['indices'], dsets['indptr']), shape=dsets['shape'])
            return dsets['genes'], dsets['gene_names'], dsets['barcodes'], matrix
        except tables.NoSuchNodeError:
            raise Exception("Genome %s does not exist in this file." % genome)
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


filtered_matrix_h5 = "/home/yuanhao/data/fn/neuron1.3m/1M_neurons_filtered_gene_bc_matrices_h5.h5"
genome = "mm10"
gene_id, gene_name, cell_id, gene_bc_sparse = get_matrix_from_h5(filtered_matrix_h5, genome)
row_attrs = {"Gene": gene_name.astype(np.str)}
col_attrs = {"CellID": cell_id.astype(np.str), "SampleID": pd.Series(cell_id.astype(np.str)).str.split("-").str.get(1).values}
loompy.create("/home/yuanhao/data/fn/neuron1.3m/1M_neurons.loom", gene_bc_sparse, row_attrs, col_attrs)




