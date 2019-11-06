import numpy as np
import pandas as pd
import h5py
import os

input_h5 = "/home/yuanhao/data/fn/neuron1.3m/1M_neurons_filt_gene_500_5000_unique_rename.loom"
hvg_pd = pd.read_csv("/home/yuanhao/data/fn/neuron1.3m/1M_neurons_filt_gene_500_5000_unique_rename_hvg.tsv", sep="\t")
dest_dir = "/home/yuanhao/data/fn/neuron1.3m/performance_test_set"
top_10k_gene = hvg_pd.index.values[:10000]
top_1k_gene = hvg_pd.index.values[:1000]
chunk_size = 10000
os.makedirs(dest_dir, exist_ok=True)


def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])


with h5py.File(input_h5, "r", libver='latest', swmr=True) as f:
    gene_name = f["row_attrs/Gene"][...].astype(np.str)
    cell_id = f["col_attrs/CellID"][...].astype(np.str)
    cell_number = f["matrix"].shape[1]
    top_10k_index, = np.where(np.isin(gene_name, top_10k_gene))
    top_1k_index, = np.where(np.isin(gene_name, top_1k_gene))
target_cell_list = [50000, 100000, 500000, cell_number, cell_number * 2]

for ii in target_cell_list:
    target_cell = ii if ii != "all" else cell_number
    chunk_target_cell = np.round(target_cell * chunk_size / cell_number).astype(np.int)
    chunk_number, chunk_residue = np.divmod(cell_number, chunk_size)
    chunk_index = np.arange(int(chunk_number), dtype=np.int32)
    out_f_10k = h5py.File("{}/1M_neurons_{}_10k.loom".format(dest_dir, human_format(ii)), "w")
    out_f_10k.create_group("row_graphs")
    out_f_10k.create_group("col_graphs")
    out_f_10k.create_group("layers")
    out_f_10k["row_attrs/Gene"] = top_10k_gene.astype(np.string_)
    out_f_10k.create_dataset("matrix", shape=(top_10k_gene.size, target_cell),
                             chunks=(top_10k_gene.size, 1), dtype=np.float32, fillvalue=0, fletcher32=False,
                             compression="gzip", shuffle=False, compression_opts=2)
    out_f_1k = h5py.File("{}/1M_neurons_{}_1k.loom".format(dest_dir, human_format(ii)), "w")
    out_f_1k.create_group("row_graphs")
    out_f_1k.create_group("col_graphs")
    out_f_1k.create_group("layers")
    out_f_1k["row_attrs/Gene"] = top_1k_gene.astype(np.string_)
    out_f_1k.create_dataset("matrix", shape=(top_1k_gene.size, target_cell),
                            chunks=(top_1k_gene.size, 1), dtype=np.float32, fillvalue=0, fletcher32=False,
                            compression="gzip", shuffle=False, compression_opts=2)
    sample_cell = 0
    index_chunk = np.arange(chunk_size)
    sample_cell_id = []
    with h5py.File(input_h5, "r", libver='latest', swmr=True) as f:
        for index, jj in enumerate(chunk_index):
            repeat = 0
            use_chunk_target_cell = chunk_target_cell.copy()
            if index == chunk_index.size - 1:
                index_chunk = np.arange(chunk_size + chunk_residue, dtype=np.int32)
                use_chunk_target_cell = target_cell - sample_cell
            this_chunk = f["matrix"][:, chunk_size * index + index_chunk]
            if use_chunk_target_cell > index_chunk.size:
                repeat += 1
                use_chunk_target_cell -= index_chunk.size
            this_sample = np.sort(np.random.choice(index_chunk, use_chunk_target_cell, replace=False))
            sample_cell += this_sample.size
            sample_cell_id.append(cell_id[chunk_size * index + this_sample].astype(np.string_))
            out_f_10k["matrix"][:, (sample_cell - this_sample.size):sample_cell] = this_chunk[top_10k_index, :][:, this_sample]
            out_f_1k["matrix"][:, (sample_cell - this_sample.size):sample_cell] = this_chunk[top_1k_index, :][:, this_sample]
            for kk in range(repeat):
                sample_cell += index_chunk.size
                sample_cell_id.append(cell_id[chunk_size * index + index_chunk].astype(np.string_))
                out_f_10k["matrix"][:, (sample_cell - index_chunk.size):sample_cell] = this_chunk[top_10k_index, :]
                out_f_1k["matrix"][:, (sample_cell - index_chunk.size):sample_cell] = this_chunk[top_1k_index, :]
            print("index: {}\tsample: {}".format(index, sample_cell))
    sample_cell_id = np.concatenate(sample_cell_id)
    out_f_10k["col_attrs/CellID"] = sample_cell_id
    out_f_1k["col_attrs/CellID"] = sample_cell_id
    out_f_10k.close()
    out_f_1k.close()

