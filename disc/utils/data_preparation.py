import h5py
import numpy as np
from collections import Counter
import os
import pandas as pd
from multiprocessing import Pool
import time


def read_loom(loom_path):
    assert os.path.exists(loom_path)
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        gene_name = f["row_attrs/Gene"][...].astype(np.str)
        assert np.max(list(Counter(gene_name).values())) == 1, "{}".format(list(filter(lambda x: x[1] > 1, Counter(gene_name).items())))
        bc_gene_mat = f["matrix"][...].transpose()
        bc_gene_mat[bc_gene_mat == -1] = np.nan
        cell_id = f["col_attrs/CellID"][...].astype(np.str) if "CellID" in f["col_attrs"].keys() else np.arange(bc_gene_mat.shape[0]).astype(np.str)
    return pd.DataFrame(bc_gene_mat, columns=gene_name, index=cell_id)


class ScanLoom:
    def __init__(self, loom_path, library_size_factor, noise_intensity=0.1, target_gene=None, min_cell=10, min_avg_exp=1,
                 z_score_library_size_factor=1000000, workers=1, scanning_batch_size=2048, log_fn=print):
        """
                Note:  If target_gene is provided, min_cell and min_avg_exp will be ignored
                output gene use the order of loom file
        """
        self.log_fn = log_fn
        self.loom_path = loom_path
        self.z_score_library_size_factor = z_score_library_size_factor
        self._scanning_batch_size = scanning_batch_size
        self.workers = workers
        assert os.path.exists(self.loom_path)
        this_time = time.time()
        self.log_fn("Scan Processes {:.0f}".format(self.workers))
        with h5py.File(self.loom_path, "r", libver='latest', swmr=True) as f:
            gene_name = f["row_attrs/Gene"][...].astype(np.str)
            assert np.max(list(Counter(gene_name).values())) == 1, "{}".format(list(filter(lambda x: x[1] > 1, Counter(gene_name).items())))
            dset = f["matrix"]
            this_shape = dset.shape
            self.gene_number = this_shape[0]
            self.cell_number = this_shape[1]
            cell_id = f["col_attrs/CellID"][...].astype(np.str) if "CellID" in f["col_attrs"].keys() else np.arange(self.cell_number).astype(np.str)
        self.task_number = np.ceil(self.cell_number / self._scanning_batch_size).astype(int)
        if target_gene is not None:
            self.calculate_mask = np.isin(gene_name, target_gene)
        else:
            self.calculate_mask = None
        #  z_* is for z-score calculation
        z_norm_sum, expressed_cell, gene_expression = self._memory_economically_scanning(self._dataset_scanning, ["plus"] * 3)
        if target_gene is None:
            self.calculate_mask = np.bitwise_and(expressed_cell >= min_cell, gene_expression > expressed_cell * min_avg_exp)
            z_norm_sum = z_norm_sum[self.calculate_mask]
        expressed_cell_calculate = expressed_cell[self.calculate_mask]
        self.z_norm_mean = np.divide(z_norm_sum, expressed_cell_calculate, out=np.zeros(self.calculate_mask.sum()), where=expressed_cell_calculate != 0)
        pre_z_norm_std, = self._memory_economically_scanning(self._calculate_norm_std, ["plus"])
        self.z_norm_std = np.sqrt(np.divide(pre_z_norm_std, expressed_cell_calculate, out=np.zeros(self.calculate_mask.sum()), where=expressed_cell_calculate != 0))
        # norm_max filtered by z-score with 1m normalization
        pre_inner_norm_max, self.zscore_cutoff, self.library_size = self._memory_economically_scanning(self._calculate_pre_inner_norm_max, ["max", "min", "append"])
        #  long vector reassignment
        self.gene_name = gene_name
        self.cell_id = cell_id
        self.expressed_cell = expressed_cell
        self.gene_expression = gene_expression
        self.noise_intensity = noise_intensity
        if library_size_factor == "median":
            self.library_size_factor = np.median(self.library_size)
        else:
            self.library_size_factor = float(library_size_factor)
        self.norm_max = np.log1p(pre_inner_norm_max * self.library_size_factor)
        #  gene_name based modification
        if target_gene is not None:
            calculate_mapping_series = pd.Series(range(self.calculate_mask.sum()), index=self.gene_name[self.calculate_mask]).reindex(target_gene[np.isin(target_gene, self.gene_name)])
            calculate_use_index = calculate_mapping_series.values.astype(np.int32)
            extra_genes = np.setdiff1d(target_gene, self.gene_name)
            self.gene_name = np.concatenate([self.gene_name, extra_genes])
            self.gene_number = self.gene_name.size
            gene_name_mapping_series = pd.Series(range(self.gene_number), index=self.gene_name).reindex(np.concatenate([target_gene, np.setdiff1d(self.gene_name, target_gene)]))
            use_index = gene_name_mapping_series.values.astype(np.int32)
            self.gene_name = self.gene_name[use_index]
            self.target_gene_mask = np.isin(self.gene_name, target_gene)
            #  re-assign all gene-related attributes
            self.expressed_cell = np.concatenate([self.expressed_cell, np.zeros(extra_genes.size)])[use_index]
            self.gene_expression = np.concatenate([self.gene_expression, np.zeros(extra_genes.size)])[use_index]
            self.z_norm_mean = np.concatenate([self.z_norm_mean[calculate_use_index], np.zeros(extra_genes.size)])
            self.z_norm_std = np.concatenate([self.z_norm_std[calculate_use_index], np.zeros(extra_genes.size)])
            self.norm_max = np.concatenate([self.norm_max[calculate_use_index], np.zeros(extra_genes.size)])
            self.zscore_cutoff = np.concatenate([self.zscore_cutoff[calculate_use_index], np.zeros(extra_genes.size)])
        else:
            self.target_gene_mask = self.calculate_mask
        #  filter
        self.target_gene = self.gene_name[self.target_gene_mask]
        self.gene_express_rate = self.expressed_cell[self.target_gene_mask] / self.cell_number
        self.norm_max = self.norm_max * (1 + noise_intensity)  # for model noise input
        self.log_fn("Scan Time: {:.2f} Seconds".format(time.time() - this_time))

    #  these functions need to return iterables
    def _dataset_scanning(self, ix):
        with h5py.File(self.loom_path, "r", libver='latest', swmr=True) as f:
            dset = f["matrix"]
            this_batch = dset[:, np.arange(ix, np.minimum(ix + self._scanning_batch_size, self.cell_number), dtype=np.int32)].astype(np.float32)
        this_library_size = np.sum(this_batch, 0)
        this_expressed_cell = (this_batch > 0).sum(1)
        this_gene_expression = this_batch.sum(1)
        use_data = this_batch if self.calculate_mask is None else this_batch[self.calculate_mask, :]
        this_z_norm_sum = np.log1p(np.divide(use_data, this_library_size, out=np.zeros_like(use_data), where=this_library_size != 0) * self.z_score_library_size_factor).sum(1)
        return this_z_norm_sum, this_expressed_cell, this_gene_expression

    def _calculate_norm_std(self, ix):
        with h5py.File(self.loom_path, "r", libver='latest', swmr=True) as f:
            dset = f["matrix"]
            this_batch = dset[:, np.arange(ix, np.minimum(ix + self._scanning_batch_size, self.cell_number), dtype=np.int32)].astype(np.float32)
        this_library_size = np.sum(this_batch, 0)
        use_data = this_batch[self.calculate_mask, :]
        this_z_norm = np.log1p(np.divide(use_data, this_library_size, out=np.zeros_like(use_data), where=this_library_size != 0) * self.z_score_library_size_factor)
        return np.sum(np.where(this_z_norm > 0, np.square(this_z_norm - np.expand_dims(self.z_norm_mean, 1)), np.zeros_like(this_z_norm)), 1),
        
    def _calculate_pre_inner_norm_max(self, ix):
        with h5py.File(self.loom_path, "r", libver='latest', swmr=True) as f:
            dset = f["matrix"]
            this_batch = dset[:, np.arange(ix, np.minimum(ix + self._scanning_batch_size, self.cell_number), dtype=np.int32)].astype(np.float32)
        this_library_size = np.sum(this_batch, 0)
        use_data = this_batch[self.calculate_mask, :]
        this_pre_norm = np.divide(use_data, this_library_size, out=np.zeros_like(use_data), where=this_library_size != 0)
        this_z_norm = np.log1p(this_pre_norm * self.z_score_library_size_factor)
        this_zscore = np.divide(this_z_norm - np.expand_dims(self.z_norm_mean, 1), np.expand_dims(self.z_norm_std, 1), out=np.zeros_like(this_z_norm), where=np.expand_dims(self.z_norm_std, 1) != 0)
        pre_inner_norm_max = np.where(this_zscore > 3, np.zeros_like(this_pre_norm), this_pre_norm).max(1)
        min_zscore = np.where(np.logical_and(use_data > 0, this_zscore >= -2.5), this_zscore, np.ones_like(this_zscore) * 10).min(1)
        return pre_inner_norm_max, min_zscore, this_library_size
    #  ####################################################

    def _memory_economically_scanning(self, worker_function, modify_type_list):
        def master_function(result, to_modify, modify_type):
            return_list = []
            for this_result, this_to_modify, this_modify_type in zip(result, to_modify, modify_type):
                if this_modify_type is "plus":
                    this_to_modify += this_result
                elif this_modify_type is "max":
                    this_to_modify = np.maximum(this_to_modify, this_result)
                elif this_modify_type is "min":
                    this_to_modify = np.minimum(this_to_modify, this_result)
                elif this_modify_type is "append":
                    this_to_modify.append(this_result)
                else:
                    raise Exception("Invalid Input")
                return_list.append(this_to_modify)
            return return_list

        to_modify_list = []
        for this_modify_type in modify_type_list:
            assert this_modify_type in ["plus", "max", "min", "append"]
            if this_modify_type in ["plus", "max", "min"]:
                if self.calculate_mask is None:
                    to_modify_list.append(np.zeros(self.gene_number))
                else:
                    to_modify_list.append(np.zeros(self.calculate_mask.sum()))
            elif this_modify_type in ["append"]:
                to_modify_list.append([])
            else:
                raise Exception("Invalid Setting")
        scan_pool = Pool(processes=self.workers)
        result_list = [scan_pool.apply_async(worker_function, args=(x * self._scanning_batch_size,)) for x in range(self.task_number)]
        scan_pool.close()
        while len(list(filter(lambda x: not x.ready(), result_list))):
            ready_index = []
            for this_index, jj in enumerate(result_list):
                if jj.ready():
                    ready_index.append(this_index)
            for this_index in ready_index[::-1]:
                ready_task = result_list.pop(this_index)
                to_modify_list = master_function(ready_task.get(), to_modify_list, modify_type_list)
        for ii, this_modify_type in enumerate(modify_type_list):
            if this_modify_type in ["append"]:
                to_modify_list[ii] = np.concatenate(to_modify_list[ii])
        return to_modify_list


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument("--loom-path", required=False, default="/home/yuanhao/data/redo/pbmc68k/fresh_68k_pbmc_donor_a_unique_rename.loom", type=str, help="loom")
    parser.add_argument("--loom-path", required=False, default="E:/DeSCI/fn/melanoma/dropseq_filt_ls.loom", type=str, help="loom")
    #parser.add_argument("--loom-path", required=False, default="/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom", type=str, help="loom")
    parser.add_argument("--workers", required=False, default=3, type=int, help="loom")
    FLAGS = vars(parser.parse_args())
    loom_path = FLAGS["loom_path"]
    data = ScanLoom(loom_path, "median", noise_intensity=0.1, workers=FLAGS["workers"])
    print("mean", np.mean(data.library_size))
    print("median", np.median(data.library_size))


