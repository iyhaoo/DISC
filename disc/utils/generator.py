import numpy as np
import h5py
from multiprocessing import Pool, Manager
import os
import time
import pandas as pd


class DataQueue:
    """
        An ultra fast generator for ultra large dataset.

        Parameters
        __________


        loom_path : str
            Input data-set path. Should be a loom-formatted file containing a gene expression matrix with genes in rows
            and cells in columns.

        target_gene : str, optional, default: None
            Only calculate specific genes, is useful in our transfer learning module.
            Output columns will follow the order of specific genes input here.

        permutation : bool, optional, default: True
            Whether use random output.

        batch_size : int, optional, default: 128
            Set batch size for output.

        chunk_number : int, optional, default: 64
            Reading chunk number, larger number for better stochasticity performance but higher memory occupation.

        chunk_size : int, optional, default: 32
            Reading chunk size, smaller number for better stochasticity performance but lower speed.

        workers : int, optional, default: 6
            Process number for reading. Larger number for better performance but higher CPU and memory occupation.

        prefetch : int, optional, default: 12
            max prefetch number for localization queue.

        refill_cutoff : int, optional, default: None
            Cutoff for starting new task. None means use automatic setting.

        log_fn : function, optional, default: print
            Logging function used for this class. Can be specified as a custom function.

        output_type : dtype, optional, default: np.float32
            Set output dtype.

        manager : manager class object, optional, default: None
            The manager class object. It's recommended only 1 manager is used in the whole program.

        debug : bool, optional, default: True
            Whether use the debug mode. Print lagged reading if lagged.
    """
    def __init__(self, loom_path, target_gene, permutation=True, batch_size=128, chunk_number=64, chunk_size=32,
                 workers=6, prefetch=12, refill_cutoff=None, log_fn=print, output_type=np.float32, manager=None,
                 debug=True):
        assert os.path.exists(loom_path)
        assert workers > 0
        if manager is None:
            manager = Manager()
        self.log_fn = log_fn
        self.output_type = output_type
        #  f open
        f = h5py.File(loom_path, "r", libver='latest', swmr=True)
        matrix_shape = f["matrix"].shape
        self.gene_number = matrix_shape[0]
        self.cell_number = matrix_shape[1]
        self.cell_id = f["col_attrs/CellID"][...].astype(np.str) if "CellID" in f["col_attrs"].keys() else np.arange(self.cell_number).astype(np.str)
        self.gene_name = f["row_attrs/Gene"][...].astype(np.str)
        self.is_overlap = self.cell_number < workers * chunk_size * chunk_number
        #  use shortcuts when reading small datasets or testing
        if self.is_overlap:
            workers = 2
            chunk_number = 1
            chunk_size = self.cell_number
            self._cache_data = f["matrix"][...]
            self.log_fn("Use small dataset read mode")
        if not permutation:
            chunk_size = chunk_number * chunk_size
            chunk_number = 1
        f.close()
        #  f close
        """
            define task queue and two queues for pre-processed data
            set up a record for running task number
            make a lock to divide epoch data
        """
        self._queue_task = manager.list()
        self._queue_t = manager.list()
        self._queue_f = manager.list()
        self._queue_local = []
        self._running_t = manager.Value("i", 0)
        self._running_f = manager.Value("i", 0)
        self._lock_t = manager.Value("i", 0)
        self._lock_f = manager.Value("i", 0)
        self.steps_per_epoch = np.ceil(self.cell_number / batch_size).astype(int)
        self._channel_remain = 0
        if prefetch is not None:
            self.prefetch = prefetch
        else:
            self.prefetch = np.ceil(chunk_number * chunk_size * workers / batch_size).astype(int)
        if refill_cutoff is not None:
            self.refill_cutoff = refill_cutoff
        else:
            self.refill_cutoff = np.ceil(chunk_number * chunk_size * (workers + 1) / batch_size).astype(int)
        self.prefetch = min(self.refill_cutoff, self.prefetch)
        self._chunk_index = np.arange(0)
        self._pool = Pool(processes=workers)
        self._pool_process = []
        #   attributes
        #   use channel t as block channel to prevent data disorder
        start_channel = False
        # tftftftf
        self._fill_queue_channel = start_channel
        # tftftftf
        self._get_queue_channel = start_channel
        self._loom_path = loom_path
        self._permutation = permutation
        self.batch_size = batch_size
        self._chunk_number = chunk_number
        self._chunk_size = chunk_size
        mapping_series = pd.Series(range(self.gene_number), index=self.gene_name).reindex(target_gene)
        self._keep_index = mapping_series.dropna().values.astype(np.int32)
        insertion_pre_index, = np.where(np.isnan(mapping_series))
        self._insertion_index = insertion_pre_index - np.arange(insertion_pre_index.size)
        self._debug = debug
        while self._init_worker():
            pass

    def _select_queue(self, channel):
        if channel:
            return self._queue_t, self._running_t, self._lock_t
        else:
            return self._queue_f,  self._running_f, self._lock_f

    def _init_worker(self):
        self._pool_process = list(filter(lambda x: not x.ready(), self._pool_process))
        if not self._chunk_index.size:
            self._fill_queue_channel = np.bitwise_not(self._fill_queue_channel)
            self._chunk_index = np.arange(np.ceil(self.cell_number / self._chunk_size).astype(np.int32))
            if self._permutation:
                self._chunk_index = np.random.permutation(self._chunk_index)
            if self._fill_queue_channel and self._running_t.value > 0:
                self._lock_t.value = 1
            if not self._fill_queue_channel and self._running_f.value > 0:
                self._lock_f.value = 1
        if (self._fill_queue_channel and self._lock_t.value) or (not self._fill_queue_channel and self._lock_f.value):
            return False
        finished_batch = len(self._queue_t) + len(self._queue_f)
        running_batch = len(self._pool_process) * self._chunk_size * self._chunk_number / self.batch_size  # an approximate estimation
        if finished_batch + running_batch > self.refill_cutoff:
            return False
        self._queue_task.append((self._chunk_index[:self._chunk_number], self._fill_queue_channel))
        self._chunk_index = self._chunk_index[self._chunk_number:]
        self._pool_process.append(self._pool.apply_async(func=self._fill_queue))
        return True

    def _fill_queue(self):
        try:
            while len(self._queue_task) == 0:
                pass
            use_chunk_index, this_queue_channel = self._queue_task.pop(0)
            this_queue, this_running, this_lock = self._select_queue(this_queue_channel)
            this_running.value += 1
            #self.log_fn("T", len(self._queue_t), self._running_t.value, self._lock_t.value)
            #self.log_fn("F", len(self._queue_f), self._running_f.value, self._lock_f.value)
            #  use memory cache to reduce disk reading operations
            if self.is_overlap:
                chunk_id = np.arange(self.cell_number)
                chunk_data = self._cache_data.copy().transpose()
            else:
                with h5py.File(self._loom_path, "r", libver='latest', swmr=True) as f:
                    dataset = f["matrix"]
                    chunk_list = []
                    chunk_id_list = []
                    #  read chunk
                    for ii in use_chunk_index:
                        chunk_list.append(dataset[:, ii * self._chunk_size:(ii + 1) * self._chunk_size])
                        chunk_id_list.append(np.arange(ii * self._chunk_size, np.minimum(self.cell_number, (ii + 1) * self._chunk_size), dtype=np.int32))
                #  reduce redundant operations
                if self._permutation:
                    chunk_id = np.concatenate(chunk_id_list)
                    chunk_data = np.hstack(chunk_list).transpose()
                else:
                    chunk_id = chunk_id_list[0]
                    chunk_data = chunk_list[0].transpose()
            #  if end
            if self._permutation:
                shuffle_index = np.random.permutation(chunk_id.size)
                chunk_id = chunk_id[shuffle_index]
                chunk_data = chunk_data[shuffle_index, :]
            chunk_data = chunk_data.astype(self.output_type)
            chunk_library_size = chunk_data.sum(1)
            chunk_data = np.insert(chunk_data[:, self._keep_index], self._insertion_index, 0, axis=1)
            for _ in np.arange(np.ceil(chunk_id.size / self.batch_size), dtype=np.int32):
                this_chunk_id = chunk_id[:self.batch_size]
                chunk_id = chunk_id[self.batch_size:]
                this_chunk_data = chunk_data[:self.batch_size, :]
                chunk_data = chunk_data[self.batch_size:, :]
                this_chunk_library_size = chunk_library_size[:self.batch_size]
                chunk_library_size = chunk_library_size[self.batch_size:]
                while len(self._queue_t) + len(self._queue_f) > self.refill_cutoff:
                    pass
                this_queue.append({"index": this_chunk_id, "data": this_chunk_data, "library_size": this_chunk_library_size})
            this_running.value -= 1
            if not this_running.value:
                this_lock.value = 0
        except Exception as e:
            self.log_fn(e)

    def _lagged_loading(self, shared_queue):
        last_time = time.time()
        while len(shared_queue) == 0:
            pass
        self._queue_local.append(shared_queue.pop(0))
        if self._debug:
            self.log_fn("Lag Load: {:.2f} seconds".format(time.time() - last_time))
        self._channel_remain -= 1

    def _queue_localize(self):
        finished_batch = len(self._queue_t) + len(self._queue_f)
        running_batch = (self._running_t.value + self._running_f.value) * self._chunk_size * self._chunk_number / self.batch_size
        if finished_batch + running_batch < self.refill_cutoff:
            self._init_worker()
        #  when an epoch is used, switch
        if self._channel_remain == 0:
            self._get_queue_channel = np.bitwise_not(self._get_queue_channel)
            self._channel_remain = self.steps_per_epoch.copy()
        this_queue, _, _ = self._select_queue(self._get_queue_channel)
        #  wait for loading
        if len(this_queue) == 0:
            self._lagged_loading(this_queue)
            return
        prefetch_number = int(min(len(this_queue), self._channel_remain, self.prefetch))
        for _ in range(prefetch_number):
            self._queue_local.append(this_queue.pop(0))
        self._channel_remain -= prefetch_number

    def __next__(self):
        if len(self._queue_local) == 0:
            self._queue_localize()
        this_dict = self._queue_local.pop(0)
        return this_dict["index"], this_dict["data"], this_dict["library_size"]

    def __getstate__(self):
        self_dict = self.__dict__.copy()
        #  remove self unpicklable contents
        del self_dict['_pool']
        del self_dict['_pool_process']
        return self_dict

    def terminate(self):
        self._pool.terminate()
        del self

    def close(self):
        self._pool.close()
        self._pool.join()


if __name__ == '__main__':
    loom_path = "E:/DeSCI/fn/melanoma/dropseq_filt_ls.loom"
    #loom_path = "/home/yuanhao/data/fn/melanoma/dropseq_filt_ls.loom"
    #loom_path = "/home/yuanhao/data/drop_seq/dropseq_filt_ls.loom"
    #loom_path = "E:/DeSCI/fn/pbmc3k/ds_to_500/downsampling_first_repeat_4/pbmc3k_filtered_ds_to_500.loom"
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        gene_name = f["row_attrs/Gene"][...].astype(np.str)
    queue = DataQueue(loom_path, gene_name, workers=3, prefetch=10, permutation=False, debug=True)
    print(queue.is_overlap)
    ii = 0
    """
    print("cell: ", queue.cell_number)
    print("epoch: ", queue.steps_per_epoch)
    index_list = []
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        for jjj in range(queue.steps_per_epoch * 10):
            index, data, ls = next(queue)
            index_list.append(index)
            if (jjj + 1) % queue.steps_per_epoch == 0:
                diff = np.setdiff1d(np.arange(queue.cell_number), np.concatenate(index_list))
                print("different", diff.size, ii)
                index_list = []
            order = np.argsort(index)
            direct_get_data = f["matrix"][:, index[order]].transpose()
            dls = direct_get_data.sum(1)
            print(dls.size)
            print(np.all(data[order, :] == direct_get_data))
            print(np.all(ls[order] == dls))
            #next(queue)
            # print(ii, index)
            # print(ii)
            ii += 1
    """
    last_time = time.time()
    while True:
        next(queue)
        # print(ii)
        ii += 1
        if time.time() - last_time > 15:
            queue.terminate()
            del queue
            print(ii)
            break
    while True:
        next(queue)
        # ii += 1
        # print(ii)





