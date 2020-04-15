import numpy as np
import pandas as pd
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool, Manager


class Evaluation:
    """
        A parallel evaluation. Calculate metrics in real time and conduct plotting and indicate optimal point.

        Parameters
        __________


        out_dir : str
            Directory to save output files from this class.

        batch_size : int
            Set batch size, indicate complete batches for evaluation.

        batch_window : int, optional, default: 10000
            Batch window size for evaluation, please note that only complete batches will be calculated.

        warm_up_cells : int, optional, default: 5000000
            Warm up cell number before determining the optimal point.

        detect_cells : int, optional, default: 250000
            Continually run cell number to make sure our optimal point is actually the optimal.

        refresh_time : float, optional, default: 0.1
            Time interval when waiting.

        log_fn : function, optional, default: print
            Logging function used for this class. Can be specified as a custom function.

        manager : manager class object, optional, default: None
            The manager class object. It's recommended only 1 manager is used in the whole program.
    """
    def __init__(self, out_dir, batch_size, batch_window=10000, warm_up_cells=5000000, detect_cells=250000, refresh_time=0.1, log_fn=print, manager=None):
        if manager is None:
            manager = Manager()
        self.pdf_file = "{}/summary.pdf".format(out_dir)
        self.summary_file = "{}/summary.tsv".format(out_dir)
        self.loss_summary = pd.DataFrame()
        self.extra_summary = pd.DataFrame()
        self.runnable = manager.Value("i", 1)
        self.batch_window = batch_window
        self.warm_up_cells = warm_up_cells
        self.optimal_point = manager.Value("i", -1)
        self.is_converge = manager.Value("i", 0)
        self.detect_cells = detect_cells
        self.refresh_time = refresh_time
        self.evaluation_list = manager.list()
        self._set_batch_size = batch_size
        self.run_cells = 0
        self._evaluated_batches = 0
        self._tmp_pd = pd.DataFrame()
        self._pool = Pool(processes=1)
        self._pool.apply_async(func=self._worker_run)
        self.log_fn = log_fn

    def __getstate__(self):
        self_dict = self.__dict__.copy()
        #  remove self unpicklable contents
        del self_dict['_pool']
        return self_dict

    def close(self):
        self.runnable.value = 0
        self._pool.close()
        self._pool.join()

    def _worker_run(self):
        try:
            next_cell_cutoff = 0
            interval_preparation = True
            while self.runnable:
                while len(self.evaluation_list) == 0:
                    time.sleep(self.refresh_time)
                    if not self.runnable.value:
                        return
                ii = self.evaluation_list.pop(0)
                if "next_cell_cutoff" in ii:
                    next_cell_cutoff = ii["next_cell_cutoff"]
                    interval_preparation = True
                    continue
                else:
                    batch_size = ii["batch_size"]
                    del ii["batch_size"]
                    self.run_cells += batch_size
                    for key, value in ii.items():
                        if interval_preparation:
                            self.loss_summary.loc[next_cell_cutoff, key] = value * batch_size
                        else:
                            self.loss_summary.loc[next_cell_cutoff, key] += value * batch_size
                    interval_preparation = False
                    if batch_size == self._set_batch_size:
                        self._tmp_pd = pd.concat([self._tmp_pd, pd.DataFrame().from_dict(ii, orient="index").reindex(list(ii.keys())).T], ignore_index=True)
                        self._evaluated_batches += 1
                if self.run_cells >= next_cell_cutoff:
                    self.loss_summary.rename({next_cell_cutoff: self.run_cells}, axis="index", inplace=True)
                    all_cell_numbers = self.loss_summary.index[::-1]
                    if all_cell_numbers.size > 1:
                        this_cell_number = all_cell_numbers[0] - all_cell_numbers[1]
                    else:
                        this_cell_number = all_cell_numbers[0]
                    self.loss_summary.loc[self.run_cells, :] /= this_cell_number
                    self.log_fn("{}".format(self.loss_summary.loc[self.run_cells, :]))
                    if self._evaluated_batches >= self.batch_window:
                        self._tmp_pd = self._tmp_pd.reindex(range(self._tmp_pd.shape[0])[-self.batch_window:])
                        self.extra_summary.loc[self.run_cells, "feature_mean"] = self._tmp_pd["feature_loss"].mean()
                        self.extra_summary.loc[self.run_cells, "feature_std"] = self._tmp_pd["feature_loss"].std()
                    else:
                        self.extra_summary.loc[self.run_cells, "feature_mean"] = np.nan
                        self.extra_summary.loc[self.run_cells, "feature_std"] = np.nan
                    self.summary = pd.concat([self.loss_summary, self.extra_summary], 1)
                    self.summary.to_csv(self.summary_file, sep="\t")
                    self._plot()
                    if self.run_cells > self.warm_up_cells:
                        if self.run_cells - self.optimal_point.value >= self.detect_cells and self.optimal_point.value != -1:
                            self.is_converge.value = 1
        except Exception as e:
            self.log_fn("Evaluation: {}".format(e))

    def _plot(self):
        run_cells_array = self.summary.index.values
        self.use_data = self.summary
        feature_loss = self.use_data["feature_loss"].values.squeeze()
        ax0_lineur = []
        ax0r_lineur = []
        xlim_max = np.clip(self.run_cells * 1.01, 1, None)
        with PdfPages(self.pdf_file) as pdf:
            if self._evaluated_batches >= self.batch_window:
                feature_convergence_estimate = self.use_data["feature_std"].values.squeeze()
                plot_mask = np.isfinite(feature_convergence_estimate)
                select_mask = run_cells_array >= self.warm_up_cells
                optimal_mask = np.logical_and(plot_mask, select_mask)
                min_feature_smooth_mean_index = run_cells_array[np.nanargmin(self.use_data["feature_mean"])]
                #  plot
                fig, ax0 = plt.subplots(nrows=1, figsize=(12, 6))
                ax0r = ax0.twinx()
                ax0r.set_ylabel("Convergence", color="red")  # we already handled the x-label with ax1
                ax0r_lineur += ax0r.plot(run_cells_array[plot_mask], feature_convergence_estimate[plot_mask], label="Feature Convergence", color="red", linestyle="--")
                ax0r_y_max = np.median([x.get_ydata() for x in ax0r_lineur]) * np.exp(1)
                ax0r_lineur += ax0r.plot([min_feature_smooth_mean_index] * 2, [0, ax0r_y_max], color="blue", label="Min Feature Loss")
                if optimal_mask.sum():
                    optimal_point = run_cells_array[optimal_mask][np.argmin(feature_convergence_estimate[optimal_mask])]
                    self.optimal_point.value = optimal_point.copy()
                    ax0r_lineur += ax0r.plot([optimal_point] * 2, [0, ax0r_y_max], color="red", label="Min Feature Convergence")
                ax0r.tick_params(axis="y", labelcolor="red")
                ax0r.set_xlim([0, xlim_max])
                ax0r.set_ylim([0, ax0r_y_max])
            else:
                fig, ax0 = plt.subplots(nrows=1, figsize=(12, 6))
            ax0.set_ylabel("Loss", color="blue")
            ax0.tick_params(axis="y", labelcolor="blue")
            ax0_lineur += ax0.plot(run_cells_array, feature_loss, label="Feature Loss", color="blue", linestyle="-.")
            ax0_y_max = np.median([x.get_ydata() for x in ax0_lineur]) * 2
            ax0.set_xlim([0, xlim_max])
            ax0.set_ylim([0, ax0_y_max])
            ax0.legend(ax0_lineur + ax0r_lineur, [l.get_label() for l in ax0_lineur + ax0r_lineur], loc="upper right")
            pdf.savefig(fig)
            plt.close()


if __name__ == "__main__":
    out_dir = "E:/DeSCI/fn/DeSCI_3/test"
    evaluator = Evaluation(out_dir, 10000, log_fn=print, batch_window=5)
    save_interval = 50000
    batch_cell = 12000
    run_cells = 0
    next_cell_cutoff = save_interval
    evaluator.evaluation_list.append({"next_cell_cutoff": next_cell_cutoff})
    jj = 0
    for ii in range(100):
        while run_cells < next_cell_cutoff:
            evaluator.evaluation_list.append({"batch_size": batch_cell,
                                              "loss": ii,
                                              "feature_loss": 3 + jj / (ii + 1e-12),
                                              "feature_loss_report": 3.5,
                                              "reconst_loss": 23,
                                              "merge_impute_loss": 3,
                                              "compression_loss": 1})
            jj += 1
            run_cells += batch_cell
        next_cell_cutoff += save_interval
        evaluator.evaluation_list.append({"next_cell_cutoff": next_cell_cutoff})
    time.sleep(5)
    print(evaluator.optimal_point.value)
    evaluator.close()
    print(evaluator.runnable.value)
    print("optimal_point: ", evaluator.optimal_point.value)

