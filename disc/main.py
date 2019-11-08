import argparse
from utils.data_preparation import *
from utils.generator import *
from utils.util import *
from model import *
from evaluation.evaluation import *


def inference(dataset, model, sess, output_dir, batch_size, manager, log_fn=print):
    with h5py.File("{}/running_info.hdf5".format(output_dir), "w") as rf:
        rf["norm_max"] = dataset.norm_max
        rf["norm_mean"] = dataset.z_norm_mean
        rf["norm_std"] = dataset.z_norm_std
        rf["target_gene"] = np.array(dataset.target_gene, dtype=h5py.special_dtype(vlen=str))
        rf["library_size_factor"] = dataset.library_size_factor
        rf["library_size"] = dataset.library_size
        rf["outlier_num"] = dataset.outlier_num
        rf["expressed_cell"] = dataset.expressed_cell
        rf["gene_expression"] = dataset.gene_expression
        rf["zscore_cutoff"] = dataset.zscore_cutoff
    feature_file = h5py.File("{}/feature.loom".format(output_dir), "w")
    feature_file.create_group("row_graphs")
    feature_file.create_group("col_graphs")
    feature_file.create_group("layers")
    feature_file["col_attrs/CellID"] = dataset.cell_id.astype(h5py.special_dtype(vlen=str))
    feature_file["row_attrs/Gene"] = np.array(["feature_{}".format(x) for x in range(model.compress_dimensions)], h5py.special_dtype(vlen=str))
    feature_matrix = feature_file.create_dataset("matrix",
                                                 shape=(model.compress_dimensions, dataset.cell_number),
                                                 chunks=(model.compress_dimensions, 1),
                                                 fletcher32=False,
                                                 dtype=np.float32)
    imputation_file = h5py.File("{}/imputation.loom".format(output_dir), "w")
    imputation_file.create_group("row_graphs")
    imputation_file.create_group("col_graphs")
    imputation_file.create_group("layers")
    imputation_file["col_attrs/CellID"] = dataset.cell_id.astype(h5py.special_dtype(vlen=str))
    imputation_file["row_attrs/Gene"] = dataset.gene_name.astype(h5py.special_dtype(vlen=str))
    imputation_matrix = imputation_file.create_dataset("matrix",
                                                       shape=(dataset.gene_number, dataset.cell_number),
                                                       chunks=(dataset.gene_number, 1),
                                                       fletcher32=False,
                                                       dtype=np.float32)
    generator = DataQueue(dataset.loom_path, dataset.gene_name, False, batch_size=batch_size, log_fn=log_fn, workers=FLAGS["generator_workers"], manager=manager)
    feed_dict = {model.norm_max: dataset.norm_max,
                 model.z_norm_mean: dataset.z_norm_mean,
                 model.z_norm_std: dataset.z_norm_std,
                 model.zscore_cutoff: dataset.zscore_cutoff,
                 model.library_size_factor: dataset.library_size_factor,
                 model.is_training: False}
    for _ in range(generator.steps_per_epoch):
        batch_index, batch_data, feed_dict[model.batch_library_size] = next(generator)
        impute_data = np.empty(batch_data.shape)
        feed_dict[model.input_raw] = batch_data[:, dataset.target_gene_mask]
        impute_data[:, dataset.target_gene_mask], output_feature = sess.run(model.output_element, feed_dict=feed_dict)
        feature_matrix[:, batch_index] = output_feature.transpose()
        imputation_matrix[:, batch_index] = np.where(np.broadcast_to(dataset.target_gene_mask, batch_data.shape), impute_data, batch_data).transpose()
    feature_file.close()
    imputation_file.close()


if __name__ == "__main__":
    last_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, type=str, help="loom")
    parser.add_argument("--out-dir", required=True, type=str, help="output folder")
    parser.add_argument("--min-expressed-cell", required=False, type=int, default=10, help="min-expressed-cell")
    parser.add_argument("--min-expressed-cell-average-expression", required=False, type=float, default=1, help="min-expressed-cell-average-expression")
    parser.add_argument("--repeats", required=False, type=int, default=1, help="repeats")
    parser.add_argument("--library-size-factor", required=False, type=str, default="1500", help="int or median")
    parser.add_argument("--depth", required=False, type=str, default="16_8_1", help="depth")
    parser.add_argument("--dimension-number", required=False, type=int, default=512, help="Dimension number")
    parser.add_argument("--memory-usage-rate", required=False, type=float, help="How many memory to use")
    parser.add_argument("--batch-size", required=False, type=int, default=128, help="Batch size")
    parser.add_argument("--compress-dimensions", required=False, type=int, default=50, help="Latent dimensions")
    parser.add_argument("--noise-intensity", required=False, type=float, default=0.1, help="noise-intensity")
    parser.add_argument("--feature-l2-factor", required=False, type=float, default=1, help="feature_l2_factor")
    parser.add_argument("--z-score-library-size-factor", required=False, type=float, default=1000000, help="z-score-library-size-factor")
    parser.add_argument("--learning-rate", required=False, type=float, default=0.001, help="learning-rate")
    parser.add_argument("--training", required=False, type=int, default=1, help="is training")
    parser.add_argument("--pretrained-model", required=False, type=str, help="pretrained model path (.pb)")
    parser.add_argument("--scan-workers", required=False, type=int, default=7, help="thread number")
    parser.add_argument("--generator-workers", required=False, type=int, default=3, help="thread number")
    parser.add_argument("--converge-number", required=False, type=int, default=10000000, help="Max cell number for training")
    parser.add_argument("--warm-up-cells", required=False, type=int, default=5000000, help="warm-up-cells")
    parser.add_argument("--save-interval", required=False, type=int, default=50000, help="save-interval")
    FLAGS = vars(parser.parse_args())
    manager = Manager()
    FLAGS["return_elements_str"] = ""
    result_dir = "{}/result".format(FLAGS["out_dir"])
    os.makedirs(result_dir, exist_ok=True)
    makeLog = MakeLogClass("{}/log.tsv".format(FLAGS["out_dir"])).make
    running_script_backup("{}/run_script".format(FLAGS["out_dir"]))
    if not FLAGS["training"]:
        assert FLAGS["pretrained_model"] is not None
        model_dir = None
    else:
        model_dir = "{}/models".format(FLAGS["out_dir"])
        os.makedirs(model_dir, exist_ok=True)
    makeLog("Batch size: {}".format(FLAGS["batch_size"]))

    if FLAGS["pretrained_model"] is not None:
        kwargs = {"target_gene": get_gene_name(FLAGS["pretrained_model"])}
    else:
        kwargs = {"min_cell": FLAGS["min_expressed_cell"],
                  "min_avg_exp": FLAGS["min_expressed_cell_average_expression"]}
    dataset = ScanLoom(loom_path=FLAGS["dataset"],
                       library_size_factor=FLAGS["library_size_factor"],
                       noise_intensity=FLAGS["noise_intensity"],
                       z_score_library_size_factor=FLAGS["z_score_library_size_factor"],
                       workers=FLAGS["scan_workers"],
                       log_fn=makeLog, **kwargs)
    #  dataset
    makeLog("Use {} cells".format(dataset.cell_number))
    makeLog("Use {} genes with min_expressed_cell of {} and min_expressed_cell_average_expression of {}".format(dataset.target_gene.size, FLAGS["min_expressed_cell"], FLAGS["min_expressed_cell_average_expression"]))
    makeLog("Use {} as library_size_factor".format(dataset.library_size_factor))
    #  model
    model = DeSCI(gene_name=dataset.target_gene,
                  depth=FLAGS["depth"],
                  repeats=FLAGS["repeats"],
                  dimension_number=FLAGS["dimension_number"],
                  compress_dimensions=FLAGS["compress_dimensions"],
                  noise_intensity=dataset.noise_intensity,
                  z_score_library_size_factor=dataset.z_score_library_size_factor,
                  log_fn=makeLog)
    makeLog("Use {} as depth".format(FLAGS["depth"]))
    makeLog("Repeats {}".format(FLAGS["repeats"]))
    feed_dict = {model.norm_max: dataset.norm_max,
                 model.z_norm_mean: dataset.z_norm_mean,
                 model.z_norm_std: dataset.z_norm_std,
                 model.zscore_cutoff: dataset.zscore_cutoff,
                 model.library_size_factor: dataset.library_size_factor}
    config = tf.ConfigProto()
    if FLAGS["memory_usage_rate"] is not None:
        config.gpu_options.per_process_gpu_memory_fraction = FLAGS["memory_usage_rate"]
    with tf.Session(config=config) as sess:
        if FLAGS["training"]:
            #  train information
            makeLog("Max cell number for training: {:.0f}".format(FLAGS["converge_number"]))
            #  make generator evaluator and training part of model
            train_generator = DataQueue(dataset.loom_path, dataset.target_gene, True, batch_size=FLAGS["batch_size"], log_fn=makeLog, workers=FLAGS["generator_workers"], manager=manager)
            evaluator = Evaluation(out_dir=FLAGS["out_dir"], batch_size=FLAGS["batch_size"], log_fn=makeLog, warm_up_cells=FLAGS["warm_up_cells"], manager=manager)
            model.training(FLAGS["learning_rate"], feature_l2_factor=FLAGS["feature_l2_factor"], push_factor=None, gene_express_rate=dataset.gene_express_rate)
            feed_dict[model.is_training] = True
            sess.run(tf.global_variables_initializer())
            #  read pre-trained model parameters if a pre-trained model is provided
            if FLAGS["pretrained_model"] is not None:
                sess.run(read_model(FLAGS["pretrained_model"], log_fn=makeLog))
            #  initiate
            runnable = True
            next_cell_cutoff = FLAGS["save_interval"]
            run_cells = 0
            #  training
            while runnable:
                evaluator.evaluation_list.append({"next_cell_cutoff": next_cell_cutoff})
                while run_cells < next_cell_cutoff:
                    _, feed_dict[model.input_raw], feed_dict[model.batch_library_size] = next(train_generator)
                    _, f_loss, run_cells = sess.run([model.train_op1, model.feature_loss, model.run_cells], feed_dict=feed_dict)
                    _, current_batch_size = sess.run([model.train_op2, model.current_batch_size], feed_dict=feed_dict)
                    #  backend evaluation
                    evaluator.evaluation_list.append({"batch_size": current_batch_size, "feature_loss": f_loss})
                #  save model
                model_save_path = "{}/model_{}_cells.pb".format(model_dir, run_cells)
                save_model(sess, model_save_path, [model.gene_name] + model.output_element + model.loss_element)
                makeLog("Run cells: {}\tOptimal Point: {}".format(run_cells, evaluator.optimal_point.value))
                makeLog("Use {:.2f} seconds\n".format(time.time() - last_time))
                last_time = time.time()
                next_cell_cutoff += FLAGS["save_interval"]
                if run_cells >= FLAGS["converge_number"] or evaluator.is_converge.value:
                    runnable = False
            train_generator.terminate()
            del train_generator
            optimal_point = evaluator.optimal_point.value
            evaluator.close()
            del evaluator
            use_pretrained_model = "{}/model_{}_cells.pb".format(model_dir, optimal_point)
        else:
            sess.run(tf.global_variables_initializer())
            use_pretrained_model = FLAGS["pretrained_model"]
        sess.run(read_model(use_pretrained_model, log_fn=makeLog))
        makeLog("Use {} for inference".format(use_pretrained_model))
        inference(dataset, model, sess, result_dir, FLAGS["batch_size"], manager, makeLog)


