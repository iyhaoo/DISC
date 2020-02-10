import argparse
from .utils.data_preparation import *
from .utils.generator import *
from .utils.utility import *
from .model import *
from .evaluation.evaluation import *


def inference(dataset, model, sess, output_dir, batch_size, workers, manager, log_fn=print):
    with h5py.File("{}/running_info.hdf5".format(output_dir), "w") as rf:
        rf["norm_max"] = dataset.norm_max
        rf["norm_mean"] = dataset.z_norm_mean
        rf["norm_std"] = dataset.z_norm_std
        rf["target_gene"] = np.array(dataset.target_gene, dtype=np.string_)
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
    feature_file["col_attrs/CellID"] = dataset.cell_id.astype(np.string_)
    feature_file["row_attrs/Gene"] = np.array(["feature_{}".format(x) for x in range(model.compress_dimensions)], np.string_)
    feature_matrix = feature_file.create_dataset("matrix",
                                                 shape=(model.compress_dimensions, dataset.cell_number),
                                                 chunks=(model.compress_dimensions, 1),
                                                 fletcher32=False,
                                                 dtype=np.float32)
    imputation_file = h5py.File("{}/imputation.loom".format(output_dir), "w")
    imputation_file.create_group("row_graphs")
    imputation_file.create_group("col_graphs")
    imputation_file.create_group("layers")
    imputation_file["col_attrs/CellID"] = dataset.cell_id.astype(np.string_)
    imputation_file["row_attrs/Gene"] = dataset.gene_name.astype(np.string_)
    imputation_matrix = imputation_file.create_dataset("matrix",
                                                       shape=(dataset.gene_number, dataset.cell_number),
                                                       chunks=(dataset.gene_number, 1),
                                                       fletcher32=False,
                                                       dtype=np.float32)
    generator = DataQueue(dataset.loom_path, dataset.gene_name, False, batch_size=batch_size, log_fn=log_fn, workers=workers, manager=manager)
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
    generator.terminate()


def main():
    last_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, type=str, help="loom")
    parser.add_argument("--out-dir", required=True, type=str, help="output folder")
    parser.add_argument("-mc", "--min-expressed-cell", required=False, type=int, default=10, help="min-expressed-cell")
    parser.add_argument("-me", "--min-expressed-cell-average-expression", required=False, type=float, default=1, help="min-expressed-cell-average-expression")
    parser.add_argument("-t", "--repeats", required=False, type=int, default=3, help="repeats")
    parser.add_argument("-sf", "--library-size-factor", required=False, type=str, default="1500", help="int or median")
    parser.add_argument("-d", "--depth", required=False, type=str, default="16_8_1", help="depth")
    parser.add_argument("-s", "--dimension-number", required=False, type=int, default=512, help="Dimension number")
    parser.add_argument("--memory-usage-rate", required=False, type=float, help="How many GPU memory to use (percentage)")
    parser.add_argument("-b", "--batch-size", required=False, type=int, default=128, help="Batch size")
    parser.add_argument("-w", "--compress-dimensions", required=False, type=int, default=50, help="Latent dimensions")
    parser.add_argument("-l", "--learning-rate", required=False, type=float, default=0.001, help="learning-rate")
    parser.add_argument("-tr", "--training", required=False, type=int, default=1, help="is training")
    parser.add_argument("--pretrained-model", required=False, type=str, help="pretrained model path (.pb)")
    parser.add_argument("--scan-workers", required=False, type=int, default=7, help="thread number")
    parser.add_argument("--generator-workers", required=False, type=int, default=3, help="thread number")
    parser.add_argument("--warm-up-cells", required=False, type=int, default=5000000, help="warm-up-cells")
    parser.add_argument("--round-number", required=False, type=int, default=5, help="round number to stop training")
    parser.add_argument("-trs", "--training-round-size", required=False, type=int, default=50000, help="training round size")
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
    makeLog("Dataset: {}".format(FLAGS["dataset"]))
    makeLog("Output Dir: {}".format(FLAGS["out_dir"]))
    makeLog("Batch size: {}".format(FLAGS["batch_size"]))
    ScanLoom_kwargs = {}
    if FLAGS["pretrained_model"] is not None:
        pretrained_gene_name, pretrained_norm_max = get_model_values_by_name(FLAGS["pretrained_model"], ["gene_name", "norm_max"])
        ScanLoom_kwargs["gene_range"] = pretrained_gene_name
    else:
        pretrained_gene_name = None
        pretrained_norm_max = None
    if FLAGS["training"]:
        ScanLoom_kwargs["min_cell"] = FLAGS["min_expressed_cell"]
        ScanLoom_kwargs["min_avg_exp"] = FLAGS["min_expressed_cell_average_expression"]
    else:
        ScanLoom_kwargs["min_cell"] = -1
        ScanLoom_kwargs["min_avg_exp"] = -1
    dataset = ScanLoom(loom_path=FLAGS["dataset"],
                       library_size_factor=FLAGS["library_size_factor"],
                       noise_intensity=0.1,
                       z_score_library_size_factor=1000000,
                       workers=FLAGS["scan_workers"],
                       log_fn=makeLog, **ScanLoom_kwargs)
    #  dataset
    makeLog("Use {} cells".format(dataset.cell_number))
    makeLog("Use {} genes with min_expressed_cell of {} and min_expressed_cell_average_expression of {}".format(dataset.target_gene.size, FLAGS["min_expressed_cell"], FLAGS["min_expressed_cell_average_expression"]))
    makeLog("Use {} as library_size_factor".format(dataset.library_size_factor))
    #  model
    input_norm_max = dataset.norm_max if pretrained_norm_max is None else pretrained_norm_max[np.isin(pretrained_gene_name, dataset.target_gene)]
    model = DISC(gene_name=dataset.target_gene,
                 norm_max=input_norm_max,
                 depth=FLAGS["depth"],
                 repeats=FLAGS["repeats"],
                 dimension_number=FLAGS["dimension_number"],
                 compress_dimensions=FLAGS["compress_dimensions"],
                 noise_intensity=dataset.noise_intensity,
                 z_score_library_size_factor=dataset.z_score_library_size_factor,
                 log_fn=makeLog)
    makeLog("Use {} as depth".format(FLAGS["depth"]))
    makeLog("Repeats {}\n".format(FLAGS["repeats"]))
    feed_dict = {model.z_norm_mean: dataset.z_norm_mean,
                 model.z_norm_std: dataset.z_norm_std,
                 model.zscore_cutoff: dataset.zscore_cutoff,
                 model.library_size_factor: dataset.library_size_factor}
    config = tf.ConfigProto()
    if FLAGS["memory_usage_rate"] is not None:
        config.gpu_options.per_process_gpu_memory_fraction = FLAGS["memory_usage_rate"]
    with tf.Session(config=config) as sess:
        if FLAGS["training"]:
            #  train information
            #  make generator evaluator and training part of model
            train_generator = DataQueue(dataset.loom_path, dataset.target_gene, True, batch_size=FLAGS["batch_size"], log_fn=makeLog, workers=FLAGS["generator_workers"], manager=manager, debug=False)
            evaluator = Evaluation(out_dir=FLAGS["out_dir"],
                                   batch_size=FLAGS["batch_size"],
                                   log_fn=makeLog,
                                   warm_up_cells=FLAGS["warm_up_cells"],
                                   detect_cells=FLAGS["round_number"] * FLAGS["training_round_size"],
                                   manager=manager)
            model.training(FLAGS["learning_rate"])
            feed_dict[model.is_training] = True
            sess.run(tf.global_variables_initializer())
            #  read pre-trained model parameters if a pre-trained model is provided
            if FLAGS["pretrained_model"] is not None:
                sess.run(read_model(FLAGS["pretrained_model"], target_gene=dataset.target_gene, log_fn=makeLog))
            #  initiate
            runnable = True
            next_cell_cutoff = FLAGS["training_round_size"]
            run_cells = 0
            #  training
            while runnable:
                evaluator.evaluation_list.append({"next_cell_cutoff": next_cell_cutoff})
                while run_cells < next_cell_cutoff:
                    _, feed_dict[model.input_raw], feed_dict[model.batch_library_size] = next(train_generator)
                    _, f_loss, run_cells = sess.run([model.train_op1, model.latent_representation_loss, model.run_cells], feed_dict=feed_dict)
                    _, current_batch_size = sess.run([model.train_op2, model.current_batch_size], feed_dict=feed_dict)
                    #  backend evaluation
                    evaluator.evaluation_list.append({"batch_size": current_batch_size, "feature_loss": f_loss})
                #  save model
                model_save_path = "{}/model_{}_cells.pb".format(model_dir, run_cells)
                save_model(sess, model_save_path, [model.gene_name] + model.output_element + model.loss_element)
                makeLog("Run cells: {}\tOptimal Point: {}".format(run_cells, evaluator.optimal_point.value))
                makeLog("Use {:.2f} seconds\n".format(time.time() - last_time))
                last_time = time.time()
                next_cell_cutoff += FLAGS["training_round_size"]
                if evaluator.is_converge.value:
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
        inference(dataset, model, sess, result_dir, FLAGS["batch_size"], FLAGS["generator_workers"], manager, makeLog)
    makeLog("View results at:\n{}".format(result_dir))


if __name__ == "__main__":
    main()

