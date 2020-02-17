import sys
import os
import shutil
import tensorflow as tf
import numpy as np
import pandas as pd


class MakeLogClass:
    def __init__(self, log_file):
        self.log_file = log_file
        if os.path.exists(self.log_file):
            os.remove(self.log_file)

    def make(self, *args):
        print(*args)
        # Write the message to the file
        with open(self.log_file, "a") as f:
            for arg in args:
                f.write("{}\r\n".format(arg))


def running_script_backup(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    script_path = sys._getframe().f_code.co_filename.replace("\\", "/").rsplit("/", 2)[0]
    target_dir = "{}/{}".format(save_dir, script_path.replace("\\", "/").rsplit("/", 1)[1])
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir, ignore_errors=True)
    shutil.copytree(script_path, target_dir)


def save_model(sess, model_save_path, output_tensor_list):
    output_graph_def = tf.compat.v1.graph_util.convert_variables_to_constants(sess, tf.get_default_graph().as_graph_def(), [x.name.rsplit(":")[0] for x in output_tensor_list])
    with tf.gfile.GFile(model_save_path, "wb") as f:
        f.write(output_graph_def.SerializeToString())


def read_model(pretrained_model, target_gene=None, log_fn=print):
    with tf.gfile.FastGFile(pretrained_model, 'rb') as f:
        graph_def = tf.GraphDef()
        graph_def.ParseFromString(f.read())
    model_gene_name = None
    if target_gene is not None:
        for ii in graph_def.node:
            if ii.op == "Const" and ii.name == "gene_name":
                model_gene_name = tf.make_ndarray(ii.attr['value'].tensor).astype(str)
        assert model_gene_name is not None
        replace_index, = np.ix_(np.isin(target_gene, model_gene_name))
        model_target_index = pd.Series(range(model_gene_name.size), index=model_gene_name).reindex(target_gene).dropna().values.astype(np.int32)
        assert model_target_index.size > 0, "Error: No intersect gene between input dataset and pretrained model."
    else:
        replace_index = None
        model_target_index = None
    assign_parameter_run_list = []
    scope_model_variables = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES)
    model_variable_name_list = [x.name.split(":")[0] for x in scope_model_variables]
    read_parameters_number = 0
    for ii in graph_def.node:
        if ii.op == "Const" and ii.name in model_variable_name_list:
            read_value = tf.make_ndarray(ii.attr['value'].tensor)
            write_tensor = scope_model_variables[model_variable_name_list.index(ii.name)]
            assign_this = True
            if read_value.shape == write_tensor.shape:
                write_value = read_value
            else:
                if target_gene is not None:
                    if ii.name in ["attention/weights_attention",
                                   "reconstructor/bias_decoder",
                                   "expression_predictor/weights_encoder",
                                   "expression_predictor/weights_psi"]:
                        write_value = write_tensor.eval()
                        write_value[replace_index] = np.take(read_value, model_target_index, axis=0)
                    elif ii.name in ["expression_predictor/phi"]:
                        write_value = write_tensor.eval()
                        write_value[:, replace_index] = np.take(read_value, model_target_index, axis=1)
                    elif ii.name in ["expression_predictor/bias_hidden_layer_1"]:
                        write_value = write_tensor.eval()
                        depth0 = int(read_value.shape[0] / model_gene_name.size)
                        merge_read_index = np.take(np.arange(depth0 * model_gene_name.size).reshape(model_gene_name.size, depth0), model_target_index, axis=0).ravel()
                        merge_write_index = np.take(np.arange(depth0 * target_gene.size).reshape(target_gene.size, depth0), replace_index, axis=0).ravel()
                        write_value[merge_write_index] = np.take(read_value, merge_read_index, axis=0)
                    elif ii.name in ["expression_predictor/weights_hidden_layer_1"]:
                        write_value = write_tensor.eval()
                        depth0 = int(read_value.shape[1] / model_gene_name.size)
                        merge_read_index = np.take(np.arange(depth0 * model_gene_name.size).reshape(model_gene_name.size, depth0), model_target_index, axis=0).ravel()
                        merge_write_index = np.take(np.arange(depth0 * target_gene.size).reshape(target_gene.size, depth0), replace_index, axis=0).ravel()
                        write_value[:, merge_write_index] = np.take(read_value, merge_read_index, axis=1)
                    else:
                        if ii.name.rsplit("_", 1)[0] in ["expression_predictor/weights_hidden_layer",
                                                         "expression_predictor/bias_hidden_layer",
                                                         "expression_predictor/weights_output_layer",
                                                         "expression_predictor/bias_output_layer"]:
                            write_value = write_tensor.eval()
                            write_value[replace_index] = np.take(read_value, model_target_index, axis=0)
                        else:
                            assign_this = False
                            write_value = None
                else:
                    assign_this = False
                    write_value = None
            if assign_this:
                    assign_parameter_run_list.append(tf.assign(write_tensor, write_value))
                    read_parameters_number += write_value.size
            else:
                log_fn("{}:\ntensor: {}\tvalue: {}".format(ii.name, write_tensor.shape, read_value.shape))
    log_fn("Read {} parameters".format(read_parameters_number))
    return assign_parameter_run_list


def get_model_values_by_name(pretrained_model, name_list):
    with tf.gfile.FastGFile(pretrained_model, 'rb') as f:
        graph_def = tf.GraphDef()
        graph_def.ParseFromString(f.read())
    output_list = [None] * len(name_list)
    for ii in graph_def.node:
        if ii.op == "Const" and ii.name in name_list:
            output_list[name_list.index(ii.name)] = tf.make_ndarray(ii.attr['value'].tensor).astype(str)
    for ii in output_list:
        assert ii is not None
    return output_list
