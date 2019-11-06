import sys
import os
import shutil
import tensorflow as tf


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


def read_model(pretrained_model, log_fn=print):
    with tf.gfile.FastGFile(pretrained_model, 'rb') as f:
        graph_def = tf.GraphDef()
        graph_def.ParseFromString(f.read())
    assign_parameter_run_list = []
    scope_model_variables = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES)
    model_variable_name_list = [x.name.split(":")[0] for x in scope_model_variables]
    read_parameters_number = 0
    for ii in graph_def.node:
        if ii.op == "Const" and ii.name in model_variable_name_list:
            this_value = tf.make_ndarray(ii.attr['value'].tensor)
            read_parameters_number += this_value.size
            this_tensor = scope_model_variables[model_variable_name_list.index(ii.name)]
            if this_value.shape == this_tensor.get_shape():
                assign_parameter_run_list.append(tf.assign(this_tensor, this_value))
            else:
                log_fn("{}:\ntensor: {}\tvalue: {}".format(ii.name, this_tensor.get_shape(), this_value.shape))
    log_fn("Read {} parameters".format(read_parameters_number))
    return assign_parameter_run_list


def get_gene_name(pretrained_model, name="gene_name"):
    with tf.gfile.FastGFile(pretrained_model, 'rb') as f:
        graph_def = tf.GraphDef()
        graph_def.ParseFromString(f.read())
    gene_name = None
    for ii in graph_def.node:
        if ii.op == "Const" and ii.name == name:
            gene_name = tf.make_ndarray(ii.attr['value'].tensor).astype(str)
    assert gene_name is not None
    return gene_name


