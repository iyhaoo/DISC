import tensorflow as tf
import numpy as np


class DISC:
    """
        An accurate and scalable imputation algorithm based on semi-supervised deep learning for single-cell
        transcriptome, which has an integrative structure of an AE and an RNN.

        Parameters
        __________


        gene_name : str
            Gene names for modelling.

        norm_max : float
            Max normalized expression for all genes specified by gene_name.
            Provided as an attribute in ScanLoom Class.

        depth : str, optional, default: "16_8_1"
            Splitting by "_", int means the channel number in each layer while the number of int means the layer number
            of the prediction matrix.

        repeats : int, optional, default: 3
            RNN Step number.

        dimension_number : int, optional, default: 512
            The projected low-dimensional length in AE.

        compress_dimensions : int, optional, default: 50
            Dimension length when compressing the latent representations over all steps.

        noise_intensity : float, optional, default: 0.1
            Uniform noise intensity.

        dropout_rate : float, optional, default: 0.5
            Dropout rate for noise target.

        output_scale_factor : int, optional, default: 2
            Output scale factor.

        z_score_library_size_factor : int, optional, default: 1000000
            Scale factor for normalization in outlier detection.

        en_de_act_fn : function, optional, default: tf.nn.tanh
            Activation function for autoencoder.

        output_activation_function : function, optional, default: tf.nn.sigmoid
            Activation function for prediction matrix and output.

        log_fn : function, optional, default: print
            Logging function used for this class. Can be specified as a custom function.
    """
    def __init__(self, gene_name, norm_max, depth="16_8_1", repeats=3, dimension_number=512, compress_dimensions=50,
                 noise_intensity=0.1, dropout_rate=0.5, output_scale_factor=2, z_score_library_size_factor=1000000,
                 en_de_act_fn=tf.nn.tanh, output_activation_function=tf.nn.sigmoid, log_fn=print):
        self.log_fn = log_fn
        gene_name_array = np.array(gene_name)
        self.gene_number = gene_name_array.size
        self.depth = [int(x) for x in depth.split("_")]
        self.repeats = repeats
        self.dimension_number = dimension_number
        self.compress_dimensions = compress_dimensions
        self.noise_intensity = noise_intensity
        self.dropout_rate = dropout_rate
        self.output_scale_factor = output_scale_factor
        self._z_score_library_size_factor = z_score_library_size_factor
        self.en_de_act_fn = en_de_act_fn
        self.output_activation_function = output_activation_function
        #  constant tensors
        self.gene_name = tf.constant(gene_name_array, name="gene_name")
        self.norm_max = tf.constant(norm_max, dtype=tf.float32, name="norm_max")
        #  feed tensors
        self.input_raw = tf.compat.v1.placeholder(dtype=tf.float32, shape=[None, self.gene_number], name="input_layer")
        self.batch_library_size = tf.compat.v1.placeholder(dtype=tf.float32, shape=[None], name="library_size")
        self.z_norm_mean = tf.compat.v1.placeholder(tf.float32, shape=[self.gene_number], name="z_norm_mean")
        self.z_norm_std = tf.compat.v1.placeholder(tf.float32, shape=[self.gene_number], name="z_norm_std")
        self.zscore_cutoff = tf.compat.v1.placeholder(tf.float32, shape=[self.gene_number], name="zscore_cutoff")
        self.library_size_factor = tf.compat.v1.placeholder(tf.float32, shape=[], name="library_size_factor")
        self.is_training = tf.compat.v1.placeholder(tf.bool, shape=[], name="is_training")
        #  general attribute tensors
        self.current_batch_size = tf.shape(self.input_raw)[0]
        self.zscore_cutoff_rescale = tf.maximum(-2.5, (tf.math.log1p(tf.math.expm1(self.zscore_cutoff * self.z_norm_std + self.z_norm_mean) / 2) - self.z_norm_mean) / self.z_norm_std)
        self.batch_library_size_expand_dims = tf.expand_dims(self.batch_library_size, 1)
        self.input_outlier_mask = self._deoutlier(self.input_raw, 3., self.z_norm_std, only_mask=True, use_right_tail=True)
        self.input_norm = self._normalization(self.input_raw * self.input_outlier_mask)
        self.known_expressed = tf.greater(self.input_norm, 0)
        self.input_norm_noise_input = self.input_norm * tf.random.uniform([self.current_batch_size, 1], 1 - noise_intensity, 1 + noise_intensity, dtype=tf.float32)
        self.input_norm_noise_target = self.input_norm * tf.random.uniform([self.current_batch_size, 1], 1 - noise_intensity, 1 + noise_intensity, dtype=tf.float32)
        with tf.compat.v1.variable_scope("expression_predictor", reuse=False):
            self._expression_predictor()
        with tf.compat.v1.variable_scope("attention", reuse=False):
            self._attention()
        with tf.compat.v1.variable_scope("imputer", reuse=False):
            self._imputer()
            self.merge_gene_features = tf.add_n([feature_ * tf.expand_dims(tf.transpose(a_t), 2) for feature_, a_t in zip(self.gene_features, self.attention_coefficients_list)])
        with tf.compat.v1.variable_scope("reconstructor", reuse=False):
            self._reconstructor()
            self.feature = tf.concat(tf.split(self.hidden_feature[0], self.repeats, 0), 1)
        with tf.compat.v1.variable_scope("compression", reuse=False):
            self._compression()
        self.output = tf.where(tf.cast(self.input_outlier_mask, tf.bool), self.merge_impute_denorm[0], self.input_raw, name="output")
        self.output_feature = self.compressed_feature[0]
        self.output_element = [self.output, self.output_feature]

    def _deoutlier(self, tensor, cutoff, std, only_mask=False, use_right_tail=False, exclude_mask=None):
        this_score = (tf.math.log1p(tf.stop_gradient(tensor) * self._z_score_library_size_factor / self.batch_library_size_expand_dims) - tf.expand_dims(self.z_norm_mean, 0)) / tf.expand_dims(std, 0)
        if use_right_tail:
            outlier_mask = tf.cast(tf.logical_not(this_score > tf.abs(cutoff)), tf.float32)
        else:
            if exclude_mask is not None:
                outlier_mask = tf.cast(tf.logical_not(tf.logical_and(this_score < -tf.abs(cutoff), tf.logical_not(exclude_mask))), tf.float32)
            else:
                outlier_mask = tf.cast(tf.logical_not(this_score < -tf.abs(cutoff)), tf.float32)
        if only_mask:
            return outlier_mask
        else:
            return tensor * outlier_mask

    def _normalization(self, tensor):
        return_tensor = tf.math.log1p(tensor * self.library_size_factor / self.batch_library_size_expand_dims) / self.norm_max
        return tf.where(tf.math.is_finite(return_tensor), return_tensor, tf.zeros_like(return_tensor))

    def _denormalization(self, tensor):
        return tf.math.expm1(tensor * self.norm_max) / self.library_size_factor * self.batch_library_size_expand_dims

    def _expression_predictor(self):
        current_input = [self.input_norm, self.input_norm_noise_input]
        self.impute_latents = []
        self.predictions = []
        self.hidden_layer = []
        self.gene_features = []
        #  learnable parameters
        self.weights_encoder = tf.compat.v1.get_variable("weights_encoder", [self.gene_number, self.dimension_number], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        self.phi = tf.compat.v1.get_variable("phi", [1, self.gene_number], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.predictor_weights_list = []
        self.predictor_bias_list = []
        self.predictor_weights_list.append(tf.compat.v1.get_variable("weights_hidden_layer_1", [self.dimension_number, self.depth[0] * self.gene_number], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025)))
        self.predictor_bias_list.append(tf.compat.v1.get_variable("bias_hidden_layer_1", [self.gene_number * self.depth[0]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        for i in range(len(self.depth) - 2):  # [0]
            self.predictor_weights_list.append(tf.compat.v1.get_variable("weights_hidden_layer_{}".format(2 + i), [self.gene_number, self.depth[i], self.depth[i + 1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)))
            self.predictor_bias_list.append(tf.compat.v1.get_variable("bias_hidden_layer_{}".format(2 + i), [self.gene_number, 1, self.depth[i + 1]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        self.predictor_weights_list.append(tf.compat.v1.get_variable("weights_output_layer_1", [self.gene_number, self.depth[-2], self.depth[-1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)))
        self.predictor_bias_list.append(tf.compat.v1.get_variable("bias_output_layer_1", [self.gene_number, 1, self.depth[-1]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        self.predictor_weights_list.append(tf.compat.v1.get_variable("weights_output_layer_2", [self.gene_number, self.depth[-2], self.depth[-1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)))
        self.predictor_bias_list.append(tf.compat.v1.get_variable("bias_output_layer_2", [self.gene_number, 1, self.depth[-1]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        self.predictor_weights_list.append(tf.compat.v1.get_variable("predictor_weights_{}".format(len(self.depth) - 1), [self.gene_number, self.depth[-2], self.depth[-1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)))
        self.predictor_bias_list.append(tf.compat.v1.get_variable("predictor_bias_{}".format(len(self.depth) - 1), [self.gene_number, 1, self.depth[-1]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        self.weights_psi = tf.compat.v1.get_variable("weights_psi", [self.gene_number, self.depth[-2], 1], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
        #  predictor
        for _ in range(self.repeats):
            #  feature refers "latent representation" of our paper
            features_stop = [tf.stop_gradient(tf.matmul(current_input_, self.weights_encoder)) for current_input_ in current_input]
            hidden_layer_1_feature = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * tf.transpose(tf.reshape(tf.matmul(self.en_de_act_fn(features_stop_), self.predictor_weights_list[0]) + self.predictor_bias_list[0], [tf.shape(features_stop_)[0], self.gene_number, self.depth[0]]), [1, 0, 2]) for features_stop_ in features_stop]
            self.gene_features.append(hidden_layer_1_feature[0])
            #  For the middle layers
            hidden_layer_feature = hidden_layer_1_feature
            for i in range(len(self.depth) - 2):  # [0]
                hidden_layer_feature = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * (tf.matmul(self.en_de_act_fn(latent_), self.predictor_weights_list[i + 1]) + self.predictor_bias_list[i + 1]) for latent_ in hidden_layer_feature]
            #  For the output layer
            output_1 = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * (tf.matmul(self.en_de_act_fn(latent_), self.predictor_weights_list[-2]) + self.predictor_bias_list[-2]) for latent_ in hidden_layer_feature]
            output_2 = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * (tf.matmul(self.en_de_act_fn(latent_), self.predictor_weights_list[-1]) + self.predictor_bias_list[-1]) for latent_ in hidden_layer_feature]
            psi = tf.sigmoid(tf.nn.selu(tf.matmul(tf.stop_gradient(hidden_layer_feature[0]), self.weights_psi)))
            output_layer_feature = [tf.transpose(tf.squeeze(self.output_scale_factor * latent_1 * psi, 2) + tf.squeeze(self.output_scale_factor * latent_2 * (1 - psi), 2)) for latent_1, latent_2 in zip(output_1, output_2)]
            self.impute_latents.append(output_layer_feature)
            self.predictions.append([self.output_activation_function(latent_) for latent_ in output_layer_feature])
            self.hidden_layer.append(tf.matmul(self.predictions[-1][0] * tf.cast(tf.logical_not(self.known_expressed), tf.float32), self.weights_encoder))
            current_input = [tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(self.predictions[-1][0])),
                             tf.where(self.known_expressed, self.input_norm_noise_input, tf.nn.dropout(tf.stop_gradient(self.predictions[-1][0]), rate=1 - self.dropout_rate))]

    def _attention(self):
        gene_features = [tf.stop_gradient(tf.nn.selu(feature_)) for feature_ in self.gene_features]
        self.weights_attention = tf.compat.v1.get_variable("weights_attention", [self.gene_number, self.depth[0], 1], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
        self.attention_coefficients_merged = tf.reshape(tf.transpose(tf.nn.softmax(tf.concat([tf.transpose(tf.matmul(feature_, self.weights_attention),[1, 0, 2]) for feature_ in gene_features], 2), 2), [2, 0, 1]), [tf.shape(gene_features[0])[1] * len(gene_features), self.gene_number])
        self.attention_coefficients_list = tf.split(self.attention_coefficients_merged, self.repeats, 0)

    def _imputer(self):
        weighted_latents = [[self.output_activation_function(l_) * a_ for l_ in latent_] for latent_, a_ in zip(self.impute_latents, self.attention_coefficients_list)]
        self.merge_impute = [tf.add_n([latent_[i] for latent_ in weighted_latents]) for i in range(len(weighted_latents[0]))]
        self.merge_impute_denorm = [self._denormalization(impute_) for impute_ in self.merge_impute]

    def _reconstructor(self):
        predictions = [tf.stop_gradient(prediction_) for prediction_ in [self.input_norm] + [predict_[0] for predict_ in self.predictions]]
        inputs = [tf.concat([tf.where(self.known_expressed, self.input_norm, prediction_) for prediction_ in predictions], axis=0),
                  tf.concat([tf.where(self.known_expressed, self.input_norm, tf.nn.dropout(prediction_, rate=1 - self.dropout_rate)) for prediction_ in predictions], axis=0)]
        self.hidden_feature = [self.en_de_act_fn(tf.matmul(input_, self.weights_encoder)) for input_ in inputs]
        self.bias_decoder = tf.compat.v1.get_variable("bias_decoder", [self.gene_number], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.split_reconst_prediction = [tf.split(self.output_activation_function(self.output_scale_factor * (self.phi + 1) * (tf.matmul(feature_, tf.transpose(self.weights_encoder)) + self.bias_decoder)), num_or_size_splits=self.repeats + 1, axis=0)[:-1] for feature_ in self.hidden_feature]
        self.reconst_prediction = [tf.add_n([predict_ * a_t for predict_, a_t in zip(reconst_predict_, self.attention_coefficients_list)]) for reconst_predict_ in self.split_reconst_prediction]

    def _compression(self):
        self.hidden_feature_compression = [tf.concat(tf.split(tf.stop_gradient(feature_), self.repeats + 1, 0), 1) for feature_ in self.hidden_feature]
        self.weight_compressor = tf.compat.v1.get_variable("weights_compressor", [self.dimension_number * (self.repeats + 1), self.compress_dimensions], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        bias_compressor = tf.compat.v1.get_variable("bias_compressor", [self.compress_dimensions], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.compressed_feature = [self.en_de_act_fn(tf.matmul(feature_, self.weight_compressor) + bias_compressor) for feature_ in self.hidden_feature_compression]
        bias_compressor_reverse = tf.compat.v1.get_variable("bias_compressor_reverse", [self.dimension_number * (self.repeats + 1)], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.reconst_feature_compression = [self.en_de_act_fn(tf.matmul(feature_, tf.transpose(self.weight_compressor)) + bias_compressor_reverse) for feature_ in self.compressed_feature]
        self.compressed_prediction = [tf.add_n(tf.split(self.output_activation_function(self.output_scale_factor * (tf.stop_gradient(self.phi) + 1) * (tf.matmul(tf.concat(tf.split(feature_, self.repeats + 1, 1)[:-1], 0), tf.stop_gradient(tf.transpose(self.weights_encoder))) + tf.stop_gradient(self.bias_decoder))) * tf.stop_gradient(self.attention_coefficients_merged), self.repeats, 0)) for feature_ in self.reconst_feature_compression]

    def training(self, learning_rate, constraint_factor=1, var_list=None):
        """
            Training function for DISC model.

            Parameters
            __________


            learning_rate : float
                Learning rate for DISC training.

            constraint_factor : int, optional, default: 1
                The factor to limits the total capacity of imputation counts.

            var_list : list, optional, default: None
                List for variables to train.
        """
        self.global_step = tf.Variable(initial_value=0, expected_shape=(), dtype=tf.int32, name="global_step", trainable=False)
        tf.compat.v1.add_to_collection(tf.compat.v1.GraphKeys.GLOBAL_STEP, self.global_step)
        self.run_cells = tf.Variable(initial_value=0, expected_shape=(), dtype=tf.int32, name="run_cells", trainable=False)
        expression_mask_1 = [tf.where(self.known_expressed, tf.zeros_like(self.input_norm), 3. * tf.ones_like(self.input_norm)) for _ in self.predictions]
        prediction_mask = tf.where(self.known_expressed, 1.5 * tf.ones_like(self.input_norm), 0.35 * tf.ones_like(self.input_norm))
        prediction_target = [tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(predict_)) for predict_ in self.split_reconst_prediction[0]]
        prediction_loss = tf.reduce_sum(tf.reduce_mean(tf.add_n([tf.square(impute_[0] - target_) * prediction_mask for impute_, target_ in zip(self.predictions, prediction_target)]), 0))
        reconst_target = tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(self.merge_impute[0]))
        reconst_mask = tf.where(self.known_expressed, 5 * tf.ones_like(self.input_norm), tf.ones_like(self.input_norm))
        #  noise target's output compared with filtered merged impute
        reconstruction_loss = tf.reduce_sum(tf.reduce_mean(tf.square(self.reconst_prediction[-1] - reconst_target) * reconst_mask, 0))
        feature_targets = [tf.zeros_like(self.hidden_layer[0])] + self.hidden_layer[:-1]
        self.latent_representation_loss = tf.add_n([tf.reduce_mean(tf.reduce_sum(tf.square(feature_ - tf.stop_gradient(target_)), 1)) for feature_, target_ in zip(self.hidden_layer, feature_targets)]) / self.repeats
        imputation_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.merge_impute[-1] - self.input_norm_noise_target) * tf.cast(self.known_expressed, tf.float32), 0))
        compression_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.compressed_prediction[0] - tf.stop_gradient(self.merge_impute[0])), 0)) + tf.reduce_sum(tf.reduce_mean(tf.square(self.reconst_feature_compression[0] - self.hidden_feature_compression[0]), 0))
        self.constraint = tf.add_n([tf.reduce_sum(tf.square(self._denormalization(impute_[0])) * expression_mask_) for impute_, expression_mask_ in zip(self.predictions, expression_mask_1)])
        #  regularizer
        regularizer1 = tf.add_n([tf.nn.l2_loss(tf.reduce_sum(tf.square(w), 0)) for w in [self.weights_encoder]])
        regularizer2 = tf.add_n([tf.nn.l2_loss(tf.reduce_sum(tf.square(w), 0)) for w in [self.predictor_weights_list[0]]])
        regularizer3 = tf.add_n([tf.nn.l2_loss(w) for w in [self.weights_attention] + [self.weights_psi]])
        regularizer4 = tf.add_n([tf.nn.l2_loss(w) for w in self.predictor_weights_list[1:]])
        regularizer5 = tf.add_n([tf.nn.l2_loss(w) for w in [self.phi]])
        regularizer6 = tf.add_n([tf.nn.l2_loss(w) for w in [self.weight_compressor]])
        optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate)
        self.loss1 = prediction_loss + reconstruction_loss + imputation_loss + \
                     0.000021 * constraint_factor * self.constraint + 0.165 * self.gene_number / 10000 * self.latent_representation_loss + \
                     0.000001 * regularizer1 + 0.000001 * regularizer2 + 0.00001 * regularizer3 + 0.000001 * regularizer4 + 0.0001 * regularizer5
        self.loss2 = compression_loss + 0.0001 * regularizer6
        self.loss_element = [self.loss1, self.loss2]
        gradients, v = zip(*optimizer.compute_gradients(self.loss1, var_list=var_list))
        gradients, _ = tf.clip_by_global_norm(gradients, 5)
        self.train_op1 = optimizer.apply_gradients(zip(gradients, v), global_step=self.global_step)
        self.train_op2 = optimizer.minimize(self.loss2)


if __name__ == '__main__':
    import h5py
    loom_path = "E:/DeSCI/fn/melanoma/dropseq_filt_ls.loom"
    # loom_path = "/home/yuanhao/data/fn/melanoma/dropseq_filt_ls.loom"
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        gene_name = f["row_attrs/Gene"][...].astype(np.str)
    model = DISC(gene_name, np.ones_like(gene_name))
    model.training(0.001, 3)
    with tf.compat.v1.Session() as sess:
        print(model.gene_name.eval())





