import tensorflow as tf
import numpy as np
import copy


class DISC:
    def __init__(self, log_fn=print):
        """
            Parameters
            __________


            log_fn : function, optional, default: print
                Logging function used for this class. Can be specified as a custom function.
        """
        self.log_fn = log_fn
        self.hyperparameter_default = {
            "model_structure": {
                "depth": "16_8_1",
                "repeat": 3,
                "dimension_number": 512,
                "compress_dimension": 50,
                "noise_intensity": 0.1,
                "dropout_rate": 0.5,
                "output_scale_factor": 2,
                "z_score_library_size_factor": 1e6,
            },
            "training": {
                "learning_rate": 0.001,
                "alpha_r": 5,
                "alpha_p1": 1.5,
                "alpha_p2": 0.35,
                "alpha_c": 3,
                "beta_1": 1,
                "beta_2": 1,
                "beta_3": 1,
                "beta_4": "1.65 * m * 1E-5",
                "beta_5": "6.3 * 1E-5",
                "beta_6": 1E-6,
                "beta_7": 1E-6,
                "beta_8": 1E-5,
                "beta_9": 1E-4,
                "beta_10": 1E-4
            }
        }
    def modeling(self, gene_name, norm_max, en_de_act_fn=tf.nn.tanh, output_activation_function=tf.nn.sigmoid, **kwargs):
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

            repeat : int, optional, default: 3
                RNN Step number.

            dimension_number : int, optional, default: 512
                The projected low-dimensional length in AE.

            compress_dimension : int, optional, default: 50
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
        """
        gene_name_array = np.array(gene_name)
        self.gene_number = gene_name_array.size
        hyperparameter_dict = self._kwarg_to_hyperparameters(kwargs, self.hyperparameter_default["model_structure"])
        self.depth = hyperparameter_dict["depth"]
        self.repeat = hyperparameter_dict["repeat"]
        self.dimension_number = hyperparameter_dict["dimension_number"]
        self.compress_dimension = hyperparameter_dict["compress_dimension"]
        noise_intensity = hyperparameter_dict["noise_intensity"]
        self.dropout_rate = hyperparameter_dict["dropout_rate"]
        self.output_scale_factor = hyperparameter_dict["output_scale_factor"]
        self._z_score_library_size_factor = hyperparameter_dict["z_score_library_size_factor"]
        self.en_de_act_fn = en_de_act_fn
        self.output_activation_function = output_activation_function
        #  constant tensors
        self.gene_name = tf.constant(gene_name_array, name="gene_name")
        self.norm_max = tf.constant(norm_max, dtype=tf.float32, name="norm_max")
        self.noise_intensity = tf.constant(noise_intensity, dtype=tf.float32, name="noise_intensity")
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
        self.batch_library_size_expand_dims = tf.expand_dims(self.batch_library_size, 1)
        self.input_not_outlier_mask = tf.logical_not(self._outlier_detector(self.input_raw, 3))
        self.input_norm = self._normalization(self.input_raw * tf.cast(self.input_not_outlier_mask, tf.float32))
        self.known_expressed = tf.greater(self.input_norm, 0)
        self.input_norm_noise_input = self.input_norm * tf.random.uniform([self.current_batch_size, 1], 1 - self.noise_intensity, 1 + self.noise_intensity, dtype=tf.float32)
        self.input_norm_noise_target = self.input_norm * tf.random.uniform([self.current_batch_size, 1], 1 - self.noise_intensity, 1 + self.noise_intensity, dtype=tf.float32)
        with tf.compat.v1.variable_scope("expression_predictor", reuse=False):
            self._ae_rnn()
        with tf.compat.v1.variable_scope("attention", reuse=False):
            self._attention()
        with tf.compat.v1.variable_scope("imputer", reuse=False):
            self._imputer()
            self.merge_gene_features = tf.add_n([feature_ * tf.expand_dims(tf.transpose(a_t), 2) for feature_, a_t in zip(self.gene_features, self.attention_coefficient_list)])
        with tf.compat.v1.variable_scope("reconstructor", reuse=False):
            self._reconstructor()
        with tf.compat.v1.variable_scope("compression", reuse=False):
            self._compression()
        self.output = tf.where(tf.cast(self.input_not_outlier_mask, tf.bool), self.merge_impute_denorm_0, self.input_raw, name="output")
        self.output_element = [self.output, self.output_feature]
        self.constraint_2 = tf.reduce_sum(self._denormalization(self.reconstructor_prediction[0]) * tf.cast(tf.logical_and(tf.logical_not(self.known_expressed), tf.greater(self.reconstructor_prediction[0], 0)), tf.float32))
        self.constraint_3 = tf.reduce_sum(self.merge_impute_denorm_0 * tf.cast(tf.logical_and(tf.logical_not(self.known_expressed), tf.greater(self.merge_impute[0], 0)), tf.float32))
        self.merge_impute_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.merge_impute_denorm_0 - self.input_raw) * tf.cast(self.known_expressed, tf.float32), 0))
        self.trainable_variables = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES)

    def _outlier_detector(self, tensor, cutoff):
        zscore = (tf.math.log1p(tf.stop_gradient(tensor) * self._z_score_library_size_factor / self.batch_library_size_expand_dims) - tf.expand_dims(self.z_norm_mean, 0)) / tf.expand_dims(self.z_norm_std, 0)
        return zscore > cutoff

    def _normalization(self, tensor):
        return_tensor = tf.math.log1p(tensor * self.library_size_factor / self.batch_library_size_expand_dims) / self.norm_max
        return tf.where(tf.math.is_finite(return_tensor), return_tensor, tf.zeros_like(return_tensor))

    def _denormalization(self, tensor):
        return tf.math.expm1(tensor * self.norm_max) / self.library_size_factor * self.batch_library_size_expand_dims

    def _ae_rnn(self):
        current_input = [self.input_norm, self.input_norm_noise_input]
        self.impute_latents = []
        self.predictions = []
        self.hidden_layer = []
        self.gene_features = []
        #  AutoEncoder
        self.weights_encoder = tf.compat.v1.get_variable("weights_encoder", [self.gene_number, self.dimension_number], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        self.phi = tf.compat.v1.get_variable("phi", [1, self.gene_number], dtype=tf.float32, initializer=tf.zeros_initializer())
        #  Prediction Matrix
        self.hidden_layer_1st_weights =tf.compat.v1.get_variable("weights_h1", [self.dimension_number, self.depth[0] * self.gene_number], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        self.hidden_layer_1st_bias = tf.compat.v1.get_variable("bias_h1", [self.gene_number * self.depth[0]], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.hidden_layer_2nd_weights_list = []
        self.hidden_layer_2nd_bias_list = []
        for ii in range(len(self.depth) - 2):  # [0]
            self.hidden_layer_2nd_weights_list.append(tf.compat.v1.get_variable("weights_h2_{}".format(2 + ii), [self.gene_number, self.depth[ii], self.depth[ii + 1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)))
            self.hidden_layer_2nd_bias_list.append(tf.compat.v1.get_variable("bias_h2_{}".format(2 + ii), [self.gene_number, 1, self.depth[ii + 1]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        self.output_layer_weights_list = [tf.compat.v1.get_variable("weights_p1", [self.gene_number, self.depth[-2], self.depth[-1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)),
                                          tf.compat.v1.get_variable("weights_p2", [self.gene_number, self.depth[-2], self.depth[-1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))]
        self.output_layer_bias_list = [tf.compat.v1.get_variable("bias_p1", [self.gene_number, 1, self.depth[-1]], dtype=tf.float32, initializer=tf.zeros_initializer()),
                                       tf.compat.v1.get_variable("bias_p2", [self.gene_number, 1, self.depth[-1]], dtype=tf.float32, initializer=tf.zeros_initializer())]
        self.psi = tf.compat.v1.get_variable("psi", [self.gene_number, self.depth[-2], 1], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
        #  RNN
        for _ in range(self.repeat):
            #  feature refers "latent representation" of our paper
            features_stop = [tf.stop_gradient(tf.matmul(x, self.weights_encoder)) for x in current_input]
            hidden_layer_1_feature = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * tf.transpose(tf.reshape(tf.matmul(self.en_de_act_fn(x), self.hidden_layer_1st_weights) + self.hidden_layer_1st_bias, [tf.shape(x)[0], self.gene_number, self.depth[0]]), [1, 0, 2]) for x in features_stop]
            self.gene_features.append(hidden_layer_1_feature[0])
            #  For the middle layers
            hidden_layer_feature = hidden_layer_1_feature
            for ii in range(len(self.depth) - 2):  # [0]
                hidden_layer_feature = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * (tf.matmul(self.en_de_act_fn(x), self.hidden_layer_2nd_weights_list[ii]) + self.hidden_layer_2nd_bias_list[ii]) for x in hidden_layer_feature]
            #  For the output layer
            output_1 = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * (tf.matmul(self.en_de_act_fn(x), self.output_layer_weights_list[0]) + self.output_layer_bias_list[0]) for x in hidden_layer_feature]
            output_2 = [tf.expand_dims(tf.transpose(self.phi + 1), 2) * (tf.matmul(self.en_de_act_fn(x), self.output_layer_weights_list[1]) + self.output_layer_bias_list[1]) for x in hidden_layer_feature]
            psi = tf.sigmoid(tf.nn.selu(tf.matmul(tf.stop_gradient(hidden_layer_feature[0]), self.psi)))
            output_layer_feature = [tf.transpose(tf.squeeze(self.output_scale_factor * x * psi, 2) + tf.squeeze(self.output_scale_factor * y * (1 - psi), 2)) for x, y in zip(output_1, output_2)]
            self.impute_latents.append(output_layer_feature)
            self.predictions.append([self.output_activation_function(x) for x in output_layer_feature])
            self.hidden_layer.append(tf.matmul(self.predictions[-1][0] * tf.cast(tf.logical_not(self.known_expressed), tf.float32), self.weights_encoder))
            #  Filter
            current_input = [tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(self.predictions[-1][0])),
                             tf.where(self.known_expressed, self.input_norm_noise_input, tf.nn.dropout(tf.stop_gradient(self.predictions[-1][0]), rate=1 - self.dropout_rate))]

    def _attention(self):
        gene_features = [tf.stop_gradient(tf.nn.selu(x)) for x in self.gene_features]
        self.weights_attention = tf.compat.v1.get_variable("weights_attention", [self.gene_number, self.depth[0], 1], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
        self.attention_coefficient_merged = tf.reshape(tf.transpose(tf.nn.softmax(tf.concat([tf.transpose(tf.matmul(x, self.weights_attention),[1, 0, 2]) for x in gene_features], 2), 2), [2, 0, 1]), [tf.shape(gene_features[0])[1] * len(gene_features), self.gene_number])
        self.attention_coefficient_list = tf.split(self.attention_coefficient_merged, self.repeat, 0)

    def _imputer(self):
        weighted_latents = [[self.output_activation_function(l_) * a_ for l_ in latent_] for latent_, a_ in zip(self.impute_latents, self.attention_coefficient_list)]
        self.merge_impute = [tf.add_n([latent_[i] for latent_ in weighted_latents]) for i in range(len(weighted_latents[0]))]
        self.merge_impute_denorm_0 = self._denormalization(self.merge_impute[0])

    def _reconstructor(self):
        predictions = [tf.stop_gradient(prediction_) for prediction_ in [self.input_norm] + [predict_[0] for predict_ in self.predictions]]
        inputs = [tf.concat([tf.where(self.known_expressed, self.input_norm, x) for x in predictions], axis=0),
                  tf.concat([tf.where(self.known_expressed, self.input_norm, tf.nn.dropout(x, rate=1 - self.dropout_rate)) for x in predictions], axis=0)]
        self.hidden_feature = [self.en_de_act_fn(tf.matmul(x, self.weights_encoder)) for x in inputs]
        self.bias_decoder = tf.compat.v1.get_variable("bias_decoder", [self.gene_number], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.split_reconstructor_prediction = [tf.split(self.output_activation_function(self.output_scale_factor * (self.phi + 1) * (tf.matmul(x, tf.transpose(self.weights_encoder)) + self.bias_decoder)), num_or_size_splits=self.repeat + 1, axis=0)[:-1] for x in self.hidden_feature]
        self.reconstructor_prediction = [tf.add_n([predict_ * a_t for predict_, a_t in zip(reconstructor_predict_, self.attention_coefficient_list)]) for reconstructor_predict_ in self.split_reconstructor_prediction]

    def _compression(self):
        self.hidden_feature_compression = [tf.concat(tf.split(tf.stop_gradient(x), self.repeat + 1, 0), 1) for x in self.hidden_feature]
        self.weights_compressor = tf.compat.v1.get_variable("weights_compressor", [self.dimension_number * (self.repeat + 1), self.compress_dimension], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        bias_compressor = tf.compat.v1.get_variable("bias_compressor", [self.compress_dimension], dtype=tf.float32, initializer=tf.zeros_initializer())
        compressed_feature = [self.en_de_act_fn(tf.matmul(x, self.weights_compressor) + bias_compressor) for x in self.hidden_feature_compression]
        self.output_feature = compressed_feature[0]
        bias_compressor_reverse = tf.compat.v1.get_variable("bias_compressor_reverse", [self.dimension_number * (self.repeat + 1)], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.compressor_reconstruction = [self.en_de_act_fn(tf.matmul(x, tf.transpose(self.weights_compressor)) + bias_compressor_reverse) for x in compressed_feature]
        self.compressor_prediction = [tf.add_n(tf.split(self.output_activation_function(self.output_scale_factor * (tf.stop_gradient(self.phi) + 1) * (tf.matmul(tf.concat(tf.split(x, self.repeat + 1, 1)[:-1], 0), tf.stop_gradient(tf.transpose(self.weights_encoder))) + tf.stop_gradient(self.bias_decoder))) * tf.stop_gradient(self.attention_coefficient_merged), self.repeat, 0)) for x in self.compressor_reconstruction]

    def training(self, **kwargs):
        """
            Training function for DISC model.

            Parameters
            __________


            learning_rate : float
                Learning rate for DISC training.

            alpha_r : float
                The weight for positive counts in reconstruction loss calculation.

            alpha_p1 : float
                The weight for positive counts in prediction loss calculation.

            alpha_p2 : float
                The weight for zero counts in prediction loss calculation.

            beta_1 : float
                The weight for imputation loss.

            beta_2 : float
                The weight for reconstruction loss.

            beta_3 : float
                The weight for prediction loss.

            beta_4 : float
                The weight for latent representation loss.

            beta_5 : float
                The weight for constraint.

            beta_6 : float
                The weight for the regulization of encoder weights and 1st hidden layer weights.

            beta_7 : float
                The weight for the regulization of 2nd hidden layer weights and output layer weights.

            beta_8 : float
                The weight for the regulization of attention weights.

            beta_9 : float
                The weight for the regulization of phi.

            beta_10 : float
                The weight for the regulization of compressor weights.

            For all hyperparameters, please define as ii = jj, in which ii is its name and jj is your defined value.
            Type str is supported, and M refers to the gene number of your input data.
            Here are default values for these hyperparameters:
            {
                "learning_rate": 0.001,
                "alpha_r": 5,
                "alpha_p1": 1.5,
                "alpha_p2": 0.35,
                "beta_1": 1,
                "beta_2": 1,
                "beta_3": 1,
                "beta_4": "1.65 * M * 1E-5",
                "beta_5": "6.3 * 1E-5",
                "beta_6": 1E-6,
                "beta_7": 1E-6,
                "beta_8": 1E-5,
                "beta_9": 1E-4,
                "beta_10": 1E-4,
            }.
            We tested the default hyperparameter set for many high-throughput single cell datasets with
            different cell numbers (thousands to millions), different platforms, different cell compositions
            and the performance was robust.
        """
        hyperparameter_dict = self._kwarg_to_hyperparameters(kwargs, self.hyperparameter_default["training"])
        self.global_step = tf.Variable(0, dtype=tf.int32, name="global_step", trainable=False)
        tf.compat.v1.add_to_collection(tf.compat.v1.GraphKeys.GLOBAL_STEP, self.global_step)
        self.run_cells = tf.Variable(0, dtype=tf.int32, name="run_cells", trainable=False)
        #  L_I
        imputation_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.merge_impute[-1] - self.input_norm_noise_target) * tf.cast(self.known_expressed, tf.float32), 0))
        #  L_R
        reconstruction_target = tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(self.merge_impute[0]))
        reconstruction_mask = tf.where(self.known_expressed,
                                       hyperparameter_dict["alpha_r"] * tf.ones_like(self.input_norm),
                                       tf.ones_like(self.input_norm))
        #  noise target's output compared with filtered merged impute
        reconstruction_loss = tf.reduce_sum(tf.reduce_mean(tf.square(self.reconstructor_prediction[-1] - reconstruction_target) * reconstruction_mask, 0))
        #  L_P
        prediction_mask = tf.where(self.known_expressed,
                                   hyperparameter_dict["alpha_p1"] * tf.ones_like(self.input_norm),
                                   hyperparameter_dict["alpha_p2"] * tf.ones_like(self.input_norm))
        prediction_target = [tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(x)) for x in self.split_reconstructor_prediction[0]]
        prediction_loss = tf.reduce_sum(tf.reduce_mean(tf.add_n([tf.square(impute_[0] - target_) * prediction_mask for impute_, target_ in zip(self.predictions, prediction_target)]), 0))
        #  L_LR
        latent_representation_targets = [tf.zeros_like(self.hidden_layer[0])] + self.hidden_layer[:-1]
        self.latent_representation_loss = tf.add_n([tf.reduce_mean(tf.reduce_sum(tf.square(feature_ - tf.stop_gradient(target_)), 1)) for feature_, target_ in zip(self.hidden_layer, latent_representation_targets)]) / self.repeat
        #  L_C
        self.constraint = tf.add_n([tf.reduce_sum(tf.square(self._denormalization(x[0])) * tf.cast(tf.logical_not(self.known_expressed), tf.float32)) for x in self.predictions])
        #  L_compressor
        self.compression_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.compressor_prediction[0] - tf.stop_gradient(self.merge_impute[0])), 0)) + tf.reduce_sum(tf.reduce_mean(tf.square(self.compressor_reconstruction[0] - self.hidden_feature_compression[0]), 0))
        #  regularizer
        regularizer1 = tf.add_n([tf.nn.l2_loss(tf.reduce_sum(tf.square(x), 0)) for x in [self.weights_encoder] + [self.hidden_layer_1st_weights]])
        regularizer2 = tf.add_n([tf.nn.l2_loss(x) for x in self.hidden_layer_2nd_weights_list + self.output_layer_weights_list])
        regularizer3 = tf.add_n([tf.nn.l2_loss(x) for x in [self.weights_attention] + [self.psi]])
        regularizer4 = tf.add_n([tf.nn.l2_loss(x) for x in [self.phi]])
        regularizer_compressor = tf.add_n([tf.nn.l2_loss(x) for x in [self.weights_compressor]])
        optimizer = tf.compat.v1.train.AdamOptimizer(hyperparameter_dict["learning_rate"])
        self.loss1 = hyperparameter_dict["beta_1"] * imputation_loss + \
                     hyperparameter_dict["beta_2"] * reconstruction_loss + \
                     hyperparameter_dict["beta_3"] * prediction_loss + \
                     hyperparameter_dict["beta_4"] * self.latent_representation_loss + \
                     hyperparameter_dict["beta_5"] * self.constraint + \
                     hyperparameter_dict["beta_6"] * regularizer1 + \
                     hyperparameter_dict["beta_7"] * regularizer2 + \
                     hyperparameter_dict["beta_8"] * regularizer3 + \
                     hyperparameter_dict["beta_9"] * regularizer4
        self.loss2 = self.compression_loss + \
                     hyperparameter_dict["beta_10"] * regularizer_compressor
        self.loss_element = [self.loss1, self.loss2]
        self.training_mask = [tf.Variable(tf.ones_like(x), trainable=False) for x in self.trainable_variables]
        self.use_training_mask = tf.Variable(False, dtype=tf.bool, trainable=False)
        gradients, v = zip(*optimizer.compute_gradients(self.loss1, var_list=self.trainable_variables))
        gradients, _ = tf.clip_by_global_norm(gradients, 5)
        gradients = [None if x is None else tf.cond(self.use_training_mask, lambda: x * y, lambda: x) for x, y in zip(gradients, self.training_mask)]
        train_op1 = optimizer.apply_gradients(zip(gradients, v), global_step=self.global_step)
        with tf.control_dependencies([tf.group(train_op1, tf.compat.v1.assign_add(self.run_cells, self.current_batch_size))]):
            self.train_op1 = tf.no_op()
        gradients, v = zip(*optimizer.compute_gradients(self.loss2))
        gradients = [None if x is None else tf.cond(self.use_training_mask, lambda: x * y, lambda: x) for x, y in zip(gradients, self.training_mask)]
        self.train_op2 = optimizer.apply_gradients(zip(gradients, v))

    def _kwarg_to_hyperparameters(self, kwargs, hyperparameter_default):
        M = m = self.gene_number
        hyperparameter_dict = copy.deepcopy(hyperparameter_default)
        is_hyperparameter_default_dict = {}
        for ii in hyperparameter_dict.keys():
            is_hyperparameter_default_dict[ii] = True
        hyperparameter_list = [x for x in hyperparameter_dict.keys()]
        for ii, jj in kwargs.items():
            if np.isin(ii, hyperparameter_list):
                is_hyperparameter_default_dict[ii] = hyperparameter_dict[ii] == jj
                hyperparameter_dict[ii] = jj
        for ii, jj in hyperparameter_dict.items():
            has_changed = False
            if ii == "depth":
                hyperparameter_dict[ii] = [int(x) for x in jj.split("_")]
                has_changed = True
            if isinstance(hyperparameter_dict[ii], str):
                hyperparameter_dict[ii] = eval(hyperparameter_dict[ii])
                has_changed = True
            if is_hyperparameter_default_dict[ii]:
                key_to_print = "{} (Default)".format(ii)
            else:
                key_to_print = ii
            if has_changed:
                value_to_print = "{} ({})".format(jj, (hyperparameter_dict[ii]))
            else:
                value_to_print = jj
            self.log_fn("{} : {}".format(key_to_print, value_to_print))

        return hyperparameter_dict

if __name__ == '__main__':
    import h5py
    loom_path = "E:/DeSCI/fn/melanoma/dropseq_filt_ls.loom"
    # loom_path = "/home/yuanhao/data/fn/melanoma/dropseq_filt_ls.loom"
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        gene_name = f["row_attrs/Gene"][...].astype(np.str)
    model = DISC()
    model.modeling(gene_name, np.ones_like(gene_name))
    model.training(learning_rate=0.001)
    with tf.compat.v1.Session() as sess:
        print(model.gene_name.eval())





