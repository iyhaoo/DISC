import tensorflow as tf
import numpy as np

#  33
class DeSCI:
    def __init__(self, gene_name, depth="16_8_1", repeats=3, dimension_number=512, compress_dimensions=50,
                 noise_intensity=0.1, dropout_rate=0.5, output_scale_factor=2, z_score_library_size_factor=1000000,
                 en_de_act_fn=tf.nn.tanh, p_act_fn=tf.nn.sigmoid, output_act_fn=tf.nn.sigmoid, log_fn=print):
        self.log_fn = log_fn
        gene_name_array = np.array(gene_name)
        self.use_gene_number = gene_name_array.size
        self.depth = [int(x) for x in depth.split("_")]
        self.repeats = repeats
        self.dimension_number = dimension_number
        self.compress_dimensions = compress_dimensions
        self.noise_intensity = noise_intensity
        self.dropout_rate = dropout_rate
        self.output_scale_factor = output_scale_factor
        self._z_score_library_size_factor = z_score_library_size_factor
        self.en_de_act_fn = en_de_act_fn
        self.p_act_fn = p_act_fn
        self.output_act_fn = output_act_fn
        #  constant tensors
        self.gene_name = tf.constant(gene_name_array, name="gene_name")
        #  feed tensors
        self.input_raw = tf.compat.v1.placeholder(dtype=tf.float32, shape=[None, self.use_gene_number], name="input_layer")
        self.batch_library_size = tf.compat.v1.placeholder(dtype=tf.float32, shape=[None], name="library_size")
        self.norm_max = tf.compat.v1.placeholder(tf.float32, shape=[self.use_gene_number], name="norm_max")
        self.z_norm_mean = tf.compat.v1.placeholder(tf.float32, shape=[self.use_gene_number], name="z_norm_mean")
        self.z_norm_std = tf.compat.v1.placeholder(tf.float32, shape=[self.use_gene_number], name="z_norm_std")
        self.zscore_cutoff = tf.compat.v1.placeholder(tf.float32, shape=[self.use_gene_number], name="zscore_cutoff")
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
        with tf.compat.v1.variable_scope("merge_imputation", reuse=False):
            self._merge_imputation()
            self.merge_gene_features = tf.add_n([feature_ * tf.expand_dims(tf.transpose(atten_), 2) for feature_, atten_ in zip(self.gene_features, tf.split(self.attention_coefficients, self.repeats, 0))])
        with tf.compat.v1.variable_scope("encoder_decoder", reuse=False):
            self._encoder_decoder()
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
        self.encoder_weights = []
        self.impute_latents = []
        self.predictions = []
        self.hidden_layer = []
        self.gene_features = []
        self.prediction_weights = []
        self.prediction_bias = []
        e_W = tf.compat.v1.get_variable("weights", [self.use_gene_number, self.dimension_number], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        self.encoder_weights.append(e_W)
        self.output_act_fn_scale_factor = tf.compat.v1.get_variable("output_act_fn_scale_factor", [1, self.use_gene_number], dtype=tf.float32, initializer=tf.zeros_initializer())
        p_W = tf.compat.v1.get_variable("hidden_weights", [self.dimension_number, self.depth[0] * self.use_gene_number], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        p_B = tf.compat.v1.get_variable("hidden_bias", [self.use_gene_number * self.depth[0]], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.prediction_weights.append(p_W)
        self.prediction_bias.append(p_B)
        for i in range(len(self.depth) - 1):  # [0, 1]
            p_W = tf.compat.v1.get_variable("predictor_weights_{}".format(i), [self.use_gene_number, self.depth[i], self.depth[i + 1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
            p_B = tf.compat.v1.get_variable("predictor_bias_{}".format(i), [self.use_gene_number, 1, self.depth[i + 1]], dtype=tf.float32, initializer=tf.zeros_initializer())
            self.prediction_weights.append(p_W)
            self.prediction_bias.append(p_B)
        self.prediction_weights.append(tf.compat.v1.get_variable("predictor_weights_{}".format(len(self.depth) - 1), [self.use_gene_number, self.depth[-2], self.depth[-1]], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25)))
        self.prediction_bias.append(tf.compat.v1.get_variable("predictor_bias_{}".format(len(self.depth) - 1), [self.use_gene_number, 1, self.depth[-1]], dtype=tf.float32, initializer=tf.zeros_initializer()))
        self.predictor_attention_weights = tf.compat.v1.get_variable("attention_weights", [self.use_gene_number, self.depth[-2], 1], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
        for _ in range(self.repeats):
            features = [tf.matmul(current_input_, e_W) for current_input_ in current_input]
            pre_features = [tf.stop_gradient(feature_) for feature_ in features]
            feature_latent = [tf.expand_dims(tf.transpose(self.output_act_fn_scale_factor + 1), 2) * tf.transpose(tf.reshape(tf.matmul(self.en_de_act_fn(pre_feature_), self.prediction_weights[0]) + self.prediction_bias[0], [tf.shape(pre_feature_)[0], self.use_gene_number, self.depth[0]]), [1, 0, 2]) for pre_feature_ in pre_features]
            self.gene_features.append(feature_latent[0])
            ###
            for i in range(len(self.depth) - 1):  # [0, 1]
                if i == len(self.depth) - 2:
                    output_1 = [tf.expand_dims(tf.transpose(self.output_act_fn_scale_factor + 1), 2) * (tf.matmul(self.en_de_act_fn(latent_), self.prediction_weights[i + 1]) + self.prediction_bias[i + 1]) for latent_ in feature_latent]
                    output_2 = [tf.expand_dims(tf.transpose(self.output_act_fn_scale_factor + 1), 2) * (tf.matmul(self.en_de_act_fn(latent_), self.prediction_weights[i + 2]) + self.prediction_bias[i + 2]) for latent_ in feature_latent]
                    attention_coefficients = tf.sigmoid(tf.nn.selu(tf.matmul(tf.stop_gradient(feature_latent[0]), self.predictor_attention_weights)))
                    feature_latent = [tf.transpose(tf.squeeze(self.output_scale_factor * latent_1 * attention_coefficients, 2) + tf.squeeze(self.output_scale_factor * latent_2 * (1 - attention_coefficients), 2)) for latent_1, latent_2 in zip(output_1, output_2)]
                else:
                    feature_latent = [tf.expand_dims(tf.transpose(self.output_act_fn_scale_factor + 1), 2) * (tf.matmul(self.en_de_act_fn(latent_), self.prediction_weights[i + 1]) + self.prediction_bias[i + 1]) for latent_ in feature_latent]
            self.impute_latents.append(feature_latent)
            self.predictions.append([self.p_act_fn(latent_) for latent_ in feature_latent])
            self.hidden_layer.append(tf.matmul(self.predictions[-1][0] * tf.cast(tf.logical_not(self.known_expressed), tf.float32), e_W))
            current_input = [tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(self.predictions[-1][0])),
                             tf.where(self.known_expressed, self.input_norm_noise_input, tf.nn.dropout(tf.stop_gradient(self.predictions[-1][0]), rate=1 - self.dropout_rate))]

    def _attention(self):
        gene_features = [tf.stop_gradient(tf.nn.selu(feature_)) for feature_ in self.gene_features]
        self.attention_weights = []
        a_W = tf.compat.v1.get_variable("attention_weights", [self.use_gene_number, self.depth[0], 1], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.25, 0.25))
        self.attention_weights.append(a_W)
        self.attention_coefficients = tf.reshape(tf.transpose(tf.nn.softmax(tf.concat([tf.transpose(tf.matmul(feature_, a_W),[1, 0, 2]) for feature_ in gene_features], 2), 2), [2, 0, 1]), [tf.shape(gene_features[0])[1] * len(gene_features), self.use_gene_number])

    def _merge_imputation(self):
        split_attention = tf.split(self.attention_coefficients, self.repeats, 0)
        weighted_latents = [[self.p_act_fn(l_) * a_ for l_ in latent_] for latent_, a_ in zip(self.impute_latents, split_attention)]
        self.merge_impute = [tf.add_n([latent_[i] for latent_ in weighted_latents]) for i in range(len(weighted_latents[0]))]
        self.merge_impute_denorm = [self._denormalization(impute_) for impute_ in self.merge_impute]

    def _encoder_decoder(self):
        predictions = [tf.stop_gradient(prediction_) for prediction_ in [self.input_norm] + [predict_[0] for predict_ in self.predictions]]
        split_attention = tf.split(self.attention_coefficients, self.repeats, 0)
        inputs = [tf.concat([tf.where(self.known_expressed, self.input_norm, prediction_) for prediction_ in predictions], axis=0),
                  tf.concat([tf.where(self.known_expressed, self.input_norm, dropout_prediction_) for dropout_prediction_ in [tf.nn.dropout(prediction_, rate=1 - self.dropout_rate) for prediction_ in predictions]], axis=0)]
        self.hidden_logit = [tf.matmul(input_, self.encoder_weights[0]) for input_ in inputs]
        self.hidden_feature = [self.en_de_act_fn(input_) for input_ in self.hidden_logit]
        self.decoder_bias = tf.compat.v1.get_variable("bias", [self.use_gene_number], dtype=tf.float32, initializer=tf.zeros_initializer())
        all_reconst_prediction = [tf.split(self.output_act_fn(self.output_scale_factor * (self.output_act_fn_scale_factor + 1) * (tf.matmul(feature_, tf.transpose(self.encoder_weights[0])) + self.decoder_bias)), num_or_size_splits=self.repeats + 1, axis=0) for feature_ in self.hidden_feature]
        self.split_reconst_prediction = [reconst_predict_[:-1] for reconst_predict_ in all_reconst_prediction]
        self.reconst_prediction = [tf.add_n([predict_ * atten_ for predict_, atten_ in zip(reconst_predict_, split_attention)]) for reconst_predict_ in self.split_reconst_prediction]

    def _compression(self):
        self.hidden_feature_compression = [tf.concat(tf.split(tf.stop_gradient(feature_), self.repeats + 1, 0), 1) for feature_ in self.hidden_feature]
        self.compress_weight = []
        c_W = tf.compat.v1.get_variable("weights", [self.dimension_number * (self.repeats + 1), self.compress_dimensions], dtype=tf.float32, initializer=tf.random_uniform_initializer(-0.025, 0.025))
        c_B = tf.compat.v1.get_variable("bias_f", [self.compress_dimensions], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.compress_weight.append(c_W)
        self.compressed_feature = [self.en_de_act_fn(tf.matmul(feature_, c_W) + c_B) for feature_ in self.hidden_feature_compression]
        rc_B = tf.compat.v1.get_variable("bias_r", [self.dimension_number * (self.repeats + 1)], dtype=tf.float32, initializer=tf.zeros_initializer())
        self.reconst_feature_compression = [self.en_de_act_fn(tf.matmul(feature_, tf.transpose(c_W)) + rc_B) for feature_ in self.compressed_feature]
        d_W = tf.stop_gradient(tf.transpose(self.encoder_weights[0]))
        d_B = tf.stop_gradient(self.decoder_bias)
        self.compressed_prediction = [tf.add_n(tf.split(self.output_act_fn(self.output_scale_factor * (tf.stop_gradient(self.output_act_fn_scale_factor) + 1) * (tf.matmul(tf.concat(tf.split(feature_, self.repeats + 1, 1)[:-1], 0), d_W) + d_B)) * tf.stop_gradient(self.attention_coefficients), num_or_size_splits=self.repeats, axis=0)) for feature_ in self.reconst_feature_compression]

    def training(self, learning_rate, push_factor=None, gene_express_rate=None):
        """
        learning_rate:
        use_push_factor: "auto", integer, float or None, optional.
                         0 or None means not use.
                         auto means model will tune it automatically.
                         integer or float means fix push_factor.
        """
        #  losses
        #  hyperparameters
        self.global_step = tf.Variable(initial_value=0, expected_shape=(), dtype=tf.int32, name="global_step", trainable=False)
        tf.compat.v1.add_to_collection(tf.compat.v1.GraphKeys.GLOBAL_STEP, self.global_step)
        self.run_cells = tf.Variable(initial_value=0, expected_shape=(), dtype=tf.int32, name="run_cells", trainable=False)
        #  for training:
        if push_factor:
            pf_initial_value = 1. if push_factor == "auto" else float(push_factor)
            self.push_factor = tf.Variable(initial_value=pf_initial_value, expected_shape=(), dtype=tf.float32, name="push_factor", trainable=False)
            expression_mask_1 = [tf.where(tf.logical_not(self.known_expressed), tf.maximum(3., self.push_factor * tf.ones_like(self.input_norm) * tf.cast(self.current_batch_size, tf.float32) * gene_express_rate), tf.zeros_like(self.input_norm)) for _ in self.predictions]
        else:
            expression_mask_1 = [tf.where(tf.logical_not(self.known_expressed), 3. * tf.ones_like(self.input_norm), tf.zeros_like(self.input_norm)) for _ in self.predictions]
        prediction_mask = tf.where(self.known_expressed, 1.5 * tf.ones_like(self.input_norm), 0.35 * tf.ones_like(self.input_norm))
        prediction_target = [tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(predict_)) for predict_ in self.split_reconst_prediction[0]]
        prediction_loss = tf.reduce_sum(tf.reduce_mean(tf.add_n([tf.square(impute_[0] - target_) * prediction_mask for impute_, target_ in zip(self.predictions, prediction_target)]), 0))
        reconst_target = tf.where(self.known_expressed, self.input_norm, tf.stop_gradient(self.merge_impute[0]))
        reconst_mask = tf.where(self.known_expressed, 5 * tf.ones_like(self.input_norm), tf.ones_like(self.input_norm))
        #  noise target's output compared with filtered merged impute
        reconst_loss = tf.reduce_sum(tf.reduce_mean(tf.square(self.reconst_prediction[-1] - reconst_target) * reconst_mask, 0))
        self.reconst_loss_report = tf.reduce_sum(tf.reduce_mean(tf.abs(self._denormalization(self.reconst_prediction[0]) - self.input_raw) * tf.cast(self.known_expressed, tf.float32), 0))
        feature_targets = [tf.zeros_like(self.hidden_layer[0])] + self.hidden_layer[:-1]
        self.feature_loss = tf.add_n([tf.reduce_mean(tf.reduce_sum(tf.square(feature_ - tf.stop_gradient(target_)), 1)) for feature_, target_ in zip(self.hidden_layer, feature_targets)]) / self.repeats
        self.feature_loss_report = tf.add_n([tf.reduce_mean(tf.reduce_sum(tf.square(feature_), 1)) for feature_ in self.hidden_layer]) / self.repeats
        merge_impute_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.merge_impute[-1] - self.input_norm_noise_target) * tf.cast(self.known_expressed, tf.float32), 0))
        self.merge_impute_loss_report = tf.reduce_sum(tf.reduce_mean(tf.abs(self.merge_impute_denorm[0] - self.input_raw) * tf.cast(self.known_expressed, tf.float32), 0))
        single_gene_expression_mask = tf.where(tf.logical_and(tf.logical_not(self.known_expressed), tf.less(self.merge_impute_denorm[0], 0.5)), tf.ones_like(self.input_norm), tf.zeros_like(self.input_norm))
        self.single_gene_loss = tf.reduce_sum(self.merge_impute_denorm[0] / self.batch_library_size_expand_dims * single_gene_expression_mask) / tf.reduce_sum(single_gene_expression_mask)
        compression_loss = tf.reduce_sum(tf.reduce_mean(tf.abs(self.compressed_prediction[0] - tf.stop_gradient(self.merge_impute[0])), 0)) + tf.reduce_sum(tf.reduce_mean(tf.square(self.reconst_feature_compression[0] - self.hidden_feature_compression[0]), 0))
        self.feature_l2 = tf.add_n([tf.reduce_sum(tf.square(self._denormalization(impute_[0])) * expression_mask_) for impute_, expression_mask_ in zip(self.predictions, expression_mask_1)])
        #  for testing:
        expression_mask_2 = tf.where(tf.logical_and(tf.logical_not(self.known_expressed), tf.greater(self.reconst_prediction[0], 0)), tf.ones_like(self.input_norm), tf.zeros_like(self.input_norm))
        expression_mask_3 = tf.where(tf.logical_and(tf.logical_not(self.known_expressed), tf.less(self._denormalization(self.merge_impute[0]), 0.5)), tf.ones_like(self.input_norm), tf.zeros_like(self.input_norm))
        self.compression_loss_report = tf.reduce_sum(tf.reduce_mean(tf.abs(self._denormalization(self.compressed_prediction[0]) - self.input_raw) * tf.cast(self.known_expressed, tf.float32), 0))
        self.feature_l2_1_report = tf.reduce_sum(tf.reduce_sum(self._denormalization(self.merge_impute[0]) * expression_mask_3, 1) / tf.where(tf.greater(tf.reduce_sum(expression_mask_3, 1), 0), tf.reduce_sum(expression_mask_3, 1), tf.reduce_sum(expression_mask_3, 1) + 1))
        self.feature_l2_2_report = tf.reduce_sum(self._denormalization(self.reconst_prediction[0]) * expression_mask_2)
        #  regularizer
        regularizer1 = tf.add_n([tf.nn.l2_loss(tf.reduce_sum(tf.square(w), 0)) for w in self.encoder_weights])
        regularizer2 = tf.add_n([tf.nn.l2_loss(tf.reduce_sum(tf.square(w), 0)) for w in [self.prediction_weights[0]]])
        regularizer3 = tf.add_n([tf.nn.l2_loss(w) for w in self.attention_weights + [self.predictor_attention_weights]])
        regularizer4 = tf.add_n([tf.nn.l2_loss(w) for w in self.prediction_weights[1:]])
        regularizer5 = tf.add_n([tf.nn.l2_loss(w) for w in [self.output_act_fn_scale_factor]])
        regularizer6 = tf.add_n([tf.nn.l2_loss(w) for w in self.compress_weight])

        self.regularizer1 = tf.add_n([tf.reduce_sum(tf.abs(w)) for w in self.encoder_weights])
        self.regularizer2 = tf.add_n([tf.reduce_sum(tf.abs(w)) for w in [self.prediction_weights[0]]])
        self.regularizer3 = tf.add_n([tf.reduce_sum(tf.abs(w)) for w in self.attention_weights])
        self.regularizer4 = tf.add_n([tf.reduce_sum(tf.abs(w)) for w in self.prediction_weights[1:]])
        self.regularizer5 = tf.add_n([tf.reduce_sum(tf.abs(w)) for w in [self.output_act_fn_scale_factor]])
        self.regularizer6 = tf.add_n([tf.reduce_sum(tf.abs(w)) for w in self.compress_weight])
        optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate)
        self.loss1 = prediction_loss + reconst_loss + merge_impute_loss + \
                     0.000021 * self.feature_l2 + 0.165 * self.use_gene_number / 10000 * self.feature_loss + \
                     0.000001 * regularizer1 + 0.000001 * regularizer2 + 0.00001 * regularizer3 + 0.000001 * regularizer4 + 0.0001 * regularizer5
        self.loss2 = compression_loss + 0.0001 * regularizer6
        self.loss_element = [self.loss1, self.loss2]
        gradients, v = zip(*optimizer.compute_gradients(self.loss1))
        gradients, _ = tf.clip_by_global_norm(gradients, 5)
        train_op1 = optimizer.apply_gradients(zip(gradients, v), global_step=self.global_step)
        with tf.control_dependencies([tf.group(train_op1, tf.compat.v1.assign_add(self.run_cells, self.current_batch_size))]):
            if push_factor is "auto":
                self.train_op1 = tf.group(tf.cond(tf.greater(self.run_cells, 50000), lambda: tf.group(tf.cond(self.single_gene_loss * 10000 <= 1, lambda: tf.compat.v1.assign_sub(self.push_factor, tf.minimum(self.push_factor, 0.001)), lambda: tf.compat.v1.assign_add(self.push_factor, 0.001))), lambda: tf.no_op()))
            else:
                self.train_op1 = tf.no_op()
        self.train_op2 = optimizer.minimize(self.loss2)


if __name__ == '__main__':
    import h5py
    loom_path = "E:/DeSCI/fn/melanoma/dropseq_filt_ls.loom"
    # loom_path = "/home/yuanhao/data/fn/melanoma/dropseq_filt_ls.loom"
    with h5py.File(loom_path, "r", libver='latest', swmr=True) as f:
        gene_name = f["row_attrs/Gene"][...].astype(np.str)
    model = DeSCI(gene_name)
    model.training(0.001, 3, np.ones_like(gene_name))
    with tf.compat.v1.Session() as sess:
        print(model.gene_name.eval())





