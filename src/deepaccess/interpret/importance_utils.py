import tensorflow as tf
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input
from tensorflow.keras import activations
import numpy as np


# saliency code modified from
# https://github.com/keras-team/keras/blob/master/examples/conv_filter_visualization.py
# All contributions by François Chollet:
# Copyright (c) 2015 - 2019, François Chollet.
# All rights reserved.

# All other contributions:
# Copyright (c) 2015 - 2019, the respective contributors.
# All rights reserved.
def saliency(
    conv_layer,
    model_file,
    input_data,
    pos_index,
    neg_index,
    n=5,
    lambda_reg=0.1,
    batch_size=100,
):
    model = load_model(model_file)
    input_l = Input(
        shape=tuple(model.layers[conv_layer].input.get_shape().as_list()[1:])
    )
    l = input_l
    N_layers = len(model.layers)

    for i in range(conv_layer, N_layers):
        l = model.layers[i](l)

    l.activation = activations.linear
    model_revised = Model(input_l, l)

    batch_grads = []

    for batch_start in range(0, input_data.shape[0], batch_size):
        if (batch_start + batch_size) > input_data.shape[0]:
            batch_in = input_data[batch_start:, :, :]
        else:
            batch_in = input_data[batch_start:(batch_start + batch_size), :, :]
        ngrads = []
        for _ in range(n):
            input_seq = tf.cast(
                batch_in + np.random.normal(
                    size=batch_in.shape, loc=0, scale=0.2
                ),
                tf.float32,
            )

            with tf.GradientTape() as tape:
                tape.watch(input_seq)
                layer_output = model_revised(input_seq, training=False)
                if len(pos_index) == 0:
                    loss = tf.reduce_mean(
                        tf.reduce_sum(
                            [-layer_output[:, ni] for ni in neg_index], axis=1
                        )
                    )
                    loss += lambda_reg * tf.abs(input_seq)
                elif len(neg_index) == 0:
                    loss = tf.reduce_mean(
                        tf.reduce_sum([layer_output[:, pi]
                                       for pi in pos_index], axis=1)
                    )
                    loss += lambda_reg * tf.abs(input_seq)
                else:
                    pos_classes = tf.reduce_sum(
                        [layer_output[:, pi] for pi in pos_index], axis=1
                    )
                    neg_classes = tf.reduce_sum(
                        [layer_output[:, ni] for ni in neg_index], axis=1
                    )
                    loss = tf.reduce_mean(pos_classes - neg_classes)
                    loss += lambda_reg * tf.abs(input_seq)

            # compute the gradient of the input picture wrt this loss
            grads = tape.gradient(loss, input_seq)

            # normalization trick: we normalize the gradient
            grads /= tf.sqrt(tf.reduce_mean(tf.square(grads))) + 1e-5
            ngrads.append(grads)

        batch_grads.append(tf.reduce_mean(ngrads, axis=0))

    del model_revised
    del model
    return np.concatenate(batch_grads, axis=0)
