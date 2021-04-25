from keras import backend as K
from keras.models import Model
from keras.layers import Input
from keras import activations
import keras
import numpy as np
from keras import activations
import tensorflow as tf #tensorflow back end assumed

# saliency code modified from
# https://github.com/keras-team/keras/blob/master/examples/conv_filter_visualization.py
#All contributions by François Chollet:
#Copyright (c) 2015 - 2019, François Chollet.
#All rights reserved.

#All other contributions:
#Copyright (c) 2015 - 2019, the respective contributors.
#All rights reserved.
def saliency(conv_layer,model_file,
             input_data,pos_index,neg_index,
             n=5,lambda_reg=0.1,batch_size=100):
    model = keras.models.load_model(model_file)
    inp = model.input                                           # input placeholder
    outputs = [layer.output for layer in model.layers]          # all layer outputs
    functors = [K.function([inp], [out]) for out in outputs]    # evaluation functions
    input_l = Input(shape=tuple(model.layers[conv_layer].input.get_shape().as_list()[1:]))
    l = input_l
    N_layers = len(model.layers)

    for i in range(conv_layer, N_layers):
        l = model.layers[i](l)
        
    l.activation = activations.linear
    model_revised = Model(input_l, l)

    batch_grads = []

    for batch_start in range(0,input_data.shape[0],batch_size):
        if (batch_start+batch_size) > input_data.shape[0]:
            batch_in = input_data[batch_start:,:,:]
        else:
            batch_in = input_data[batch_start:(batch_start+batch_size),:,:]
        ngrads = []
        for _ in range(n):
            input_seq = tf.Variable(batch_in+np.random.normal(size=batch_in.shape,loc=0,scale=0.2),dtype=tf.float32)
            layer_output = model_revised(input_seq, training=False)
    
            if len(pos_index) == 0:
                loss = K.mean(K.sum([-layer_output[:,ni] for ni in neg_index],axis=1)) + lambda_reg*K.abs(input_seq)
            elif len(neg_index) == 0:
                loss = K.mean(K.sum([layer_output[:,pi] for pi in pos_index],axis=1)) + lambda_reg*K.abs(input_seq)
            else:
                loss = K.mean(K.sum([layer_output[:,pi] for pi in pos_index],axis=1)-K.sum([layer_output[:,ni] for ni in neg_index],axis=1)) + lambda_reg*K.abs(input_seq)

            with tf.GradientTape() as tape:
                # compute the gradient of the input picture wrt this loss
                grads = tape.gradient(loss, input_seq)
                grads = tf.convert_to_tensor(grads,dtype=tf.float32)
                # normalization trick: we normalize the gradient
                print(type(grads))
                grads /= (K.sqrt(K.mean(K.square(grads))) + 1e-5)
                ngrads.append(grads)
                
        batch_grads.append(np.mean(ngrads,axis=0))

    del model_revised
    del model
    return np.concatenate(batch_grads,axis=0)
