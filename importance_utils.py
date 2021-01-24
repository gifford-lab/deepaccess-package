from keras import backend as K
from keras.models import Model
from keras.layers import Input
from keras import activations
import keras
import numpy as np
from keras import activations

# saliency code modified from
# https://github.com/keras-team/keras/blob/master/examples/conv_filter_visualization.py
#All contributions by François Chollet:
#Copyright (c) 2015 - 2019, François Chollet.
#All rights reserved.

#All other contributions:
#Copyright (c) 2015 - 2019, the respective contributors.
#All rights reserved.
def saliency(conv_layer,model_file,input_data,pos_index,neg_index,n=5,lambda_reg=0.1,batch_size=100):
    K.set_learning_phase(0)
    model = keras.models.load_model(model_file)
    print(model.summary())
    inp = model.input                                           # input placeholder
    outputs = [layer.output for layer in model.layers]          # all layer outputs
    functors = [K.function([inp], [out]) for out in outputs]    # evaluation functions
    input_l = Input(shape=tuple(model.layers[conv_layer].input.get_shape().as_list()[1:]))
    l = input_l
    for i in range(conv_layer, len(model.layers)):
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
            layer_outs = [func([batch_in + np.random.normal(size=batch_in.shape)]) for func in functors]
            layer_output = model_revised.layers[-1].get_output_at(1)
            input_seq = model_revised.input
            if len(pos_index) == 0:
                loss = K.mean(K.sum([-layer_output[:,ni] for ni in neg_index],axis=1)) + lambda_reg*K.abs(input_seq)
            elif len(neg_index) == 0:
                loss = K.mean(K.sum([layer_output[:,pi] for pi in pos_index],axis=1)) + lambda_reg*K.abs(input_seq)
            else:
                loss = K.mean(K.sum([layer_output[:,pi] for pi in pos_index],axis=1)-K.sum([layer_output[:,ni] for ni in neg_index],axis=1)) + lambda_reg*K.abs(input_seq)
            # compute the gradient of the input picture wrt this loss
            grads = K.gradients(loss, input_seq)[0]
        
            # normalization trick: we normalize the gradient
            grads /= (K.sqrt(K.mean(K.square(grads))) + 1e-5)

            # this function returns the loss and grads given the input picture
            iterate = K.function([input_seq], [loss, grads])
            if conv_layer == 0:
                loss,gradient = iterate([batch_in+np.random.normal(size=batch_in.shape,loc=0,scale=0.2)])
            else:
                loss,gradient = iterate(layer_outs[conv_layer-1])
            ngrads.append(gradient)
        batch_grads.append(np.mean(ngrads,axis=0))

    del model_revised
    del model
    return np.concatenate(batch_grads,axis=0)

def trace_to_conv_layer(conv_layer,model_file,class_index,seed_seqs,lambda_reg):
    K.set_learning_phase(0)
    model = keras.models.load_model(model_file)
    inp = model.input                                           # input placeholder
    outputs = [layer.output for layer in model.layers]          # all layer outputs
    functors = [K.function([inp], [out]) for out in outputs]    # evaluation functions
    test = seed_seqs
    layer_outs = [func([test]) for func in functors]
    for i in range(0,conv_layer):
        model.layers.pop()
    input_l = Input(shape=tuple(model.layers[conv_layer].input.get_shape().as_list()[1:]))
    seqs_shape = tuple([seed_seqs.shape[0]]+model.layers[conv_layer].input.get_shape().as_list()[1:])
    l = input_l
    for i in range(conv_layer, len(model.layers)):
        l = model.layers[i](l)
    model_revised = Model(input_l, l)
    layer_output = model_revised.layers[-1].get_output_at(1)
    input_seq = model_revised.input
    if class_index==0: other_index=1
    if class_index==1: other_index=0
    loss = K.mean(layer_output[:,class_index] - layer_output[:,other_index]) - lambda_reg*K.pow(K.sum(K.pow(K.abs(input_seq), 2.0)), 1. / 2.0)

    # compute the gradient of the input picture wrt this loss
    grads = K.gradients(loss, input_seq)[0]
    
    # normalization trick: we normalize the gradient
    grads /= (K.sqrt(K.mean(K.square(grads))) + 1e-5)

    # this function returns the loss and grads given the input picture
    iterate = K.function([input_seq], [loss, grads])

    changes = np.zeros(seqs_shape)
    scores = []
    scores_new = []
    for j in range(seed_seqs.shape[0]):
        input_seq_data = np.copy(np.reshape(layer_outs[conv_layer-1][0][j,:,:],(1,seqs_shape[1],seqs_shape[2])))
        #og_data = np.copy(input_seq_data)
        scores.append(model_revised.predict(input_seq_data))
        step=0.1
        #grads = []
        # run gradient ascent for 20 steps
        #for i in range(100):
        loss_value, grads_value = iterate([input_seq_data])
        changes[j,:,:] = grads_value 
        input_seq_data += grads_value * step
        #changes[j,:,:] = input_seq_data
        
        scores_new.append(model_revised.predict(input_seq_data))
        
    return {"grads":changes,"scores":scores,"scores_grads":scores_new}
