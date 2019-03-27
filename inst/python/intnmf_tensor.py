#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
#import numpy as np
import os
import tensorflow as tf
#import argparse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "1"  # or any {'0', '1', '2'}


"""
Created on Thu Mar 21 2019

@author: Andres Quintero
"""

def iNMF(matrices, rank, iterations, lamb, sparcity=0, nrep=200, steps=100, stop_threshold=40):
    
    nviews = len(matrices)
    
    Xs_tensor = []
    Xs_shapes = []
    for i in range(nviews):
        Xs_shapes.append(matrices[i].shape)
        Xs_tensor.append(tf.constant(matrices[i], name = ("X" + str(i)), dtype=tf.float32))

    
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                             Initialize W matrices                         ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    initializer = tf.random_uniform_initializer(minval=0, maxval=2)
    
    Ws_tensor = [tf.Variable(initializer(shape=[Xs_shapes[i][0], rank]), 
                             name=("W" + str(i))) for i in range(nviews)]
        
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                             Initialize H matrix                           ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    
    H = tf.Variable(initializer(shape=[rank, Xs_shapes[0][1]]), name="H")
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                  Initialize view specific H matrices                      ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    Hviews_tensor = [tf.Variable(initializer(shape=[rank, Xs_shapes[i][1]]), 
                                 name= ("Hview" + str(i))) for i in range(nviews)]
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                Save initial max exposures in H matrices                   ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    Hs_summed = [tf.math.add(H, Hv) for Hv in Hviews_tensor]
    oldExposures = tf.concat([tf.math.argmax(Hsummed, axis=0) for Hsummed in Hs_summed ], 0)
    const = 0
    
    for inner in range(iterations):
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                              Update H                                     ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        with tf.name_scope('update_H') as scope:
            num = tf.reduce_sum([tf.matmul(Xs_tensor[j], Ws_tensor[j], transpose_a=True) for j in range(nviews)], 0)
            den = tf.reduce_sum([(tf.matmul(tf.reduce_sum([H, Hviews_tensor[j]], 0), tf.matmul(Ws_tensor[j], Ws_tensor[j], transpose_a=True), transpose_a=True)) for j in range(nviews)], 0)
            H_new = tf.math.multiply(H,  tf.transpose(tf.math.divide(num, den)))
            
            H_new  = tf.where(tf.math.is_nan(H_new ), tf.zeros_like(H_new ), H_new )
            update_H = H.assign(H_new)
            
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                                 Update W                                  ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##                        
        for i in range(nviews):
            X, W, Hv = Xs_tensor[i], Ws_tensor[i], Hviews_tensor[i]
            
            #save W for calculating delta with the updated W
            W_old = tf.Variable(W, name = ("W_old" + str(i)))
            save_W = W_old.assign(W)
            
            with tf.name_scope(('update_W' + str(i))) as scope:
                #update operation for W (after updating H)
                HHv     = tf.reduce_sum([H, Hv], 0)
                HHvHHvW = tf.matmul(tf.matmul(HHv, HHv, transpose_b=True), W, transpose_b=True)
                HvHvW   = tf.matmul(tf.matmul(Hv, Hv, transpose_b=True), W, transpose_b=True) + lamb
                denH = tf.reduce_sum([HHvHHvW, HvHvW], 0) + sparcity
                HHvX_denH =tf.math.divide(tf.matmul(HHv,  X, transpose_b=True), denH)
                HHvX_denH = tf.where(tf.math.is_nan(HHvX_denH), tf.zeros_like(HHvX_denH), HHvX_denH)
    
                W_new = W * tf.transpose(HHvX_denH)             
                W_new  = tf.where(tf.math.is_nan(W_new ), tf.zeros_like(W_new ), W_new )
                update_W = W.assign(W_new)
    
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                    Update view Specific Hs                                ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
                
        count_Hvs_updated = []
        for i in range(nviews):
            count_Hv_updated_tmp = tf.Variable(0, name = ("updateoldHview" + str(i)), dtype=tf.int32)
            count_Hvs_updated.append(count_Hv_updated_tmp )
        
                
        for i in range(nviews):
            X, W, Hv = Xs_tensor[i], Ws_tensor[i], Hviews_tensor[i]
            
            with tf.name_scope(('update_View_H' + str(i))) as scope:
                #update operation for view specific H (after updating H and Ws)
                
                WW = tf.matmul(W, W, transpose_a=True)
                HHvWW = tf.matmul(tf.reduce_sum([H, Hv], 0), WW, transpose_a=True)
                lambHvWW = tf.matmul((Hv * lamb), WW, transpose_a=True)
                denHv = tf.reduce_sum([HHvWW, lambHvWW], 0) 
                XH_denHv = tf.math.divide(tf.matmul(X, W, transpose_a=True), denHv)
                XH_denHv = tf.where(tf.math.is_nan(XH_denHv), tf.zeros_like(XH_denHv), XH_denHv)
    
                Hv_new = Hv * tf.transpose(XH_denHv)
                Hv_new  = tf.where(tf.math.is_nan(Hv_new ), tf.zeros_like(Hv_new ), Hv_new )
                update_Hv = Hv.assign(Hv_new)
    
        Hs_summed = [tf.math.add(H, Hv) for Hv in Hviews_tensor]
        newExposures = tf.concat([tf.math.argmax(Hsummed, axis=0) for Hsummed in Hs_summed ], 0)
        
        if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
            oldExposures = newExposures
            const = 0
        else:
            const += 1
            if const == stop_threshold:
                break
        
    frobNorm = []
    for i in range(nviews):
        HHv = tf.math.add(H, Hviews_tensor[i])
        fb = tf.linalg.norm(Xs_tensor[i] - tf.matmul(Ws_tensor[i], HHv)) / tf.linalg.norm(Xs_tensor[i])
        frobNorm.append(fb.numpy())
    
    Ws_num = [Wi.numpy() for Wi in Ws_tensor]
    H_num  = H.numpy()
    Hv_num = [Hv.numpy() for Hv in Hviews_tensor] 
        
    return Ws_num, H_num, Hv_num, inner+1, frobNorm 
