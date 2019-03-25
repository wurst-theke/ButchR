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


##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                Define Function to run NMF in tensorflow 2                 ##
##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

def NMF(matrix,r,iterations, stop_threshold=40):
    
    # NMF in tensorflow
    n = matrix.shape[0]
    m = matrix.shape[1]     
    X = tf.constant(matrix, name = ("X"), dtype=tf.float32)

    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                   Initialize W and H matrices                         ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    initializer = tf.random_uniform_initializer(minval=0, maxval=1)
    
    H = tf.Variable(initializer(shape=[r, m]), name="H")
    W = tf.Variable(initializer(shape=[n, r]), name="W")
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                          Run NMF                                      ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

    const = 0
    oldExposures = tf.math.argmax(H, axis=0)

    for i in range(iterations):
        # update H
        WTX  = tf.matmul(W, X, transpose_a=True)
        WTW  = tf.matmul(W, W, transpose_a=True)
        WTWH = tf.matmul(WTW, H, transpose_a=True)
        newH = tf.math.divide(tf.math.multiply(H, WTX), WTWH)
        newH = tf.where(tf.math.is_nan( newH ), tf.zeros_like( newH ), newH )
        update_H = H.assign(newH)

        # update W
        XHT  = tf.matmul(X, H, transpose_b=True)
        WH   = tf.matmul(W, H)
        WHHT = tf.matmul(WH, H, transpose_b=True)
        newW = tf.math.divide(tf.math.multiply(W, XHT), WHHT)
        newW = tf.where(tf.math.is_nan( newW ), tf.zeros_like( newW ), newW )
        update_W = W.assign(newW)
        
        # Evaluate convergence
        
        newExpo = tf.math.argmax(H, axis=0)
        if tf.reduce_all(tf.math.equal(oldExposures, newExpo)).__invert__():
            oldExposures = newExpo
            const = 0
        else:
            const += 1
            #print(f'new eval diff, {const}')
            if const == stop_threshold:
                #print(f"NMF converged after {i} iterations")
                break
    
    frobNorm = tf.linalg.norm(X - tf.matmul(W, H)) / tf.linalg.norm(X)
    return W.numpy(),H.numpy(), i, frobNorm.numpy()
