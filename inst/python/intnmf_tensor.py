#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import os
import tensorflow as tf
#import argparse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "1"  # or any {'0', '1', '2'}


"""
Created on Thu Mar 21 2019

@author: Andres Quintero
"""


def iNMF(matrix_list, rank, iterations, L, Sp=0, stop_threshold=40):
    matrix_list = [np.transpose(Xv) for Xv in matrix_list]
    
    K = len(matrix_list)
    N = matrix_list[0].shape[0]
    Ms = [Xv.shape[1] for Xv in matrix_list]
    
    # X matrices to tensor constant
    Xs = [tf.constant(matrix_list[i], name = ("X" + str(i)), dtype=tf.float32) for i in range(K)]
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                         Initialize W matrices                         ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    initializer = tf.random_uniform_initializer(minval=0, maxval=2)
    
    Ws = [tf.Variable(initializer(shape=[rank, Ms[i]]),
                      name=("W" + str(i))) for i in range(K)]
        
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                             Initialize H matrix                       ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    
    H = tf.Variable(initializer(shape=[N, rank]), name="H")
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                  Initialize view specific H matrices                  ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    Hvs = [tf.Variable(initializer(shape=[N, rank]),
                       name= ("Hview" + str(i))) for i in range(K)]

    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##            Save initial max exposures in H matrices                   ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    Hts = [tf.add(H, Hv) for Hv in Hvs]
    oldExposures = tf.concat([tf.math.argmax(Ht, axis=1) for Ht in Hts], 0)
    const = 0
    
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                       Start matrix factorization                      ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    for inner in range(iterations):
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                          Update H                                 ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        num   = tf.reduce_sum([tf.matmul(Xs[i], Ws[i], transpose_b=True) for i in range(K)], 0)
        den_K = []
        for i in range(K):
            Ht = tf.reduce_sum([H, Hvs[i]], 0)
            WW = tf.matmul(Ws[i], Ws[i], transpose_b=True)
            den_K.append(tf.matmul(Ht, WW))
        den    = tf.reduce_sum(den_K, 0)
        H_new  = tf.multiply(H,  tf.divide(num, den))
        H_new  = tf.where(tf.math.is_nan(H_new), tf.zeros_like(H_new), H_new)
        H.assign(H_new)
        
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                            Update Ws                              ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##                        

        for i in range(K):
            Ht     = tf.reduce_sum([H, Hvs[i]], 0)
            HtHt   = tf.matmul(Ht, Ht, transpose_a=True)
            
            HvHv   = tf.matmul(Hvs[i], Hvs[i], transpose_a=True)
            HsHs   = tf.reduce_sum([HtHt, tf.multiply(L,  HvHv)], 0)            
            den    = tf.matmul(HsHs, Ws[i]) + Sp
            HX_den = tf.divide(tf.matmul(Ht, Xs[i], transpose_a=True), den)
            W_new  = tf.multiply(Ws[i], HX_den)
            W_new  = tf.where(tf.math.is_nan(W_new), tf.zeros_like(W_new), W_new)
            Ws[i].assign(W_new)
            
        #    WV = W+Vs[i]
        #    den = (WV.T.dot(WV) + L*Vs[i].T.dot(Vs[i])).dot(Hs[i]) + Sp
        #    Hs[i] *= (WV.T.dot(Xs[i]))/den
    
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                    Update view Specific Hs                        ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        for i in range(K):
            Ht     = tf.reduce_sum([H, tf.multiply((1+L), Hvs[i])], 0)
            WW     = tf.matmul(Ws[i], Ws[i], transpose_b=True)
            den    = tf.matmul(Ht, WW)
            XW     = tf.matmul(Xs[i], Ws[i], transpose_b=True)
            XW_den = tf.divide(XW, den)
            Hv_new = tf.multiply(Hvs[i], XW_den)
            Hv_new = tf.where(tf.math.is_nan(Hv_new), tf.zeros_like(Hv_new), Hv_new)
            Hvs[i].assign(Hv_new)
    
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                    Evaluate Convergence                           ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        Hts = [tf.add(H, Hv) for Hv in Hvs]
        newExposures = tf.concat([tf.math.argmax(Ht, axis=1) for Ht in Hts], 0)
        
        if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
            oldExposures = newExposures
            const = 0
        else:
            const += 1
            if const == stop_threshold:
                break
            
    frobNorm = []
    for i in range(K):
        Ht = tf.add(H, Hvs[i])
        fb = tf.linalg.norm(Xs[i] - tf.matmul(Ht, Ws[i])) / tf.linalg.norm(Xs[i])
        frobNorm.append(fb.numpy())
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##             Convert to numpy, transpose and return                    ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    
    Ws_num = [Wi.numpy().T for Wi in Ws]
    H_num  = H.numpy().T
    Hvs_num = [Hv.numpy().T for Hv in Hvs] 
        
    return Ws_num, H_num, Hvs_num, inner+1, frobNorm 
