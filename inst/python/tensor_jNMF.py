#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 11:36:20 2020

@author: andresq
"""

from __future__ import division
import numpy as np
import os
import tensorflow as tf

os.environ['TF_CPP_MIN_LOG_LEVEL'] = "1"  # or any {'0', '1', '2'}




##%%
#
###–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
###                       Reading matrices to factorize                       ##
###–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##os.chdir("/Users/andresq/phd/main_project/src/")
#
## Toy RNAseq and Toy ATACseq
##filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/toy_rna_norm_mat.csv", "/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/toy_atac_norm_mat.csv"]
## Toy RNAseq
##filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/toy_rna_norm_mat.csv"]
## Toy RNAseq and Toy RNAseq
##filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/toy_rna_norm_mat.csv", "/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/toy_rna_norm_mat.csv"]
## RNAseq and ATACseq
#filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/rna_norm_mat.csv", "/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/atac_norm_mat.csv"]
## RNAseq 
##filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/rna_norm_mat.csv"]
## RNAseq and RNAseq duplicated
##filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/rna_norm_mat.csv", "/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/rna_norm_mat.csv"]
## 3 views RNAseq and ATACseq
##filenames.append("/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/toy_rna_norm_mat.csv")
##filenames.extend(filenames)
#
#Xs_list = [np.loadtxt(path) for path in filenames]
#Xs_list 

#%%
#Xs_list = [np.transpose(X)for X in Xs_list]


#%%

#[X.shape for X in Xs_list]


#%%


def jNMF_tensor_py(matrix_list, rank, iterations, Sp, stop_threshold=40):
    matrix_list = [np.transpose(Xv) for Xv in matrix_list]

    #K = len(matrix_list)
    Nviews = len(matrix_list)
    N = matrix_list[0].shape[0]
    Ms = [Xv.shape[1] for Xv in matrix_list]

    # X matrices to tensor constant
    Xs = [tf.constant(matrix_list[i], name = ("X" + str(i)), dtype=tf.float32) for i in range(Nviews)]

    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                         Initialize W matrices                         ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    initializer = tf.random_uniform_initializer(minval=0, maxval=2)

    Ws = [tf.Variable(initializer(shape=[rank, Ms[i]]),
                      name=("W" + str(i))) for i in range(Nviews)]


    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                             Initialize H matrix                       ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

    H = tf.Variable(initializer(shape=[N, rank]), name="H")

    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                  Initialize view specific H matrices                  ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    #Hvs = [tf.Variable(initializer(shape=[N, rank]),
    #                   name= ("Hview" + str(i))) for i in range(K)]

    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##            Save initial max exposures in H matrices                   ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    #Hts = [tf.add(H, Hv) for Hv in Hvs]
    #oldExposures = tf.concat([tf.math.argmax(Ht, axis=1) for Ht in Hts], 0)
    oldExposures = tf.math.argmax(H, axis=1)
    const = 0        

    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##                       Start matrix factorization                      ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    for inner in range(iterations):
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                          Update H                                 ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        num   = tf.reduce_sum([tf.matmul(Xs[i], Ws[i], transpose_b=True) for i in range(Nviews)], 0)
        den_K = tf.reduce_sum([tf.matmul(Ws[i], Ws[i], transpose_b=True) for i in range(Nviews)], 0)
        den    = tf.matmul(H, den_K, transpose_b=True)
        H_new  = tf.multiply(H,  tf.divide(num, den))
        H_new  = tf.where(tf.math.is_nan(H_new), tf.zeros_like(H_new), H_new)
        H.assign(H_new)

        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                            Update Ws                              ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##                        
        HtH   = tf.matmul(H, H, transpose_a=True)
        for i in range(Nviews):
            den    = tf.matmul(HtH, Ws[i]) + Sp            
            HX_den = tf.divide(tf.matmul(H, Xs[i], transpose_a=True), den)
            W_new  = tf.multiply(Ws[i], HX_den)
            W_new  = tf.where(tf.math.is_nan(W_new), tf.zeros_like(W_new), W_new)
            Ws[i].assign(W_new)


        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
        ##                    Evaluate Convergence                           ##
        ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##        
        newExposures = tf.math.argmax(H, axis=1)

        if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
            oldExposures = newExposures
            const = 0
        else:
            const += 1
            if const == stop_threshold:
                break

    frobNorm = []
    for i in range(Nviews):
        #Ht = tf.add(H, Hvs[i])
        fb = tf.linalg.norm(Xs[i] - tf.matmul(H, Ws[i])) / tf.linalg.norm(Xs[i])
        frobNorm.append(fb.numpy())
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
    ##             Convert to numpy, transpose and return                    ##
    ##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

    Ws_num = [Wi.numpy().T for Wi in Ws]
    H_num  = H.numpy().T

    return Ws_num, H_num, inner+1, frobNorm

#%%

#jNMF_tensor(matrix_list = Xs_list,
#            rank           = 10,
#            iterations     = 1000,
#            Sp             = 0, 
#            stop_threshold = 40)    



#%%


#
#iNMF(matrix_list = Xs_list,
#     rank           = 10,
#     iterations     = 1000,
#     L              = 1,
#     Sp             = 0, 
#     stop_threshold = 40)    


#%%
