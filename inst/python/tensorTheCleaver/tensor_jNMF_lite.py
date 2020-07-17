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




#%%

##---------------------------------------------------------------------------##
##                       Reading matrices to factorize                       ##
##---------------------------------------------------------------------------##
# RNAseq and ATACseq
#filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/rna_norm_mat.csv", "/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/atac_norm_mat.csv"]
#
#Xs_list = [np.loadtxt(path) for path in filenames]
#Xs_list 

#%%
#Xs_list = [np.transpose(X)for X in Xs_list]


#%%

#[X.shape for X in Xs_list]


#%%

#@tf.function
##---------------------------------------------------------------------------##
##                               Update H                                    ##
##---------------------------------------------------------------------------##
def jNMF_uptade_H(Xs, Ws, H, Nviews):
    num   = tf.reduce_sum([tf.matmul(Xs[i], Ws[i], transpose_b=True) for i in range(Nviews)], 0)
    den_K = tf.reduce_sum([tf.matmul(Ws[i], Ws[i], transpose_b=True) for i in range(Nviews)], 0)
    den    = tf.matmul(H, den_K, transpose_b=True)
    H_new  = tf.multiply(H,  tf.divide(num, den))
    H_new  = tf.where(tf.math.is_nan(H_new), tf.zeros_like(H_new), H_new)
    #H.assign(H_new)
    return H_new

##---------------------------------------------------------------------------##
##                            Update Ws                                      ##
##---------------------------------------------------------------------------##                        
def jNMF_uptade_Ws(Xs, Ws, H, Nviews, Sp):
    HtH   = tf.matmul(H, H, transpose_a=True)
    for i in range(Nviews):
        den    = tf.matmul(HtH, Ws[i]) + Sp            
        HX_den = tf.divide(tf.matmul(H, Xs[i], transpose_a=True), den)
        W_new  = tf.multiply(Ws[i], HX_den)
        W_new  = tf.where(tf.math.is_nan(W_new), tf.zeros_like(W_new), W_new)
        Ws[i].assign(W_new)
    return Ws

##---------------------------------------------------------------------------##
##                   define objective function                               ##
##---------------------------------------------------------------------------##
def nmf_obj_eval(Xs, Ws, H, Nviews):
    # Frobenius term
    frob_c = []
    for i in range(Nviews):
        frob_c.append(tf.linalg.norm(Xs[i] - tf.matmul(H, Ws[i])) / tf.linalg.norm(Xs[i]))
    frob_c = tf.reduce_sum(frob_c)
    
    return frob_c 
        

def jNMF_tensor_py(matrix_list, rank, n_initializations, iterations, Sp, stop_threshold=40):
    matrix_list = [np.transpose(Xv) for Xv in matrix_list]

    #K = len(matrix_list)
    Nviews = len(matrix_list)
    N = matrix_list[0].shape[0]
    Ms = [Xv.shape[1] for Xv in matrix_list]
    
    # Initialization Metrics 
    frobNorm = []
    iter_to_conv = []
    W_eval = []
    
    # X matrices to tensor constant
    Xs = [tf.constant(matrix_list[i], name = ("X" + str(i)), dtype=tf.float32) for i in range(Nviews)]

    ##-----------------------------------------------------------------------##
    ##                              N inits                                  ##
    ##-----------------------------------------------------------------------##
    # cycle through n initializations and choose best factorization
    for init_n in range(n_initializations):
    
        ##-------------------------------------------------------------------##
        ##                     Initialize W matrices                         ##
        ##-------------------------------------------------------------------##
        initializer = tf.random_uniform_initializer(minval=0, maxval=2)
    
        Ws = [tf.Variable(initializer(shape=[rank, Ms[i]]),
                          name=("W" + str(i))) for i in range(Nviews)]
    
    
        ##-------------------------------------------------------------------##
        ##                         Initialize H matrix                       ##
        ##-------------------------------------------------------------------##    
        H = tf.Variable(initializer(shape=[N, rank]), name="H")
    
        ##-------------------------------------------------------------------##
        ##        Save initial max exposures in H matrices                   ##
        ##-------------------------------------------------------------------##
        oldExposures = tf.math.argmax(H, axis=1)
        const = 0        

        ##-------------------------------------------------------------------##
        ##                   Start matrix factorization                      ##
        ##-------------------------------------------------------------------##
        for inner in range(iterations):
            ##---------------------------------------------------------------##
            ##                          Update H                             ##
            ##---------------------------------------------------------------##
            H = jNMF_uptade_H(Xs, Ws, H, Nviews)
    
            ##---------------------------------------------------------------##
            ##                            Update Ws                          ##
            ##---------------------------------------------------------------##   
            Ws = jNMF_uptade_Ws(Xs, Ws, H, Nviews, Sp)
    
            ##---------------------------------------------------------------##
            ##                    Evaluate Convergence                       ##
            ##---------------------------------------------------------------##        
            newExposures = tf.math.argmax(H, axis=1)
    
            if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
                oldExposures = newExposures
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    break
        
        ##-------------------------------------------------------------------##
        ##                     Initialize frob error                         ##
        ##-------------------------------------------------------------------##
        if init_n == 0 :
            Best_frob = nmf_obj_eval(Xs, Ws, H, Nviews)
            Best_H    = H
            Best_Ws   = Ws
            
        #Hvs = [tf.Variable(initializer(shape=[N, rank]),
        #                   name= ("Hview" + str(i))) for i in range(K)]
        
        ##-------------------------------------------------------------------##
        ##         Evaluate if best factorization initialization             ##
        ##-------------------------------------------------------------------##
        frobInit = nmf_obj_eval(Xs, Ws, H, Nviews)
        #        frobNorm_init = []
        #        for i in range(Nviews):
        #            #Ht = tf.add(H, Hvs[i])
        #            fb = tf.linalg.norm(Xs[i] - tf.matmul(H, Ws[i])) / tf.linalg.norm(Xs[i])
        #            frobNorm_init.append(fb)
        #        frobNorm_init = tf.reduce_sum(frobNorm_init)
        frobNorm.append(frobInit)
        iter_to_conv.append(inner+1)
        W_eval.append(tf.concat(Ws, 1))
        
        if frobInit < Best_frob :
            #print('is less')
            Best_frob = frobInit
            Best_H    = H
            Best_Ws   = Ws
        #x = nmf_obj_eval(Xs, Ws, H, Nviews)
        #print("Best frob:", x.numpy())
        #print("Current frob", frobInit.numpy())
        #fb = tf.reduce_sum(fb, 0)
    
    ##-----------------------------------------------------------------------##
    ##             Convert to numpy, transpose and return                    ##
    ##-----------------------------------------------------------------------##
    Ws_num = [Wi.numpy().T for Wi in Best_Ws]
    H_num  = Best_H.numpy().T

    frobNorm = [i.numpy() for i in frobNorm]
    W_eval_num  = [i.numpy().T for i in W_eval]

    return Ws_num, H_num, iter_to_conv, frobNorm, W_eval_num
    

#@tf.function
#def jNMF_run_one_init(Xs, rank, Ms, N, Nviews, Sp, iterations, stop_threshold):
#    ##-------------------------------------------------------------------##
#    ##                     Initialize W matrices                         ##
#    ##-------------------------------------------------------------------##
#    initializer = tf.random_uniform_initializer(minval=0, maxval=2)
#
#    Ws = [tf.Variable(initializer(shape=[rank, Ms[i]]),
#                      name=("W" + str(i))) for i in range(Nviews)]
#    ##-------------------------------------------------------------------##
#    ##                         Initialize H matrix                       ##
#    ##-------------------------------------------------------------------##
#    H = tf.Variable(initializer(shape=[N, rank]), name="H")
#
#    ##-------------------------------------------------------------------##
#    ##        Save initial max exposures in H matrices                   ##
#    ##-------------------------------------------------------------------##
#    oldExposures = tf.math.argmax(H, axis=1)
#    const = 0        
#
#    ##-------------------------------------------------------------------##
#    ##                   Start matrix factorization                      ##
#    ##-------------------------------------------------------------------##
#    for inner in range(iterations):
#        ##---------------------------------------------------------------##
#        ##                          Update H                             ##
#        ##---------------------------------------------------------------##
#        #H.assign(jNMF_uptade_H(Xs, Ws, H, Nviews))
#        H = jNMF_uptade_H(Xs, Ws, H, Nviews)
#        ##---------------------------------------------------------------##
#        ##                            Update Ws                          ##
#        ##---------------------------------------------------------------##   
#        Ws = jNMF_uptade_Ws(Xs, Ws, H, Nviews, Sp)
#        ##---------------------------------------------------------------##
#        ##                    Evaluate Convergence                       ##
#        ##---------------------------------------------------------------##        
#        newExposures = tf.math.argmax(H, axis=1)
#
#        if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
#            oldExposures = newExposures
#            const = 0
#        else:
#            const += 1
#            if const == stop_threshold:
#                break
#    return Ws, H, inner+1
#    
#    
#    
#    
#def jNMF_tensor_py(matrix_list, rank, n_initializations, iterations, Sp, stop_threshold=40):
#    matrix_list = [np.transpose(Xv) for Xv in matrix_list]
#
#    #K = len(matrix_list)
#    Nviews = len(matrix_list)
#    N = matrix_list[0].shape[0]
#    Ms = [Xv.shape[1] for Xv in matrix_list]
#    frobNorm = []
#    iter_to_conv = []
#    
#    # X matrices to tensor constant
#    Xs = [tf.constant(matrix_list[i], name = ("X" + str(i)), dtype=tf.float32) for i in range(Nviews)]
#
#    ##-----------------------------------------------------------------------##
#    ##                              N inits                                  ##
#    ##-----------------------------------------------------------------------##
#    # cycle through n initializations and choose best factorization
#    for init_n in range(n_initializations):
#        Ws, H, inner = jNMF_run_one_init(Xs, rank, Ms, N, Nviews, Sp, iterations, stop_threshold)
#        
#        ##-------------------------------------------------------------------##
#        ##                     Initialize frob error                         ##
#        ##-------------------------------------------------------------------##
#        if init_n == 0 :
#            Best_frob = nmf_obj_eval(Xs, Ws, H, Nviews)
#            Best_H    = H
#            Best_Ws   = Ws
#    
#        ##-------------------------------------------------------------------##
#        ##         Evaluate if best factorization initialization             ##
#        ##-------------------------------------------------------------------##
#        frobInit = nmf_obj_eval(Xs, Ws, H, Nviews)
#        frobNorm.append(frobInit)
#        iter_to_conv.append(inner)
#        #        if frobInit < Best_frob :
#        #            print('is less')
#        #            Best_frob = frobInit
#        #            Best_H    = H
#        #            Best_Ws   = Ws
#        #        x = nmf_obj_eval(Xs, Ws, H, Nviews)
#        #        print("Best frob:", x.numpy())
#        print("Current frob", frobInit.numpy())
#        #fb = tf.reduce_sum(fb, 0)
#    
#    ##-----------------------------------------------------------------------##
#    ##             Convert to numpy, transpose and return                    ##
#    ##-----------------------------------------------------------------------##
#    frobNorm = [i.numpy() for i in frobNorm]
#    #frobNorm.numpy()
#    #    Ws_num = [Wi.numpy().T for Wi in Ws]
#    #    H_num  = H.numpy().T
#    Ws_num = [Wi.numpy().T for Wi in Best_Ws]
#    H_num  = Best_H.numpy().T
#
#    return Ws_num, H_num, iter_to_conv, frobNorm

#%%

#jNMF_tensor_py(matrix_list = Xs_list,
#            rank              = 10,
#            n_initializations = 10,
#            iterations        = 1000,
#            Sp                = 0, 
#            stop_threshold    = 10)    



#%%


#
#iNMF(matrix_list = Xs_list,
#     rank           = 10,
#     iterations     = 1000,
#     L              = 1,
#     Sp             = 0, 
#     stop_threshold = 40)    


#%%
