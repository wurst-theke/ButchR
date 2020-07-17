#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import os
import tensorflow as tf
# import argparse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "1"  # or any {'0', '1', '2'}


"""
Created on Thu Mar 21 2019

@author: Andres Quintero
"""

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
def iNMF_uptade_H(Xs, Ws, H, Hvs, Nviews):
    num   = tf.reduce_sum([tf.matmul(Xs[i], Ws[i], transpose_b=True) for i in range(Nviews)], 0)
    den_K = []
    for i in range(Nviews):
        Ht = tf.reduce_sum([H, Hvs[i]], 0)
        WW = tf.matmul(Ws[i], Ws[i], transpose_b=True)
        den_K.append(tf.matmul(Ht, WW))
    den    = tf.reduce_sum(den_K, 0)
    H_new  = tf.multiply(H,  tf.divide(num, den))
    H_new  = tf.where(tf.math.is_nan(H_new), tf.zeros_like(H_new), H_new)
    return H_new

##---------------------------------------------------------------------------##
##                    Update view specficic Ws                               ##
##---------------------------------------------------------------------------##                        
def iNMF_uptade_Ws(Xs, Ws, H, Hvs, Nviews, lamb, Sp):
    for i in range(Nviews):
        Ht     = tf.reduce_sum([H, Hvs[i]], 0)
        HtHt   = tf.matmul(Ht, Ht, transpose_a=True)
    
        HvHv   = tf.matmul(Hvs[i], Hvs[i], transpose_a=True)
        HsHs   = tf.reduce_sum([HtHt, tf.multiply(lamb,  HvHv)], 0)
        den    = tf.matmul(HsHs, Ws[i]) + Sp
        HX_den = tf.divide(tf.matmul(Ht, Xs[i], transpose_a=True), den)
        W_new  = tf.multiply(Ws[i], HX_den)
        W_new  = tf.where(tf.math.is_nan(W_new), tf.zeros_like(W_new), W_new)
        Ws[i].assign(W_new)
    return Ws

##---------------------------------------------------------------------------##
##                    Update view specficic Hs                               ##
##---------------------------------------------------------------------------##                        
def iNMF_uptade_Hs(Xs, Ws, H, Hvs, Nviews, lamb, Sp):
    for i in range(Nviews):
        Ht     = tf.reduce_sum([H, tf.multiply((1 + lamb), Hvs[i])], 0)
        WW     = tf.matmul(Ws[i], Ws[i], transpose_b=True)
        den    = tf.matmul(Ht, WW)
        XW     = tf.matmul(Xs[i], Ws[i], transpose_b=True)
        XW_den = tf.divide(XW, den)
        Hv_new = tf.multiply(Hvs[i], XW_den)
        Hv_new = tf.where(tf.math.is_nan(Hv_new), tf.zeros_like(Hv_new), Hv_new)
        Hvs[i].assign(Hv_new)
    return Hvs

##---------------------------------------------------------------------------##
##                   define objective function                               ##
##---------------------------------------------------------------------------##
def inmf_obj_eval(Xs, Ws, H, Hvs, Nviews, Sp, lamb):
    # Frobenius term
    #    frob_c = []
    #    for i in range(Nviews):
    #        frob_c.append(tf.linalg.norm(Xs[i] - tf.matmul(H, Ws[i])) / tf.linalg.norm(Xs[i]))
    #    frob_c = tf.reduce_sum(frob_c)
    frob_c   = []
    pen_c    = []
    sparse_c = []
    for i in range(Nviews):
        Ht = tf.add(H, Hvs[i])
        #frob_c.append(tf.linalg.norm(Xs[i] - tf.matmul(Ht, Ws[i])) / tf.linalg.norm(Xs[i]))
        frob_ci = tf.linalg.norm(Xs[i] - tf.matmul(Ht, Ws[i]))
        frob_c.append(tf.math.square(frob_ci))
        
        pen_ci = tf.math.square(tf.linalg.norm(tf.matmul(Hvs[i], Ws[i])))
        pen_c.append(lamb * pen_ci)
        
        sparse_c.append(Sp *tf.reduce_sum(Ws[i]))
        
    frob_c   = tf.reduce_sum(frob_c)
    pen_c    = tf.reduce_sum(pen_c)
    sparse_c = tf.reduce_sum(sparse_c)
    #return frob_c 
    return frob_c + pen_c + sparse_c

##---------------------------------------------------------------------------##
##                   define objective function                               ##
##---------------------------------------------------------------------------##
def inmf_max_exp(H, Hvs):
    Hts = [tf.add(H, Hv) for Hv in Hvs]
    max_exposures = tf.concat([tf.math.argmax(Ht, axis=1) for Ht in Hts], 0)
    return max_exposures 

##---------------------------------------------------------------------------##
##                        define main function                               ##
##---------------------------------------------------------------------------##
def iNMF_tensor_py(matrix_list, rank, n_initializations, iterations, Sp, 
                   stop_threshold=40, lamb = 10, **kwargs):
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
        ##                     Initialize H matrix                           ##
        ##-------------------------------------------------------------------##    
        H = tf.Variable(initializer(shape=[N, rank]), name="H")        
        ##-------------------------------------------------------------------##
        ##               Initialize view specific H matrices                 ##
        ##-------------------------------------------------------------------##
        Hvs = [tf.Variable(initializer(shape=[N, rank]),
                           name= ("Hview" + str(i))) for i in range(Nviews)]    
        ##-------------------------------------------------------------------##
        ##        Save initial max exposures in H matrices                   ##
        ##-------------------------------------------------------------------##
        oldExposures = inmf_max_exp(H, Hvs)
        const = 0       



        ##-------------------------------------------------------------------##
        ##                   Start matrix factorization                      ##
        ##-------------------------------------------------------------------##
        for inner in range(iterations):
            ##---------------------------------------------------------------##
            ##                          Update H                             ##
            ##---------------------------------------------------------------##
            H = iNMF_uptade_H(Xs, Ws, H, Hvs, Nviews)
            ##---------------------------------------------------------------##
            ##                   Update view specficic Ws                    ##
            ##---------------------------------------------------------------##   
            Ws = iNMF_uptade_Ws(Xs, Ws, H, Hvs, Nviews, lamb, Sp)
            ##---------------------------------------------------------------##
            ##                    Update view specficic Hs                   ##
            ##---------------------------------------------------------------##                        
            Hvs = iNMF_uptade_Hs(Xs, Ws, H, Hvs, Nviews, lamb, Sp)

    
            ##---------------------------------------------------------------##
            ##                    Evaluate Convergence                       ##
            ##---------------------------------------------------------------##        
            newExposures = inmf_max_exp(H, Hvs)
    
            if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
                oldExposures = newExposures
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    break
            
            #print("Best frob:", inmf_obj_eval(Xs, Ws, H, Hvs, Nviews, Sp, lamb).numpy())
        
        ##-------------------------------------------------------------------##
        ##                     Initialize frob error                         ##
        ##-------------------------------------------------------------------##
        if init_n == 0 :
            Best_frob = inmf_obj_eval(Xs, Ws, H, Hvs, Nviews, Sp, lamb)
            Best_H    = H
            Best_Ws   = Ws
            Best_Hvs   = Hvs
            
        #Hvs = [tf.Variable(initializer(shape=[N, rank]),
        #                   name= ("Hview" + str(i))) for i in range(K)]
        
        ##-------------------------------------------------------------------##
        ##         Evaluate if best factorization initialization             ##
        ##-------------------------------------------------------------------##
        frobInit = inmf_obj_eval(Xs, Ws, H, Hvs, Nviews, Sp, lamb)
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
        x = inmf_obj_eval(Xs, Ws, H, Hvs, Nviews, Sp, lamb)
        #print("Best frob:", x.numpy())
        #print("Current frob", frobInit.numpy())
    
    
    ##-----------------------------------------------------------------------##
    ##             Convert to numpy, transpose and return                    ##
    ##-----------------------------------------------------------------------##
    Ws_num  = [Wi.numpy().T for Wi in Best_Ws]
    H_num   = Best_H.numpy().T
    Hvs_num = [Hvi.numpy().T for Hvi in Best_Hvs]


    frobNorm = [i.numpy() for i in frobNorm]
    W_eval_num  = [i.numpy().T for i in W_eval]

    return Ws_num, H_num, Hvs_num, iter_to_conv, frobNorm, W_eval_num

    



#%%

#iNMF_tensor_py(matrix_list = Xs_list,
#            rank              = 10,
#            n_initializations = 10,
#            iterations        = 1000,
#            Sp                = 0, 
#            stop_threshold    = 5)    



#%%


#
#iNMF(matrix_list = Xs_list,
#     rank           = 10,
#     iterations     = 1000,
#     L              = 1,
#     Sp             = 0, 
#     stop_threshold = 40)    


#%%
