#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#from __future__ import division
##import numpy as np
#import os
import tensorflow as tf
##import argparse
#os.environ['TF_CPP_MIN_LOG_LEVEL'] = "1"  # or any {'0', '1', '2'}


"""
Created on Thu Jan  9 11:36:20 2020

@author: Andres Quintero
"""


##---------------------------------------------------------------------------##
##                Define Function to run NMF in tensorflow 2                 ##
##---------------------------------------------------------------------------##
def NMF_tensor_py(matrix, rank, n_initializations, iterations, stop_threshold=40, 
                  n_neighbors = 4, alpha = 0.1, lamb = 10, graph = None, **kwargs):
    
    
    # NMF in tensorflow
    n = matrix.shape[0]
    m = matrix.shape[1]     
    X = tf.constant(matrix, name = ("X"), dtype=tf.float32)
    
    # Initialization Metrics 
    frobNorm = []
    iter_to_conv = []
    W_eval = []
    
    ##-----------------------------------------------------------------------##
    ##                   Compute nearest neighbors graph                     ##
    ##-----------------------------------------------------------------------##
    if graph is None:
        r = tf.reduce_sum(tf.transpose(X)*tf.transpose(X), 1)
        # turn r into column vector
        r = tf.reshape(r, [-1, 1])
        D = r - 2*tf.matmul(X, X, transpose_a=True) + tf.transpose(r)
        #print(D)
    else:
        #print(graph)
        D = graph
    # Find top n_neighbors neighbors
    b = tf.nn.top_k(D, (m-n_neighbors-1))
    kth = tf.reduce_min(b.values, axis=1)
    #kth = tf.broadcast_to(kth , D.shape)
    
    # Binary matrix with nearest neighbors across colummns
    G = tf.less(D, kth)
    # Binary matrix with nearest neighbors across rows and colummns
    G = tf.math.logical_or(G, tf.transpose(G))
    G = tf.cast(G, tf.float32)
    #print(G)
    
    # diagonal matrix whose entries are column (or row, since W is symmetric) sums of W
    D = tf.reduce_sum(G, axis=0)
    D = tf.linalg.tensor_diag(D)
    
    L = D - G 
    #    print(D.numpy()[0:5, 0:5] )
    #    print(G.numpy()[0:5, 0:5] )
    #    print(L.numpy()[0:5, 0:5] )
    #print(G)
    #print('nearest neighbors graph G completed')

    ##-----------------------------------------------------------------------##
    ##                   define objective function                           ##
    ##-----------------------------------------------------------------------##
    def nmf_obj_eval(X, W, H, L, alpha, lamb):
        # Frobenius term
        frob_c   = tf.math.square(tf.linalg.norm(X - tf.matmul(W, H), 
                                                 ord = 'fro', axis=[-2,-1]))
        #frob_c   = tf.linalg.norm(X - tf.matmul(W, H), ord = 'fro', axis=[-2,-1]) /  tf.linalg.norm(X, ord = 'fro', axis=[-2,-1])
        # Graph term
        graph_c  = tf.matmul(tf.matmul(H, L), H, transpose_b=True)
        graph_c  = lamb * tf.linalg.trace(graph_c)
        # Sparsity constrain
        sparse_c = alpha * tf.reduce_sum(tf.linalg.norm(H, ord = 1, axis = 1))
        #print(frob_c)
        #print(graph_c)
        #print(sparse_c)
        return frob_c + graph_c + sparse_c
        #return frob_c 
        

    #print(top2)
    #x
    ##-----------------------------------------------------------------------##
    ##                              N inits                                  ##
    ##-----------------------------------------------------------------------##
    # cycle through n initializations and choose best factorization
    for init_n in range(n_initializations):
        
        
        ##-------------------------------------------------------------------##
        ##                  Initialize W and H matrices                      ##
        ##-------------------------------------------------------------------##
        initializer = tf.random_uniform_initializer(minval=0, maxval=1)
        
        H = tf.Variable(initializer(shape=[rank, m]), name="H")
        W = tf.Variable(initializer(shape=[n, rank]), name="W")
    
    
        ##-------------------------------------------------------------------##
        ##        Save initial max exposures in H matrices                   ##
        ##-------------------------------------------------------------------##
        #Hts = [tf.add(H, Hv) for Hv in Hvs]
        #oldExposures = tf.concat([tf.math.argmax(Ht, axis=1) for Ht in Hts], 0)
        oldExposures = tf.math.argmax(H, axis=0)
        const = 0            
    
        ##-------------------------------------------------------------------##
        ##                          Run NMF                                  ##
        ##-------------------------------------------------------------------##
        for inner in range(iterations):
            ##---------------------------------------------------------------##
            ##                          Update H                             ##
            ##---------------------------------------------------------------##
            XtW   = tf.matmul(X, W, transpose_a=True)
            GHt   = lamb * tf.matmul(G, H, transpose_b=True)
            H_num = (2 * (XtW + GHt)) - alpha
            
            WtW   = tf.matmul(W, W, transpose_a=True)
            HWtW  = tf.matmul(H, WtW, transpose_a=True)
            DHt   = lamb * tf.matmul(D, H, transpose_b=True)
            H_den = 2 * (HWtW + DHt)
            
            newH = tf.math.multiply(H, tf.transpose((tf.math.divide(H_num, H_den))))
            newH = tf.where(tf.math.is_nan( newH ), tf.zeros_like( newH ), newH )
            H.assign(newH)            
    
            ##---------------------------------------------------------------##
            ##                          Update W                             ##
            ##---------------------------------------------------------------##
            XHt   = tf.matmul(X, H, transpose_b=True)
            WH    = tf.matmul(W, H)
            WHHt  = tf.matmul(WH, H, transpose_b=True)
      
            newW  = tf.math.multiply(W, tf.math.divide(XHt, WHHt))
            newW  = tf.where(tf.math.is_nan( newW ), tf.zeros_like( newW ), newW )
            #print(newW.numpy()[1:5,1:5] )            
            
            W_cor = tf.reduce_sum(newW, axis=0)
            newW  = tf.math.divide(newW, W_cor)
            W.assign(newW)

            
            ##---------------------------------------------------------------##
            ##                    Evaluate Convergence                       ##
            ##---------------------------------------------------------------##
            #print(nmf_obj_eval(X, W, H, L, alpha, lamb).numpy())
            newExposures = tf.math.argmax(H, axis=0)
            if tf.reduce_all(tf.math.equal(oldExposures, newExposures)).__invert__():
                oldExposures = newExposures
                const = 0
            else:
                const += 1
                #print(f'new eval diff, {const}')
                if const == stop_threshold:
                    #print(f"NMF converged after {i} iterations")
                    break
        
        ##-------------------------------------------------------------------##
        ##                     Initialize frob error                         ##
        ##-------------------------------------------------------------------##
        if init_n == 0 :    
            Best_frob = nmf_obj_eval(X, W, H, L, alpha, lamb)
            #Best_frob = tf.linalg.norm(X - tf.matmul(W, H)) / tf.linalg.norm(X)
            Best_H    = H
            Best_W    = W
        
        
        ##-------------------------------------------------------------------##
        ##         Evaluate if best factorization initialization             ##
        ##-------------------------------------------------------------------##
        frobInit = nmf_obj_eval(X, W, H, L, alpha, lamb)
        #frobInit = tf.linalg.norm(X - tf.matmul(W, H)) / tf.linalg.norm(X)
        # Append to list of initialization metrics
        frobNorm.append(frobInit)
        iter_to_conv.append(inner+1)
        W_eval.append(W)
        
        if frobInit < Best_frob :
            #print('is less')
            Best_frob = frobInit
            Best_H    = H
            Best_W    = W
        #x = nmf_obj_eval(X, W, H, L, alpha, lamb)
        #print("Best frob:", Best_frob.numpy())
        #print("Current frob", frobInit.numpy())
        #fb = tf.reduce_sum(fb, 0)
    
    ##-----------------------------------------------------------------------##
    ##             Convert to numpy, transpose and return                    ##
    ##-----------------------------------------------------------------------##
    W_num  = Best_W.numpy()
    H_num  = Best_H.numpy()

    frobNorm    = [i.numpy() for i in frobNorm]
    W_eval_num  = [i.numpy() for i in W_eval]
    
    return W_num, H_num, iter_to_conv, frobNorm, W_eval_num
    #return W_num, H_num, iter_to_conv, frobNorm


#%%


#NMF_tensor_py(matrix = Xs_list[0],
#            rank              = 10,
#            n_initializations = 1,
#            iterations        = 10,   
#            stop_threshold    = 40,
#            n_neighbors = 4,
#            more =1)    


#%%

#Xs_list[0].T

#%%
#import numpy as np
#filenames = ["/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/rna_norm_mat.csv", "/Users/andresq/phd/main_project/test_cases/buenrostro_AML/data/multiview/norm_matrices/atac_norm_mat.csv"]
#Xs_list = [np.loadtxt(path) for path in filenames]
#Xs_list 
#
#
#filenames = "/Users/andresq/phd/main_project/test_cases/scCAT_HumanEmbryo/results/integrative/regulatoryRelationships_1000Kb/HumanEmbryo_regulatoryRelationships_1000Kb.csv"
#X = np.loadtxt(filenames)
#X
#
