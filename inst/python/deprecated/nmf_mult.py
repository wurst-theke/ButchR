#-*- coding: utf-8 -*-
# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

import numpy as np
import cudamat as cm
import argparse

def NMF(X,r,iterations,H=None,W=None, niter_test_conv=10, stop_threshold=40):
    # nmf in cudamat with preallocated matrices

    rng = np.random
    n = np.size(X,0)
    m = np.size(X,1)
    if(H is None):
        H = rng.random((r,m)).astype(np.float32)
    if(W is None):
        W = rng.random((n,r)).astype(np.float32)

    H_gpu = cm.CUDAMatrix(H)
    W_gpu = cm.CUDAMatrix(W)
    X_gpu = cm.CUDAMatrix(X)
    WTW_gpu = cm.empty((r, r))
    WTWH_gpu = cm.empty(H.shape)
    WTX_gpu = cm.empty(H.shape)
    XHT_gpu = cm.empty(W.shape)
    WH_gpu = cm.empty(X.shape)
    WHHT_gpu = cm.empty(W.shape)

    #frobNorm = float('inf')
    const = 0
    oldExposures = np.argmax(H, axis=0)

    for i in range(0,iterations):
        # update H
        cm.dot(W_gpu.T, X_gpu, target=WTX_gpu)
        cm.dot(W_gpu.T, W_gpu, target=WTW_gpu)
        cm.dot(WTW_gpu, H_gpu, target=WTWH_gpu)
        H_gpu.mult(WTX_gpu).divide(WTWH_gpu)

        # update W
        cm.dot(X_gpu, H_gpu.T, target=XHT_gpu)
        cm.dot(W_gpu, H_gpu, target=WH_gpu)
        cm.dot(WH_gpu, H_gpu.T, target=WHHT_gpu)
        W_gpu.mult(XHT_gpu).divide(WHHT_gpu)

        if i % niter_test_conv == 0:
            newExpo = np.argmax(H_gpu.asarray(), axis=0)
            if (oldExposures != newExpo).any():
                oldExposures = newExpo
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    print "NMF converged after %i iterations" % i
                    break

    return W_gpu.asarray(),H_gpu.asarray()

def NMFsparse(X, r, sparseH=0, sparseW=0, iterations=1000, H=None, W=None, niter_test_conv=10, stop_threshold=40):
    rng = np.random
    n = np.size(X,0)
    m = np.size(X,1)
    if(H is None):
        H = rng.random((r,m)).astype(np.float32)
    if(W is None):
        W = rng.random((n,r)).astype(np.float32)

    H_gpu = cm.CUDAMatrix(H)
    W_gpu = cm.CUDAMatrix(W)
    X_gpu = cm.CUDAMatrix(X)
    XHT_gpu = cm.empty(W.shape)
    WH_gpu = cm.empty(X.shape)
    WHHT_gpu = cm.empty(W.shape)
    Wrowsum_gpu = cm.empty([1, r])
    WTWH_gpu = cm.empty(H.shape)
    WTX_gpu = cm.empty(H.shape)

    ### helpers as sum is slower than cm.dot
    nones_gpu = cm.CUDAMatrix(np.ones([1,n]))

    const = 0
    oldExposures = np.argmax(H, axis=0)

    for i in range(0,iterations):

        # update W
        cm.dot(W_gpu, H_gpu, target=WH_gpu)
        cm.dot(X_gpu, H_gpu.T, target=XHT_gpu)
        cm.dot(WH_gpu, H_gpu.T, target=WHHT_gpu)
        WHHT_gpu.add(sparseW)
        W_gpu.mult(XHT_gpu).divide(WHHT_gpu)

        # normalize W
        cm.dot(nones_gpu, W_gpu, target=Wrowsum_gpu) # slower correct version: W_gpu.sum(0, target=rowsum_gpu)
        W_gpu.div_by_row(Wrowsum_gpu)

        # update H
        cm.dot(W_gpu.T, X_gpu, target=WTX_gpu)
        cm.dot(W_gpu.T, WH_gpu, target=WTWH_gpu)
        WTWH_gpu.add(sparseH)
        H_gpu.mult(WTX_gpu).divide(WTWH_gpu)

        if i % niter_test_conv == 0:
            newExpo = np.argmax(H_gpu.asarray(), axis=0)
            if (oldExposures != newExpo).any():
                oldExposures = newExpo
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    print "NMF converged after %i iterations" % i
                    break

    return W_gpu.asarray(),H_gpu.asarray()

def NMFaffine(X, r, sparseH=0, sparseW=0, iterations=1000, H=None, W=None, niter_test_conv=10, stop_threshold=40):
    rng = np.random
    n = np.size(X,0)
    m = np.size(X,1)
    if(H is None):
        H = rng.random((r,m)).astype(np.float32)
    if(W is None):
        W = rng.random((n,r)).astype(np.float32)

    H_gpu = cm.CUDAMatrix(H)
    W_gpu = cm.CUDAMatrix(W)
    X_gpu = cm.CUDAMatrix(X)
    W0_gpu = cm.CUDAMatrix(X.mean(1)[:,None])
    XHT_gpu = cm.empty(W.shape)
    WH_gpu = cm.empty(X.shape)
    WHHT_gpu = cm.empty(W.shape)
    rowsum_gpu = cm.empty([1, r])
    Xcolsum_gpu = cm.empty([n, 1])
    WHcolsum_gpu = cm.empty([n, 1])
    WTWH_gpu = cm.empty(H.shape)
    WTX_gpu = cm.empty(H.shape)

    ### helpers as sum is slower than cm.dot
    nones_gpu = cm.CUDAMatrix(np.ones([1,n]))
    mones_gpu = cm.CUDAMatrix(np.ones([m,1]))

    const = 0
    oldExposures = np.argmax(H, axis=0)

    X_gpu.sum(1, target=Xcolsum_gpu)
    for i in range(0,iterations):
        # update W
        cm.dot(W_gpu, H_gpu, target=WH_gpu)
        WH_gpu.add_col_vec(W0_gpu)
        cm.dot(X_gpu, H_gpu.T, target=XHT_gpu)
        cm.dot(WH_gpu, H_gpu.T, target=WHHT_gpu)
        WHHT_gpu.add(sparseW)
        W_gpu.mult(XHT_gpu).divide(WHHT_gpu)

        # normalize W
        cm.dot(nones_gpu, W_gpu, target=rowsum_gpu) # slower correct version: W_gpu.sum(0, target=rowsum_gpu)
        W_gpu.div_by_row(rowsum_gpu)

        # update H
        cm.dot(W_gpu.T, X_gpu, target=WTX_gpu)
        cm.dot(W_gpu.T, WH_gpu, target=WTWH_gpu)
        WTWH_gpu.add(sparseH)
        H_gpu.mult(WTX_gpu).divide(WTWH_gpu)

        # update W0
        cm.dot(WH_gpu, mones_gpu, target=WHcolsum_gpu) # slower correct version: WH_gpu.sum(1, target=WHcolsum_gpu)
        W0_gpu.mult(Xcolsum_gpu).divide(WHcolsum_gpu)

        if i % niter_test_conv == 0:
            newExpo = np.argmax(H_gpu.asarray(), axis=0)
            if (oldExposures != newExpo).any():
                oldExposures = newExpo
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    print "NMF converged after %i iterations" % i
                    break

    return W_gpu.asarray(),H_gpu.asarray(),W0_gpu.asarray()

def NMFsemi(X, r, iterations=1000, G=None, niter_test_conv=10, stop_threshold=40):

    n = np.size(X,0)
    m = np.size(X,1)

    if G is None:
        G = np.random.random((m,r)).astype(np.float32)
    elif G.strides[1] > G.strides[0]:
        # this is a check wether the data array is correct and not a transpose
        # of some other array that is hard to process (due to strides problem)
        G = G.copy()

    # allocate the matrices on the GPU
    G_gpu = cm.CUDAMatrix(G)
    F_gpu = cm.empty((n,r))
    X_gpu = cm.CUDAMatrix(X)
    GTG_gpu = cm.empty((r,r))
    GTGpinv_gpu = cm.empty((r,r))
    XG_gpu = cm.empty((n,r))
    XTF_gpu = cm.empty((m,r))
    FTF_gpu = cm.empty((r,r))
    XTFgreater_gpu = cm.empty((m,r))
    FTFgreater_gpu = cm.empty((r,r))
    XTFpos_gpu = cm.empty((m,r))
    XTFneg_gpu = cm.empty((m,r))
    FTFpos_gpu = cm.empty((r,r))
    FTFneg_gpu = cm.empty((r,r))
    GFTFneg_gpu = cm.empty((m,r))
    GFTFpos_gpu = cm.empty((m,r))

    const = 0
    oldExposures = np.argmax(G, axis=0)

    for i in range(iterations):
        # F = XG(G.T G)^-1
        cm.dot(G_gpu.T, G_gpu, target=GTG_gpu)
        try:
            GTGpinv_gpu = cm.CUDAMatrix(np.linalg.inv(GTG_gpu.asarray()))
        except LinAlgError:
            GTGpinv_gpu = cm.CUDAMatrix(np.linalg.pinv(GTG_gpu.asarray()))
        cm.dot(X_gpu, G_gpu, target=XG_gpu)
        cm.dot(XG_gpu, GTGpinv_gpu, target=F_gpu)

        # preparation and calculation of the matrix separations
        cm.dot(X_gpu.T, F_gpu, target=XTF_gpu)
        cm.dot(F_gpu.T, F_gpu, target=FTF_gpu)

        XTF_gpu.greater_than(0, target=XTFgreater_gpu)
        XTF_gpu.mult(XTFgreater_gpu, target=XTFpos_gpu)
        XTFpos_gpu.subtract(XTF_gpu, target=XTFneg_gpu)

        FTF_gpu.greater_than(0, target=FTFgreater_gpu)
        FTF_gpu.mult(FTFgreater_gpu, target=FTFpos_gpu)
        FTFpos_gpu.subtract(FTF_gpu, target=FTFneg_gpu)

        # compute the G update
        cm.dot(G_gpu, FTFpos_gpu, target=GFTFpos_gpu)
        cm.dot(G_gpu, FTFneg_gpu, target=GFTFneg_gpu)

        XTFpos_gpu.add(GFTFneg_gpu)
        XTFneg_gpu.add(GFTFpos_gpu)
        XTFpos_gpu.divide(XTFneg_gpu)
        cm.sqrt(XTFpos_gpu)

        G_gpu.mult(XTFpos_gpu)

        if i % niter_test_conv == 0:
            newExpo = np.argmax(G_gpu.asarray(), axis=0)
            if (oldExposures != newExpo).any():
                oldExposures = newExpo
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    print "NMF converged after %i iterations" % i
                    break

    return F_gpu.asarray(), G_gpu.asarray().T

def NMFconvex(X, r, iterations=1000, G=None, niter_test_conv=10, stop_threshold=40):

    n = np.size(X,0)
    m = np.size(X,1)

    if G is None: # implemnt k means initialization
        G = np.random.random((m,r)).astype(np.float32)
        W = np.random.random((m,r)).astype(np.float32)
    else:
        G += 0.2
    Wi = np.dot(G, np.linalg.inv(np.dot(G.T,G)))
    Wipos = (np.abs(Wi) + Wi) / 2
    W = Wipos + 0.2 * np.sum(np.abs(Wi)) / Wi.size

    G_gpu = cm.CUDAMatrix(G)
    W_gpu = cm.CUDAMatrix(W)
    X_gpu = cm.CUDAMatrix(X)

    XTX_gpu= cm.dot(X_gpu.T, X_gpu)
    XTXpos_gpu = cm.empty((m,m))
    XTX_gpu.greater_than(0, target=XTXpos_gpu)
    XTXpos_gpu.mult(XTX_gpu)
    XTXneg_gpu = cm.empty((m,m))
    XTXpos_gpu.subtract(XTX_gpu, target=XTXneg_gpu)

    XTXnegW_gpu = cm.empty((m,r))
    XTXposW_gpu = cm.empty((m,r))
    GWT_gpu = cm.empty((m,m))
    update1_gpu = cm.empty((m,r))
    update2_gpu = cm.empty((m,r))

    GTG_gpu = cm.empty((r,r))
    XTXnegG_gpu = cm.empty((m,r))
    XTXposG_gpu = cm.empty((m,r))

    const = 0
    oldExposures = np.argmax(G, axis=1)

    for i in range(0, iterations):

        cm.dot(XTXneg_gpu, W_gpu, target=XTXnegW_gpu)
        cm.dot(XTXpos_gpu, W_gpu, target=XTXposW_gpu)

        # Update G
        cm.dot(G_gpu, W_gpu.T, target=GWT_gpu)
        # G *= np.sqrt((XTXposW + np.dot(GWT, XTXnegW))/(XTXnegW+np.dot(GWT, XTXposW)))
        cm.dot(GWT_gpu, XTXnegW_gpu, target=update1_gpu)
        cm.dot(GWT_gpu, XTXposW_gpu, target=update2_gpu)
        update1_gpu.add(XTXposW_gpu)
        update2_gpu.add(XTXnegW_gpu)
        update1_gpu.divide(update2_gpu)
        cm.sqrt(update1_gpu)
        G_gpu.mult(update1_gpu)

        # Update W
        cm.dot(G_gpu.T, G_gpu, target=GTG_gpu)
        #W *= np.sqrt((np.dot(XTXpos, G) + np.dot(XTXnegW, GTG))/(np.dot(XTXneg, G) + np.dot(XTXposW, GTG)))
        cm.dot(XTXpos_gpu, G_gpu, target=XTXposG_gpu)
        cm.dot(XTXneg_gpu, G_gpu, target=XTXnegG_gpu)
        cm.dot(XTXnegW_gpu, GTG_gpu, target=update1_gpu)
        cm.dot(XTXposW_gpu, GTG_gpu, target=update2_gpu)
        update1_gpu.add(XTXposG_gpu)
        update2_gpu.add(XTXnegG_gpu)
        update1_gpu.divide(update2_gpu)
        cm.sqrt(update1_gpu)
        W_gpu.mult(update1_gpu)

        if i % niter_test_conv == 0:
            newExpo = np.argmax(G_gpu.asarray(), axis=1)
            if (oldExposures != newExpo).any():
                oldExposures = newExpo
                const = 0
            else:
                const += 1
                if const == stop_threshold:
                    print "NMF converged after %i iterations" % i
                    break

    return cm.dot(X_gpu, W_gpu).asarray(), G_gpu.asarray().T

if __name__ == "__main__":

    import time
    t0 = time.time()

    print "NMF-CUDA from a python script via cudamat"

    parser = argparse.ArgumentParser(description='A python script for running NMF on Cuda',
                                     epilog='Dependencies: numpy, cudamat',
                                     usage='prog <filename> [options]')
    parser.add_argument("filename", default=None, help="The file of the input matrix X")
    parser.add_argument("-k", "-K", dest="rank", type=int, default=2, help="Factorization Rank (default: 2)")
    parser.add_argument("-i", "-I", dest="iter", type=int, default=2000, help="Maximum number of iterations")
    parser.add_argument("-j", "-J", dest="niter_test_conv", type=int, default=10000, help="Perform a convergence test each <niter_test_conv> iterations (default: 10000). If this value is greater than <nIters> (see '-i' option), no test is performed")
    parser.add_argument("-t", "-T", dest="stop_threshold", type=int, default=40, help="When matrix H has not changed on the last <stop_threshold> times that the convergence test has been performed, it is considered that the algorithm has converged to a solution and stops it.")
    parser.add_argument("-s", "-S", dest="type", default=None, type=str, help="Type of NMF, P for sparse, S for semi, C for convex and A for affine, else normal NMF")
    parser.add_argument("-ho", "-HO", dest="sparseH", default=0, type=float, help="Sparseness parameter of H matrix for sparse and affine NMF")
    parser.add_argument("-wo", "-WO", dest="sparseW", default=0, type=float, help="Sparseness parameter of W matrix for sparse and affine NMF")
    parser.add_argument("-g", "-G", dest="gpuID", default=0, type=int, help="ID of the GPU, if multiple GPUs are available")
    parser.add_argument("-e", "-E", dest="encoding", default="txt", help="File encoding of input and output matrices")
    parser.add_argument("-sets", dest="sets", type=bool, default=False, help="Binary whether seed from random state should be set")
    parser.add_argument("-sv", dest="seedv", type=int, default=0, help="Binary whether seed from random state should be set")

    args = parser.parse_args()

    if args.encoding == "npy":
        X = np.load(args.filename)
    else:
        X = np.loadtxt(args.filename)

    X = X.astype(np.float32)

    cm.cuda_set_device(args.gpuID)
    cm.cublas_init()

    if args.sets:
        np.random.seed(seed=args.seedv)

    if args.type == "P":
        print "Running sparse NMF with sparseness constraints %f for H and %f for W" % (args.sparseH, args.sparseW)
        W, H = NMFsparse(X, args.rank, args.sparseH, args.sparseH, args.iter, niter_test_conv=args.niter_test_conv, stop_threshold=args.stop_threshold)
        frobNorm = np.linalg.norm(X-np.dot(W,H)) / np.linalg.norm(X)
    elif args.type == "A":
        print "Running affine NMF"
        W, H, W0 = NMFaffine(X, args.rank, args.sparseH, args.sparseH, args.iter, niter_test_conv=args.niter_test_conv, stop_threshold=args.stop_threshold)
        np.savetxt(args.filename + "_W0.txt", W0)
        frobNorm = np.linalg.norm(X-np.dot(W,H)-W0) / np.linalg.norm(X)
    elif args.type == "S":
        print "Running semi NMF"
        W, H = NMFsemi(X, args.rank, args.iter, niter_test_conv=args.niter_test_conv, stop_threshold=args.stop_threshold)
        frobNorm = np.linalg.norm(X-np.dot(W,H)) / np.linalg.norm(X)
    elif args.type == "C":
        print "Running convex NMF"
        W, H = NMFconvex(X, args.rank, args.iter, niter_test_conv=args.niter_test_conv, stop_threshold=args.stop_threshold)
        frobNorm = np.linalg.norm(X-np.dot(W,H)) / np.linalg.norm(X)
    else:
        print "Running default NMF"
        W, H = NMF(X, args.rank, args.iter, niter_test_conv=args.niter_test_conv, stop_threshold=args.stop_threshold)
        frobNorm = np.linalg.norm(X-np.dot(W,H)) / np.linalg.norm(X)

    if args.encoding == "npy":
        # the RcppCNPy library needs arrays of float64 and with correct strides
        np.save(args.filename + "_H.npy", H.astype(np.float64).copy())
        np.save(args.filename + "_W.npy", W.astype(np.float64).copy())
    else:
        np.savetxt(args.filename + "_H.txt", H)
        np.savetxt(args.filename + "_W.txt", W)

    print "Distance: " + str(frobNorm)
    t1 = time.time()
    print "Time take by NMF-cuda: ", t1-t0

