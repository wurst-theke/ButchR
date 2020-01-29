
#environment(myn_compute_OptKStats_NMF) <- asNamespace("Bratwurst")
#
#
# jnmf_exp <- run_joinNMF_tensor(norm_mat_list,
#                                ranks  = 2:5,
#                                n_initializations     = 3,
#                                iterations            = 10^4,
#                                convergence_threshold = 2,
#                                Sp = 0)
# jnmf_exp@OptKStats
