# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#' Translate from Bratwurst to YAPSA
#'
#' Bratwurst (this package, wrapping NMF on GPU functionality) and YAPSA are
#' complementary. Both do signature analysis, but Bratwurst implements
#' unsupervised analyses (based on NMF), YAPSA runs supervised analyses (based
#' on LCD deconvolutions). Both have their own conventions and historically
#' evolved data structures. This function translates from the Bratwurst world
#' to the YAPSA world.
#'
#' @param nmf.exp Input data of class nmfExperiment
#'
#' @return
#'  A list of lists:
#'  outer list: different factorization ranks,
#'  inner list: different layers of information per factorization rank, e.g.
#'  exposures, norm_exposures, signatures or out_sig_ind_df
#'
#' @export
#'
#' @seealso \code{\link{normalizeW}}
#'
#' @examples
#'  NULL
#'
translateBratwurstToYAPSA <- function(nmf.exp, normalize = TRUE){
  if(normalize){
    norm.nmf.exp <- normalizeW(nmf.exp)
  } else {
    norm.nmf.exp <- nmf.exp
  }
  kList <- as.numeric(as.character(names(WMatrix(nmf.exp))))
  BratwurstListsList <-
    lapply(kList, function(current_k){
      temp_sig_df <- as.data.frame(WMatrix(norm.nmf.exp, k = current_k))
      rownames(temp_sig_df) <- rownames(norm.nmf.exp)
      colnames(temp_sig_df) <- paste0("S", seq(dim(temp_sig_df)[2]))
      temp_expo_df <- as.data.frame(HMatrix(norm.nmf.exp, k = current_k))
      colnames(temp_expo_df) <- colnames(norm.nmf.exp)
      rownames(temp_expo_df) <- paste0("S", seq(dim(temp_expo_df)[1]))
      temp_normExpo_df <- normalize_df_per_dim(temp_expo_df, in_dimension = 2)
      temp_sigInd_df <-
        data.frame(sig = colnames(temp_sig_df),
                   index = seq(dim(temp_sig_df)[2]),
                   colour = rainbow(dim(temp_sig_df)[2]),
                   process = rep("unknown", dim(temp_sig_df)[2]))
      return(list(exposures = temp_expo_df,
                  norm_exposures = temp_normExpo_df,
                  signatures = temp_sig_df,
                  out_sig_ind_df = temp_sigInd_df))
    })
  names(BratwurstListsList) <- kList
  return(BratwurstListsList)
}
