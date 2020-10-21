#------------------------------------------------------------------------------#
#                      NMF-GPU plot generation - FUNCTIONS                     #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                           Optimal K metrics                                  #
#------------------------------------------------------------------------------#
#' Plots optimal K metrics
#'
#' For every factorization rank the Frobenius error,
#' coefficient variation of Frobenius error,
#' sum Silhouette Width, mean Silhouette width,
#' cophenetic coefficient and mean Amari distance is shown.
#'
#' @param nmf_exp an object of class ButchR_NMF, ButchR_joinNMF, or
#' ButchR_integrativeNMF.
#' @param plot_vars a character vector with the ids of the metrics to display.
#' possible values are: FrobError, FrobError_min, FrobError_mean, FrobError_cv,
#' FrobError_sd, sumSilWidth, meanSilWidth, copheneticCoeff, and meanAmariDist.
#' Default value: c("FrobError", "FrobError_cv", "sumSilWidth", "meanSilWidth",
#' "copheneticCoeff", "meanAmariDist").
#'
#' @return a ggplot figure with the values for the selected factorization metrics
#' @import ggplot2 dplyr
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- run_NMF_tensor(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 2)
#' gg_plotKStats(nmf_exp)
#' }
gg_plotKStats <- function(nmf_exp,
                          plot_vars = c("FrobError", "FrobError_cv",
                                        "sumSilWidth", "meanSilWidth",
                                        "copheneticCoeff", "meanAmariDist")) {

  if (!class(nmf_exp) %in% c("ButchR_NMF", "ButchR_joinNMF", "ButchR_integrativeNMF")) {
    stop("\ngg_plotKStats is only supported for objects of class: \n",
         "ButchR_NMF, ButchR_joinNMF, or ButchR_integrativeNMF\n")
  }
  if (!is.character(plot_vars)) {
    stop("\nplot_vars: Expecting a character vector with IDs of the
         factorization metrics to visualize, e.g.:\n",
         "c('FrobError', 'FrobError_cv', 'sumSilWidth', 'meanSilWidth',
         'copheneticCoeff', 'meanAmariDist')\n")
  }
  if (!all(plot_vars %in% c("FrobError", "FrobError_min", "FrobError_mean",
                           "FrobError_cv", "FrobError_sd", "sumSilWidth",
                           "meanSilWidth", "copheneticCoeff", "meanAmariDist"))) {
    stop("\nPossible factorization metrics are: FrobError, FrobError_min,
    FrobError_mean, FrobError_cv, FrobError_sd, sumSilWidth, meanSilWidth,
    copheneticCoeff, and meanAmariDist.\n")
  }




  # Frobenius error for all initializations
  frob_df <- nmf_exp@FrobError %>%
    tidyr::pivot_longer(everything(), names_to = "k", values_to = "Stat") %>%
    dplyr::mutate(Metric = "FrobError") %>%
    dplyr::mutate(k = as.numeric(sub("^k", "", .data$k)))
  # Optimal factorization rank metrix coputed by compute_OptKStats_NMF
  metrics_df <- nmf_exp@OptKStats[,-1] %>%
    tidyr::pivot_longer(names_to = "Metric", values_to = "Stat", -.data$k)
  # Plot and highlight optK
  dplyr::bind_rows(frob_df, metrics_df) %>%
    dplyr::filter(.data$Metric %in% plot_vars) %>%
    dplyr::mutate(Metric = factor(.data$Metric, levels = unique(.data$Metric))) %>%
    ggplot(aes(x = .data$k, y = .data$Stat)) +
    geom_vline(xintercept = nmf_exp@OptK, color = "firebrick") +
    geom_point() +
    facet_wrap(.~Metric, scales = "free") +
    theme_bw()
}



#------------------------------------------------------------------------------#
#                                  Riverplot                                   #
#------------------------------------------------------------------------------#
#' @rdname generateRiverplot-methods
#' @aliases generateRiverplot,ANY,ANY-method
#' @import nnls riverplot grDevices
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- run_NMF_tensor(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 2)
#' plot(generateRiverplot(nmf_exp))
#' plot(generateRiverplot(nmf_exp, ranks = 2:5))
#' }
setMethod("generateRiverplot",
          "ButchR_NMF",
          function(nmf_exp, edges.cutoff = 0,
                   useH=FALSE, color=TRUE,
                   ranks=NULL) {

            #------------------------------------------------------------------#
            #                      Check if ranks exists                       #
            #------------------------------------------------------------------#
            if (!is.null(ranks)) {
              if (!is.numeric(ranks) | !length(ranks) >=2) {
                stop("Expecting a numeric vector with two or more values",
                     "\nPlease select from ranks = ",
                     paste0(nmf_exp@OptKStats$k, collapse = ","))
              }
              if (!all(ranks %in% nmf_exp@OptKStats$k)) {
                stop("No factorization present for all k in the range k =",
                     paste0(ranks, collapse = ","),
                     "\nPlease select from ranks = ",
                     paste0(nmf_exp@OptKStats$k, collapse = ","))
              }
            }

            #------------------------------------------------------------------#
            #                      Retrieve list of matrices                   #
            #------------------------------------------------------------------#
            if(useH) {
              W_list <- lapply(HMatrix(nmf_exp), t)
            } else {
              W_list <- WMatrix(nmf_exp)
            }


            #------------------------------------------------------------------#
            #                  Build data frame with node names                #
            #------------------------------------------------------------------#
            if (is.null(ranks)) {
              ranks <- nmf_exp@OptKStats$k
            } else {
              idx <- as.character(nmf_exp@OptKStats$rank_id[nmf_exp@OptKStats$k %in% ranks])
              W_list <- W_list[idx]
            }
            nodes <- do.call(rbind, lapply(ranks, function(i) {
              data.frame(ID = sapply(1:i, function(n) paste(i, "_S", n, sep = "")),
                         x = i)
            }))
            #------------------------------------------------------------------#
            #          Build data frame with edges values - NNLS               #
            #------------------------------------------------------------------#
            edges <- do.call(rbind, lapply(1:(length(ranks)-1), function(i){
              k     <- ranks[i]
              kplus <- ranks[i+1]
              df <- data.frame(
                N1 = rep(sapply(1:k, function(n) paste(k, "_S", n, sep = "")),
                         each = kplus),
                N2 = rep(sapply(1:kplus, function(n) paste(kplus, "_S", n, sep = "")),
                         time = k),
                Value = as.vector(t(nnls_sol(W_list[[i]], W_list[[i + 1]])))
              )
              df$ID <- paste(df$N1, df$N2)
              return(df)
            }))

            edges_dfList <- split(edges, f = edges$N2)
            edges_vecList <- lapply(edges_dfList, function(current_df){
              norm_vec <- current_df$Value / sum(current_df$Value)
              names(norm_vec) <- current_df$ID
              return(norm_vec)
            })
            edges_vec <- do.call(c, edges_vecList)
            names(edges_vec) <- unlist(lapply(edges_vecList, names))
            edges$rescaled <- edges_vec[as.character(edges$ID)]
            edges <- edges[edges$Value > edges.cutoff, ]
            nodes$y <- yNodeCoords(nodes, edges, rank_flag = TRUE,
                                   start_coords = c(-0.5, 0.5),
                                   edgeWeightColName = "rescaled",
                                   scale_fun = function(x){return(x^(1 / 3))})
            edges <- reorderEdges(nodes, edges)
            if (color){
              pca <- stats::prcomp(t(do.call(cbind, W_list)))
              pca <- apply(pca$x, 2, function(r) {
                r <- r - min(r)
                return(r / max(r))
              })
              col <- rgb(pca[, 1], pca[, 2], pca[, 3], 0.9)
              nodes$col <- col
            }
            ret <- makeRiver(nodes = nodes, edges = edges)
            return(ret)
          }
)


#------------------------------------------------------------------------------#
#  Riverplot - similarities between signatures at different ranks - FUNCTIONS  #
#------------------------------------------------------------------------------#
#' Find non-negative exposures of one matrix in another matrix
#'
#' @param B of A * X = B
#' @param A of A * X = B
#'
#' @return X of A * X = B
#'
#' @import nnls
#'
nnls_sol <- function(B, A) {
  X <- matrix(0, nrow = ncol(B), ncol = ncol(A))
  for(i in 1:ncol(B))
    X[i, ] <- stats::coef(nnls(A, B[, i]))
  X
}

#' Order the riverplot nodes to minimize crossings
#'
#' @param nodes riverplot nodes
#' @param edges riverplot nodes
#' @param rank_flag default FALSE
#' @param scale_factor default 1
#' @param start_coords default c(1,2)
#' @param edgeWeightColName "rescaled"
#' @param scale_fun function
#'
#' @return node_ypos
#'
yNodeCoords <- function(nodes, edges,
                        rank_flag = FALSE,
                        scale_factor = 1,
                        start_coords = c(1, 2),
                        edgeWeightColName = "rescaled",
                        scale_fun = function(x) { return(x) }){
  ranks <- unique(nodes[, 2])
  node_ypos <- rep(1, nrow(nodes))
  names(node_ypos) <- nodes[, 1]
  node_ypos[c(1, 2)] <- start_coords
  for(current_rank in ranks[-1]){
    #for(current_rank in c(3:10)){
    my_nodes <- as.character(nodes[which(nodes[, 2] == current_rank), 1])
    yCoords <- lapply(my_nodes, function(current_node){
      current_edges <- edges[which(as.character(edges[, 2]) == current_node), ]
      yCoord <- sum(scale_fun(current_edges[, edgeWeightColName]) *
                      node_ypos[as.character(current_edges[, 1])])
      return(yCoord)
    })
    my_factor <- scale_factor * (current_rank / (current_rank - 1))
    if(rank_flag) {
      node_ypos[my_nodes] <- rank(unlist(yCoords), ties.method = "first")
      if(any(start_coords < 0)){
        node_ypos[my_nodes] <-
          rank(unlist(yCoords), ties.method = "first") -
          0.5 * (current_rank + 1)
      }
    } else {
      node_ypos[my_nodes] <- my_factor * unlist(yCoords)
    }
  }
  return(node_ypos)
}

#' Order the riverplot edges to prevent crossings from a single node
#'
#' @param nodes riverplot nodes
#' @param edges riverplot nodes#'
#' @return edges
#'
reorderEdges <- function(nodes, edges){
  node_ypos <- nodes$y
  names(node_ypos) <- nodes$ID
  tempSum <- cumsum(unique(nodes[, 2]))
  offsetClasses <- c(0, tempSum[1:length(tempSum) - 1])
  offsets <- rep(offsetClasses, times = unique(nodes[, 2]))
  node_ranks <- node_ypos + offsets
  edgesOrder <- order(node_ranks[as.character(edges$N1)],
                      node_ranks[as.character(edges$N2)])
  return(edges[edgesOrder, ])
}



#' Relabel and recolour a riverplot
#'
#' @param in_riverplot riverplot object
#' @param in_list list of new colors and labels to use
#'
#' @return relabeled riverplots
#' @export
#'
relabelRiverplot <- function(in_riverplot, in_list){
  in_riverplot$nodes$labels <-
    in_list$name_vector[as.character(in_riverplot$nodes$ID)]
  tempVec <-lapply(
    strsplit(in_list$sig_names[as.character(in_riverplot$nodes$ID)], split = "_"),
    function(x) utils::head(x, n = 1))
  tempVec <- unlist(tempVec)
  in_riverplot$nodes$col <-
    as.character(in_list$col_vector[as.character(tempVec)])
  temp_list <-
    lapply(seq_len(dim(in_riverplot$nodes)[1]), function(current_nodeInd){
      in_riverplot$styles[[current_nodeInd]]$col <-
        as.character(in_riverplot$nodes$col[current_nodeInd])
      return(in_riverplot$styles[[current_nodeInd]])
    })
  names(temp_list) <- names(in_riverplot$styles)
  in_riverplot$styles <- temp_list
  return(in_riverplot)
}

#------------------------------------------------------------------------------#
#                            Recovery plots for matrix                         #
#------------------------------------------------------------------------------#
#' @rdname recovery_plot-methods
#' @aliases recovery_plot,ANY,ANY-method
#' @import ggplot2 dplyr cowplot
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' data(leukemia)
#' leukemia_nmf_exp <- run_NMF_tensor(X = leukemia$matrix,
#'                                    ranks = 2:4,
#'                                    method = "NMF",
#'                                    n_initializations = 10,
#'                                    extract_features = TRUE)
#' recovery_plot(HMatrix(leukemia_nmf_exp, k = 4),
#'               leukemia$annotation$ALL_AML)
#' }
setMethod("recovery_plot",
          "matrix",
          function(x, annot){
            h <- x
            # Add sig IDs if missing
            if (is.null(rownames(h))) {
              rownames(h) <- paste0('Sig ',1:nrow(h))
            }
            # check i annot lenght == h

            if (is.factor(annot) | is.character((annot))) {
              # annot_list <- list(main_annot = as.factor(annot))
              annot_factor <- as.factor(annot)
              if (is.null(names(annot))) {
                names(annot_factor) <- colnames(h)
              }
            } else {
              stop("Not a valid annotation input")
            }
            n_samples  <- ncol(h)
            #which.a = annotID
            #annot.factor <- annot[,annotID]


            ## -------------------------------------------------------------------##
            ##                        Find ranks                                  ##
            ##--------------------------------------------------------------------##
            # cycle annot levels
            lIds <- stats::setNames(levels(annot_factor), levels(annot_factor))
            ALL_RNKS <-  lapply(lIds,function(l) {
              # cycle h matrix rows and find ranks
              lapply(stats::setNames(1:nrow(h), rownames(h)),function(i) {
                exp   <-  sort(h[i,],decreasing=TRUE) # sorted exposure
                i_rnk <-  match(names(annot_factor)[annot_factor==l], names(exp))
                sort(i_rnk[!is.na(i_rnk)]) # keep steps/ranks
              })
              #print(RNKS)
              #return(RNKS)
            })
            # ALL_RNKS

            ## -------------------------------------------------------------------##
            ##                  Find AUC and P-value                              ##
            ##--------------------------------------------------------------------##
            AUC_singleannot <- lapply(ALL_RNKS,function(r) {
              # AUC random set
              AUC_RAND <- do.call("rbind",lapply(r, function(x) {
                l = lapply(1:500,function(i) {
                  sample(1:n_samples, length(x))
                })
                aux = auc(l, max = n_samples)
                return(c(mean(aux), stats::sd(aux)))
              }))

              # AUC
              #AUC <-  lapply(ALL_RNKS, auc, max = n_samples)
              AUC <-  auc(r, max = n_samples)
              #print(AUC)

              # Find P - value
              AUC_df <- data.frame(AUC_RAND, AUC)
              colnames(AUC_df) = c('mean','sd','val')
              AUC_df <- AUC_df %>%
                tibble::rownames_to_column("SignatureID") %>%
                mutate(z = (.data$val - .data$mean)/.data$sd) %>%
                mutate(p = ifelse(.data$z>0,
                                  stats::pnorm(.data$z, lower.tail=FALSE),
                                  stats::pnorm(.data$z)))

              #Return randon and AUC - P-val
              return(AUC_df)
            })
            # AUC_allannot <- bind_rows(AUC_singleannot, .id = "Annotation_level") %>%
            #   mutate(Annotation = "main_annot")
            AUC_allannot <- bind_rows(AUC_singleannot, .id = "Annotation_level")

            # Add min and max to rank, for step plot
            # cycle all annots
            # cycle annot levels
            ALL_RNKS <- lapply(ALL_RNKS, function(x){
              # cycle h matrix rows and find ranks
              lapply(x, function(xi) c(0, xi, n_samples))
            })



            ALL_RNKS_df <- bind_rows(ALL_RNKS, .id = "Annotation_level") %>%
              tidyr::pivot_longer(-c("Annotation_level"),  names_to = "SignatureID", values_to = "Rank") %>%
              left_join(AUC_allannot, by = c("Annotation_level", "SignatureID"))


            #return(ALL_RNKS_df)

            gg_recov <- ALL_RNKS_df %>%
              group_by(.data$Annotation_level, .data$SignatureID ) %>%
              mutate(Frequency = c(seq(0, 1, length.out = n()-1), 1)) %>% # all y axis step
              mutate(issignif = .data$p < 0.05) %>%

              ggplot(aes(x = .data$Rank, y = .data$Frequency, color = .data$SignatureID,
                         linetype = .data$issignif, size = .data$issignif)) +

              # geom_step(data = function(x){x %>% filter(!issignif)}, size  = 0.5) +
              # geom_step(data = function(x){x %>% filter(issignif)}, size  = 1.5) +
              geom_step() +

              geom_abline(intercept = 0, slope = 1/n_samples) +
              facet_wrap(.~Annotation_level) +
              # chance line style
              scale_linetype_manual(name = c("Significant p-val<0.05"),
                                    values = c("TRUE" = 1, "FALSE" = 2)) +
              scale_size_manual(name = c("Significant p-val<0.05"),
                                values = c("TRUE" = 1, "FALSE" = 0.5)) +
              #theme_bw() +
              theme_cowplot() +
              panel_border(color = "grey40", size = 1, linetype = 1,
                           remove = FALSE)


            #return(ALL_RNKS_df)
            return(gg_recov)
          }
)


#' Helper function to estimate AUC in the recovery plots
#'
#' @param rnk.list list of ranks for a particular annotation
#' @param max maximum rank
#' @return AUC
#'
auc <- function(rnk.list, max = NULL) {
  aux <-  sapply(rnk.list,function(rnk) {
    if (is.null(max)) {
      max <-  max(rnk)
    }
    rnk <-  sort(rnk)
    X <-  0
    i <-  1
    ngenes <-  length(rnk)
    while ((rnk[i] <= max) && (i <= length(rnk))) {
      X <-  X + max -rnk[i]
      i <-  i+1
    }

    rauc <-  X/(i-1)/max
    rauc
  })
  return(aux)
}