#------------------------------------------------------------------------------#
#                      NMF-GPU plot generation - FUNCTIONS                     #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                           Optimal K metrics                                  #
#------------------------------------------------------------------------------#
#' Plots optimal K metrics
#'
#' For every factorization rank the Frobinius error,
#' coefficient variation of Frobinius error,
#' sum Silhouette Width, mean Silhouette width,
#' cophenetic coefficient and mean Amari distance is shown
#'
#' @param nmf_exp A nmfExperiment_lite object
#' @param plot_vars character - ids of the metrics to display
#'
#' @return a ggplot figure with the values for the six factorization metrics
#' @import ggplot2 dplyr
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' data("leukemia")
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 2)
#' gg_plotKStats(nmf_exp)
#' }
gg_plotKStats <- function(nmf_exp,
                          plot_vars = c("FrobError", "cv", "sumSilWidth",
                                        "meanSilWidth", "copheneticCoeff",
                                        "meanAmariDist")) {
  frob_df <- nmf_exp@FrobError %>%
    tidyr::pivot_longer(everything(), names_to = "k", values_to = "Stat") %>%
    dplyr::mutate(Metric = "FrobError") %>%
    dplyr::mutate(k = as.numeric(sub("^k", "", .data$k)))

  metrics_df <- nmf_exp@OptKStats[,-1] %>%
    tidyr::pivot_longer(names_to = "Metric", values_to = "Stat", -.data$k)

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
#' nmf_exp <- runNMFtensor_lite(leukemia$matrix, ranks = 2:10,
#' method = "NMF",
#' n_initializations = 2)
#' plot(generateRiverplot(nmf_exp))
#' plot(generateRiverplot(nmf_exp, ranks = 2:5))
#' }
setMethod("generateRiverplot",
          "nmfExperiment_lite",
          function(nmf_exp, edges.cutoff = 0,
                   useH=FALSE, color=TRUE,
                   ranks=NULL) {
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