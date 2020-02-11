# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#==============================================================================#
#                      NMF-GPU plot generation - FUNCTIONS                     #
#==============================================================================#
#' Science like figures - ggplot2 theme
#'
#' @return
#'
#' @import ggplot2
#' @export
#'
#' @examples
science_theme <- function() {
  science_theme <-
    theme(panel.grid.major = element_line(size = 0.5, color = "grey"),
          axis.line = element_line(size = 0.7, color = "black"),
          text = element_text(size = 14))
  science_theme <- science_theme +
    background_grid(major = "xy", minor = "none")
  return(science_theme)
}

#' Universal ploting function for Optimal K Statistics
#'
#' @param nmf.exp
#' @param plot.vars quality metrics to be displayed
#' @param optK optimal factorization rank, highlighted in the plot if non-NULL
#'
#' @return
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
plotKStats <- function(nmf.exp,
                       plot.vars = c("FrobError", "cv", "sumSilWidth",
                                  "meanSilWidth", "copheneticCoeff",
                                  "meanAmariDist"),
                       optK = NULL) {
  frobError.df <- melt(as.data.frame(FrobError(nmf.exp)))
  frobError.df[, 1] <- factor(gsub("X", "", frobError.df[, 1]))
  frobError.df <- data.frame("k" = frobError.df[, 1],
                             "variable" = "FrobError",
                             "value" = frobError.df[, 2])
  optKStats.df <- melt(as.data.frame(OptKStats(nmf.exp)), id.vars = "k")
  meanError.df <- optKStats.df[optKStats.df$variable == "mean", ]
  meanError.df$variable <- unique(frobError.df$variable)
  #plot.vars <- as.character(unique(optKStats.df$variable)[-1:-3])
  optKStats.df <- optKStats.df[optKStats.df$variable %in% plot.vars, ]
  optKStats.df <- rbind(frobError.df, optKStats.df)
  optKStats.df$k <- as.numeric(as.character(optKStats.df$k))
  gg.optK <- ggplot() + geom_point(data = optKStats.df, aes(x = k, y = value),
                                   col = "black", size = 0.75)
  gg.optK <- gg.optK + geom_point(data = meanError.df,
                                  aes(x = k, y = value), col = "red")
  gg.optK <- gg.optK + facet_wrap(~variable, scales = "free_y")
  gg.optK <- gg.optK + xlab("K") + ylab("") + theme_bw() + science_theme()
  gg.optK <- gg.optK + theme(strip.background = element_rect(fill = "white"))
  if(!is.null(optK)){
    if(min(optKStats.df$k) <= optK & optK <= max(optKStats.df$k)){
      rect_df <- data.frame(xmin = optK - 0.5, xmax = optK + 0.5,
                            ymin = -Inf, ymax = Inf)
      gg.optK <- gg.optK +
        geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                      ymax = ymax),
                  fill = "grey50", alpha = 0.3, inherit.aes = FALSE)
    }
  }
  return(gg.optK)
}

#' Ploting function for Rank vs. Frobenius Error
#'
#' @param nmf.exp
#'
#' @return
#'
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
plotRankedFrobErrors <- function(nmf.exp) {
  # Prepare frob error data for plotting.
  frobErrorRanks.df <- apply(FrobError(nmf.exp), 2, function(matrix.col) {
    ranks <- order(matrix.col, decreasing = T)
    rank.df <- data.frame(rank = 1:length(ranks),
                          froberror = matrix.col[ranks])
    return(rank.df)
  })
  frobErrorRanks.df <- melt(frobErrorRanks.df, value.name = "frobError",
                            id.vars = "rank")
  grid.labels <- paste("K =", frobErrorRanks.df$L1)
  frobErrorRanks.df$L1 <- factor(grid.labels, levels = unique(grid.labels))
  # Plot ranked frob error
  gg.rankedFrobError <- ggplot(frobErrorRanks.df, aes(x = rank, y = frobError))
  gg.rankedFrobError <-
    gg.rankedFrobError + geom_line(col = "black", size = 0.75)
  gg.rankedFrobError <- gg.rankedFrobError + facet_wrap(~L1, scales = "free")
  gg.rankedFrobError <-
    gg.rankedFrobError + xlab("Error rank") + ylab("Frobenius error")
  gg.rankedFrobError <- gg.rankedFrobError + theme_bw() + theme_cowplot()
  gg.rankedFrobError <- gg.rankedFrobError +
    theme(strip.background = element_rect(fill = "white"))
  return(gg.rankedFrobError)
}

#' Generate ColorMap for ComplexHeatmaps
#'
#' @param matrix.list
#'
#' @return
#'
#' @importFrom circlize colorRamp2
#' @export
#'
#' @examples
getColorMap <- function(matrix.list) {
  matrix.max <- quantile(unlist(matrix.list), 0.95)
  matrix.min <- quantile(unlist(matrix.list), 0.05)
  if(matrix.min > 0) {
    col.map <- colorRamp2(c(0, matrix.max), c("white", "red"))
  } else {
    col.map <- colorRamp2(c(matrix.min, 0, matrix.max),
                          c("blue", "white", "red"))
  }
  return(col.map)
}

#' Computes color annotation for colData from NMF-experiment class
#'
#' @param col.data
#'
#' @return
#'
#' @import RColorBrewer
#' @export
#'
#' @examples
getColorAnno <- function(col.data) {
  # Define own palate
  own.cols <- c("red", "darkblue", "orange", "darkgreen")
  # Define usable qualitative colour palates from RColorBreweer
  qual.colPals <- sort(c("Dark2" = 8, "Paired" = 12, "Pastel1" = 9,
                         "Pastel2" = 8, "Set1" = 9, "Set2" = 8, "Set3" = 12,
                         "own" = length(own.cols)))  # "Accent" = 8,
  # Build color annotation for complexHeatmap
  anno.colours <- lapply(seq(2, ncol(col.data)), function(i) {
    n.levels <- length(unique(col.data[, i]))
    col.pal <- names(qual.colPals)[which(qual.colPals > n.levels)[i - 1]]
    if( col.pal != "own" ) {
      anno.col <- brewer.pal(qual.colPals[col.pal], col.pal)
      anno.col <- anno.col[1:n.levels]
    } else {
      anno.col <- own.cols[1:n.levels]
    }
    names(anno.col) <- unique(col.data[, i])
    return(anno.col)
  })
  names(anno.colours) <- colnames(col.data)[2:ncol(col.data)]
  return(anno.colours)
}

#' Plot H-Matrix as heatmap for a given factorization rank K
#'
#' @param nmf.exp
#' @param k
#'
#' @return
#'
#' @import ComplexHeatmap
#' @export
#'
#' @examples
plotHMatrix <- function(nmf.exp, k = 2){
  # Get matrix
  m <- HMatrix(nmf.exp, k = k)
  # Get signature names
  signature.names <- getSignatureNames(nmf.exp, k)
  signature.names <- paste("Signature:", signature.names, sep = "\n")
  # Define colormap for heatmaps
  col.map <- getColorMap(m)
  # Create heatmap annotation
  n.cols <- ncol(colData(nmf.exp))
  colnames(m) <- colData(nmf.exp)[, 1]
  rownames(m) <- signature.names
  col.anno <- getColorAnno(colData(nmf.exp))
  heat.anno <- HeatmapAnnotation(colData(nmf.exp)[, 2:n.cols], col = col.anno)
  H.heatmap <- Heatmap(matrix = m,
                       col = col.map,
                       cluster_rows = F,
                       top_annotation = heat.anno,
                       clustering_distance_rows = "pearson",
                       heatmap_legend_param = list(color_bar = "continuous"),
                       column_title = sprintf("H-Matrix for K = %s", k))
  return(H.heatmap)
}


#' OLD VERSION MIGHT BE REMOVABLE
#' Plot H or W-Matrix as heatmaps for a given range of K's
#'
#' @param matrix.list
#' @param titles
#' @param trans
#' @param cluster_columns
#' @param show.row
#'
#' @return
#'
#' @import ComplexHeatmap
#' @export
#'
#' @examples
plotHeatmap4MatrixList <- function(nmf.exp, H = T, W = F, titles = NULL,
                                   trans = F, cluster_columns = F,
                                   show.row = F) {
  # Get Matrix List
  matrix.list <- HMatrix(nmf.exp)
  if(trans) matrix.list <- lapply(matrix.list, t)
  # Define colormap for heatmaps
  col.map <- getColorMap(matrix.list)
  # Create annotations
  heat.anno <- HeatmapAnnotation(colData(nmf.exp))
  # Plot heatmaps
  heatmap.list <- Heatmap(matrix.list[[1]],
                          col = col.map,
                          cluster_columns = F,
                          show_row_names = show.row,
                          clustering_distance_rows = "pearson",
                          heatmap_legend_param = list(color_bar = "continuous"),
                          #column_title = colData(nmf.exp))
                          column_title = titles[1])
  for (i in 2:length(matrix.list)){
    heatmap.list <- heatmap.list + Heatmap(matrix.list[[i]],
                                           col = col.map,
                                           cluster_columns = F,
                                           show_row_names = show.row,
                                           clustering_distance_rows = "pearson",
                                           show_heatmap_legend = FALSE,
                                           column_title = titles[i])
  }
  return(heatmap.list)
}


#' Title
#'
#' @param nmf.exp
#' @param sig.combs
#' @param col
#'
#' @return
#'
#' @export
#' @import reshape2
#' @import ggplot2
#' @import cowplot
#'
#' @examples
plotSignatureFeatures <- function(nmf.exp, sig.combs = T, col = "blue") {
  # Get number of features per signature combination and create a barplot.
  n.peaks <- as.data.frame(table(FeatureStats(nmf.exp)[, 1]))
  colnames(n.peaks) <- c("sigCombId", "value")
  # Build annotation matrix and order it
  anno.matrix <-
    as.data.frame(do.call(rbind, strsplit(as.character(n.peaks$sigCombId),
                                          split = "")))
  cluster.matrix <- apply(anno.matrix, 1, as.numeric)
  # If only signatures should be plotted, filter here.
  if(!sig.combs) {
    i.keep <- which(apply(cluster.matrix, 2, sum) == 1)
    anno.matrix <- anno.matrix[i.keep, ]
    cluster.matrix <- cluster.matrix[, i.keep]
    n.peaks <- n.peaks[i.keep, ]
  }
  i.order <- unlist(orderBinary(cluster.matrix))
  anno.matrix <- t(anno.matrix)[, i.order]
  colnames(anno.matrix) <- 1:ncol(anno.matrix)
  rownames(anno.matrix) <- signature.names
  anno.matrix <- melt(anno.matrix)
  # Reorder number of features data.frame
  n.peaks$sigCombId <- factor(n.peaks$sigCombId,
                              levels = n.peaks$sigCombId[i.order])
  # Prepare heatmap annotation plot for signature feature plot
  gg.heat <- ggplot(anno.matrix,
                    aes_string(x = "Var2", y = "Var1", fill = "value"))
  gg.heat <- gg.heat + geom_tile(col = "black", size = 0.5)
  gg.heat <- gg.heat + scale_fill_manual(values = c("white", "black"))
  gg.heat <- gg.heat + ylab("Signatures") + xlab ("Signature combinations")
  gg.heat <- gg.heat + scale_x_discrete(expand = c(10^(-2), 10^(-2)))
  gg.heat <- gg.heat + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL)
  gg.heat <- gg.heat + theme_bw() + theme_cowplot() +
    theme(legend.position = "none", axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.line = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank())
  # Prepare barplot for signature feature plot
  gg.bar <- ggplot(n.peaks, aes_string(x = "sigCombId", y = "value"))#,
                                     # fill = clusterId))
  gg.bar <- gg.bar + geom_bar(colour = "black", fill = col, stat = "identity",
                              position = position_dodge())
  gg.bar <- gg.bar + xlab("") + ylab("#Features")
  gg.bar <- gg.bar + theme_bw() + theme_cowplot()
  gg.bar <- gg.bar + scale_x_discrete(expand = c(10^(-2), 10^(-2)))
  gg.bar <- gg.bar + labs(x = NULL)
  gg.bar <- gg.bar + theme(axis.ticks.x = element_blank(),
                           axis.text.x = element_blank())
  # Combine annotation heatmap and barplot
  gg.feature <- plot_grid(gg.bar, gg.heat, nrow = 2, align = "v")
  return(gg.feature)
}

#=============================================================================#
#  Riverplot - similarities between signatures at different ranks - FUNCTIONS #
#=============================================================================#
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
    X[i, ] <- coef(nnls(A, B[, i]))
  X
}

#' Order the riverplot nodes to minimize crossings
#'
#' @param nodes
#' @param edges
#' @param rank_flag
#' @param scale_factor
#' @param start_coords
#' @param edgeWeightColName
#' @param sacle_fun
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
#' @param nodes
#' @param edges
#'
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

#' NMF riverplot
#'
#' Generate a riverplot object that displays the similarities between
#' signatures at different factorization ranks
#'
#' @param nmf_exp an nmfExperiment or a nmfExperiment_lite object
#' @param edges.cutoff cutoff until which similarities are displayed
#' @param useH whether to relate signatures (FALSE) or exposures (TRUE)
#' @param color whether to colorize the nodes based on PCA of the signatures
#'
#' @return ret riverplot object
#'
#' @import riverplot
#' @export
#'
#' @docType methods
#' @rdname generateRiverplot-methods
#'
#' @examples plt <- generateRiverplot(nmf_exp, edges.cutoff = 0.2)
#' plot(plt, plot_area = 1, yscale = 0.6, nodewidth = 0.5)
setGeneric("generateRiverplot", function(nmf_exp, ...) standardGeneric("generateRiverplot"))

#' @rdname generateRiverplot-methods
#' @aliases generateRiverplot,ANY,ANY-method
#'
#' @export
#'
setMethod("generateRiverplot",
          "nmfExperiment",
          function(nmf_exp, edges.cutoff = 0, useH=FALSE,
                   color=TRUE) {
            if(useH)
              W.list <- lapply(HMatrix(nmf_exp), t)
            else
              W.list <- WMatrix(nmf_exp)
            rank.max <- max(as.numeric(names(W.list)))
            rank.min <- min(as.numeric(names(W.list)))
            nodes <- do.call(rbind, lapply(rank.min:rank.max, function(i) {
              data.frame(ID = sapply(1:i, function(n) paste(i, "_S", n, sep = "")),
                         x = i)
            }))
            edges <- do.call(rbind, lapply(rank.min:(rank.max - 1), function(i){
              df <- data.frame(
                N1 = rep(sapply(1:i, function(n) paste(i, "_S", n, sep = "")),
                         each = i + 1),
                N2 = rep(sapply(1:(i + 1), function(n) paste(i + 1, "_S", n, sep = "")),
                         time = i),
                Value = as.vector(t(nnls_sol(W.list[[as.character(i)]], W.list[[as.character(i + 1)]])))
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
              pca <- prcomp(t(do.call(cbind, W.list)))
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


#' Preparation for relabelling and recolouring of a riverplot
#'
#' @param in_nmf.exp
#' @param in_signatures_df
#' @param in_sigInd_df
#' @param in_normalize
#' @param in_optK If non-NULL, specifies the optimal factorization rank.
#'  The signatures at this rank are the reference signatures, all other
#'  signatures will get the suffix "-like"
#'
#' @return
#' @export
#'
#' @examples
attributeComparisonSignatures <- function(in_nmf.exp, in_signatures_df,
                                          in_sigInd_df, in_normalize = TRUE,
                                          in_optK = NULL){
  my_NMFlistsList <- translateBratwurstToYAPSA(in_nmf.exp,
                                               normalize = in_normalize)
  new_NMFlistsList <-
    deriveSigInd_df(in_signatures_df, my_NMFlistsList, in_sigInd_df)
  sigNameVecList <-
    lapply(names(new_NMFlistsList), function(current_rank) {
      current_list <- new_NMFlistsList[[current_rank]]
      current_df <- current_list$out_sig_ind
      currentVec <- as.character(current_df$match)
      names(currentVec) <- paste0(current_rank, "_",
                                  as.character(current_df$sig))
      return(currentVec)
    })
  sigNameVec <- do.call(c, sigNameVecList)
  compositeNameVec <- paste0(names(sigNameVec), "\n", sigNameVec)
  names(compositeNameVec) <- names(sigNameVec)
  if(!is.null(in_optK)){
    other_ind <- which(unlist(lapply(
      compositeNameVec, function(x)
        strsplit(x, split = "_")[[1]][1])) != in_optK)
    compositeNameVec[other_ind] <- paste0(compositeNameVec[other_ind], "-\nlike")
  }
  sigColVector <- in_sigInd_df$colour
  names(sigColVector) <- in_sigInd_df$sig
  return(list(sig_names = sigNameVec,
              name_vector = compositeNameVec,
              col_vector = sigColVector))
}

#' Relabel and recolour a riverplot
#'
#' @param in_riverplot
#' @param in_list
#'
#' @return
#' @export
#'
#' @examples
relabelRiverplot <- function(in_riverplot, in_list){
  in_riverplot$nodes$labels <-
    in_list$name_vector[as.character(in_riverplot$nodes$ID)]
  tempVec <-unlist(lapply(strsplit(
    in_list$sig_names[as.character(in_riverplot$nodes$ID)], split = "_"),
    head, 1))
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
