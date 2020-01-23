#==============================================================================#
#                      NMF-GPU plot generation - FUNCTIONS                     #
#==============================================================================#

#' @rdname generateRiverplot-methods
#' @aliases generateRiverplot,ANY,ANY-method
#'
#' @importFrom YAPSA normalize_df_per_dim
#' @export
#'
#' @examples
#' plot(generateRiverplot(nmf_exp))
#' plot(generateRiverplot(nmf_exp, ranks = 2:5))
setMethod("generateRiverplot",
          "nmfExperiment_lite",
          function(nmf_exp, edges.cutoff = 0,
                   useH=FALSE, color=TRUE,
                   ranks=NULL) {
            #----------------------------------------------------------------------------#
            #                      Retrieve list of matrices                             #
            #----------------------------------------------------------------------------#
            if(useH)
              W_list <- lapply(HMatrix(nmf_exp), t)
            else
              W_list <- WMatrix(nmf_exp)

            #----------------------------------------------------------------------------#
            #                  Build data frame with node names                          #
            #----------------------------------------------------------------------------#
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
            #----------------------------------------------------------------------------#
            #          Build data frame with edges values - NNLS                         #
            #----------------------------------------------------------------------------#
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
              pca <- prcomp(t(do.call(cbind, W_list)))
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