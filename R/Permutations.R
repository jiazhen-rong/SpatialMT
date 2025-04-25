#' Perform Celltype Test Using Nearest Neighbors for Each Spot
#'
#' This function permutes the VAF and celltype proportion relationships
#' and then perform spot level celltype proportion test with neighbors.
#'
#' @param celltypes A string vector of unique cell types.
#' @param voi A string vector of variants of interest.
#' @param N A matrix of total coverage counts (MT variant x cell).
#' @param vaf A matrix of variant allele frequencies (MT variant x cell).
#' @param Ws A matrix of cell type weights (cell x celltype).
#' @param spatial_coords_temp A matrix or data frame of spatial coordinates (cell x 2) with row names as cell IDs.
#' @param test_type A string specifying the regression test: "linear" or "weighted".
#' @param permute_num Number of permutations for the weighted test (default 1000).
#' @param k_neighbors Number of nearest neighbors to use for each spot (default = 100).
#' @param method P-value correction method. Options are "FDR", "FWER", or "Raw" (default = "Raw").
#'
#' @return A list with one element per spot (named by the spot ID). Each element is a list containing:
#'   \item{intercept}{A matrix of intercept estimates (rows = variants, columns = cell types).}
#'   \item{coef}{A matrix of regression coefficients.}
#'   \item{pval}{A matrix of raw p-values.}
#'   \item{adjusted_pval}{A matrix of adjusted p-values (if method != "Raw"); otherwise, same as raw.}
#'
#' @examples
#' \dontrun{
#'   # results_knn <- celltype_test_knn_permute(celltypes, voi, N, vaf, Ws, spatial_coords_temp,
#'   #                                  test_type = "linear", permute_num = 1000,
#'   #                                  k_neighbors = 100, method = "FDR")
#' }
#'
#' @export
celltype_test_knn_permute <- function(celltypes, vars, N_voi, vaf, Ws, spatial_coords_temp,
                                      test_type = c("linear", "weighted"),
                                      permute_num = 1000, k_neighbors = 100,
                                      method = "Raw", sample_idx_local=NULL,disease_celltype="BE",
                                      ratio_threshold=0.03,
                                      vaf_cellprop=F,exclude_plot_idx=NULL,
                                      p_thresh=0.05,max_log10p_cap=6,coef_plot_option="negative",plot=T) {

  # Define Continuous Color Scale
  if(coef_plot_option=="grey"){
    pval_palette <- scale_color_gradientn(
      colors = c("darkgrey","grey","red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(0, -log10(p_thresh), max_log10p_cap)),  # Map p-values to colors
      limits = c(0, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }else if(coef_plot_option=="negative"){
    pval_palette <- scale_color_gradientn(
      colors = c("blue","grey", "red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(-max_log10p_cap, 0, max_log10p_cap)),  # Map p-values to colors
      limits = c(-max_log10p_cap, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }

  # Find common cells across datasets
  common_cells <- intersect(intersect(colnames(N_voi), rownames(Ws)), rownames(spatial_coords_temp))
  if(length(common_cells) == 0) stop("No common cells found across datasets.")

  N_voi <- N_voi[, common_cells, drop = FALSE]
  vaf <- vaf[, common_cells, drop = FALSE]
  Ws <- Ws[common_cells, , drop = FALSE]
  spatial_coords_temp = spatial_coords_temp[common_cells,]

  # break the relationship between vaf (X,N remains the same) and Ws
  sample_idx = sample(dim(N_voi)[2])
  N_voi_temp = N_voi[,sample_idx,drop=F];colnames(N_voi_temp) = colnames(N_voi);
  N_voi = N_voi_temp
  vaf_temp <- vaf[, sample_idx, drop = FALSE];colnames(vaf_temp) = colnames(vaf);
  vaf = vaf_temp

  spatial_coords_temp <- as.data.frame(spatial_coords_temp[common_cells, ])
  colnames(spatial_coords_temp) = c("X","Y")

  # Use FNN::get.knn to obtain indices of k nearest neighbors for each cell
  knn_res <- FNN::get.knn(as.matrix(spatial_coords_temp[, c("X", "Y")]), k = k_neighbors)

  results_list <- list()  # will store results for each spot
  all_pvals <- c()        # to collect p-values for multiple testing correction

  # Loop over each spot (each row in spatial_coords_temp)
  for (i in seq_len(nrow(spatial_coords_temp))) {
    if(i%%1000==0 ){
      print(paste0(i,"/", nrow(spatial_coords_temp)," spots being processed."))
    }
    spot_id <- rownames(spatial_coords_temp)[i]
    neighbor_idx <- knn_res$nn.index[i, ]  # indices of nearest neighbors
    neighbor_ids <- rownames(spatial_coords_temp)[neighbor_idx]
    test_ids = c(neighbor_ids,spot_id ) # include the spot for test

    # Subset matrices for these neighbors
    N_neighbors <- N_voi_temp[, test_ids, drop = FALSE]
    vaf_neighbors <- vaf_temp[, test_ids, drop = FALSE]
    Ws_neighbors <- Ws[test_ids, , drop = FALSE]

    # Initialize matrices for this spot's regression results
    intercept_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    coef_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    pval_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    rownames(intercept_mat) <- vars
    colnames(intercept_mat) <- celltypes
    rownames(coef_mat) <- vars
    colnames(coef_mat) <- celltypes
    rownames(pval_mat) <- vars
    colnames(pval_mat) <- celltypes

    # Loop over each variant
    for (j in seq_along(vars)) {
      var <- vars[j]
      N_j <- N_neighbors[j, ]
      vaf_j <- vaf_neighbors[j, ]
      Ws_j = Ws_neighbors

      if(sum(vaf_j > 0)/length(vaf_j) < ratio_threshold){ # excluding outliers
        intercept_mat[j, ] <- 0
        coef_mat[j, ] <- 0
        pval_mat[j, ] <- 1
        next
      }

      # If a set of control spots is desired, determine sample_idx_local
      if(is.null(sample_idx_local) && !is.null(disease_celltype)) {
        sample_idx_local <- intersect(
          rownames(Ws)[(Ws[, disease_celltype] < 0.3) & !is.na(Ws[, disease_celltype])],
          colnames(N_voi)[N_voi[var, ] > 1]
        )
      }

      #print(paste0(length(sample_idx_local), " paired control spots."))
      # Append control spots if available
      if(length(sample_idx_local) > 0){
        N_j <- append(N_j, N_voi[j, sample_idx_local])
        vaf_j <- append(vaf_j, vaf[j, sample_idx_local])
        Ws_j = rbind(Ws_neighbors,Ws[sample_idx_local,])
      }

      for (k in seq_along(celltypes)) {
        data_df <- data.frame(vaf_j = vaf_j,
                              N_j = N_j,
                              W_sk = Ws_j[, celltypes[k]])

        # Check for variation in predictors
        if (length(unique(data_df$W_sk)) == 1 || length(unique(data_df$vaf_j)) == 1) {
          next
        }

        # Run regression (linear or weighted)
        if (test_type == "linear") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data_df), error = function(e) NULL)
        } else if (test_type == "weighted") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data_df, weights = sqrt(N_j)), error = function(e) NULL)
        }

        if (is.null(fit)) next

        res <- summary(fit)
        if (nrow(res$coefficients) < 2) next

        intercept_mat[j, k] <- res$coefficients[1, 1]
        coef_mat[j, k] <- res$coefficients[2, 1]
        pval_mat[j, k] <- res$coefficients[2, 4]

        all_pvals <- c(all_pvals, pval_mat[j, k])
      } # end celltype loop

      # plot vaf vs celltype diagnostic for grid
      if(vaf_cellprop){
        intercept=intercept_mat[j,]
        coef=coef_mat[j,]
        pval=pval_mat[j,]
        #print(pval)
        temp_plot_df = as.data.frame(cbind(spatial_coords_temp[c(test_ids,sample_idx_local),c("X","Y")],
                                           vaf_j,N_j,Ws_j))
        colnames(temp_plot_df) = c("X","Y",var,paste0(var,"_cov"),celltypes)
        grid_vaf_cellprop_plot[[grid]] = plot_vaf_cellprop(j,
                                                           vaf[,c(test_ids,sample_idx_local),drop=F],Ws_j,temp_plot_df,
                                                           intercept,coef,pval,permute,return_plot=T)
      }
    } # end variant loop

    results_list[[spot_id]] <- list(intercept = intercept_mat,
                                    coef = coef_mat,
                                    pval = pval_mat)
  } # end spot loop
  # Plot celltype diagnostic plots
  if(plot){
    plot_list  = list()
    # Generate Plots for Each celltype
    for (celltype in celltypes) {
      # For each spot (each element in results_list), extract the adjusted p-value and coefficient for the given variant and cell type.
      spot_values <- sapply(names(results_list), function(spot_id) {
        res <- results_list[[spot_id]]
        #if (is.null(res)) return(NA)
        #if (!(variant %in% rownames(res$adjusted_pval))) return(NA)
        #if (!(ct %in% colnames(res$adjusted_pval))) return(NA)

        p_val <- res$pval[, celltype]
        coef_val <- res$coef[, celltype]
        # Compute signed -log10(adjusted p)
        signed_logp <- sign(coef_val) * (-log10(p_val))
        return(signed_logp)

      })
      spatial_coords_temp$signed_logP = spot_values
      spatial_coords_temp$capped_signed_logP =  pmax(pmin(spot_values, max_log10p_cap),
                                                     -max_log10p_cap)
      # exclude spots with low diseased celltype proportions
      final_plot_temp = spatial_coords_temp
      final_plot_temp[exclude_plot_idx,"capped_signed_logP"] = 0

      p <- ggplot(final_plot_temp, aes(x = X, y = Y,color = capped_signed_logP)) +
        geom_point(size=0.1) +
        pval_palette +
        labs(title = paste(var, "-", celltype, "-", " P value"),
             x = "X", y = "Y", fill = "-log10(P-value)") +
        theme_minimal() +
        coord_fixed()  # Keep aspect ratio square

      # Store plot in list only if it has data
      plot_list[[paste(var, celltype, sep = "_")]] <- p
    }

    # Arrange Plots in a Grid
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ceiling(sqrt(length(celltypes))))  # Adjust ncol if needed
  }else{
    combined_plot=NULL
  }
  return(list(results=results_list,
              spatial_coords = spatial_coords_temp,
              combined_plot=combined_plot))
}


#' Compute a Simple Path for the spots Using MST and DFS.
#'
#' This function computes a traversal path through spatial coordinates using a
#' Minimum Spanning Tree (MST) and Depth-First Search (DFS). The path is not
#' necessarily the shortest but provides a simple and efficient ordering.
#'
#' @param spatial_coords A matrix or data frame with two columns representing
#'        the spatial coordinates (e.g., X and Y positions).
#' @param plot Logical (default: TRUE). If TRUE, plots the path over the spatial coordinates.
#'
#' @return A vector of indices representing the order in which the coordinates are visited.
#'
#' @examples
#' spatial_coords <- matrix(runif(20), ncol=2)
#' path_order <- mst_dfs_path(spatial_coords, plot=TRUE)
#'
#' @import igraph
#' @export
UpdatePath <- function(spatial_coords_temp,common_cells,plot=T) {

  if(length(common_cells) == 0) stop("No common cells found across datasets.")
  spatial_coords_temp = spatial_coords_temp[common_cells,]

  # Create distance matrix
  dist_matrix <- as.matrix(dist(spatial_coords_temp))

  # Convert to graph and compute MST
  g <- graph_from_adjacency_matrix(dist_matrix, mode = "undirected", weighted = TRUE)
  mst <- mst(g)

  # Perform DFS traversal from the first node
  dfs_result <- dfs(mst, root = 1, unreachable = FALSE)
  path_order <- dfs_result$order
  path_order_barcodes = rownames(spatial_coords_temp)[path_order]
  # Reorder spatial coordinates based on path
  if(plot){
    # Plot the path
    ordered_df <- as.data.frame(spatial_coords_temp[path_order, ])
    plot(spatial_coords_temp, pch = 16, main = "MST-DFS Path")
    lines(ordered_df[, 1], ordered_df[, 2], col = "red")  # Adjusted for generic column names
  }
  return(path_order_barcodes)  # Return the ordered coordinates
}


#' Perform Celltype Test Using Nearest Neighbors for Each Spot.
#' Sequentially Update regression result with biglm to speed up.
#'
#' This function permutes the VAF and celltype proportion relationships
#' and then perform spot level celltype proportion test with neighbors.
#'
#' @param celltypes A string vector of unique cell types.
#' @param voi A string vector of variants of interest.
#' @param N A matrix of total coverage counts (MT variant x cell).
#' @param vaf A matrix of variant allele frequencies (MT variant x cell).
#' @param Ws A matrix of cell type weights (cell x celltype).
#' @param spatial_coords_temp A matrix or data frame of spatial coordinates (cell x 2) with row names as cell IDs.
#' @param test_type A string specifying the regression test: "linear" or "weighted".
#' @param permute_num Number of permutations for the weighted test (default 1000).
#' @param k_neighbors Number of nearest neighbors to use for each spot (default = 100).
#' @param method P-value correction method. Options are "FDR", "FWER", or "Raw" (default = "Raw").
#'
#' @return A list with one element per spot (named by the spot ID). Each element is a list containing:
#'   \item{intercept}{A matrix of intercept estimates (rows = variants, columns = cell types).}
#'   \item{coef}{A matrix of regression coefficients.}
#'   \item{pval}{A matrix of raw p-values.}
#'   \item{adjusted_pval}{A matrix of adjusted p-values (if method != "Raw"); otherwise, same as raw.}
#'
#' @examples
#' \dontrun{
#'   # results_knn <- celltype_test_knn_permute(celltypes, voi, N, vaf, Ws, spatial_coords_temp,
#'   #                                  test_type = "linear", permute_num = 1000,
#'   #                                  k_neighbors = 100, method = "FDR")
#' }
#'
#' @import biglm
#' @export
celltype_test_knn_permute_fastupdate <- function(celltypes, vars, N_voi, vaf, Ws, spatial_coords_temp,
                                      path_order_barcodes,
                                      test_type = c("linear", "weighted"),
                                      permute_num = 1000, k_neighbors = 100,
                                      method = "Raw", sample_idx=NULL,disease_celltype="BE",
                                      ratio_threshold=0.03,
                                      vaf_cellprop=F,exclude_plot_idx=NULL,
                                      p_thresh=0.05,max_log10p_cap=6,coef_plot_option="negative",plot=T) {

  # Define Continuous Color Scale
  if(coef_plot_option=="grey"){
    pval_palette <- scale_color_gradientn(
      colors = c("darkgrey","grey","red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(0, -log10(p_thresh), max_log10p_cap)),  # Map p-values to colors
      limits = c(0, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }else if(coef_plot_option=="negative"){
    pval_palette <- scale_color_gradientn(
      colors = c("blue","grey", "red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(-max_log10p_cap, 0, max_log10p_cap)),  # Map p-values to colors
      limits = c(-max_log10p_cap, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }

  # Find common cells across datasets
  common_cells <- intersect(intersect(colnames(N_voi), rownames(Ws)), rownames(spatial_coords_temp))
  if(length(common_cells) == 0) stop("No common cells found across datasets.")

  N_voi <- N_voi[, common_cells, drop = FALSE]
  vaf <- vaf[, common_cells, drop = FALSE]
  Ws <- Ws[common_cells, , drop = FALSE]

  # break the relationship between vaf (X,N remains the same) and Ws
  sample_idx = sample(dim(N_voi)[2])
  N_voi_temp = N_voi[,sample_idx,drop=F];colnames(N_voi_temp) = colnames(N_voi);
  N_voi = N_voi_temp
  vaf_temp <- vaf[, sample_idx, drop = FALSE];colnames(vaf_temp) = colnames(vaf);
  vaf = vaf_temp

  spatial_coords_temp <- as.data.frame(spatial_coords_temp[common_cells, ])
  colnames(spatial_coords_temp) = c("X","Y")

  # Use FNN::get.knn to obtain indices of k nearest neighbors for each cell
  knn_res <- FNN::get.knn(as.matrix(spatial_coords_temp[, c("X", "Y")]), k = k_neighbors)

  results_list <- list()  # will store results for each spot
  all_pvals <- c()        # to collect p-values for multiple testing correction
  previous_models <- list()  # Store previous models for incremental updates

  # Loop over each spot (each row in spatial_coords_temp)
  for (i in seq_along(path_order_barcodes)) {
    if(i%%1000==0){
      print(paste0(i,"/", nrow(spatial_coords_temp)," spots being processed."))
    }
    spot_id <- path_order_barcodes[i]
    spot_idx = which(rownames(spatial_coords_temp) == spot_id)
    neighbor_idx <- knn_res$nn.index[spot_idx, ]  # indices of nearest neighbors
    neighbor_ids <- rownames(spatial_coords_temp)[neighbor_idx]
    test_ids = c(neighbor_ids,spot_id) # include the spot for test

    # Subset matrices for these neighbors
    N_neighbors <- N_voi_temp[, test_ids, drop = FALSE]
    vaf_neighbors <- vaf_temp[, test_ids, drop = FALSE]
    Ws_neighbors <- Ws[test_ids, , drop = FALSE]

    # Initialize matrices for this spot's regression results
    intercept_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    coef_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    pval_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    rownames(intercept_mat) <- vars
    colnames(intercept_mat) <- celltypes
    rownames(coef_mat) <- vars
    colnames(coef_mat) <- celltypes
    rownames(pval_mat) <- vars
    colnames(pval_mat) <- celltypes

    # Loop over each variant
    for (j in seq_along(vars)) {
      var <- vars[j]
      N_j <- N_neighbors[j, ]
      vaf_j <- vaf_neighbors[j, ]
      Ws_j = Ws_neighbors

      if(sum(vaf_j > 0)/length(vaf_j) < ratio_threshold){ # excluding outliers
        intercept_mat[j, ] <- 0
        coef_mat[j, ] <- 0
        pval_mat[j, ] <- 1
        next
      }

      # If a set of control spots is desired, determine sample_idx_local
      if(is.null(sample_idx) && !is.null(disease_celltype)) {
        sample_idx_local <- intersect(
          rownames(Ws)[(Ws[, disease_celltype] < 0.3) & !is.na(Ws[, disease_celltype])],
          colnames(N_voi)[N_voi[var, ] > 1]
        )
      } else {
        sample_idx_local <- sample_idx
      }

      #print(paste0(length(sample_idx_local), " paired control spots."))
      # Append control spots if available
      if(length(sample_idx_local) > 0){
        N_j <- append(N_j, N_voi[j, sample_idx_local])
        vaf_j <- append(vaf_j, vaf[j, sample_idx_local])
        Ws_j = rbind(Ws_neighbors,Ws[sample_idx_local,])
      }

      for (k in seq_along(celltypes)) {
        data_df <- data.frame(vaf_j = vaf_j,
                              N_j = N_j,
                              W_sk = Ws_j[, celltypes[k]])

        # Check for variation in predictors
        if (length(unique(data_df$W_sk)) == 1 || length(unique(data_df$vaf_j)) == 1) {
          next
        }

        if (i == 1 || is.null(previous_models[[var]][[celltype]])) {
          # First time fitting this variant-celltype pair, use biglm()
          formula <- as.formula("vaf_j ~ W_sk")
          if (test_type == "linear") {
            fit <- tryCatch(biglm(formula, data = data_df), error = function(e) NULL)
          } else if (test_type == "weighted") {
            fit <- tryCatch(biglm(formula, data = data_df, weights = sqrt(N_j)), error = function(e) NULL)
          }
        } else {
          # Incrementally update the model
          fit <- tryCatch(update(previous_models[[var]][[celltype]], data = data_df), error = function(e) NULL)
        }


        if (is.null(fit)) next

        res <- summary(fit)
        if (nrow(res$mat) < 2) next

        intercept_mat[j, k] <- res$mat[1, 1]
        coef_mat[j, k] <- res$mat[2, 1]
        pval_mat[j, k] <- res$mat[2, 5]

        all_pvals <- c(all_pvals, pval_mat[j, k])
      } # end celltype loop

      # plot vaf vs celltype diagnostic for grid
      if(vaf_cellprop){
        intercept=intercept_mat[j,]
        coef=coef_mat[j,]
        pval=pval_mat[j,]
        #print(pval)
        temp_plot_df = as.data.frame(cbind(spatial_coords_temp[c(test_ids,sample_idx_local),c("X","Y")],
                                           vaf_j,N_j,Ws_j))
        colnames(temp_plot_df) = c("X","Y",var,paste0(var,"_cov"),celltypes)
        grid_vaf_cellprop_plot[[grid]] = plot_vaf_cellprop(j,
                                                           vaf[,c(test_ids,sample_idx_local),drop=F],Ws_j,temp_plot_df,
                                                           intercept,coef,pval,permute,return_plot=T)
      }
    } # end variant loop

    results_list[[spot_id]] <- list(intercept = intercept_mat,
                                    coef = coef_mat,
                                    pval = pval_mat)
  } # end spot loop
  # Plot celltype diagnostic plots
  if(plot){
    plot_list  = list()
    # Generate Plots for Each celltype
    for (celltype in celltypes) {
      # For each spot (each element in results_list), extract the adjusted p-value and coefficient for the given variant and cell type.
      spot_values <- sapply(names(results_list), function(spot_id) {
        res <- results_list[[spot_id]]
        #if (is.null(res)) return(NA)
        #if (!(variant %in% rownames(res$adjusted_pval))) return(NA)
        #if (!(ct %in% colnames(res$adjusted_pval))) return(NA)

        p_val <- res$pval[, celltype]
        coef_val <- res$coef[, celltype]
        # Compute signed -log10(adjusted p)
        signed_logp <- sign(coef_val) * (-log10(p_val))
        return(signed_logp)

      })
      spatial_coords_temp$signed_logP = spot_values
      spatial_coords_temp$capped_signed_logP =  pmax(pmin(spot_values, max_log10p_cap),
                                                     -max_log10p_cap)
      # exclude spots with low diseased celltype proportions
      final_plot_temp = spatial_coords_temp
      final_plot_temp[exclude_plot_idx,"capped_signed_logP"] = 0

      p <- ggplot(final_plot_temp, aes(x = X, y = Y,color = capped_signed_logP)) +
        geom_point(size=0.1) +
        pval_palette +
        labs(title = paste(var, "-", celltype, "-", " P value"),
             x = "X", y = "Y", fill = "-log10(P-value)") +
        theme_minimal() +
        coord_fixed()  # Keep aspect ratio square

      # Store plot in list only if it has data
      plot_list[[paste(var, celltype, sep = "_")]] <- p
    }

    # Arrange Plots in a Grid
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ceiling(sqrt(length(celltypes))))  # Adjust ncol if needed
  }else{
    combined_plot=NULL
  }
  return(list(results=results_list,
              spatial_coords = spatial_coords_temp,
              combined_plot=combined_plot))
}



