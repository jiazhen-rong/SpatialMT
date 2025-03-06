#' Plot Simple Model Bayes Factor
#'
#' This function computes and visualizes the Bayes Factor (BF) for a given variant
#' in spatial mitochondtial data. The BF is calculated using a binomial
#' likelihood ratio test, comparing a baseline probability (`p0`) estimated from
#' all cells to the maximum likelihood estimate (`p1`) for each spatial spot.
#'
#' @param var Character string, the variant of interest.
#' @param N Matrix or data frame, the total counts for each spot.
#' @param af.dm Matrix or data frame, the allele frequency data.
#' @param spatial_coords Data frame or matrix with spatial coordinates (SPATIAL_1, SPATIAL_2).
#' @param BF_cap Numeric, maximum value to cap the Bayes Factor at (default: 3).
#'
#' @return A `ggplot2` object showing the spatial distribution of the Bayes Factor.
#' @export
#'
#' @examples
#' PlotBayesFactor("3054_G>C", N, af.dm, spatial_coords.plot=T)
PlotBayesFactor <- function(var, N, af.dm, spatial_coords,BF_cap=3,plot=F) {
  # Compute baseline probability from across all cells
  p0 <- mean(af.dm[var, ], na.rm = TRUE)
  X_var = round(N[var, ]*af.dm[var,])
  # Compute Bayes Factor
  BF <- log10(dbinom(X_var, size = N[var, ], prob = af.dm[var, ])) -
    log10(dbinom(X_var, size = N[var, ], prob = p0))

  if(plot){
    # Create data frame for plotting
    plot_df <- as.data.frame(cbind(spatial_coords, BF[rownames(spatial_coords)]))
    colnames(plot_df) <- c("SPATIAL_1", "SPATIAL_2", "BF")
    # Cap BF values at 3
    plot_df$BF_capped <- pmin(plot_df$BF, BF_cap)
    # Generate plot
    p<- ggplot(plot_df, aes(x = SPATIAL_1, y = SPATIAL_2)) +
      # First layer: All points in grey
      geom_point(aes(color = BF_capped), size = 0.5) +
      # Second layer: Highlight BF > 1 with larger points
      geom_point(data = subset(plot_df, BF > 1),
                 aes(color = BF_capped), size = 1) +
      #scale_color_gradientn(colors = c("grey", "red")) +
      scale_color_gradientn(colors = c("blue", "grey", "red"),
                            values = c(0, 0.5, 1),  # Ensure 0 is at 0.5 position
                            limits = c(-3, 3)) +
      scale_size_continuous(range = c(0, 5)) +  # Adjust size range
      theme_minimal() +
      ggtitle(paste0("Bayes Factor ", var)) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
  return(BF)
}


# Transformation function mapping raw BF to [0, 1]
transformBF <- function(BF, k = 2, BF_mid = 0) {
  # BF_mid is the midpoint: BF == BF_mid gives a score of 0.5.
  score <- 1 / (1 + exp(-k * (BF - BF_mid)))
  score
}

#' Plot Heatmap of Bayes Factor for Variants of Interest
#'
#' This function computes and visualizes the Bayes Factor (BF) for a list of variants
#' in a spatial transcriptomics dataset. The BF is calculated using a binomial likelihood
#' ratio test, comparing a baseline probability (`p0`) estimated from all cells to the
#' maximum likelihood estimate (`p1`) for each spatial spot.
#'
#' To order the variants by selection process as MAESTER:
#' variants are first ordered by hierarchical clustering and
#' the spots are iteratively re‐ordered such that those with higher BF in each variant group
#' are placed together.
#'
#' @param voi Character vector, list of variants of interest.
#' @param N Matrix (voi x spot), total counts per variant per spot.
#' @param af.dm Matrix (voi x spot), allele frequency data per variant per spot.
#' @param spatial_coords Data frame (spot x 2), spatial coordinates (SPATIAL_1, SPATIAL_2).
#' @param BF_cap Numeric, maximum value to cap the Bayes Factor (default: 3).
#' @param ngroups Numeric, number of groups to cut the variant dendrogram (default: number of variants).
#'
#' @return A heatmap showing the Bayes Factor for each variant.
#' @export
#' @import pheatmap ComplexHeatmap circlize
#' @examples
#' PlotBayesFactorHeatmap(c("3054_G>C", "2021_A>T"), N, af.dm, spatial_coords)
PlotBayesFactorHeatmap <- function(voi, N, af.dm, spatial_coords, BF_cap = 3,
                                   cor_order=T,ngroups = length(voi),method=c("BF","sigmoid")) {

  if(method == "BF"){
    # Compute BF for each variant
    BF_matrix <- sapply(voi, function(var) {
      p0 <- mean(af.dm[var, ], na.rm = TRUE)
      X_var = round(N[var, ]*af.dm[var,])
      BF <- log10(dbinom(X_var, size = N[var, ], prob = af.dm[var, ])) -
        log10(dbinom(X_var, size = N[var, ], prob = p0))
      BF
    })
    BF_df <- as.data.frame(BF_matrix)
    rownames(BF_df) <- colnames(N)  # Use spot names as rownames
  }else if(method == "sigmoid"){
    # compute sigmoid function so final value between -1 to 1.
    BF_matrix <- sapply(voi, function(var) {
      BF = rep(NA,dim(N)[2]); names(BF) = colnames(N)
      p0 <- mean(af.dm[var, ], na.rm = TRUE)
      # Case 1: N=0, X=0
      BF[N[var,] ==0] = 0
      X_var = round(N[var, ]*af.dm[var,])
      # Case 2: X>0, N>0, positive evidence
      spot_idx = ((X_var > 0) & (N[var,]>0))
      BF[spot_idx] <- log10(dbinom(X_var[spot_idx],
                                   size = N[var, spot_idx],
                                   prob = af.dm[var, spot_idx])) -
        log10(dbinom(X_var[spot_idx],
                     size = N[var, spot_idx],
                     prob = p0))
      # the sigmoid transformation convert to 0~1
      # X=0, N > 0
      # the sigmoid transformation convert to -1~0
      # Case 3:
      neg_spot_idx = ((X_var == 0) & (N[var,]>0))
      BF[neg_spot_idx ] <- -log10(N[var,neg_spot_idx] + 1)
      BF
    })
    BF_df <- as.data.frame(BF_matrix)
    rownames(BF_df) <- colnames(N)  # Use spot names as rownames
  }

  # Capping BF values for visualization
  BF_df_capped <- as.data.frame(apply(BF_df, 2, function(x) pmin(x, BF_cap)))

  if(cor_order){ # Ordering as MAESTER, mimicing selection process
    # Compute correlation among variants
    cor_mat <- cor(BF_df, use = "pairwise.complete.obs")
    var_clust <- hclust(as.dist(1 - cor_mat))
    #plot(var_clust$height, ylim = c(0, max(var_clust$height)))
    hm1 <- Heatmap(cor_mat,
                   col = colorRamp2(c(-1,0,1), c("blue", "#DDDDDD", "red")),
                   cluster_columns = var_clust,
                   cluster_rows = var_clust,
                   row_split = switch(ngroups < length(voi), ngroups),
                   column_split = switch(ngroups < length(voi), ngroups),
                   show_row_dend = F, # without this the visualizationn does not complete
                   show_column_dend = F, # without this the visualizationn does not complete
                   row_gap = unit(0.5, "mm"),
                   column_gap = unit(0.5, "mm"),
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   row_title_gp = gpar(fontsize = 10),
                   width = unit(100, "mm"),
                   height = unit(100, "mm"),
                   column_title = ngroups)
    print(hm1)

    # Decide groups of the variants
    groups <- cutree(var_clust, k = min(ngroups, length(voi)))
    ordered_variants <- var_clust$labels[var_clust$order]
    group_list <- split(ordered_variants, groups[ordered_variants])
    group_order <- unique(groups[ordered_variants])
    group_list <- group_list[as.character(group_order)]

    # Subset and order columns accordingly
    BF_ordered <- BF_df_capped[, ordered_variants, drop = FALSE]

    # --- Spot (row) ordering ---
    # Now iterate over the groups in reverse order.
    # For each group, if the group has a single variant, order spots by descending BF of that variant;
    # if multiple variants, order spots by descending row sum across the group.
    for (group_variants in rev(group_list)) {
      if (length(group_variants) == 1) {
        BF_ordered <- BF_ordered[order(-BF_ordered[[group_variants]]), , drop = FALSE]
      } else {
        BF_ordered <- BF_ordered[order(-rowSums(BF_ordered[, group_variants, drop = FALSE])), , drop = FALSE]
      }
    }

    # --- Plotting ---
    # Now plot the ordered matrix without further clustering.
    # Create color breaks. Here we use 101 breaks so that our 100 colors map well.
    my_breaks <- seq(-BF_cap, BF_cap, length.out = 100)

    # Create a palette for the negative side (blue to grey) and for the positive side (grey to red)
    neg_colors <- colorRampPalette(c("blue", "grey"))(50)
    pos_colors <- colorRampPalette(c("grey", "red"))(50)
    my_colors <- c(neg_colors, pos_colors)

    hm2=pheatmap(
      t(as.matrix(BF_ordered)),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      labels_col = F,
      color = my_colors,
      main = "Bayes Factor Heatmap",
      angle_col = 90
    )
    print(hm2)

  }else{
    # Perform hierarchical clustering on spots
    dist_mat <- dist(BF_df)
    cluster_order <- hclust(dist_mat)$order
    ordered_spots <- rownames(BF_df)[cluster_order]

    p2 <- pheatmap(
      BF_df_capped[ordered_spots, ],  # Reorder spots based on clustering
      cluster_rows = FALSE,           # Already clustered, no need for pheatmap clustering
      cluster_cols = TRUE,            # Cluster variants (columns)
      labels_row=F,
      color = my_colors,
      breaks = my_breaks,
      main = "Bayes Factor Heatmap",
      angle_col = 90
    )
    return(p2)
  }
}


#' Plot Heatmap of VAF for Variants of Interest
#'
#' This function visualizes the variant allele frequency (VAF) for a list of variants
#' in a spatial transcriptomics dataset. The VAF values are extracted from `af.dm`
#' and ordered similarly to the BF heatmap example: variants are clustered and grouped,
#' and spots are iteratively re‐ordered by descending VAF values within each variant group.
#'
#' @param voi Character vector, list of variants of interest.
#' @param af.dm Matrix (variants x spot), allele frequency data per variant per spot.
#' @param spatial_coords Data frame (spot x 2), spatial coordinates (SPATIAL_1, SPATIAL_2).
#' @param VAF_cap Numeric, maximum value to cap the VAF for visualization (default: 100).
#' @param ngroups Numeric, number of groups to cut the variant dendrogram (default: number of variants).
#'
#' @return A heatmap showing the VAF for each variant across spots.
#' @export
#' @import pheatmap
#' @examples
#' PlotVAFHeatmap(c("3054_G>C", "2021_A>T"), af.dm, spatial_coords)
PlotVAFHeatmap <- function(voi, af.dm, spatial_coords, VAF_cap = 100, ngroups = length(voi)) {
  # Subset af.dm for the variants of interest.
  # Here, af.dm is assumed to have rows as variants and columns as spots.
  VAF_matrix <- af.dm[voi, ]

  # Transpose so that rows represent spots and columns represent variants.
  VAF_df <- as.data.frame(t(VAF_matrix))

  # Cap VAF values (if necessary) for visualization.
  VAF_df_capped <- as.data.frame(lapply(VAF_df, function(x) pmin(x, VAF_cap)))

  # --- Variant (column) ordering ---
  # Compute the correlation matrix among variants using the capped VAF values.
  cor_mat <- cor(VAF_df_capped, use = "pairwise.complete.obs")
  var_clust <- hclust(as.dist(1 - cor_mat))

  hm1 <- Heatmap(cor_mat,
                 col = colorRamp2(c(-1,0,1), c("blue", "#DDDDDD", "red")),
                 cluster_columns = var_clust,
                 cluster_rows = var_clust,
                 row_split = switch(ngroups < length(voi), ngroups),
                 column_split = switch(ngroups < length(voi), ngroups),
                 show_row_dend = F, # without this the visualizationn does not complete
                 show_column_dend = F, # without this the visualizationn does not complete
                 row_gap = unit(0.5, "mm"),
                 column_gap = unit(0.5, "mm"),
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 row_title_gp = gpar(fontsize = 10),
                 width = unit(100, "mm"),
                 height = unit(100, "mm"),
                 column_title = ngroups)
  print(hm1)


  # Cut the dendrogram into groups (if ngroups is less than the number of variants,
  # some variants will be grouped together).
  groups <- cutree(var_clust, k = min(ngroups, length(voi)))

  # Order variants based on the clustering order.
  ordered_variants <- var_clust$labels[var_clust$order]

  # Create a list of variant groups.
  group_list <- split(ordered_variants, groups[ordered_variants])
  group_order <- unique(groups[ordered_variants])
  group_list <- group_list[as.character(group_order)]

  # Reorder columns (variants) of the VAF matrix accordingly.
  VAF_ordered <- VAF_df_capped[, ordered_variants, drop = FALSE]

  # --- Spot (row) ordering ---
  # Iterate over variant groups in reverse order. For each group, sort spots by:
  # - A single variant: descending VAF for that variant.
  # - Multiple variants: descending sum of VAFs across the variants in that group.
  for (group_variants in rev(group_list)) {
    if (length(group_variants) == 1) {
      VAF_ordered <- VAF_ordered[order(-VAF_ordered[[group_variants]]), , drop = FALSE]
    } else {
      VAF_ordered <- VAF_ordered[order(-rowSums(VAF_ordered[, group_variants, drop = FALSE])), , drop = FALSE]
    }
  }

  # --- Plotting ---
  pheatmap(
    t(as.matrix(VAF_ordered)),
    cluster_rows = FALSE,    # Use custom row ordering.
    cluster_cols = FALSE,    # Use custom column ordering.
    labels_col = FALSE,
    color = colorRampPalette(c("white", "blue"))(50),
    main = "VAF Heatmap",
    angle_col = 90
  )
}
