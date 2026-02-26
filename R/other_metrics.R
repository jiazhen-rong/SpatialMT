#' Pairwise Spatial LISI-like Mixing Metric
#'
#' Computes a pairwise spatial mixing score between all combinations of
#' variants of interest (VOI). For each variant pair, the function evaluates
#' whether spatial k-nearest neighbors of variant-positive spots contain
#' signals from one or both variants.
#'
#' The returned metric summarizes:
#' - var_12: Mixing score from variant 1 perspective
#' - var_21: Mixing score from variant 2 perspective
#' - avg: Average of both directions
#'
#' @param voi Character vector of variant IDs to analyze.
#' @param N Numeric matrix of coverage counts (variants x spots).
#' @param af.dm Numeric matrix of allele frequencies (variants x spots).
#' @param spatial_coords Numeric matrix of spatial coordinates
#'   (rows = spots, columns = spatial dimensions).
#' @param coverage_thresh Minimum coverage required to include a spot.
#'   Default is 0.
#' @param k Number of nearest neighbors for kNN graph. Default is 50.
#'
#' @return A named list where each element corresponds to a variant pair
#'   and contains a numeric vector with:
#'   \describe{
#'     \item{var_12}{Mean mixing score from variant 1 perspective}
#'     \item{var_21}{Mean mixing score from variant 2 perspective}
#'     \item{avg}{Average of the two directional scores}
#'   }
#'
#' @importFrom FNN get.knn
#' @export
pairwise_LISI <- function(voi=NULL,N=NULL,af.dm=NULL,spatial_coords=NULL,
                          coverage_thresh=0,k = 50){
  var_pairs = combn(voi,2,simplify = FALSE)
  LISI_score = list()
  for(pair in var_pairs){
    print(pair)
    var_1 = pair[1]
    var_2 = pair[2]
    var_1_spot_to_keep = intersect(colnames(N)[(N[var_1,] > coverage_thresh) &
                                                 (af.dm[var_1,]>0)],
                                   rownames(spatial_coords))
    var_2_spot_to_keep = intersect(colnames(N)[(N[var_2,] > coverage_thresh) &
                                                 (af.dm[var_2,]>0)],
                                   rownames(spatial_coords))
    knn_res <- get.knn(spatial_coords, k = k)
    neighbor_list <- apply(knn_res$nn.index, 1, function(idxs) {
      rownames(spatial_coords)[idxs]
    })
    colnames(neighbor_list) <- rownames(spatial_coords)
    C_s_12_list = c()
    for(s in var_1_spot_to_keep){
      s_nb=neighbor_list[,s]
      if(any(af.dm[var_1,s_nb]>0) & any(af.dm[var_2,s_nb]>0)){
        c_s=2
      }else if (any(af.dm[var_1,s_nb]>0) & all(af.dm[var_2,s_nb]==0)){
        c_s = 1
      }else if (any(af.dm[var_2,s_nb]>0) & all(af.dm[var_1,s_nb]==0)){
        c_s = 1
      }else{
        c_s=0
      }
      C_s_12_list <- append(C_s_12_list,c_s)
    }
    C_s_21_list = c()
    for(s in var_2_spot_to_keep){
      s_nb=neighbor_list[,s]
      if(any(af.dm[var_1,s_nb]>0) & any(af.dm[var_2,s_nb]>0)){
        c_s=2
      }else if (any(af.dm[var_1,s_nb]>0) & all(af.dm[var_2,s_nb]==0)){
        c_s = 1
      }else if (any(af.dm[var_2,s_nb]>0) & all(af.dm[var_1,s_nb]==0)){
        c_s = 1
      }else{
        c_s=0
      }
      C_s_21_list <- append(C_s_21_list,c_s)
    }
    avg_Cs = mean(c(C_s_12_list,C_s_21_list),na.rm=T)
    LISI_metric = c(mean(C_s_12_list),mean(C_s_21_list),avg_Cs)
    names(LISI_metric) = c("var_12","var_21","avg")
    LISI_score[[paste(pair,collapse=",")]] <- LISI_metric
  }
  return(LISI_score)
}


#' Spatial Self-Mixing LISI
#'
#' Computes a per-variant spatial mixing score using k-nearest neighbors (kNN).
#' For each variant, spots passing a coverage threshold are split into
#' variant-positive (VAF > \code{vaf_thresh}) and variant-negative (the remaining
#' eligible spots). For each variant-positive spot, the function checks whether
#' its kNN neighborhood contains both positive and negative spots.
#'
#' The per-spot score is:
#' \itemize{
#'   \item 2 if both positive and negative neighbors are present
#'   \item 1 if only one group is present
#'   \item 0 if no eligible neighbors are present
#' }
#' The reported score for a variant is the mean over variant-positive spots.
#'
#' @param voi Character vector of variant IDs to analyze.
#' @param N Numeric matrix of coverage counts (variants x spots).
#' @param af.dm Numeric matrix of allele frequencies (variants x spots).
#' @param spatial_coords Numeric matrix/data.frame of spatial coordinates
#'   with rownames as spot IDs (rows = spots, cols = spatial dimensions).
#' @param coverage_thresh Minimum coverage required to include a spot.
#'   Default is 0.
#' @param vaf_thresh Minimum allele frequency threshold to define variant-positive
#'   spots. Default is 0.25.
#' @param k Number of nearest neighbors used to build the kNN graph. Default is 50.
#'
#' @return Named numeric vector of LISI-like scores for each variant in \code{voi}.
#'   Variants with no positive or no negative spots among eligible spots return \code{NA}.
#'
#' @importFrom FNN get.knn
#' @export
LISI_self <- function(voi,N,af.dm,spatial_coords,
                      coverage_thresh = 0,vaf_thresh = 0.25,k = 50) {
  LISI_scores <- c()
  # keep only spots that exist in spatial coords 
  # spots_all <- rownames(spatial_coords)
  # 
  # N <- N[, intersect(colnames(N), spots_all), drop = FALSE]
  # af.dm <- af.dm[, intersect(colnames(af.dm), spots_all), drop = FALSE]
  
  # kNN on all spatial coords 
  knn_res <- get.knn(spatial_coords, k = k)
  neighbor_list <- apply(knn_res$nn.index, 1, function(idxs) {
    rownames(spatial_coords)[idxs]
  })
  colnames(neighbor_list) <- rownames(spatial_coords)
  
  for (var in voi) {
    # candidate spots that pass coverage and exist in coords
    eligible <- intersect(
      colnames(N)[N[var, ] > coverage_thresh],
      rownames(spatial_coords)
    )
    
    # define variant-positive and non-variant groups among eligible
    var_pos <- intersect(eligible, colnames(af.dm)[af.dm[var, ] > vaf_thresh])
    var_neg <- setdiff(eligible, var_pos)  # includes 0 <= VAF <= thresh
    
    # handle edge cases
    if (length(var_pos) == 0 || length(var_neg) == 0) {
      LISI_scores <- append(LISI_scores,NA)
      next
    }
    
    C_s <- numeric(length(var_pos))
    names(C_s) <- var_pos
    
    for (s in var_pos) {
      s_nb <- neighbor_list[, s]
      s_nb <- intersect(s_nb, eligible)  # only count neighbors that are eligible
      
      has_pos <- any(s_nb %in% var_pos)
      has_neg <- any(s_nb %in% var_neg)
      
      C_s[s] <- if (has_pos && has_neg) 2 else if (has_pos || has_neg) 1 else 0
    }
    
    LISI_scores <- append(LISI_scores, mean(C_s, na.rm = TRUE))
  }
  names(LISI_scores)<- voi
  return(LISI_scores)
}



#' Global Moran's I
#'
#' Computes global Moran's I for each variant using allele frequency values
#' across spots, with spots filtered by a minimum coverage threshold.
#' Spatial neighbors are defined via a k-nearest neighbor (kNN) graph built
#' from \code{spatial_coords}, and weights are row-standardized (\code{style = "W"}).
#'
#' In addition to the variants in \code{voi}, the function computes a
#' "control" Moran's I by permuting (\code{sample}) the allele frequencies of
#' the first variant in \code{voi} while keeping coverage and coordinates fixed.
#'
#' @param voi Character vector of variant IDs to analyze. The first element is
#'   also used to generate the control.
#' @param N Numeric matrix of coverage counts (variants x spots).
#' @param af.dm Numeric matrix of allele frequencies (variants x spots).
#' @param spatial_coords Numeric matrix/data.frame of spatial coordinates
#'   with rownames as spot IDs (rows = spots, cols = spatial dimensions).
#' @param coverage_thresh Minimum coverage required to include a spot.
#'   Default is 1.
#' @param k Number of nearest neighbors used to build the kNN graph. Default is 50.
#' @param seed Integer seed for reproducible control permutation. Default is 42.
#'
#' @return A matrix with rows \code{c(voi, "control")} and columns:
#'   \code{mi_list} (Moran's I statistic) and \code{p_list} (p-value).
#'
#' @details
#' Uses \code{spdep::knearneigh}, \code{spdep::knn2nb}, \code{spdep::nb2listw}
#' and \code{spdep::moran.test} with \code{zero.policy = TRUE}.
#'
#' @importFrom spdep knearneigh knn2nb nb2listw moran.test
#' @export
moran_I_knn <- function(voi=NULL,N=NULL,af.dm=NULL,spatial_coords=NULL,
                        coverage_thresh=1,k=50,seed=42){
  mi_list = c()
  p_list = c()
  for(var in voi){
    # filter spots by coverage & nonzero VAF
    keep <- which((N[var, ] > coverage_thresh))#& 
    #(af.dm[var, ] > 0))
    spot_ids <- intersect(colnames(N)[keep],rownames(spatial_coords))
    coords    <- spatial_coords[spot_ids, , drop=FALSE]
    values    <- af.dm[var, spot_ids]
    
    # build k‐NN neighbor list
    knn       <- spdep::knearneigh(coords, k = k)
    nb        <- spdep::knn2nb(knn)
    
    # convert to spatial weights (row‐standardized)
    lw        <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # compute global Moran's I
    mi_test   <- spdep::moran.test(values, lw, zero.policy = TRUE)
    mi_stat   <- mi_test$estimate[["Moran I statistic"]]
    p_value   <- mi_test$p.value
    
    cat(sprintf("Variant %s: Moran's I = %.4f (p = %.3g)\n", 
                var, mi_stat, p_value))
    mi_list <- append(mi_list,mi_stat)
    p_list <- append(p_list,p_value)
  } 
  # compare with control
  base_var = voi[1]
  set.seed(seed)
  ctrl_af <- sample(af.dm[base_var, ]) 
  names(ctrl_af) <- colnames(af.dm)
  ctrl_N  <- N[base_var, ]
  keep <- which((ctrl_N> coverage_thresh))
  spot_ids <- intersect(names(ctrl_N)[keep],rownames(spatial_coords))
  coords    <- spatial_coords[spot_ids, , drop=FALSE]
  values    <- ctrl_af[spot_ids]
  
  # build k‐NN neighbor list
  knn       <- spdep::knearneigh(coords, k = k)
  nb        <- spdep::knn2nb(knn)
  
  # convert to spatial weights (row‐standardized)
  lw        <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # compute global Moran's I
  mi_test   <- spdep::moran.test(values, lw, zero.policy = TRUE)
  mi_stat   <- mi_test$estimate[["Moran I statistic"]]
  p_value   <- mi_test$p.value
  
  cat(sprintf("Control Variant: Moran's I = %.4f (p = %.3g)\n", 
              mi_stat, p_value))
  
  mi_list <- append(mi_list,mi_stat)
  p_list <- append(p_list,p_value)
  res = cbind(mi_list=mi_list,p_list=p_list);rownames(res) = c(voi,"control")
  
  return(res)
}