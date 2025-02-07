#' Regression Test for Identifying Significant Mutations and Carrier Celltypes.
#'
#' `celltype_test()` implements the statistical tests for identifying significant celltypes given spatial decomposition and varaint information.
#'
#' @param celltypes A string vector containing all unique celltypes.
#' @param voi A string vector containing names of variant of interest. Usually a subset of variants that went through basic filtering.
#' @param N A MT variant x cell total coverage matrix.
#' @param X A MT variant x cell alternative allele count matrix. Optional if vaf matrix is provided.
#' @param vaf A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as af.dm matrix.
#' @param Ws A cell x celltype weight matrix. Representing spatial decomposition results from tools like RCTD.
#' @param spatial_coords A cell x spatial coordinates matrix. Representing the actual spot/cell location in space.
#' @param test_type Type of Celltype test. Options from c("linear","weighted").
#' @param permute_num A numeric number representing permuation number for each variant.
#' @param save_path path for output results
#' @param figure_height output pdf figure height
#' @param figure_width output pdf figure width
#'
#' @param ... Other options used to control matching behavior between duplicate
#'   strings. Passed on to [stringi::stri_opts_collator()].
#' @returns A list object containing fitted parameters (intercept, coefficient and p-values).
#' @import Matrix
#'
#' @examples
#' celltype_test(celltypes=c("BE","DYSP","Normal"), voi=voi_list,
#' N=N,vaf=af.dm, test_type="linear",permute_num=1000)
#'
#' @export
#'
celltype_test <- function(celltypes=NULL,voi=NULL,N=NULL,vaf=NULL,X=NULL,Ws=NULL,spatial_coords=NULL,
                          test_type=c("linear","weighted"),permute_num=1000,plot=F,save_path=NULL,
                          figure_height=10,figure_width=10){
  # check for input and conditions
  if(test_type == "linear"){
    message("Performing Linear Rgression ...")
  }else if(test_type == "weighted"){
    message("Performing Linear Rgression weighted by coverage...")
  }else{
    error("Not a valid test. Please choose from documented tests.")
  }

  if(is.null(vaf)){
    if(!is.null(N)){
      vaf_j = X[j,]/N[j,]
    }else{
      error("Missing alternative allelle count matrix N or frequency matrix vaf.")
    }
  }
  # check if celltypes aligned with weight matrix
  if(is.null(celltypes)){
    message("Missing celltypes...Set as the column names of celltype weight matrix Ws.")
    celltypes = colnames(Ws)
  }else if(any(celltypes != colnames(Ws))){
    error("Celltypes and the column names of celltype weight matrix Ws not matched...Please check the variables celltypes or Ws.")
  }
  if(is.null(save_path)){
    save_path = "output/"
  }else{
    dir.create(save_path)
  }

  # create plotting data frame
  if(plot==T){
    intersect_bc=intersect(intersect(colnames(N),rownames(Ws)),rownames(spatial_coords))
    plot_df=cbind(spatial_coords[intersect_bc,],
                  Ws[intersect_bc,],
                  as.data.frame(t(vaf[,intersect_bc])))
    # add in coverage of variant of interest
    cov = N[,intersect_bc];rownames(cov) =  paste0(rownames(N),"_cov")
    plot_df = cbind(plot_df,t(cov))
  }

  # for each variant j, fit a Xs/Ns vs Ws Linear regression
  intercept_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  coef_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  pval_df=Matrix(NA,nrow=length(voi),ncol=length(celltypes))

  if(plot==T){
    pdf(paste0(save_path,"/",test_type,"_regression_anova_p_VAF_vs_celltype_plot.pdf"),
        height=figure_height,width=figure_width)
  }
  for(j in 1:length(voi)){
    var = voi[j]
    print(paste0(j,": ", var))
    N_j = N[j,]
    vaf_j = vaf[j,]


    for(k in 1:length(celltypes)){ # for each celltype

      data = data.frame(vaf_j=vaf_j,N_j,W_sk=Ws[,k])
      cell_ratio_bins = seq(from=0,to=1,by=0.05)

      if(test_type == "linear"){
        fit=lm(vaf_j ~ W_sk,data=data)# weights equals to inverse of variance (1/sigma_i^2)
      }else if(test_type == "weighted"){
        fit=lm(vaf_j ~ W_sk,data=data, weights=sqrt(N_j))# weights equals to inverse of variance (1/sigma_i^2)
      }
      res=summary(fit)
      intercept_df[j,k]=res$coefficients[1,1]
      coef_df[j,k]= res$coefficients[2,1]

      if(test_type == "linear"){
        pval_df[j,k] = res$coefficients[2,4]
      }else if(test_type == "weighted"){
        # doing permutations and see where the data locates
        set.seed(42)
        coef_list <- c()
        #  decouple VAF with Ws
        for(l in 1:permute_num){
          idx=sample(1:length(vaf_j),size=length(vaf_j))
          vaf_shuffled = vaf_j[idx]
          #vaf_shuffled = pmin(X_shuffled/(N_j+1e-6),1) # cap max VAF to be by 1
          fit=lm(vaf_shuffled ~ Ws[,k],weights=sqrt(N_j[idx]))# weights equals to inverse of variance (1/sigma_i^2)
          res_shuf=summary(fit)
          coef_list <- append(coef_list,res_shuf$coefficients[2,1])
        }
        pval_df[j,k] = 1- sum(res$coefficients[2,1] >= coef_list)/length(coef_list)
      }
    }
    # for variant j, plot out the fitted values on the diagnostic plots
    if(plot==T){
      intercept=intercept_df[j,]
      coef=coef_df[j,]
      pval=pval_df[j,]
      plot_vaf_cellprop(j,af.dm,Ws,plot_df,intercept,coef,pval)
    }
  }
  if(plot==T){
    dev.off()
    message("Diagnostic plots generated.")
  }
  rownames(intercept_df)=voi
  rownames(coef_df)=voi
  rownames(pval_df)=voi
  colnames(intercept_df) = colnames(Ws)
  colnames(coef_df) = colnames(Ws)
  colnames(pval_df) = colnames(Ws)

  return(res=list(intercept=intercept_df,coef=coef_df,pval=pval_df))
}


#' Power Analysis for Identifying Negative Mutations and Carrier Celltypes.
#' Simulate alternative allele count values Xi through Binomial distribution Xi ~ Binom(Ni,beta_0 + beta_1 * Wik)
#'
#' @param beta A numeric value. Equivalent as beta_1 in the simulation.
#' @param Nj A total coverage counts for a given variant j
#' @param X A MT variant x cell alternative allele count matrix. Optional if vaf matrix is provided.
#' @param vaf A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as af.dm matrix.
#' @param Ws A cell x celltype weight matrix. Representing spatial decomposition results from tools like RCTD.
#' @param spatial_coords A cell x spatial coordinates matrix. Representing the actual spot/cell location in space.
#' @param test_type Type of Celltype test. Options from c("linear","weighted").
#' @param permute_num A numeric number representing permuation number for each variant.
#'
calc_power <- function(beta,Nj,Wj,vaf_j,alpha,n_sim,verbose=T,report_frac=3){
  rejections = 0
  p_vals = c()
  for(i in 1:n_sim){
    if(verbose){
      if(i%%(round(n_sim/report_frac))){
       print(paste0("sim:",i,"/",n_sim))
      }
    }
    X_sim = rbinom(n=length(Nj),size=Nj,prob=beta*Wj)
    data = data.frame(vaf_j=X_sim/(Nj+0.001),Nj,Wj=Wj)
    fit=lm(vaf_j ~ Wj,data=data)
    res=summary(fit)
    p_vals = append(p_vals,res$coefficients[2,4])
  }
  rejections  = sum(p_vals <= alpha,na.rm=T)
  power = rejections/(n_sim)
  return(power)
}

#' Power Analysis for Identifying Negative Mutations and Carrier Celltypes.
#' for all variants of interest & all celltypes
#' set permute_num=20 for now, shall be 1000
#' @export
power_analysis_all <- function(voi=NULL,celltypes=NULL,Ws=NULL,N=NULL,vaf=NULL,X=NULL,
                               sample_num=100,alpha=0.05,n_sim=10,beta_threshold=0.5,plot=T,save_path=NULL,
                               figure_height=8,figure_width=8){
  if(is.null(save_path)){
    save_path = "output/"
    message("Setting output folder to be 'outoput'")
    dir.create(save_path)

  }else{
    dir.create(save_path)
  }
  message("Running Power Analysis for Negative Celltypes and Mutations")
  intersect_bc = intersect(colnames(N),rownames(Ws))
  beta_list = c()

  if(plot){
    pdf( paste0(save_path,"/example_power_analysis.pdf"),height=figure_height,width=figure_width)
  }
  for(j in 1:length(voi)){
    var = voi[j];
    print(paste0("voi ",j,":",var))
    N_j = N[j,intersect_bc]
    vaf_j = vaf[j,intersect_bc]
    if(plot){
      # Dynamically set the plotting grid to fit the number of cell types
      par(mfrow = c(ceiling(sqrt(length(celltypes))), ceiling(sqrt(length(celltypes)))))
    }
    for(k in 1:length(celltypes)){
      sample_idx = sample(intersect_bc,sample_num)
      Nj = N_j[sample_idx] # sampled coverage
      Wj = Ws[sample_idx,k] # celltype ratio
      beta_grid <- seq(0.1,1,by=0.2)

      power_results = sapply(beta_grid,function(beta){
        #print(beta)
        calc_power(beta,Nj,Wj,vaf_j,alpha,n_sim,verbose=F)
      })
      names(power_results) = beta_grid
      beta_list[paste0(voi[j],"_",celltypes[k])] <- power_results

      if(plot){
        plot(beta_grid,power_results,type="l",col="blue",lwd=2,
             xlab="Beta",ylab="Power",main=paste0(var,", ",celltypes[k]),
             ylim=c(-0.1,1),xlim=c(0,1)
        )
        points(c(0.1,0.5,1),power_results[as.character(c(0.1,0.5,1))],pch=16)
        abline(h=0.8,col="orange",lty=2)
      }
    }
  }
  if(plot){
    dev.off()
  }
  return(beta_list)

}


# Grid based test - Scan statistics
#' Perform Celltype Test Within Spatial Grids
#'
#' This function partitions the spatial coordinates into a 20x20 grid and performs
#' the linear regression test within each grid to identify significant cell types.
#'
#' @param celltypes A string vector of unique cell types.
#' @param voi A string vector of variants of interest.
#' @param N A matrix of total coverage counts (MT variant x cell).
#' @param vaf A matrix of variant allele frequencies (MT variant x cell).
#' @param Ws A matrix of cell type weights (cell x celltype).
#' @param spatial_coords A matrix of spatial coordinates (cell x 2).
#' @param test_type A string specifying test type: "linear" or "weighted".
#' @param permute_num Number of permutations for weighted test.
#' @param grid_size Number of bins per axis (default = 20 for 20x20 grid).
#'
#' @return A nested list of p-values for each grid and variant.
#' @export
celltype_test_grid <- function(celltypes, voi, N, vaf, Ws, spatial_coords,
                               test_type = c("linear", "weighted"),
                               permute_num = 1000, grid_size = 20,verbose=F) {

  # Ensure spatial_coords is a matrix/dataframe
  spatial_coords <- as.data.frame(spatial_coords)
  colnames(spatial_coords) <- c("X", "Y") # Ensure correct column names

  # Create Grid Labels
  spatial_coords$X_bin <- cut(spatial_coords$X, breaks = grid_size, labels = FALSE)
  spatial_coords$Y_bin <- cut(spatial_coords$Y, breaks = grid_size, labels = FALSE)

  # Create Unique Grid IDs
  spatial_coords$Grid_ID <- paste0("Grid_", spatial_coords$X_bin, "_", spatial_coords$Y_bin)

  # Find common cells across datasets
  intersect_bc <- intersect(intersect(colnames(N), rownames(Ws)), rownames(spatial_coords))

  # Merge Spatial Data with Weights and VAF
  plot_df <- cbind(spatial_coords[intersect_bc, ],
                   Ws[intersect_bc, ],
                   as.data.frame(t(vaf[, intersect_bc])))

  # Add Coverage Data
  cov <- N[, intersect_bc]
  rownames(cov) <- paste0(rownames(N), "_cov")
  plot_df <- cbind(plot_df, t(cov))

  # Initialize Output
  results <- list()

  # Loop through each grid
  for (grid in unique(spatial_coords$Grid_ID)) {
    if(verbose){
      message(paste("Processing:", grid))
    }

    # Subset data for the current grid
    grid_cells <- spatial_coords$Grid_ID == grid
    if (sum(grid_cells) == 0) next  # Skip empty grids

    N_grid <- N[, grid_cells, drop = FALSE]
    vaf_grid <- vaf[, grid_cells, drop = FALSE]
    Ws_grid <- Ws[grid_cells, , drop = FALSE]

    # Initialize Data Storage
    intercept_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
    coef_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
    pval_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))

    rownames(intercept_df) <- voi
    colnames(intercept_df) <- celltypes
    rownames(coef_df) <- voi
    colnames(coef_df) <- celltypes
    rownames(pval_df) <- voi
    colnames(pval_df) <- celltypes

    # Loop through each variant
    for (j in 1:length(voi)) {
      var <- voi[j]
      if(verbose){
        message(paste("  Variant:", var))
      }

      N_j <- N_grid[j, ]
      vaf_j <- vaf_grid[j, ]

      for (k in 1:length(celltypes)) { # Loop through each cell type
        data <- data.frame(vaf_j = vaf_j, N_j = N_j, W_sk = Ws_grid[, k])

        # Check if there is enough variation in the predictor
        if (length(unique(data$W_sk)) == 1 || length(unique(data$vaf_j)) == 1) {
          message(paste("Skipping grid", grid, "for", var, "-", celltypes[k], "due to insufficient variation"))
          next
        }

        # Perform Linear Regression Safely
        if (test_type == "linear") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data), error = function(e) NULL)
        } else if (test_type == "weighted") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data, weights = sqrt(N_j)), error = function(e) NULL)
        }

        # Skip if fit fails
        if (is.null(fit)) {
          message(paste("Skipping grid", grid, "for", var, "-", celltypes[k], "due to regression failure"))
          next
        }

        # Store Results
        res <- summary(fit)

        # Ensure coefficients exist
        if (nrow(res$coefficients) < 2) {
          message(paste("Skipping grid", grid, "for", var, "-", celltypes[k], "due to missing coefficients"))
          next
        }

        intercept_df[j, k] <- res$coefficients[1, 1]
        coef_df[j, k] <- res$coefficients[2, 1]
        pval_df[j, k] <- res$coefficients[2, 4]
      }
    }

    # Store grid results
    results[[grid]] <- list(intercept = intercept_df, coef = coef_df, pval = pval_df)
  }

  return(results)
}

#' Plot P-Values from Grid-Based Regression
#'
#' This function plots the mean p-value per grid after running the `celltype_test_grid` function.
#'
#' @param results_grid A list containing p-values for each spatial grid.
#' @param grid_size Number of bins per axis (default = 20 for 20x20 grid).
#'
#' @return A `ggplot2` object showing p-value distribution across the spatial grid.
#' @export
#'
plot_pval_grid <- function(results_grid, var_name=NULL,grid_size = 20) {

  # Extract Grid IDs
  grid_ids <- names(results_grid)

  # Convert Grid IDs into Numeric X and Y Coordinates
  grid_data <- data.frame(
    Grid_ID = grid_ids,
    X_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 2)),
    Y_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 3)),
    Mean_Pval = NA
  )

  # Compute Mean P-value for Each Grid (Averaged Across Variants & Celltypes)
  for (i in seq_along(grid_ids)) {
    grid <- grid_ids[i]
    pvals <- results_grid[[grid]]$pval
    grid_data$Mean_Pval[i] <- mean(pvals, na.rm = TRUE)  # Average p-values
  }

  # Define Log-Scaled P-value
  grid_data$Log_Pval <- -log10(grid_data$Mean_Pval)  # Transform p-values

  # Avoid infinite values: Cap log-transformed values at a max threshold
  grid_data$Log_Pval[is.infinite(grid_data$Log_Pval)] <- max(grid_data$Log_Pval[!is.infinite(grid_data$Log_Pval)], na.rm = TRUE)

  # Define Custom Color Scale
  pval_palette <- scale_fill_gradientn(
    colors = c("red", "grey"),  # Red for significant, grey for non-significant
    values = scales::rescale(c(-log10(0.05), max(grid_data$Log_Pval, na.rm = TRUE))),  # Scale range
    limits = c(0, max(grid_data$Log_Pval, na.rm = TRUE)),  # Ensures proper mapping
    na.value = "white"  # Assigns white for missing values
  )

  # Plot Log-Scaled P-value Distribution Across the Grid
  pval_plot <- ggplot(grid_data, aes(x = X_bin, y = Y_bin, fill = Log_Pval)) +
    geom_tile(color = "white") +  # Grid Cells
    pval_palette +
    labs(title = paste0(var_name,"Log-Scaled P-value Distribution Across Spatial Grid"),
         x = "Grid X", y = "Grid Y", fill = "-log10(P-value)") +
    theme_minimal() +
    coord_fixed()  # Keep aspect ratio square

  pval_plot  # Display the plot


  return(pval_plot)
}

#' Plot P-Values from Grid-Based Regression for Each Cell Type
#'
#' This function plots the -log10 p-values for each cell type and variant
#' within the spatial grid, using a continuous color scale.
#'
#' @param results_grid A list containing p-values for each spatial grid.
#' @param grid_size Number of bins per axis (default = 20 for 20x20 grid).
#'
#' @return A `ggplot2` object showing p-value distributions for each cell type.
#' @export
#'
plot_pval_grid_per_celltype <- function(results_grid, grid_size = 20) {

  # Extract Grid IDs
  grid_ids <- names(results_grid)
  all_celltypes <- colnames(results_grid[[grid_ids[1]]]$pval)  # Extract cell types

  # Initialize List to Store Plots
  plot_list <- list()

  # Loop Through Each Cell Type and Create a Plot
  for (celltype in all_celltypes) {
    grid_data <- data.frame(
      Grid_ID = rep(grid_ids, each = length(results_grid[[1]]$pval[, celltype])),
      Variant = rep(rownames(results_grid[[1]]$pval), times = length(grid_ids)),
      X_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 2)),
      Y_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 3)),
      Pval = unlist(lapply(results_grid, function(res) res$pval[, celltype]))
    )

    # Drop missing values (NA p-values)
    grid_data <- grid_data %>% filter(!is.na(Pval))

    # If no valid data, skip this celltype
    if (nrow(grid_data) == 0) {
      message(paste("Skipping cell type", celltype, "due to missing data"))
      next
    }

    # Transform p-values: -log10 and cap at 4
    grid_data$Log_Pval <- -log10(grid_data$Pval)
    grid_data$Log_Pval[grid_data$Log_Pval > 4] <- 4  # Cap at 4

    # Define Continuous Color Scale
    pval_palette <- scale_fill_gradientn(
      colors = c("grey", "red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(0, -log10(0.05), 4)),  # Map p-values to colors
      limits = c(0, 4),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )

    # Generate Plots for Each Variant
    for (variant in unique(grid_data$Variant)) {
      df_variant <- grid_data %>% filter(Variant == variant)

      # Skip empty data
      if (nrow(df_variant) == 0) {
        message(paste("Skipping", variant, "for", celltype, "due to missing data"))
        next
      }

      p <- ggplot(df_variant, aes(x = X_bin, y = Y_bin, fill = Log_Pval)) +
        geom_tile(color = "white") +  # Grid Cells
        pval_palette +
        labs(title = paste(variant, "-", celltype),
             x = "Grid X", y = "Grid Y", fill = "-log10(P-value)") +
        theme_minimal() +
        coord_fixed()  # Keep aspect ratio square

      # Store plot in list only if it has data
      plot_list[[paste(variant, celltype, sep = "_")]] <- p
    }
  }

  # Ensure at least one plot exists
  if (length(plot_list) == 0) {
    message("No valid plots available")
    return(NULL)
  }

  # Arrange Plots in a Grid
  combined_plot <- plot_grid(plotlist = plot_list, ncol = ceiling(sqrt(length(all_celltypes))))  # Adjust ncol if needed

  return(combined_plot)
}
