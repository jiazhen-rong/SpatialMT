## create synthetic location 
# Function to create synthetic grid
simulate_clone_location <- function(grid_size = 30,clones,plot=T,
                                    clone_colors=NULL,save=F,save_path=NULL,
                                    height=10,width=10) {
  x_coord <- rep(1:grid_size, each = grid_size)
  y_coord <- rep(1:grid_size, grid_size)
  plot_df <- data.frame(x_coord = x_coord, y_coord = y_coord, label = "Other")
  
  for (clone in clones) {
    x_start <- clone$x_start
    y_start <- clone$y_start
    width <- clone$width
    height <- clone$height
    label <- paste0("Clone_", clone$id)
    
    plot_df[(plot_df$x_coord %in% x_start:(x_start + width - 1)) &
              (plot_df$y_coord %in% y_start:(y_start + height - 1)), "label"] <- label
  }
  
  # plot the simulated clones in space
  p=ggplot(plot_df, aes(x = x_coord, y = y_coord, color = label)) +
    geom_point(size=2) + theme_bw() +
    scale_color_manual(values =clone_colors) +
    coord_fixed()+
    ggtitle("Simulated Clones") 
  if(plot){
    print(p)
  }
  if(save){
    saveRDS(plot_df,paste0(save_path,"/spatial_coords.rds"))
    pdf(paste0(save_path,"/spatial_locations.pdf"),height=height,width=width)
    print(p)
    dev.off()
  }
  
  return(plot_df)
}


# Function to simulate VAF (X) matrix
simulate_vaf <- function(plot_df, grid_size = 30, var_num = 5, noise_var_num=15,
                         N = 10, clone_num = 3, clone_size = 50, 
                         vaf = 0.5, background_vaf = 0.1,
                         save=F,save_path=NULL) {
  vaf_mtx <- matrix(0, nrow = (clone_num * var_num+noise_var_num), ncol = grid_size^2)
  rownames(vaf_mtx) <- paste0("Var_", 1:((clone_num) * var_num+noise_var_num))
  # simulate real signals in each clone
  for (j in 1:clone_num) {
    for (i in ((j - 1) * var_num + 1):(j * var_num)) {
      vaf_mtx[i, plot_df$label == paste0("Clone_", j)] <- rbinom(sum(plot_df$label == paste0("Clone_", j)), size = N, prob = vaf)
    }
  }
  # simulate noisy signal for all cells/spots
  for (i in ((clone_num * var_num) + 1):(((clone_num * var_num) + noise_var_num))) {
    vaf_mtx[i, ] <- rbinom(n = grid_size^2, size = N, prob = background_vaf)
  }
  if(save){
    saveRDS(vaf_mtx,paste0(save_path,"/vaf_",vaf,"_N_",N,".rds")) 
  }
  return(vaf_mtx)
}

plot_vaf<- function(vaf_mtx,plot_df,vaf=NA,N=NA,clone_num=3,plot=T,save=F,save_path=NULL,
                    width=10,height=10){
  clone_idx = unlist(lapply(1:clone_num,function(j){
    which(plot_df$label==paste0("Clone_",j))}))
  pdf(paste0(save_path,"/vaf_",vaf,"_N_",N,"_Spot_Var_Heatmap.pdf"),
      width=width,height=height)
  # plot clone cells x all variants 
  pheatmap(vaf_mtx[,clone_idx],cluster_rows=F,cluster_cols =F,
           main="Clone Spots")
  # plot all cells x all variants
  pheatmap(vaf_mtx,cluster_rows=F,main="All Spots")
  dev.off()
}


library(truncnorm)
library(grid)
library(gridExtra)
# Simulate Cell-type Ratio - spot x celltype
simulate_celltype_ratio <- function(plot_df, 
                                    celltypes=NULL,
                                    celltype_ratio=0.75,
                                    sd=0.1,
                                    clone_num=3,
                                    grid_size=30,
                                    save=F,save_path=NULL,plot=T,
                                    distribution=c("TruncN","Beta"),
                                    cell_range=c("Clone","Custom"),
                                    celltype_idx=NULL) {
  
  Ws <- matrix(0, nrow = grid_size^2,ncol=length(celltypes))
  colnames(Ws) = celltypes
  
  # For Clone_j cells/spots
  for(j in 1:length(celltypes)){
    if(cell_range == "Clone"){
      if(j != length(celltypes)){
        idx_clone_j <- which(plot_df[,"label"] == paste0("Clone_",j))
      }else{
        idx_clone_j <- which(plot_df[,"label"] == "Other")
      }
    }else if(cell_range == "Custom"){
      idx_clone_j = celltype_idx[[j]]
    }
    if(distribution == "TruncN"){
      result <-  rtruncnorm(length(idx_clone_j), a = 0, b = 1, 
                            mean = celltype_ratio, sd = sd)
    }else if(distribution == "Beta"){ # from beta distribution
      mu=celltype_ratio
      a=mu*(mu*(1-mu)/sd^2-1)
      b=(1-mu)*(mu*(1-mu)/sd^2-1)
      result <- rbeta(n=length(idx_clone_j),shape1=a,shape2=b)
    }
    # assign cell-type j; then split by rest 
    Ws[idx_clone_j,j] <- result
    #Ws[idx_clone_j,-j] <- matrix((1 - result) / (length(celltypes) - 1),
    #                             nrow = length(idx_clone_j),
    #                             ncol = length(celltypes) - 1)
    
    # fill rest with normal cells
    Ws[idx_clone_j,length(celltypes)] = 1-result
  }
  if(save){
    saveRDS(Ws,paste0(save_path,"/sim_celltype_ratio_",
                      celltype_ratio,"_",cell_range,".rds"))
  }
  if(plot){
    temp_df = as.data.frame(cbind(plot_df,Ws))
    g_list = list()
    for(celltype in celltypes){
      p <- ggplot(temp_df,aes(x=x_coord,y=y_coord,color=.data[[celltype]]))+
        geom_point() +
        theme_minimal() + 
        scale_color_viridis_c(option = "D", limits = c(0, 1)) +
        #scale_color_gradient(low = "lightblue", high = "darkblue") +
        theme_bw() +
        labs(x = "X", y = "Y", title = paste0("Simulated Celltype ",celltype)) +
        theme(plot.title = element_text(hjust = 0.5))
      g_list[[celltype]] <- p
    }
    pdf(paste0(save_path,"/simulated_celltype_spatial_ratio_",
               celltype_ratio,"_",cell_range,".pdf"),width=10,height=10)
    grid.draw(grid.arrange(grobs=g_list,ncol=ceiling(sqrt(length(celltypes)))))
    dev.off()
  }
  return(Ws)
}

library(dplyr)
# for purly synthetic data
calc_metric<- function(pred_var,ground_truth,clone_num){
  class_labels=c(paste0("Clone_",1:clone_num),"Noise")
  # Convert to factor for consistency
  pred_var <- factor(pred_var,levels=class_labels)
  ground_truth <- factor(ground_truth,levels=class_labels)
  
  # Define positive and negative (Noise) classes
  positive_classes <- paste0("Clone_",1:clone_num)
  
  # Create binary vectors: 1 for positive, 0 for negative
  true_labels <- ifelse(ground_truth %in% positive_classes, 1, 0)
  pred_labels <- ifelse(pred_var %in% positive_classes, 1, 0)
  
  # Compute confusion matrix components
  TP <- sum(pred_labels == 1 & true_labels == 1)  # True Positives
  TN <- sum(pred_labels == 0 & true_labels == 0)  # True Negatives
  FP <- sum(pred_labels == 1 & true_labels == 0)  # False Positives (Type I Error)
  FN <- sum(pred_labels == 0 & true_labels == 1)  # False Negatives (Type II Error)
  
  # Compute metrics
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), NA)
  type_I_error <- ifelse((FP + TN) > 0, FP / (FP + TN), NA)  # False Positive Rate
  type_II_error <- ifelse((FN + TP) > 0, FN / (FN + TP), NA)  # False Negative Rate
  
  # Compute confusion matrix (for multi-class)
  cm <- table(ground_truth, pred_var)
  
  # Initialize storage
  class_accuracies <- numeric(clone_num+1)
  class_f1_scores <- numeric(clone_num+1)
  class_support <- numeric(clone_num+1)  # Number of samples per class
  
  # Compute per-class metrics
  for (i in seq_along(class_labels)){
    class <- class_labels[i]
    
    TP <- cm[class, class] %>% ifelse(is.na(.), 0, .)  # True Positives
    FN <- sum(cm[class, ]) - TP                      # False Negatives
    FP <- sum(cm[, class]) - TP                      # False Positives
    TN <- sum(cm) - (TP + FN + FP)                   # True Negatives
    
    # Class-wise accuracy
    class_accuracies[i] <- ifelse(sum(cm[class, ]) > 0, TP / sum(cm[class, ]), NA)
    
    # Precision, Recall, and F1-score
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
    
    class_f1_scores[i] <- f1_score
    class_support[i] <- sum(cm[class, ])  # Number of samples in true class
  }
  
  # Overall metrics
  overall_accuracy <- sum(diag(cm)) / sum(cm)  # (Sum of diagonal / total samples)
  avg_class_accuracy <- mean(class_accuracies)  # Mean per-class accuracy
  weighted_f1 <- sum(class_f1_scores * class_support) / sum(class_support)  # Weighted F1-score
  
  # Return as a named list
  return(list(
    overall_accuracy = overall_accuracy,
    avg_class_accuracy = avg_class_accuracy,
    weighted_f1_score = weighted_f1,
    #class_accuracies = setNames(class_accuracies, class_labels),
    type_I_error=type_I_error,
    type_II_error=type_II_error
  ))
}


# for data synthesized from real data
calc_metric_realsim <- function(pred_var,ground_truth,celltypes,sig_celltypes=c("BE")){
  #class_labels=c(paste0("Clone_",1:clone_num),"Noise")
  class_labels=c(celltypes,"Noise")
  # Convert to factor for consistency
  pred_var <- factor(pred_var,levels=class_labels)
  ground_truth <- factor(ground_truth,levels=class_labels)
  
  # Define positive and negative (Noise) classes
  positive_classes <- sig_celltypes
  
  # Create binary vectors: 1 for positive, 0 for negative
  true_labels <- ifelse(ground_truth %in% positive_classes, 1, 0)
  pred_labels <- ifelse(pred_var %in% positive_classes, 1, 0)
  
  # Compute confusion matrix components
  TP <- sum(pred_labels == 1 & true_labels == 1)  # True Positives
  TN <- sum(pred_labels == 0 & true_labels == 0)  # True Negatives
  FP <- sum(pred_labels == 1 & true_labels == 0)  # False Positives (Type I Error)
  FN <- sum(pred_labels == 0 & true_labels == 1)  # False Negatives (Type II Error)
  
  # Compute metrics
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), NA)
  type_I_error <- ifelse((FP + TN) > 0, FP / (FP + TN), NA)  # False Positive Rate
  type_II_error <- ifelse((FN + TP) > 0, FN / (FN + TP), NA)  # False Negative Rate
  
  # Compute confusion matrix (for multi-class)
  cm <- table(ground_truth, pred_var)
  
  # Initialize storage
  class_accuracies <- numeric(clone_num+1)
  class_f1_scores <- numeric(clone_num+1)
  class_support <- numeric(clone_num+1)  # Number of samples per class
  
  # Compute per-class metrics
  for (i in seq_along(class_labels)){
    class <- class_labels[i]
    
    TP <- cm[class, class] %>% ifelse(is.na(.), 0, .)  # True Positives
    FN <- sum(cm[class, ]) - TP                      # False Negatives
    FP <- sum(cm[, class]) - TP                      # False Positives
    TN <- sum(cm) - (TP + FN + FP)                   # True Negatives
    
    # Class-wise accuracy
    class_accuracies[i] <- ifelse(sum(cm[class, ]) > 0, TP / sum(cm[class, ]), NA)
    
    # Precision, Recall, and F1-score
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
    
    class_f1_scores[i] <- f1_score
    class_support[i] <- sum(cm[class, ])  # Number of samples in true class
  }
  
  # Overall metrics
  overall_accuracy <- sum(diag(cm)) / sum(cm)  # (Sum of diagonal / total samples)
  avg_class_accuracy <- mean(class_accuracies* class_support,na.rm=T)/sum(class_support)   # Mean per-class accuracy
  weighted_f1 <- sum(class_f1_scores * class_support) / sum(class_support)  # Weighted F1-score
  
  # Return as a named list
  return(list(
    overall_accuracy = overall_accuracy,
    avg_class_accuracy = avg_class_accuracy,
    weighted_f1_score = weighted_f1,
    type_I_error=type_I_error,
    type_II_error=type_II_error
  ))
}


# Function to simulate VAF (X) matrix
simulate_vaf_from_read_data <- function(temp_plot_df,vaf_mu=0.5,vaf_sd=0.2,
                                        noise_mu=0.1,noise_sd=0.2,
                                        N_voi,clone_size,coherence,
                                        var_num=1,clone_num=3,noise_var_num=1,
                                        plot=F) {
  clone_col = paste0("CloneSize_",clone_size,"_Coherence_",coherence)
  clone_id = temp_plot_df[,clone_col]
  
  # Check clone assignment column exists
  if (!clone_col %in% colnames(temp_plot_df)) {
    stop(paste("Column", clone_col, "not found in temp_plot_df"))
  }
  
  # Extract clone membership per spot
  n_spots <- nrow(temp_plot_df)
  n_mean <- ceiling(apply(N_voi,2,mean))#apply(N_voi,2,mean)
  # Simulate read depth per spot: ceiling of average depth per variant
  read_depths <- rpois(n_spots, lambda = n_mean)
  read_depths[read_depths == 0] <- 1  # avoid zero
  
  vaf_mtx <- matrix(0, nrow = (clone_num * var_num+noise_var_num), ncol = n_spots)
  rownames(vaf_mtx) <- paste0("Var_", 1:((clone_num) * var_num+noise_var_num))
  colnames(vaf_mtx) <- rownames(temp_plot_df)
  
  # Simulate variant reads per spot and per variant
  # For each variant
  for (j in 1:(clone_num * var_num+noise_var_num)) {
    if(j <= clone_num * var_num){
      k = ceiling(j/var_num)
      idx_clone_j <- which(clone_id == paste0(k))
        
      # simulate vaf and X from beta distributions
      mu=vaf_mu;sd=vaf_sd
      a=mu*(mu*(1-mu)/sd^2-1)
      b=(1-mu)*(mu*(1-mu)/sd^2-1)
      sim_vaf <- rbeta(n=length(idx_clone_j),shape1=a,shape2=b)
      
      # Simulate observed alt read count from binomial
      X_j <- rbinom(length(idx_clone_j), size = read_depths[idx_clone_j], prob = sim_vaf)
      vaf_mtx[j, idx_clone_j] <- X_j
    }else{ # noise in all spots
      idx_clone_j <- 1:n_spots #which(is.na(clone_id))
      mu=noise_mu;sd=noise_sd
      a=mu*(mu*(1-mu)/sd^2-1)
      b=(1-mu)*(mu*(1-mu)/sd^2-1)
      sim_vaf <- rbeta(n=length(idx_clone_j),shape1=a,shape2=b)
      # Simulate observed alt read count from binomial
      X_j <- rbinom(length(idx_clone_j), size = read_depths[idx_clone_j], prob = sim_vaf)
      vaf_mtx[j, idx_clone_j] <- X_j
    }
  }
  if(plot){
    order_id = unlist(lapply(1:clone_num,function(clone_l){which(clone_id == paste0(clone_l))}))
    pheatmap(vaf_mtx[,order_id],show_colnames = F,cluster_row=F,cluster_col=F)
    pheatmap(vaf_mtx[,order_id],show_colnames = F,cluster_row=T,cluster_col=T)
    
  }
  return(list(X_mtx=vaf_mtx,N = read_depths))
}

calc_spot_level_metrics <- function(pred_sig_spots, real_sig_spots) {
  # Ensure inputs are logical vectors
  pred_sig_spots <- as.logical(pred_sig_spots)
  real_sig_spots <- as.logical(real_sig_spots)
  
  # Confusion matrix components
  TP <- sum(pred_sig_spots & real_sig_spots)
  TN <- sum(!pred_sig_spots & !real_sig_spots)
  FP <- sum(pred_sig_spots & !real_sig_spots)
  FN <- sum(!pred_sig_spots & real_sig_spots)
  
  # Calculate metrics
  accuracy    <- (TP + TN) / length(pred_sig_spots)
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall      <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1_score    <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 
    2 * precision * recall / (precision + recall) else NA
  type_I_error <- if ((FP + TN) > 0) FP / (FP + TN) else NA
  type_II_error <- if ((TP + FN) > 0) FN / (TP + FN) else NA
  
  # Return the metrics as a list
  return(list(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    type_I_error = type_I_error,
    type_II_error = type_II_error
  ))
}
spot_fdr <- function(p_df) {
  # Flatten the entire p-value matrix into a vector
  p_vector <- as.vector(as.matrix(p_df))
  # Adjust p-values globally using the BH method
  p_adj <- p.adjust(p_vector, method = "BH")
  # Reshape the adjusted p-values back to the original matrix dimensions
  p_adj_mat <- matrix(p_adj, nrow = nrow(p_df), ncol = ncol(p_df))
  rownames(p_adj_mat) <- rownames(p_df)
  colnames(p_adj_mat) <- colnames(p_df)
  as.data.frame(p_adj_mat)
}


