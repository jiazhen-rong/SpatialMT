library(Seurat)
library(ggplot2)
library(Matrix)
library(spacexr)
library(gridExtra)

setwd("~/nzhanglab/project/jrong/mito_LT/")
source("scripts/our_model/example_data/210215_FunctionsGeneral.R") # from MAESTER paper
# spatial locations & varied clone size, coherence

# load spatial locations
temp_plot_df = readRDS("scripts/simulation_experiment/results/plot_df_withClones_A200S6.rds")

#### Load spatial MAESTER data
slide_name = "a200_s6"
# (1) mitochondrial data
maegatk.rse = readRDS("data/Sydney_Bracht/A200_s6/MAESTER/MAESTER_subsets/maegatk_mr1.rds") # master result object
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
Ws = readRDS("~/nzhanglab/project/jrong/mito_LT/results/Sydney_Bracht/sample_level/a200_s6_high_cov/RCTD/Major/RCTD_res.rds")
Ws = as.data.frame(normalize_weights(Ws@results$weights))

temp_plot_df = cbind(temp_plot_df,Ws[rownames(temp_plot_df),])
voi = read.csv("data/Sydney_Bracht/A200_s6/MAESTER/200_s6_voi_be25vaf.tsv",header=F)[,1]
N_voi = N[voi,]

# simulate vaf and N; do not change Ws.
vaf_list = c(0.1,0.25,0.5,0.75,0.9)
coherence_list = c(0,1,2)
clone_size_list = c(5,25,50,100)
# each with 1 variant?
for(vaf in vaf_list){
  for(clone_size in clone_size_list){
    for(coherence in coherence_list){
      X_mtx <- simulate_vaf_from_read_data(temp_plot_df,
                                           vaf_mu=vaf,vaf_sd=0.2,
                                           noise_mu=0.1,noise_sd=0.2,
                                           N_voi,clone_size,coherence,
                                           var_num=1,clone_num=3,noise_var_num=1,
                                           plot=F)
      saveRDS(X_mtx,paste0("scripts/simulation_experiment/data/sim_real/vaf_",
                                   vaf,"_clonesize_",clone_size,"_coherence_",coherence,".rds"))
    }
  }
}

# simulate cell-type ratio or use real distribution?
Ws = Ws[rownames(temp_plot_df),]

# bulk + spot level significance test
celltypes= colnames(Ws)
voi = paste0("Var_",1:4)
spatial_coords = read.csv("data/Sydney_Bracht/A200_higher_cov_20231125/a200_s6/a200_s6_MatchedBeadLocation.csv",header=T,row.names = 1)
spatial_coords = spatial_coords[colnames(af.dm),]
colnames(spatial_coords) = c("X","Y")

library(SpatialMT)
for(vaf in vaf_list){
  print(paste0("vaf:",vaf))
  for(clone_size in clone_size_list){
    print(paste0("clone_size:",clone_size))
    for(coherence in coherence_list){
      print(paste0("coherence:",coherence))
      # load simulated data
      sim_data = readRDS(paste0("scripts/simulation_experiment/data/sim_real/vaf_",
                                vaf,"_clonesize_",clone_size,"_coherence_",
                                coherence,".rds"))
      X_sim = sim_data$X_mtx
      N_sim = Matrix(rep(sim_data$N,dim(X_sim)[1]),byrow=T,nrow=dim(X_sim)[1])
      colnames(N_sim) = colnames(X_sim);rownames(N_sim) = rownames(X_sim)
      af.dm = X_sim/N_sim
      save_path=paste0("scripts/simulation_experiment/results/sim_real/",
             "vaf_",vaf,"_clonesize_",clone_size,"_coherence_",coherence,"/")
      # # bulk level test
      res_lg = celltype_test(celltypes = celltypes,voi=voi,
                             N=N_sim,
                             vaf=af.dm,
                             Ws=Ws,spatial_coords = spatial_coords,
                             test_type = "linear",plot=T,
                             save_path=save_path,
                             method="FDR")
      saveRDS(res_lg,paste0(save_path,"/res_lg.rds"))
      
      # spot-level significance test
      pdf(paste0(save_path,"/spot_level_significance.pdf"),width=15,height=15)
      for(var in voi){
        #colnames(N_voi) = sapply(strsplit(colnames(N_voi),"-1"),"[[",1)
        res=celltype_test_knn(celltypes, vars=var, N_sim[var,,drop=F], 
                              vaf=af.dm[var,,drop=F], Ws, spatial_coords, 
                          test_type = "linear", permute_num = 1000, 
                          k_neighbors = 100, method = "Raw", 
                          sample_idx = NULL, disease_celltype = "BE", 
                          ratio_threshold = 0.03, vaf_cellprop = F, 
                          exclude_plot_idx = NULL, 
                          p_thresh = 0.05, max_log10p_cap = 6, 
                          coef_plot_option = "negative")
        saveRDS(res,paste0(save_path,"/",var,"res_spot_level.rds"))
        print(res$combined_plot)
      }
      dev.off()
    }
  }
}

# caluclate bulk-level metric
library(dplyr)
clone_num=3;
var_num=1; noise_var_num=1;
vaf_list = c(0.1,0.25,0.5,0.75,0.9)
coherence_list = c(0,1,2)
clone_size_list = c(5,25,50,100)

# ground truth label
ground_truth = c(rep(paste0("BE"),each=var_num*clone_num),
                 rep("Noise",noise_var_num))
names(ground_truth) = paste("Var",1:length(ground_truth))

metric_list = c();
# change celltype ratio 
alpha_threshold=0.05;
celltypes = colnames(Ws)
for(clone_size in clone_size_list){
  for(vaf in vaf_list){
    for(coherence in coherence_list){
        save_path=paste0("scripts/simulation_experiment/results/sim_real/",
                                   "vaf_",vaf,"_clonesize_",clone_size,
                         "_coherence_",coherence,"/")
        results = readRDS(paste0(save_path,"/res_lg.rds"))
        adjusted_pval <- results$adjusted_pval
        coef_mat <- results$coef
        pred_var = rep("Noise",length(ground_truth))
        names(pred_var) = names(ground_truth)
        
        # should only be significant for BE cell-types
        for (i in 1:length(celltypes)){
          celltype_i  = celltypes[i]
          sig_idx <- which(adjusted_pval[, celltype_i] < alpha_threshold &
                             coef_mat[, celltype_i] > 0)
          pred_var[sig_idx] = paste0(celltype_i)
        }
        
        metric = c(clone_size,vaf,coherence,
                   unlist(calc_metric_realsim(pred_var,ground_truth,
                                              celltypes,sig_celltypes = c("BE"))))
        metric_list = rbind(metric_list,metric);
    }
  }
}
colnames(metric_list) = c("clone_size","vaf","coherence",
                          "overall_acc","avg_acc","weighted_f1",
                          "T1E","T2E")

metric_list = as.data.frame(metric_list)
saveRDS(metric_list,"scripts/simulation_experiment/results/sim_real/metric_list.rds")

library(ggplot2)
library(reshape2)
# Convert str to numeric values
metric_cols <- c("clone_size","vaf","coherence",
                 "overall_acc", "avg_acc", "weighted_f1", "T1E", "T2E")
metric_list[metric_cols] <- lapply(metric_list[metric_cols], as.numeric)

# Check structure
str(metric_list)  # Ensure all relevant columns are numeric

# Melt data for easier plotting
df_long <- melt(metric_list, id.vars = metric_cols ,
                measure.vars = c("overall_acc", "avg_acc", "weighted_f1", "T1E", "T2E"),
                variable.name = "Metric", value.name = "Value")
pdf("scripts/simulation_experiment/results/sim_real/clone_size_metric.pdf",width=10,height=10)
p<-ggplot(df_long, aes(x = as.factor(clone_size), y = Value, 
                       color = as.factor(vaf),
                       group = interaction(Metric, vaf))) +
  geom_line(size = 1,aes(alpha=0.3)) +
  geom_point(size = 2) +
  facet_grid(coherence~ Metric, scales = "free_y") +
  labs(x = "clone_size", y = "Metric Value",color = "VAF") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_color_viridis_d()
print(p)
dev.off()

pdf("scripts/simulation_experiment/results/sim_real/vaf_metric.pdf",width=10,height=10)
p<-ggplot(df_long, aes(x = as.factor(vaf), y = Value, 
                       color = as.factor(clone_size),
                       group = interaction(Metric, clone_size))) +
  geom_line(size = 1,aes(alpha=0.3)) +
  geom_point(size = 2) +
  facet_grid(coherence~ Metric, scales = "free_y") +
  labs(x = "vaf", y = "Metric Value",color = "clone_size") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_color_viridis_d()
print(p)
dev.off()


pdf("scripts/simulation_experiment/results/sim_real/coherence_metric.pdf",width=10,height=10)
p<-ggplot(df_long, aes(x = as.factor(coherence), y = Value, 
                       color = as.factor(clone_size),
                       group = interaction(Metric, clone_size))) +
  geom_line(size = 1,aes(alpha=0.3)) +
  geom_point(size = 2) +
  facet_grid(vaf ~ Metric, scales = "free_y") +
  labs(x = "Coherence", y = "Metric Value",color = "clone_size") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_color_viridis_d()
print(p)
dev.off()


# Spot-level significance
# caluclate bulk-level metric
clone_num=3;
var_num=1; noise_var_num=1;
vaf_list = c(0.1,0.25,0.5,0.75,0.9)
coherence_list = c(0,1,2)
clone_size_list = c(5,25,50,100)

metric_list = c();
vars = c(paste0("Var_",1:clone_num),"Noise")
# change celltype ratio 
alpha_threshold=0.05;
for(vaf in vaf_list){
  print(paste0("vaf: ",vaf))
  for(clone_size in clone_size_list){
    print(paste0("clone size: ",clone_size))
    for(coherence in coherence_list){
      print(paste0("coherence: ",coherence))
      for(j in 1:length(vars)){
        var = vars[j]
        save_path=paste0("scripts/simulation_experiment/results/sim_real/",
                         "vaf_",vaf,"_clonesize_",clone_size,
                         "_coherence_",coherence)
        results = readRDS(paste0(save_path,"/Var_",j,"res_spot_level.rds"))
        spot_results <- results$results
        coef_list <- t(sapply(spot_results, function(barcode_res) {
          barcode_res$coef
        }));colnames(coef_list) = celltypes 
        pval_list <- t(sapply(spot_results, function(barcode_res) {
          barcode_res$pval
        }));colnames(pval_list) = celltypes 
        pval_fdr_list <- spot_fdr(pval_list)
        # significant for BE cell-types
        pred_sig_spots = ((coef_list[,1] > 0) & (pval_fdr_list[,1] < 0.05))
        real_sig_spots = (temp_plot_df[names(pred_sig_spots),
                         paste0("CloneSize_",clone_size,"_Coherence_",coherence)] == j)
        names(real_sig_spots) = names(pred_sig_spots)
        real_sig_spots[is.na(real_sig_spots)] = F
        metrics <- c(clone_size,vaf,coherence,var,
                     unlist(calc_spot_level_metrics(pred_sig_spots, real_sig_spots)))
        # for DEBUG
        #metrics[is.na(metrics)] = 0
        metric_list <- rbind(metric_list,metrics)
        
      }
    }
  }
}

colnames(metric_list) = c("clone_size","vaf","coherence","Var",
                          "TP","TN","FP","FN",
                          "acc","precision","recall","f1",
                          "T1E","T2E")

metric_list = as.data.frame(metric_list)
saveRDS(metric_list,"scripts/simulation_experiment/results/sim_real/spot_level_metric_list.rds")

metric_list[is.na(metric_list)] = 0  

library(dplyr)
library(reshape2)
library(ggplot2)
# 1. Filter out the noise variant.
metric_list_real <- metric_list %>% filter(Var != "Noise")
#metric_list_real <- metric_list %>% subset(Var == "Var_1")

# Convert str to numeric values
metric_cols <- c("clone_size","vaf","coherence",
                 "TP","TN","FP","FN",
                 "acc","precision","recall","f1",
                 "T1E","T2E")
metric_list_real[metric_cols] <- lapply(metric_list_real[metric_cols], as.numeric)

# Check structure
str(metric_list_real)  # Ensure all relevant columns are numeric

# 2. Average across the 3 real variants for each condition combination.
# Here we assume the metric columns are stored as characters; if not, you may skip the as.numeric() conversion.
agg_metrics <- metric_list_real %>%
  group_by(clone_size, vaf, coherence) %>%
  summarise(
    acc       = mean(as.numeric(acc)),
    precision = mean(as.numeric(precision)),
    recall    = mean(as.numeric(recall)),
    f1        = mean(as.numeric(f1)),
    T1E       = mean(as.numeric(T1E)),
    T2E       = mean(as.numeric(T2E))
  ) %>%
  ungroup()

# 3. Reshape the aggregated data into long format.
# This will create columns: clone_size, vaf, coherence, Metric, and Value.
metrics_long <- melt(agg_metrics,
                     id.vars = c("clone_size", "vaf", "coherence"),
                     variable.name = "Metric",
                     value.name = "Value")

# 4. Create the line plot.
# X axis: clone_size, Y axis: Value; color by coherence (converted to factor),
# and facet by vaf (rows) and Metric (columns).
pdf("scripts/simulation_experiment/results/sim_real/spot_clone_size_metric.pdf",width=15,height=15)
p1 <- ggplot(metrics_long, aes(x = clone_size, y = Value, 
                              color = as.factor(coherence), 
                              group = as.factor(coherence))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ vaf + Metric, scales = "free_y") +
  labs(x = "Clone Size", y = "Metric Value", color = "Coherence",
       title = "Spot-Level Metrics (Averaged over Real Variants)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
        strip.text = element_text(face = "bold"))
print(p1)

metrics_long$clone_size <- factor(metrics_long$clone_size,levels=clone_size_list)
p2<- ggplot(metrics_long, aes(x = vaf, y = Value, 
                               color = as.factor(coherence), 
                               group = as.factor(coherence))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ clone_size + Metric, scales = "free_y",
             nrow = 4, ncol = 6, 
             drop = TRUE) +
  #facet_grid(clone_size ~ Metric, scales = "free_y") +
  labs(x = "VAF", y = "Metric Value", color = "Coherence",
       title = "Spot-Level Metrics (Averaged over Real Variants)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
        strip.text = element_text(face = "bold"))
print(p2)

p3<- ggplot(metrics_long, aes(x = coherence, y = Value, 
                              color = as.factor(clone_size), 
                              group = as.factor(clone_size))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ vaf + Metric, scales = "free_y") +
  labs(x = "Coherence", y = "Metric Value", color = "clone_size",
       title = "Spot-Level Metrics (Averaged over Real Variants)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
        strip.text = element_text(face = "bold"))

print(p3)

dev.off()

