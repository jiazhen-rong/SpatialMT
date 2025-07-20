library(SpatialMT)
library(Seurat)
library(ggplot2)
library(Matrix)
library(spacexr)
library(gridExtra)
library(pheatmap)

setwd("~/nzhanglab/project/jrong/mito_LT/scripts/simulation_experiment/")
source("../our_model/example_data/210215_FunctionsGeneral.R") # from MAESTER paper

# Load Simulated Data
clone_size=10#50;
clone_num=3;
vaf_list = c(0.1,0.25,0.5,0.75,0.9); 
n_list = c(2,5,8,10); 
celltype_ratio_list=c(0.1,0.25,0.5,0.75);
var_num=5; noise_var_num=15;
cell_range="Custom"
path = paste0("data/clone_",clone_num,"_size_",clone_size,"/")
for(vaf in vaf_list){
  print(paste0("vaf: ",vaf))
  for(n in n_list){
    print(paste0("N: ",n))
    for(celltype_ratio in celltype_ratio_list){
      print(paste0("Ratio: ",celltype_ratio," Range: ",cell_range))
      X = readRDS(paste0(path,"/vaf_",vaf,"_N_",n,".rds"));colnames(X) = paste0("BC_",1:ncol(X))
      Ws = readRDS(paste0(path,"/sim_celltype_ratio_",
                          celltype_ratio,"_",cell_range,".rds"))
      rownames(Ws) = paste0("BC_",1:ncol(X)); 
      celltypes = colnames(Ws)
      N = Matrix(n,nrow=nrow(X),ncol=ncol(X))
      colnames(N) = paste0("BC_",1:ncol(X)); rownames(N) = rownames(X)
      af.dm = as.matrix(X/N)
      colnames(af.dm) = paste0("BC_",1:ncol(X)); rownames(af.dm) = rownames(X)
      voi = rownames(X)
      spatial_coords = readRDS(paste0(path,"/spatial_coords.rds"))
      rownames(spatial_coords) = paste0("BC_",1:ncol(X))
      
      # using SpatialMT package
      save_path=paste0("results/","/clone_",clone_num,"_size_",clone_size,
                       "/vaf_",vaf,"_N_",n,"_ratio_",
                       celltype_ratio,"_range_",cell_range,"/")
      dir.create(save_path,recursive = T)
      res_lg = celltype_test(celltypes = celltypes,voi=voi,N=N,
                             vaf=af.dm,Ws=Ws,spatial_coords = spatial_coords,
                             test_type = "linear",plot=T,
                             save_path=save_path,method="FDR",
                             figure_height = 20, figure_width = 20)
      saveRDS(res_lg,paste0(save_path,"/res_lg.rds"))
      
      sig_var = filter_significant_variants(res_lg,alpha_threshold = 0.05,
                                            output_file=paste0(save_path,"sig_vars_fdr.txt"))
    }
  }
}

# Load Simulated Data
clone_size=50;
clone_num=3;
vaf_list = c(0.1,0.25,0.5,0.75,0.9); 
n_list = c(2,5,8,10); 
celltype_ratio_list=c(0.1,0.25,0.5,0.75);
var_num=5; noise_var_num=15;
cell_range="Custom"
path = paste0("data/clone_",clone_num,"_size_",clone_size,"/")
for(vaf in vaf_list){
  print(paste0("vaf: ",vaf))
  for(n in n_list){
    print(paste0("N: ",n))
    for(celltype_ratio in celltype_ratio_list){
      print(paste0("Ratio: ",celltype_ratio," Range: ",cell_range))
      X = readRDS(paste0(path,"/vaf_",vaf,"_N_",n,".rds"));colnames(X) = paste0("BC_",1:ncol(X))
      Ws = readRDS(paste0(path,"/sim_celltype_ratio_",
                          celltype_ratio,"_",cell_range,".rds"))
      rownames(Ws) = paste0("BC_",1:ncol(X)); 
      celltypes = colnames(Ws)
      N = Matrix(n,nrow=nrow(X),ncol=ncol(X))
      colnames(N) = paste0("BC_",1:ncol(X)); rownames(N) = rownames(X)
      af.dm = as.matrix(X/N)
      colnames(af.dm) = paste0("BC_",1:ncol(X)); rownames(af.dm) = rownames(X)
      voi = rownames(X)
      spatial_coords = readRDS(paste0(path,"/spatial_coords.rds"))
      rownames(spatial_coords) = paste0("BC_",1:ncol(X))
      
      # using SpatialMT package
      save_path=paste0("results/","/clone_",clone_num,"_size_",clone_size,
                       "/vaf_",vaf,"_N_",n,"_ratio_",
                       celltype_ratio,"_range_",cell_range,"/")
      dir.create(save_path,recursive = T)
      res_lg = celltype_test(celltypes = celltypes,voi=voi,N=N,
                             vaf=af.dm,Ws=Ws,spatial_coords = spatial_coords,
                             test_type = "linear",plot=T,
                             save_path=save_path,method="FDR",
                             figure_height = 20, figure_width = 20)
      saveRDS(res_lg,paste0(save_path,"/res_lg.rds"))
      
      sig_var = filter_significant_variants(res_lg,alpha_threshold = 0.05,
                                            output_file=paste0(save_path,"sig_vars_fdr.txt"))
    }
  }
}

clone_size=25;
clone_num=3;
vaf_list = c(0.1,0.25,0.5,0.75,0.9); 
n_list = c(2,5,8,10); 
celltype_ratio_list=c(0.1,0.25,0.5,0.75);
var_num=5; noise_var_num=15;
cell_range="Custom"
path = paste0("data/clone_",clone_num,"_size_",clone_size,"/")

for(vaf in vaf_list){
  print(paste0("vaf: ",vaf))
  for(n in n_list){
    print(paste0("N: ",n))
    for(celltype_ratio in celltype_ratio_list){
      print(paste0("Ratio: ",celltype_ratio," Range: ",cell_range))
      X = readRDS(paste0(path,"/vaf_",vaf,"_N_",n,".rds"));colnames(X) = paste0("BC_",1:ncol(X))
      Ws = readRDS(paste0(path,"/sim_celltype_ratio_",
                          celltype_ratio,"_",cell_range,".rds"))
      rownames(Ws) = paste0("BC_",1:ncol(X)); 
      celltypes = colnames(Ws)
      N = Matrix(n,nrow=nrow(X),ncol=ncol(X))
      colnames(N) = paste0("BC_",1:ncol(X)); rownames(N) = rownames(X)
      af.dm = as.matrix(X/N)
      colnames(af.dm) = paste0("BC_",1:ncol(X)); rownames(af.dm) = rownames(X)
      voi = rownames(X)
      spatial_coords = readRDS(paste0(path,"/spatial_coords.rds"))
      rownames(spatial_coords) = paste0("BC_",1:ncol(X))
      
      # using SpatialMT package
      save_path=paste0("results/","/clone_",clone_num,"_size_",clone_size,
                       "/vaf_",vaf,"_N_",n,"_ratio_",
                       celltype_ratio,"_range_",cell_range,"/")
      dir.create(save_path,recursive = T)
      res_lg = celltype_test(celltypes = celltypes,voi=voi,N=N,
                             vaf=af.dm,Ws=Ws,spatial_coords = spatial_coords,
                             test_type = "linear",plot=T,
                             save_path=save_path,method="FDR",
                             figure_height = 20, figure_width = 20)
      saveRDS(res_lg,paste0(save_path,"/res_lg.rds"))
      
      sig_var = filter_significant_variants(res_lg,alpha_threshold = 0.05,
                                            output_file=paste0(save_path,"sig_vars_fdr.txt"))
    }
  }
}

# Calculate Acc of bulk level test
clone_size=10#50;
clone_num=3;
vaf_list = c(0.1,0.25,0.5,0.75,0.9); 
n_list = c(2,5,8,10); 
celltype_ratio_list=c(0.1,0.25,0.5,0.75);
var_num=5; noise_var_num=15;
cell_range="Custom"

# ground truth label
ground_truth = c(rep(paste0("Clone_",1:clone_num),each=var_num),
                 rep("Noise",noise_var_num))
names(ground_truth) = paste("Var",1:length(ground_truth))

metric_list = c();
# change celltype ratio 
alpha_threshold=0.05;
for(clone_size in c(10,25,50)){
  for(vaf in vaf_list){
    for(n in n_list){
      for(cr in celltype_ratio_list){
        save_path=paste0("results/","/clone_",clone_num,"_size_",clone_size,
                         "/vaf_",vaf,"_N_",n,"_ratio_",
                         cr,"_range_",cell_range,"/")
        results = readRDS(paste0(save_path,"/res_lg.rds"))
        adjusted_pval <- results$adjusted_pval
        coef_mat <- results$coef
        pred_var = rep("Noise",length(ground_truth))
        names(pred_var) = names(ground_truth)
        
        for (i in 1:(length(celltypes)-1)){
          celltype_i  = celltypes[i]
          sig_idx <- which(adjusted_pval[, celltype_i] < alpha_threshold & 
                             coef_mat[, celltype_i] > 0)
          pred_var[sig_idx] = paste0("Clone_",i)
        }
        metric = c(clone_size,n,vaf,cr,cell_range,
                   unlist(calc_metric(pred_var,ground_truth,clone_num)))
        metric_list = rbind(metric_list,metric);
      }
    }
  }
}
colnames(metric_list) = c("clone_size","N","vaf","celltype_ratio",
                          "cell_range","overall_acc","avg_acc","weighted_f1",
                          "T1E","T2E")
                          
metric_list = as.data.frame(metric_list)
saveRDS(metric_list,"results/metric_list.rds")

library(ggplot2)
library(reshape2)
metric_list <- readRDS("results/metric_list.rds")
# Convert str to numeric values
metric_cols <- c("clone_size","N","vaf","celltype_ratio",
                 "overall_acc", "avg_acc", "weighted_f1", "T1E", "T2E")
metric_list[metric_cols] <- lapply(metric_list[metric_cols], as.numeric)

# Check structure
str(metric_list)  # Ensure all relevant columns are numeric

# Melt data for easier plotting
df_long <- melt(metric_list, id.vars = c("clone_size", "N", "vaf", "celltype_ratio", "cell_range"),
                measure.vars = c("overall_acc", "avg_acc", "weighted_f1", "T1E", "T2E"),
                variable.name = "Metric", value.name = "Value")

df_long_clone_25 <- df_long[df_long$clone_size == 25,]
df_long_clone_50 <- df_long[df_long$clone_size == 50,]
df_long_clone_10 <- df_long[df_long$clone_size == 10,]
# Facetted Line Plots
pdf("results/grid_bulk_metrics.pdf",width=15,height=15)
p1 <- ggplot(df_long_clone_10, aes(x = as.numeric(N), y = Value, 
                          color = as.factor(celltype_ratio), 
                          group = as.factor(celltype_ratio))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(vaf ~ Metric, scales = "free_y") +
  scale_color_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "N", y = "Metric Value") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2 <- ggplot(df_long_clone_25, aes(x = as.numeric(N), y = Value, 
                                   color = as.factor(celltype_ratio), 
                                   group = as.factor(celltype_ratio))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(vaf ~ Metric, scales = "free_y") +
  scale_color_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "N", y = "Metric Value") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p3 <- ggplot(df_long_clone_50, aes(x = as.numeric(N), y = Value, 
                                   color = as.factor(celltype_ratio), 
                                   group = as.factor(celltype_ratio))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(vaf ~ Metric, scales = "free_y") +
  scale_color_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "N", y = "Metric Value") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Keep only weighted_f1
df_f1 <- metric_list %>%
  select(clone_size, N, vaf, celltype_ratio, weighted_f1)

p4 <- ggplot(df_f1, aes(x = as.numeric(clone_size), y = weighted_f1, 
                             color = as.factor(celltype_ratio), 
                             group = as.factor(celltype_ratio))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(vaf ~ N) +
  scale_color_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "Clone Size", y = "Weighted F1 Score") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p4)

p5 <- ggplot(df_f1, aes(x = as.numeric(clone_size), y = weighted_f1, 
                        color = as.factor(N), 
                        group = as.factor(N))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(celltype_ratio ~ vaf) +
  scale_color_brewer(palette = "Set2", name = "N") +
  labs(x = "Clone Size", y = "Weighted F1 Score") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p5)

p6 <- ggplot(df_f1, aes(x = as.numeric(clone_size), y = weighted_f1, 
                        color = as.factor(vaf), 
                        group = as.factor(vaf))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(celltype_ratio ~ N) +
  scale_color_brewer(palette = "Set2", name = "vaf") +
  labs(x = "Clone_size", y = "Weighted F1 Score") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p7 <- ggplot(df_f1, aes(x = as.numeric(vaf), y = weighted_f1, 
                        color = as.factor(clone_size), 
                        group = as.factor(clone_size))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(celltype_ratio ~ N) +
  scale_color_brewer(palette = "Set2", name = "vaf") +
  labs(x = "VAF", y = "Weighted F1 Score") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf("results/grid_bulk_metrics.pdf",width=15,height=15)
print(p1);print(p2);print(p3);print(p4);print(p5);print(p6);print(p7)
dev.off()

##################################3

p2 <- ggplot(df_long, aes(x = as.numeric(N), y = Value, 
                          color = as.factor(clone_size), 
                          group = as.factor(clone_size))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_grid(vaf ~ Metric, scales = "free_y") +
  scale_color_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "N", y = "Metric Value") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p1)

# Create faceted bar plot
# Faceted plot: N on x-axis, facet by celltype_ratio (rows) and vaf (columns)
pdf("results/clone_3_size_10/metric_line_plots.pdf",width=15,height=10)
p<-ggplot(df_long, aes(x = as.factor(N), y = Value, fill = as.factor(celltype_ratio))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(vaf ~ Metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "N", y = "Metric Value") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p)
dev.off()

# Keep only weighted_f1
df_f1 <- metric_list %>%
  select(clone_size, N, vaf, celltype_ratio, weighted_f1)

pdf("results/clone_metric.pdf",width=15,height=10)
# Plot
p2<-ggplot(df_f1, aes(x = as.factor(clone_size), y = weighted_f1, fill = as.factor(celltype_ratio))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(vaf ~ N) +
  scale_fill_brewer(palette = "Set2", name = "Celltype Ratio") +
  labs(x = "Clone Size", y = "Weighted F1 Score") +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(color = "grey50", fill = "grey95", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p2)
dev.off()

