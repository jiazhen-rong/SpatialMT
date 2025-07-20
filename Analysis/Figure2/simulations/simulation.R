library(Seurat)
library(ggplot2)
library(Matrix)
library(spacexr)
library(gridExtra)
library(pheatmap)

setwd("~/nzhanglab/project/jrong/mito_LT/scripts/simulation_experiment/")
source("../our_model/example_data/210215_FunctionsGeneral.R") # from MAESTER paper
source("simulation_functions.R")

# 3 Clones with Size 50 
# Define clone locations
clones <- list(
  list(id = 1, x_start = 3, y_start = 6, width = 10, height = 5),
  list(id = 2, x_start = 21, y_start = 11, width = 10, height = 5),
  list(id = 3, x_start = 11, y_start = 21, width = 10, height = 5)
)
clone_colors = c("orange", "red", "blue", "lightgrey")
plot_df = simulate_clone_location(grid_size=30,clones,
                                  plot=F,clone_colors =clone_colors,
                                  save=T,
                                  save_path="data/clone_3_size_50/",
                                  width=10,height=10)
for(vf in c(0.1,0.25,0.75,0.9)){
  for(n in c(2,5,8,10)){
  X_mtx <- simulate_vaf(plot_df,grid_size = 30, 
                          var_num = 5, noise_var_num=15,
                          N = n, clone_num = 3, clone_size = 50, 
                          vaf = vf, background_vaf = 0.1,
                          save=T,save_path="data/clone_3_size_50/")
  }
}


celltypes = c("C1","C2","C3","C4")
for(celltype_ratio in c(0.25,0.26,0.3,0.4,0.5,0.75)){
  print(celltype_ratio)
  Ws <- simulate_celltype_ratio(plot_df=data.frame(plot_df), 
                                grid_size = 30, 
                                celltypes = celltypes,
                                celltype_ratio = celltype_ratio, sd = 0.2,
                                save=T,save_path="data/clone_3_size_50/",plot=T,
                                distribution="Beta")
}

