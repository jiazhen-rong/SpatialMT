library(Seurat)
library(ggplot2)
library(Matrix)
library(spacexr)
library(gridExtra)
library(pheatmap)

setwd("~/nzhanglab/project/jrong/mito_LT/scripts/simulation_experiment/")
source("../our_model/example_data/210215_FunctionsGeneral.R") # from MAESTER paper
source("simulation_functions.R")

# Define clone locations - clone size 10 
clone_size=10
clones <- list(
  list(id = 1, x_start = 3, y_start = 6, width = 5, height = 2),
  list(id = 2, x_start = 21, y_start = 11, width = 5, height = 2),
  list(id = 3, x_start = 11, y_start = 21, width = 5, height = 2)
)
clone_colors = c("orange", "red", "blue", "lightgrey")
plot_df = simulate_clone_location(grid_size=30,clones,
                                  plot=F,clone_colors =clone_colors,
                                  save=T,
                                  save_path="data/clone_3_size_10/",
                                  width=10,height=10)
for(vf in c(0.1,0.25,0.5,0.75,0.9)){
  for(n in c(2,5,8,10)){
    X_mtx <- simulate_vaf(plot_df,grid_size = 30, 
                            var_num = 5, noise_var_num=15,
                            N = n, clone_num = 3, clone_size = 50, 
                            vaf = vf, background_vaf = 0.1,
                            save=T,save_path="data/clone_3_size_10/")
  }
}


celltypes = c("C1","C2","C3","C4")
c1_idx=which((plot_df$x_coord >= 1) & (plot_df$x_coord <= 20) & 
               (plot_df$y_coord >= 1) & (plot_df$y_coord <= 10))
c2_idx=which((plot_df$x_coord >= 11) & (plot_df$x_coord <= 30) & 
               (plot_df$y_coord >= 11) & (plot_df$y_coord <= 20))
c3_idx=which((plot_df$x_coord >= 1) & (plot_df$x_coord <= 20) & 
               (plot_df$y_coord >= 21) & (plot_df$y_coord <= 30))
c4_idx = setdiff(1:(grid_size^2),c(c1_idx,c2_idx,c3_idx))

celltype_idx = list(
  "1"= c1_idx,
  "2" = c2_idx,
  "3" = c3_idx,
  "4" = c4_idx
)
for(celltype_ratio in c(0.1,0.25,0.5,0.75)){
  print(celltype_ratio)
  Ws <- simulate_celltype_ratio(plot_df=data.frame(plot_df), 
                                grid_size = 30, 
                                celltypes = celltypes,
                                celltype_ratio = celltype_ratio, sd = 0.2,
                                save=T,save_path="data/clone_3_size_10/",plot=T,
                                distribution="Beta",
                                cell_range="Clone")
}

for(celltype_ratio in c(0.1,0.25,0.5,0.75)){
  print(celltype_ratio)
  Ws <- simulate_celltype_ratio(plot_df=data.frame(plot_df), 
                                grid_size = 30, 
                                celltypes = celltypes,
                                celltype_ratio = celltype_ratio, sd = 0.2,
                                save=T,save_path="data/clone_3_size_10/",plot=T,
                                distribution="Beta",
                                cell_range="Custom",
                                celltype_idx=celltype_idx)
}

# Define clone locations - clone size 10 
clone_size=25;
clones <- list(
  list(id = 1, x_start = 3, y_start = 3, width = 5, height = 5),
  list(id = 2, x_start = 21, y_start = 11, width = 5, height = 5),
  list(id = 3, x_start = 11, y_start = 21, width = 5, height = 5)
)
clone_colors = c("orange", "red", "blue", "lightgrey")
plot_df = simulate_clone_location(grid_size=30,clones,
                                  plot=F,clone_colors =clone_colors,
                                  save=T,
                                  save_path="data/clone_3_size_25/",
                                  width=10,height=10)
for(vf in c(0.1,0.25,0.5,0.75,0.9)){
  for(n in c(2,5,8,10)){
    X_mtx <- simulate_vaf(plot_df,grid_size = 30, 
                          var_num = 5, noise_var_num=15,
                          N = n, clone_num = 3, clone_size = 50, 
                          vaf = vf, background_vaf = 0.1,
                          save=T,save_path="data/clone_3_size_25/")
  }
}
