
setwd("~/nzhanglab/project/jrong/mito_LT/scripts/our_model/")
source("example_data/210215_FunctionsGeneral.R") # from MAESTER paper

# define LISI and Moran's I function
library(FNN)
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

# coverage & neighbor based moran's I
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

### load data
# 1.BE as beta test
# load transcriptomics/seeker data
seu <-readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/A200_s6/MAESTER/MAESTER_subsets/a200_s6_final.rds")
# load MT data
maegatk.rse = readRDS("../../data/Sydney_Bracht/A200_s6/MAESTER/MAESTER_subsets/maegatk_mr1.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
voi = c("3054_G>C","15777_G>C","3071_T>C")
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]
spatial_coords = seu@reductions$SPATIAL@cell.embeddings

A200_s6_sec1_LISI = pairwise_LISI(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                                  coverage_thresh=0,k = 50)
A200_s6_sec1_MI = moran_I_knn(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                              coverage_thresh=1,k=50,seed=42)

# 2.BE sec 2
seu <-readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/A200_s6_sec2/A200_s6_2_maester/analysis/A200s6_2_detailedannotations_seurat_final.rds")
# load MT data
maegatk.rse = readRDS("../../data/Sydney_Bracht/A200_s6_sec2/A200_s6_2_maester/analysis/A200s6_2_maegatk_final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
spatial_coords = seu@reductions$SPATIAL@cell.embeddings
# no pairwise, only MI
A200_s6_sec2_MI = moran_I_knn(voi=c("3054_G>C"),N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                              coverage_thresh=1,k=50,seed=42)

# 3. Small Bowel Disease
# load transcriptomics/seeker data
seu <-readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/EFF3/MAESTER/run1run2_comb/rds/EEF3_seurat_final.rds")
# load MT data
maegatk.rse <- readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/EFF3/MAESTER/run1run2_comb/rds/EEF3_maegatk_final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
# load variant of interest
voi = c("9811_G>A","9797_T>C")
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]
spatial_coords = seu@reductions$SPATIAL@cell.embeddings
EEF3_LISI = pairwise_LISI(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                                  coverage_thresh=0,k = 30)
EEF3_MI = moran_I_knn(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                              coverage_thresh=1,k=30,seed=42)

# 4. A210
# load transcriptomics/seeker data
seu <-readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/A210/combined/analysis/seurat_final.rds")
# load MT data
maegatk.rse <- readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/A210/combined/analysis/maegatk_final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
# load variant of interest
#voi = read.table(paste0(data_path,slide_name,"_voi_be25vaf.tsv"),sep="\t")[,1]
voi = c("2245_A>G","11969_G>A","3866_T>C")
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]
spatial_coords = seu@reductions$SPATIAL@cell.embeddings

A210_LISI = pairwise_LISI(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                          coverage_thresh=0,k = 30)
A210_MI = moran_I_knn(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                      coverage_thresh=1,k=30,seed=42)

# 5. CRC
seu <-readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/CRC076_C/seeker/run2/Seurat/CRC076C_seurat_final.rds")
# load MT data
maegatk.rse <- readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/CRC076_C/seeker/run2/Seurat/CRC076C_maegatk_final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
# load variant of interest
voi = c("4239_C>A","9982_G>T","5380_T>A","3957_A>C",
            "7411_A>C","11892_A>C")
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]
spatial_coords = seu@reductions$SPATIAL@cell.embeddings

CRC_LISI = pairwise_LISI(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                          coverage_thresh=0,k = 30)
CRC_MI = moran_I_knn(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                      coverage_thresh=1,k=30,seed=42)

# 5. cellline mixture
maegatk.rse = readRDS("../../data/Sydney_Bracht/mouse_cellline_mix/MAESTER/MT-R_maegatk_final.rds") # master result object
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
voi = c("3010_G>A","9545_A>G","8002_C>T","11152_T>C","10176_G>A","13500_T>C")
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]
spatial_obj <- readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/mouse_cellline_mix/MAESTER/MT-R_seurat_final.rds")
spatial_obj <- subset(spatial_obj, nFeature_RNA > 20)
spatial_coords = as.data.frame(spatial_obj@reductions$SPATIAL@cell.embeddings)

CDX_LISI = pairwise_LISI(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                         coverage_thresh=0,k = 50)
CDX_MI = moran_I_knn(voi=voi,N=N,af.dm=af.dm,spatial_coords=spatial_coords,
                     coverage_thresh=1,k=50,seed=42)

# Save as table for future look up
LISI_list = list(A200_s6_sec1_LISI, A210_LISI, EEF3_LISI, CDX_LISI, CRC_LISI)
MI_list = list(A200_s6_sec1_MI, A200_s6_sec2_MI, A210_MI, 
               EEF3_MI, CDX_MI, CRC_MI)

names(LISI_list) <- c("A200_s6_sec1", "A210", "EEF3", "CDX", "CRC")

# Initialize an empty list to collect results
lisi_avg_list <- list()

# Loop through each sample's LISI result
for (sample_name in names(LISI_list)) {
  sample_lisi <- LISI_list[[sample_name]]
  
  # For each variant pair, extract the "avg" value
  if (length(sample_lisi) > 0) {
    for (pair in names(sample_lisi)) {
      avg_value <- sample_lisi[[pair]]["avg"]
      lisi_avg_list[[length(lisi_avg_list) + 1]] <- data.frame(
        sample = sample_name,
        variant_pair = pair,
        avg_lisi = avg_value
      )
    }
  }
}

# Combine all into one data frame
lisi_avg_df <- do.call(rbind, lisi_avg_list)

# Write to txt and csv
write.table(lisi_avg_df, file = "LISI_allpairs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(lisi_avg_df, file = "LISI_allpairs.csv", row.names = FALSE)


sample_names <- c("A200_s6_sec1", "A200_s6_sec2", "A210", "EEF3", "CDX", "CRC")
names(MI_list) <- sample_names
# Combine all Moran's I data.frames into one
mi_combined_df <- do.call(rbind, lapply(names(MI_list), function(sample_name) {
  df <- as.data.frame(MI_list[[sample_name]])
  df$variant <- paste0(sample_name,"_",rownames(df))
  df$sample <- sample_name
  rownames(df) <- NULL
  return(df)
}))
# Rename and reorder columns
colnames(mi_combined_df)[colnames(mi_combined_df) == "mi_list"] <- "mi_value"
colnames(mi_combined_df)[colnames(mi_combined_df) == "p_list"] <- "p_value"
mi_combined_df <- mi_combined_df[, c("sample", "variant", "mi_value", "p_value")]

# Export
write.table(mi_combined_df, "moran_i_combined.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(mi_combined_df, "moran_i_combined.csv", row.names = FALSE)
