library(grid)
library(gridExtra)
library(spacexr)
library(CARD)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)

# load data (X - alternate allele counts, N - total counts, Ws - celltype)
setwd("~/Documents/GitHub/SpatialMT/examples/")
source("utility_prev_literature/210215_FunctionsGeneral.R") # from MAESTER paper
#source("../diagnostic_plot.R")

# (1) mitochondrial data
#slide_name = "a200_s6"
maegatk.rse = readRDS("maegatk_mr1.rds") # master result object
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
# load variant of interest
voi = read.table("200_s6_voi_be25vaf.tsv",sep="\t")[,1]
if(slide_name=="a200_s6"){
  voi <- append(voi,"3071_T>C")
}
#subset data to the variant of interest
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]
# remove "-1" in the spot barcodes
colnames(af.dm) = sapply(strsplit(colnames(af.dm),"-"),"[[",1)
colnames(N) = sapply(strsplit(colnames(N),"-"),"[[",1)
# spatial coordinates
spatial_coords = read.csv(paste0(slide_name,"_MatchedBeadLocation.csv"))
rownames(spatial_coords) =spatial_coords[,1]; spatial_coords[,1] <- NULL
# Load celltype ratios
# RCTD ratio
rctd_ratio_major = readRDS("rctd_ratio_major.rds")
rctd_ratio_be = readRDS("rctd_ratio_be.rds")
# normalize the weights to sum to 1 in each spot
rctd_ratio_major = as.data.frame(as.matrix(normalize_weights(rctd_ratio_major)))
rctd_ratio_be = as.data.frame(as.matrix(normalize_weights(rctd_ratio_be)))

# (2) Determine celltype by linear regression
# choose a celltype ratio to work on
Ws = rctd_ratio_major[colnames(af.dm),]
celltypes= colnames(Ws)

# choose one weight/proportions for plotting
norm_weights = rctd_ratio_major
plot_df = cbind(spatial_coords[rownames(norm_weights),],norm_weights)
plot_df = cbind(plot_df[colnames(af.dm),],as.data.frame(t(af.dm*0.01)))
# also add in coverage
cov = N;rownames(cov) =  paste0(rownames(N),"_cov")
plot_df = cbind(plot_df,t(cov))

# for each variant j, fit a Xs/Nx vs Ws Linear regression
intercept_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
coef_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
pval_df=Matrix(NA,nrow=length(voi),ncol=length(celltypes))

pdf("linear_regression_anova_p_VAF_vs_celltype_plot.pdf",height=10,width=10)
for(j in 1:length(voi)){
  print(j)
  # LG
  var = voi[j]
  # weighted LG
  N_j = N[j,]
  vaf_j = af.dm[j,] # X_j/N_j
  for(k in 1:length(celltypes)){ # kx2
    #print(celltypes[k])
    #plot(Ws[,k],vaf_j,color)
    data = data.frame(vaf_j=vaf_j*0.01,N_j,W_sk=Ws[,k])
    # method 1: LG
    cell_ratio_bins = seq(from=0,to=1,by=0.05)
    fit=lm(vaf_j ~ W_sk,data=data,)# weights equals to inverse of variance (1/sigma_i^2)
    res=summary(fit)
    intercept_df[j,k]=res$coefficients[1,1]
    coef_df[j,k]= res$coefficients[2,1]
    pval_df[j,k] = res$coefficients[2,4]
  }
  # for variant j, plot out the fitted values on the diagnostic plots
  intercept=intercept_df[j,]
  coef=coef_df[j,]
  pval=pval_df[j,]
  plot_vaf_cellprop(j,af.dm,Ws,plot_df,intercept,coef,pval)
}
dev.off()

# weighted LG
intercept_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
coef_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
pval_df=Matrix(NA,nrow=length(voi),ncol=length(celltypes))
pdf("cell_mix/weighted_LG_permute_p_VAF_vs_celltype_plot.pdf",height=10,width=10)
for(j in 1:length(voi)){
  print(j)
  # LG
  var = voi[j]
  # weighted LG
  N_j = N[j,]
  vaf_j = af.dm[j,] # X_j/N_j
  for(k in 1:length(celltypes)){ # kx2
    #print(celltypes[k])
    #plot(Ws[,k],vaf_j,color)
    data = data.frame(vaf_j=vaf_j*0.01,N_j,W_sk=Ws[,k])
    # method 1: LG
    cell_ratio_bins = seq(from=0,to=1,by=0.05)
    fit=lm(vaf_j ~ W_sk,data=data,weights=sqrt(N_j))# weights equals to inverse of variance (1/sigma_i^2)
    res=summary(fit)
    intercept_df[j,k]=res$coefficients[1,1]
    coef_df[j,k]= res$coefficients[2,1]
    # doing 1000 permutations and see where the data locates
    set.seed(42)
    coef_list <- c()
    #  decouple VAF with Ws
    for(l in 1:1000){
      idx=sample(1:length(vaf_j),size=length(vaf_j))
      vaf_shuffled = vaf_j[idx] *0.01
      #vaf_shuffled = pmin(X_shuffled/(N_j+1e-6),1) # cap max VAF to be by 1
      fit=lm(vaf_shuffled ~ Ws[,k],weights=sqrt(N_j[idx]))# weights equals to inverse of variance (1/sigma_i^2)
      res_shuf=summary(fit)
      coef_list <- append(coef_list,res_shuf$coefficients[2,1])
    }
    pval_df[j,k] = 1- sum(res$coefficients[2,1] >= coef_list)/length(coef_list)
  }
  # for variant j, plot out the fitted values on the diagnostic plots
  intercept=intercept_df[j,]
  coef=coef_df[j,]
  pval=pval_df[j,]
  plot_vaf_cellprop(j,af.dm,Ws,plot_df,intercept,coef,pval,permute=T)
}
dev.off()



# Moran's I
spot_distances <- distances(spatial_coords[,1:2])
spot_dists_inv <- 1/(1+spot_distances)
#diag(spot_dists_inv) <- 0
moran_i_list <- list()
match_idx <- match(colnames(N),rownames(spatial_coords))
spot_dists_inv = spot_dists_inv[match_idx,match_idx]
diag(spot_dists_inv) <- 0

for(j in 1:length(voi)){
  var = voi[j]
  print(paste0(j,": ",var))
  #N_j = N[j,intersect_bc]
  #vaf_j = af.dm[j,intersect_bc] # X_j/N_j
  N_j = N[j,]
  vaf_j = af.dm[j,] # X_j/N_j
  res <- Moran.I(vaf_j, spot_dists_inv, alternative="greater")
  moran_i_list[[var]] <- res
}
saveRDS(moran_i_list,"voi_moran_i_list.rds")

mI <- unlist(lapply(voi,function(x){(moran_i_list[[x]])$observed}))
names(mI) <- voi
sort(mI,decreasing = T)[1:5]
#15642_T>A     14005_T>A      8588_T>G      6991_T>A     13025_A>C
#-0.0004701676 -0.0004945259 -0.0005914952 -0.0007012713 -0.0007636173

png("moranI_voi.png",width=500,height=500)
hist(mI,main="a200_S6 Moran's I")
dev.off()




