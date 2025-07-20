library(Seurat)
library(ggplot2)
library(Matrix)
library(spacexr)
library(gridExtra)


setwd("~/nzhanglab/project/jrong/mito_LT/")
source("scripts/our_model/example_data/210215_FunctionsGeneral.R") # from MAESTER paper

# load bulk level frequency
bulk_vaf = read.table("data/Sydney_Bracht/A200_higher_cov_20231125/TWIST/gatk_output/filtered/vcf_files/mutect2_a200_s6_mito.vcf.gz")
bulk_vaf = bulk_vaf[(unlist(lapply(bulk_vaf[,4],nchar))==1) & (unlist(lapply(bulk_vaf[,5],nchar))==1),]
real_vaf = paste0(bulk_vaf[,2],"_",bulk_vaf[,4],">",bulk_vaf[,5])
temp_vaf = as.numeric(sapply(strsplit(bulk_vaf[,10],":"),"[[",3))

png("results/Sydney_Bracht/Simulation/bulk_hist.png",width=500,height=300)
hist(temp_vaf,breaks=100,main="Bulk VAF in Real A200_S6")
dev.off()

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

# remove "-1" in the spot barcodes
colnames(af.dm) = sapply(strsplit(colnames(af.dm),"-"),"[[",1)
colnames(N) = sapply(strsplit(colnames(N),"-"),"[[",1)

# spatial coordinates 
#spatialRNA_rds = readRDS(paste0("data/Sydney_Bracht/A200_higher_cov_20231125/a200_s6/a200_s6_seurat.rds"))
#spatialRNA_rds = subset(spatialRNA_rds,nCount_RNA >=30)
#SpatialDimPlot(spatialRNA_rds)

# take intersections - 6003 out of 6008

## Diagnostic plot
# N + X + vaf distribution
N_real = as.data.frame(t(N[real_vaf,]))
af.dm_real = as.data.frame(t(af.dm[real_vaf,]))
hist(colSums(N_real))
sum(colSums(N_real) > 5)
N_voi = as.data.frame(t(N[voi[,1],]))
af.dm_voi = as.data.frame(t(af.dm[voi[,1],]))

#ggplot(N_real,aes(x=`64_C>T`)) + geom_density(alpha=0.5)
png("results/Sydney_Bracht/Simulation/real_coverage.png",width=500,height=300)
plot(density(pmin(N_real[,1],10)),main="Coverage Distribution",
     ylim=c(0,10),xlim=c(-1,11))
for(i in 1:(length(real_vaf)-1)){
  lines(density(pmin(N_real[,i],10))  )
}
dev.off()

png("results/Sydney_Bracht/Simulation/real_vaf.png",width=500,height=300)
plot(density(af.dm_real[,1]),main="VAF Distribution",ylim=c(0,0.5),xlim=c(-1,100))
for(i in 1:(length(real_vaf)-1)){
  lines(density(af.dm_real[,i]))  
}
dev.off()

png("results/Sydney_Bracht/Simulation/real_vaf.png",width=500,height=300)
plot(density(af.dm_real[,1]),main="VAF Distribution",ylim=c(0,10))
for(i in 1:(length(real_vaf)-1)){
  lines(density(af.dm_real[,i]))  
}
dev.off()

voi = read.csv("data/Sydney_Bracht/A200_s6/MAESTER/200_s6_voi_be25vaf.tsv")
intersect(voi[,1],real_vaf)

# scatter plot - VAF vs VAF, color by coverage
plot_df = data.frame(bulk_vaf=temp_vaf,
                     MAESTER_vaf=colMeans(af.dm_real)/100,
                     MAESTER_coverage_log = log10(colSums(N_real)+1))

png("results/Sydney_Bracht/Simulation/Scatter_plot_VAF.png",width=600,height=500)
ggplot(plot_df,aes(x=bulk_vaf,y=MAESTER_vaf,color=MAESTER_coverage_log))+ geom_point(size=5,alpha=0.5)+
  theme_bw() + 
  scale_color_continuous(low="grey",high="red")+
  xlab("Bulk VAF from TWIST") + ylab("Average MAESTER VAF") +
  ggtitle("Variants from TWIST")
dev.off()


# Simulate coordinates
x_coord = rep(c(1:30),each=30);y_coord=rep(c(1:30),30)
plot_df = data.frame(x_coord=x_coord,y_coord=y_coord)

plot_df$label = "Other"
plot_df[(plot_df$x_coord %in% c(5:10)) & (plot_df$y_coord %in% c(5:10)),"label"] = "Clone_1" 
plot_df[(plot_df$x_coord %in% c(20:25)) & (plot_df$y_coord %in% c(10:15)),"label"] = "Clone_2" 
plot_df[(plot_df$x_coord %in% c(10:15)) & (plot_df$y_coord %in% c(20:25)),"label"] = "Clone_3" 

png("results/Sydney_Bracht/Simulation/Simulated_Clones.png",width=500,height=500)
ggplot(plot_df,aes(x=x_coord,y=y_coord,color=label)) + geom_point() + theme_bw()+
  scale_color_manual(values = c("orange","red","blue","lightgrey"))
dev.off()

png("results/Sydney_Bracht/Simulation/Simulated_Clones.png",width=500,height=500)
plot_df_random = plot_df
plot_df_random[,c(1,2)] = plot_df[sample(1:50),c(1,2)]
ggplot(plot_df_random,aes(x=x_coord,y=y_coord,color=label)) + geom_point() + theme_bw()+
  scale_color_manual(values = c("orange","red","blue","lightgrey"))
dev.off()

# For each clone, simulate 5 vairant
library(Matrix)
library(pheatmap)
vaf_mtx = matrix(0,nrow=3*5,ncol=dim(plot_df)[1])
rownames(vaf_mtx) = paste0("Var_",c(1:15))

for(i in 1:5){
  vaf_mtx[i,plot_df$label == "Clone_1"] = rbinom(size=20,
                              n=sum(plot_df$label == "Clone_1"),prob=0.1)
}
for(i in 5:10){
  vaf_mtx[i,plot_df$label == "Clone_2"] = rbinom(size=20,
                              n=sum(plot_df$label == "Clone_2"),prob=0.1)
}
for(i in 10:15){
  vaf_mtx[i,plot_df$label == "Clone_3"] = rbinom(size=20,
                                                 n=sum(plot_df$label == "Clone_3"),prob=0.1)
}
png("results/Sydney_Bracht/Simulation/Simulated_Clones_heatmap.png")
pheatmap(vaf_mtx[,plot_df$label != "Other"])
dev.off()

png("results/Sydney_Bracht/Simulation/Simulated_Clones_heatmap_allcells.png")
pheatmap(vaf_mtx)
dev.off()

####### Simulate from real distributions - high quality spots, 6008 spots
seu <- readRDS("data/Sydney_Bracht/A200_higher_cov_20231125/a200_s6/a200_s6_seurat.rds")
spatial_coords = read.csv("data/Sydney_Bracht/A200_higher_cov_20231125/a200_s6/a200_s6_MatchedBeadLocation.csv")
rownames(spatial_coords) = spatial_coords[,1]
spatial_coords = spatial_coords[colnames(af.dm),]

# load RCTD celltype decompositions
rctd_res = readRDS("./results/Sydney_Bracht/sample_level/a200_s6_high_cov/RCTD/Major/RCTD_res.rds")
sample_weights = as.matrix(normalize_weights(rctd_res@results$weights))

celltype0 = data.frame(barcodes=rownames(spatial_coords),celltype="Rest") 
rownames(celltype0) = celltype0[,1]
celltype0[intersect(rownames(sample_weights),rownames(celltype0)),"celltype"] = 
  colnames(sample_weights)[apply(sample_weights[intersect(rownames(sample_weights),rownames(celltype0)),],
                                 1,function(x){which.max(x)})]

plot_df = cbind(spatial_coords,celltype0[rownames(spatial_coords),2])
colnames(plot_df) = c("barcodes","spatial_1","spatial_2","celltype_rctd")


# bottom
plot_df[(-plot_df$spatial_1 < -4900) & (-plot_df$spatial_1 > -5000) &
          (-plot_df$spatial_2 < -4250) & (-plot_df$spatial_2 > -4300) &
          (plot_df$celltype_rctd =="BE"),"barcodes"]
# upper left
plot_df[(-plot_df$spatial_1 < -1750) & (-plot_df$spatial_1 > -1850) &
          (-plot_df$spatial_2 < -4750) & (-plot_df$spatial_2 > -4800) &
          (plot_df$celltype_rctd =="BE"),"barcodes"]
# upper right
plot_df[(-plot_df$spatial_1 < -2100) & (-plot_df$spatial_1 > -2200) &
          (-plot_df$spatial_2 < -3250) & (-plot_df$spatial_2 > -3500) &
          (plot_df$celltype_rctd =="BE"),"barcodes"]

center1 = "ACAGCGCGGGCGAC"
center2 = "AGCAGTCTACGGGG"
center3 = "AATGAACGATATAT"

p1 <- ggplot(plot_df,aes(x=-spatial_2,y=-spatial_1,color=celltype_rctd)) +
  geom_point(size=0.1) + theme_bw() + ggrastr::geom_point_rast(raster.dpi = 100) + 
  geom_point(data =plot_df[c(center1,center2,center3),],aes(x=-spatial_2,y=-spatial_1),size=3,color="red") + 
  coord_fixed(ratio=1) + ggtitle("Three Artifical Clone Center")

library(distances)
library(igraph)

# clone size
clone_size_list = c(5,25,50,100)
coherence_list = c(0,1,2)
spatial_coord = spatial_coords[,c(2,3)]
dist <- as.matrix(distances::distances(as.matrix(spatial_coord)))
rownames(dist) <- rownames(spatial_coord);colnames(dist) <- rownames(spatial_coord)
celltype = celltype0[,2];names(celltype) = celltype0[,1]

clone_size=5;center=center1;clone_type="BE";coherence=1;neighbor_max=6;

#####
# This function finds clone neighbor given spatial coordinates and clone center
# input params:
# @spatial_coords: spot x spatial coords dataframe
# @ dist: Eucledian distance matrix (spot x spot). If NULL, the distance matrix will be auto-calculated.
# @center: (x,y) of clone center or barcode name
# @clone_size: size of the clones to be simulated
# @coherence: coherence of clones. 0 - nearest neighbors; 1- skip 1 nearest neighbors; 2- skip 2 nearest nrighbors
# @ neighbor_max: when coherence > 0, maximum neighbors allowed for each degree of nearest neighbors
# @celltype: a named vector contianing celltype of each spot. If NULL, 
#            directly choose nearest neighbors spots.
# @clone_type: a string of celltype of the spots that you want the clone to be generated in.
# output:
#   nn - clone members, including the input clone center 
# import distances,igraph
#####
find_clone_neighbor <- function(spatial_coord,dist=NULL,center,clone_size=50,coherence=0,neighbor_max=6,
                                celltype=NULL,clone_type="BE"){
  # calculate spot spatial distances to each other
  if(is.null(dist)){
    dist <- as.matrix(distances(as.matrix(spatial_coord)))
    rownames(dist) <- rownames(spatial_coord);colnames(dist) <- rownames(spatial_coord)
  }
  # subset if celltype is given
  if(!is.null(celltype)){
    dist <- dist[celltype ==clone_type,celltype ==clone_type]
  }
  
  # coherence = 0 - direct nearest neighbors
  # find nearest neighbors of each clone
  if(coherence == 0){
    nn = names(sort(dist[center,])[1:clone_size])
  }else{ # coherence > 0
    nn = c(center)
    history_nn = c(center)
    spots = center
    i = 1
    while(length(nn) < clone_size){
      if(i==1){ # first layer
        cur_layer_nn = names(sort(dist[spots,])[2:neighbor_max])
      }else{
        #print("Yes")
        cur_layer_nn = unique(as.vector(sapply(spots,function(x){
          #print(x)
          names(sort(dist[x,])[2:(neighbor_max+1)])})))
        cur_layer_nn = setdiff(cur_layer_nn,history_nn)
      }
      
      if(i%%(coherence+1) == 0){
        #print(i)
        if(length(unique(append(setdiff(cur_layer_nn,history_nn), nn))) > clone_size){
          cur_layer_nn = sample(cur_layer_nn,size=(clone_size-length(nn)))
        }
        nn = unique(append(setdiff(cur_layer_nn,history_nn), nn))
      }
      spots = setdiff(cur_layer_nn,history_nn)
      history_nn = unique(append(cur_layer_nn, history_nn))
      i=i+1
    }
  }
  return(nn)
}

# Simualte spatial coordinates of the clones
nn_list = list()
temp_plot_df = plot_df
centers = c(center1,center2,center3)
for(clone_size in clone_size_list){
    for(coherence in coherence_list){
      temp_name=paste0("CloneSize_",clone_size,"_Coherence_",coherence)
      temp_plot_df[,temp_name] = NA
      for(j in 1:length(centers)){
        center = centers[j]
        nn=find_clone_neighbor(spatial_coord,dist=dist,center=center,clone_size=clone_size,
                             coherence=coherence,neighbor_max=6,
                          celltype=celltype,clone_type="BE")
        nn_list[[paste0("CloneSize_",clone_size)]][[paste0("Coherence_",coherence)]][[paste0("Center_",center)]] = nn
        temp_plot_df[nn,temp_name] = paste0(j)
      }
  }
}
saveRDS(nn_list,"scripts/simulation_experiment/results/nnlist_A200S6.rds")

g_list = list()
for(coherence in coherence_list){
  for(clone_size in clone_size_list){
    for(j in 1:length(centers)){
      temp_name = paste0("CloneSize_",clone_size,"_Coherence_",coherence)
      print(temp_name)
      temp_df = temp_plot_df[!is.na(temp_plot_df[,temp_name]),]
      
      p <- ggplot(temp_plot_df,aes(x=-spatial_2,y=-spatial_1))+
        geom_point(size=0.1,color="Grey") + theme_bw()+ coord_fixed(ratio=1)+
        geom_point(data=temp_df,aes(x=-spatial_2,y=-spatial_1,color=.data[[temp_name]],),size=0.5) + 
        ggtitle(temp_name) + xlab("") + ylab("")
      g_list[[temp_name]] <- p
    }
  }
}
pdf("scripts/simulation_experiment/results/plot_df_withClones_A200S6.pdf",height=15,width=20)
gridExtra::grid.arrange(grobs = g_list,ncol=4,padding = unit(0, "line"))
dev.off()

saveRDS(temp_plot_df,"scripts/simulation_experiment/results/plot_df_withClones_A200S6.rds")

temp_plot_df = readRDS("scripts/simulation_experiment/results/plot_df_withClones_A200S6.rds")

## Now, simulate the VAF distribution
informative_vars = 15
var_sim = sample(real_vaf,size=informative_vars *length(centers))
vaf_sim = matrix(0,nrow=dim(N)[2], ncol=informative_vars *length(centers))
N_sim =  matrix(0,nrow=dim(N)[2], ncol=informative_vars *length(centers))
rownames(vaf_sim) = colnames(N);colnames(vaf_sim) = paste0(var_sim)
rownames(N_sim) = colnames(N);colnames(N_sim) = paste0(var_sim)
clone_size=50;coherence=1
temp_clone = temp_plot_df[,paste0("CloneSize_",clone_size,"_Coherence_",coherence)]
temp_clone[is.na(temp_clone)] = 0
p_j_list = c()
for(j in 1:informative_vars){
  vaf_j = var_sim[j]
  p_j = sample(temp_vaf,size=1)
  p_j_list = append(p_j,p_j_list)
  N_j = N[vaf_j,]*2
  N_j[N_j == 0] = rpois(n=sum(N_j == 0),lambda=ceiling(mean(N_j)))
  N_sim[,j] = N_j
  vaf_sim[temp_clone == 1,j] = rbinom(size=N_j,n=clone_size,prob=p_j)/(N_j[temp_clone == 1]+0.01)
  
}
for(j in (informative_vars+1):(2*informative_vars)){
  vaf_j = var_sim[j]
  p_j = sample(temp_vaf,size=1)
  p_j_list = append(p_j,p_j_list)
  N_j = N[vaf_j,]*2
  N_j[N_j == 0] = rpois(n=sum(N_j == 0),lambda=ceiling(mean(N_j)))
  N_sim[,j] = N_j
  vaf_sim[temp_clone == 2,j] = rbinom(size=N_j,n=clone_size,prob=p_j)/(N_j[temp_clone == 2]+0.01)
}
for(j in (2*informative_vars+1):(3*informative_vars)){
  vaf_j = var_sim[j]
  p_j = sample(temp_vaf,size=1)
  p_j_list = append(p_j,p_j_list)
  N_j = N[vaf_j,]*2
  N_j[N_j == 0] = rpois(n=sum(N_j == 0),lambda=ceiling(mean(N_j)))
  N_sim[,j] = N_j
  vaf_sim[temp_clone == 3,j] = rbinom(size=N_j,n=clone_size,prob=p_j)/(N_j[temp_clone == 3]+0.01)
}

save(var_sim,vaf_sim,N_sim,p_j_list,clone_size,file = paste0("scripts/simulation_experiment/results/TWIST_clone_",clone_size,"_coherence_",coherence,".RData"))
#png("results/Sydney_Bracht/Simulation/Simulated_Clones_beta.png",width=800,height=1000)

pdf("scripts/simulation_experiment/results/Simulated_Clones_TWIST.pdf",width=8,height=5)
hist(p_j_list,main="Histogram of Sampled P_j from TWIST",breaks=50)
anno_row=as.data.frame(temp_plot_df[temp_clone !=0,paste0("CloneSize_",clone_size,"_Coherence_",coherence)])
rownames(anno) = rownames(temp_plot_df[temp_clone !=0,]);colnames(anno)="Clone"
anno_col = as.data.frame(rep(c(1,2,3),each=informative_vars));
rownames(anno_col) = var_sim;colnames(anno_col)="Clone"
pheatmap(vaf_sim[temp_clone != 0,],annotation_row = anno,annotation_col=anno_col,show_rownames=F,main="VAF")
pheatmap(N_sim[temp_clone != 0,],annotation_row = anno, annotation_col=anno_col,show_rownames=F,main="Coverage N")
#print(p1)
#print(p2)
#grid.arrange(p1[[4]],p2[[4]])
dev.off()

## Now, simulate the VAF distribution
informative_vars = 15
var_sim = sample(real_vaf,size=informative_vars *length(centers))
vaf_sim = matrix(0,nrow=dim(N)[2], ncol=informative_vars *length(centers))
N_sim =  matrix(0,nrow=dim(N)[2], ncol=informative_vars *length(centers))
rownames(vaf_sim) = colnames(N);colnames(vaf_sim) = paste0(var_sim)
rownames(N_sim) = colnames(N);colnames(N_sim) = paste0(var_sim)
clone_size=50;coherence=1
temp_clone = temp_plot_df[,paste0("CloneSize_",clone_size,"_Coherence_",coherence)]
temp_clone[is.na(temp_clone)] = 0
p_j_list = c()
for(j in 1:informative_vars){
  vaf_j = var_sim[j]
  p_j = rbeta(shape1 = 0.5,shape2=0.5,n=1)
  p_j_list <- append(p_j,p_j_list )
  N_j = N[vaf_j,]*2
  N_j[N_j == 0] = rpois(n=sum(N_j == 0),lambda=ceiling(mean(N_j)))
  N_sim[,j] = N_j
  vaf_sim[temp_clone == 1,j] = rbinom(size=N_j[temp_clone == 1],n=clone_size,prob=p_j)/(N_j[temp_clone == 1]+0.01)
  
}
for(j in (informative_vars+1):(2*informative_vars)){
  vaf_j = var_sim[j]
  p_j = rbeta(shape1 = 0.5,shape2=0.5,n=1)
  p_j_list <- append(p_j,p_j_list )
  N_j = N[vaf_j,]*2
  N_j[N_j == 0] = rpois(n=sum(N_j == 0),lambda=ceiling(mean(N_j)))
  N_sim[,j] = N_j
  vaf_sim[temp_clone == 2,j] = rbinom(size=N_j[temp_clone == 2],n=clone_size,prob=p_j)/(N_j[temp_clone == 2]+0.01)
}

for(j in (2*informative_vars+1):(3*informative_vars)){
  vaf_j = var_sim[j]
  p_j = rbeta(shape1 = 0.5,shape2=0.5,n=1)
  p_j_list <- append(p_j,p_j_list )
  N_j = N[vaf_j,]*2
  N_j[N_j == 0] = rpois(n=sum(N_j == 0),lambda=ceiling(mean(N_j)))
  N_sim[,j] = N_j
  vaf_sim[temp_clone == 3,j] = rbinom(size=N_j[temp_clone == 3],n=clone_size,prob=p_j)/(N_j[temp_clone == 3]+0.01)
}

save(var_sim,vaf_sim,N_sim,p_j_list,clone_size,file = paste0("scripts/simulation_experiment/results/beta_distribution_clone_",clone_size,"_coherence_",coherence,".RData"))
#png("results/Sydney_Bracht/Simulation/Simulated_Clones_beta.png",width=800,height=1000)
pdf("scripts/simulation_experiment/results/Simulated_Clones_beta.pdf",width=8,height=5)
hist(p_j_list,main="Histogram of Sampled P_j from Beta Distribution(a=0.5,beta=0.5)",breaks=20)
anno_row=as.data.frame(temp_plot_df[temp_clone !=0,paste0("CloneSize_",clone_size,"_Coherence_",coherence)])
rownames(anno) = rownames(temp_plot_df[temp_clone !=0,]);colnames(anno)="Clone"
anno_col = as.data.frame(rep(c(1,2,3),each=informative_vars));
rownames(anno_col) = var_sim;colnames(anno_col)="Clone"
pheatmap(vaf_sim[temp_clone != 0,],annotation_row = anno,annotation_col=anno_col,show_rownames=F,main="VAF")
pheatmap(N_sim[temp_clone != 0,],annotation_row = anno, annotation_col=anno_col,show_rownames=F,main="Coverage N")
#print(p1)
#print(p2)
#grid.arrange(p1[[4]],p2[[4]])
dev.off()
