library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(Matrix)
library(ggrastr)
library(RColorBrewer)
library(tidyr)

rm(list=ls())
gc()

#MAESTER source code from Miller et al. 
source("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/210215_FunctionsGeneral.R")

## RCTD  ---------------------------------------------------------------------------------------------------------------------
#reference generation 

#loading in data
seu <- Read10X('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/output/CRC076/filtered_feature_bc_matrix')
seu <- CreateSeuratObject(seu)
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/maegatk.rds')
#using paired MAESTER data to isolate high quality cells, reduce ambients from snRNA-seq

# Generate cell IDs that will match Maegatk
seu$cellMerge <- colnames(seu)

# Only keep cells with a cellMerge id that occurs once, then intersect with Maester data
cells.ch <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell %>% unname
cells.ch <- cells.ch[cells.ch %in% colnames(maegatk.rse)]

# Subset for common cells
seu <- subset(seu, subset = cellMerge %in% cells.ch)
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL


seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000) 
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 4000)

# Calculate cell cycle scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu <- CellCycleScoring(seu, 
                        s.features = s.genes, 
                        g2m.features = g2m.genes, 
                        set.ident = FALSE)

#visualization optional"
RidgePlot(seu, features = c("S.Score", "G2M.Score"), ncol = 2)

all.genes <- rownames(seu)
seu <- ScaleData(seu, 
                 features = all.genes, 
                 vars.to.regress = c("S.Score", "G2M.Score"))

seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = FALSE)
set.seed(101)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.1)
seu <- RunUMAP(seu, dims = 1:20, min.dist = .2, n_neighbors = 50, spread = 2)
DimPlot(seu, reduction = 'umap')

# You can compare the clustering with and without cell cycle influence
DimPlot(seu, reduction = 'umap', group.by = "Phase")

seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]

markers <- FindAllMarkers(seu, min.pct = .25)
write.csv(markers, '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/analysis/markers.csv')

cluster.ids <-  c("Colorectal cancer", #epithelial
                  "Tumor-associated stroma", 
                  "Immune", #myeloid
                  "Endothelial", 
                  "Fibroblasts", #CAFs?
                  "Immune", #B-cells
                  "Immune") #dendritic cells

names(cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, cluster.ids)

cluster_colors <- c(
  "Colorectal cancer" = "#C6D3BD",  
  "Tumor-associated stroma" = "#751385",   
  "Endothelial" = "red",                         
  "Immune" = "#19368a",  
  "Fibroblasts" = "gold"
)
plt <- DimPlot(seu, reduction = 'umap', cols = cluster_colors)

pdf('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/analysis/UMAP_annotation.pdf', width = 8, height = 8)
print(plt)
dev.off()
saveRDS(seu, file = "CRC076_sn_seurat_final.rds")


library(spacexr)
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/rctd_ref')
counts <- GetAssayData(object = seu, assay = "RNA", slot = "counts")
annotations <- seu@active.ident
write.csv(annotations,'CRC076_sn_celltype_annotations.csv')

counts <- as.data.frame(counts) # converting count matrix into data.frame as read.csv would do

meta_data <- read.csv("CRC076_sn_celltype_annotations.csv") # load in meta_data (barcodes, clusters)
cell_types <- meta_data$x; names(cell_types) <- colnames(counts) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type

### Create the Reference object
reference <- Reference(counts, cell_types,require_int = F)
saveRDS(reference, 'CRC076_sn_ref_RCTD.rds')

rm(reference)
rm(counts)

# Subset and order MAEGATK RSE object
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/analysis')
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/maegatk.rds')
maegatk.rse <- maegatk.rse[,colnames(seu)]

saveRDS(maegatk.rse, file = "CRC076_sn_maegatk_final.rds")
## -------------------------------------------------------------------------------------------------------------------

## Loading in raw spatial data ------------------------------------------------------------------------------------------------------------------
seu <- readRDS(file = "/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/seeker/run2/CRC076_C/CRC076_C_seurat.rds")

## Spatial transcriptomics QC  ------------------------------------------------------------------------------------------------------------------
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/seeker/run2/CRC076_C/Seurat/QC')

SpatialFeaturePlot(seu, 
                   features = "log_umi", 
                   pt.size = 1,
                   min.cutoff = "q1",
                   max.cutoff = "q99") +
  scale_fill_gradientn(colors = c("#efefef", "#CBC3E3", "#301934"))+
  labs(title = "UMI Counts per Cell") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right", 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.line = element_line(color = "black")
  )
ggsave("UMI_feature_plot.pdf", width = 10, height = 8)

plot_data <- data.frame(
  log10_UMIs = seu@meta.data[["log_umi"]]
)

p <- ggplot(plot_data, aes(x = log10_UMIs)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = NA, alpha = 0.7) +
  labs(x = expression(log[10]~"UMIs per feature"),
       y = "Number of features",
       title = "Distribution of UMIs per Feature") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"), 
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    axis.line = element_line(color = "black")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_x_continuous(limits = c(0, NA)) + 
  ylim(0,2000) + 
  xlim(0,5) 

ggsave("UMI_distribution_histogram.pdf", p, width = 10, height = 8)

## Filtering by RNA count and gene count for high quality cells
Idents(seu) <- seu$orig.ident
p1 <- VlnPlot(seu, features = c('nCount_RNA', 'nFeature_RNA', 'log_umi'))
pdf('QC_violinplot_CRC076.pdf')
print(p1)
dev.off()

seu[["percent.mt.human"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
p2 <- VlnPlot(seu, feature = c('percent.mt.human'))
pdf('mito_violinplot_CRC076.pdf')
print(p2)
dev.off()


## RCTD cell type annotation ---------------------------------------------------------------------------------------------------------------------
analysis_folder <-'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/RCTD'

args=c("CRC076_C","/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/seeker/run2/CRC076_C",'colon')
slide_name <- args[1] 
print(slide_name)
data_path <- args[2]
print(data_path)
ref_data= args[3] 
print(ref_data)

# load reference data
if(ref_data=="colon"){
  reference <-readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/rctd_ref/CRC076_sn_ref_RCTD.rds') 
}

# output directory
dir_path <- analysis_folder
print(dir_path)
dir.create(dir_path)

# load sample spatial data
count_path=data_path
mtx <- readMM(paste0(count_path,"/",slide_name,"_MoleculesPerMatchedBead.mtx"))
barcodes = as.data.frame(read.table(paste0(count_path,"/",slide_name,"_barcodes.tsv")))
features=read.table(paste0(count_path,"/",slide_name,"_genes.tsv"),stringsAsFactors = F)

rownames(mtx) <- features[,1]
colnames(mtx) <- barcodes[,1]
print("mtx:")
print(dim(mtx))

nUMI <- colSums(mtx) # In this case, total counts per pixel is nUMI

if(file.exists(paste0(count_path,"/",slide_name,"_coordinates.csv"))){
  coords = read.table(paste0(count_path,"/",slide_name,"_coordinates.csv"),sep=",",header=T)
  colnames(coords) <-c("barcodes","xcoord","ycoord")
  rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
}else{
  coords = read.table(paste0(count_path,"/",slide_name,"_MatchedBeadLocation.csv"),sep=",",header=T)
  colnames(coords) <-c("barcodes","xcoord","ycoord")
  rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
  coords = coords[colnames(mtx),]
}
print("spot:")
print(dim(coords))

x = coords$xcoord 
y = coords$ycoord

### Create SpatialRNA object
puck <- SpatialRNA(coords, mtx, nUMI)

## Examine SpatialRNA object (optional)
print(dim(puck@counts)) # observe Digital Gene Expression matrix
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/RCTD')
pdf(paste0(dir_path,"diagnostic_stats.pdf"),width=4,height=4)
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

dev.off()

myRCTD <- create.RCTD(puck, reference, max_cores = 1, UMI_min = 25)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
results$weights_doublet = normalize_weights(results$weights_doublet) 

cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- dir_path

# Plotting
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plots the number of confident pixels of each cell type in 'doublet_mode' (meaning two cell types max should be assigned to each bead)
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
# save the weights
saveRDS(myRCTD,paste0(dir_path,"/RCTD_res.rds"))

#for celltypeID
data <- data.frame(results[["results_df"]],results[["weights_doublet"]]) 
data <- data %>% rename(first_type.1 = 'type1_score', second_type.1 = 'type2_score')
data <- data %>% filter(spot_class != 'reject') %>% select(-spot_class)

data$cell_type <- ifelse(data$type1_score > 0.7, as.character(data$first_type),
                         ifelse(data$type2_score > 0.7, as.character(data$second_type),
                                "Undetermined"))

annotations <- data %>% select(first_type)
seu <- seu[,colnames(seu) %in% rownames(data)]
seu$annotation <- annotations
Idents(seu) <- seu$annotation

cluster_colors <- c(
  "Colorectal cancer" = "#C6D3BD",  
  "Tumor-associated stroma" = "#751385",   
  "Endothelial" = "red",                         
  "Immune" = "#19368a",  
  "Fibroblasts" = "gold"
)
spatial_plt <- SpatialDimPlot(seu, cols = cluster_colors, pt.size.factor = 1.3)
spatial_plt

pdf('spatial_RCTD_annotations.pdf')
print(spatial_plt)
dev.off()

saveRDS(seu, 'CRC076_seurat_RCTD.rds')
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## MAESTER preprocessing  ------------------------------------------------------------------------------------------------------------------
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/rds')

# Loading in data
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/maegatk_output/maegatk.rds')

# Only keep cells with a cellMerge id that occurs once, then intersect with Maester data
seu$cellMerge <- paste0(cutf(colnames(seu), d = "-", f = 1), "-1")
cells.ch <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell %>% unname
cells.ch <- cells.ch[cells.ch %in% colnames(maegatk.rse)]

seu <- subset(seu, subset = cellMerge %in% cells.ch)
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL

seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]

seu$SPATIAL_1 <-seu@reductions$SPATIAL@cell.embeddings[,1]
seu$SPATIAL_2 <-seu@reductions$SPATIAL@cell.embeddings[,2]

maegatk.rse <- maegatk.rse[,colnames(seu)]

saveRDS(seu, 'CRC076C_seurat_final.rds' )
saveRDS(maegatk.rse, 'CRC076C_maegatk_final.rds')

# Prepare allele frequency matrix
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/pipeline_output')
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100
write.csv(af.dm, file = "af_dm.csv")

# Extract cell IDs
cells.tib <- as_tibble(seu@meta.data, rownames = "cell")
CellSubset.ls <- list(unionCells = cells.tib$cell)

# Get the mean allele frequency and coverage for every cell subset
mean_af.ls <- lapply(CellSubset.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubset.ls, function(x) rowMeans(assays(maegatk.rse)[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset
quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)

library(stringr)
start_time <- Sys.time()
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubset.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
Sys.time() - start_time

# Get the mean quality for each variant. This can take a few hours.
start_time <- Sys.time()
qual.num <- sapply(rownames(af.dm), function(x) {
  pos <- as.numeric( cutf(x, d = "_") )
  message(pos)
  mut <- cutf(x, d = ">", f = 2)
  # Only use cells in which the base was sequenced. Use reverse only because that's how we amplify transcripts.
  covered <- assays(maegatk.rse)[[str_c(mut, "_counts_rev")]][pos,] > 0
  # Get the mean quality for this call
  qual <- mean( assays(maegatk.rse)[[str_c(mut, "_qual_rev")]][pos, covered] )
  return(qual)
})
Sys.time() - start_time

library(tibble)
library(readr)
# Collect all information in a tibble
vars.tib <- as_tibble(do.call(cbind, c(mean_af.ls, mean_cov.ls, unlist(quantiles.ls, recursive = F))), rownames = "var")
vars.tib <- add_column(vars.tib, quality = qual.num, .before = 2)

# Save for fast loading next time
write_tsv(vars.tib, "CRC076C_allvariants.txt")
vars.tib <- read.table('CRC076C_allvariants.txt', header = TRUE)

# Parameters to select clones with high VAF
voi.ch <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 20) %>%
  filter(q99.unionCells > 20) %>% .$var
voi.ch

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

write.table(voi.ch, file = "CRC076_voi_vaf.tsv", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)

TWIST_variants <- c("9982_G>T","4239_C>A", "5380_T>A","3957_A>C", "7411_A>C",  "11892_A>C")

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Lineage AF plots  ------------------------------------------------------------------------------------------------------------------
for (v in TWIST_variants) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  cells.tib$af_voi[cells.tib$cov_voi < 3] <- NA
  
  p <- ggplot(data = cells.tib, aes(x = SPATIAL_1, y = SPATIAL_2, color = af_voi))+
    geom_point(size = .25) +  
    ggtitle(v) + 
    scale_color_gradientn(colors = c("#c8dfea", "#608da2", "darkblue"), limits = c(0.01, 100), na.value = "#ECECEC") + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) + 
    theme_bw() +
    theme(panel.grid = element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5)) +   
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank()) 
  
  name <- paste0(v, ".pdf")
  ggsave(name, plot = p, width = 7, height = 7)
  
  cells.tib$af_voi <- NULL
  cells.tib$cov_voi <- NULL
}
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Mitochondrial QC: coverage  ----------------------------------------------------------------------
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/coverage')
# Set y axis parameters
ymax <- 200

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

# Plot coverage along chrM
base.tib <- tibble(base = 1:16569, depth = rowMeans(assays(maegatk.rse)[["coverage"]]))

coverage <- 
  base.tib %>% ggplot() +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-ymax*0.2, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
  ylab("Mean coverage per cell") + xlab("Position along chrM") +  ggtitle('CRC076, colorectal cancer') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))

ggsave('CRC076_MAESTER_coverage.png', plot = coverage)

#comparing to Curio Seeker data: 

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/curio_compare')
maegatk.curio <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/curio_compare/CRC076_curio.mr1/maegatk.rds')
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/single_nuclei_references/MAESTER/output/CRC076/maegatk.rds')

#subsetting so the two have the same cells
int <- intersect(colnames(maegatk.curio), colnames(maegatk.rse))
maegatk.curio <- maegatk.curio[, int]
maegatk.rse <- maegatk.rse[, int]

#make sure theyre the same
dim(maegatk.curio)
dim(maegatk.rse)

# Plot mean coverage per position per cell ----------------------------------------------------------------------
# Set y axis parameters
ymax <- 1000000

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

# Plot coverage along chrM
base.tib1 <- tibble(base = 1:16569, depth = rowSums(assays(maegatk.curio)[["coverage"]]))
base.tib2 <- tibble(base = 1:16569, depth = rowSums(assays(maegatk.rse)[["coverage"]]))

coverage_curio <- 
  base.tib1 %>% ggplot() +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "gold", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  ylab("Mean coverage per cell") + xlab("Position along chrM") + ggtitle('MT before enrichment') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))

coverage_curio
pdf('Mito_coverage_curio.pdf')
print(coverage_curio)
dev.off()

coverage_maester <- 
  base.tib2 %>% ggplot() +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  ylab("Mean coverage per cell") + xlab("Position along chrM") + ggtitle('MT after enrichment') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))

coverage_maester
pdf('mito_coverage_MAESTER.pdf')
print(coverage_maester)
dev.off()

combined_data <- rbind(transform(base.tib1, dataset = "Slide-seq"), 
                       transform(base.tib2, dataset = "MAESTER"))

ratio_plot<-
  ggplot(combined_data) +
  geom_boxplot(aes(y = depth, fill = dataset)) + 
  scale_y_continuous(trans = "log10") + 
  ylab('Mean coverage per cell') + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("#64b53b", "gold"))

ggsave("Mito_boxplot.pdf",ratio_plot)

combined_coverage <- 
  base.tib1 %>% ggplot() +
  geom_bar(data = base.tib2,
           aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), 
           stat = "identity", fill = "#64b53b", width = 1) +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), 
           stat = "identity", fill = "gold", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  ylab("Mean coverage per cell") + 
  xlab("Position along chrM") + 
  ggtitle('Coverage comparison: Before and After enrichment') + 
  theme_classic() +
  theme(aspect.ratio = 0.5,
        plot.title = element_text(hjust = 0.5))

combined_coverage
pdf('Mito_coverage_compare.pdf')
print(combined_coverage)
dev.off()
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## mtDNA sequencing confirmation of variants -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/TWIST_validation/cancer')

cancer <- read.csv("/Users/sydneybracht/Downloads/mitochondrial_base_counts.csv")
mutations <- c(rownames(af.dm))
reference <- character(length(mutations))

for (i in 1:length(mutations)) {
  # Split the string at the underscore
  parts <- strsplit(mutations[i], "_", fixed = TRUE)[[1]]
  
  # If there's a second part, get its first character
  if (length(parts) >= 2) {
    reference[i] <- substr(parts[2], 1, 1)
  } else {
    reference[i] <- NA
  }
}

# Extract positions from mutation names
positions <- numeric(length(mutations))
for (i in 1:length(mutations)) {
  parts <- strsplit(mutations[i], "_", fixed = TRUE)[[1]]
  if (length(parts) >= 1) {
    positions[i] <- as.numeric(parts[1])
  } else {
    positions[i] <- NA
  }
}

variants <- character(length(mutations))
for (i in 1:length(mutations)) {
  parts1 <- strsplit(mutations[i], "_", fixed = TRUE)[[1]]
  if (length(parts1) >= 2) {
    after_underscore <- parts1[2]
    parts2 <- strsplit(after_underscore, ">", fixed = TRUE)[[1]]
    if (length(parts2) >= 2) {
      variants[i] <- parts2[2]
    } else {
      variants[i] <- NA
    }
  } else {
    variants[i] <- NA
  }
}

known_variants <- data.frame(
  Position = positions,
  RefNucleotide = reference,
  VarNucleotide = variants,
  stringsAsFactors = FALSE
)

# Convert the data from wide to long format for the percent columns
allele_freq <- cancer %>%
  select(Position, A_percent, C_percent, G_percent, T_percent) %>%
  pivot_longer(
    cols = c(A_percent, C_percent, G_percent, T_percent),
    names_to = "Nucleotide",
    values_to = "Percent"
  ) %>%
  # Clean up the nucleotide names by removing the "_percent" suffix
  mutate(Nucleotide = gsub("_percent", "", Nucleotide))

# Join with the known variants data to identify reference vs variant alleles
variant_data <- allele_freq %>%
  inner_join(known_variants, by = "Position") %>%
  mutate(
    Type = case_when(
      Nucleotide == RefNucleotide ~ "Reference",
      Nucleotide == VarNucleotide ~ "Variant",
      TRUE ~ "Other"
    )
  )


# Extract only the variant alleles for our analysis
variant_alleles <- variant_data %>%
  filter(Type == "Variant") %>%
  # For very small percentages, add a small value to make them visible on log scale
  mutate(Percent = ifelse(Percent == 0, 0.01, Percent))

p1 <- ggplot(variant_alleles, aes(x = Percent)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(
    title = "Distribution of Variant Allele Frequencies in Mitochondrial Genome",
    x = "Allele Frequency (%)",
    y = "Log(count)"
  ) +
  scale_y_continuous(trans = "log10") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

pdf('logcount_allelefreq.pdf')
print(p1)
dev.off()

p2 <- ggplot(variant_alleles, aes(x = Percent)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black") +
  scale_x_log10(
    breaks = c(0.1, 0.5, 1, 5, 10, 25, 50, 100),
    labels = c("0.1", "0.5", "1", "5", "10", "25", "50", "100")
  ) +
  annotation_logticks(sides = "b") +
  labs(
    title = "Log-Scale Distribution of Variant Allele Frequencies",
    x = "Allele Frequency (%) - Log Scale",
    y = "Count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

pdf('count_logallelefreq.pdf')
print(p2)
dev.off()

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Plots from Caleb Lareau's critcism of ReDeem -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Parameters to select clones with low thresholds to plot noise vs. real
voi.ch <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 5) %>%
  filter(q99.unionCells >= 5) %>% .$var

# plot 1: # of cells with the mutation vs. # of reads per molecule per cell thats part of the lineage, coloring the variants of interest. 
af.dm_mini <- af.dm[voi.ch,]
#First, from af.dm calculating how many cells there are with a given mutation (threshold is AF > 20%)

cells_per_lineage <- data.frame(
  Variant = rownames(af.dm_mini),
  NumCells = rowSums(af.dm_mini > 15, na.rm = TRUE)
)

#subsetting to only plot variants that have any cells.
has_cells <- cells_per_lineage %>% filter(cells_per_lineage$NumCells > 1)
af.dm_mini <- af.dm_mini[rownames(has_cells),]
af.dm_mini[af.dm_mini <10 ] <- NA #unreliable low VAF

#getting number of reads per cell for each variant.
coverage_info <- as_tibble(colnames(af.dm_mini)) %>% rename(value = 'barcode')
for (v in rownames(af.dm_mini)) {
  message(v)
  col_name <- v
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),coverage_info$barcode]
}

# If that cell is NA in the allele frequency matrix (aka, not part of the lineage), the coverage becomes NA. not included in the average coverage calcualtions. 
for (v in rownames(af.dm_mini)) {
  # Skip the 'value' or 'barcode' column in coverage_info
  if (v %in% colnames(coverage_info) && v != 'barcode') {
    # Get the NAs mask from af.dm for this variant
    na_mask <- is.na(af.dm_mini[v, ])
    
    # Apply the mask to coverage_info - set to NA where af.dm is NA
    coverage_info[[v]][na_mask] <- NA
  }
}
coverage_info <- coverage_info %>% select(-barcode)

af.dm_mini <- af.dm_mini[colnames(coverage_info),]

all_variants <- rownames(af.dm_mini)
spatial_variants <- TWIST_variants

plot_data <- data.frame(
  Variant = colnames(coverage_info),
  meanCoverage = colMeans(coverage_info, na.rm = TRUE), 
  cellsperlienage = rowSums(af.dm_mini > 15, na.rm = TRUE)
)

variants <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 20) %>%
  filter(q99.unionCells > 20) %>% .$var

plot_data$annotation <- "Detected"
plot_data$annotation[plot_data$Variant %in% variants] <- "Signfificant by maegatk"
plot_data$annotation[plot_data$Variant %in% spatial_variants] <- "Spatial Variant"

# Convert annotation to factor with specific levels for proper ordering in legend
plot_data$annotation <- factor(plot_data$annotation, levels = c("Detected", "Signfificant by maegatk", "Spatial Variant"))

# Then use that column for coloring
p1 <- ggplot(plot_data, aes(x = cellsperlienage, y = meanCoverage, color = annotation)) +
  geom_point(size = 1.25) + 
  geom_hline(yintercept = 3, linetype = 2) +
  scale_x_log10() + 
  scale_y_log10() + 
  xlab('Cells with variant') + 
  ylab('Mean reads per cell with variant') +
  scale_color_manual(values = c("Detected" = "lightgrey", "Signfificant by maegatk" = "blue", "Spatial Variant" = "red")) +
  theme_classic()
p1

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/artifact_analysis')
pdf('lineagesize_vs_coverage.pdf')
print(p1)
dev.off()

## Making a graph of mismatch vs. position in the read (is there end of sequence bias?)
counts_per_position <- read.csv('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/artifact_analysis/position_counts_simple.csv')
counts_per_position$Position <- (counts_per_position$Position) / 250 #getting fraction for plotting
positions <- data.frame(Variant = spatial_variants, position = c('226.0', '230.0', "234.8","180.0","227.0", "208.0"))

# Convert positions$position to numeric for comparison
numeric_positions <- as.numeric(as.character(positions$position)) / 250

counts_per_position$spatial <- FALSE

# Find the closest positions with some tolerance
tolerance <- 0.01  # Adjust this based on your data precision
for (pos in numeric_positions) {
  # Find positions within tolerance
  matches <- which(abs(counts_per_position$Position - pos) < tolerance)
  if (length(matches) > 0) {
    counts_per_position$spatial[matches] <- TRUE
  }
}

# Create the plot with conditional coloring
p <- ggplot(counts_per_position, aes(x = Position, y = Count, fill = spatial)) +
  geom_bar(stat = "identity", width = 0.015) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  scale_y_log10() +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "red"), guide = "none") +
  theme_classic() +
  labs(
    x = "Position in Read",
    y = "Count of Alternative Allele"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )
p

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/artifact_analysis')
pdf('position_of_mismatches.pdf')
print(p)
dev.off()

# Where are the variants compared to primer binding sites?
library(ggforce) # For curved text and arc drawing

# Define the mitochondrial genome size
mito_size <- 16569

# Define primer sites (from your list)
primer_sites <- c(400, 953, 637, 1848, 3679, 7561, 10264, 10496, 11900, 13926, 14263, 15643, 
                  3777, 5145, 6957, 7852, 8766, 9535, 10127, 11684, 13758, 14492, 15432, 
                  3537, 4923, 6742, 7609, 8541, 9316, 11451, 13515, 14664, 15260, 6563, 
                  9847, 11654, 7921, 10114, 10132, 8815, 10884, 4833, 6324, 11223, 13069, 
                  14937, 2524, 6124, 10994, 12831, 14789, 2110, 5910, 10761, 12601, 16791, 
                  2360, 13472, 9858, 9921, 12725, 7899, 2323, 1895, 8071)

# Define mitochondrial genome features (approximate positions based on the human mitochondrial genome)
# You may need to adjust these to match your specific mitochondrial genome
features <- data.frame(
  name = c("D-loop", "RNR1", "RNR2", "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"),
  start = c(16024, 648, 1671, 3307, 4470, 5904, 7586, 8366, 8527, 9207, 10059, 10470, 10760, 12337, 14149, 14747),
  end = c(576, 1601, 3229, 4262, 5511, 7445, 8269, 8572, 9207, 9990, 10404, 10766, 12137, 14148, 14673, 15887),
  type = c("Control", "rRNA", "rRNA", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding", "Coding")
)

# Define tRNA positions (these are approximate based on human mtDNA)
trnas <- data.frame(
  name = paste0("tRNA-", 1:22),
  position = c(577, 1602, 3230, 4263, 4329, 4400, 5512, 5587, 5761, 5826, 6900, 7446, 7518, 8295, 9991, 10405, 12138, 12207, 14674, 14742, 15888, 15953),
  type = "tRNA"
)

# Create a data frame for the primer sites
primers_df <- data.frame(
  position = primer_sites,
  type = "Primer"
)

# Define variant sites 
variants_raw <- TWIST_variants

# Extract positions from variant notation
variant_positions <- as.numeric(sapply(strsplit(variants_raw, "_"), function(x) x[1]))

# Extract mutation types
mutation_types <- sapply(strsplit(variants_raw, "_"), function(x) x[2])

# Create a data frame for variants
variants_df <- data.frame(
  position = variant_positions,
  mutation = mutation_types,
  type = "Variant"
)

# Function to convert positions to radians for plotting in a circle
pos_to_rad <- function(pos, genome_size = mito_size) {
  # Convert to radians, starting from the top (π/2) and going clockwise
  rads <- 2 * pi * (pos / genome_size)
  return(pi/2 - rads)  # Adjust to start from top and go clockwise
}

# Create data frame for the outer circle (genome)
circle_points <- data.frame(
  angle = seq(0, 2*pi, length.out = 500),
  radius = 1
)
circle_x <- circle_points$radius * cos(circle_points$angle)
circle_y <- circle_points$radius * sin(circle_points$angle)
circle_df <- data.frame(x = circle_x, y = circle_y)

# Create data frame for feature arcs
features_plot <- data.frame()
for(i in 1:nrow(features)) {
  # Handle features that cross the origin (like D-loop)
  if(features$start[i] > features$end[i]) {
    # First part: from start to end of genome
    start_angle <- pos_to_rad(features$start[i])
    end_angle <- pos_to_rad(mito_size)
    
    # Add points for the first arc
    arc_angles1 <- seq(start_angle, end_angle, length.out = 100)
    radius <- 0.8  # Slightly inside the main circle
    x1 <- radius * cos(arc_angles1)
    y1 <- radius * sin(arc_angles1)
    temp_df1 <- data.frame(
      x = x1, y = y1, 
      name = features$name[i],
      type = features$type[i]
    )
    
    # Second part: from start of genome to end position
    start_angle <- pos_to_rad(1)
    end_angle <- pos_to_rad(features$end[i])
    
    # Add points for the second arc
    arc_angles2 <- seq(start_angle, end_angle, length.out = 100)
    x2 <- radius * cos(arc_angles2)
    y2 <- radius * sin(arc_angles2)
    temp_df2 <- data.frame(
      x = x2, y = y2, 
      name = features$name[i],
      type = features$type[i]
    )
    
    features_plot <- rbind(features_plot, temp_df1, temp_df2)
  } else {
    # Regular feature that doesn't cross the origin
    start_angle <- pos_to_rad(features$start[i])
    end_angle <- pos_to_rad(features$end[i])
    
    # Add points for the arc
    arc_angles <- seq(start_angle, end_angle, length.out = 100)
    radius <- 0.8  # Slightly inside the main circle
    x <- radius * cos(arc_angles)
    y <- radius * sin(arc_angles)
    temp_df <- data.frame(
      x = x, y = y, 
      name = features$name[i],
      type = features$type[i]
    )
    
    features_plot <- rbind(features_plot, temp_df)
  }
}

# Convert tRNA positions to coordinates
trna_coords <- data.frame(
  x = 0.85 * cos(pos_to_rad(trnas$position)),
  y = 0.85 * sin(pos_to_rad(trnas$position)),
  name = trnas$name,
  type = trnas$type
)

# Convert primer site positions to coordinates for the outer spikes
primer_coords <- data.frame(
  x_inner = 1.05 * cos(pos_to_rad(primers_df$position)),
  y_inner = 1.05 * sin(pos_to_rad(primers_df$position)),
  x_outer = 1.15 * cos(pos_to_rad(primers_df$position)),
  y_outer = 1.15 * sin(pos_to_rad(primers_df$position)),
  name = paste0("Primer_", primers_df$position),
  type = primers_df$type
)

# Convert variant positions to coordinates for outer dashes (placed under the primer dashes)
variant_coords <- data.frame(
  x_inner = 1.05 * cos(pos_to_rad(variants_df$position)),
  y_inner = 1.05 * sin(pos_to_rad(variants_df$position)),
  x_outer = 1.15 * cos(pos_to_rad(variants_df$position)),
  y_outer = 1.15 * sin(pos_to_rad(variants_df$position)),
  name = paste0("Variant_", variants_df$position),
  mutation = variants_df$mutation,
  type = variants_df$type
)

# Create labels for genes at the middle of their arc
gene_labels <- data.frame()
for(i in 1:nrow(features)) {
  # Handle features that cross the origin
  if(features$start[i] > features$end[i]) {
    # For features crossing the origin, place label in the middle of the largest segment
    if((mito_size - features$start[i]) > features$end[i]) {
      # Label in the first segment
      mid_pos <- features$start[i] + (mito_size - features$start[i])/2
    } else {
      # Label in the second segment
      mid_pos <- features$end[i]/2
    }
  } else {
    # Regular feature
    mid_pos <- features$start[i] + (features$end[i] - features$start[i])/2
  }
  
  angle <- pos_to_rad(mid_pos)
  radius <- 0.8  # Place labels directly in the feature arcs
  
  # Calculate text angle to make it readable
  text_angle <- (angle * 180/pi) + 90
  # Adjust text angle to always be readable (not upside-down)
  if (text_angle > 90 && text_angle < 270) {
    text_angle <- text_angle + 180
  }
  
  temp_label <- data.frame(
    x = radius * cos(angle),
    y = radius * sin(angle),
    angle = text_angle, # Adjusted for text orientation
    name = features$name[i],
    type = features$type[i]
  )
  
  gene_labels <- rbind(gene_labels, temp_label)
}

# Create the plot
p <- ggplot() +
  # Add the main circle
  geom_path(data = circle_df, aes(x = x, y = y), color = "black", size = 0.5) +
  
  # Add the feature arcs with proper coloring
  geom_path(data = features_plot, 
            aes(x = x, y = y, group = name, color = type), 
            size = 8) +
  
  # Add tRNA markers
  geom_point(data = trna_coords, 
             aes(x = x, y = y, color = type), 
             size = 3, shape = 18) +
  
  # Add variant site markers first (so they appear under the primer markers)
  geom_segment(data = variant_coords,
               aes(x = x_inner, y = y_inner, xend = x_outer, yend = y_outer),
               color = "red", size = 1) +
  
  # Add primer site spikes on top
  geom_segment(data = primer_coords,
               aes(x = x_inner, y = y_inner, xend = x_outer, yend = y_outer),
               color = "green4", size = 1) +
  
  # Add gene labels
  geom_text(data = gene_labels,
            aes(x = x, y = y, label = name, angle = angle),
            size = 3, hjust = 0.5, vjust = 0.5, color = "black", fontface = "bold") +
  
  # Color scheme similar to the reference image
  scale_color_manual(values = c(
    "Coding" = "grey70",
    "rRNA" = "lightblue",
    "tRNA" = "magenta",
    "Control" = "gold",
    "Primer" = "green4"
  ), guide = "none") +  # Hide the default legend
  
  # Set theme for circular plot
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(ratio = 1) +
  labs(title = "Mitochondrial Genome Map")

# Legend modifications
p <- p + 
  guides(color = guide_legend(
    title = "Features",
    override.aes = list(
      shape = c(NA, NA, 18, NA),
      linetype = c(1, 1, 0, 0)
    )
  ))

# Create custom legend
legend_entries <- data.frame(
  x = rep(1.5, 6),
  y = seq(-1.1, -1.6, length.out = 6),
  type = c("Coding", "Control", "rRNA", "tRNA", "Primer", "Variant")
)

# Match colors with the main plot
legend_colors <- c(
  "Coding" = "grey70",
  "Control" = "gold", 
  "rRNA" = "lightblue",
  "tRNA" = "magenta",
  "Primer" = "green4",
  "Variant" = "red"
)

# Add squares for feature types
p <- p +
  geom_rect(data = legend_entries[1:4,],
            aes(xmin = x - 0.1, xmax = x + 0.1, 
                ymin = y - 0.03, ymax = y + 0.03,
                fill = type),
            color = "black") +
  geom_segment(data = legend_entries[5,],
               aes(x = x - 0.1, y = y, xend = x + 0.1, yend = y),
               color = "green4", size = 1) +
  geom_segment(data = legend_entries[6,],
               aes(x = x - 0.1, y = y, xend = x + 0.1, yend = y),
               color = "red", size = 1) +
  geom_text(data = legend_entries,
            aes(x = x + 0.15, y = y, label = type),
            hjust = 0, size = 3) +
  scale_fill_manual(values = legend_colors, guide = "none")

print(p)

# To save the plot
ggsave("mitochondrial_genome_map.pdf", p, width = 10, height = 10)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Pairwise scatterplots of AF and coverage for variants of interest.
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB042/seeker/CRC076/CRC076_C/MAESTER/pipeline_output')
## VAF plots for all pairs of variants colored by coverage conditionally

# Define your variants
variants_of_interest <- c('4239_C>A', '3957_A>C', '7411_A>C', '5380_T>A', '11892_A>C', '9982_G>T')

# Get coverage for all variants
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'barcode')

for (v in variants_of_interest) {
  message(v)
  col_name <- v
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric(cutf(v, d = "_")),coverage_info$barcode]
}

# Generate all unique pairs
variant_pairs <- combn(variants_of_interest, 2, simplify = FALSE)

# Create plots for each pair
plot_list <- list()

for (i in seq_along(variant_pairs)) {
  var1 <- variant_pairs[[i]][1]
  var2 <- variant_pairs[[i]][2]
  
  message(paste("Creating plot for", var1, "vs", var2))
  
  # Create plot data
  plot_data <- data.frame(
    af_1 = af.dm[var1,], 
    af_2 = af.dm[var2,], 
    coverage_var1 = coverage_info[[var1]],
    coverage_var2 = coverage_info[[var2]]
  )
  
  # Create grouping variable 
  plot_data$group <- ifelse(plot_data$af_2 > plot_data$af_1, "above", "below")
  
  # Split data for different color palettes
  data_above <- plot_data[plot_data$group == "above", ]
  data_below <- plot_data[plot_data$group == "below", ]
  
  # Create the plot with two different color palettes
  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.7) +
    # Points above y=x line (var1 coverage) - using magma palette
    geom_point(data = data_above, aes(x = af_1, y = af_2, color = coverage_var1, size = coverage_var1)) +
    scale_color_viridis_c(option = "magma", direction = -1, name = paste(var1, "Coverage")) +
    # Add second layer for points below y=x line
    ggnewscale::new_scale_color() +
    # Points below y=x line (var2 coverage) - using viridis palette  
    geom_point(data = data_below, aes(x = af_1, y = af_2, color = coverage_var2, size = coverage_var2)) +
    scale_color_viridis_c(option = "viridis", direction = -1, name = paste(var2, "Coverage")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    xlab(var1) + 
    ylab(var2) +
    labs(subtitle = paste("Points above y=x:", var1, "coverage (magma), below y=x:", var2, "coverage (viridis)"))
  
  # Store plot in list
  plot_list[[i]] <- p
  
  # Also display the plot
  print(p)
  
  # Save individual plot
  filename <- paste0('AF_compare_', gsub(">", "_", var1), '_vs_', gsub(">", "_", var2), '.pdf')
  pdf(filename, width = 7, height = 7)
  print(p)
  dev.off()
}

# Create simplified plots for grid (keep x/y labels and legend, remove main title, square aspect ratio)
plot_list_grid <- list()

## VAF plots for all pairs of variants colored by coverage conditionally

# Define your variants
variants_of_interest <- c('4239_C>A', '3957_A>C', '7411_A>C', '5380_T>A', '11892_A>C', '9982_G>T')

# Get coverage for all variants
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'barcode')

for (v in variants_of_interest) {
  message(v)
  col_name <- v
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric(cutf(v, d = "_")),coverage_info$barcode]
}

# Generate all unique pairs
variant_pairs <- combn(variants_of_interest, 2, simplify = FALSE)

# Create plots for each pair
plot_list <- list()

for (i in seq_along(variant_pairs)) {
  var1 <- variant_pairs[[i]][1]
  var2 <- variant_pairs[[i]][2]
  
  message(paste("Creating plot for", var1, "vs", var2))
  
  # Create plot data
  plot_data <- data.frame(
    af_1 = af.dm[var1,], 
    af_2 = af.dm[var2,], 
    coverage_var1 = coverage_info[[var1]],
    coverage_var2 = coverage_info[[var2]]
  )
  
  # Create grouping variable 
  plot_data$group <- ifelse(plot_data$af_2 > plot_data$af_1, "above", "below")
  
  # Split data for different color palettes
  data_above <- plot_data[plot_data$group == "above", ]
  data_below <- plot_data[plot_data$group == "below", ]
  
  # Create the plot with two different color palettes
  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.7) +
    # Points above y=x line (var1 coverage) - using magma palette
    geom_point(data = data_above, aes(x = af_1, y = af_2, color = coverage_var1, size = coverage_var1)) +
    scale_color_viridis_c(option = "magma", direction = -1, name = paste(var1, "Coverage")) +
    # Add second layer for points below y=x line
    ggnewscale::new_scale_color() +
    # Points below y=x line (var2 coverage) - using viridis palette  
    geom_point(data = data_below, aes(x = af_1, y = af_2, color = coverage_var2, size = coverage_var2)) +
    scale_color_viridis_c(option = "viridis", direction = -1, name = paste(var2, "Coverage")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    xlab(var1) + 
    ylab(var2) +
    coord_fixed(ratio = 1) 
  # Store plot in list
  plot_list[[i]] <- p
  
  # Also display the plot
  print(p)
  
  # Save individual plot as square
  filename <- paste0('AF_compare_', gsub(">", "_", var1), '_vs_', gsub(">", "_", var2), '.pdf')
  pdf(filename, width = 7, height = 7)
  print(p)
  dev.off()
}

# Create gridded version with 2 plots per row, multiple pages
library(gridExtra)
plots_per_page <- 4  # 2 rows x 2 columns = 4 plots per page
n_plots <- length(plot_list)
n_pages <- ceiling(n_plots / plots_per_page)

pdf('AF_compare_all_pairs_grid.pdf', width = 14, height = 14)
for (page in 1:n_pages) {
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(page * plots_per_page, n_plots)
  
  plots_for_page <- plot_list[start_idx:end_idx]
  
  # Use grid.arrange with just the available plots
  do.call(grid.arrange, c(plots_for_page, list(ncol = 2)))
}
dev.off()

