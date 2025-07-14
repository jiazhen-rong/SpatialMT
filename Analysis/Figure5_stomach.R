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

rm(list=ls())
gc()

#MAESTER source code from Miller et al. 
source("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/210215_FunctionsGeneral.R")

#Loading in raw data
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/A210_combined_final/maegatk.rds')
seu <- readRDS(file = '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/Run2/curio_seeker_output/A210_seurat.rds')

## Spatial transcriptomics QC  ------------------------------------------------------------------------------------------------------------------
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/Run2/QC')

SpatialFeaturePlot(seu, 
                   features = "log10_numReads", 
                   pt.size = 1,
                   min.cutoff = "q1",
                   max.cutoff = "q99") +
  scale_fill_gradientn(colors = c("#efefef", "#CBC3E3", "#301934"))+
  labs(title = "Reads per Cell") +
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
  log10_UMIs = seu@meta.data[["log10_numReads"]]
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

ggsave("count_distribution_histogram.pdf", p, width = 10, height = 8)

#how many mitochondrial transcripts do I get? 
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
mean(seu[["percent.mt"]]$percent.mt) #3.375829

mt_genes <- rownames(seu)[grep("^MT-", rownames(seu))]
mt_counts <- Matrix::colSums(seu@assays$RNA@counts[mt_genes, ])
seu$mt_umi_counts <- mt_counts
mean(mt_counts) # 5.99 mitochondrial transcripts per cell
## ------------------------------------------------------------------------------------------------------------------

## RCTD  ---------------------------------------------------------------------------------------------------------------------

#reference generation 
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/reference_dataset')
all_gier <- readRDS('seurat_final.rds') #dataset from Gier et al.

#keep only GC, IM, VC, FB identities
all_gier <- subset(all_gier, ident = 'BE', invert = TRUE)
all_gier <- subset(all_gier, ident = 'NS', invert = TRUE)
all_gier <- subset(all_gier, ident = 'SQ', invert = TRUE)
gastric_cardia <- subset(all_gier, ident = 'TC', invert = TRUE)
rm(all_gier) #memory reasons

gastric_cardia <- FindNeighbors(gastric_cardia, dims = 1:15)
gastric_cardia <- FindClusters(gastric_cardia, resolution = 0.2)
gastric_cardia <- RunUMAP(gastric_cardia, dims = 1:15,spread = 2, min.dist = .1)
DimPlot(gastric_cardia, reduction = "umap", label = TRUE)

stomach_cell_types_genes <- list(
  "Parietal Cells" = c("ATP4A", "ATP4B", "CKB", "CLIC6"),
  "Chief Cells" = c("PGA5", "CHIA", "PGC", "ITIH5", "LIPF"),
  "Enteroendocrine Cells" = c("CHGA", "CHGB", "NEUROD1"),
  "Neck Mucous Cells" = c("TFF2", "PGC", "MUC6","AQP5","FMOD"),
  "Proliferating cells" = c("MKI67", "STMN1"),
  "Pit (Folveolar) Cells" = c("PHGR1","CAPN8", "MUC5AC", "GKN1", "GKN2","TFF1"),
  "Stem/Progenitor Cells" = c("LGR5", "LEFTY1", "OLFM4", "CD44"),
  "Enterocytes" = c("ALPI", "APOA1", "VIL1", "FABP1", "FABP2", "CDH17", "SI", "MUC13", "KRT20")
)

new.cluster.ids <- c("Pit",
                     "Neck", 
                     "Immune", 
                     "Vascular", 
                     "Fibroblast", 
                     "Immune",
                     "Immune",
                     "Parietal",
                     "Enteroendocrine",
                     "Cycling",
                     "Chief",
                     "Immune",
                     "Vascular")

names(new.cluster.ids) <- levels(gastric_cardia)
gastric_cardia <- RenameIdents(gastric_cardia, new.cluster.ids)

pdf('annotations_GC.pdf')
print(DimPlot(gastric_cardia))
dev.off()

saveRDS(gastric_cardia, 'gastric_reference_final.rds')


### Create the Reference object
library(spacexr)
counts <- GetAssayData(object = gastric_cardia, assay = "RNA", slot = "counts")
annotations <- gastric_cardia@active.ident
write.csv(annotations,'GC_celltype_annotations.csv')
counts <- as.data.frame(counts) # converting count matrix into data.frame as read.csv would do
meta_data <- read.csv("GC_celltype_annotations.csv") # load in meta_data (barcodes, clusters)
cell_types <- meta_data$x; names(cell_types) <- colnames(counts) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
reference <- Reference(counts, cell_types,require_int = F)
rm(counts)
saveRDS(reference, 'GC_ref_RCTD.rds')

## running RCTD on  sample:
analysis_folder <-'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/Run2/RCTD'

args=c("A210",'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/Run2/curio_seeker_output','stomach')
slide_name <- args[1] 
print(slide_name)
data_path <- args[2]
print(data_path)
ref_data= args[3] 
print(ref_data)

# load reference data
if(ref_data=="stomach"){
  reference <-readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/reference_dataset/GC_ref_RCTD.rds') 
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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/RCTDresults/A210/RCTD')
pdf(paste0(dir_path,"diagnostic_stats.pdf"),width=4,height=4)
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

dev.off()

myRCTD <- create.RCTD(puck, reference, max_cores = 1, UMI_min = 60)
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

data$cell_type <- ifelse(data$type1_score > 0.6, as.character(data$first_type),
                         ifelse(data$type2_score > 0.6, as.character(data$second_type),
                                "Undetermined"))

annotations <- data %>% select(first_type)
seu <- seu[,colnames(seu) %in% rownames(data)]
seu$annotation <- annotations
Idents(seu) <- seu$annotation

cluster_colors <- c("Chief" = "#66CCEE",
                    "Cycling" = "#CCBB44",
                    "Enteroendocrine" = "#4477AA", 
                    "Neck" = "#228833",
                    "Parietal" = "#AA3377", 
                    "Pit" = "#EE6677",
                    "Fibroblast" = "#8AC3C2", 
                    "Vascular" = "red", 
                    "Immune" = "orange")

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/curio_seeker/A210/Run2/RCTD')
pdf('spatial_RCTD_annotations.pdf')
print(SpatialDimPlot(seu, cols = cluster_colors, pt.size.factor = 1.3))
dev.off()
saveRDS(seu, 'A210_seurat_RCTD.rds')
## ------------------------------------------------------------------------------------------------------------------

## Mitochondrial genotpying  ------------------------------------------------------------------------------------------------------

## Only necessary if barcode in seu does NOT contain '-1' at the end.
seu$cellMerge <- paste0(cutf(colnames(seu), d = "-", f = 1), "-1")

# Only keep cells with a cellMerge id that occurs once
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

# saving final objects:
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis')
saveRDS(seu, 'seurat_final.rds' )
saveRDS(maegatk.rse, 'maegatk_final.rds' )

# Prepare allele frequency matrix
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

# Get the mean quality for each variant
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
write_tsv(vars.tib, "A210_allvariants.txt")
vars.tib <- read_tsv("A210_allvariants.txt")

# Parameters to select clones with high VAF
voi.ch <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 25) %>%
  filter(q99.unionCells> 25) %>% .$var
voi.ch

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

write.table(voi.ch, file = "A210_voi_be25vaf.tsv", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)

#spatial lineage plots
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis/lineage')
for (v in voi.ch) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  cells.tib$af_voi[cells.tib$cov_voi < 3] <- NA
  cells.tib$af_voi[cells.tib$af_voi < 35] <- NA
  
  p <- ggplot(data = cells.tib %>% 
                arrange(!is.na(af_voi)), # This makes NA values plot first (bottom)
              aes(x = SPATIAL_1, y = SPATIAL_2, color = af_voi)) +
    geom_point(aes(size = !is.na(af_voi)), 
               show.legend = FALSE) +
    scale_size_manual(values = c("FALSE" = 0.1, "TRUE" = 0.7)) +
    ggtitle(v) + 
    scale_color_gradientn(colors = c("#bed8e6", "#608da2", "darkblue"), 
                          limits = c(1, 100), 
                          na.value = "#ECECEC") + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) + 
    theme_bw() +
    theme(panel.grid = element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5)) +   
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank())
  name <- paste0(v, ".pdf")
  ggsave(name, plot = p)
  
  cells.tib$af_voi <- NULL
  cells.tib$cov_voi <- NULL
}

#coloring by coverage 
for (v in voi.ch) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  cells.tib$af_voi[cells.tib$cov_voi < 3] <- NA
  cells.tib$af_voi[cells.tib$af_voi < 35] <- NA
  
  p <- ggplot(data = cells.tib,
              aes(x = SPATIAL_1, y = SPATIAL_2, color = cov_voi)) +
    geom_point(size = 1) + 
    ggtitle(v) + 
    scale_color_gradientn(colors = c("#bed8e6", "#608da2", "darkblue"), 
                          limits = c(1, 700), 
                          na.value = "#ECECEC") + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) + 
    theme_bw() +
    theme(panel.grid = element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5)) +   
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank())
  
  name <- paste0(v, "_cov.pdf")
  ggsave(name, plot = p)
  
  cells.tib$af_voi <- NULL
  cells.tib$cov_voi <- NULL
}
## ------------------------------------------------------------------------------------------------------------------

## cell type composition of lineages  ------------------------------------------------------------------------------------------------------

# loading in significant cells from statistical tests: 
sig_cells <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/sig_spots/A210_var_spot_list.rds')
lineage_3866 <- sig_cells[["3866_T>C"]][["fdr"]][["sig_spots"]][["GC"]]
subset_3866 <- seu[,lineage_3866]
#percent of total cells
dim(subset_3866)[2]/dim(seu)[2] * 100

lineage_11969 <- sig_cells[["11969_G>A"]][["fdr"]][["sig_spots"]][["GC"]]
subset_11969  <- seu[,lineage_11969]
#percent of total cells
dim(subset_11969)[2]/dim(seu)[2] * 100

lineage_2245 <- sig_cells[["2245_A>G"]][["fdr"]][["sig_spots"]][["GC"]]
subset_2245  <- seu[,lineage_2245]
#percent of total cells
dim(subset_2245)[2]/dim(seu)[2] * 100


Idents(subset_3866) <- subset_3866$annotation
Idents(subset_11969) <- subset_11969$annotation
Idents(subset_2245) <- subset_2245$annotation

# Get the annotation table
table(subset_3866$annotation)
table(subset_11969$annotation)
table(subset_2245$annotation)
table(seu$annotation)

# Calculate percentages
annotation_counts <- table(subset_2245$annotation)
annotation_percentages <- prop.table(annotation_counts) * 100
plot_data <- data.frame(
  Annotation = names(annotation_percentages),
  Percentage = as.numeric(annotation_percentages)
)

# Create 100% stacked bar plot
ggplot(plot_data, aes(x = factor(1), y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Percentage Contribution of Each Annotation",
       y = "Percentage",
       x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Set3")
ggsave('lineage_2445_detailed_annotations.pdf')

# Calculate percentages
annotation_counts <- table(subset_11969$annotation)
annotation_percentages <- prop.table(annotation_counts) * 100
plot_data <- data.frame(
  Annotation = names(annotation_percentages),
  Percentage = as.numeric(annotation_percentages)
)

# Create 100% stacked bar plot
ggplot(plot_data, aes(x = factor(1), y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Percentage Contribution of Each Annotation",
       y = "Percentage",
       x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Set3")
ggsave('lineage_11969_detailed_annotations.pdf')

# Calculate percentages
annotation_counts <- table(subset_3866$annotation)
annotation_percentages <- prop.table(annotation_counts) * 100
plot_data <- data.frame(
  Annotation = names(annotation_percentages),
  Percentage = as.numeric(annotation_percentages)
)

# Create 100% stacked bar plot
ggplot(plot_data, aes(x = factor(1), y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Percentage Contribution of Each Annotation",
       y = "Percentage",
       x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Set3")
ggsave('lineage_3866_detailed_annotations.pdf')

# Calculate percentages
annotation_counts <- table(seu$annotation)
annotation_percentages <- prop.table(annotation_counts) * 100
plot_data <- data.frame(
  Annotation = names(annotation_percentages),
  Percentage = as.numeric(annotation_percentages)
)

# Create 100% stacked bar plot
ggplot(plot_data, aes(x = factor(1), y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Percentage Contribution of Each Annotation",
       y = "Percentage",
       x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Set3")
ggsave('all_detailed_annotations.pdf')
## ------------------------------------------------------------------------------------------------------------------

#differential gene expression between clusters and cells in a lineage -------------------------------------------------------------

# Broad identities:
tissue_types <- c("Chief" = "GC",
                    "Cycling" = "GC",
                    "Enteroendocrine" = "GC", 
                    "Neck" = "GC",
                    "Parietal" = "GC", 
                    "Pit" = "GC",
                    "Fibroblast" = "NonEpi", 
                    "Vascular" = "NonEpi", 
                    "Immune" = "NonEpi")

seu <- RenameIdents(seu, tissue_types)
seu$tissue <- Idents(seu)

all_gastric_cardia <- colnames(seu[,seu$tissue =='GC'])
cells_in_lineages <- c(colnames(subset_3866), colnames(subset_11969))
nonlineage_gc <- intersect(cells_in_lineages,all_gastric_cardia)

seu_subset <- seu[,seu$tissue =='GC']
seu_subset$lineage_id <- "NotLineage"
lineage_cells <- colnames(seu_subset) %in% cells_in_lineages
seu_subset$lineage_id[lineage_cells] <- "Lineage"

# Set identities to your grouping variable
Idents(seu_subset) <- "lineage_id"

de_results <- FindMarkers(seu_subset, 
                          ident.1 = "Lineage",  
                          ident.2 = "NotLineage", 
                          min.pct = 0.1,       
                          logfc.threshold = 0,  
                          test.use = "wilcox") 


volcano_plot <- EnhancedVolcano(de_results,
                                lab = rownames(de_results),
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                pCutoff = 0.05,
                                FCcutoff = 0.5,
                                pointSize = 2.0,
                                labSize = 3.0,
                                title = 'Group1 vs Group2',
                                subtitle = 'Differential Expression',
                                # Only label top significant genes to avoid crowding
                                selectLab = rownames(de_results)[which(
                                  de_results$p_val_adj < 0.01 & 
                                    abs(de_results$avg_log2FC) > 1)],
                                # Color scheme
                                col = c("gray", "gray", "red", "red"),
                                colAlpha = 1.0,
                                # Add more details
                                caption = "FC cutoff: 0.5, p-value cutoff: 0.05",
                                legendPosition = "right",
                                drawConnectors = FALSE,
                                widthConnectors = 0.5,
                                # Make plot square
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE) + 
  theme(aspect.ratio = 1)  # This makes it square


setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis')
pdf('diff_gene_exp_lineages.pdf')
print(volcano_plot)
dev.off()
## ------------------------------------------------------------------------------------------------------------------

# Mitochondrial genome coverage comparison ------------------------------------------------------------------------------------------------------
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/A210_combined_final/maegatk.rds')
maegatk.curio <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/curio_compare/A210_curio.mr1/final/maegatk.rds')

# Plot coverage along chrM
base.tib1 <- tibble(base = 1:16569, depth = rowSums(assays(maegatk.curio)[["coverage"]]))
base.tib2 <- tibble(base = 1:16569, depth = rowSums(assays(maegatk.rse)[["coverage"]]))

base.tib1.1 <- tibble(base = 1:16569, depth = rowMeans(assays(maegatk.curio)[["coverage"]]))
base.tib2.1 <- tibble(base = 1:16569, depth = rowMeans(assays(maegatk.rse)[["coverage"]]))
mean(base.tib1.1$depth) # .33 transcripts per bead
mean(base.tib2.1$depth) # 5.49 transcripts per bead

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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis')
pdf('A210_boxplot.pdf')
print(ratio_plot)
dev.off()
## ------------------------------------------------------------------------------------------------------------------

## Do variants co-occur in cells? 
## VAF for 3866_T>C vs 11969_G>A colored by coverage: 

variants_of_interest <- c('3866_T>C', '11969_G>A','2245_A>G')
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'barcode')

for (v in variants_of_interest) {
  message(v)
  col_name <- v
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),coverage_info$barcode]
  coverage_info[col_name] 
}

plot_data <- data.frame(
  af_1 = af.dm['3866_T>C',], 
  af2 = af.dm['11969_G>A',],
  coverage_3866 = coverage_info$`3866_T>C`,
  coverage_11969 = coverage_info$`11969_G>A`
)

# Create grouping variable 
plot_data$group <- ifelse(plot_data$af2 > plot_data$af_1, "above", "below")

# Split data for different color palettes
data_above <- plot_data[plot_data$group == "above", ]
data_below <- plot_data[plot_data$group == "below", ]

p1 <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 1) +
  # Points above y=x line (3866 coverage) - using magma palette
  geom_point(data = data_above, aes(x = af_1, y = af2, color = coverage_3866,
                                    size = coverage_3866, alpha =1)) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "3866_T>C Coverage") +
  # Add second layer for points below y=x line
  ggnewscale::new_scale_color() +
  # Points below y=x line (11969 coverage) - using viridis palette  
  geom_point(data = data_below, aes(x = af_1, y = af2, color = coverage_11969,
                                    size = coverage_11969, alpha = 1)) +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "11969_G>A Coverage") +
  coord_fixed() +  # Makes the plot square
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

pdf('AF_compare_spatiallineages_3866_11969.pdf')
print(p1)
dev.off()

cells_3866 <- sum(plot_data$af_1 > 0.5 & plot_data$coverage_3866 > 5, =na.rm = TRUE)
print(paste("Cells with 3866_T>C AF > 50% and coverage > 10:", cells_3866))

# For variant 11969_G>A  
cells_11969 <- sum(plot_data$af2 > 0.5 & plot_data$coverage_11969 > 5, na.rm = TRUE)
print(paste("Cells with 11969_G>A AF > 50% and coverage > 10:", cells_11969))


# Plot 2: 3866 vs 2245
plot_data_2 <- data.frame(
  af_1 = af.dm['3866_T>C',], 
  af2 = af.dm['2245_A>G',],
  coverage_3866 = coverage_info$`3866_T>C`,
  coverage_2245 = coverage_info$`2245_A>G`
)

# Create grouping variable 
plot_data_2$group <- ifelse(plot_data_2$af2 > plot_data_2$af_1, "above", "below")

data_above_2 <- plot_data_2[plot_data_2$group == "above", ]
data_below_2 <- plot_data_2[plot_data_2$group == "below", ]

p2 <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 1) +
  # Points above y=x line (3866 coverage) - using magma palette
  geom_point(data = data_above_2, aes(x = af_1, y = af2, color = coverage_3866,
                                      size = coverage_3866, alpha = 1.0)) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "3866_T>C Coverage") +
  # Add second layer for points below y=x line
  ggnewscale::new_scale_color() +
  # Points below y=x line (2245 coverage) - using viridis palette  
  geom_point(data = data_below_2, aes(x = af_1, y = af2, color = coverage_2245,
                                      size = coverage_2245, alpha = 1.0)) +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "2245_A>G Coverage") +
  coord_fixed() +  # Makes the plot square
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  xlab('3866_T>C') + 
  ylab('2245_A>G')

p2

pdf('AF_compare_spatiallineages_3866_2245.pdf')
print(p2)
dev.off()

# Plot 3: 11969 vs 2245
plot_data_3 <- data.frame(
  af_1 = af.dm['11969_G>A',], 
  af2 = af.dm['2245_A>G',],
  coverage_11969 = coverage_info$`11969_G>A`,
  coverage_2245 = coverage_info$`2245_A>G`
)

# Create grouping variable 
plot_data_3$group <- ifelse(plot_data_3$af2 > plot_data_3$af_1, "above", "below")

data_above_3 <- plot_data_3[plot_data_3$group == "above", ]
data_below_3 <- plot_data_3[plot_data_3$group == "below", ]

p3 <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 1) +
  # Points above y=x line (11969 coverage) - using magma palette
  geom_point(data = data_above_3, aes(x = af_1, y = af2, color = coverage_11969,
                                      size = coverage_11969, alpha = 1.0)) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "11969_G>A Coverage") +
  # Add second layer for points below y=x line
  ggnewscale::new_scale_color() +
  # Points below y=x line (2245 coverage) - using viridis palette  
  geom_point(data = data_below_3, aes(x = af_1, y = af2, color = coverage_2245,
                                      size = coverage_2245, alpha = 1.0)) +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "2245_A>G Coverage") +
  coord_fixed() +  # Makes the plot square
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  xlab('11969_G>A') + 
  ylab('2245_A>G')

pdf('AF_compare_spatiallineages_11969_2245.pdf')
print(p3)
dev.off()


## relationships between lineages: 
plot_data_comprehensive <- data.frame(
  af_3866 = af.dm['3866_T>C',],   
  af_2245 = af.dm['2245_A>G',],
  af_11969 = af.dm['11969_G>A',],
  coverage_3866 = coverage_info$`3866_T>C`,
  coverage_2245 = coverage_info$`2245_A>G`,
  coverage_11969 = coverage_info$`11969_G>A`
)

## For 2245 and 3866: 
# Analysis 1: Spots with 2245_A>G (AF > 25%, cov > 5) but 0 VAF for 3866_T>C
spots_2245_positive_3866_zero <- sum(
  plot_data_comprehensive$af_2245 > 25 & 
    plot_data_comprehensive$af_3866 < 25 &
   plot_data_comprehensive$coverage_3866 > 3,
  na.rm = TRUE
)

print(paste("Spots with 2245_A>G (AF > 25%, cov > 5) but 0 VAF for 3866_T>C:", spots_2245_positive_3866_zero))

# Analysis 2: Spots with both variants (high confidence)
spots_both_variants <- sum(
  plot_data_comprehensive$af_2245 > 25 & 
    plot_data_comprehensive$coverage_2245 > 3 & 
    plot_data_comprehensive$af_3866 > 25 & 
    plot_data_comprehensive$coverage_3866 > 3, 
  na.rm = TRUE
)

print(paste("Spots with both 2245_A>G and 3866_T>C (AF > 25%, cov > 5):", spots_both_variants))

# Analysis 3: Spots with 3866_T>C but not 2245_A>G
spots_3866_only <- sum(
  plot_data_comprehensive$af_3866 > 25 & 
    plot_data_comprehensive$af_2245 < 25 & 
    plot_data_comprehensive$coverage_2245 > 5, 
  na.rm = TRUE
)

print(paste("Spots with 3866_T>C (AF > 25%, cov > 5) but not 2245_A>G (AF < 25%, cov > 5):", spots_3866_only))

## for 11969 and 2245
# Analysis 1: Spots with 2245_A>G (AF > 25%, cov > 5) but 0 VAF for 11969_G>A
spots_2245_positive_11969_zero <- sum(
  plot_data_comprehensive$af_2245 > 25 & 
    plot_data_comprehensive$coverage_11969 > 5 & 
    plot_data_comprehensive$af_11969 < 25, 
  na.rm = TRUE
)

print(paste("Spots with 2245_A>G (AF > 25%, cov > 5) but 0 VAF for 11969_G>A:", spots_2245_positive_11969_zero))

# Analysis 2: Spots with both 11969_G>A and 2245_A>G (high confidence)
spots_both_11969_2245 <- sum(
  plot_data_comprehensive$af_2245 > 25 & 
    plot_data_comprehensive$coverage_2245 > 5 & 
    plot_data_comprehensive$af_11969 > 25 & 
    plot_data_comprehensive$coverage_11969 > 5, 
  na.rm = TRUE
)

print(paste("Spots with both 2245_A>G and 11969_G>A (AF > 25%, cov > 5):", spots_both_11969_2245))

# Analysis 3: Spots with 11969_G>A but not 2245_A>G
spots_11969_only <- sum(
  plot_data_comprehensive$af_11969 > 25 & 
    plot_data_comprehensive$af_2245 < 25 & 
    plot_data_comprehensive$coverage_2245 > 5, 
  na.rm = TRUE
)

print(paste("Spots with 11969_G>A (AF > 25%, cov > 5) but not 2245_A>G (AF < 25%, cov > 5):", spots_11969_only))

## For 11696 and 3866: 
# Analysis 1: Spots with 3866_T>C (AF > 25%, cov > 5) but 0 VAF for 11969_G>A
spots_3866_positive_11969_zero <- sum(
  plot_data_comprehensive$af_3866 > 25 & 
    plot_data_comprehensive$af_11969 < 25 & 
    plot_data_comprehensive$coverage_11969 > 5,
  na.rm = TRUE
)

print(paste("Spots with 3866_T>C (AF > 25%, cov > 5) but 0 VAF for 11969_G>A:", spots_3866_positive_11969_zero))

# Analysis 2: Spots with both 11969_G>A and 3866_T>C (high confidence)
spots_both_11969_3866 <- sum(
  plot_data_comprehensive$af_3866 > 25 & 
    plot_data_comprehensive$coverage_3866 > 5 & 
    plot_data_comprehensive$af_11969 > 25 & 
    plot_data_comprehensive$coverage_11969 > 5, 
  na.rm = TRUE
)

print(paste("Spots with both 3866_T>C and 11969_G>A (AF > 25%, cov > 5):", spots_both_11969_3866))

# Analysis 3: Spots with 11969_G>A but not 3866_T>C
spots_11969_only_vs_3866 <- sum(
  plot_data_comprehensive$af_11969 > 25 & 
    plot_data_comprehensive$af_3866 < 25 & 
    plot_data_comprehensive$coverage_3866 > 5, 
  na.rm = TRUE
)

print(paste("Spots with 11969_G>A (AF > 25%, cov > 5) but not 3866_T>C (AF < 25%, cov > 5):", spots_11969_only_vs_3866))



## ------------------------------------------------------------------------------------------------------------------

## Plots from Caleb Lareau's critcism of ReDeem: 

# plot 1: # of cells with the mutation vs. # of reads per molecule per cell thats part of the lineage, coloring the variants of interest. 

# Parameters to select all clones detected by maegatk
voi.ch <- vars.tib %>% filter(mean_cov.unionCells > 1,quality > 15) %>%
  filter(q99.unionCells> 15) %>% .$var

af.dm_mini <- af.dm[voi.ch,]

#First, from af.dm calculating how many cells there are with a given mutation (threshold is AF > 20%)
cells_per_lineage <- data.frame(
  Variant = rownames(af.dm_mini),
  NumCells = rowSums(af.dm_mini > 15, na.rm = TRUE)
)

#subsetting to only plot variants that have any cells.
has_cells <- cells_per_lineage %>% filter(cells_per_lineage$NumCells > 1)
af.dm_mini <- af.dm_mini[rownames(has_cells),]
af.dm_mini[af.dm_mini < 10 ] <- NA #unreliable low VAF

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
spatial_variants <- c('3688_T>C', '11969_G>A', '2245_A>G')

plot_data <- data.frame(
  Variant = colnames(coverage_info),
  meanCoverage = colMeans(coverage_info, na.rm = TRUE), 
  cellsperlienage = rowSums(af.dm_mini > 15, na.rm = TRUE)
)

variants <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 21) %>%
  filter(q99.unionCells > 23) %>% .$var


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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis/artifact')
pdf('lineagesize_vs_coverage.pdf')
print(p1)
dev.off()

## Making a graph of mismatch vs. position in the read (is there end of sequence bias?)
counts_per_position <- read.csv('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis/artifact/position_summary_simple.csv')
counts_per_position$Position <- (counts_per_position$Position) / 250 #getting fraction for plotting

positions <- data.frame(Variant = spatial_variants, position = c('89', '69', '135'))

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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB040/MAESTER/MAESTER_a210/combined/analysis/artifact')
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
variants_raw <- c("3866_T>C", "11969_G>A", '2245_A>G')

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


