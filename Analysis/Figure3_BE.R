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
library(circlize)
library(ggrastr)
library(scCustomize)
library(stringr)
library(Cairo)

rm(list=ls())
gc()

#MAESTER source code from Miller et al. 
source("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/210215_FunctionsGeneral.R")

# RCTD  ---------------------------------------------------------------------------------------------------------------------
#reference generation 
#reference dataset from Rodrigo to make RCTD reference out of: 
be_sc <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Gier_reviews/2023_revision/seurat_final.rds')
DimPlot(be_sc)

#removing some of the identities to make a cleaner object
be_sc <- subset(be_sc, ident = 'GC')
saveRDS(be_sc, '/Users/sydneybracht/Downloads/GC_only_singlcell.rds')

cluster_colors <- c("BE" = "#E00081",
                    "FB" = "#8AC3C2", 
                    "VC" = "#A2B7D7", 
                    "SQ" = "#EA685F", 
                    "IM" = "#538E8E")

plt <- DimPlot(be_sc, cols = cluster_colors, raster  = FALSE)
pdf('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/RCTD/UMAP.pdf')
print(plt)
dev.off()

# Downsample to 10,000 cells per identity. 
number_cells <- 10000
cell_types <- unique(be_sc$type)
selected_cells <- c()
set.seed(42)
cells_per_type <- 10000

# Sample cells from each cell type
for (ct in cell_types) {
  # Get cells of this type
  cells_of_type <- rownames(be_sc@meta.data[be_sc$type == ct, ])
  
  # Sample cells (or take all if fewer than requested)
  n_cells <- min(cells_per_type, length(cells_of_type))
  sampled_cells <- sample(cells_of_type, n_cells, replace = FALSE)
  
  # Add to the list of selected cells
  selected_cells <- c(selected_cells, sampled_cells)
}

# Subset the Seurat object with selected cells
downsampled_seurat <- subset(be_sc, cells = selected_cells)

# Verify the number of cells per cell type
print(table(downsampled_seurat$type))
rm(be_sc)

### Create the Reference object
library(spacexr)
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/RCTD')
counts <- GetAssayData(object = be_sc, assay = "RNA", slot = "counts")
annotations <- be_sc@active.ident
write.csv(annotations,'BE_celltype_annotations.csv')

rm(downsampled_seurat)
counts <- as.data.frame(counts) # converting count matrix into data.frame as read.csv would do

meta_data <- read.csv("BE_celltype_annotations.csv") # load in meta_data (barcodes, clusters)
cell_types <- meta_data$x; names(cell_types) <- colnames(counts) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
reference <- Reference(counts, cell_types,require_int = F)
rm(counts)
saveRDS(reference, 'BE_ref_celltypes_RCTD.rds')
rm(reference)

## running RCTD on  sample:
seu <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/curio_seeker_output/a200_s6_seurat.rds')

analysis_folder = "/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/RCTD"
args=c("a200_s6","/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/curio_seeker_output", 'subset')
slide_name <- args[1] 
print(slide_name)
data_path <- args[2]
print(data_path)
ref_data= args[3] # choose from "Major" or "subset"
print(ref_data)

# load reference data
if(ref_data=="subset"){
  reference <-readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/RCTD/major/BE_ref_reduced_RCTD.rds') 
}else if(ref_data=="detailed"){
  reference <-readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB034/curio_seeker/RCTD Deconvolution/BE_reference/be_ref_RCTD.rds') 
}

# output directory
dir_path <- '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/RCTD'

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

### Create SpatialRNA object
puck <- SpatialRNA(coords, mtx, nUMI)

## Examine SpatialRNA object (optional)
print(dim(puck@counts)) # observe Digital Gene Expression matri
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

pdf(paste0(dir_path,"diagnostic_stats.pdf"),width=4,height=4)
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

dev.off()

myRCTD <- create.RCTD(puck, reference, max_cores = 1,UMI_min = 50)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- dir_path
dir.create(resultsdir)
#> Warning in dir.create(resultsdir): 'RCTD_Plots' already exists
# Plotting

# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

# save the weights
saveRDS(myRCTD,paste0(dir_path,"/RCTD_res.rds"))

data <- data.frame(results[["results_df"]],results[["weights_doublet"]]) 
data <- data %>% rename(first_type.1 = 'type1_score', second_type.1 = 'type2_score')
data <- data %>% select(-first_class, -second_class, -conv_all, -conv_doublet, -min_score, -singlet_score)
data <- data %>% filter(spot_class != 'reject') %>% select(-spot_class)

data$cell_type <- ifelse(data$type1_score > 0.70, as.character(data$first_type),
                         ifelse(data$type2_score > 0.70, as.character(data$second_type),
                                "Undetermined"))

data <- data %>% filter(cell_type != 'Undetermined')
annotations <- data %>% select(cell_type)

seu <-seu[,colnames(seu) %in% rownames(data)]
seu$annotation <- annotations
Idents(seu) <- seu$annotation

cluster_colors_subset <- c("BE" = "#E00081",
                    "FB" = "#8AC3C2", 
                    "VC" = "#A2B7D7", 
                    "SQ" = "#EA685F", 
                    "IM" = "#538E8E")

cluster_colors_detailed <- c("Ciliated" = "#9FD9C9",
                    "Colonocyte" = "#F9F7BE", 
                    "Crypt" = "#C5C2E2", 
                    "Cycling" = "#E78D7B", 
                    "Enterocyte" = "#8CBADC",
                    "Enteroendocrine" = "#F0C389",
                    "Foveolar" = "#C1E09F",
                    "Goblet" = "#F3C5D3",
                    "Intermediate" = "#6B5EA1",
                    "Neck" = "#BB85B4")

plt <- SpatialDimPlot(seu, cols = cluster_colors)
pdf('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/RCTD/RCTD_spatial_plot.pdf')
print(plt)
dev.off()

## Mitochondrial genotpying  ------------------------------------------------------------------------------------------------------

maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/comined_runs/maegatk_output/maegatk.rds')
## Only necessary if barcode in seu does NOT contain '-1' at the end.
seu$cellMerge <- paste0(cutf(colnames(seu), d = "-", f = 1), "-1")

# Only keep cells with a cellMerge id that occurs once, then intersect with Maester data
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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/')
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
write_tsv(vars.tib, "A200s6_variants.txt")
vars.tib <- read.table('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/a200_s6_allvariants.txt', header = TRUE)

# Parameters to select clones with high VAF
voi.ch <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 25) %>%
  filter(q99.unionCells> 25) %>% .$var
voi.ch

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

write.table(voi.ch, file = "A200s6_voi_be25vaf.tsv", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)

#spatial lineage plots
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/lineageplots_newcolor')
voi.ch <- c('15777_G>C', '3071_T>C', '3054_G>C')
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
    scale_color_gradientn(colors = c("white", "#608da2", "darkblue"), 
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

#Combo plot 
variant_data <- list()

# Process each variant
for (v in voi.ch) {
  temp_data <- cells.tib
  
  # Add info of variant of interest
  temp_data$af_voi <- af.dm[v,temp_data$cell]
  temp_data$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric(cutf(v, d = "_")), temp_data$cell]
  temp_data$af_voi[temp_data$cov_voi < 3] <- NA
  temp_data$af_voi[temp_data$af_voi < 35] <- NA
  
  # Add variant name
  temp_data$variant <- v
  
  # Store in list
  variant_data[[v]] <- temp_data
}

# Combine all variant data
combined_data <- bind_rows(variant_data)

# Create the plot
p <- ggplot(data = combined_data %>% 
              arrange(!is.na(af_voi)), # This makes NA values plot first (bottom)
            aes(x = SPATIAL_1, y = SPATIAL_2)) +
  # Layer for NA values (grey dots)
  geom_point(data = combined_data %>% filter(is.na(af_voi)),
             size = 0.1, color = "#ECECEC") +
  # Layer for non-NA values (colored by variant)
  geom_point(data = combined_data %>% filter(!is.na(af_voi)),
             aes(color = variant), size = 0.7) +
  scale_color_manual(values = c("3071_T>C" = "#6B251E", 
                                "3054_G>C" = "#1D702E", 
                                "15777_G>C" = "#1E1F6A")) +
  ggtitle("Multiple Variants") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

# Save the plot
ggsave("combined_variants.pdf", plot = p)

# bulk AF for variants of interest: 
mean(af.dm['15777_G>C',])
mean(af.dm['3071_T>C',])
mean(af.dm['3054_G>C',])
## ------------------------------------------------------------------------------------------------------------------

# Coverage comparison before and after enrichment
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/curio_compare')
maegatk.rse  <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/comined_runs/maegatk_output/maegatk.rds')
maegatk.curio <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/curio_compare/output/final/maegatk.rds')

seu <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/curio/a200_s6/curio_seeker_output/a200_s6_seurat.rds')

colnames(seu) <- paste0(colnames(seu), '-1')

#subseting maegatk objects to only contain cells found in the Seeker data
int <- intersect(colnames(maegatk.rse), colnames(seu))
maegatk.rse <- maegatk.rse[, int]

int2 <- intersect(colnames(maegatk.curio), colnames(seu))
maegatk.curio <- maegatk.curio[, int2]

#subsetting so the two have the same cells
int3 <- intersect(colnames(maegatk.curio), colnames(maegatk.rse))
maegatk.curio <- maegatk.curio[, int3]
maegatk.rse <- maegatk.rse[, int3]
seu <- seu[, colnames(maegatk.rse)]

#make sure theyre the same
dim(maegatk.curio)
dim(maegatk.rse)

# Plot mean coverage per position per cell ----------------------------------------------------------------------
# Set y axis parameters
ymax <- 100000

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

a200s6_coverage_curio <- 
  base.tib1 %>% ggplot() +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-ymax*0.2, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
  ylab("Reads per position") + xlab("Position along chrM") + ggtitle('BE Biopsy before enrichment') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))

a200s6_coverage_curio
pdf('A200s6_coverage_curio.pdf')
print(a200s6_coverage_curio)
dev.off()

a200s6_coverage_maester <- 
  base.tib2 %>% ggplot() +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  ylab("Mean coverage per cell") + xlab("Position along chrM") + ggtitle('HGD Biopsy after enrichment') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))

a200s6_coverage_maester
pdf('A200s6_coverage_maester.pdf')
print(a200s6_coverage_maester)
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

pdf('A200s6_boxplot.pdf')
print(ratio_plot)
dev.off()

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
  ggtitle('A200s6 comparison: Before and After enrichment') + 
  theme_classic() +
  theme(aspect.ratio = 0.5,
        plot.title = element_text(hjust = 0.5))

combined_coverage

# Save to PDF
pdf('A200s6_coverage_combined.pdf')
print(combined_coverage)
dev.off()
## ------------------------------------------------------------------------------------------------------------------

## cell type composition of lineages  ------------------------------------------------------------------------------------------------------

# loading in significant cells from statistical tests: 
sig_cells <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_var_spot_list.rds')
lineage_15777 <- sig_cells[["15777_G>C"]][["permute"]][["sig_spots"]][["BE"]]
subset_15777 <- seu[,lineage_15777]

lineage_3071 <- sig_cells[["3071_T>C"]][["permute"]][["sig_spots"]][["BE"]]
subset_3071  <- seu[,lineage_3071]

lineage_3054 <- sig_cells[["3054_G>C"]][["permute"]][["sig_spots"]][["BE"]]
subset_3054  <- seu[,lineage_3054]

cells_in_lineages <- c(colnames(subset_3054), colnames(subset_3071), colnames(subset_15777))

subset_3054 <-subset_3054[,colnames(subset_3054) %in% rownames(data)]
subset_3054$annotation <- annotations
Idents(subset_3054) <- subset_3054$annotation

subset_3071 <-subset_3071[,colnames(subset_3071) %in% rownames(data)]
subset_3071$annotation <- annotations
Idents(subset_3071) <- subset_3071$annotation

subset_15777 <-subset_15777[,colnames(subset_15777) %in% rownames(data)]
subset_15777$annotation <- annotations
Idents(subset_15777) <- subset_15777$annotation

# Get the annotation table
table(subset_15777$annotation)
table(subset_3071$annotation)
table(subset_3054$annotation)
table(seu$annotation)

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

## Plots from Caleb Lareau's critcism of ReDeem: 

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
spatial_variants <- c('15777_G>C', '3071_T>C', '3054_G>C')

plot_data <- data.frame(
  Variant = colnames(coverage_info),
  meanCoverage = colMeans(coverage_info, na.rm = TRUE), 
  cellsperlienage = rowSums(af.dm_mini > 15, na.rm = TRUE)
)

variants <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 21) %>%
  filter(q99.unionCells > 23) %>% .$var

plot_data$annotation <- "Detected"
plot_data$annotation[plot_data$Variant %in% variants] <- "Significant by maegatk"
plot_data$annotation[plot_data$Variant %in% spatial_variants] <- "Spatial Variant"


# Convert annotation to factor with specific levels for proper ordering in legend
plot_data$annotation <- factor(plot_data$annotation, levels = c("Detected", "Significant by maegatk", "Spatial Variant"))

# Then use that column for coloring
p1 <- ggplot(plot_data, aes(x = cellsperlienage, y = meanCoverage, color = annotation)) +
  geom_point(size = 1.25) + 
  geom_hline(yintercept = 3, linetype = 2) +
  scale_x_log10() + 
  scale_y_log10() + 
  xlab('Cells with variant') + 
  ylab('Mean reads per cell with variant') +
  scale_color_manual(values = c("Detected" = "lightgrey", "Significant by maegatk" = "blue", "Spatial Variant" = "red")) +
  theme_classic()
p1

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/not_artifact_plot')
pdf('lineagesize_vs_coverage.pdf')
print(p1)
dev.off()

## Making a graph of mismatch vs. position in the read (is there end of sequence bias?)

counts_per_position <- read.csv('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/not_artifact_plot/position_counts_simple.csv')
counts_per_position$Position <- (counts_per_position$Position) / 250 #getting fraction for plotting

positions <- data.frame(Variant = spatial_variants, position = c('134.0', '120.0', "69.0"))

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
variants_raw <- voi.ch

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
