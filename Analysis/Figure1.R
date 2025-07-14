library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(SummarizedExperiment)
library(Matrix)
library(ggrastr)
library(RColorBrewer)
library(scDblFinder)
library(SoupX)
library(DropletUtils)

rm(list=ls())
gc()

setwd("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat")

### Mouse tumor single cell RNA-seq data ------------------------------------------------------------------------------------------------------------------------
seu.mt <- Read10X(
  '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/raw_data/outs/MT_adjusted_counts_rho=0.1',
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE)
seu.mt <- CreateSeuratObject(seu.mt)

## subseting Seurat objects based on sufficient RNA content.
VlnPlot(seu.mt, feature = c('nCount_RNA', 'nFeature_RNA'))
seu.mt <- subset(seu.mt, subset = nFeature_RNA > 250 & nFeature_RNA > 300)

## subseting Seurat objects based on low MT content. 
seu.mt[["percent.mt.human"]] <- PercentageFeatureSet(seu.mt, pattern = "^GRCh38-MT-")
seu.mt[["percent.mt.mouse"]] <- PercentageFeatureSet(seu.mt, pattern = "^GRCm39-mt-")
VlnPlot(seu.mt, feature = c('percent.mt.human', 'percent.mt.mouse'))
seu.mt <- subset(seu.mt, subset = percent.mt.human < 10)
seu.mt <- subset(seu.mt, subset = percent.mt.mouse < 10)

# preprocessing 
seu.mt <- NormalizeData(seu.mt, normalization.method = "LogNormalize", scale.factor = 10000) 
seu.mt <- FindVariableFeatures(seu.mt, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu.mt)
seu.mt <- ScaleData(seu.mt, features = all.genes)
seu.mt <- RunPCA(seu.mt, features = VariableFeatures(object = seu.mt), verbose = FALSE)

#doublet removal
sce <- as.SingleCellExperiment(seu.mt)
sce <- scDblFinder(sce, nfeatures = 3000,includePCs = 1:20)
doublet_status <- sce$scDblFinder.class
rm(sce)

seu.mt <- AddMetaData(object = seu.mt, metadata = doublet_status, col.name = "doublets")
seu.mt <- subset(seu.mt, subset = doublets == "singlet")
saveRDS(seu.mt, '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/rds/MT_filtered.rds')

ElbowPlot(seu.mt)
set.seed(101)
seu.mt <- FindNeighbors(seu.mt, dims = 1:10)
seu.mt <- FindClusters(seu.mt, resolution = 0.1)
seu.mt <- RunUMAP(seu.mt, dims = 1:10, spread = 1, min.dist = 0.3)
plt <- DimPlot(seu.mt, reduction = 'umap', label = TRUE)

pdf('UMAP.pdf')
print(plt)
dev.off()

markers <- FindAllMarkers(seu.mt, min.pct = .25)
write.csv(markers,'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat/markers.csv')

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat/feature_plots')
#QC UMAP plots
counts <- FeaturePlot(seu.mt, features = 'nCount_RNA')
pdf('UMAP_counts.pdf')
print(counts)
dev.off()

genes <- FeaturePlot(seu.mt, features = 'nFeature_RNA')
pdf('UMAP_features.pdf')
print(genes)
dev.off()

## genes from literature that should be differentially expressed in the two cell lines:
# Gene markers for SW1990: CALB1, DHRS9, GNGT1
# Gene markers for CAPAN2: NALF1, PRKG1, PLCH1

SW1990 <- FeaturePlot(seu.mt, features = c('GRCh38-CALB1', 'GRCh38-DHRS9','GRCh38-GNGT1'))
pdf('UMAP_SWI_markers.pdf')
print(SW1990)
dev.off()
#clusters 0 and 3

capan2 <- FeaturePlot(seu.mt, features = c('GRCh38-NALF1', 'GRCh38-PRKG1','GRCh38-PLCH1'))
pdf('UMAP_CAPAN2_markers.pdf')
print(capan2)
dev.off()
#cluster 1 and 2 

SM <- FeaturePlot(seu.mt, features = c('GRCm39-Actb')) #smooth muscle 
pdf('UMAP_smoothmuscle_markers.pdf')
print(SM)
dev.off()
#cluster 4/5

IM <- FeaturePlot(seu.mt, features = c('GRCm39-Il1b')) #Immune
pdf('UMAP_immune_markers.pdf')
print(IM)
dev.off()
#cluster 4/5

## species annotation
species <- c('Human', 'Human', 'Human','Human','Mouse', 'Mouse',"Mouse")
names(species) <- levels(seu.mt)
seu.mt <- RenameIdents(object = seu.mt, species)

plt <- DimPlot(seu.mt)
pdf('UMAP_speciesannotation.pdf')
print(plt)
dev.off()

## cell type annotations
seu.mt <- FindClusters(seu.mt, resolution = 0.1)
cluster_ids <- c('SW1990', 'CAPAN2', 'CAPAN2', 'SW1990','IM-TME', 'IM-TME',"IM-TME")
names(cluster_ids) <- levels(seu.mt)
seu.mt <- RenameIdents(object = seu.mt, cluster_ids)

# confirming DEGs with that of literature 
FindMarkers(seu.mt, ident.1 = 'CAPAN2', ident.2 = 'SW1990')
#DEGs from confidently annotated clusters.
# Gene markers for SW1990: CALB1, CEACAM6, GNGT1, MUC5AC, VIM
# Gene markers for CAPAN2: NALF1, PRKG1, PLCH1, CDKN2A

cluster_colors <- c("SW1990" = "#D55E00", 
                    "CAPAN2" = "#009E73", 
                    "IM-TME" = "#E1BE6A")

plt <- DimPlot(seu.mt, cols = cluster_colors)
pdf('UMAP_mousetumor_annotated.pdf')
print(plt)
dev.off()

saveRDS(seu.mt,'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat/rds/MT_final.rds')

## Make reference genome with this data for RCTD: 
library(spacexr)

#Do I need to downsample (RCTD takes maximum 10,000 cells per cluster)
table(seu.mt@active.ident)

setwd("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat/rctd_ref")
counts <- GetAssayData(object = seu.mt, assay = "RNA", slot = "counts")
counts <- as.data.frame(counts) # converting count matrix into data.frame as read.csv would do

cluster_ids <- c('SW1990', 'CAPAN2','IM-TME')
names(cluster_ids) <- levels(seu.mt)
seu.mt <- RenameIdents(object = seu.mt, cluster_ids)

annotations <- seu.mt@active.ident
write.csv(annotations,'MT_celltype_annotations.csv')

meta_data <- read.csv("MT_celltype_annotations.csv") # load in meta_data (barcodes, clusters)
cell_types <- meta_data$x; names(cell_types) <- colnames(counts) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type

### Create the Reference object
mouse_reference <- Reference(counts, cell_types,require_int = F)
saveRDS(mouse_reference, 'PDAC_mouse_ref_RCTD.rds')

## Single cell mitochondrial genotyping -------------------------------------------------------------

#preprocessing 
#MAESTER source code from Miller et al. 
source("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/210215_FunctionsGeneral.R")

#loading in data
seu <- readRDS(file = '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat/rds/MT_final.rds')
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/MAESTER/MT_sc_final/maegatk.rds')

# Generate cell IDs that will match Maegatk
seu$cellMerge <- colnames(seu)

# Only keep cells with a cellMerge id that occurs once, then intersect with Maester data
cells.ch <- tibble(cell = seu$cellMerge) %>% group_by(cell) %>% filter(n()==1) %>% .$cell %>% unname
cells.ch <- cells.ch[cells.ch %in% colnames(maegatk.rse)]

# Subset for common cells
seu <- subset(seu, subset = cellMerge %in% cells.ch)
seu <- RenameCells(seu, new.names = as.character(seu$cellMerge))
seu$cellMerge <- NULL

seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]

saveRDS(seu, file = "MT_sc_final.rds")

# Subset and order RSE object
maegatk.rse <- maegatk.rse[,colnames(seu)]
saveRDS(maegatk.rse, file = "MT_maegatk_final.rds")

# Prepare allele frequency matrix
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))*100
write.csv(af.dm, file = "af_dm.csv")

# Check (should be TRUE)
all(colnames(af.dm) == colnames(seu))

# Extract Seurat metadata
cells.tib <- as_tibble(seu@meta.data, rownames = "cell")

CellSubsets.ls <- list(unionCells = cells.tib$cell)

lengths(CellSubsets.ls)
# -----------------------------------------------------------------------------------------------------------------------------------------
# Get the mean allele frequency and coverage for every cell subset
mean_af.ls <- lapply(CellSubsets.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(CellSubsets.ls, function(x) rowMeans(assays(maegatk.rse)[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset
quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)

library(stringr)
start_time <- Sys.time()
quantiles.ls <- lapply(quantiles, function(x) lapply(CellSubsets.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
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
write_tsv(vars.tib, "MT_sc_allvariants.txt")
vars.tib <- read_tsv('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/MAESTER/MT_sc_final/analysis/MT_sc_allvariants.txt')

# Parameters to select clones with high VAF
voi.ch <- vars.tib %>% filter(mean_cov.unionCells >= 1, quality >= 27) %>% 
  filter(q99.unionCells > 25) %>% .$var

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )

write.table(voi.ch, file = "MT_sc_25vaf.tsv", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)
# -----------------------------------------------------------------------------------------------------------------------------------------
## lineage plots
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/MAESTER/MT_sc_final/analysis/lineage')
for (v in voi.ch) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  cells.tib$af_voi[cells.tib$cov_voi < 5] <- NA
  
  p <- ggplot(data = cells.tib, aes(x = UMAP_1, y = UMAP_2, color = af_voi))+
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
# -----------------------------------------------------------------------------------------------------------------------------------------
#QC metrics -- Violin plot of AF for two variants of interest: 

# Create a dataframe with just the column of interest
seu_2 <- subset(seu, ident = 'IM-TME', invert = TRUE)
af_voi_3010 <- af.dm['3010_G>A',colnames(seu_2)]
af_voi_9545 <- af.dm['9545_A>G',colnames(seu_2)]

plot_data <- data.frame(
  cell = colnames(seu_2), 
  type = Idents(seu_2),
  values_3010= af_voi_3010, 
  values_9545 = af_voi_9545
)

# Set custom colors for each cell type
violin_colors <- c("CAPAN2" = "#133C5C", "SW1990" = "#804D60")
dot_colors <- c("CAPAN2" = "#1E206B", "SW1990" = "#391F2A")

ggplot(plot_data, aes(x = type, y = values_9545, fill = type)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(aes(color = type), width = 0.1, height = 0, size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Violin Plot of 9545_A>G",
       x = "Type",
       y = "Values") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  # Apply the custom colors
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = dot_colors)

ggsave("violin_plot_9545_jitter.pdf", width = 6, height = 8, dpi = 300)


# Plot coverage per position 
ymax <- 1000

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
  ylab("Mean coverage per cell") + xlab("Position along chrM") +  ggtitle('Mouse tumor, PDAC') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))


#pseudobulked expression for these two variants.

annotations <- data.frame(cells = colnames(seu), annotation = Idents(seu))
rm(seu)

variants_of_interest <- c('3010_G>A')
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'cells')
for (v in variants_of_interest) {
  message(v)
  col_name <- 'coverage'
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),coverage_info$cells]
}

for (v in variants_of_interest) {
  message(v)
  col_name <- 'allele_freq'
  coverage_info[col_name] <- af.dm[v,]
}

merged_df <- merge(annotations, coverage_info, by = "cells")



pseudobulk_weighted <- merged_df %>%
  group_by(annotation) %>%  # Replace 'cluster' with your cluster column name
  summarise(
    pseudobulk_af = sum(allele_freq * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    total_coverage = sum(coverage, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  filter(total_coverage > 0)  # Remove clusters with no coverage

p1 <- ggplot(pseudobulk_weighted, aes(x = annotation, y = pseudobulk_af)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Pseudobulk Allele Frequency by Cell Type",
    subtitle = "Weighted by coverage",
    x = "Cell Type",
    y = "Pseudobulk Allele Frequency"
  )
p1

pdf('pseudobulked_af_3010.pdf')
print(p1)
dev.off()

variants_of_interest <- c('9545_A>G')
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'cells')
for (v in variants_of_interest) {
  message(v)
  col_name <- 'coverage'
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),coverage_info$cells]
}

for (v in variants_of_interest) {
  message(v)
  col_name <- 'allele_freq'
  coverage_info[col_name] <- af.dm[v,]
}

merged_df <- merge(annotations, coverage_info, by = "cells")

pseudobulk_weighted <- merged_df %>%
  group_by(annotation) %>%  
  summarise(
    pseudobulk_af = sum(allele_freq * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    total_coverage = sum(coverage, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  filter(total_coverage > 0)

p1 <- ggplot(pseudobulk_weighted, aes(x = annotation, y = pseudobulk_af)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Pseudobulk Allele Frequency by Cell Type",
    subtitle = "Weighted by coverage",
    x = "Cell Type",
    y = "Pseudobulk Allele Frequency"
  )
p1

pdf('pseudobulked_af_9545.pdf')
print(p1)
dev.off()

# -----------------------------------------------------------------------------------------------------------------------------------------
##spatial mouse data ----------------------------------------------------------------------------------------------------------------------------

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb')
source("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB025/A200/MAESTER/a200_s6_1/analysis/210215_FunctionsGeneral.R")

seu <- readRDS(file = "/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/output/MT_seurat.rds")

## QC 
## Filtering first by log(umi) to remove empty beads
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

ggplot(plot_data, aes(x = log10_UMIs)) +
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

## Filtering by RNA count and gene count for high quality cells
Idents(seu) <- seu$orig.ident
VlnPlot(seu, features = c('nCount_RNA', 'nFeature_RNA', 'log_umi'))
seu[["percent.mt.human"]] <- PercentageFeatureSet(seu, pattern = "^GRCh38-MT-")
seu[["percent.mt.mouse"]] <- PercentageFeatureSet(seu, pattern = "^mm10---mt-")
VlnPlot(seu, feature = c('percent.mt.human', 'percent.mt.mouse'))

seu <- subset(seu, subset = percent.mt.human < 10 & percent.mt.human < 10)
seu <- subset(seu, log_umi > 1.5 & nCount_RNA > 250 & nFeature_RNA > 100)
SpatialDimPlot(seu)

#PCA             
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 5000)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))

#UMAP embedding
set.seed(101) #for reproducibility
seu <- RunUMAP(seu, dims = 1:25, n.neighbors = 25, min.dist = .1, spread = 3)
seu <- FindNeighbors(seu, dims = 1:25)
seu <- FindClusters(seu, resolution = .25)
DimPlot(seu)

## BANKSY embeddings

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 8000)
banksy_seu <- RunPCA(seu, assay = "SCT", features = VariableFeatures(object = seu), verbose = FALSE)

set.seed(101) #for reproducibility
banksy_seu <- FindNeighbors(banksy_seu, reduction = "pca", dims = 1:25)
banksy_seu <- FindClusters(banksy_seu, resolution = .25, verbose = FALSE)

banksy_seu <- RunUMAP(banksy_seu, reduction = "pca", dims = 1:25, n.neighbors = 100, min.dist = .1, spread = 3)

spatial_plt <- SpatialDimPlot(seu)
pdf('spatial_leiden_annotations.pdf')
print(spatial_plt)
dev.off()

plt <- DimPlot(seu)
pdf('UMAP_leiden_annotations.pdf')
print(plt)
dev.off()

library(spacexr)
#RCTD
analysis_folder <-'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD'

args=c("MT",'/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/output','mouse')
slide_name <- args[1] 
print(slide_name)
data_path <- args[2]
print(data_path)
ref_data= args[3] 
print(ref_data)

# load reference data
if(ref_data=="mouse"){
  reference <-readRDS("/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/10x/MT/Seurat/rctd_ref/PDAC_mouse_ref_RCTD.rds") 
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

#change mm10 to GRCm39 so naming convention matches
features$V1 <- gsub('mm10---', 'GRCm39-', features$V1 )
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
print(dim(puck@counts)) # observe Digital Gene Expression matri
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

pdf(paste0(dir_path,"diagnostic_stats.pdf"),width=4,height=4)
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI

plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

dev.off()

myRCTD <- create.RCTD(puck, reference, max_cores = 1,UMI_min = 100)
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
data <- data %>% select(-first_class, -second_class, -conv_all, -conv_doublet, -min_score, -singlet_score)
data <- data %>% filter(spot_class != 'reject') %>% select(-spot_class)

data$cell_type <- ifelse(data$type1_score > 0.90, as.character(data$first_type),
                         ifelse(data$type2_score > 0.90, as.character(data$second_type),
                                "Undetermined"))

data <- data %>% filter(cell_type != 'Undetermined')
annotations <- data %>% select(cell_type)

seu <-seu[,colnames(seu) %in% rownames(data)]
seu$annotation <- annotations
Idents(seu) <- seu$annotation

SpatialDimPlot(seu)

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD')
cluster_colors <- c("SWI1990" = "#D55E00", 
                    "CAPAN2" = "#009E73", 
                    "IM-TME" = "#E1BE6A")


spatial_plt <- SpatialDimPlot(seu, group.by = 'annotation', cols = cluster_colors)
pdf('spatial_RCTD_annotations.pdf')
print(spatial_plt)
dev.off()

write.csv(annotations, '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD/celltype_IDs.csv' )

#annotations <- read.csv('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD/celltype_IDs.csv')
#rownames(annotations) <- annotations$X
#annotations <- annotations %>% select(-X)

seu <- seu[,colnames(seu) %in% rownames(annotations)]
seu$annotation <- annotations
banksy_seu <- banksy_seu[,colnames(banksy_seu) %in% rownames(annotations)]
banksy_seu$annotation <- annotations

pdf('UMAP_RCTD_annotations.pdf')
print(DimPlot(banksy_seu, group.by = 'annotation', cols = cluster_colors))
dev.off()

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD/')
saveRDS(banksy_seu, 'banksy_MT_RCTD.rds')
saveRDS(seu, 'seurat_MT_RCTD.rds')

## MAETSER preprocessing ---------------------------------------------------------------------------------------------------------------
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER')
#Loading in data
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/maester_output/maester/maegatk.rds')
seu <- readRDS(file = '/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD/banksy_MT_RCTD.rds')

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
seu <- seu[,colnames(maegatk.rse)]

#adding into metadata percent mouse genes and percent human genes for later plotting. 
calculate_species_mixture <- function(seu, 
                                      human_prefix = "^GRCh38-", 
                                      mouse_prefix = "^mm10---") {
  # Get raw count matrix
  counts <- GetAssayData(seu, slot = "counts")
  
  # Identify human and mouse genes
  human_genes <- grep(human_prefix, rownames(counts), value = TRUE)
  mouse_genes <- grep(mouse_prefix, rownames(counts), value = TRUE)
  
  # Calculate number of genes expressed per cell (>0 counts)
  human_expressed <- Matrix::colSums(counts[human_genes,] > 0)
  mouse_expressed <- Matrix::colSums(counts[mouse_genes,] > 0)
  
  # Calculate total genes expressed per cell
  total_expressed <- Matrix::colSums(counts > 0)
  
  # Calculate percentages
  percent_human <- (human_expressed / total_expressed) * 100
  percent_mouse <- (mouse_expressed / total_expressed) * 100
  
  # Create dataframe with results
  results <- data.frame(
    cell_id = colnames(seu),
    total_genes_expressed = total_expressed,
    human_genes_expressed = human_expressed,
    mouse_genes_expressed = mouse_expressed,
    percent_human = round(percent_human, 2),
    percent_mouse = round(percent_mouse, 2)
  )
  
  # Add results to Seurat object metadata
  seu$percent_human <- results$percent_human
  seu$percent_mouse <- results$percent_mouse
  seu$total_genes_expressed <- results$total_genes_expressed
  
  # Return both the results dataframe and updated Seurat object
  return(list(
    results = results,
    seurat_object = seu
  ))
}

results <- calculate_species_mixture(seu)
results <- data.frame(cellID = results[["results"]][["cell_id"]], mouse = results[["results"]][["percent_mouse"]], human = results[["results"]][["percent_human"]])

seu$percent_mouse <- results$mouse
seu$percent_human <- results$human

saveRDS(seu, file = "MT_curio_final.rds")
saveRDS(maegatk.rse, file = "MT_curio_maegatk_final.rds")


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Mitochondrial genotyping -------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER')
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

#what is mean mito coverage per cell now? 
base.tib <- tibble(depth = colMeans(assays(maegatk.rse)[["coverage"]]))
mean(base.tib$depth)


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
write_tsv(vars.tib, "MT_allvariants.txt")
vars.tib <- read.table('MT_allvariants.txt', header = TRUE)

# Parameters to select clones with high VAF
voi.ch <- vars.tib %>% filter(mean_cov.unionCells >= 1,quality >= 21) %>%
  filter(q99.unionCells> 20) %>% .$var
voi.ch

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TWIST MT DNA sequencing validation ------------------------------------------------------------------------------------------------------------------------------------------------

## Loading in variants detected in whole MT genome sequencing for two two cell lines: 
swi1990 <- read.table('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB045/analysis/genotypes/SWI1990_unique.txt')[,2]
capan2 <- read.table('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB045/analysis/genotypes/CAPAN2_unique.txt')[,2]
twist_indeces <- c(capan2, swi1990)

# Function to remove last 4 characters, so that I can match the indeces from TWIST and MAESTER for keeping.
remove_last_4 <- function(x) {
  substr(x, 1, nchar(x) - 4)
}
maester_indeces <- as.numeric(sapply(voi.ch, remove_last_4))
matching_indices <- swi1990[swi1990 %in% maester_indeces]
matching_positions <- which(maester_indeces %in% matching_indices)

# Use these positions to select from the original voi.ch
var_of_interest <- c('3010_G>A','8002_C>T', '13500_T>C', '9545_A>G', '10176_G>A', '11152_T>C')

write.table(voi.ch, file = "MT_voi_vaf.tsv", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# lineage and coverage plots (spatial and UMAP) ------------------------------------------------------------------------------------------------------------------------------------------------------

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/lineage/umap')

for (v in var_of_interest) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  
  #coverage and AF filtering
  cells.tib$af_voi[cells.tib$cov_voi < 12]<- NA
  cells.tib$af_voi[cells.tib$af_voi < 30] <- NA
  
  #filtering out cells when plotting if there are less than 50% human genes 
  cells.tib$af_voi[cells.tib$percent_mouse > 50]<- NA
  
  p <- ggplot(data = cells.tib %>% 
                arrange(!is.na(af_voi)), # This makes NA values plot first (bottom)
              aes(x = UMAP_1, y = UMAP_2, color = af_voi)) +
    geom_point(aes(size = !is.na(af_voi)), 
               show.legend = FALSE) +
    scale_size_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/lineage/spatial/')
for (v in var_of_interest) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  cells.tib$af_voi[cells.tib$cov_voi < 12] <- NA
  cells.tib$af_voi[cells.tib$af_voi < 30] <- NA
  
  p <- ggplot(data = cells.tib %>% 
                arrange(!is.na(af_voi)), # This makes NA values plot first (bottom)
              aes(x = SPATIAL_1, y = SPATIAL_2, color = af_voi)) +
    geom_point(aes(size = !is.na(af_voi)), 
               show.legend = FALSE) +
    scale_size_manual(values = c("FALSE" = 0.1, "TRUE" = 0.9)) +
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

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/lineage/coverage')
for (v in var_of_interest) {
  message(v)
  
  # Add info of variant of interest
  cells.tib$af_voi <- af.dm[v,cells.tib$cell]
  cells.tib$cov_voi <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),cells.tib$cell]
  cells.tib$cov_voi[cells.tib$cov_voi < 10] <- NA
  
  p <- ggplot(data = cells.tib %>% 
                arrange(!is.na(af_voi)), # This makes NA values plot first (bottom)
              aes(x = SPATIAL_1, y = SPATIAL_2, color = cov_voi)) +
    geom_point(aes(size = !is.na(cov_voi)), 
               show.legend = FALSE) +
    scale_size_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
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
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# coverage per position ------------------------------------------------------------------------------------------------------------------------------------------------------
# Set y axis parameters
ymax <- 200

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

# Plot coverage along chrMe
base.tib <- tibble(base = 1:16569, depth = rowMeans(assays(maegatk.rse)[["coverage"]]))


coverage <- 
  base.tib %>% ggplot() +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)), stat = "identity", fill = "#64b53b", width = 1) + 
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = ycoord-ymax*0.2, label = cutf(Names, d = "\\.", f = 2)), size = 3) +
  ylab("Mean coverage per cell") + xlab("Position along chrM") +  ggtitle('Mouse tumor - right (spatial)') + 
  theme_classic() +
  theme(aspect.ratio = 0.5) + theme(plot.title = element_text(hjust = 0.5))

ggsave('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/MAESTER/lineage/MT_MAESTER_coverage.png', plot = coverage)

## comparing between Curio Seeker (unenriched) and SUMMIT (enriched)
setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/coverage_compare')

#reloading in data
maegatk.curio <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/maester_output/curio_sequencing/maegatk.rds')
maegatk.rse <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/maester_output/maester/maegatk.rds')
seu <- readRDS('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/seeker_comb/RCTD/banksy_MT_RCTD.rds')

seu <- subset(seu, ident = '2', invert = TRUE) #removing mouse from analyses

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
pdf('MT_coverage_curio.pdf')
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
pdf('MT_coverage_MAESTER.pdf')
print(coverage_maester)
dev.off()

#Boxplot 
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

ggsave("MT_boxplot.pdf",ratio_plot)

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
pdf('MT_coverage_compare.pdf')
print(combined_coverage)
dev.off()

## per cell enrichment?


# For maegatk.curio
total_coverage_curio <- data.frame(reads_per_cell = colSums(assays(maegatk.curio)[["coverage"]]))
mean(total_coverage_curio$reads_per_cell)

total_coverage_curio <- data.frame(reads_per_locus = rowSums(assays(maegatk.curio)[["coverage"]]))
mean(total_coverage_curio$reads_per_locus)


# For maegatk.rse (SUMMIT)
total_coverage_SUMMIT <- data.frame(all_reads = colSums(assays(maegatk.rse)[["coverage"]]))
mean(total_coverage_SUMMIT$all_reads)

## per locus enrichment?
total_coverage_curio <- data.frame(all_reads = rowSums(assays(maegatk.curio)[["coverage"]]))
mean(total_coverage_curio$all_reads)

total_coverage_SUMMIT <- data.frame(all_reads = rowSums(assays(maegatk.rse)[["coverage"]]))
mean(total_coverage_SUMMIT$all_reads)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Mito QC post processing ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# plotting RNA transcripts detected vs MAESTER coverage -----------------------------------------------------------------
transcripts <- data.frame(RNA_counts = seu$nCount_RNA, MAESTER = (colSums(assays(maegatk.rse)[["coverage"]])))

p <- ggplot(transcripts, aes(x = RNA_counts, y = MAESTER))+
  geom_point()  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "RNA transcripts detected", 
       y = "mtRNA transcripts detected")

pdf('RNA_mtRNA_scatterplot.pdf')
print(p)
dev.off()

## Violin plots to compare AF detection between spatial and sc:  
seu_2 <- subset(seu, ident = 'Mouse TME', invert = TRUE)
af_voi_3010 <- af.dm['3010_G>A',colnames(seu_2)]
af_voi_9545 <- af.dm['9545_A>G',colnames(seu_2)]

coverage_3010 <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf('3010_G>A', d = "_") ),colnames(seu_2)]
coverage_9545 <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf('9545_A>G', d = "_") ),colnames(seu_2)]

plot_data <- data.frame(
  cell = colnames(seu_2), 
  type = Idents(seu_2),
  values_3010= af_voi_3010, 
  values_9545 = af_voi_9545, 
  coverage_3010 = coverage_3010, 
  coverage_9545 = coverage_9545
)

#filtering to only plot cells that would pass earlier filtration: minimally 3 reads for either variant. 
filtered_plot_data <- plot_data %>%
  filter(coverage_3010 >= 3 | coverage_9545 >= 3)

# Create the violin plot
violin_colors <- c("Capan2" = "#133C5C", "SW1990" = "#804D60")
dot_colors <- c("Capan2" = "#1E206B", "SW1990" = "#391F2A")

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER/lineage')
ggplot(filtered_plot_data, aes(x = type, y = values_3010, fill = type)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(aes(color = type), width = 0.1, height = 0, size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Violin Plot of 3010_A>C",
       x = "Cell line",
       y = "VAF") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = dot_colors)
ggsave("violin_plot_3010_jitter.pdf", width = 6, height = 8, dpi = 300)

ggplot(filtered_plot_data, aes(x = type, y = values_9545, fill = type)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(aes(color = type), width = 0.1, height = 0, size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Violin Plot of 9545_A>G",
       x = "Cell line",
       y = "VAF") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = dot_colors)

ggsave("violin_plot_9545_jitter.pdf", width = 6, height = 8, dpi = 300)

## pseudobulked expression for these two variants.

annotations <- data.frame(cells = colnames(seu), annotation = Idents(seu))
rm(seu)

variants_of_interest <- c('3010_G>A')
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'cells')
for (v in variants_of_interest) {
  message(v)
  col_name <- 'coverage'
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),coverage_info$cells]
}

for (v in variants_of_interest) {
  message(v)
  col_name <- 'allele_freq'
  coverage_info[col_name] <- af.dm[v,]
}

merged_df <- merge(annotations, coverage_info, by = "cells")


pseudobulk_weighted <- merged_df %>%
  group_by(annotation) %>%  # Replace 'cluster' with your cluster column name
  summarise(
    pseudobulk_af = sum(allele_freq * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    total_coverage = sum(coverage, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  filter(total_coverage > 0)  # Remove clusters with no coverage

p1 <- ggplot(pseudobulk_weighted, aes(x = annotation, y = pseudobulk_af)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Pseudobulk Allele Frequency by Cell Type",
    subtitle = "Weighted by coverage",
    x = "Cell Type",
    y = "Pseudobulk Allele Frequency"
  )
p1

setwd('/Users/sydneybracht/Library/CloudStorage/GoogleDrive-sbracht@sydshafferlab.com/My Drive/Experiments/SB044/Curio/MAESTER')

pdf('pseudobulked_af_3010.pdf')
print(p1)
dev.off()

variants_of_interest <- c('9545_A>G')
coverage_info <- as_tibble(colnames(af.dm)) %>% rename(value = 'cells')
for (v in variants_of_interest) {
  message(v)
  col_name <- 'coverage'
  coverage_info[col_name] <- assays(maegatk.rse)[["coverage"]][as.numeric( cutf(v, d = "_") ),coverage_info$cells]
}

for (v in variants_of_interest) {
  message(v)
  col_name <- 'allele_freq'
  coverage_info[col_name] <- af.dm[v,]
}

merged_df <- merge(annotations, coverage_info, by = "cells")

pseudobulk_weighted <- merged_df %>%
  group_by(annotation) %>%  
  summarise(
    pseudobulk_af = sum(allele_freq * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    total_coverage = sum(coverage, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  filter(total_coverage > 0)

p1 <- ggplot(pseudobulk_weighted, aes(x = annotation, y = pseudobulk_af)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Pseudobulk Allele Frequency by Cell Type",
    subtitle = "Weighted by coverage",
    x = "Cell Type",
    y = "Pseudobulk Allele Frequency"
  )
p1

pdf('pseudobulked_af_9545.pdf')
print(p1)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------ƒs




