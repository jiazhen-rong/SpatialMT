#' Calculate Differentially Expressed Genes Between Two Groups of Cells
#'
#' This function performs differential expression analysis between two groups of cells,
#' defined by their cell barcodes. It creates a temporary metadata column in the Seurat object
#' to assign cells to either "Lineage1" or "Lineage2" based on the provided barcode lists,
#' subsets the object to only include those cells, and then runs Seurat's \code{FindMarkers()}
#' to identify differentially expressed genes.
#'
#' @param seu A Seurat object containing gene expression data.
#' @param lineage1_cells A character vector of cell barcodes belonging to the first group.
#' @param lineage2_cells A character vector of cell barcodes belonging to the second group.
#' @param assay A character string specifying which assay to use for DE analysis (default "RNA").
#' @param ... Additional arguments passed to \code{FindMarkers()} (e.g., test.use, logfc.threshold, etc.).
#'
#' @return A data frame containing the differentially expressed genes between the two groups.
#'
#' @examples
#' \dontrun{
#'   # Suppose you have two sets of cell barcodes from your AF.dm/N analysis:
#'   degs <- calculate_DEGs_from_cell_barcodes(seu, 
#'             lineage1_cells = c("Cell1", "Cell2", "Cell3"),
#'             lineage2_cells = c("Cell10", "Cell11", "Cell12"),
#'             assay = "RNA",
#'             test.use = "wilcox", logfc.threshold = 0.25)
#'   head(degs)
#' }
#'
#' @export
calculate_DEGs_from_cell_barcodes <- function(seu, lineage1_cells, lineage2_cells, assay = "RNA", ...) {
  # Create a temporary metadata column to hold lineage assignments
  all_cells <- colnames(seu)
  lineage_temp <- rep(NA, length(all_cells))
  names(lineage_temp) <- all_cells
  
  lineage_temp[all_cells %in% lineage1_cells] <- "Lineage1"
  lineage_temp[all_cells %in% lineage2_cells] <- "Lineage2"
  
  seu$lineage_temp <- lineage_temp
  
  # Subset Seurat object to only include cells that are assigned to a lineage
  # Identify the cells to keep (those with a non-NA assignment)
  cells_to_keep <- names(lineage_temp)[!is.na(lineage_temp)]
  # Subset the Seurat object by cell names
  seu_subset <- subset(seu, cells = cells_to_keep)
  
  # Set the identities to the temporary lineage assignment
  Idents(seu_subset) <- seu_subset$lineage_temp
  
  # Run differential expression analysis
  degs <- FindMarkers(seu_subset, ident.1 = "Lineage1", ident.2 = "Lineage2", assay = assay)
  
  return(degs)
}


#' Calculate and Plot Gene Set Enrichment Between Two Lineages
#'
#' This function calculates gene set enrichment scores between two lineages in a Seurat object.
#' It takes two sets of cell barcodes (one for each lineage) and a named list of gene sets (e.g., cancer hallmark gene sets).
#' The function computes module scores using Seurat's \code{AddModuleScore} for each lineage, averages the module scores,
#' and then plots a heatmap of the average scores with gene sets as rows and lineages as columns.
#'
#' @param seu A Seurat object containing gene expression data.
#' @param lineage1_cells A character vector of cell barcodes for lineage 1.
#' @param lineage2_cells A character vector of cell barcodes for lineage 2.
#' @param gene_programs A named list of gene sets (e.g., hallmark gene sets) where each element is a vector of gene symbols.
#' @param assay A character string specifying the assay to use for module scoring (default "RNA").
#' @param name_prefix A character string for the prefix used by \code{AddModuleScore} (default "Hallmark").
#' @param ... Additional arguments passed to \code{AddModuleScore} (e.g., nbin).
#'
#' @return A list with two elements:
#'   \item{avg_scores}{A matrix of average module scores (gene programs as rows, lineages as columns).}
#'   \item{heatmap}{A pheatmap object displaying the enrichment heatmap.}
#'
#' @examples
#' \dontrun{
#'   # Assuming 'seu' is your Seurat object and hallmark_list is your named list of gene sets:
#'   res <- calculate_gene_set_enrichment_between_lineages(
#'            seu, 
#'            lineage1_cells = l1_spots, 
#'            lineage2_cells = l2_spots, 
#'            gene_programs = hallmark_list, 
#'            assay = "RNA"
#'          )
#'   res$heatmap
#' }
#'
#' @export
calculate_gene_set_enrichment_between_lineages <- function(seu, lineage1_cells, lineage2_cells, 
                                                           gene_programs, assay = "RNA", 
                                                           name_prefix = "Hallmark", ...) {
  # Subset the Seurat object for each lineage
  seu_l1 <- subset(seu, cells = lineage1_cells)
  seu_l2 <- subset(seu, cells = lineage2_cells)
  seu_l3 <- subset(seu, cells = lineage3_cells)
  
  # Ensure each gene set in gene_programs is a character vector
  #gene_programs <- lapply(gene_programs, as.character)
  
  # Calculate module scores for gene programs in each subset.
  # AddModuleScore creates columns like "Hallmark1", "Hallmark2", etc.
  seu_l1 <- AddModuleScore(seu_l1, features = gene_programs, name = names(gene_programs), assay = assay,nbin=10)
  seu_l2 <- AddModuleScore(seu_l2, features = gene_programs, name = names(gene_programs), assay = assay,nbin=10)
  seu_l3 <- AddModuleScore(seu_l3, features = gene_programs, name = names(gene_programs), assay = assay,nbin=10)
  
  # Determine the number of gene sets
  n_sets <- length(gene_programs)
  score_names <- paste0(names(gene_programs), 1:n_sets)#paste0(name_prefix, 1:n_sets)
  
  # Fetch the module scores from each subset
  scores_l1 <- FetchData(seu_l1, vars = score_names)
  scores_l2 <- FetchData(seu_l2, vars = score_names)
  scores_l3 <- FetchData(seu_l3, vars = score_names)
  
  # Compute average module scores for each gene set in each lineage
  avg_l1 <- colMeans(scores_l1, na.rm = TRUE)
  avg_l2 <- colMeans(scores_l2, na.rm = TRUE)
  avg_l3 <- colMeans(scores_l3, na.rm = TRUE)
  
  # Combine into a matrix: rows = gene sets, columns = two lineages
  avg_scores <- rbind(Lineage1 = avg_l1, Lineage2 = avg_l2 )#, lineage3=avg_l3)
  avg_scores <- t(avg_scores)  # Now rows = gene sets
  
  # Replace row names with the original gene set names if available
  if (!is.null(names(gene_programs))) {
    rownames(avg_scores) <- names(gene_programs)
  }
  
  # Load pheatmap package and plot the heatmap with a blue-grey-red palette
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install the 'pheatmap' package.")
  }
  library(pheatmap)
  colnames(avg_scores) = c("15777_G>C","3054_G>C")
  heatmap_obj <- pheatmap(avg_scores, 
                          color = colorRampPalette(c("blue", "grey", "red"))(100),
                          cluster_rows = TRUE, 
                          cluster_cols = FALSE,
                          cellwidth = 15,      # Ensure cells are square
                          cellheight = 10,
                          show_rownames = TRUE,
                          show_colnames = TRUE,
                          main = "Gene Set Enrichment Between Lineages")
  
  return(list(avg_scores = avg_scores, heatmap = heatmap_obj))
}