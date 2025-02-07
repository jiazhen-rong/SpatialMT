#  helper functions of the spatial mito model

# plotting VAF vs celltype proportions, colored by coverage
plot_vaf_cellprop <- function(i,af.dm,norm_weights,plot_df,intercept,coef,pval,permute=F){
  if(permute==T){
    p_str = "permute_p"
  }else{
    p_str = "ANOVA_p"
  }
  # specify a variant
  var = rownames(af.dm)[i]
  p_list = list()
  for(k in 1:length(colnames(norm_weights))){
    celltype=colnames(norm_weights)[k]
    fit_df=data.frame(x=seq(0,1,0.01),y=intercept[k]+coef[k]*seq(0,1,0.01))
    p_list[[celltype]]  = local({
      celltype=celltype
      var=var
      var_cov=paste0(var,"_cov")
      temp=plot_df[order(plot_df[,var],decreasing=F),]
      #ggplot(temp %>% arrange(get(var_cov)),
      #aes(x=get(celltype),y=get(var),color=get(var_cov))) + #,size=get(var_cov))) +
      ggplot(temp %>% arrange(get(var_cov))) +
        geom_point(size=1,alpha=0.5,aes(x=get(celltype),y=get(var),color=get(var_cov))) +
        geom_line(data=fit_df,aes(x=x,y=y)) +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =  element_text(angle = 90,vjust=0.5,hjust=1))+
        scale_color_continuous(low="grey",high="red") +
        ggtitle(paste0(celltype," ",p_str,": ",formatC(pval[k],format="e",digits=2))) +
        xlab("Celltype Proportions") +  ylab("VAF") +
        guides(color=guide_legend(title="Coverage"))#,size=guide_legend(title="Coverage"))
    })

  }
  # to speed up, can choose to save in pdf
  grid.arrange(grobs=p_list,top=var)
}

#' Plot Variant Diagnostic Plots
#'
#' Generates diagnostic plots for single-cell mitochondrial variant analysis.
#' This function produces three types of plots for each cell type:
#' (1) X vs N (Alternative Allele Count vs Total Count)
#' (2) VAF vs N (Variant Allele Frequency vs Total Count)
#' (3) Histogram of the square root of VAF
#' Additionally, it includes combined plots for all cell types.
#'
#' @param X A named vector or matrix containing the alternative allele counts per cell.
#' @param N A named vector or matrix containing the total coverage per cell.
#' @param af_matrix A named vector or matrix containing the variant allele frequencies per cell.
#' @param cell_label A named factor or character vector assigning cell types to each barcode.
#' @param variant_name A string specifying the variant being analyzed (e.g., "3054_G>C").
#' @param colors A named vector of colors for each cell type.
#' @param shapes A named vector of shapes for each cell type.
#'
#' @return A `ggplot2` object with a grid of diagnostic plots for each cell type and combined plots.
#'
#' @details
#'
#' - The last row of the grid includes combined plots across all cell types.
#' - If a cell type has no matching barcodes, it is skipped.
#'
#' @examples
#' # Define Colors and Shapes
#' cell_types <- c("BE", "SQ", "NE", "Rest")
#' colors <- setNames(c("red", "chartreuse4", "black", "orange"), cell_types)
#' shapes <- setNames(c(16, 17, 18, 19), cell_types)  # Different point shapes
#'
#' # Run function on an example variant
#' plot_variant_diagnostics(X["3054_G>C", ], N["3054_G>C", ], af.dm["3054_G>C", ],
#'                          cell_label, "3054_G>C", colors, shapes)
#'
#' @export

plot_variant_diagnostics <- function(X, N, af_matrix, cell_label, variant_name, colors,shapes) {
  unique_types <- unique(cell_label)
  plot_list <- list()

  # Ensure barcodes match across inputs
  common_cells <- intersect(names(cell_label), names(N))
  print(length(common_cells))
  # Handle both matrix and vector cases
  if (is.matrix(N)) {
    N <- N[, common_cells, drop = FALSE]
    X <- X[, common_cells, drop = FALSE]
    af_matrix <- af_matrix[, common_cells, drop = FALSE]
  } else {
    N <- N[common_cells]
    X <- X[common_cells]
    af_matrix <- af_matrix[common_cells]
  }
  cell_label <- cell_label[common_cells]

  # Get Global Axis Limits
  max_N <- max(N, na.rm = TRUE)  # Max total count
  max_X <- max(X, na.rm = TRUE)  # Max alternative allele count
  max_AF <- max(af_matrix, na.rm = TRUE)  # Max VAF
  min_N <- min(N, na.rm = TRUE)

  # Generate plots for each cell type
  for (cell_type in unique_types) {
    # Ensure subsetting is correct
    matching_cells <- which(cell_label == cell_type)
    # If no matching cells, skip to avoid errors
    if (length(matching_cells) == 0) {
      message(paste("Skipping", cell_type, "because no matching cells found."))
      next
    }

    df <- data.frame(
      Total_Count = N[matching_cells],
      Alt_Allele_Count = jitter(X[matching_cells], 0.25),
      VAF = af_matrix[matching_cells],
      Cell_Type = cell_type
    )

    # X vs N Plot
    p1 <- ggplot(df, aes(x = Total_Count, y = Alt_Allele_Count,shape = Cell_Type, color = Cell_Type)) +
      geom_point(color = colors[cell_type]) +
      labs(title = paste(variant_name, "-", cell_type), x = "Total Count", y = "Alt Allele Count") +
      scale_shape_manual(values = shapes) +
      theme_minimal() +
      xlim(min_N, max_N) + ylim(0, max_X)

    # VAF vs N Plot
    p2 <- ggplot(df, aes(x = Total_Count, y = VAF,shape = Cell_Type, color = Cell_Type)) +
      geom_point(color = colors[cell_type]) +
      labs(title = paste(variant_name, "-", cell_type), x = "Total Count", y = "VAF") +
      scale_shape_manual(values = shapes) +
      theme_minimal() +
      ylim(c(0, 1)) +
      xlim(min_N, max_N) + ylim(0, 1)

    # Histogram of VAF
    p3 <- ggplot(df, aes(x = VAF)) +
      geom_histogram(fill = colors[cell_type], bins = 100, alpha = 0.7, boundary = 0) +
      scale_y_continuous(trans = "sqrt") +
      labs(title = paste(variant_name, "-", cell_type), x = "VAF", y = "SQRT of Count") +
      theme_minimal() +
      xlim(c(0, 1))

    # Store plots in list
    plot_list[[paste(cell_type, "XvsN")]] <- p1
    plot_list[[paste(cell_type, "VAFvsN")]] <- p2
    plot_list[[paste(cell_type, "HistVAF")]] <- p3
  }


  # **Final Row: Combined Plot for All Cell Types**
  df_all <- data.frame(
    Total_Count = N,
    Alt_Allele_Count = jitter(X, 0.25),
    VAF = af_matrix,
    Cell_Type = cell_label
  )

  p1_all <- ggplot(df_all, aes(x = Total_Count, y = Alt_Allele_Count, shape = Cell_Type, color = Cell_Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(title = paste(variant_name, "- All Cell Types"), x = "Total Count", y = "Alt Allele Count") +
    theme_minimal() +
    xlim(0, max_N) + ylim(0, max_X)

  p2_all <- ggplot(df_all, aes(x = Total_Count, y = VAF, shape = Cell_Type, color = Cell_Type)) +
    geom_point(size = 2) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(title = paste(variant_name, "- All Cell Types"), x = "Total Count", y = "VAF") +
    theme_minimal() +
    xlim(0, max_N) + ylim(0, 1)

  plot_list[["All_XvsN"]] <- p1_all
  plot_list[["All_VAFvsN"]] <- p2_all

  # Arrange all plots in a grid
  if (length(plot_list) > 0) {
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 3)
  } else {
    message("No valid plots to display.")
    combined_plot <- NULL
  }

  return(combined_plot)
}


#' Plot VAF, Alternative Allele Counts, and Total Counts in Spatial Coordinates
#'
#' This function visualizes spatial distributions of Variant Allele Frequency (VAF),
#' Alternative Allele Counts (X), and Total Counts (N) across spots in spatial transcriptomics data.
#'
#' @param X A named vector of alternative allele counts per spot.
#' @param N A named vector of total coverage counts per spot.
#' @param af_matrix A named vector of variant allele frequencies per spot.
#' @param spatial_coords A matrix or data frame with spot coordinates (rows = spots, columns = X and Y positions).
#' @param variant_name A string specifying the variant being analyzed (e.g., "3054_G>C").
#'
#' @return A `ggplot2` object showing three spatial plots: VAF, X, and N.
#'
#' @examples
#' plot_spatial_vaf(X["3054_G>C", ], N["3054_G>C", ], af.dm["3054_G>C", ],
#'                  spatial_coords, "3054_G>C")
#'
#' @export
plot_spatial_vaf <- function(X, N, af_matrix, spatial_coords, variant_name) {

  # Ensure spot names are aligned across all data inputs
  common_spots <- intersect(names(X), rownames(spatial_coords))
  X <- X[common_spots]
  N <- N[common_spots]
  af_matrix <- af_matrix[common_spots]
  spatial_coords <- spatial_coords[common_spots, ]

  # Convert data to a plot-friendly format
  df <- data.frame(
    X_Coord = spatial_coords[, 1],
    Y_Coord = spatial_coords[, 2],
    VAF = af_matrix,
    Alt_Allele_Count = X,
    Total_Count = N
  )

  # Define gradient color palettes
  vaf_palette <- scale_color_gradient(low = "grey", high = "red", limits = c(0, 1))
  count_palette <- scale_color_gradient(low = "grey", high = "blue")

  # Helper function for dual-layer plotting
  highlight_plot <- function(data, value_col, title, color_palette) {
    ggplot(data, aes(x = X_Coord, y = Y_Coord)) +
      geom_point(color = "grey", size = 0.5) +  # All points in grey
      geom_point(data = subset(data, get(value_col) > 0), aes(color = get(value_col)), size = 1) +  # Highlighted points
      color_palette +
      labs(title = title, x = "X Coord", y = "Y Coord") +
      theme_minimal()
  }

  # Generate plots
  p1 <- highlight_plot(df, "VAF", paste("VAF of", variant_name), vaf_palette)
  p2 <- highlight_plot(df, "Alt_Allele_Count", paste("Alt Allele Count of", variant_name), count_palette)
  p3 <- highlight_plot(df, "Total_Count", paste("Total Count of", variant_name), count_palette)

  # Arrange plots in a grid
  combined_plot <- plot_grid(p1, p2, p3, ncol = 3)

  return(combined_plot)
}

