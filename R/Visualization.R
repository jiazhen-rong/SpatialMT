#' Plot variant-by-celltype significance as grouped bars
#'
#' Creates a grouped bar plot of -log10(p) per (Variant, Celltype), using
#' adjusted p-values from SpatialMT celltype association results.
#'
#' By default, all celltypes present in `res$adjusted_pval` are plotted.
#' If `celltypes` is provided, only those columns are plotted.
#'
#' @param res A result list from the `celltype_test()` function:
#'   \itemize{
#'     \item \code{adjusted_pval}: matrix/data.frame of p-values (variants x celltypes)
#'     \item \code{coef}: matrix/data.frame of coefficients (variants x celltypes) used to optionally mask negatives
#'   }
#' @param celltypes Character vector of celltypes to plot. If \code{NULL}, plot all.
#' @param keep_positive_only Logical; if TRUE (default), any entries with coef < 0 have p set to 1 (logp=0).
#' @param order_by Celltype name used to order variants (descending -log10(p) in that celltype).
#'   If NULL, variants keep their original order.
#' @param p_threshold Numeric p-value threshold line (default 0.05).
#' @param logy Logical; if TRUE, use log10 y-scale with special handling of zeros (default TRUE).
#' @param zero_floor Numeric; when \code{logy=TRUE}, values with logp==0 are plotted at this floor (default 1e-2).
#' @param fill_values Optional named vector of colors for celltypes (names must match celltypes).
#'   If NULL, ggplot default palette is used.
#' @param base_size Base font size for theme_minimal (default 13).
#' @param width,height PDF width/height if saving.
#' @param outfile If provided, save a PDF to this path. If NULL, no file is saved.
#' @param draw_axis_break Logical; if TRUE and \code{logy=TRUE}, draw a small '//' mark to indicate axis break.
#'
#' @return A ggplot object 
#' @import dplyr,tidyr,ggplot2
#' @export
plot_lineage_significance <- function(
    res,
    celltypes = NULL,
    keep_positive_only = TRUE,
    order_by = NULL,
    p_threshold = 0.05,
    logy = TRUE,
    zero_floor = 1e-2,
    fill_values = NULL,
    base_size = 13,
    width = 10,
    height = 5,
    outfile = NULL,
    draw_axis_break = TRUE,
    return_plot=T,
    title=NULL
) {
  
  if(is.null(res$adjusted_pval)){
    stop("No p values detected in the input.")
  }
  
  pval_filter=res_lg$adjusted_pval
  pval_filter[res_lg$coef <0] = 1
  pval_df <- as.data.frame(as.matrix(pval_filter))
  pval_df$Variant <- rownames(pval_filter)
  
  # Long format
  df_long <- pval_df %>%
    tidyr::pivot_longer(
      cols = -Variant,
      names_to = "Celltype",
      values_to = "pval"
    )
  
  # Subset celltypes if requested
  if (!is.null(celltypes)) {
    df_long <- df_long %>%
      filter(Celltype %in% celltypes)
  }
  
  # Summarize (kept for safety, but usually already one value each)
  df_summary <- df_long %>%
    dplyr::group_by(.data$Variant, .data$Celltype) |>
    dplyr::summarise(logp = -log10(.data$pval), .groups = "drop")
  
  # Order variants by chosen celltype
  if (!is.null(order_by)) {
    if (!(order_by %in% df_summary$Celltype)) {
      warning(sprintf("order_by='%s' not found in plotted celltypes; keeping original order.", order_by))
    } else {
      ord <- df_summary %>%
        dplyr::filter(.data$Celltype == order_by) %>%
        dplyr::arrange(dplyr::desc(.data$logp)) %>%
        dplyr::pull(.data$Variant)
      
      # Keep only variants present
      df_summary$Variant <- factor(df_summary$Variant, levels = unique(ord))
    }
  }
  
  # Consider for zeros for log y-scale plotting
  if (logy) {
    df_summary <- df_summary %>%
      dplyr::mutate(
        logp_plot = ifelse(.data$logp == 0, zero_floor, .data$logp),
        is_zero = .data$logp == 0,
        logp_label = ifelse(is_zero, "0", sprintf("%.2f", logp)),# label 0 as "0", others as "%.2f"
         signif_label = case_when(
          logp > 3 ~ "***",
          logp > 2 ~ "**",
          logp > -log10(0.05) ~ "*",
          TRUE ~ ""
        )
      )
    y_aes <- ggplot2::aes(y = .data$logp_plot)
  } else {
    df_summary$logp_plot <- df_summary$logp
    y_aes <- ggplot2::aes(y = .data$logp)
  }
    
  # Build plot
  p <- ggplot(df_summary,aes(x = .data$Variant, fill = .data$Celltype)) +
    geom_bar(
      y_aes,
      stat = "identity",
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.7
    ) +
    labs(x = "Variant", y = expression(-log[10](p)), fill = NULL) +
    ggtitle(title) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.line.x = element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,  face = "bold", size = base_size + 2)
    ) + 
    # significance stars
    geom_text(
      aes(label = signif_label, y = logp_plot + 0.3),
      position = position_dodge(width = 0.7),
      size = 5
    ) +
    # Threshold line
    geom_hline(
      yintercept = -log10(p_threshold),
      linetype = "dashed",
      color = "black"
    ) 
  
  # Fill override
  if (!is.null(fill_values)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_values)
  }
  
  # Log y scale
  if (logy) {
    p <- p + ggplot2::scale_y_log10(
      limits = c(zero_floor, NA),
      breaks = c(zero_floor, 1e-1, 1e0, 1e1, 1e2),
      labels = c("0", "1e-1", "1", "10", "100")
    ) +
      ggplot2::geom_hline(yintercept = 1, color = "black", linewidth = 0.5)
  }
  
  # Save or return
  if (!is.null(outfile)) {
    pdf(outfile, width = width, height = height)
    if (logy && draw_axis_break) {
      # draw with grid to add the // mark
      grid::grid.newpage()
      grid::grid.draw(ggplot2::ggplotGrob(p))
      grid::grid.rect(
        x = grid::unit(0.15, "npc"),
        y = grid::unit(0.3, "npc"),
        width = grid::unit(0.008, "npc"),
        height = grid::unit(0.05, "npc"),
        gp = grid::gpar(col = NA, fill = "white")
      )
      grid::grid.text(
        "//",
        x = grid::unit(0.15, "npc"),
        y = grid::unit(0.3, "npc"),
        gp = grid::gpar(fontsize = 16, fontface = "bold")
      )
    } else {
      print(p)
    }
    dev.off()
  }
  if(return_plot){
    return(p)
  }
  #invisible(p)
}