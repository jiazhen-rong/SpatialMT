#' Plot global variant-by-celltype significance as grouped bars
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
  #invisible(p)
  if(return_plot){
    return(p)
  }
}

#' Plot power curves for multiple cell types across lineages
#'
#' This function extracts power vectors from a named list (beta_list) whose names
#' end with one of the specified cell type suffixes (e.g., "_BE", "_IM", ...),
#' converts them into a long data frame, and draws power curves vs. effect size
#' faceted by lineage.
#'
#' @param beta_list Named list. Each entry is a numeric vector of power values,
#'   aligned with `effect_sizes`. Names must follow: "<lineage>_<celltype>".
#' @param effect_sizes Numeric vector of effect sizes (x-axis).
#' @param celltype_pattern Regex for cell type suffixes to include.
#'   Default matches patterns like "BE|IM|SQ|FB|VC".
#' @param desired_order Optional character vector giving the facet order for lineages.
#'   If NULL, ordering is whatever appears in the data.
#' @param celltype_colors Named character vector mapping celltype -> hex color.
#'   If NULL, uses the default palette for BE/SQ/IM/FB/VC.
#' @param facet_ncol Number of columns in facet_wrap. If NULL, ggplot decides.
#' @param title Plot title.
#' @param save_pdf Logical; if TRUE, save a PDF.
#' @param save_path Directory to save the PDF when save_pdf=TRUE.
#' @param file_name Output PDF filename when save_pdf=TRUE.
#' @param width,height PDF dimensions in inches.
#' @param return_plot Return the plot
#'
#' @return 
#'   - plot: a ggplot object
#'
#' @examples
#' res <- plot_power_curves(beta_list, effect_sizes,
#'   desired_order = c("3054_G>C","3071_T>C","15777_G>C"),
#'   save_pdf = TRUE, save_path = save_path
#' )
#' @import dplyr,tidyr,ggplot2
#' @export
plot_power_curves <- function(
    beta_list,
    effect_sizes,
    celltype_pattern = "BE|IM|SQ|FB|VC",
    desired_order = NULL,
    color_map = NULL,
    facet_ncol = NULL,
    title = "Power Analysis",
    save_pdf = FALSE,
    save_path = NULL,
    file_name =NULL,
    width = 8,
    height = 5,
    return_plot=T
) {
  
  # 1) Select entries by suffix
  suffix_regex <- paste0("(", celltype_pattern, ")$")
  selected_names <- grep(suffix_regex, names(beta_list), value = TRUE)
  if (length(selected_names) == 0) {
    stop("No entries in beta_list matched suffix pattern: ", celltype_pattern)
  }
  
  # 2) Build long data frame
  power_df <- do.call(rbind, lapply(selected_names, function(name) {
    power_values <- beta_list[[name]]
    
    if (!is.numeric(power_values)) {
      stop("beta_list[['", name, "']] is not numeric.")
    }
    if (length(power_values) != length(effect_sizes)) {
      stop("Length mismatch for '", name, "': power_values has length ",
           length(power_values), " but effect_sizes has length ", length(effect_sizes), ".")
    }
    lineage <- sub(paste0("_", suffix_regex), "", name)
    celltype <- sub(paste0("^.*_",suffix_regex), "\\1", name)
    
    data.frame(
      lineage = lineage,
      celltype = celltype,
      effect_size = effect_sizes,
      power = power_values,
      stringsAsFactors = FALSE
    )
  }))
  
  # 3) Optional lineage ordering
  if (!is.null(desired_order)) {
    power_df$lineage <- factor(power_df$lineage, levels = desired_order)
  }
  
  # 4) Plot
  p <- ggplot(power_df, aes(x = effect_size, y = power,
      color = celltype, shape = celltype, linetype = celltype
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    facet_wrap(~ lineage, ncol = facet_ncol) 
  
    if (is.null(celltype_colors)) {
    } else {
      p <- p + scale_color_manual(values = celltype_colors)
    }
  
    p <- p + labs(
      x = "Effect Size", y = "Power", color = "Cell Type",
      title = title
    ) +
    ylim(c(0, 1)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(
        fill = "grey95", color = "grey40", linewidth = 0.7
      ),
      strip.text = element_text(face = "bold", color = "black")
    ) +
    guides(
      color = guide_legend(title = "Cell Type"),
      shape = guide_legend(title = "Cell Type"),
      linetype = guide_legend(title = "Cell Type")
    )
  
  if (save_pdf) {
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    out_file <- file.path(save_path, file_name)
    pdf(out_file, width = width, height = height)
    print(p)
    dev.off()
  }
    
  if(return_plot){
    return(p)
  }
}

