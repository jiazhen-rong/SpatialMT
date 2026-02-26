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
    dplyr::group_by(.data$Variant, .data$Celltype) %>%
    dplyr::summarise(pval=.data$pval, logp = -log10(.data$pval), .groups = "drop")

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
  #if (logy) {
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
  if (logy) {
    y_aes <-  aes(y = .data$logp_plot)
  }else{
    y_aes <-  aes(y = .data$logp)
    # }
  }
  #}
  # } else {
  #   df_summary$logp_plot <- df_summary$logp
  #   y_aes <-  aes(y = .data$logp)
  # }

  # Build plot
  p <- ggplot(df_summary,aes(x = .data$Variant, fill = .data$Celltype)) +
    geom_bar(
      y_aes,
      stat = "identity",
      position =  position_dodge(width = 0.7),
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
      plot.title =   element_text(
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
    p <- p +   scale_fill_manual(values = fill_values)
  }

  # Log y scale
  if (logy) {
    p <- p +   scale_y_log10(
      limits = c(zero_floor, NA),
      breaks = c(zero_floor, 1e-1, 1e0, 1e1, 1e2),
      labels = c("0", "1e-1", "1", "10", "100")
    ) +
        geom_hline(yintercept = 1, color = "black", linewidth = 0.5)
  }

  # Save or return
  if (!is.null(outfile)) {
    pdf(outfile, width = width, height = height)
    if (logy && draw_axis_break) {
      # draw with grid to add the // mark
      grid::grid.newpage()
      grid::grid.draw(  ggplotGrob(p))
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
    celltype_colors = NULL,
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


plot_localized_test <- function(
    voi,
    save_path,
    Ws,
    celltypes = NULL,             # default: all Ws colnames
    signed = TRUE,
    adjust_method = "BH",
    fdr_thresh = 0.05,
    celltype_colors = NULL,
    file_prefix = "spatial_maxlogp",
    width_all = 8, height_all = 4
) {
  library(dplyr); library(ggplot2); library(patchwork)

  stopifnot(!is.null(voi), length(voi) > 0)
  stopifnot(!is.null(save_path))
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

  # Default celltypes
  all_ct <- colnames(Ws)
  if (is.null(celltypes)) {
    celltypes <- all_ct
  }

  # Check provided celltypes exist in Ws
  if (!all(celltypes %in% colnames(Ws))) {
    stop("Some provided celltypes are not columns of Ws: ",
         paste(setdiff(celltypes, colnames(Ws)), collapse = ", "))
  }

  # container
  var_spot_list <- list()

  for (var in voi) {
    f <- file.path(save_path, paste0(var, "_spot_test.rds"))
    if (!file.exists(f)) {
      stop("Missing file: ", f, "\nMake sure your spot test files are saved as <save_path>/<VAR>_spot_test.rds")
    }
    var_spot_res <- readRDS(f)

    # Build signed pval matrix: rows = celltypes, cols = barcodes
    var_spot_pval <- sapply(names(var_spot_res$results), function(bc) {
      coef <-var_spot_res$results[[bc]]$coef
      pval <- var_spot_res$results[[bc]]$pval

      if (signed) {
        pval_signed <- pval
        pval_signed[!is.na(coef) & coef < 0] <- 1
        return(pval_signed)
      } else {
        return(pval)
      }
    })

    var_spot_list[[var]] <- list()
    var_spot_list[[var]]$var_spot_pval <- var_spot_pval

    # per-variant FDR adjustment per celltype (across spots)
    var_spot_list[[var]]$fdr <- list()
    var_spot_list[[var]]$fdr$thresh <- fdr_thresh
    var_spot_list[[var]]$fdr$adjusted_pval <- t(apply(var_spot_pval, 1, function(p) p.adjust(p, method = adjust_method)))

    # count significant spots per celltype (after FDR)
    spot_sig_num_fdr <- sapply(1:length(all_ct), function(i) {
      sum(var_spot_list[[var]]$fdr$adjusted_pval[i, ] <= fdr_thresh, na.rm = TRUE)
    })
    print(length(spot_sig_num_fdr))
    names(spot_sig_num_fdr) <- all_ct
    var_spot_list[[var]]$fdr$spot_sig_num <- spot_sig_num_fdr

    # store list of significant spots per celltype
    var_spot_list[[var]]$fdr$sig_spots <- setNames(vector("list", length(celltypes)), celltypes)
    for (i in seq_along(all_ct)) {
      celltype_i <- all_ct[i]
      sig_spots_fdr <- colnames(var_spot_pval)[var_spot_list[[var]]$fdr$adjusted_pval[i, ] <= fdr_thresh]
      var_spot_list[[var]]$fdr$sig_spots[[celltype_i]] <- sig_spots_fdr
    }
  }

  # Aggregate per-celltype
  df_min_p <- lapply(voi, function(var) {
    M <- var_spot_list[[var]]$var_spot_pval
    M_adj <- var_spot_list[[var]]$fdr$adjusted_pval

    data.frame(
      variant = var,
      celltype = all_ct,
      min_pval = apply(M, 1, min, na.rm = TRUE),
      min_adjusted_pval = apply(M_adj, 1, min, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  df_min_p =do.call(rbind,df_min_p)

  df_min_p <- df_min_p %>%
    dplyr::mutate(
      variant = factor(variant, levels = voi),
      neglog10_p = -log10(min_adjusted_pval),
      logp_plot = ifelse(is.infinite(neglog10_p) | neglog10_p == 0, 1e-2, neglog10_p),
      signif_label = dplyr::case_when(
        neglog10_p > 3 ~ "***",
        neglog10_p > 2 ~ "**",
        neglog10_p > -log10(fdr_thresh) ~ "*",
        TRUE ~ ""
      )
    )

  # after computing full df_min_p (with all celltypes), create subset
  if (!is.null(celltypes)) {
    # keep only requested subtypes that exist
    keep_ct <- intersect(celltypes, unique(df_min_p$celltype))
    df_min_p <- df_min_p %>% dplyr::filter(celltype %in% keep_ct)
  }
  # and then use df_min_p_sub for plotting instead of df_min_p

  # plotting helper: all celltypes per variant
  max_plot_p = max(df_min_p$logp_plot)
  make_variant_plot <- function(df, var_name) {
    subdf <- df %>% dplyr::filter(variant == var_name)
    p <- ggplot(subdf, aes(x = factor(celltype, levels = celltypes), y = logp_plot, fill = celltype)) +
         geom_col(width = 0.7) +
         scale_y_log10(limits = c(1e-2, max_plot_p ),
                             breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2),
                             labels = c("0", "1e-1", "1", "10", "100")) #+
    if(!is.null(celltype_colors)){p <- p+ scale_fill_manual(values = celltype_colors[celltypes], drop = FALSE)}
    p <- p +
       geom_text( aes(label = signif_label, y = logp_plot + 0.3), size = 4) +
       geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black", size = 0.4) +
       geom_hline(yintercept = 1, color = "black", linewidth = 0.6)+
       labs(title = var_name, x = NULL, y = expression(-log[10](p))) +
       theme_minimal(base_size = 13) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1),
             axis.line.y =    element_line(color = "black"),
             panel.grid.major =    element_blank(),
             panel.grid.minor =    element_blank(),
             legend.position = "bottom",
             axis.line = element_line(color = "black")
             )
    p
  }

  plot_all_obj <- NULL
  plots <- lapply(unique(df_min_p$variant), function(v) make_variant_plot(df_min_p, v))
  plot_all_obj <- patchwork::wrap_plots(plots, ncol = length(plots)) +
  patchwork::plot_annotation(title = "Minimum spot p-value (per celltype)")
  ggsave(filename = file.path(save_path, paste0(file_prefix, "_celltypes.pdf")),
         plot = plot_all_obj, width = width_all, height = height_all)

  # Save aggregated data and var_spot_list
  saveRDS(df_min_p, file.path(save_path, paste0(file_prefix, "_df_min_p.rds")))
  saveRDS(var_spot_list, file.path(save_path, paste0(file_prefix, "_var_spot_list.rds")))

  return(list(df_min_p = df_min_p,
              var_spot_list = var_spot_list,
              plot_all = plot_all_obj))
}
