library(grid)
library(gridExtra)
library(spacexr)
library(CARD)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(SpatialMT)
library(patchwork)

setwd("~/nzhanglab/project/jrong/mito_LT/scripts/our_model/")
source("example_data/210215_FunctionsGeneral.R") # from MAESTER paper
source("diagnostic_plot.R")

# load mito & RCTD results from spatial data
# (1) mitochondrial data
maegatk.rse = readRDS("../../data/Sydney_Bracht/mouse_cellline_mix/MAESTER/MT-R_maegatk_final.rds") # master result object
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)

# load variant of interest
voi = read.table("../../data/Sydney_Bracht/mouse_cellline_mix/MAESTER/MT_voi_vaf.tsv",sep="\t")[,1]
voi = c("3010_G>A","9545_A>G","8002_C>T","11152_T>C","10176_G>A","13500_T>C")
#voi = c("3010_G>A",rownames(N)[grepl("10176",rownames(N))],
#        rownames(N)[grepl("11152",rownames(N))])# from TWIST
#voi=c("3010_G>A","9545_A>G")
#subset data to the variant of interest
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]

# remove "-1" in the spot barcodes
#colnames(af.dm) = sapply(strsplit(colnames(af.dm),"-"),"[[",1)
#colnames(N) = sapply(strsplit(colnames(N),"-"),"[[",1)

# spatial coordinates 
spatial_obj <- readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/mouse_cellline_mix/MAESTER/MT-R_seurat_final.rds")
spatial_obj <- subset(spatial_obj, nFeature_RNA > 20)
spatial_coords = as.data.frame(spatial_obj@reductions$SPATIAL@cell.embeddings)

# Load celltype ratios
# RCTD ratio
rctd_res = readRDS("../../results/Sydney_Bracht/mouse_cellline_mix/spatial/RCTD/DBremoved/RCTD_res.rds")
Ws = as.data.frame(as.matrix(normalize_weights(rctd_res@results$weights))) 
celltypes= colnames(Ws)

library(ggnewscale)
temp_df= Ws#as.data.frame(cbind(Ws,t(N[c("3010_G>A","9545_A>G"),])))
celltpyes =colnames(Ws)

# Create a new column 'color_cat' based on the conditions.
temp_df$color_cat <- with(temp_df, ifelse(SWI1990 <= 0.05 & CAPAN2 <= 0.05, "None",
                                          ifelse(SWI1990 <= 0.05, "CAPAN2",
                                                 ifelse(CAPAN2 <= 0.05, "SW1990", "Mixed"))))
temp_df$color_cat2 <- with(temp_df, ifelse(SWI1990 <= 0.05 & EMT <= 0.05, "None",
                                          ifelse(SWI1990 <= 0.05, "EMT",
                                                 ifelse(EMT <= 0.05, "SW1990", "Mixed"))))
temp_df$color_cat3 <- with(temp_df, ifelse(CAPAN2 <= 0.05 & EMT <= 0.05, "None",
                                           ifelse(CAPAN2 <= 0.05, "EMT",
                                                  ifelse(EMT <= 0.05, "CAPAN2", "Mixed"))))
table(temp_df$color_cat)
table(temp_df$color_cat2)
table(temp_df$color_cat3)

# Create the scatter plot.
p1 <- ggplot(temp_df, aes(x = SWI1990, y = CAPAN2, color = color_cat)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("SW1990" = "#e56126",#"#956173",
                                "CAPAN2" = "#9e2064",#"#114D6E",
                                "None" = "black",
                                "Mixed" = "grey"),
                     name = "Celltype") +
  labs(x = "SW1990", 
       y = "CAPAN2", 
       title = "Celltype Ratio Scatter Plot") +
  # Add a vertical dashed line at x=0.05 (SWI threshold)
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#9e2064", size = 1) +
  # Add a horizontal dashed line at y=0.05 (CAPAN2 threshold)
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#e56126", size = 1) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")  # this keeps axis lines
        )
p2 <- ggplot(temp_df, aes(x = SWI1990, y = EMT, color = color_cat2)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("SW1990" = "#e56126",
                                "EMT" = "#d7df27",#"#BA952F",
                                "None" = "black",
                                "Mixed" = "grey"),
                     name = "Celltype") +
  labs(x = "SW1990", 
       y = "EMT", 
       title = "Celltype Ratio Scatter Plot") +
  # Add a vertical dashed line at x=0.05 (SWI threshold)
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#e56126", size = 1) +
  # Add a horizontal dashed line at y=0.05 (CAPAN2 threshold)
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d7df27", size = 1) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
        )  # this keeps axis lines)
p3 <- ggplot(temp_df, aes(x = EMT, y = CAPAN2, color = color_cat3)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("EMT" = "#d7df27",
                                "CAPAN2" = "#9e2064",
                                "None" = "black",
                                "Mixed" = "grey"),
                     name = "Celltype") +
  labs(x = "EMT", 
       y = "CAPAN2", 
       title = "Celltype Ratio Scatter Plot") +
  # Add a vertical dashed line at x=0.05 (SWI threshold)
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#9e2064", size = 1) +
  # Add a horizontal dashed line at y=0.05 (CAPAN2 threshold)
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d7df27", size = 1) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black")  # this keeps axis lines
        )

pdf("cell_mix/celltype_ratio_scatter.pdf",height=5,width=15)
print(p1+p2+p3)
dev.off()

###### Linear Regression
save_path = "cell_mix/test_package/"
# overall celltype test for voi
res_lg = celltype_test(celltypes = celltypes,voi=voi,N=N,
                       vaf=af.dm,Ws=Ws,spatial_coords = spatial_coords,
                       test_type = "linear",plot=T,save_path=save_path,method="FDR")
saveRDS(res_lg,paste0(save_path,"/res_lg.rds"))

res_lg = readRDS(paste0(save_path,"/res_lg.rds"))

# visualize
library(tidyverse)
pval_filter=res_lg$adjusted_pval
pval_filter[res_lg$coef <0] = 1
pval_df <- as.data.frame(pval_filter)
pval_df$Variant <- rownames(pval_filter)

library(tidyr)
# Pivot to long format
df_long <- pval_df %>%
  pivot_longer(-Variant, names_to = "Celltype", values_to = "pval")

# Compute group label
df_long <- df_long %>%
  mutate(Group = case_when(
    Celltype == "SWI1990" ~ "SWI1990",
    Celltype == "CAPAN2" ~ "CAPAN2",
    TRUE ~ "Others"
  ))

# Average others per variant and combine
df_summary <- df_long %>%
  group_by(Variant, Group) %>%
  summarize(logp = -log10(mean(pval)), .groups = "drop") %>%
  mutate(Group = factor(Group, levels = c("SWI1990", "CAPAN2", "Others")))

# Optional: order variants by SWI1990 signal
df_summary$Variant <- factor(df_summary$Variant,
                             levels = df_summary %>%
                               filter(Group == "SWI1990") %>%
                               arrange(desc(logp)) %>%
                               pull(Variant))

pdf("cell_mix/test_package/lineage_significance.pdf",width=5,height=5)
ggplot(df_summary, aes(x = Variant, y = logp, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("SWI1990" = "#956173", "CAPAN2" = "#114D6E",
                               "Others" = "grey")) +
  labs(x = "Variant", y = expression(-log[10](p)), fill = NULL) +
  coord_cartesian(ylim = c(-0.1, NA)) +  # ensures y starts from 0
  theme_minimal(base_size = 13) +
  #scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = "grey50")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# Replace 0 with a small value for plotting, keep original for labels
df_summary <- df_summary %>%
  mutate(
    logp_plot = ifelse(logp == 0, 1e-2, logp),
    is_zero = logp == 0,
    logp_label = ifelse(is_zero, "0", sprintf("%.2f", logp))  # label 0 as "0", others as "%.2f"
  )

# Plot
p <- ggplot(df_summary, aes(x = Variant, y = logp_plot, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  # Annotate values
  #geom_text(aes(label = logp_label),
  #          position = position_dodge(width = 0.7),
  #          vjust = -0.3, size = 3.5) +
  
  # Styling
  scale_fill_manual(values = c("SWI1990" = "#956173", 
                               "CAPAN2" = "#114D6E", 
                               "Others" = "grey")) +
  
  # Y-axis log scale with custom breaks and labels
  scale_y_log10(
    limits = c(1e-2, NA),
    breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2),
    labels = c("0", "1e-1", "1", "10", "100")
  ) +

  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "Variant", y = expression(-log[10](p)), fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.y = element_line(color = "black"),
    panel.grid.minor.y = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = "grey50")
  )

# Draw base plot
grid.newpage()
grid.draw(ggplotGrob(p))
grid.rect(x = unit(0.15, "npc"), 
          y = unit(0.3, "npc"), 
          width = unit(0.008, "npc"), 
          height = unit(0.05, "npc"),
          gp = gpar(col = NA, fill = "white"))
# Draw zigzag "//" on y-axis to simulate break
grid.text("//",
          x = unit(0.15, "npc"),  # position right over y-axis line
          y = unit(0.3, "npc"),   # visually between 0.005 and 0.01 on log scale
          gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()

sig_var = filter_significant_variants(res_lg,alpha_threshold = 0.05,
                                      output_file="cell_mix/test_package/sig_vars_fdr.txt")
  #as.data.frame((res_lg$coef > 0) & (res_lg$pval <0.05))
# beta test for voi
beta_list = power_analysis_all(voi=voi,celltypes = celltypes,Ws=Ws,
                               N=N,vaf=af.dm,X=NULL,
                               sample_num=100,alpha=0.05,n_sim=10,
                               beta_threshold = 0.5,plot=T,save_path = save_path)
saveRDS(beta_list,paste0(save_path,"/beta_list.rds"))
# ways of visualizing beta values
effect_sizes = seq(0.1, 1, by = 0.2)
# Filter to keep only the 3 desired cell types
selected_names <- grep("SWI1990|CAPAN2|EMT", names(beta_list), value = TRUE)
# Convert to data frame
power_df <- do.call(rbind, lapply(selected_names, function(name) {
  power_values <- beta_list[[name]]
  lineage <- sub("_(SWI1990|CAPAN2|EMT)$", "", name)
  celltype <- sub("^.*_(SWI1990|CAPAN2|EMT)$", "\\1", name)
  data.frame(
    lineage = lineage,
    celltype = celltype,
    effect_size = effect_sizes,
    power = power_values,
    stringsAsFactors = FALSE
  )
}))
power_df_sub=power_df[(power_df$lineage %in% c("3010_G>A","9545_A>G"))
                      & (power_df$celltype %in% c("SWI1990","CAPAN2","EMT")),]

power_df2 <- subset(power_df, celltype != "EMT")
desired_order <- voi
power_df2$lineage <- factor(power_df2$lineage, levels = desired_order)
p2<-ggplot(power_df2, aes(x = effect_size, y = power, 
                     color = celltype, shape = celltype, linetype = celltype)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  facet_wrap(~ lineage) +
  #geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("SWI1990" = "#956173", 
                                "CAPAN2" = "#114D6E", 
                                "EMT" = "grey")) +
  scale_shape_manual(values = c("SWI1990" = 16, "CAPAN2" = 17, "EMT" = 15)) +
  scale_linetype_manual(values = c("SWI1990" = "solid", "CAPAN2" = "dashed", "EMT" = "dotdash")) +
  labs(
    x = "Effect Size", y = "Power", color = "Cell Type",
    title = "Power Analysis"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
    #strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.7),
    #strip.text = element_text(face = "bold", color = "black")
    panel.grid.major    = element_blank(),     # remove major grid
    panel.grid.minor    = element_blank(),     # remove minor grid
    strip.background    = element_blank(),     # remove facet header box
    #axis.line.x        = element_line(color = "black"),   # keep x-axis line
    #axis.line.y        = element_line(color = "black")  
    #panel.border        = element_rect(color = "black", fill = NA),
  )
pdf(paste0(save_path,"/power_visualization.pdf"),width=8,height=5)
#print(p1)
print(p2)
dev.off()

# plot side by side
df_summary <- df_summary %>%
  mutate(
    Group = recode_factor(Group,
                          "SWI1990"   = "SWI1990",
                          "CAPAN2"    = "CAPAN2",
                          "Others"    = "Mouse TME"),
    # replace zeros with 1e-2 for plotting on log scale
    logp_plot = ifelse(logp == 0, 1e-2, logp)
  )

make_variant_panel <- function(v) {
  # significance barplot
  # prepare the threshold
  p_thresh <- -log10(0.05)
  p_sig <- df_summary %>%
    filter(Variant == v) %>%
    ggplot(aes(x = Group, y = logp_plot, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = 1, color = "black", size = 0.5) +
    geom_hline(yintercept = p_thresh, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c(
      "SWI1990"   = "#956173",
      "CAPAN2"    = "#114D6E",
      "Mouse TME" = "#ba952f"
    )) +
    scale_y_log10(
      limits = c(1e-2, max(df_summary$logp_plot)),
      breaks = c(1e-2, 1e-1, 1, 10, 100),
      labels = c("0", "1e-1", "1", "10", "100")
    ) +
    labs(title=v,
         x = "",
         y = "-log10p") +#expression(-log[10](p))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),#element_line(color = "grey80"),
      panel.grid.minor   = element_blank(),
      axis.line          = element_blank(),  
      axis.line.y.left   = element_line(color = "black"),
      axis.ticks         = element_line(color = "black"),
      axis.text.x        = element_text(angle = 45, hjust = 1),
      legend.position    = "none"
    )
  
  # power curve (unchanged)
  p_pow <- power_df2 %>%
    filter(lineage == v) %>%
    ggplot(aes(x = effect_size, y = power,
               color = celltype, shape = celltype, linetype = celltype)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    ylim(c(0,1))+
    scale_color_manual(values = c(
      "SWI1990" = "#956173",
      "CAPAN2"  = "#114D6E",
      "EMT"     = "#ba952f"
    )) +
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
    labs(x = "Effect Size", y = "Power") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      legend.position   = "none"
    )
  
  # stack & title
  (p_sig / p_pow) + plot_annotation(title = v)
}

# Build and save 2×3
panels    <- lapply(voi, make_variant_panel)
final_plt <- wrap_plots(panels, ncol = 6)
pdf("cell_mix/test_package/linegaes_sig_power_stack.pdf",width=15,height=5)
print(final_plt)
dev.off()
###### Diagnostic Plots
X= round(N*af.dm)
cell_label = spatial_obj@active.ident
cell_types= c("SWI1990","CAPAN2","IM-SM-ECM")
colors = setNames(c("#D46127","#059e74","#E1BD6A"),cell_types)
shapes = setNames(c(1,2,3),cell_types)
sizes = setNames(c(2,2,1),cell_types)

library(cowplot)
library(grid)
library(gridExtra)
pdf(paste0(save_path,"/diagnostic_plots.pdf"),height=10,width=10)
for(var in voi){
  p1<- plot_variant_diagnostics(X[var,],N[var,],af.dm[var,],
                                cell_label,var,colors,shapes)
  print(p1)
}
dev.off()
pdf(paste0(save_path,"/spatial_vaf_diagnostic.pdf"),height=5,width=20)
for(var in voi){
  p2 <- plot_spatial_vaf(X[var,],N[var,],af.dm[var,],spatial_coords,var)
  print(p2)
}
dev.off()

# Grid based Scan statistics
voi_sub = c("3010_G>A","9545_A>G","8002_C>T","11152_T>C","10176_G>A","13500_T>C")
# selecting spots as control
spatial_coords=spatial_coords
voi=voi
pdf(paste0(save_path,"/grid_diagnostics.pdf"),height=20,width=20)
for(var in voi_sub){
  sample_idx = intersect(rownames(Ws)[(Ws$SWI1990 < 0.2)&(Ws$CAPAN2 < 0.2)],
                         colnames(N)[N[var,] > 4])
  results_grid <-celltype_test_grid(
    celltypes=colnames(Ws),
    voi=c(var),
    N=N[var,,drop=F],
    vaf=af.dm[var,,drop=F],
    Ws=Ws,
    spatial_coords = spatial_coords,
    test_type="linear",
    grid_size=10,
    verbose=T,
    vaf_cellprop=F,
    #disease_celltype = "SWI1990",
    sample_idx = sample_idx,
    method="FDR"
  )
  p1=plot_pval_grid_per_celltype(results_grid$results,
                                 coef_plot_option="negative",method="FDR")
  p2=plot_pval_grid_per_celltype(results_grid$results,
                                 coef_plot_option="grey",method="FDR")
  print(p1)
  print(p2)
  # debug
  #grid.draw(results_grid_5$vaf_cellprop_plots[["Grid_1_1"]])
}
dev.off()

