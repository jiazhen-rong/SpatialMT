library(ggplot2)
library(ggnewscale)  
library(patchwork)

setwd("~/nzhanglab/project/jrong/mito_LT/scripts/simulation_experiment/")
source("simulation_functions.R")

# 1. synthetic grids
clone_size=50
clones <- list(
  list(id = 1, x_start = 3, y_start = 6, width = 10, height = 5),
  list(id = 2, x_start = 21, y_start = 11, width = 10, height = 5),
  list(id = 3, x_start = 11, y_start = 21, width = 10, height = 5)
)
#clone_colors = c("orange", "red", "blue", "lightgrey")
clone_colors = c("#e56126", "#9e2064", "#d7df27", "lightgrey")
plot_df = simulate_clone_location(grid_size=30,clones,
                                  plot=F,clone_colors =clone_colors,
                                  save=T,
                                  save_path="data/clone_3_size_50/",
                                  width=10,height=10)

celltype_ratio=0.5;cell_range="Custom"
save_path = "data/clone_3_size_50/"
W <- readRDS(paste0(save_path,"/sim_celltype_ratio_",
                  celltype_ratio,"_",cell_range,".rds"))

plot_df <- cbind(plot_df, as.data.frame(W))

# simulated celltype
p1 <- ggplot() +
  # 1) “Other” spots: fill by C4 (white→dark grey), no outline
  geom_point(
    data   = subset(plot_df, label == "Other"),
    aes(
      x     = x_coord,
      y     = y_coord,
      fill  = C4
    ),
    shape  = 21,
    size   = 1.5,
    stroke = 0        # removes the border entirely
  ) +
  scale_fill_gradient(
    name     = "Other",
    low      = "white",
    high     = "darkgrey",
    limits   = c(0, 1),
    oob      = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 4)
  ) +
  
  # 2) Clone_1: fill by C1 (white→orange), no outline
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, label == "Clone_1"),
    aes(
      x     = x_coord,
      y     = y_coord,
      fill  = C1
    ),
    shape  = 21,
    size   = 1.7,
    stroke = 0        # removes the border
  ) +
  scale_fill_gradient(
    name     = "C1",
    low      = "white",
    high     = "#e56126",
    limits   = c(0, 1),
    oob      = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 1)
  ) +
  
  # 3) Clone_2: fill by C2 (white→red), no outline
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, label == "Clone_2"),
    aes(
      x     = x_coord,
      y     = y_coord,
      fill  = C2
    ),
    shape  = 21,
    size   = 1.7,
    stroke = 0        # removes the border
  ) +
  scale_fill_gradient(
    name     = "C2",
    low      = "white",
    high     = "#9e2064",
    limits   = c(0, 1),
    oob      = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 2)
  ) +
  
  # 4) Clone_3: fill by C3 (white→blue), no outline
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, label == "Clone_3"),
    aes(
      x     = x_coord,
      y     = y_coord,
      fill  = C3
    ),
    shape  = 21,
    size   = 1.7,
    stroke = 0        # removes the border
  ) +
  scale_fill_gradient(
    name     = "C3",
    low      = "white",
    high     = "#d7df27",
    limits   = c(0, 1),
    oob      = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 3)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1) +
  labs(
    x     = "",
    y     = "",
    title = "Celltype Ratio"
  )

# plot spatial vaf
vf=0.5
n=5
X_mtx <- t(readRDS(paste0("data/clone_3_size_50/vaf_",vf,"_N_",n,".rds")))
vaf = X_mtx/n
plot_df$var1 = vaf[,1]
plot_df$var2 = vaf[,6]
plot_df$var3 = vaf[,11]

p2<-ggplot() +
  # 1) Background: only “Other” spots in light grey
  geom_point(
    data   = subset(plot_df, label == "Other"),
    aes(
      x = x_coord,
      y = y_coord
    ),
    shape  = 21,
    size   = 1.5,
    fill   = "lightgrey",
    color  = NA     # no border
  ) +
  
  # 2) Clone_1: fill by var1 (white→orange)
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, label == "Clone_1"),
    aes(
      x    = x_coord,
      y    = y_coord,
      fill = var1
    ),
    shape  = 21,
    size   = 1.7,
    stroke = 0       # no border
  ) +
  scale_fill_gradient(
    name   = "C1",
    low    = "white",
    high   = "blue",
    limits = c(0, 1),
    oob    = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 1)
  ) +
  
  # 3) Clone_2: fill by var2 (white→red)
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, label == "Clone_2"),
    aes(
      x    = x_coord,
      y    = y_coord,
      fill = var2
    ),
    shape  = 21,
    size   = 1.7,
    stroke = 0
  ) +
  scale_fill_gradient(
    name   = "C2",
    low    = "white",
    high   = "#00fff0",
    limits = c(0, 1),
    oob    = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 2)
  ) +
  
  # 4) Clone_3: fill by var3 (white→blue)
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, label == "Clone_3"),
    aes(
      x    = x_coord,
      y    = y_coord,
      fill = var3
    ),
    shape  = 21,
    size   = 1.7,
    stroke = 0
  ) +
  scale_fill_gradient(
    name   = "C3",
    low    = "white",
    high   = "skyblue",
    limits = c(0, 1),
    oob    = scales::squish,
    na.value = NA,
    guide = guide_colorbar(order = 3)
  ) +
  
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x     = "",
    y     = "",
    title = "Spatial VAF"
  )

pdf("results/fig2_syngrid.pdf",width=10,height=5)
#fig <- (p1 + p2) + plot_layout(guides = "collect")
#fig <- fig & theme(legend.position  = "bottom",legend.direction = "horizontal")
p1v <- p1 + theme(
  legend.position  = "bottom",
  legend.direction = "vertical",
  legend.key.height = unit(0.5, "cm")
)
p2v <- p2 + theme(
  legend.position  = "bottom",
  legend.direction = "vertical",
  legend.key.height  = unit(0.5, "cm")
)

# combine without merging legends:
fig_separate <- p1v + p2v + plot_layout(ncol = 2, guides = "keep")

print(fig_separate)
dev.off()

# Seperately plot the clones
library(ggplot2)

make_spatial_vaf_plot <- function(df, var, clone_color, title) {
  ggplot() +
    # background “Other” spots
    geom_point(
      data   = subset(df, label == "Other"),
      aes(x = x_coord, y = y_coord),
      shape  = 21, size = 1.5,
      fill   = "lightgrey",
      color  = NA
    ) +
    # only the relevant clone’s variant spots
    geom_point(
      data   = subset(df, label == paste0("Clone_", which(c("var1","var2","var3")==var))),
      aes_string(x = "x_coord", y = "y_coord", fill = var),
      shape  = 21, size = 1.7, stroke = 0
    ) +
    scale_fill_gradient(
      name   = title,
      low    = "white",
      high   = clone_color,
      limits = c(0, 1),
      oob    = scales::squish,
      na.value = NA
    ) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.direction = "horizontal"
    ) +
    labs(x = NULL, y = NULL)
}

p_var1 <- make_spatial_vaf_plot(plot_df, "var1", "blue", "Variant 1 (C1)")
p_var2 <- make_spatial_vaf_plot(plot_df, "var2", "#00fff0",    "Variant 2 (C2)")
p_var3 <- make_spatial_vaf_plot(plot_df, "var3", "skyblue",   "Variant 3 (C3)")

# Arrange them side-by-side (or in any layout you like):
library(patchwork)
fig <- (p_var1 | p_var2 | p_var3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",  # stack the separate legends one above the next
    legend.direction = "horizontal"  # arrange each legend’s keys in a column
  )

pdf("results/fig2_syngrid.pdf",width=10,height=5)
print(fig)
dev.off()


## Plot synthetic data based on SUMMIT
# Simulate coordinates
plot_df <- readRDS("results/plot_df_withClones_A200S6.rds")
# 3) Load the simulated VAF/N matrices for clone_size = 100, coherence = 0 (i.e. sparsity = 0)
data = readRDS("data/sim_real/vaf_0.5_clonesize_100_coherence_0.rds")
vaf_sim = t(data$X_mtx/(data$N+0.01))
plot_df$var1 <- vaf_sim[, 1]  
plot_df$var2 <- vaf_sim[, 2]   
plot_df$var3 <- vaf_sim[, 3] 
plot_df$carrier = "Non-carrier"
plot_df$carrier[plot_df$celltype_rctd=="BE"] = "Carrier"


# 6) Example: Plot the spatial VAF for “clone1_var1” as a continuous heatmap:
library(ggplot2)

# ── (A) First, mark carriers vs. non‐carriers in plot_df
plot_df$carrier <- "Non‐carrier"
plot_df$carrier[plot_df$celltype_rctd == "BE"] <- "Carrier"
rctd_res = readRDS("../../results/Sydney_Bracht/sample_level/a200_s6_high_cov/RCTD/Major/RCTD_res.rds")
Ws = as.data.frame(as.matrix(normalize_weights(rctd_res@results$weights)))
plot_df$BE_ratio = NA;
plot_df$BE_ratio = Ws[rownames(plot_df),"BE"]

p3 <- ggplot() +
  # 1) draw Non‐carriers in solid grey
  geom_point(
    data = subset(plot_df, carrier == "Non‐carrier"),
    aes(x = -spatial_2, y = -spatial_1),
    shape  = 21,
    size   = 0.5,
    fill   = "grey",
    color  = NA    # no outline
  ) +
  
  # 2) draw Carriers and map fill = BE_ratio
  geom_point(
    data = plot_df,#subset(plot_df, carrier == "Carrier"),
    aes(x = -spatial_2, y = -spatial_1, fill = BE_ratio),
    shape  = 21,
    size   = 0.5,
    stroke = 0     # remove any border
  ) +
  scale_fill_gradient(
    name  = "BE ratio",
    low   = "grey",
    high  = "red",
    limits = c(0, 1),
    oob   = scales::squish,
    na.value = NA
  ) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(    panel.border     = element_rect(color="black", fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5))+
  labs(title = "Carrier Ratio")

library(ggplot2)
library(ggnewscale)   # for new_scale_fill()

# ── Plot all three variants on the same spatial canvas:
p4<- ggplot() +
  # 1) Draw spots where var1, var2, var3 are all zero in light grey
  geom_point(
    data = subset(plot_df, var1 == 0 & var2 == 0 & var3 == 0),
    aes(x = -spatial_2, y = -spatial_1),
    shape  = 21,
    size   = 0.5,
    fill   = "lightgrey",
    #color  = NA,
  ) +
  
  # 2) Overlay spots with var1 > 0 (white→red)
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, var1 > 0),
    aes(
      x    = -spatial_2,
      y    = -spatial_1,
      fill = var1
    ),
    shape  = 21,
    size   = 0.7,
    stroke = 0
  ) +
  scale_fill_gradient(
    name   = "Variant 1 VAF",
    low    = "white",
    high   = "blue",
    limits = c(0, 1),
    oob    = scales::squish,
    na.value = NA,
    guide=guide_colorbar(order=1)
  ) +
  
  # 3) Overlay spots with var2 > 0 (white→green)
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, var2 > 0),
    aes(
      x    = -spatial_2,
      y    = -spatial_1,
      fill = var2
    ),
    shape  = 21,
    size   = 0.5,
    stroke = 0
  ) +
  scale_fill_gradient(
    name   = "Variant 2 VAF",
    low    = "white",
    high   = "#00fff0",
    limits = c(0, 1),
    oob    = scales::squish,
    na.value = NA,
    guide=guide_colorbar(order=2)
  ) +
  
  # 4) Overlay spots with var3 > 0 (white→blue)
  new_scale_fill() +
  geom_point(
    data   = subset(plot_df, var3 > 0),
    aes(
      x    = -spatial_2,
      y    = -spatial_1,
      fill = var3
    ),
    shape  = 21,
    size   = 0.5,
    stroke = 0
  ) +
  scale_fill_gradient(
    name   = "Variant 3 VAF",
    low    = "white",
    high   = "skyblue",
    limits = c(0, 1),
    oob    = scales::squish,
    na.value = NA,
    guide=guide_colorbar(order=3)
  ) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(    panel.border     = element_rect(color="black", fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5))+
  labs(title = "Spatial VAF")

# 3) Bottom‐aligned vertical legends for each plot
common_legend_theme <- theme(
  legend.position   = "bottom",
  legend.direction  = "horizontal",  # keys left→right
  legend.box        = "vertical",    # stack multiple legends top→bottom
  legend.title      = element_text(size = 10),
  legend.text       = element_text(size = 8),
  legend.key.height = unit(0.4, "cm"),
  plot.title        = element_text(hjust = 0.5)
)

p1v <- p1   + common_legend_theme
p2v <- p2   + common_legend_theme
p3v <- p3 + common_legend_theme
p4v <- p4 + common_legend_theme

# 4) Arrange in 2×2, each with its own legend:
update_geom_defaults("point", list(shape = 21, stroke = 0))
library(grid)
library(cowplot)
library(patchwork)

fig_final <- (p1v + p2v) /
  (p3v + p4v) +
  plot_layout(    ncol   = 2,
                  nrow   = 2,
                  guides = "keep",
                  align  = "hv",
                  heights = c(1, 1)    # both rows get equal share of vertical space
  )

common_margin <- theme(
  plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
)
p1r <- p1v + common_margin
p2r <- p2v + common_margin
p3r <- p3v + common_margin
p4r <- p4v + common_margin

# now grid them and align both horizontally & vertically
fig_aligned <- plot_grid(
  p1r, p2r,
  p3r, p4r,
  ncol   = 2,
  align  = "hv",
  axis   = "tblr"   # align top, bottom, left & right axes
)
pdf("results/fig2_syndata_stacked.pdf",width=10,height=10)
print(fig_aligned)
dev.off()

#pdf("scripts/simulation_experiment/results/fig2_syn_realdata.pdf",width=10,height=5)
#print(fig_final)
#dev.off()
make_spatial_vaf_plot <- function(df, var, clone_color, title) {
  ggplot() +
    # background “Other” spots
    geom_point(
      data = subset(df, var1 == 0 & var2 == 0 & var3 == 0),
      aes(x = -spatial_2, y = -spatial_1),
      shape  = 21, size = 1.5,
      fill   = "lightgrey",
      color  = NA
    ) +
    # only the relevant clone’s variant spots
    geom_point(
      data   = subset(df, 
                      label == paste0("Clone_", which(c("var1","var2","var3")==var))),
      aes_string(x = "x_coord", y = "y_coord", fill = var),
      shape  = 21, size = 1.7, stroke = 0
    ) +
    scale_fill_gradient(
      name   = title,
      low    = "white",
      high   = clone_color,
      limits = c(0, 1),
      oob    = scales::squish,
      na.value = NA
    ) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.direction = "horizontal"
    ) +
    labs(x = NULL, y = NULL)
}


p_var1 <- make_spatial_vaf_plot(plot_df, "var1", "blue", "Var1")
p_var2 <- make_spatial_vaf_plot(plot_df, "var2", "#00fff0","Var2")
p_var3 <- make_spatial_vaf_plot(plot_df, "var3", "skyblue","Var3")

# Arrange them side-by-side (or in any layout you like):
library(patchwork)
fig <- (p_var1 | p_var2 | p_var3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",  # stack the separate legends one above the next
    legend.direction = "horizontal"  # arrange each legend’s keys in a column
  )

pdf("results/fig2_realsummit.pdf",width=10,height=5)
print(fig)
dev.off()

##### 
# load metrics
grid_metric = readRDS("results/metric_list.rds")

df_plot <- grid_metric %>%
  filter(vaf == 0.25, celltype_ratio == 0.25) %>%
  mutate(
    clone_size    = as.numeric(as.character(clone_size)),
    N             = as.numeric(as.character(N)),
    weighted_f1   = as.numeric(as.character(weighted_f1))
  )

# Now draw clone_size vs. weighted_f1, coloring by N (coverage) with an orange gradient:
p1<-ggplot(df_plot, aes(x = clone_size, y = weighted_f1, group = factor(N), color = N)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_gradient(
    name   = "Coverage (N)",
    low    = "grey",
    high   = "black",#"orange",
    limits = range(df_plot$N),
    oob    = scales::squish
  ) +
  labs(
    x     = "Clone Size",
    y     = "Weighted F1 Score",
    title = "F1 Score vs. Clone Size"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

grid_metric_clean <- grid_metric %>%
  mutate(
    clone_size     = as.numeric(as.character(clone_size)),
    N              = as.numeric(as.character(N)),
    vaf            = as.numeric(as.character(vaf)),
    celltype_ratio = as.numeric(as.character(celltype_ratio)),
    weighted_f1    = as.numeric(as.character(weighted_f1))
  )
df_plot <- grid_metric_clean %>%
  filter(celltype_ratio == 0.25, N == 5)

# Plot VAF (x) vs. weighted_f1 (y), colored by clone_size using an orange gradient
p2 <- ggplot(df_plot, aes(x = vaf, y = weighted_f1, color = clone_size, group = clone_size)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_gradient(
    name   = "Clone Size",
    low    = "grey",
    high   = "black",#"blue",
    limits = range(df_plot$clone_size),
    oob    = scales::squish
  ) +
  labs(
    x     = "VAF",
    y     = "Weighted F1 Score",
    title = "F1 vs. VAF"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p1<- p1+ 
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black")
p2<- p2 + 
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black")

pdf("results/grid_metric.pdf",width=10,height=5)
print(p1+p2)
dev.off()

############
real_metric = readRDS("results/sim_real/metric_list.rds")
# 1) Ensure real_metric columns are numeric
real_metric_clean <- real_metric %>%
  mutate(
    clone_size   = as.numeric(as.character(clone_size)),
    vaf          = as.numeric(as.character(vaf)),
    coherence    = as.numeric(as.character(coherence)),
    weighted_f1  = as.numeric(as.character(weighted_f1))
  )

# 2) Filter for vaf == 0.25
df3 <- real_metric_clean %>%
  filter(vaf == 0.5)

# 3) Plot p3: x = clone_size, y = weighted_f1, color = coherence (sparsity)
p3 <- ggplot(df3, aes(x = clone_size, y = weighted_f1,
                      color = coherence, group = coherence)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_gradient(
    name = "Sparsity\n(coherence)",
    low  = "lightgrey",
    high = "black",#"red",
    limits = range(df3$coherence),
    oob    = scales::squish
  ) +
  labs(
    x     = "Clone Size",
    y     = "Weighted F1 Score",
    title = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+ 
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black")

# 4) Plot p4: x = coherence, y = weighted_f1, color = clone_size
df4 <- real_metric_clean %>%
  filter(coherence == 0)

# df4 <- real_metric_clean %>%
#   filter(vaf == 0.25)
p4 <- ggplot(df4, aes(x = vaf, y = weighted_f1,
                      color = clone_size, group = clone_size)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_gradient(
    name = "Clone Size",
    low  = "grey",
    high = "black",#"blue",
    limits = range(df3$clone_size),
    oob    = scales::squish
  ) +
  labs(
    x     = "VAF",
    y     = "Weighted F1 Score",
    title = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+ 
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black")

# 5) Display the plots
library(gridExtra)
library(patchwork)
pdf("results/fig2b_metrics.pdf",width=10,height=10)
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
(p1 + p2) / (p3 + p4)
dev.off()


#(6) plot minimum p-value & spot level significance
library(dplyr)
library(ggplot2)
library(patchwork)

make_minp_barplot <- function(var_spot_pval, celltypes, title){
  df_min <- data.frame(
    celltype = celltypes,
    min_pval = apply(var_spot_pval, 1, min, na.rm = TRUE)
  ) %>%
    mutate(
      neglog10 = -log10(min_pval),
      neglog10 = ifelse(neglog10==0, 1e-2, neglog10),
    )
  
  ggplot(df_min, aes(x = celltype, y = neglog10, fill = celltype)) +
    geom_col() +
    scale_y_log10(
      name   = "-log10 pmin",
      limits = c(1e-2, NA),
      breaks = c(1e-2, 1e-1, 1, 10, 100),
      labels = c("0",   "1e-1",  "1",  "10",  ">100"),
      expand = c(0, 0)
    ) +
    geom_hline(yintercept = 1, color = "black", size = 0.5) +
    labs(title = title, x = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.line.y       = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# helper to make a spot-levell significance map:
make_spatial_sigmap <- function(plot_df, pvals, title){
  df <- plot_df %>%
    mutate(
      pval   = pvals,
      logp   = -log10(pval),
      logp   = ifelse(is.infinite(logp), NA, logp)
    )
  
  ggplot(df, aes(x = -spatial_2, y = -spatial_1, color = logp)) +
    geom_point(size = 1) +
    scale_color_viridis_c(
      name   = expression(-log[10](p)),
      option = "D",
      na.value = "grey80"
    ) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(hjust = 0.5)
    )
}


# ──────────────────────────────────────────────────────────────────────────────
# 1) GRID data
# ──────────────────────────────────────────────────────────────────────────────
vf <- 0.5; n <- 5
# (a) read in the spatial coords & labels you generated for the grid
clone_size=50
clones <- list(
  list(id = 1, x_start = 3, y_start = 6, width = 10, height = 5),
  list(id = 2, x_start = 21, y_start = 11, width = 10, height = 5),
  list(id = 3, x_start = 11, y_start = 21, width = 10, height = 5)
)
clone_colors = c("orange", "red", "blue", "lightgrey")
plot_df_grid= simulate_clone_location(grid_size=30,clones,
                                  plot=F,clone_colors =clone_colors,
                                  save=T,
                                  save_path="data/clone_3_size_50/",
                                  width=10,height=10)
colnames(plot_df_grid) = c("spatial_1","spatial_2","label")
# (b) load =per‐spot test results
X_mtx = readRDS("data/clone_3_size_50/vaf_0.5_N_5.rds")
colnames(X_mtx) = c(1:900)
N = X_mtx; N[,]=5; af.dm = X_mtx/N;
celltype_ratio=0.5;cell_range="Custom"
save_path = "data/clone_3_size_50/"
Ws <- readRDS(paste0(save_path,"/sim_celltype_ratio_",
                    celltype_ratio,"_",cell_range,".rds"))
colnames(Ws) = c("Clone_1","Clone_2","Clone_3","Other");
rownames(Ws) = c(1:900)
var="Var_1"; save_path ="results/clone_3_size_50/vaf_0.5_N_10_ratio_0.5_range_Custom/"
spatial_coords=plot_df_grid
celltypes=unique(plot_df_grid$label)
library(SpatialMT)
sample_idx=rownames(Ws)[Ws[,1] <0.25]
exclude_plot_idx = rownames(Ws)[Ws[,1] <0.25]
result_var=celltype_test_knn(celltypes,vars=c(var),
                             N_voi=N[var,,drop=F],vaf=af.dm[var,,drop=F],
                             Ws,spatial_coords,
                             test_type="linear",k_neighbors=50,
                             method="Raw",sample_idx=sample_idx,
                             disease_celltype = "Clone_1",
                             ratio_threshold=0.02,exclude_plot_idx = exclude_plot_idx,
                             vaf_cellprop = F,coef_plot_option = "negative")
var_spot_res_grid <- result_var
# var_spot_pval: a matrix [celltype × spot] of signed p-values
var_spot_pval_grid <- sapply(
  names(var_spot_res_grid$results),
  function(bc){
    out <- var_spot_res_grid$results[[bc]]$pval
    # if you used signed_pval logic:
    coef <- var_spot_res_grid$results[[bc]]$coef
    out[coef < 0] <- 1
    out
  }
)
rownames(var_spot_pval_grid) <- names(var_spot_res_grid$results[[1]]$pval)  
celltypes_grid <- sort(unique(plot_df_grid$label))

# (c) make the two plots
bar_grid <- make_minp_barplot(
  var_spot_pval = var_spot_pval_grid,
  celltypes     = celltypes_grid,
  title         = "Min p‐value - Grid"
)+  
  scale_fill_manual(
    # force Clone_1→orange, Clone_2→red, Clone_3→blue, Other→lightgrey
    values = c(
      "Clone_1" = "orange",
      "Clone_2" = "red",
      "Clone_3" = "blue",
      "Other"   = "lightgrey"
    )
  )

map_grid <- make_spatial_sigmap(
  plot_df = plot_df_grid,
  pvals    = colMeans(var_spot_pval_grid),  # or pick one variant’s p‐vals
  title    = "Spatial −log10 pmin"
)
##############
voi=var
var_spot_list=list();var_spot_list[["Var_1"]]=var_spot_pval_grid
celltypes=c("Clone_1","Clone_2","Clone_3","Other")


# ──────────────────────────────────────────────────────────────────────────────
# 2) SUMMIT (“real‐like”) synthetic data
# ──────────────────────────────────────────────────────────────────────────────
# (a) read in plot_df for summit synthetic
plot_df_summit <- readRDS("results/plot_df_withClones_A200S6.rds")

# (b) load the spot‐test results you ran for summit data:
var_spot_res_sum <- readRDS("results/sim_real/vaf_0.5_clonesize_100_coherence_0/Var_1res_spot_level.rds")
var_spot_pval_sum <- sapply(
  names(var_spot_res_sum$results),
  function(bc){
    out <- var_spot_res_sum$results[[bc]]$pval
    coef <- var_spot_res_sum$results[[bc]]$coef
    out[coef < 0] <- 1
    out
  }
)
celltypes_sum <- colnames(readRDS(
  "../../results/Sydney_Bracht/sample_level/a200_s6_high_cov/RCTD/Major/RCTD_res.rds"  
)@results$weights)

# (c) make the two plots
bar_summit <- make_minp_barplot(
  var_spot_pval = var_spot_pval_sum,
  celltypes     = celltypes_sum,
  title         = "Min p‐value in Space"
)

df_BE_FB <- data.frame(
  celltype = celltypes_sum,
  min_pval = apply(var_spot_pval_sum, 1, min, na.rm = TRUE)
) %>%
  mutate(
    neglog10 = -log10(min_pval),
    neglog10 = ifelse(neglog10==0, 1e-2, neglog10),
  )%>%
  # floor zeros at 1e-2 for log‐scale
  #mutate(min_pval = pmin(min_pval, 100)) %>%
  filter(celltype %in% c("BE", "FB")) 

bar_summit <- ggplot(df_BE_FB, aes(x = celltype, y = neglog10, fill = celltype)) +
  geom_col() +
  scale_y_log10(
    name   = "-log10 pmin",
    limits = c(1e-2, NA),
    breaks = c(1e-2, 1e-1, 1, 10, 100),
    labels = c("0", "1e-1", "1", "10", ">100"),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      BE = "red",
      FB = "lightgrey"
    )
  ) +
  geom_hline(yintercept = 1, color = "black", size = 0.5) +
  labs(
    title = "Summit Min p-value",
    x     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x      = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5),
    axis.line.y       = element_line(color = "black"),
  )


a=plot_df_summit[names(var_spot_res_sum$results),]
map_summit <- make_spatial_sigmap(
  plot_df = a,
  pvals    = colMeans(var_spot_pval_sum),
  title    = "Summit Spatial Significance"
)


# ──────────────────────────────────────────────────────────────────────────────
# 3) Arrange bar significance & plots
# ──────────────────────────────────────────────────────────────────────────────
common_legend <- theme(
  legend.position   = "bottom",
  legend.direction  = "horizontal",
  plot.title        = element_text(hjust = 0.5)
)

p1v <- bar_grid   + common_legend
p2v <- map_grid   + common_legend
p3v <- bar_summit + common_legend
p4v <- map_summit + common_legend

final_fig <- (p1v + p2v) /
  (p3v + p4v) +
  plot_layout(guides = "keep", widths = c(1,1), heights = c(1,1))
final_stacked <- p1v / p3v +
  plot_layout(guides = "keep", heights = c(1, 1))
pdf("results/fig2f_stacked.pdf", width = 4, height = 10)
print(final_stacked)
dev.off()
