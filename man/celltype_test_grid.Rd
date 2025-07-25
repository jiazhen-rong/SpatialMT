% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carrier_celltype_determination.R
\name{celltype_test_grid}
\alias{celltype_test_grid}
\title{Perform Celltype Test Within Spatial Grids}
\usage{
celltype_test_grid(
  celltypes,
  voi,
  N_voi,
  vaf,
  Ws,
  spatial_coords,
  test_type = c("linear", "weighted"),
  permute_num = 1000,
  grid_size = 20,
  verbose = F,
  vaf_cellprop = T,
  power = T,
  disease_celltype = NULL,
  method = "Raw",
  sample_idx = NULL,
  min_cell_grid = 20
)
}
\arguments{
\item{celltypes}{A string vector of unique cell types.}

\item{voi}{A string vector of variants of interest.}

\item{vaf}{A matrix of variant allele frequencies (MT variant x cell).}

\item{Ws}{A matrix of cell type weights (cell x celltype).}

\item{spatial_coords}{A matrix of spatial coordinates (cell x 2).}

\item{test_type}{A string specifying test type: "linear" or "weighted".}

\item{permute_num}{Number of permutations for weighted test.}

\item{grid_size}{Number of bins per axis (default = 20 for 20x20 grid).}

\item{method}{P-value correction methods. Options from c("FDR","FWER","Raw").}

\item{sample_idx}{Spots as control for each grid's regression. Usually low diseased celltype proportions and not low coverage spots.}

\item{N}{A matrix of total coverage counts (MT variant x cell).}
}
\value{
A nested list of p-values for each grid and variant.
}
\description{
This function partitions the spatial coordinates into a 20x20 grid and performs
the linear regression test within each grid to identify significant cell types.
}
