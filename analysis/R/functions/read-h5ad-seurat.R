#!/usr/bin/env Rscript
# Kamil Slowikowski <kslowikowski@mgb.org>
# 2025-05-05
#
# This script shows an example of how to read the ircolitis .h5ad data files into Seurat.
#

library(glue)
library(stringr)
library(Matrix)
library(dplyr)
library(rhdf5)
library(Seurat)

base_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206299/suppl"
gz_file <- "GSE206299_ircolitis-tissue-cd8.h5ad.gz"
h5_file <- str_remove(gz_file, ".gz$")

# Download the .h5ad.gz file
if (!file.exists(h5_file)) {
  download.file(glue("{base_url}/{gz_file}"), destfile=gz_file)
  # Decompress the file
  system(glue("gzip -d {gz_file}"))
}

# Inspect the contents of the file
h5_table <- h5ls(h5_file)
table(h5_table$group)
##        /                 /obs    /obs/__categories                 /uns
##        4                  112                   64                    5
##  /uns/de             /uns/knn   /uns/knn/nn_method /uns/knn/simil_cells
##       11                    2                    2                    3
## /uns/pca                 /var                   /X
##        7                   10                    3


# Let's make a function to read this kind of h5ad file
read_h5ad <- function(my_h5_file) {
  # Read the obs table
  obs <- h5read(my_h5_file, "/obs")
  for (my_col in names(obs)) {
    if (my_col %in% c("_index", "__categories")) {
      next
    }
    obs[[my_col]] <- obs[[my_col]]
    if (my_col %in% names(obs[["__categories"]])) {
      ix <- obs[[my_col]] >= 0
      itoa <- obs[["__categories"]][[my_col]]
      obs[[my_col]][ix] <- itoa[1 + obs[[my_col]][ix]]
    }
  }
  obs$`__categories` <- NULL
  obs <- as.data.frame(obs)
  # Read gene names
  var <- as.data.frame(h5read(my_h5_file, "/var"))
  colnames(var)[1] <- "ensembl_gene"
  # Read the counts matrix
  h5 <- rhdf5::h5read(my_h5_file, "X", , bit64conversion = 'double')
  h5$shape <- c(nrow(obs), nrow(var))
  counts <- Matrix::sparseMatrix(
    dims   = h5$shape,
    i      = as.numeric(h5$indices),
    p      = as.numeric(h5$indptr),
    x      = as.numeric(h5$data),
    index1 = FALSE
  )
  counts <- t(counts)
  colnames(counts) <- obs$cell
  rownames(counts) <- var$ensembl_gene
  return(list(counts = counts, var = var, obs = obs))
}

ad <- read_h5ad(h5_file)
# Seurat does not like to accept `ad$counts`, so we need to make a new variable?!
my_counts <- ad$counts
my_meta <- ad$obs
rownames(my_meta) <- my_meta$cell

# Make a new Seurat object
seurat_object <- CreateSeuratObject(counts = my_counts, meta.data = my_meta)
## 28165 features across 25341 samples within 1 assay
## Active assay: RNA (28165 features, 0 variable features)

