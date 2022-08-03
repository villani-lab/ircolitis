#!/usr/bin/env Rscript

library(Matrix)
library(RColorBrewer)
# library(cowplot)
library(data.table)
# library(doParallel)
# library(foreach)
library(ggbeeswarm)
library(naturalsort)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(ggstance)
library(glue)
library(harmony)
library(janitor)
library(magrittr)
library(pals)
library(patchwork)
library(pbapply)
# library(princurve)
library(readxl)
library(rhdf5)
library(scattermore)
library(scales)
library(scico)
library(shadowtext)
library(sitools)
library(tidygraph)
library(tidyverse)
library(uwot)
# library(waffle)
# Bioconductor
# library(monocle)
# library(SingleR)
library(fgsea)
library(qs)
# My functions

source("R/functions/helpers.R")
source("R/functions/read-h5-files.R")
source("R/functions/write-mtx.R")
source("R/functions/do-analysis.R")
source("R/plot-analysis.R")
source("R/functions/write-de.R")
source("R/colors-Luoma2020.R")
source("R/functions/theme-kamil.R")
theme_set(theme_kamil)

conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Global parameters
########################################################################

# Quality control parameters
# The minimum number of genes detected in a cell
min_genes <- 500
# The maximum percent of reads from a cell that are assigned to MT genes
max_mito <- 30
min_mito <- 0.1
# The minimum percent of cells in which a gene is detected
# min_pct <- 0.5
# The minimum number of cells detected in which a gene is detected
# min_cells <- 100

# Define important gene sets
########################################################################
features <- data.table::fread(
  Sys.glob("data/Luoma2020/geo/*/features.tsv.gz")[1],
  header = FALSE
)
colnames(features) <- c("ensembl_id", "symbol", "type")
features <- unique(features[,c("ensembl_id","symbol")])
ensembl_to_symbol <- unlist(split(features$symbol, features$ensembl_id))
#
# MT
mito_genes <- names(ensembl_to_symbol)[
  stringr::str_detect(ensembl_to_symbol, "^MT-")
]
length(mito_genes)
# TCR
tcr_genes <- janitor::clean_names(data.table::fread(
  "data/genenames.org/tcr.tsv"
))
ix <- str_detect(tcr_genes$approved_name, "(joining|variable|diversity)")
# tcr_genes$approved_symbol[ix]
tcr_genes <- tcr_genes$ensembl_gene_id[ix] %>% discard(function(x) x == "")
length(tcr_genes)
substr(sort(ensembl_to_symbol[tcr_genes]), 1, 4) %>% unique
x <- ensembl_to_symbol[rownames(counts)]
tcr_genes2 <- x[str_detect(x, "^TR(A|B|D|G)(V|D|J)")] %>% names
length(intersect(tcr_genes, tcr_genes2)) / length(tcr_genes2)
tcr_genes <- union(tcr_genes, tcr_genes2)
tcr_genes <- intersect(tcr_genes, rownames(counts))
rm(tcr_genes2)
# BCR
bcr_genes <- janitor::clean_names(data.table::fread(
  "data/genenames.org/immunoglobulins.tsv",
))
ix <- str_detect(bcr_genes$approved_name, "(joining|variable|diversity)")
# bcr_genes$approved_symbol[ix]
bcr_genes <- bcr_genes$ensembl_gene_id[ix] %>% discard(function(x) x == "")
length(bcr_genes)
# substr(sort(ensembl_to_symbol[bcr_genes]), 1, 4) %>% unique
x <- ensembl_to_symbol[rownames(counts)]
bcr_genes2 <- x[str_detect(x, "^IG(H|K|L)(V|D|J)")] %>% names
length(intersect(bcr_genes, bcr_genes2)) / length(bcr_genes2)
bcr_genes <- union(bcr_genes, bcr_genes2)
bcr_genes <- intersect(bcr_genes, rownames(counts))
rm(bcr_genes2)
print_status(glue("{length(mito_genes)} MT genes"))
print_status(glue("{length(tcr_genes)} TCR genes"))
print_status(glue("{length(bcr_genes)} BCR genes"))



cache_luoma <- "cache/Luoma2020-meta-counts.qs"
if (file.exists(cache_luoma)) {
  x <- qread(cache_luoma)
  meta <- x$meta
  counts <- x$counts
  rm(x)
} else {
  # Metadata
  ########################################################################

  gse_id <- "GSE144469"
  gsets <- GEOquery::getGEO(
    GEO       = gse_id,
    destdir   = "data/Luoma2020",
    GSEMatrix = TRUE,
    getGPL    = TRUE
  )

  # Chemistry
  # scRNA-seq (10X Genomics 5' V1 GEX)

  # Get a dataframe describing the columns of the expression matrix.
  x <- Biobase::pData(gsets[[2]])
  sample_info <- Biobase::pData(gsets[[3]])
  my_cols <- intersect(colnames(x), colnames(sample_info))
  sample_info <- sample_info[,my_cols]
  x <- x[,my_cols]
  sample_info <- rbind(sample_info, x)
  ix <- which(
    sapply(
      my_cols,
      function(x) {
        length(unique(sample_info[[x]]))
      }
    ) > 1
  )
  sample_info <- sample_info[,ix]
  sample_info <- janitor::clean_names(sample_info)
  sample_info <- sample_info %>%
    dplyr::select(!starts_with("data_processing")) %>%
    dplyr::select(!starts_with("supplementary_file")) %>%
    dplyr::select(!starts_with("characteristics_ch")) %>%
    dplyr::select(!starts_with("x10x")) %>%
    dplyr::select(!starts_with("relation"))
  sample_info$treatment_protocol_ch1 <- NULL
  sample_info$extract_protocol_ch1_1 <- NULL
  sample_info$donor <- str_match(sample_info$title, "(\\S+) sample")[,2]
  sample_info$type <- str_split_fixed(sample_info$title, "sample ", 2)[,2]
  sample_info$channel <- sprintf("%s-%s",
    sample_info$donor,
    str_split_fixed(sample_info$type, " ", 2)[,1]
  )
  colnames(sample_info) <- str_replace(colnames(sample_info), "_ch1", "")
  #
  sample_info$case <- ifelse(
    str_detect(sample_info$colon_status, "CPI colitis"),
    "Case",
    "Control"
  )
  #
  class_key <- c(
    "+CPI colitis"    = "CPI colitis",
    "+CPI no colitis" = "CPI no colitis",
    "Normal"          = "Normal"
  )
  sample_info <- sample_info %>%
    mutate(class = recode(colon_status, !!!class_key))
  class_short_key <- c(
    "+CPI colitis"    = "C",
    "+CPI no colitis" = "NC",
    "Normal"          = "N"
  )
  sample_info <- sample_info %>%
    mutate(class_short = recode(colon_status, !!!class_short_key))
  sample_info$class2 <- sample_info$tissue_region
  sample_info$facs_sorting <- str_split_fixed(sample_info$channel, "-", 2)[,2]
  # sample_info$case <- sample_info$class_short

  print_status(glue("{length(unique(sample_info$donor))} donors"))


  # Expression data
  ########################################################################
  print_status("Reading expression data")
  expression_data_file <- "cache/Luoma2020_expression_data.rda"
  if (!file.exists(expression_data_file)) {
    mtx_files <- Sys.glob("data/Luoma2020/geo/*/matrix.mtx.gz")
    #
    mats <- pbapply::pblapply(mtx_files, function(mtx_file) {
      sample <- mtx_file %>% dirname %>% basename
      mat <- Matrix::readMM(mtx_file)
      mat_rownames <- data.table::fread(
        glue("data/Luoma2020/geo/{sample}/features.tsv.gz"),
        header = FALSE
      )
      # rownames(mat) <- sprintf("%s|%s", mat_rownames$V1, mat_rownames$V2)
      rownames(mat) <- mat_rownames$V1
      colnames(mat) <- sprintf("%s|%s", sample, data.table::fread(
        glue("data/Luoma2020/geo/{sample}/barcodes.tsv.gz"),
        header = FALSE
      )$V1)
      mat
    })
    stopifnot(length(unique(sapply(mats, nrow))) == 1)
    counts <- do.call(cbind, mats)
    #
    genes <- data.table::fread(
      Sys.glob("data/Luoma2020/geo/*/features.tsv.gz")[1],
      header = FALSE
    )[,1:2]
    colnames(genes) <- c("ensembl_id", "symbol")
    #
    ensembl_to_symbol <- unlist(with(genes, split(symbol, ensembl_id)))
    #
    save(
      list = c("counts", "genes", "ensembl_to_symbol"),
      file = expression_data_file
    )
  } else {
    load(expression_data_file)
  }
  ensembl_to_symbol[1:5]

  sparsity <- 1 - length(counts@x) / counts@Dim[1] / counts@Dim[2]
  n_genes <- nrow(counts)
  n_cells <- ncol(counts)
  print_status(glue(
    "Unfiltered counts matrix has {comma(n_genes)} genes and {comma(n_cells)} cells"
  ))
  print_status(glue(
    "Unfiltered counts matrix is {signif(100 * sparsity, 3)}% sparse"
  ))


  # Cell information
  ########################################################################
  print_status("Cell information")

  meta <- data.table(
    cell = colnames(counts)
  )
  meta$channel  <- str_split_fixed(colnames(counts), "\\|", 2)[,1]
  # meta$donor <- str_split_fixed(meta$channel, "-", 2)[,1]
  # meta$sort <- str_split_fixed(meta$channel, "-", 2)[,2]
  meta <- dplyr::left_join(
    x  = meta,
    y  = sample_info %>%
      dplyr::filter(library_strategy == "RNA-Seq") %>%
      dplyr::select(
        channel, case, class, class2, class_short, donor,
        facs_sorting
      ),
    by = "channel"
  )
  stopifnot(nrow(meta) == ncol(counts))
  stopifnot(
    setdiff(
      unique(sample_info$channel), unique(meta$channel)
    ) == character(0)
  )
  meta$n_counts    <- Matrix::colSums(counts)
  meta$n_features  <- Matrix::colSums(counts > 0)
  meta$mito_counts <- colSums(counts[mito_genes,])
  meta$mito_pct    <- 100 * meta$mito_counts / meta$n_counts
  meta <- meta %>%
    mutate(pass_qc = mito_pct < max_mito & n_features > min_genes)
  meta$chemistry <- "5p"
  meta[1:5,]
  print_status(
    sprintf(
      "%s (%s%%) of cells pass QC with n_features > %s and mito_pct < %s%%",
      comma(sum(meta$pass_qc)), signif(100 * sum(meta$pass_qc) / ncol(counts), 3),
      min_genes, max_mito
    )
  )

  n_between <- sum(meta$mito_pct > min_mito & meta$mito_pct < max_mito)
  p <- ggplot(meta) +
    aes(n_features, mito_pct) +
    stat_binhex(bins = 101) +
    geom_hline(yintercept = c(max_mito, min_mito), size = 0.3, linetype = 2) +
    # geom_vline(xintercept = 500, size = 0.3) +
    scale_x_log10(labels = label_number_si()) +
    scale_y_sqrt(breaks = c(1, 10, 20, 40, 80, 100)) +
    annotate(
      geom = "text",
      x = 1, y = max_mito - 5,
      label = glue("{signif(100 * n_between / nrow(meta), 3)}%"),
      vjust = 1,
      hjust = 0,
      size = 5
    ) +
    # scale_y_continuous(
    #   breaks = pretty_breaks(5),
    #   expand = expansion(mult = c(0.1, 0.1))
    # ) +
    # scale_fill_scico(breaks = log_breaks(5), trans = "log10") +
    scale_fill_gradientn(
      colors = scico::scico(20)[3:20],
      breaks = log_breaks(5),
      trans = "log10"
    ) +
    annotation_logticks(sides = "b", base = 10, color = "grey20", size = 0.4) +
    labs(
      title = "Genes and mito. reads for each cell",
      x = "Number of genes", y = "Percent mito. reads"
    ) +
    guides(fill = guide_colorbar(title = "Cells", barheight = 7))
  my_ggsave(
    "raw-genes-vs-mito",
    out_dir = "results/Luoma2020",
    type = "pdf",
    plot = p,
    scale = 0.33, width = 16, height = 9, units = "in", dpi = 300
  )

  p <- ggplot(meta) +
    aes(n_features, mito_pct) +
    stat_binhex(bins = 101) +
    # geom_hline(yintercept = c(max_mito, min_mito), size = 0.3, linetype = 2) +
    # geom_vline(xintercept = 500, size = 0.3) +
    scale_x_log10(labels = label_number_si()) +
    scale_y_sqrt(breaks = c(1, 10, 20, 40, 80, 100)) +
    scale_fill_gradientn(
      colors = scico::scico(20)[3:20],
      breaks = log_breaks(5),
      trans = "log10"
    ) +
    annotation_logticks(sides = "b", base = 10, color = "grey20", size = 0.4) +
    labs(
      title = "Genes and mito. reads for each cell",
      x = "Number of genes", y = "Percent mito. reads"
    ) +
    guides(fill = guide_colorbar(title = "Cells", barheight = 7))
  my_ggsave(
    "raw-genes-vs-mito-by-channel",
    out_dir = "results/Luoma2020",
    type = "pdf",
    plot = p + facet_wrap(~ channel),
    scale = 1.67, width = 16, height = 9, units = "in", dpi = 300
  )


  # Add TCR information
  ########################################################################

  # tcr_files <- Sys.glob(
  #   "data/Luoma2020/geo/*contig_annotations.csv.gz"
  # )
  # length(tcr_files)
  # source("R/load-tcr.R")
  # tcr <- load_tcr(tcr_files)

  tcr_file <- "data/Luoma2020/GSE144469_TCR_filtered_contig_annotations_all.csv.gz"
  source("R/functions/load-tcr.R")
  tcr <- load_tcr(tcr_file)
  tcr$channel <- sprintf("%s-CD3", str_split_fixed(tcr$cell, "-", 2)[,2])
  tcr$cell <- sprintf(
    "%s|%s-1",
    tcr$channel,
    str_split_fixed(tcr$barcode, "-", 2)[,1]
  )
  tcr$donor <- str_split_fixed(tcr$channel, "-", 2)[,1]

  # Join
  meta <- dplyr::left_join(
    x = meta,
    y = tcr %>% select(-channel, -donor),
    by = "cell"
  )


  # Exclude cells with zero or too much MT
  ########################################################################

  ix_mito_pct <- with(meta, mito_pct > min_mito & mito_pct < max_mito)
  meta   <- meta[ix_mito_pct,]
  counts <- counts[,ix_mito_pct]
  counts <- counts[Matrix::rowSums(counts) > 0,]


  # Exclude poor quality cells before further analysis
  ########################################################################

  # ix_cells <- meta$is_cell == "cell" & meta$n_features > min_genes
  ix_cells <- with(meta, n_features > min_genes)
  sum(ix_cells)
  counts <- counts[,ix_cells]
  meta   <- meta[ix_cells,]
  counts <- counts[Matrix::rowSums(counts) > 0,]


  # Exclude poor channels
  ########################################################################

  x <- (
    meta %>%
    # filter(! cell %in% unwanted_cells) %>%
    group_by(channel) %>%
    summarize(genes = median(n_features)) %>%
    arrange(genes)
  )$channel
  meta$channel2 <- factor(meta$channel, x)
  p <- ggplot(meta) +
    aes(x = n_features, fill = class_short) +
    geom_histogram(bins = 50)  +
    scale_x_log10(name = "Detected Genes", labels = label_number_si()) +
    scale_fill_manual(name = "Class", values = my_colors$class_short) +
    scale_y_continuous(expand = expansion(0.15)) +
    annotation_logticks(sides = "b") +
    facet_wrap(~ channel2, ncol = 3)
  my_ggsave(
    "histogram-channel-cells",
    out_dir = "results/Luoma2020",
    plot = p,
    type = "pdf",
    scale = 1.0, width = 16.5, height = 30, units = "in", dpi = 300,
    limitsize = FALSE
  )

  x <- meta %>%
    dplyr::group_by(channel, case) %>%
    dplyr::summarize(n_cells = n())
  p <- ggplot(x) +
    aes(x = n_cells, y = reorder(channel, n_cells), fill = case) +
    geom_colh() +
    scale_fill_manual(values = pals::okabe()[2:3]) +
    scale_x_continuous(trans = "log10") +
    # scale_y_discrete(expand = c(0.035, 0)) +
    annotation_logticks(side = "b") +
    # scale_x_continuous(labels = scales::comma) +
    labs(
      y = NULL, x = "Cells", fill = NULL,
      title = glue("Cells with >{min_genes} genes for {nrow(x)} channels")
    ) +
    theme(panel.grid.major.x = element_line(size = 1))
  fig_height <- 0.2 * nrow(x)
  my_ggsave(
    "bars-channel",
    out_dir = "results/Luoma2020",
    plot = p,
    type = "pdf",
    scale = 1.8, width = 6.5, height = fig_height, units = "in", dpi = 300
  )

  x <- meta %>%
    dplyr::group_by(donor, case) %>%
    dplyr::summarize(n_cells = n())
  p <- ggplot(x) +
    aes(x = n_cells, y = reorder(donor, n_cells), fill = case) +
    geom_colh() +
    scale_fill_manual(values = pals::okabe()[2:3]) +
    scale_x_continuous(trans = "log10") +
    # scale_y_discrete(expand = c(0.035, 0)) +
    annotation_logticks(side = "b") +
    # scale_x_continuous(labels = scales::comma) +
    labs(
      y = NULL, x = "Cells", fill = NULL,
      title = glue("Cells with >{min_genes} genes for {nrow(x)} donors")
    ) +
    theme(panel.grid.major.x = element_line(size = 1))
  fig_height <- 0.2 * nrow(x)
  my_ggsave(
    "bars-donor",
    out_dir = "results/Luoma2020",
    plot = p,
    type = "pdf",
    scale = 1.8, width = 6.5, height = fig_height, units = "in", dpi = 300
  )

  x <- meta %>%
    # filter(! cell %in% unwanted_cells) %>%
    group_by(channel, case, class_short) %>%
    summarize(
      cells = length(cell),
      genes = median(n_features)
    )
  p <- ggplot(x) +
    aes(x = cells, y = genes, fill = class_short) +
    geom_vline(xintercept = 100, linetype = 2) +
    geom_point(data = x, shape = 21, size = 5, stroke = 0.2) +
    # geom_hline(yintercept = 400, linetype = 2) +
    # geom_text_repel(
    #   data = subset(x, cells < 500 | genes < 500),
    #   mapping = aes(label = channel)
    # ) +
    scale_x_log10(expand = expansion(0.1)) +
    # scale_y_log10(expand = expansion(0.1), breaks = log_breaks(8)) +
    scale_fill_manual(values = my_colors$class_short) +
    annotation_logticks(sides = "b") +
    labs(
      title = glue("Cells and median genes for {nrow(x)} channels"),
      x = "Cells",
      y = "Median genes per cell",
      fill = NULL
    )
  my_ggsave(
    "scatter-channel-cells-vs-median-genes",
    out_dir = "results/Luoma2020",
    plot = p,
    type = "pdf",
    scale = 1.8, width = 3.5, height = 2.5, units = "in", dpi = 300
  )
  # Require a median of 400 genes per cell
  include_channels <- x$channel[
    x$cells > 100
  ]

  ix_channel <- meta$channel %in% include_channels
  meta <- meta[ix_channel,]
  counts <- counts[,ix_channel]
  stopifnot(all(meta$cell == colnames(counts)))

  qsave(list(meta = meta, counts = counts), cache_luoma)
}


# TCRdiv and alpha diversity
########################################################################

tcr_donors <- sort(unique(
  meta$donor[meta$channel %in% include_channels & meta$class_short != "CDMC"]
))
#
d <- tcr %>%
  dplyr::filter(donor %in% tcr_donors) %>%
  dplyr::group_by(donor) %>%
  dplyr::count(
    TRBV,
    TRAV,
    TRB_cdr3_trim,
    TRA_cdr3_trim
  ) %>%
  dplyr::arrange(-n) %>%
  dplyr::mutate(rank = seq(length(n)), percent = 100 * n / sum(n))
d <- inner_join(d, meta %>% select(donor, case) %>% unique, by = "donor")
#
p1 <- ggplot(d) +
  aes(x = rank, y = percent, group = donor, color = case) +
  geom_line(size = 1) +
  scale_color_manual(values = my_colors$case) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(sides = "bl") +
  labs(x = "Rank", y = "Percent")
# p1
#
x <- tcr %>%
  filter(donor %in% tcr_donors) %>%
  mutate(
    SAMPLE = donor,
    CLONE = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim)
  )
# x <- tcr %>% filter(donor %in% unique(obs$donor)) %>%
#   mutate(SAMPLE = donor, CLONE = TRB_cdr3_trim)
sample_curve_file <- "cache/Luoma2020_tcr_sample_curve.rds"
if (!file.exists(sample_curve_file)) {
  sample_curve <- alakazam::alphaDiversity(
    x,
    group = "SAMPLE",
    clone = "CLONE",
    min_q = 0,
    max_q = 4,
    step_q = 0.1,
    ci = 0.95,
    nboot = 200
  )
  saveRDS(sample_curve, sample_curve_file)
} else {
  sample_curve <- readRDS(sample_curve_file)
}
#
d_div <- sample_curve@diversity
d_div <- left_join(d_div, meta %>% select(donor, case) %>% unique, by = c("SAMPLE" = "donor"))
#
p2 <- ggplot(d_div) +
  aes(Q, D, color = case, group = SAMPLE) +
  geom_line(size = 1) +
  scale_color_manual(values = my_colors$case, name = NULL)
p <- (
    p1 + theme(legend.position = "none") +
      labs(title = "Rank abundance")
  ) + (
    p2 + labs(title = "Sample diversity")
  )
# p
x <- left_join(x, meta %>% select(donor, case) %>% unique, by = c("donor"))
d_div_counts <- x %>%
  group_by(donor, case) %>%
  summarise(n = length(CLONE), n_unique = length(unique(CLONE)))
p3 <- ggplot(d_div_counts) +
  geom_colh(
    mapping = aes(x = n, y = reorder(donor, n_unique), fill = case),
    color = "white", size = 0.5, alpha = 0.7
  ) +
  geom_colh(
    mapping = aes(x = n_unique, y = reorder(donor, n_unique), fill = case),
    color = "white", size = 0.5
  ) +
  scale_x_continuous(labels = comma) +
  labs(title = "Number of TCRs", x = "TCRs (Unique | Total)", y = NULL) +
  scale_fill_manual(values = my_colors$case, name = NULL) +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank()
  )
p <- (
    p1 + theme(legend.position = "none") +
      labs(title = "Rank abundance")
  ) + (
    p2 + labs(title = "Sample diversity") +
      theme(legend.position = "none")
  ) + p3
# p
my_ggsave(
  "tcr-diversity",
  out_dir = "results/Luoma2020",
  plot = p,
  type = "pdf",
  scale = 1, width = 14, height = 6, units = "in", dpi = 300
)


# Analysis for each channel separately
########################################################################
if (FALSE) {

my_channels <- sort(unique(meta$channel))
#
for (my_channel in my_channels) {
  keep_cells <- which(meta$channel == my_channel)
  analysis_name <- my_channel
  a1_file <- as.character(glue("results/Luoma2020/per-channel/{analysis_name}/{analysis_name}.rds"))
  file.exists(a1_file)
  dir.create(dirname(a1_file), showWarnings = FALSE, recursive = TRUE)
  try({
    if (!file.exists(a1_file)) {
      source("R/do-analysis.R")
      a1 <- do_analysis(
        obs           = meta[keep_cells,],
        counts        = counts[,keep_cells],
        exclude_genes = union(tcr_genes, bcr_genes),
        mito_genes    = mito_genes,
        min_percent   = 100 * (50 / ncol(counts[,keep_cells])),
        loess_span    = 0.01,
        n_pcs         = 20,
        n_harmony     = 0,
        harmony_vars  = c("channel"),
        n_knn         = 30,
        leiden_res    = 1.3,
        leiden_iter   = 20,
        umap_spread   = 1,
        umap_min_dist = 0.05
      )
      print_status(glue("Writing {a1_file}"))
      saveRDS(a1, a1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {a1_file}"))
      a1 <- readRDS(a1_file)
      print_status(glue("done"))
    }
    #
    a1_cells_file <- as.character(glue("results/Luoma2020/per-channel/{analysis_name}/cells.tsv.gz"))
    write_tsv(a1$obs, a1_cells_file)
    #
    a1$de$symbol <- ensembl_to_symbol[a1$de$feature]
    a1$de %<>% select(feature, symbol, everything())
    a1_de_file <- as.character(glue("results/Luoma2020/per-channel/{analysis_name}/de.xlsx"))
    write_de(a1$de, a1_de_file, n = 1000)
    #
    assign(analysis_name, a1)
    source("R/plot-analysis.R")
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = glue("results/Luoma2020/per-channel/{analysis_name}/figures"),
      rowname_key   = ensembl_to_symbol,
      exclude_genes = union(tcr_genes, bcr_genes)
    )
    #
    source("R/cellbrowser-colitis.R")
    make_cellbrowser(
      analysis_name = analysis_name,
      out_dir = glue("results/Luoma2020/per-channel/{analysis_name}/figures"),
      ensembl_to_symbol
    )
    #
    rm(list = c("a1", analysis_name))
  })
}

}

# analysis
########################################################################

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6546699/
# Meg: ITGA2B, PF4, VWF
# E: CA1, HBB, KLF1, TFR2
# DC: CCR2, IRF8, MPEG1
# G: ELANE, MPO, LYZ, CSF1R, CTSG, PRTN3, AZU1
# Ly1: RGS1, NPTX2, DDIT4, ID2
# Ly2: DNTT, RAG1, RAG2 HSC: CRHBP, HLF, DUSP1
#
# And for the Lin−CD34/CD164 data set:
# E: KLF1, CA1
# Meg: ITGA2B, PLEK
# BEM: CLC, CPA3, HDC
# Ly: DNTT, CD79A, VPREB1
# DC: IRF8, SPIB, IGKC
# M: LYZ, MS4A6A, ANXA2
# N: ELANE
# HSC: HLF, ADGRG6, CRHBP, PCDH9

min_genes <- 500
n_pcs <- 20
keep_cells <- with(meta, which(class_short != "CDMC"))
length(keep_cells)
# length(keep_cells) / nrow(meta)
source("R/do-analysis.R")
analysis_name <- as.character(glue("a1_min_genes{min_genes}_n_pcs{n_pcs}"))
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/data/{analysis_name}.rds"))
file.exists(a1_file)
dir.create(dirname(a1_file), showWarnings = FALSE, recursive = TRUE)
if (!file.exists(a1_file)) {
  assign(analysis_name, do_analysis(
    obs           = meta[keep_cells,],
    counts        = counts[,keep_cells],
    exclude_genes = union(tcr_genes, bcr_genes),
    mito_genes    = mito_genes,
    min_percent   = 100 * (300 / ncol(counts[,keep_cells])),
    loess_span    = 0.01,
    n_pcs         = n_pcs,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 30,
    leiden_res    = 1.3,
    leiden_iter   = 20,
    umap_spread   = 1,
    umap_min_dist = 0.05
  ))
  a1 <- get(analysis_name)
  print_status(glue("Writing {a1_file}"))
  saveRDS(a1, a1_file)
  print_status(glue("done"))
# }
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- readRDS(a1_file)
  a1$obs$case <- "Case"
  a1$obs$case[str_detect(a1$obs$channel, "^CT")] <- "Control"
  a1$obs$case[str_detect(a1$obs$channel, "^NC")] <- "Control"
  print_status(glue("done"))
}
a1$de$symbol <- ensembl_to_symbol[a1$de$feature]
a1$de %<>% select(feature, symbol, everything())
assign(analysis_name, a1)
source("R/plot-analysis.R")
plot_analysis(
  analysis_name = analysis_name,
  out_dir = glue("results/Luoma2020/{analysis_name}/figures"),
  rowname_key = ensembl_to_symbol,
  exclude_genes = union(tcr_genes, bcr_genes)
)
source("R/cellbrowser-colitis.R")
make_cellbrowser(
  analysis_name = analysis_name,
  out_dir = glue("results/Luoma2020/{analysis_name}/figures"),
  cb_dir = glue("cellbrowser/Luoma2020/{analysis_name}"),
  ensembl_to_symbol
)
# a1_exclude_cells <- a1$obs$cell[
#   a1$obs$leiden %in% c("21", "20", "25", "22", "18", "11", "23", "24")
# ]
rm(list = c("a1", analysis_name))



# Only the CD45 channels
########################################################################
min_genes <- 500
n_pcs <- 20
keep_cells <- which(str_detect(meta$channel, "CD45"))
length(keep_cells)
# length(keep_cells) / nrow(meta)
source("R/do-analysis.R")
analysis_name <- as.character(glue("a1_CD45_min_genes{min_genes}_n_pcs{n_pcs}"))
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/data/{analysis_name}.rds"))
file.exists(a1_file)
dir.create(dirname(a1_file), showWarnings = FALSE, recursive = TRUE)
if (!file.exists(a1_file)) {
  assign(analysis_name, do_analysis(
    obs           = meta[keep_cells,],
    counts        = counts[,keep_cells],
    exclude_genes = union(tcr_genes, bcr_genes),
    mito_genes    = mito_genes,
    min_percent   = 100 * (300 / ncol(counts[,keep_cells])),
    loess_span    = 0.01,
    n_pcs         = n_pcs,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 30,
    leiden_res    = 1.3,
    leiden_iter   = 20,
    umap_spread   = 1,
    umap_min_dist = 0.05
  ))
  a1 <- get(analysis_name)
  print_status(glue("Writing {a1_file}"))
  saveRDS(a1, a1_file)
  print_status(glue("done"))
# }
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- readRDS(a1_file)
  print_status(glue("done"))
}
a1$de$symbol <- ensembl_to_symbol[a1$de$feature]
a1$de %<>% select(feature, symbol, everything())
assign(analysis_name, a1)
source("R/plot-analysis.R")
plot_analysis(
  analysis_name = analysis_name,
  out_dir = glue("results/Luoma2020/{analysis_name}/figures"),
  rowname_key = ensembl_to_symbol,
  exclude_genes = union(tcr_genes, bcr_genes)
)
source("R/cellbrowser-colitis.R")
make_cellbrowser(
  analysis_name = analysis_name,
  out_dir = glue("results/Luoma2020/{analysis_name}/figures"),
  cb_dir = glue("cellbrowser/Luoma2020/{analysis_name}"),
  ensembl_to_symbol
)
# a1_exclude_cells <- a1$obs$cell[
#   a1$obs$leiden %in% c("21", "20", "25", "22", "18", "11", "23", "24")
# ]
rm(list = c("a1", analysis_name))


# Analysis 2021-09-21
########################################################################

stopifnot(nrow(meta) == ncol(counts))
stopifnot(all(meta$cell == colnames(counts)))

analysis_name <- "luoma_a4"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
file.exists(a1_file)
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
ix_keep <- rep(TRUE, ncol(counts))
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(counts[,ix_keep]),
    n_cells       = ncol(counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 40,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  a1 <- run_analysis(
    obs           = meta[ix_keep,],
    counts        = counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})
assign(analysis_name, a1)
# rm(list = c("a1", analysis_name))
source("R/functions/make-cellbrowser.R")
make_cellbrowser(
  analysis_name     = analysis_name,
  out_dir           = file.path(dirname(a1_file), "figures"),
  cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)


analysis_name <- "luoma_cd3_a5"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
file.exists(a1_file)
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
ix_keep <- meta$facs_sorting == "CD3"
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(counts[,ix_keep]),
    n_cells       = ncol(counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 40,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  a1 <- run_analysis(
    obs           = meta[ix_keep,],
    counts        = counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})
assign(analysis_name, a1)
# rm(list = c("a1", analysis_name))
source("R/functions/make-cellbrowser.R")
make_cellbrowser(
  analysis_name     = analysis_name,
  out_dir           = file.path(dirname(a1_file), "figures"),
  cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)



########################################################################
########################################################################


# Round 0 - CD45 channels
########################################################################
analysis_name <- "luoma_cd45_a5"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
file.exists(a1_file)
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
ix_keep <- meta$facs_sorting == "CD45"
sum(ix_keep)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(counts[,ix_keep]),
    n_cells       = ncol(counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  a1 <- run_analysis(
    obs           = meta[ix_keep,],
    counts        = counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})

# assign(analysis_name, a1)
# # rm(list = c("a1", analysis_name))
# source("R/functions/make-cellbrowser.R")
# make_cellbrowser(
#   analysis_name     = analysis_name,
#   out_dir           = file.path(dirname(a1_file), "figures"),
#   cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
#   ensembl_to_symbol = ensembl_to_symbol
# )


# Round 1
########################################################################
# luoma_cd45_a5
# Keep: 3, 4, 5, 7, 8, 12, 16, 17, 18, 19, 20, 22
# Discard: 1, 2, 6, 9, 10, 11, 13, 14, 15, 21, 23
analysis_name <- "luoma_cd45_a5"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
table(a0$obs$leiden0.933)
cells_keep <- a0$obs$cell[a0$obs$leiden0.933 %in% c(3, 4, 5, 7, 8, 12, 16, 17, 18, 19, 20, 22)]
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell1"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}

source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})

# assign(analysis_name, a1)
# # rm(list = c("a1", analysis_name))
# source("R/functions/make-cellbrowser.R")
# make_cellbrowser(
#   analysis_name     = analysis_name,
#   out_dir           = file.path(dirname(a1_file), "figures"),
#   cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
#   ensembl_to_symbol = ensembl_to_symbol
# )

# Round 2
########################################################################
# luoma_cd45_a5_tcell1
# Discard: 11, 12, 16, 18, 19, 20, 21
analysis_name <- "luoma_cd45_a5_tcell1"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
table(a0$obs$leiden0.933)
cells_keep <- a0$obs$cell[!a0$obs$leiden0.933 %in% c(11, 12, 16, 18, 19, 20, 21)]
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})

# assign(analysis_name, a1)
# # rm(list = c("a1", analysis_name))
# source("R/functions/make-cellbrowser.R")
# make_cellbrowser(
#   analysis_name     = analysis_name,
#   out_dir           = file.path(dirname(a1_file), "figures"),
#   cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
#   ensembl_to_symbol = ensembl_to_symbol
# )


# Round 3 - CD4
########################################################################
# luoma_cd45_a5_tcell2_cd4_1
# Discard: 10, 11, 14, 15, 16
# Assign to CD4: 1, 4, 6, 8, 12
# Assign to CD8: 2, 3, 5, 7, 9, 13, 17
analysis_name <- "luoma_cd45_a5_tcell2"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
table(a0$obs$leiden0.933)
cells_keep <- a0$obs$cell[!a0$obs$leiden0.933 %in% c(10, 11, 14, 15, 16)]
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
ix_cd4 <- a0$obs$leiden0.933 %in% c(1, 4, 6, 8, 12)
ix_keep <- ix_keep & ix_cd4
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd4_1"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})

# assign(analysis_name, a1)
# # rm(list = c("a1", analysis_name))
# source("R/functions/make-cellbrowser.R")
# make_cellbrowser(
#   analysis_name     = analysis_name,
#   out_dir           = file.path(dirname(a1_file), "figures"),
#   cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
#   ensembl_to_symbol = ensembl_to_symbol
# )


# Round 4 - CD8
########################################################################
# luoma_cd45_a5_tcell2_cd8_1
# Discard: 10, 11, 14, 15, 16
# Assign to CD4: 1, 4, 6, 8, 12
# Assign to CD8: 2, 3, 5, 7, 9, 13, 17
analysis_name <- "luoma_cd45_a5_tcell2"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
table(a0$obs$leiden0.933)
cells_keep <- a0$obs$cell[!a0$obs$leiden0.933 %in% c(10, 11, 14, 15, 16)]
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
ix_cd8 <- a0$obs$leiden0.933 %in% c(2, 3, 5, 7, 9, 13, 17)
ix_keep <- ix_keep & ix_cd8
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd8_1"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})

# Round 4 - CD4
########################################################################
# luoma_cd45_a5_tcell2_cd4_1
# Assign to CD8: 7, 11
analysis_name <- "luoma_cd45_a5_tcell2_cd4_1"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
table(a0$obs$leiden0.933)
cells_keep <- a0$obs$cell[!a0$obs$leiden0.933 %in% c(7, 11)]
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd4_2"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})

# Round 4 - CD8
########################################################################
# luoma_cd45_a5_tcell2
analysis_name <- "luoma_cd45_a5_tcell2"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
a_cd8 <- qread("results/Luoma2020/luoma_cd45_a5_tcell2_cd8_1/luoma_cd45_a5_tcell2_cd8_1.qs")
a_cd4 <- qread("results/Luoma2020/luoma_cd45_a5_tcell2_cd4_1/luoma_cd45_a5_tcell2_cd4_1.qs")
table(a_cd8$obs$leiden0.933)
table(a_cd4$obs$leiden0.933)
cells_keep <- unique(c(
  a_cd8$obs$cell,
  a_cd4$obs$cell[a_cd4$obs$leiden0.933 %in% c(7, 11)]
))
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd8_2"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})


# Round 5 - CD4
########################################################################
# Assign to CD8: 9
analysis_name <- "luoma_cd45_a5_tcell2_cd4_2"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
table(a0$obs$leiden0.933)
cells_keep <- a0$obs$cell[!a0$obs$leiden0.933 %in% c(9)]
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd4_3"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})
assign(analysis_name, a1)
source("R/functions/make-cellbrowser.R")
make_cellbrowser(
  analysis_name     = analysis_name,
  out_dir           = file.path(dirname(a1_file), "figures"),
  cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)

# Round 5 - CD8
########################################################################
# luoma_cd45_a5_tcell2
analysis_name <- "luoma_cd45_a5_tcell2"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
a_cd8 <- qread("results/Luoma2020/luoma_cd45_a5_tcell2_cd8_2/luoma_cd45_a5_tcell2_cd8_2.qs")
a_cd4 <- qread("results/Luoma2020/luoma_cd45_a5_tcell2_cd4_2/luoma_cd45_a5_tcell2_cd4_2.qs")
table(a_cd8$obs$leiden0.933)
table(a_cd4$obs$leiden0.933)
cells_keep <- unique(c(
  a_cd8$obs$cell,
  a_cd4$obs$cell[a_cd4$obs$leiden0.933 %in% c(9)]
))
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd8_3"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})
assign(analysis_name, a1)
source("R/functions/make-cellbrowser.R")
make_cellbrowser(
  analysis_name     = analysis_name,
  out_dir           = file.path(dirname(a1_file), "figures"),
  cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)

# Round 6 - CD8
########################################################################
analysis_name <- "luoma_cd45_a5_tcell2_cd8_3"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)
#
cells_keep <- unique(c(
  a0$obs$cell[!a0$obs$leiden0.933 %in% c(4, 9)]
))
ix_keep <- a0$obs$cell %in% cells_keep
sum(ix_keep)
#
analysis_name <- "luoma_cd45_a5_tcell2_cd8_4"
a1_file <- as.character(glue("results/Luoma2020/{analysis_name}/{analysis_name}.qs"))
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(a0$counts[,ix_keep]),
    n_cells       = ncol(a0$counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(a0$counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 25,
    harmony_vars  = c("channel"),
    n_knn         = 50,
    leiden_res    = seq(0.5, 1.8, length.out = 10),
    leiden_iter   = 10,
    umap_spread   = 1,
    umap_min_dist = 0.25,
    log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
  )
  a1_params_file <- glue("{dirname(a1_file)}/params.json")
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
  #
  keep_cols <- colnames(a0$obs)
  keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|UMAP|leiden|cluster)")]
  #
  a1 <- run_analysis(
    obs           = a0$obs[ix_keep, keep_cols, with = FALSE],
    counts        = a0$counts[,ix_keep],
    exclude_genes = c(mito_genes, bcr_genes),
    mito_genes    = mito_genes,
    params        = params
  )
  a1$params <- params
  print_status(glue("Writing {a1_file}"))
  qsave(a1, a1_file)
  print_status(glue("done"))
} else {
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
}
source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = file.path(dirname(a1_file), "figures"),
  rowname_key   = ensembl_to_symbol,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
  do_pb         = FALSE
)
})
assign(analysis_name, a1)
source("R/functions/make-cellbrowser.R")
make_cellbrowser(
  analysis_name     = analysis_name,
  out_dir           = file.path(dirname(a1_file), "figures"),
  cb_dir            = glue("cellbrowser/colitis/Luoma2020/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)


# Supplemental tables
########################################################################
analyses <- c(
  "luoma_cd45_a5_tcell2_cd8_3",
  "luoma_cd45_a5_tcell2_cd8_4",
  "luoma_cd45_a5_tcell2_cd4_3"
)

retval <- rbindlist(lapply(analyses, function(analysis_name) {
	retval <- fread(as.character(glue(
    "results/a20/{analysis_name}/figures/de-case-vs-control/de_summary_case-vs-control.tsv"
  )))
  # retval$analysis <- analysis_name_pub[analysis_name]
  retval$analysis <- analysis_name
  retval
})) %>% relocate(analysis)
fwrite(retval, "paper/luoma_genes_case-vs-control.tsv", sep = "\t")

de_case <- rbindlist(lapply(analyses, function(analysis_name) {
  x <- fread(glue(
    "results/a20/{analysis_name}/figures/de-case-vs-control/de_donor_case-vs-control.tsv.gz"
  ))
  # x$analysis <- analysis_name_pub[analysis_name]
  x$analysis <- analysis_name
  x$cluster <- "all cells"
  #
  y <- fread(glue(
    "results/a20/{analysis_name}/figures/de-case-vs-control/de_case-vs-control.tsv.gz"
  ))
  # y$analysis <- analysis_name_pub[analysis_name]
  y$analysis <- analysis_name
  y$GeneName <- NULL
  x <- rbind(x, y)
})) %>% relocate(analysis, cluster, Gene, percent)
de_case <- clean_names(de_case)
de_case <- de_case %>% mutate_if(is.numeric, signif, 3) # saves 50% file size
fwrite(de_case, "paper/luoma_de-case-vs-control.tsv.gz", sep = "\t")


# Direct comparison
########################################################################
my_pairs <- list(
  "cd8_nocycling" = c(
    "Tissue CD8 T cells",
    "luoma_cd45_a5_tcell2_cd8_4",
    c(
      "CXCL13", "IL26", "IL17A", "CTLA4", "IFNG", "IL21", "IL10",
      "LINC02195", "SPINK2", "SPON1", "IL22", "MKI67"
    )
  ),
  "cd8_cycling" = c(
    "Tissue CD8 T cells",
    "luoma_cd45_a5_tcell2_cd8_3",
    c(
      "CXCL13", "IL26", "IL17A", "CTLA4", "IFNG", "IL21", "IL10",
      "LINC02195", "SPINK2", "SPON1", "IL22", "MKI67"
    )
  ),
  "cd4" = c(
    "Tissue CD4 T cells",
    "luoma_cd45_a5_tcell2_cd4_3",
    c(
      "CXCL13", "DUSP4", "HLA-DRA", "EBI3", "JAKMIP1", "ETV7", "LAIR2",
      "IL1R2", "HAVCR2", "IL7", "SPINK2", "MKI67", "IL17A"
    )
  )
)
my_pair <- "cd4"

d_both <- rbindlist(list(
  fread("paper/luoma_de-case-vs-control.tsv.gz"),
  fread("paper/de-case-vs-control.tsv.gz")
))

for (my_pair in names(my_pairs)) {
  message(my_pair)
  d <- d_both
  d %>% count(analysis)
  #
  x1 <- my_pairs[[my_pair]][1]
  x2 <- my_pairs[[my_pair]][2]
  # x1 <- "Tissue CD8 T cells"
  # x2 <- "luoma_cd45_a5_tcell2_cd8_4"
  d %<>% filter(analysis %in% c(x1, x2))
  d %<>% filter(cluster == "all cells")
  d %>% count(analysis)
  #
  d <- inner_join(
    x = d %>% filter(analysis == x1),
    y = d %>% filter(analysis == x2),
    by = c("ensembl_id", "gene")
  )
  head(d)
  my_genes <- tail(my_pairs[[my_pair]], -2)
  my_cor <- cor.test(d$log_fc.x, d$log_fc.y, method = "pearson")
  xrng <- range(d$log_fc.x, drop.na = TRUE)
  yrng <- range(d$log_fc.y, drop.na = TRUE)
  my_d <- d %>% filter(gene %in% my_genes)
  #
  p <- ggplot() +
    geom_scattermore(
      data = d,
      mapping = aes(x = log_fc.x, y = log_fc.y),
      pixels = c(800, 800),
      pointsize = 1,
      alpha = 0.5
    ) +
    # geom_hex(
    #   data = d,
    #   mapping = aes(x = log_fc.x, y = log_fc.y),
    #   bins = 201
    # ) +
    # scale_fill_gradientn(colors = scico::scico(11)[2:11]) +
    geom_abline(
      intercept = 0, slope = c(-1, 0, 1), linetype = 1, size = 0.2,
      alpha = 0.2
    ) +
    geom_vline(xintercept = 0, size = 0.2, alpha = 0.2) +
    annotate(
      geom = "text",
      label = glue(
        "{comma(nrow(d))} genes\nPearson's r = {signif(my_cor$estimate, 1)}"
      ),
      x = xrng[1],
      y = yrng[2],
      hjust = "inward", vjust = "inward"
    ) +
    geom_linerange(
      data = my_d,
      mapping = aes(xmin = ci_l.x, xmax = ci_r.x, y = log_fc.y),
      size = 0.2, alpha = 0.2
    ) +
    geom_linerange(
      data = my_d,
      mapping = aes(ymin = ci_l.y, ymax = ci_r.y, x = log_fc.x),
      size = 0.2, alpha = 0.2
    ) +
    geom_point(
      data = my_d,
      mapping = aes(x = log_fc.x, y = log_fc.y),
      size = 0.5, color = "red"
    ) +
    geom_text_repel(
      data = my_d,
      mapping = aes(x = log_fc.x, y = log_fc.y, label = gene),
      fontface = "italic"
    ) +
    scale_x_continuous(
      breaks = seq(-10, 10, by = 2),
      labels = function(x) fractional::fractional(2 ^ x)
      # expand = expansion(mult = 0.1)
    ) +
    scale_y_continuous(
      breaks = seq(-10, 10, by = 2),
      labels = function(x) fractional::fractional(2 ^ x)
      # expand = expansion(mult = 0.1)
    ) +
    labs(
      x = "Our data",
      y = "Luoma 2020",
      title = "Fold-changes for case vs control"
    )
  my_ggsave(
    glue("compare-logfc-case-vs-control-{my_pair}"),
    out_dir = "results/Luoma2020",
    type = "pdf",
    plot = p,
    scale = 1, width = 5, height = 5, units = "in", dpi = 300
  )
  #
  my_ens <- (
    d %>% mutate(score = (log_fc.x + log_fc.y) / 2) %>%
    arrange(-score) %>%
    filter(gene %in% my_genes)
    # head(20)
  )$ensembl_id
  #
  d2 <- d_both
  d2 %>% count(analysis)
  #
  x1 <- my_pairs[[my_pair]][1]
  x2 <- my_pairs[[my_pair]][2]
  # x1 <- "Tissue CD8 T cells"
  # x2 <- "luoma_cd45_a5_tcell2_cd8_4"
  d2 %<>% filter(analysis %in% c(x1, x2))
  d2 %<>% filter(cluster == "all cells")
  d2 %>% count(analysis)
  d2 %<>% filter(ensembl_id %in% my_ens)
  d2 %<>% mutate(ensembl_id = factor(ensembl_id, rev(my_ens)))
  #
  xrng <- range(c(d2$ci_l, d2$ci_r))
  p <- ggplot(d2) +
    aes(group = analysis) +
    ggforestplot::geom_stripes(
      mapping = aes(y = ensembl_id)
    ) +
    geom_vline(
      # xintercept = seq(floor(xrng[1]), ceiling(xrng[2]), by = 2),
      xintercept = seq(log2(1/16), log2(256), by = 2),
      color = "white", size = 0.5
    ) +
    geom_linerange(
      mapping = aes(xmin = ci_l, xmax = ci_r, y = ensembl_id, color = analysis),
      size = 0.4,
      position = position_dodge(width = 0.5)
    ) +
    geom_point(
      mapping = aes(x = log_fc, y = ensembl_id, color = analysis),
      size = 1.8,
      position = position_dodge(width = 0.5)
    ) +
    scale_color_manual(values = pals::okabe()[c(7,1)]) +
    scale_y_discrete(
      labels = function(x) ensembl_to_symbol[x]
    ) +
    scale_x_continuous(
      breaks = seq(-10, 10, by = 2),
      labels = function(x) fractional::fractional(2 ^ x)
      # expand = expansion(mult = 0.1)
    ) +
    # facet_row(vars(analysis)) +
    labs(
      x = NULL, y = NULL,
      # title = glue("Fold change for Case vs Control {my_pair}"),
      title = "Fold change for Case vs Control",
      subtitle = my_pair
    ) +
    theme(
      axis.text.y = element_text(face = 'italic'),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    )
  my_ggsave(
    glue("lanes-logfc-case-vs-control-{my_pair}"),
    out_dir = "results/Luoma2020",
    type = "pdf",
    plot = p,
    scale = 1, width = 5, height = 6.5, units = "in", dpi = 300
  )
}

# Annotation with SingleR
########################################################################

analysis_name0 <- "luoma_cd45_a5_tcell2_cd8_3"
a0_file <- as.character(glue("results/Luoma2020/{analysis_name0}/{analysis_name0}.qs"))
stopifnot(file.exists(a0_file))
a0 <- qread(a0_file)

analysis_name1 <- "a12_4_4_t4_cd8_1_2"
a1_file <- as.character(glue("results/a20/{analysis_name1}/data/{analysis_name1}.qs"))
stopifnot(file.exists(a1_file))
a1 <- qread(a1_file)

# pacman::p_install(SingleR)
library(SingleR)

a0$log2cpm <- do_log2cpm(a0$counts, total = median(colSums(counts)))
a1$log2cpm <- do_log2cpm(a1$counts, total = median(colSums(counts)))

if (!"pred" %in% names(a0)) {
  # 40 minutes
  system.time({
  a0$pred <- SingleR(
    test      = a0$log2cpm,
    ref       = a1$log2cpm,
    labels    = a1$obs$leiden,
    de.method = "wilcox"
  )
  table(a0$pred$labels)
  qsave(a0, a0_file)
  })
}

a0$obs$leiden_pred <- a0$pred$labels

d <- a0$obs %>% dplyr::count(leiden, leiden_pred) %>%
  mutate(
    leiden = naturalfactor(leiden),
    leiden_pred = naturalfactor(leiden_pred)
  )
d <- d %>% group_by(leiden) %>% mutate(pct_leiden = n / sum(n))
d <- d %>% group_by(leiden_pred) %>% mutate(pct_leiden_pred = n / sum(n))

p <- ggplot(d) +
  aes(x = leiden, y = leiden_pred, fill = 100 * pct_leiden) +
  geom_tile() +
  scale_fill_gradientn(
    name = "Percent of Luoma cluster",
    colors = scico::scico(n = 9, palette = "grayC")
  ) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  geom_point(
    data = d %>% group_by(leiden) %>% top_n(n = 1, wt = pct_leiden),
    size = 1, color = "white"
  ) +
  labs(
    x = "Luoma clusters",
    y = "Our clusters"
  ) +
  theme(legend.position = "bottom")
my_ggsave(
  glue("singler-{analysis_name0}-{analysis_name1}"),
  out_dir = "results/Luoma2020",
  type = "pdf",
  plot = p,
  scale = 1, width = 7, height = 5, units = "in", dpi = 300
)

cluster_colors <- mpn65
names(cluster_colors) <- seq_along(mpn65)
p1 <- plot_hexmix(
  x = a0$obs$UMAP1,
  y = a0$obs$UMAP2,
  group = a0$obs$leiden_pred,
  group_colors = cluster_colors,
  bins = 301
) +
labs(
  title = glue(
    "{length(unique(a0$obs$leiden_pred))} clusters of {comma(nrow(a0$obs))} cells from {length(unique(a0$obs$donor))} donors"
  )
) +
guides(fill = guide_legend()) +
theme(legend.position = "right")
my_ggsave(
  glue("umap-clusters-hex-leiden_pred-{analysis_name0}"),
  out_dir = "results/Luoma2020",
  plot = p1,
  type = "pdf",
  scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
)

p1 <- plot_scattermore(
  x = -a0$obs$UMAP1,
  y = -a0$obs$UMAP2,
  group = a0$obs$leiden_pred,
  group_colors = cluster_colors,
  pixels = 800,
  alpha = 0.5
)
my_ggsave(
  glue("umap-clusters-leiden_pred-{analysis_name0}"),
  out_dir = "results/Luoma2020",
  plot = p1,
  type = "pdf",
  scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
)

# Annotation with clustifyr
########################################################################

features <- data.table::fread(
  Sys.glob("data/Luoma2020/geo/*/features.tsv.gz")[1],
  header = FALSE
)
colnames(features) <- c("ensembl_id", "symbol", "type")
features <- unique(features[,c("ensembl_id","symbol")])
ensembl_to_symbol <- unlist(split(features$symbol, features$ensembl_id))


my_pairs <- list(
  c(
    analysis_name0 = "luoma_cd45_a5_tcell2_cd8_3",
    analysis_name1 = "a12_4_4_t4_cd8_1_2"
  ),
  c(
    analysis_name0 = "luoma_cd45_a5_tcell2_cd4_3",
    analysis_name1 = "a12_4_4_t4_cd4_2_2"
  )
)

for (my_pair in my_pairs) {
  analysis_name0 <- my_pair[["analysis_name0"]]
  analysis_name1 <- my_pair[["analysis_name1"]]

  # analysis_name0 <- "luoma_cd45_a5_tcell2_cd8_3"
  a0_file <- as.character(glue("results/Luoma2020/{analysis_name0}/{analysis_name0}.qs"))
  stopifnot(file.exists(a0_file))
  a0 <- qread(a0_file)

  # analysis_name1 <- "a12_4_4_t4_cd8_1_2"
  a1_file <- as.character(glue("results/a20/{analysis_name1}/data/{analysis_name1}.qs"))
  stopifnot(file.exists(a1_file))
  a1 <- qread(a1_file)

  out_dir <- glue("results/Luoma2020/{analysis_name0}/figures/classify-{analysis_name1}")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  a0$log2cpm <- do_log2cpm(a0$counts, total = median(Matrix::colSums(a0$counts)))
  a1$log2cpm <- do_log2cpm(a1$counts, total = median(Matrix::colSums(a1$counts)))

  y <- with(a1$obs, model.matrix(~ 0 + factor(leiden)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/") # means
  pb <- as(a1$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "factor\\(leiden\\)", "")
  pb <- as.matrix(pb)
  #
  pb_meta <- as_tibble(data.frame(cluster = colnames(pb)))
  stopifnot(nrow(pb_meta) == ncol(pb))
  #
  #pb_ova <- rbindlist(
  #  pblapply(sort(unique(pb_meta$cluster)), function(this_cluster) {
  #    pb_meta$x <- pb_meta$cluster == this_cluster
  #    des1 <- with(pb_meta, model.matrix(~ x))
  #    ob <- pb[rowMeans(pb) > 0.5,]
  #    fit1 <- lmFit(object = ob, design = des1)
  #    fit1 <- eBayes(fit1)
  #    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  #    res <- topTable(fit1, coef = 2, number = 1e6)
  #    res$coef <- this_cluster
  #    res$ensembl_id <- rownames(res)
  #    return(res)
  #  })
  #)
  ##
  #my_genes <- unique(
  #  (
  #    pb_ova %>% filter(adj.P.Val < 0.05, abs(logFC) > log2(1.5))
  #  )$ensembl_id
  #)
  #length(my_genes)

  a1$de <- presto::wilcoxauc(a1$log2cpm, a1$obs$leiden)

  my_genes <- unique(
    (
      a1$de %>% group_by(group) %>% top_n(n = 100, wt = auc)
    )$feature
  )

  # calculate correlation
  res <- clustifyr::clustify(
    input       = a0$log2cpm,
    metadata    = as.character(a0$obs$leiden),
    ref_mat     = as.matrix(pb),
    query_genes = my_genes
  )

  fwrite(
    clustifyr::cor_to_call(res),
    file.path(out_dir, glue("clustifyr-calls-{analysis_name0}-{analysis_name1}.tsv")),
    sep = "\t"
  )

  d <- clustifyr::cor_to_call(res) %>%
    dplyr::rename(ours = "type", luoma = "cluster")
  # d$luoma <- naturalfactor(d$luoma)
  # d$ours <- naturalfactor(d$ours)
  d$id <- seq(nrow(d))
  d <- d %>% pivot_longer(cols = c("luoma", "ours"))
  d <- d %>% mutate(value = naturalfactor(value))

  p1 <- ggplot(d) +
    aes(x = name, y = id, label = value, color = value) +
    geom_point(size = 8) +
    scale_color_manual(values = mpn65) +
    geom_text(size = 5, color = "black", nudge_x = ifelse(d$name == "luoma", -1, 1)) +
    scale_x_discrete(expand = c(1, 1)) +
    theme_void() +
    theme(legend.position = "none")
  # p2 <- ggplot(d %>% select(r, id) %>% unique) +
  #   aes(x = r, y = id) +
  #   geom_colh()
  my_ggsave(
    glue("clustifyr-table-{analysis_name0}-{analysis_name1}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p1,
    scale = 1, width = 2, height = 5, units = "in", dpi = 300
  )

  # # calculate correlation
  # res_cell <- clustifyr::clustify(
  #   input       = a0$log2cpm,
  #   metadata    = as.character(a0$obs$leiden),
  #   ref_mat     = as.matrix(pb),
  #   query_genes = my_genes,
  #   per_cell    = TRUE
  # )
  # res_cell_cor <- clustifyr::cor_to_call(res_cell)
  # res_cell_cor <- left_join(
  #   x = res_cell_cor,
  #   y = a0$obs %>% select(cell, leiden),
  #   by = c("cluster" = "cell")
  # )

  # d <- res_cell_cor %>% ungroup %>% dplyr::count(type, leiden)
  # colnames(d) <- c("luoma", "ours", "n")
  # d <- d %>% mutate(ours = as.character(ours)) %>% pivot_longer(cols = c("luoma", "ours"))

  # library(ggalluvial)
  # p <- d %>%
  #   mutate(
  #     luoma = (as.factor(sprintf("L%s", luoma))),
  #     ours = (as.factor(ours))
  #   ) %>%
  #   ggplot(aes(y = n, axis1 = luoma, axis2 = ours)) +
  #   geom_alluvium(aes(fill = ours), aes.bind="flows", width = 1/12) +
  #   geom_stratum(width = 1/4, fill = "white", color = "black") +
  #   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  #   scale_x_discrete(limits = c("Luoma", "Ours"),
  #                    expand = c(.05, .05)) +
  #   scale_fill_manual(values = mpn65) +
  #   labs(y = "Cases") +
  #   theme_minimal() +
  #   theme(legend.position = "none")
  # my_ggsave(
  #   glue("clustifyr-alluvial-{analysis_name0}-{analysis_name1}"),
  #   out_dir = out_dir,
  #   type = "pdf",
  #   plot = p,
  #   scale = 1, width = 7, height = 5, units = "in", dpi = 300
  # )


  d <- as.data.frame(summary(as(res, "dgCMatrix")))
  colnames(d) <- c("leiden", "leiden_pred", "correlation")
  d <- d %>% mutate(leiden = naturalfactor(leiden), leiden_pred = naturalfactor(leiden_pred))
  p <- ggplot(d) +
    aes(x = leiden, y = leiden_pred, fill = correlation) +
    geom_tile() +
    scale_fill_gradientn(
      # colors = scico::scico(n = 9, palette = "bilbao"),
      colors = tail(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")), 6),
      # limits = c(0, 1)
      breaks = pretty_breaks(3)
    ) +
    guides(fill = guide_colorbar(barwidth = 15)) +
    geom_point(
      data = d %>% group_by(leiden) %>% top_n(n = 1, wt = correlation),
      size = 1, color = "white"
    ) +
    labs(
      x = "Luoma clusters",
      y = "Our clusters",
      title = glue("clustifyr {analysis_name0}, ref = {analysis_name1}")
    ) +
    theme(legend.position = "bottom")
  my_ggsave(
    glue("clustifyr-{analysis_name0}-{analysis_name1}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1, width = 7, height = 5, units = "in", dpi = 300
  )

  d <- clustifyr::cor_to_call(res)
  to_clustifyr <- unlist(split(d$type, d$cluster))

  a0$obs$leiden_clustifyr <- unname(to_clustifyr[as.character(a0$obs$leiden)])
  a0$obs$leiden_clustifyr <- naturalfactor(as.character(a0$obs$leiden_clustifyr))

  fwrite(
    a0$obs %>% dplyr::count(leiden, leiden_clustifyr) %>% arrange(-n),
    file.path(out_dir, glue("clustifyr-{analysis_name0}-{analysis_name1}.tsv")),
    sep = "\t"
  )

  p1 <- plot_scattermore(
    x = a0$obs$UMAP1,
    y = a0$obs$UMAP2,
    group = a0$obs$leiden_clustifyr,
    group_colors =  cluster_colors[levels(a0$obs$leiden_clustifyr)],
    pixels = 800,
    alpha = 0.5
  )
  my_ggsave(
    glue("umap-clusters-leiden_clustifyr-{analysis_name0}"),
    out_dir = out_dir,
    plot = p1,
    type = "pdf",
    scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
  )

}

# Annotation with logistic regression
########################################################################

predict_logistic <- function(input_data, input_labels, query_data) {
  stopifnot(length(input_labels) == ncol(input_data))
  common_features <- intersect(rownames(input_data), rownames(query_data))
  stopifnot(length(common_features) > 5)
  input_data <- as.data.frame(t(input_data[common_features,]))
  input_data$label <- as.character(input_labels)
  query_data <- as.data.frame(t(query_data[common_features,]))
  model <- nnet::multinom(label ~., data = input_data)
  pred <- predict(model, query_data)
  return(list(model = model, pred = pred))
}

my_genes <- unique(
  (
    a1$de %>% group_by(group) %>% top_n(n = 5, wt = abs(auc - 0.5))
  )$feature
)
length(my_genes)

pred1 <- predict_logistic(
  input_data   = scale_data(a1$log2cpm[my_genes,]),
  input_labels = a1$obs$leiden,
  query_data   = scale_data(a0$log2cpm[my_genes,])
)

a0$obs$leiden_logistic <- pred1$pred

d <- a0$obs %>% dplyr::count(leiden, leiden_logistic) %>%
  mutate(
    leiden = naturalfactor(leiden),
    leiden_logistic = naturalfactor(leiden_logistic)
  )
d <- d %>% group_by(leiden) %>% mutate(pct_leiden = n / sum(n))
d <- d %>% group_by(leiden_logistic) %>% mutate(pct_leiden_logistic = n / sum(n))
p <- ggplot(d) +
  aes(x = leiden, y = leiden_logistic, fill = 100 * pct_leiden) +
  geom_tile() +
  scale_fill_gradientn(
    name = "Percent of Luoma cluster",
    colors = scico::scico(n = 9, palette = "grayC")
  ) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  geom_point(
    data = d %>% group_by(leiden) %>% top_n(n = 1, wt = pct_leiden),
    size = 1, color = "white"
  ) +
  labs(
    x = "Luoma clusters",
    y = "Our clusters",
    title = "Predictions from logistic regression"
  ) +
  theme(legend.position = "bottom")
my_ggsave(
  glue("logistic-{analysis_name0}-{analysis_name1}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1, width = 7, height = 5, units = "in", dpi = 300
)

a0$obs$leiden_logistic <- naturalfactor(as.character(a0$obs$leiden_logistic))
p1 <- plot_scattermore(
  x = -a0$obs$UMAP1,
  y = -a0$obs$UMAP2,
  group = a0$obs$leiden_logistic,
  group_colors =  cluster_colors[levels(a0$obs$leiden_logistic)],
  pixels = 800,
  alpha = 0.5
)
my_ggsave(
  glue("umap-clusters-leiden_logistic-{analysis_name0}"),
  out_dir = out_dir,
  plot = p1,
  type = "pdf",
  scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
)




