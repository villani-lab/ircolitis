#!/usr/bin/env Rscript
# Kamil Slowikowski
# 2021-05-12
#

# libraries {{{
library(pacman)
pacman::p_load(
  jsonlite,
  Matrix,
  data.table,
  ggbeeswarm,
  ggforce,
  ggplot2,
  ggrepel,
  ggraph,
  ggridges,
  ggstance,
  glue,
  harmony,
  magrittr,
  naturalsort,
  pals,
  patchwork,
  parallel,
  pbapply,
  pheatmap,
  princurve,
  progressr,
  qs,
  readxl,
  rhdf5,
  scales,
  scattermore,
  scico,
  shadowtext,
  sitools,
  tidygraph,
  tidyverse,
  uwot
)
source("R/functions/helpers.R")
source("R/functions/read-h5-files.R")
source("R/functions/theme-kamil.R")
source("R/functions/read-sparse-csv.R")
source("R/functions/mpn65.R")
source("R/functions/do-demux.R")
source("R/functions/do-analysis.R")
source("R/functions/write-de.R")
theme_set(theme_kamil)
# }}}

# Globals {{{
source("R/mt-tcr-bcr-genes.R") # Gene sets
out_dir <- "results/blood"
dir.create(out_dir, showWarnings = FALSE)
# g_file_hashtags <- "analysis/terra/blood_2020-11-30/sample-hashtags.tsv"
g_file_hashtags <- "/projects/irae_blood/cellranger_output/colitis_hashing_info.csv"
x <- fread(g_file_hashtags)
x$sample <- str_remove(x$sample, "_\\d+$")
hashtag_to_donor <- x$sample
names(hashtag_to_donor) <- sprintf("%s_%s", x$batch, x$hash)
rm(x)
# }}}

# Quality control parameters {{{
# The minimum number of genes detected in a cell
min_genes <- 500
# The maximum percent of reads from a cell that are assigned to MT genes
max_mito <- 20
min_mito <- 0.1
# The minimum percent of cells in which a gene is detected
# min_pct <- 0.5
# The minimum number of cells detected in which a gene is detected
# min_cells <- 100
# }}}

# Sample information {{{
sample_info <- readRDS("cache/sample_info.rds")
donor_info <- sample_info[,colnames(sample_info)[sapply(colnames(sample_info), function(col) {
  n1 <- length(unique(sample_info[["donor"]]))
  n2 <- nrow(unique(sample_info[,c("donor", col)]))
  n1 == n2
})]] %>% unique
# }}}

# read RNA-seq expression data for blood single cell immune cells {{{
print_status("Reading expression data")
my_h5_files <- Sys.glob(
  glue(
    "/projects/irae_blood/cellranger_output/C*_gex/raw_feature_bc_matrix.h5"
  )
)
print_status(glue("Reading {length(my_h5_files)} h5 files"))
counts_file <- "results/blood/counts-raw.qs"
if (!file.exists(counts_file)) {
  counts <- read_h5_files(my_h5_files, min_var = 10)
  qsave(counts, counts_file)
} else {
  counts <- qread(counts_file)
}
colnames(counts) <- str_remove(colnames(counts), "-1$")
#
genes <- h5read(my_h5_files[1], "matrix/features")
genes <- tibble(ensembl_id = genes$id, symbol = genes$name)
#
ensembl_to_symbol <- unlist(with(genes, split(symbol, ensembl_id)))
# }}}

# Quality control for raw data {{{
obs <- data.table(
  cell = colnames(counts)
)
obs$batch <- str_split_fixed(colnames(counts), "\\|", 2)[,1]
obs$batch <- naturalfactor(obs$batch)
stopifnot(nrow(obs) == ncol(counts))
obs$n_counts    <- Matrix::colSums(counts)
obs$n_features  <- Matrix::colSums(counts > 0)
obs$mito_counts <- colSums(counts[mito_genes,])
obs$mito_pct    <- 100 * obs$mito_counts / obs$n_counts
obs <- obs %>%
  mutate(pass_qc = mito_pct < max_mito & mito_pct > min_mito & n_features > min_genes)
obs[1:5,]
print_status(
  sprintf(
    "%s (%s%%) of cells pass QC with n_features > %s and mito_pct < %s%%",
    comma(sum(obs$pass_qc)), signif(100 * sum(obs$pass_qc) / ncol(counts), 3),
    min_genes, max_mito
  )
)

p <- ggplot(obs) +
  aes(n_features, mito_pct) +
  stat_binhex(bins = 101) +
  geom_hline(yintercept = c(max_mito, min_mito), size = 0.3, linetype = 2) +
  geom_vline(xintercept = min_genes, size = 0.3, linetype = 2) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_continuous(breaks = c(0, max_mito, 50, 100)) +
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
  guides(fill = guide_colorbar(title = "Cells", barheight = 7)) +
  facet_wrap(~ batch)
my_ggsave(
  "genes-vs-mito-raw",
  out_dir = glue("{out_dir}/figures"),
  type = "pdf",
  plot = p,
  scale = 1, width = 16, height = 9, units = "in", dpi = 300
)


# Discard cells that don't pass QC
########################################################################

stopifnot(all(colnames(counts) == obs$cell))
counts <- counts[,obs$pass_qc]
obs <- obs[obs$pass_qc,]

obs$cell <- str_remove(obs$cell, "_gex")
obs$batch <- str_remove(obs$batch, "_gex")
colnames(counts) <- str_remove(colnames(counts), "_gex")
stopifnot(all(colnames(counts) == obs$cell))

p <- ggplot(obs) +
  aes(n_features, mito_pct) +
  stat_binhex(bins = 101) +
  geom_hline(yintercept = c(max_mito, min_mito), size = 0.3, linetype = 2) +
  geom_vline(xintercept = min_genes, size = 0.3, linetype = 2) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_continuous(breaks = c(0, max_mito, 50, 100)) +
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
  guides(fill = guide_colorbar(title = "Cells", barheight = 7)) +
  facet_wrap(~ batch)
my_ggsave(
  "genes-vs-mito-filtered",
  out_dir = glue("{out_dir}/figures"),
  type = "pdf",
  plot = p,
  scale = 1, width = 16, height = 9, units = "in", dpi = 300
)

# souporcell
########################################################################

normalize_soup <- function(x) {
  x <- as.data.table(x)
  my_cols <- c("log_prob_singleton", "log_prob_doublet")
  mat <- as.matrix(x[,..my_cols])
  mat <- 1 / (mat / rowSums(mat))
  mat <- mat / rowSums(mat)
  x[["prob_singleton"]] <- mat[,1]
  x[["prob_doublet"]] <- mat[,2]
  #
  my_cols <- colnames(sc)[stringr::str_detect(colnames(sc), "^cluster")]
  mat <- as.matrix(x[, ..my_cols])# %>% exp
  mat <- 1 / ( mat / rowSums(mat) )
  mat <- mat / rowSums(mat)
  for (i in 0:7) {
    x[[sprintf("prob%s", i)]] <- mat[,i+1]
  }
  return(x)
}

# 1st round of souporcell on raw cells from cellranger output
# sc_files <- Sys.glob("analysis/terra/blood_souporcell_2021-08-12/demux_output/C*/clusters.tsv")
sc_files <- Sys.glob("/projects/irae_blood/demux_output/C*/clusters.tsv")
read_sc <- function(sc_file) {
  x <- fread(sc_file)
  x$batch <- str_remove(basename(dirname(sc_file)), "_gex")
  x
}
sc <- rbindlist(lapply(sc_files, read_sc), fill = TRUE)
sc <- normalize_soup(sc)
sc$cell <- sprintf("%s|%s", sc$batch, sc$barcode)
sc$cell <- stringr::str_remove(sc$cell, "-1$")

qsave(sc, "results/blood/souporcell.qs")

intersect(obs$cell, sc$cell) %>% length

sc_count <- sc[,.(count = .N), by = .(batch, status, assignment)]
sc_count$batch <- naturalfactor(sc_count$batch)

setdiff(sc_count$batch, obs$batch)

sc_x <- sc[
  status == "singlet" &
  cell %in% obs$cell
][
  ,.(count = .N), by = .(batch, status, assignment)
]
sc_x$batch <- naturalfactor(sc_x$batch)

p <- ggplot(sc_x) +
  aes(y = assignment, x = count, label = count) +
  geom_colh() +
  geom_text(
    hjust = -0.1,
    size = 5
  ) +
  facet_wrap(~ batch) +
  scale_x_continuous(labels = label_number_si()) +
  expand_limits(x = max(sc_x$count) * 1.5) +
  labs(
    title = glue(
      "Demultiplexing results (souporcell)\n{length(unique(sc_x$batch))} batches, {comma(sum(sc_x$count))} cells"
    ),
    x = "Cells",
    y = NULL
  )
my_ggsave(
  "souporcell-by-batch",
  out_dir = glue("{out_dir}/figures/demux"),
  type = "pdf",
  plot = p,
  scale = 2, width = 5, height = 5, units = "in", dpi = 300
)
#
p <- ggplot(sc_x) +
  aes(x = count) +
  geom_histogram(bins = 50) +
  # geom_vline(xintercept = median(sc_x$count)) +
  facet_wrap(~ batch) +
  # annotate(
  #   geom = "text",
  #   x = median(sc_x$count), y = Inf,
  #   vjust = 1.1, hjust = -0.1,
  #   label = median(sc_x$count),
  #   size = 5
  # ) +
  scale_x_continuous(labels = label_number_si()) +
  labs(
    title = glue(
      "Demultiplexing results (souporcell)\n{length(unique(sc_x$batch))} batches, {comma(sum(sc_x$count))} cells"
    ),
    x = "Cells per sample"
  )
my_ggsave(
  "souporcell-histogram",
  out_dir = glue("{out_dir}/figures/demux"),
  type = "pdf",
  plot = p,
  scale = 2, width = 5, height = 3, units = "in", dpi = 300
)


# demuxEM
########################################################################

adt_files <- Sys.glob(
  glue(
    "/projects/irae_blood/cellranger_output/C*_fbc/*_fbc.csv"
  )
)
adt_batches <- basename(dirname(adt_files))

read_hto <- function(adt_file) {
  hto <- read_sparse_csv(adt_file)
  hto <- hto[stringr::str_detect(rownames(hto), "^HT_"),]
  hto[,Matrix::colSums(hto) > 0]
}
#
call_hto <- function(adt_file) {
  batch <- basename(dirname(adt_file))
  hto <- read_hto(adt_file)
  res <- data.frame(
    # batch = batch,
    cell = glue("{batch}|{colnames(hto)}"),
    demux = do_demux(hto)
  )
}
dm_file <- glue("{out_dir}/demuxem.qs")
if (!file.exists(dm_file)) {
  dm <- data.table::rbindlist(pblapply(adt_files, call_hto))
  dm$cell <- str_remove(dm$cell, "_fbc")
  dm$cell <- str_remove(dm$cell, "_5p")
  # dm$cell <- str_replace(dm$cell, "\\|", "_gex|")
  dm$batch <- str_remove(dm$cell, "\\|.+$")
  dm$batch <- str_remove(dm$batch, "_gex")
  dm$cell <- str_remove(dm$cell, "_gex")
  qsave(dm, dm_file)
} else {
  dm <- qread(dm_file)
}

dm_count <- dm[,.(count = .N), by = .(batch, demux)]
dm_count$batch <- naturalfactor(dm_count$batch)

# Upset comparison
########################################################################

# ix <- obs$batch == "C3_CD45B"

# Pass QC cells, only for demuxEM batches
s_passqc <- obs$cell#[ix]
#
# Singlets from Souporcell, only for demuxEM batches
s_sc     <- sc$cell[
  sc$cell %in% s_passqc &
  # sc$batch %in% unique(obs$batch) &
  sc$status == "singlet"
]
#
# Singlets from demuxEM
s_dm     <- dm$cell[
  sc$cell %in% s_passqc &
  # dm$batch %in% unique(obs$batch) &
  !(dm$demux %in% c("Negative", "Doublet"))
]
#
length(s_passqc)
length(intersect(s_passqc, s_sc))
length(intersect(s_passqc, s_dm))

x <- data.frame(
  # status = c(
  #   "passqc",
  #   "souporcell",
  #   "demuxem",
  #   "passqc and souporcell",
  #   "passqc and demuxem",
  #   "passqc and souporcell and demuxem"
  # ),
  id = 1:6,
  passqc = c(
    TRUE, FALSE, FALSE, TRUE, TRUE, TRUE
  ),
  souporcell = c(
    FALSE, TRUE, FALSE, TRUE, FALSE, TRUE
  ),
  demuxem = c(
    FALSE, FALSE, TRUE, FALSE, TRUE, TRUE
  ),
  cells = c(
    length(s_passqc),
    length(s_sc),
    length(s_dm),
    length(intersect(s_passqc, s_sc)),
    length(intersect(s_passqc, s_dm)),
    length(intersect(intersect(s_passqc, s_dm), s_sc))
  )
)
x <- x %>% pivot_longer(c(-id, -cells))
x$name <- str_remove(as.character(x$name), "\n")
x$name <- factor(x$name, rev(c("passqc", "souporcell", "demuxem")))

p1 <- ggplot(x %>% select(id, cells) %>% unique) +
  geom_col(aes(x = id, y = cells)) +
  geom_text(
    aes(x = id, y = cells, label = signif(cells / 1e3, 3)),
    vjust = -0.1
  ) +
  labs(x = NULL, y = "Cells (Thousands)") +
  scale_y_continuous(
    labels = function(x) x / 1e3
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p2 <- ggplot(x) +
  geom_point(
    aes(x = id, y = name, fill = value),
    shape = 21, stroke = 0.2, size = 5
  ) +
  scale_fill_manual(
    values = c("white", "grey40"), guide = "none"
  ) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(expand = expansion(0, 0.8)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 0)
  )
p <- p1 / p2 + plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = glue("Cell counts for {length(unique(dm$batch))} batches")
  )
my_ggsave(
  "upset-passqc-souporcell-demux",
  out_dir = glue("{out_dir}/figures/demux"),
  type = "pdf",
  plot = p,
  scale = 1, width = 5, height = 5, units = "in", dpi = 300
)


# demuxEM figures
########################################################################

dm_x <- dm[
  !demux %in% c("Doublet", "Negative") &
  cell %in% obs$cell
][
  ,.(count = .N), by = .(batch, demux)
]
dm_x$batch <- naturalfactor(dm_x$batch)

#
p <- ggplot(dm_x) +
  aes(y = demux, x = count, label = count) +
  geom_colh() +
  geom_text(
    hjust = -0.1,
    size = 5
  ) +
  facet_wrap(~ batch) +
  scale_x_continuous(labels = label_number_si()) +
  expand_limits(x = max(dm_x$count) * 1.5) +
  labs(
    title = glue(
      "Demultiplexing results (demuxEM)\n{length(unique(dm_x$batch))} batches, {comma(sum(dm_x$count))} cells"
    ),
    x = "Cells",
    y = NULL
  )
my_ggsave(
  "demuxem-by-batch",
  out_dir = glue("{out_dir}/figures/demux"),
  type = "pdf",
  plot = p,
  scale = 2, width = 5, height = 5, units = "in", dpi = 300
)
#
p <- ggplot(dm_x) +
  aes(x = count) +
  geom_histogram(bins = 50) +
  # geom_vline(xintercept = median(dm_x$count)) +
  facet_wrap(~ batch) +
  # annotate(
  #   geom = "text",
  #   x = median(dm_x$count), y = Inf,
  #   vjust = 1.1, hjust = -0.1,
  #   label = median(dm_x$count),
  #   size = 5
  # ) +
  scale_x_continuous(labels = label_number_si()) +
  labs(
    title = glue(
      "Demultiplexing results (demuxEM)\n{length(unique(dm_x$batch))} batches, {comma(sum(dm_x$count))} cells"
    ),
    x = "Cells per sample"
  )
my_ggsave(
  "demuxem-histogram",
  out_dir = glue("{out_dir}/figures/demux"),
  type = "pdf",
  plot = p,
  scale = 2, width = 5, height = 3, units = "in", dpi = 300
)

my_batches <- unique(dm$batch)
if (FALSE) {
  # demuxEM UMAP
  for (my_batch in my_batches) {
    message(my_batch)
    x <- dm[batch == my_batch]
    #
    adt_file <- Sys.glob(
      glue(
        "/projects/irae_blood/cellranger_output/{my_batch}_fbc/*_fbc.csv"
      )
    )
    read_hto <- function(adt_file) {
      hto <- read_sparse_csv(adt_file)
      hto <- hto[stringr::str_detect(rownames(hto), "^HT_"),]
      hto[,Matrix::colSums(hto) > 10]
    }
    hto <- read_hto(adt_file)
    #
    hto <- t(hto)
    rownames(hto) <- sprintf("%s|%s", my_batch, rownames(hto))
    ix <- rownames(hto) %in% x$cell
    stopifnot(sum(ix) > 0)
    hto <- hto[ix,]
    #
    hto <- hto / rowSums(hto)
    hto <- as.matrix(hto)
    #
    x <- x[x$cell %in% rownames(hto),]
    stopifnot(all(x$cell == rownames(hto)))
    #
    ix_pos <- x$demux != "Negative"
    x <- x[ix_pos,]
    hto <- hto[ix_pos,]
    #
    ix_passqc <- x$cell %in% obs$cell
    x <- x[ix_passqc,]
    hto <- hto[ix_passqc,]
    #
    x_umap <- uwot::umap(
      X = log2(1 + hto),
      n_threads = 8
    ) %>% as.data.table
    x_umap$assignment <- x$demux
    x_umap$assignment.x <- x$demux
    x_umap$cell <- x$cell
    #
    p <- ggplot(x_umap[sample(nrow(x_umap)),]) +
      geom_scattermore(
        aes(x = V1, y = V2, color = append_n(assignment.x)),
        pixels = c(900, 900), alpha = 0.3,
        pointsize = 2.5
      ) +
      scale_color_manual(values = c("grey70", mpn65), name = "Cells") +
      guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
      labs(
        title = glue("{my_batch} - UMAP on HT proportions"),
        x = "UMAP1", y = "UMAP2"
      )
    my_ggsave(
      slug    = glue("{my_batch}_hto-umap"),
      out_dir = glue("{out_dir}/figures/demux"),
      plot    = p,
      type    = "pdf",
      scale   = 1.0, width = 7, height = 4, units = "in", dpi = 300
    )
  }
}


# Merge demuxEM and Souporcell
########################################################################

x_d <- dm[!demux %in% c("Negative", "Doublet")]
x_s <- sc[status == "singlet"]

b1 <- naturalsort(names(table(x_d$batch)))
b2 <- naturalsort(names(table(x_s$batch)))
b3 <- naturalsort(names(table(obs$batch)))
# setdiff(b3, b1)
# setdiff(b3, b2)

all(obs$cell %in% x_s$cell)
all(obs$cell %in% x_d$cell)

common_cells <- intersect(intersect(x_d$cell, x_s$cell), obs$cell)
x_d <- x_d[x_d$cell %in% common_cells]
x_s <- x_s[x_s$cell %in% common_cells]
# obs <- obs[obs$cell %in% common_cells]

m <- merge.data.table(x_d, x_s, by = "cell", all = TRUE)
m <- data.table(table(m$demux, m$assignment, m$batch.y))
colnames(m) <- c("demux", "soup", "batch", "n")

m <- m %>%
  group_by(batch, demux) %>%
  mutate(
    n_demux = sum(n),
    pct_demux = 100 * n / (sum(n) + 1)
  )
m <- m %>%
  group_by(batch, soup) %>%
  mutate(
    n_soup = sum(n),
    pct_soup = 100 * n / (sum(n) + 1)
  ) %>%
  mutate(id = glue("{batch}_{demux}"))

dir.create(glue("{out_dir}/demux-match"), showWarnings = FALSE)
fwrite(m, glue("{out_dir}/demux-match/merge-demux-souporcell.tsv"), sep = "\t")

length(unique(m$id))

m_id <- m %>% group_by(batch, soup) %>%
  top_n(n = 1, wt = pct_soup + pct_demux) %>%
  # filter(pct_soup > 50 & pct_demux > 50) %>%
  mutate(id = glue("{batch}_{demux}"))
m_id$keep <- with(m_id, pct_soup > 50 & pct_demux > 50)
p <- ggplot() +
  geom_hline(yintercept = 50, size = 0.2, linetype = 2) +
  geom_vline(xintercept = 50, size = 0.2, linetype = 2) +
  geom_point(
    data = m_id,
    mapping = aes(pct_demux, pct_soup, color = append_n(as.character(keep)))
  ) +
  scale_color_manual(
    values = c("grey60", "grey20"),
    name = "Keep"
  ) +
  labs(
    x = "demuxEM", y = "Souporcell",
    title = "Matches between Souporcell and demuxEM",
    subtitle = glue(
      "Percent of cells in the best match, total of {nrow(m_id)} matches"
    )
  )
my_ggsave(
  slug    = "scatter-pct_demux-pct_soup",
  out_dir = glue("{out_dir}/figures/demux"),
  plot    = p,
  type    = "pdf",
  scale   = 1.0, width = 6, height = 4, units = "in", dpi = 300
)

# How many cells are discarded?
m_id %<>% mutate(batch_soup = glue("{batch}_{soup}"))
sc %<>% mutate(batch_soup = glue("{batch}_{assignment}"))

sum(
  sc$status == "singlet" &
  sc$cell %in% obs$cell
)
sum(
  sc$status == "singlet" &
  sc$cell %in% obs$cell &
  sc$batch_soup %in% m_id$batch_soup[m_id$keep]
)

p <- ggplot() +
  geom_tile(
    data = m,
    mapping = aes(x = soup, y = demux, fill = n)
  ) +
  geom_point(
    data = m_id,
    mapping = aes(x = soup, y = demux),
    color = "white", size = 1
  ) +
  scale_fill_gradientn(
    colors = scico::scico(n = 100, palette = "batlow"),
    na.value = "grey10",
    # na.value = "black",
    trans = "log10",
    name = "Cells"
  ) +
  guides(fill = guide_colorbar(barheight = 25)) +
  facet_wrap(~ naturalfactor(batch)) +
  labs(
    x = "Souporcell",
    y = "demuxEM",
    title = "Matching singlets between demuxEM and Souporcell",
    subtitle = glue(
      "Number of cells in each match, {length(unique(m$batch))} batches, {comma(sum(m$n))} cells"
    )
  )
my_ggsave(
  slug    = "matching-singlets",
  out_dir = glue("{out_dir}/figures/demux"),
  plot    = p,
  type    = "pdf",
  scale   = 1.0, width = 12, height = 9, units = "in", dpi = 300
)




# Label the cells with the demultiplexed sample
########################################################################
if (!"demux" %in% colnames(obs)) {
  obs <- left_join(obs, dm, by = c("cell", "batch"))
}
ix <- which(!obs$demux %in% c("Doublet", "Negative"))
stopifnot(all(colnames(counts) == obs$cell))
obs <- obs[ix,]
counts <- counts[,ix]
ix <- !is.na(obs$demux)
obs <- obs[ix,]
counts <- counts[,ix]
stopifnot(all(colnames(counts) == obs$cell))
obs$batch_demux <- sprintf("%s_%s", str_split_fixed(obs$batch, "_", 2)[,1], obs$demux)
obs$donor <- hashtag_to_donor[obs$batch_demux]

obs <- left_join(obs, donor_info, by = "donor")

source("R/load-sample-therapy.R")
if (!"drug" %in% colnames(obs)) {
  obs <- left_join(obs, sample_therapy, by = "donor")
}
stopifnot(all(colnames(counts) == obs$cell))

obs %>% count(donor, case, drug, demux)

# }}}

# Run analyses {{{

########################################################################
analysis_name <- "blood1"
a1_file <- as.character(glue("results/blood/{analysis_name}/{analysis_name}.qs"))
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
    n_harmony     = 0,
    # harmony_vars  = c("channel"),
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
    obs           = obs[ix_keep,],
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
  out_dir       = glue("results/blood/{analysis_name}/figures"),
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
  out_dir           = glue("results/blood/{analysis_name}/figures"),
  cb_dir            = glue("cellbrowser/colitis/blood/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)



########################################################################
analysis_name <- "blood2"
a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
file.exists(a1_file)
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
source("R/functions/do-analysis.R")
if (!file.exists(a1_file)) {
  ix_keep <- rep(TRUE, ncol(counts))
  #
  params <- list(
    analysis_name = analysis_name,
    n_genes       = nrow(counts[,ix_keep]),
    n_cells       = ncol(counts[,ix_keep]),
    min_percent   = 100 * (50 / ncol(counts[,ix_keep])),
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 50,
    n_harmony     = 30,
    harmony_vars  = c("batch"),
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
    obs           = obs[ix_keep,],
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

if (FALSE) {
  # Fix the donor labels
  table(a1$obs$donor)
  a1$obs$batch_demux <- sprintf("%s_%s", str_split_fixed(a1$obs$batch, "_", 2)[,1], a1$obs$demux)
  table(a1$obs$batch_demux)
  a1$obs$donor <- hashtag_to_donor[a1$obs$batch_demux]
  length(unique(a1$obs$donor))
  #
  # SIC_100_443 (from Pool C2) & SIC_43_441 (from Pool C3) should be INCLUDED in the initial QC, clustering and subclustering, but then EXCLUDED
  a1$obs$donor_batch <- sprintf("%s_%s", a1$obs$donor, a1$obs$batch)
  #
  # Reload donor info after fixing donors
  source("R/load-sample-therapy.R")
  keep_cols <- c(
    setdiff(
      x = colnames(a1$obs),
      y = c(colnames(sample_therapy), colnames(donor_info))
    ),
    "donor"
  )
  a1$obs <- a1$obs[,..keep_cols]
  a1$obs <- left_join(a1$obs, donor_info, by = "donor")
  if (!"drug" %in% colnames(a1$obs)) {
    a1$obs <- left_join(a1$obs, sample_therapy, by = "donor")
  }
  stopifnot(all(colnames(a1$counts) == a1$obs$cell))
  # a1$obs %>% count(case, donor)
  a1$obs$leiden <- a1$obs$leiden0.933
  qsave(a1, a1_file)
}

a1$obs %>% count(donor, case, drug) %>% arrange(-n)

source("R/colors-int.R")
source("R/plot-analysis.R")
a1$obs$leiden <- a1$obs$leiden0.933
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = glue("results/a20/{analysis_name}/figures"),
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
  out_dir           = glue("results/a20/{analysis_name}/figures"),
  cb_dir            = glue("cellbrowser/colitis/a20/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)


# Round 1
########################################################################
blood2_celltypes <- list(
  "myeloid" = c(1, 7, 13, 15, 16, 17),
  "tcell" = c(2, 3, 4, 5, 6, 9, 10, 11, 12, 16),
  "bcell" = c(8)
)
#
analysis_name <- "blood2"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))

celltypes <- names(blood2_celltypes)

my_celltype <- "bcell"

for (my_celltype in celltypes) {

ix_keep <- which(a1$obs$leiden0.933 %in% blood2_celltypes[[my_celltype]])
length(ix_keep)

analysis_name <- sprintf("blood2_%s1", my_celltype)
c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
if (!file.exists(c1_file)) {
  dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
  #
  source("R/functions/do-analysis.R")
  if (!file.exists(c1_file)) {
    #
    params <- list(
      analysis_name = analysis_name,
      n_genes       = nrow(a1$counts[,ix_keep]),
      n_cells       = ncol(a1$counts[,ix_keep]),
      min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
      loess_span    = 0.01,
      n_pcs         = 'mcv',
      max_pcs       = 50,
      n_harmony     = 30,
      harmony_vars  = c("batch"),
      n_knn         = 50,
      leiden_res    = seq(0.5, 1.8, length.out = 10),
      leiden_iter   = 10,
      umap_spread   = 1,
      umap_min_dist = 0.25,
      log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
    )
    c1_params_file <- glue("{dirname(c1_file)}/params.json")
    writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
    #
    keep_cols <- colnames(a1$obs)
    keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
    c1 <- run_analysis(
      obs           = a1$obs[ix_keep,..keep_cols],
      counts        = a1$counts[,ix_keep],
      exclude_genes = c(mito_genes, bcr_genes),
      mito_genes    = mito_genes,
      params        = params
    )
    c1$params <- params
    print_status(glue("Writing {c1_file}"))
    qsave(c1, c1_file)
    print_status(glue("done"))
  } else {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
  }
  #
  source("R/colors-int.R")
  source("R/plot-analysis.R")
  c1$obs$leiden <- c1$obs$leiden0.933
  assign(analysis_name, c1)
  try({
  plot_analysis(
    analysis_name = analysis_name,
    out_dir       = glue("results/a20/{analysis_name}/figures"),
    rowname_key   = ensembl_to_symbol,
    exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
    do_pb         = FALSE
  )
  })
}

}

# # blood2_bcell1
# c1$obs %>% count(leiden, donor, batch) %>% filter(donor == "SIC_140")
# c1$obs %>% count(leiden)
# (1010 + 960) / 2020

# Round 2
########################################################################
cell_ids <- list("tcell" = c(), "bcell" = c(), "myeloid" = c())
xs <- lapply(c("tcell", "bcell", "myeloid"), function(my_celltype) {
  analysis_name <- sprintf("blood2_%s1", my_celltype)
  qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))
})
names(xs) <- c("tcell", "bcell", "myeloid")

cell_ids[["tcell"]] <- unique(c(
  with(xs[["tcell"]]$obs, cell[leiden0.933 %in% 1:10]),
  with(xs[["bcell"]]$obs, cell[leiden0.933 %in% c(5, 6)]),
  with(xs[["myeloid"]]$obs, cell[leiden0.933 %in% c(9, 10)])
))
cell_ids[["bcell"]] <- unique(c(
  with(xs[["tcell"]]$obs, cell[leiden0.933 %in% c(12)]),
  with(xs[["bcell"]]$obs, cell[leiden0.933 %in% c(1,2,3,4,8)]),
  with(xs[["myeloid"]]$obs, cell[leiden0.933 %in% c(12)])
))
cell_ids[["myeloid"]] <- unique(c(
  with(xs[["tcell"]]$obs, cell[leiden0.933 %in% c(11, 13)]),
  with(xs[["bcell"]]$obs, cell[leiden0.933 %in% c(7)]),
  with(xs[["myeloid"]]$obs, cell[leiden0.933 %in% c(1,2,3,4,5,6,7,8,11,13,14)])
))
sapply(cell_ids, length)


for (my_celltype in celltypes) {

ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
length(ix_keep)

analysis_name <- sprintf("blood2_%s2", my_celltype)
c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
# if (!file.exists(c1_file)) {
  dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
  #
  source("R/functions/do-analysis.R")
  if (!file.exists(c1_file)) {
    #
    params <- list(
      analysis_name = analysis_name,
      n_genes       = nrow(a1$counts[,ix_keep]),
      n_cells       = ncol(a1$counts[,ix_keep]),
      min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
      loess_span    = 0.01,
      n_pcs         = 'mcv',
      max_pcs       = 50,
      n_harmony     = 30,
      harmony_vars  = c("batch"),
      n_knn         = 50,
      leiden_res    = seq(0.5, 1.8, length.out = 10),
      leiden_iter   = 10,
      umap_spread   = 1,
      umap_min_dist = 0.25,
      log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
    )
    c1_params_file <- glue("{dirname(c1_file)}/params.json")
    writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
    #
    keep_cols <- colnames(a1$obs)
    keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
    c1 <- run_analysis(
      obs           = a1$obs[ix_keep,..keep_cols],
      counts        = a1$counts[,ix_keep],
      exclude_genes = c(mito_genes, bcr_genes),
      mito_genes    = mito_genes,
      params        = params
    )
    c1$params <- params
    print_status(glue("Writing {c1_file}"))
    qsave(c1, c1_file)
    print_status(glue("done"))
  } else {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
  }
  #
  source("R/colors-int.R")
  source("R/plot-analysis.R")
  # c1$obs$leiden <- c1$obs$leiden0.933
  c1$obs$leiden <- c1$obs$leiden1.51
  assign(analysis_name, c1)
  try({
  plot_analysis(
    analysis_name = analysis_name,
    out_dir       = glue("results/a20/{analysis_name}/figures"),
    rowname_key   = ensembl_to_symbol,
    exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
    do_pb         = FALSE
  )
  })
# }

}


# Round 3
########################################################################
xs <- lapply(c("tcell", "bcell", "myeloid"), function(my_celltype) {
  analysis_name <- sprintf("blood2_%s2", my_celltype)
  qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))
})
names(xs) <- c("tcell", "bcell", "myeloid")

cell_ids <- list()
cell_ids[["tcell"]] <- unique(c(
  with(
    xs[["tcell"]]$obs,
    cell[
      leiden1.51 %in% c(1:12) &
      !(batch %in% c("C3_CD45B"))
    ]
  )
))
cell_ids[["bcell"]] <- unique(c(
  with(
    xs[["bcell"]]$obs,
    cell[
      leiden1.51 %in% c(1:8) &
      !(batch %in% c("C1_CD8A", "C1_CD8B", "C2_CD8A", "C2_CD8B", "C3_CD8A", "C3_CD8B"))
    ]
  )
))
cell_ids[["myeloid"]] <- unique(c(
  with(
    xs[["myeloid"]]$obs,
    cell[
      leiden1.51 %in% c(1:16, 20) &
      !(batch %in% c("C1_CD8A", "C1_CD8B", "C2_CD8A", "C2_CD8B", "C3_CD8A", "C3_CD8B"))
    ]
  )
))
sapply(cell_ids, length)


for (my_celltype in celltypes) {

ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
length(ix_keep)

analysis_name <- sprintf("blood2_%s3", my_celltype)
c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
# if (!file.exists(c1_file)) {
  dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
  #
  source("R/functions/do-analysis.R")
  if (!file.exists(c1_file)) {
    #
    params <- list(
      analysis_name = analysis_name,
      n_genes       = nrow(a1$counts[,ix_keep]),
      n_cells       = ncol(a1$counts[,ix_keep]),
      min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
      loess_span    = 0.01,
      n_pcs         = 'mcv',
      max_pcs       = 50,
      n_harmony     = 30,
      harmony_vars  = c("batch"),
      n_knn         = 50,
      leiden_res    = seq(0.5, 1.8, length.out = 10),
      leiden_iter   = 10,
      umap_spread   = 1,
      umap_min_dist = 0.25,
      log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
    )
    c1_params_file <- glue("{dirname(c1_file)}/params.json")
    writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
    #
    keep_cols <- colnames(a1$obs)
    keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
    c1 <- run_analysis(
      obs           = a1$obs[ix_keep,..keep_cols],
      counts        = a1$counts[,ix_keep],
      exclude_genes = c(mito_genes, bcr_genes),
      mito_genes    = mito_genes,
      params        = params
    )
    c1$params <- params
    print_status(glue("Writing {c1_file}"))
    qsave(c1, c1_file)
    print_status(glue("done"))
  } else {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
  }
  #
  source("R/colors-int.R")
  source("R/plot-analysis.R")
  # c1$obs$leiden <- c1$obs$leiden0.933
  c1$obs$leiden <- c1$obs$leiden1.51
  assign(analysis_name, c1)
  try({
  plot_analysis(
    analysis_name = analysis_name,
    out_dir       = glue("results/a20/{analysis_name}/figures"),
    rowname_key   = ensembl_to_symbol,
    exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
    do_pb         = FALSE
  )
  })
# }

}


# Round 4
########################################################################
xs <- lapply(c("tcell", "bcell", "myeloid"), function(my_celltype) {
  analysis_name <- sprintf("blood2_%s3", my_celltype)
  qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))
})
names(xs) <- c("tcell", "bcell", "myeloid")

cell_ids <- list()
cell_ids[["tcell"]] <- unique(c(
  with(
    xs[["tcell"]]$obs,
    cell[
      leiden1.51 %in% c(1:13, 15) &
      !(batch %in% c("C3_CD45B"))
    ]
  )
))
cell_ids[["bcell"]] <- unique(c(
  with(
    xs[["bcell"]]$obs,
    cell[
      !(batch %in% c("C1_CD8A", "C1_CD8B", "C2_CD8A", "C2_CD8B", "C3_CD8A", "C3_CD8B"))
    ]
  )
))
cell_ids[["myeloid"]] <- unique(c(
  with(
    xs[["myeloid"]]$obs,
    cell[
      !(batch %in% c("C1_CD8A", "C1_CD8B", "C2_CD8A", "C2_CD8B", "C3_CD8A", "C3_CD8B"))
    ]
  )
))
sapply(cell_ids, length)

celltypes <- names(cell_ids)

for (my_celltype in celltypes) {

ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
length(ix_keep)

analysis_name <- sprintf("blood2_%s4", my_celltype)
c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
# if (!file.exists(c1_file)) {
  dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
  #
  source("R/functions/do-analysis.R")
  if (!file.exists(c1_file)) {
    #
    params <- list(
      analysis_name = analysis_name,
      n_genes       = nrow(a1$counts[,ix_keep]),
      n_cells       = ncol(a1$counts[,ix_keep]),
      min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
      loess_span    = 0.01,
      n_pcs         = 'mcv',
      max_pcs       = 50,
      n_harmony     = 30,
      harmony_vars  = c("batch"),
      n_knn         = 50,
      leiden_res    = seq(0.5, 1.8, length.out = 10),
      leiden_iter   = 10,
      umap_spread   = 1,
      umap_min_dist = 0.25,
      log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
    )
    c1_params_file <- glue("{dirname(c1_file)}/params.json")
    writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
    #
    keep_cols <- colnames(a1$obs)
    keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
    c1 <- run_analysis(
      obs           = a1$obs[ix_keep,..keep_cols],
      counts        = a1$counts[,ix_keep],
      exclude_genes = c(mito_genes, bcr_genes),
      mito_genes    = mito_genes,
      params        = params
    )
    c1$params <- params
    print_status(glue("Writing {c1_file}"))
    qsave(c1, c1_file)
    print_status(glue("done"))
  } else {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
  }
  #
  source("R/colors-int.R")
  source("R/plot-analysis.R")
  # c1$obs$leiden <- c1$obs$leiden0.933
  c1$obs$leiden <- c1$obs$leiden1.51
  assign(analysis_name, c1)
  try({
  plot_analysis(
    analysis_name = analysis_name,
    out_dir       = glue("results/a20/{analysis_name}/figures"),
    rowname_key   = ensembl_to_symbol,
    exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
    do_pb         = FALSE
  )
  })
# }

}


# Round 5
########################################################################
xs <- lapply(c("tcell", "bcell", "myeloid"), function(my_celltype) {
  analysis_name <- sprintf("blood2_%s4", my_celltype)
  qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))
})
names(xs) <- c("tcell", "bcell", "myeloid")

cell_ids <- list()
cell_ids[["tcell"]] <- unique(c(
  with(
    xs[["tcell"]]$obs,
    cell[
      leiden1.51 %in% c(1:15)
    ]
  )
))
cell_ids[["myeloid"]] <- unique(c(
  with(
    xs[["myeloid"]]$obs,
    cell[
      leiden1.51 %in% c(1:15)
    ]
  )
))
sapply(cell_ids, length)

for (my_celltype in names(cell_ids)) {

  ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
  length(ix_keep)

  analysis_name <- sprintf("blood2_%s5", my_celltype)
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  # if (!file.exists(c1_file)) {
    dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
    #
    source("R/functions/do-analysis.R")
    if (!file.exists(c1_file)) {
      #
      params <- list(
        analysis_name = analysis_name,
        n_genes       = nrow(a1$counts[,ix_keep]),
        n_cells       = ncol(a1$counts[,ix_keep]),
        min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
        loess_span    = 0.01,
        n_pcs         = 'mcv',
        max_pcs       = 50,
        n_harmony     = 30,
        harmony_vars  = c("batch"),
        n_knn         = 50,
        leiden_res    = seq(0.5, 1.8, length.out = 10),
        leiden_iter   = 10,
        umap_spread   = 1,
        umap_min_dist = 0.25,
        log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
      )
      c1_params_file <- glue("{dirname(c1_file)}/params.json")
      writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
      #
      keep_cols <- colnames(a1$obs)
      keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
      c1 <- run_analysis(
        obs           = a1$obs[ix_keep,..keep_cols],
        counts        = a1$counts[,ix_keep],
        exclude_genes = c(mito_genes, bcr_genes),
        mito_genes    = mito_genes,
        params        = params
      )
      c1$params <- params
      print_status(glue("Writing {c1_file}"))
      qsave(c1, c1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {c1_file}"))
      c1 <- qread(c1_file)
      print_status(glue("done"))
    }
    #
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # c1$obs$leiden <- c1$obs$leiden0.933
    my_leiden <- "leiden1.51"
    if (my_celltype == "myeloid") {
      # my_leiden <- "leiden0.933"
      my_leiden <- "leiden0.789"
      # with(c1$obs, table(leiden0.644, leiden0.933))
    } else if (my_celltype == "tcell") {
      my_leiden <- "leiden1.08"
    }
    c1$obs$leiden <- c1$obs[[my_leiden]]
    assign(analysis_name, c1)
    try({
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = glue("results/a20/{analysis_name}/figures_{my_leiden}"),
      rowname_key   = ensembl_to_symbol,
      exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
      do_pb         = FALSE
    )
    })
  # }

}


# Round 6: Attempt to split CD4 and CD8
########################################################################

analysis_name <- "blood2_tcell5"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))
#
cell_ids <- list(
  "cd4" = unique(c(
    with(
      a1$obs,
      cell[
        leiden1.66 %in% c(1, 5, 12)
      ]
    )
  )),
  "cd8" = unique(c(
    with(
      a1$obs,
      cell[
        !leiden1.66 %in% c(1, 5, 12)
      ]
    )
  ))
)
sapply(cell_ids, length)

for (my_celltype in names(cell_ids)) {

  ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
  length(ix_keep)

  analysis_name <- sprintf("blood2_tcell5_%s_1", my_celltype)
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  # if (!file.exists(c1_file)) {
    dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
    #
    source("R/functions/do-analysis.R")
    if (!file.exists(c1_file)) {
      #
      params <- list(
        analysis_name = analysis_name,
        n_genes       = nrow(a1$counts[,ix_keep]),
        n_cells       = ncol(a1$counts[,ix_keep]),
        min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
        loess_span    = 0.01,
        n_pcs         = 'mcv',
        max_pcs       = 50,
        n_harmony     = 30,
        harmony_vars  = c("batch"),
        n_knn         = 50,
        leiden_res    = seq(0.5, 1.8, length.out = 10),
        leiden_iter   = 10,
        umap_spread   = 1,
        umap_min_dist = 0.25,
        log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
      )
      c1_params_file <- glue("{dirname(c1_file)}/params.json")
      writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
      #
      keep_cols <- colnames(a1$obs)
      keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
      c1 <- run_analysis(
        obs           = a1$obs[ix_keep,..keep_cols],
        counts        = a1$counts[,ix_keep],
        exclude_genes = c(mito_genes, bcr_genes),
        mito_genes    = mito_genes,
        params        = params
      )
      c1$params <- params
      print_status(glue("Writing {c1_file}"))
      qsave(c1, c1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {c1_file}"))
      c1 <- qread(c1_file)
      print_status(glue("done"))
    }
    #
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # c1$obs$leiden <- c1$obs$leiden0.933
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    out_dir <- glue("results/a20/{analysis_name}/figures_{my_leiden}")
    assign(analysis_name, c1)
    try({
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = out_dir,
      rowname_key   = ensembl_to_symbol,
      exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
      do_pb         = FALSE
    )
    })
  # }

    c1$adt_clr <- a1$adt_clr[,ix_keep]

    # x1 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD8 ADT_C0046",]),
    #   by = "cell"
    # )
    # x2 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD4 ADT_C0072",]),
    #   by = "cell"
    # )

    if ("adt_clr" %in% names(c1)) {
      for (my_protein in c("CD8 ADT_C0046", "CD4 ADT_C0072")) {
          x <- left_join(
            x = c1$obs %>% select(cell, UMAP1, UMAP2),
            y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
            by = "cell"
          )
          x$protein[is.na(x$protein)] <- 0
          p <- plot_hexgene(
            x            = x$UMAP1,
            y            = x$UMAP2,
            z            = x$protein,
            # z            = x1$protein / x2$protein,
            bins         = 201,
            palette      = "davos",
            direction    = -1,
            use_quantile = FALSE,
            legend       = FALSE,
            italic       = FALSE
          ) +
          guides(
            fill = guide_colorbar(
              direction = "horizontal",
              title.position = "top",
              title = "CLR",
              barwidth = 15
            )
          ) +
          labs(subtitle = my_protein)
          my_ggsave(
            glue("umap-{safe(my_protein)}"),
            # glue("umap-CD8-to-CD4"),
            out_dir = glue("{out_dir}/adt-umap"),
            type = "pdf",
            plot = p,
            scale = 1, width = 4, height = 5, units = "in", dpi = 300
          )
      }
    }

}

# Round 7: Continue to split CD4 and CD8
########################################################################

analysis_name <- "blood2_tcell5"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))

#
cell_ids <- list(
  "cd4" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_1/data/blood2_tcell5_cd4_1.qs")))$obs,
      cell[
        leiden1.51 %in% c(1, 2, 3, 4, 6, 7, 8, 9, 10, 12)
      ]
    )
  )),
  "cd8" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_1/data/blood2_tcell5_cd4_1.qs")))$obs,
      cell[
        leiden1.51 %in% c(5, 11)
      ]
    ),
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd8_1/data/blood2_tcell5_cd8_1.qs")))$obs,
      cell
    )
  ))
)
sapply(cell_ids, length)

for (my_celltype in names(cell_ids)) {

  ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
  length(ix_keep)

  analysis_name <- sprintf("blood2_tcell5_%s_2", my_celltype)
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  # if (!file.exists(c1_file)) {
    dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
    #
    source("R/functions/do-analysis.R")
    if (!file.exists(c1_file)) {
      #
      params <- list(
        analysis_name = analysis_name,
        n_genes       = nrow(a1$counts[,ix_keep]),
        n_cells       = ncol(a1$counts[,ix_keep]),
        min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
        loess_span    = 0.01,
        n_pcs         = 'mcv',
        max_pcs       = 50,
        n_harmony     = 30,
        harmony_vars  = c("batch"),
        n_knn         = 50,
        leiden_res    = seq(0.5, 1.8, length.out = 10),
        leiden_iter   = 10,
        umap_spread   = 1,
        umap_min_dist = 0.25,
        log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
      )
      c1_params_file <- glue("{dirname(c1_file)}/params.json")
      writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
      #
      keep_cols <- colnames(a1$obs)
      keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
      c1 <- run_analysis(
        obs           = a1$obs[ix_keep,..keep_cols],
        counts        = a1$counts[,ix_keep],
        exclude_genes = c(mito_genes, bcr_genes),
        mito_genes    = mito_genes,
        params        = params
      )
      c1$params <- params
      print_status(glue("Writing {c1_file}"))
      qsave(c1, c1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {c1_file}"))
      c1 <- qread(c1_file)
      print_status(glue("done"))
    }
    #
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # c1$obs$leiden <- c1$obs$leiden0.933
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    out_dir <- glue("results/a20/{analysis_name}/figures_{my_leiden}")
    assign(analysis_name, c1)
    try({
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = out_dir,
      rowname_key   = ensembl_to_symbol,
      exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
      do_pb         = FALSE
    )
    })
  # }

    c1$adt_clr <- a1$adt_clr[,ix_keep]

    # x1 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD8 ADT_C0046",]),
    #   by = "cell"
    # )
    # x2 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD4 ADT_C0072",]),
    #   by = "cell"
    # )

    if ("adt_clr" %in% names(c1)) {
      for (my_protein in c("CD8 ADT_C0046", "CD4 ADT_C0072")) {
          x <- left_join(
            x = c1$obs %>% select(cell, UMAP1, UMAP2),
            y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
            by = "cell"
          )
          x$protein[is.na(x$protein)] <- 0
          p <- plot_hexgene(
            x            = x$UMAP1,
            y            = x$UMAP2,
            z            = x$protein,
            # z            = x1$protein / x2$protein,
            bins         = 201,
            palette      = "davos",
            direction    = -1,
            use_quantile = FALSE,
            legend       = FALSE,
            italic       = FALSE
          ) +
          guides(
            fill = guide_colorbar(
              direction = "horizontal",
              title.position = "top",
              title = "CLR",
              barwidth = 15
            )
          ) +
          labs(subtitle = my_protein)
          my_ggsave(
            glue("umap-{safe(my_protein)}"),
            # glue("umap-CD8-to-CD4"),
            out_dir = glue("{out_dir}/adt-umap"),
            type = "pdf",
            plot = p,
            scale = 1, width = 4, height = 5, units = "in", dpi = 300
          )
      }
    }

}

# Round 8: Continue to split CD4 and CD8
########################################################################

analysis_name <- "blood2_tcell5"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))

#
cell_ids <- list(
  "cd4" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_2/data/blood2_tcell5_cd4_2.qs")))$obs,
      cell[
        leiden1.51 %in% c(1, 2, 3, 4, 7, 8, 9, 10, 11, 13)
      ]
    )
  )),
  "cd8" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_2/data/blood2_tcell5_cd4_2.qs")))$obs,
      cell[
        leiden1.51 %in% c(5, 6, 12, 14)
      ]
    ),
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd8_2/data/blood2_tcell5_cd8_2.qs")))$obs,
      cell
    )
  ))
)
sapply(cell_ids, length)

for (my_celltype in names(cell_ids)) {

  ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
  length(ix_keep)

  analysis_name <- sprintf("blood2_tcell5_%s_3", my_celltype)
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  # if (!file.exists(c1_file)) {
    dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
    #
    source("R/functions/do-analysis.R")
    if (!file.exists(c1_file)) {
      #
      params <- list(
        analysis_name = analysis_name,
        n_genes       = nrow(a1$counts[,ix_keep]),
        n_cells       = ncol(a1$counts[,ix_keep]),
        min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
        loess_span    = 0.01,
        n_pcs         = 'mcv',
        max_pcs       = 50,
        n_harmony     = 30,
        harmony_vars  = c("batch"),
        n_knn         = 50,
        leiden_res    = seq(0.5, 1.8, length.out = 10),
        leiden_iter   = 10,
        umap_spread   = 1,
        umap_min_dist = 0.25,
        log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
      )
      c1_params_file <- glue("{dirname(c1_file)}/params.json")
      writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
      #
      keep_cols <- colnames(a1$obs)
      keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
      c1 <- run_analysis(
        obs           = a1$obs[ix_keep,..keep_cols],
        counts        = a1$counts[,ix_keep],
        exclude_genes = c(mito_genes, bcr_genes),
        mito_genes    = mito_genes,
        params        = params
      )
      c1$params <- params
      print_status(glue("Writing {c1_file}"))
      qsave(c1, c1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {c1_file}"))
      c1 <- qread(c1_file)
      print_status(glue("done"))
    }
    #
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # c1$obs$leiden <- c1$obs$leiden0.933
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    out_dir <- glue("results/a20/{analysis_name}/figures_{my_leiden}")
    assign(analysis_name, c1)
    try({
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = out_dir,
      rowname_key   = ensembl_to_symbol,
      exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
      do_pb         = FALSE
    )
    })
  # }

    c1$adt_clr <- a1$adt_clr[,ix_keep]

    # x1 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD8 ADT_C0046",]),
    #   by = "cell"
    # )
    # x2 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD4 ADT_C0072",]),
    #   by = "cell"
    # )

    if ("adt_clr" %in% names(c1)) {
      for (my_protein in c("CD8 ADT_C0046", "CD4 ADT_C0072")) {
          x <- left_join(
            x = c1$obs %>% select(cell, UMAP1, UMAP2),
            y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
            by = "cell"
          )
          x$protein[is.na(x$protein)] <- 0
          p <- plot_hexgene(
            x            = x$UMAP1,
            y            = x$UMAP2,
            z            = x$protein,
            # z            = x1$protein / x2$protein,
            bins         = 201,
            palette      = "davos",
            direction    = -1,
            use_quantile = FALSE,
            legend       = FALSE,
            italic       = FALSE
          ) +
          guides(
            fill = guide_colorbar(
              direction = "horizontal",
              title.position = "top",
              title = "CLR",
              barwidth = 15
            )
          ) +
          labs(subtitle = my_protein)
          my_ggsave(
            glue("umap-{safe(my_protein)}"),
            # glue("umap-CD8-to-CD4"),
            out_dir = glue("{out_dir}/adt-umap"),
            type = "pdf",
            plot = p,
            scale = 1, width = 4, height = 5, units = "in", dpi = 300
          )
      }
    }

}

for (analysis_name in c("blood2_tcell5_cd8_3", "blood2_tcell5_cd4_3")) {
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    assign(analysis_name, c1)
    source("R/functions/make-cellbrowser.R")
    make_cellbrowser(
      analysis_name     = analysis_name,
      out_dir           = glue("results/a20/{analysis_name}/figures_{my_leiden}"),
      cb_dir            = glue("cellbrowser/colitis/a20/{analysis_name}"),
      ensembl_to_symbol = ensembl_to_symbol
    )
  }
}


# Round 9: Continue to split CD4 and CD8
########################################################################

analysis_name <- "blood2_tcell5"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))

#
cell_ids <- list(
  "cd4" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_3/data/blood2_tcell5_cd4_3.qs")))$obs,
      cell[
        !leiden1.51 %in% c(1)
      ]
    ),
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd8_3/data/blood2_tcell5_cd8_3.qs")))$obs,
      cell[
        leiden1.51 %in% c(13)
      ]
    )
  )),
  "cd8" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_3/data/blood2_tcell5_cd4_3.qs")))$obs,
      cell[
        leiden1.51 %in% c(1)
      ]
    ),
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd8_3/data/blood2_tcell5_cd8_3.qs")))$obs,
      cell[
        !leiden1.51 %in% c(13)
      ]
    )
  ))
)
sapply(cell_ids, length)

for (my_celltype in names(cell_ids)) {

  ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
  length(ix_keep)

  analysis_name <- sprintf("blood2_tcell5_%s_4", my_celltype)
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  # if (!file.exists(c1_file)) {
    dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
    #
    source("R/functions/do-analysis.R")
    if (!file.exists(c1_file)) {
      #
      params <- list(
        analysis_name = analysis_name,
        n_genes       = nrow(a1$counts[,ix_keep]),
        n_cells       = ncol(a1$counts[,ix_keep]),
        min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
        loess_span    = 0.01,
        n_pcs         = 'mcv',
        max_pcs       = 50,
        n_harmony     = 30,
        harmony_vars  = c("batch"),
        n_knn         = 50,
        leiden_res    = seq(0.5, 1.8, length.out = 10),
        leiden_iter   = 10,
        umap_spread   = 1,
        umap_min_dist = 0.25,
        log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
      )
      c1_params_file <- glue("{dirname(c1_file)}/params.json")
      writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
      #
      keep_cols <- colnames(a1$obs)
      keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
      c1 <- run_analysis(
        obs           = a1$obs[ix_keep,..keep_cols],
        counts        = a1$counts[,ix_keep],
        exclude_genes = c(mito_genes, bcr_genes),
        mito_genes    = mito_genes,
        params        = params
      )
      c1$params <- params
      print_status(glue("Writing {c1_file}"))
      qsave(c1, c1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {c1_file}"))
      c1 <- qread(c1_file)
      print_status(glue("done"))
    }
    #
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # c1$obs$leiden <- c1$obs$leiden0.933
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    out_dir <- glue("results/a20/{analysis_name}/figures_{my_leiden}")
    assign(analysis_name, c1)
    try({
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = out_dir,
      rowname_key   = ensembl_to_symbol,
      exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
      do_pb         = FALSE
    )
    })
  # }

    c1$adt_clr <- a1$adt_clr[,ix_keep]

    # x1 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD8 ADT_C0046",]),
    #   by = "cell"
    # )
    # x2 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD4 ADT_C0072",]),
    #   by = "cell"
    # )

    if ("adt_clr" %in% names(c1)) {
      for (my_protein in c("CD8 ADT_C0046", "CD4 ADT_C0072")) {
          x <- left_join(
            x = c1$obs %>% select(cell, UMAP1, UMAP2),
            y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
            by = "cell"
          )
          x$protein[is.na(x$protein)] <- 0
          p <- plot_hexgene(
            x            = x$UMAP1,
            y            = x$UMAP2,
            z            = x$protein,
            # z            = x1$protein / x2$protein,
            bins         = 201,
            palette      = "davos",
            direction    = -1,
            use_quantile = FALSE,
            legend       = FALSE,
            italic       = FALSE
          ) +
          guides(
            fill = guide_colorbar(
              direction = "horizontal",
              title.position = "top",
              title = "CLR",
              barwidth = 15
            )
          ) +
          labs(subtitle = my_protein)
          my_ggsave(
            glue("umap-{safe(my_protein)}"),
            # glue("umap-CD8-to-CD4"),
            out_dir = glue("{out_dir}/adt-umap"),
            type = "pdf",
            plot = p,
            scale = 1, width = 4, height = 5, units = "in", dpi = 300
          )
      }
    }

}

for (analysis_name in c("blood2_tcell5_cd8_4", "blood2_tcell5_cd4_4")) {
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    assign(analysis_name, c1)
    source("R/functions/make-cellbrowser.R")
    make_cellbrowser(
      analysis_name     = analysis_name,
      out_dir           = glue("results/a20/{analysis_name}/figures_{my_leiden}"),
      cb_dir            = glue("cellbrowser/colitis/a20/{analysis_name}"),
      ensembl_to_symbol = ensembl_to_symbol
    )
  }
}


# Round 10: Continue to split CD4 and CD8
########################################################################

analysis_name <- "blood2_tcell5"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))

#
cell_ids <- list(
  "cd4" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_4/data/blood2_tcell5_cd4_4.qs")))$obs,
      cell[
        !leiden1.51 %in% c(12)
      ]
    )
  )),
  "cd8" = unique(c(
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd4_4/data/blood2_tcell5_cd4_4.qs")))$obs,
      cell[
        leiden1.51 %in% c(12)
      ]
    ),
    with(
      qread(as.character(glue("results/a20/blood2_tcell5_cd8_4/data/blood2_tcell5_cd8_4.qs")))$obs,
      cell
    )
  ))
)
sapply(cell_ids, length)

for (my_celltype in names(cell_ids)) {

  ix_keep <- which(a1$obs$cell %in% cell_ids[[my_celltype]])
  length(ix_keep)

  analysis_name <- sprintf("blood2_tcell5_%s_5", my_celltype)
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  # if (!file.exists(c1_file)) {
    dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
    #
    source("R/functions/do-analysis.R")
    if (!file.exists(c1_file)) {
      #
      params <- list(
        analysis_name = analysis_name,
        n_genes       = nrow(a1$counts[,ix_keep]),
        n_cells       = ncol(a1$counts[,ix_keep]),
        min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
        loess_span    = 0.01,
        n_pcs         = 'mcv',
        max_pcs       = 50,
        n_harmony     = 30,
        harmony_vars  = c("batch"),
        n_knn         = 50,
        leiden_res    = seq(0.5, 1.8, length.out = 10),
        leiden_iter   = 10,
        umap_spread   = 1,
        umap_min_dist = 0.25,
        log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
      )
      c1_params_file <- glue("{dirname(c1_file)}/params.json")
      writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
      #
      keep_cols <- colnames(a1$obs)
      keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
      c1 <- run_analysis(
        obs           = a1$obs[ix_keep,..keep_cols],
        counts        = a1$counts[,ix_keep],
        exclude_genes = c(mito_genes, bcr_genes),
        mito_genes    = mito_genes,
        params        = params
      )
      c1$params <- params
      if ("adt_clr" %in% names(a1)) {
        c1$adt_counts <- a1$adt_counts[,ix_keep]
        c1$adt_clr <- a1$adt_clr[,ix_keep]
      }
      print_status(glue("Writing {c1_file}"))
      qsave(c1, c1_file)
      print_status(glue("done"))
    } else {
      print_status(glue("Reading {c1_file}"))
      c1 <- qread(c1_file)
      print_status(glue("done"))
    }
    #
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # c1$obs$leiden <- c1$obs$leiden0.933
    my_leiden <- "leiden1.51"
    c1$obs$leiden <- c1$obs[[my_leiden]]
    out_dir <- glue("results/a20/{analysis_name}/figures_{my_leiden}")
    assign(analysis_name, c1)
    try({
    plot_analysis(
      analysis_name = analysis_name,
      out_dir       = out_dir,
      rowname_key   = ensembl_to_symbol,
      exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
      do_pb         = FALSE
    )
    })
  # }

    # x1 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD8 ADT_C0046",]),
    #   by = "cell"
    # )
    # x2 <- left_join(
    #   x = c1$obs %>% select(cell, UMAP1, UMAP2),
    #   y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr["CD4 ADT_C0072",]),
    #   by = "cell"
    # )

    if ("adt_clr" %in% names(c1)) {
      for (my_protein in c("CD8 ADT_C0046", "CD4 ADT_C0072")) {
          x <- left_join(
            x = c1$obs %>% select(cell, UMAP1, UMAP2),
            y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
            by = "cell"
          )
          x$protein[is.na(x$protein)] <- 0
          p <- plot_hexgene(
            x            = x$UMAP1,
            y            = x$UMAP2,
            z            = x$protein,
            # z            = x1$protein / x2$protein,
            bins         = 201,
            palette      = "davos",
            direction    = -1,
            use_quantile = FALSE,
            legend       = FALSE,
            italic       = FALSE
          ) +
          guides(
            fill = guide_colorbar(
              direction = "horizontal",
              title.position = "top",
              title = "CLR",
              barwidth = 15
            )
          ) +
          labs(subtitle = my_protein)
          my_ggsave(
            glue("umap-{safe(my_protein)}"),
            # glue("umap-CD8-to-CD4"),
            out_dir = glue("{out_dir}/adt-umap"),
            type = "pdf",
            plot = p,
            scale = 1, width = 4, height = 5, units = "in", dpi = 300
          )
      }
    }

}

# for (analysis_name in c("blood2_tcell5_cd8_5", "blood2_tcell5_cd4_5")) {
for (analysis_name in c("blood2_tcell5_cd8_5", "blood2_tcell5_cd4_5", "blood2_myeloid5", "blood2_bcell5")) {
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    # my_leiden <- "leiden1.51"
    # c1$obs$leiden <- c1$obs[[my_leiden]]
    c1$obs$cluster <- c1$obs$leiden
    assign(analysis_name, c1)
    out_dir <- glue("results/a20/{analysis_name}/figures")
    source("R/colors-int.R")
    source("R/plot-analysis.R")
    # try({
    #   plot_analysis(
    #     analysis_name = analysis_name,
    #     out_dir       = out_dir,
    #     rowname_key   = ensembl_to_symbol,
    #     exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
    #     do_pb         = FALSE
    #   )
    # })
    source("R/functions/make-cellbrowser.R")
    make_cellbrowser(
      analysis_name     = analysis_name,
      # out_dir           = glue("results/a20/{analysis_name}/figures_{my_leiden}"),
      out_dir           = out_dir,
      cb_dir            = glue("cellbrowser/colitis/a20/{analysis_name}"),
      ensembl_to_symbol = ensembl_to_symbol
    )
  }
}


# B cells without SIC_140
########################################################################

analysis_name <- "blood2_bcell4"
a1 <- qread(as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs")))

cell_ids <- with(a1$obs, cell[donor != "SIC_140"])
length(cell_ids)

{
  ix_keep <- which(a1$obs$cell %in% cell_ids)
  length(ix_keep)
  #
  analysis_name <- sprintf("blood2_bcell5")
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  file.exists(c1_file)
  dir.create(dirname(c1_file), recursive = TRUE, showWarnings = FALSE)
  #
  source("R/functions/do-analysis.R")
  if (!file.exists(c1_file)) {
    #
    params <- list(
      analysis_name = analysis_name,
      n_genes       = nrow(a1$counts[,ix_keep]),
      n_cells       = ncol(a1$counts[,ix_keep]),
      min_percent   = 100 * (50 / ncol(a1$counts[,ix_keep])),
      loess_span    = 0.01,
      n_pcs         = 'mcv',
      max_pcs       = 50,
      n_harmony     = 30,
      harmony_vars  = c("batch"),
      n_knn         = 50,
      leiden_res    = seq(0.5, 1.8, length.out = 10),
      leiden_iter   = 10,
      umap_spread   = 1,
      umap_min_dist = 0.25,
      log_file      = as.character(glue("{dirname(c1_file)}/analysis.log"))
    )
    c1_params_file <- glue("{dirname(c1_file)}/params.json")
    writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), c1_params_file) 
    #
    keep_cols <- colnames(a1$obs)
    keep_cols <- keep_cols[!str_detect(keep_cols, "^(PC|leiden|UMAP)")]
    c1 <- run_analysis(
      obs           = a1$obs[ix_keep,..keep_cols],
      counts        = a1$counts[,ix_keep],
      exclude_genes = c(mito_genes, bcr_genes),
      mito_genes    = mito_genes,
      params        = params
    )
    c1$params <- params
    c1$adt_clr <- a1$adt_clr[,ix_keep]
    print_status(glue("Writing {c1_file}"))
    qsave(c1, c1_file)
    print_status(glue("done"))
  } else {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
  }
  #
  source("R/colors-int.R")
  source("R/plot-analysis.R")
  # c1$obs$leiden <- c1$obs$leiden0.933
  my_leiden <- "leiden1.51"
  c1$obs$leiden <- c1$obs[[my_leiden]]
  out_dir <- glue("results/a20/{analysis_name}/figures_{my_leiden}")
  assign(analysis_name, c1)
  try({
  plot_analysis(
    analysis_name = analysis_name,
    out_dir       = out_dir,
    rowname_key   = ensembl_to_symbol,
    exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
    do_pb         = FALSE
  )
  })
  #
  source("R/functions/make-cellbrowser.R")
  make_cellbrowser(
    analysis_name     = analysis_name,
    out_dir           = glue("results/a20/{analysis_name}/figures_{my_leiden}"),
    cb_dir            = glue("cellbrowser/colitis/a20/{analysis_name}"),
    ensembl_to_symbol = ensembl_to_symbol
  )
  #
  if ("adt_clr" %in% names(c1)) {
    for (my_protein in c("CD8 ADT_C0046", "CD4 ADT_C0072")) {
        x <- left_join(
          x = c1$obs %>% select(cell, UMAP1, UMAP2),
          y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
          by = "cell"
        )
        x$protein[is.na(x$protein)] <- 0
        p <- plot_hexgene(
          x            = x$UMAP1,
          y            = x$UMAP2,
          z            = x$protein,
          # z            = x1$protein / x2$protein,
          bins         = 201,
          palette      = "davos",
          direction    = -1,
          use_quantile = FALSE,
          legend       = FALSE,
          italic       = FALSE
        ) +
        guides(
          fill = guide_colorbar(
            direction = "horizontal",
            title.position = "top",
            title = "CLR",
            barwidth = 15
          )
        ) +
        labs(subtitle = my_protein)
        my_ggsave(
          glue("umap-{safe(my_protein)}"),
          # glue("umap-CD8-to-CD4"),
          out_dir = glue("{out_dir}/adt-umap"),
          type = "pdf",
          plot = p,
          scale = 1, width = 4, height = 5, units = "in", dpi = 300
        )
    }
  }
}


# Add VDJ data to the existing objects
########################################################################
for (analysis_name in c("blood2_bcell4", "blood2_myeloid5", "blood2_tcell5")) {
  #
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    #
    if (!"has_tcr" %in% colnames(c1$obs)) {
      tcr_files <- Sys.glob(
        "/projects/irae_blood/cellranger_output/C*_tcr/all_contig_annotations.csv"
      )
      source("R/functions/load-tcr.R")
      tcr <- load_tcr(tcr_files)
      tcr$cell <- str_replace(tcr$cell, "_tcr", "_gex")
      tcr$cell <- str_remove(tcr$cell, "-1")
      tcr$channel <- str_remove(tcr$channel, "_tcr")
      #
      if (length(intersect(tcr$cell, c1$obs$cell)) > 1) {
        c1$obs <- dplyr::left_join(
          x = c1$obs,
          y = tcr %>% select(-channel),
          by = "cell"
        )
        tcr_cols <- c("TRBV", "TRAV", "TRB_cdr3", "TRA_cdr3")
        c1$obs$has_tcr <- with(c1$obs,
          !is.na(TRB_cdr3) & !is.na(TRA_cdr3)
        )
      }
    }
    if (!"has_bcr" %in% colnames(c1$obs)) {
      bcr_files <- Sys.glob(
        "/projects/irae_blood/cellranger_output/C*_bcr/all_contig_annotations.csv"
      )
      source("R/functions/load-bcr.R")
      bcr <- load_bcr(bcr_files)
      bcr$cell <- str_replace(bcr$cell, "_bcr", "_gex")
      bcr$cell <- str_remove(bcr$cell, "-1")
      bcr$channel <- str_remove(bcr$channel, "_bcr")
      #
      if (length(intersect(bcr$cell, c1$obs$cell)) > 1) {
        c1$obs <- dplyr::left_join(
          x = c1$obs,
          y = bcr %>% select(-channel, -barcode),
          by = "cell"
        )
        c1$obs$has_bcr <- with(c1$obs,
          !is.na(IGH_cdr3) & (
            !is.na(IGL_cdr3) | !is.na(IGK_cdr3)
          )
        )
      }
    }
    stopifnot(all(c1$obs$cell == colnames(c1$counts)))
    qsave(c1, c1_file)
  }
}

for (analysis_name in c("a12_4_4_t4_cd8_1_2", "a12_4_4_t4_cd4_2_2", "a12_4_4_b5_1_3")) {
  #
  t1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(t1_file)) {
    print_status(glue("Reading {t1_file}"))
    t1 <- qread(t1_file)
    print_status(glue("done"))
    #
    if (!"has_tcr" %in% colnames(t1$obs)) {
      tcr_files <- Sys.glob(
        "analysis/terra/cellranger-vdj_2020-07-27/output/*_TCR*/all_contig_annotations.csv"
      )
      source("R/functions/load-tcr.R")
      tcr <- load_tcr(tcr_files)
      tcr$cell <- str_replace(tcr$cell, "_TCR", "_GEX")
      tcr$channel <- str_remove(tcr$channel, "_tcr")
      #
      if (length(intersect(tcr$cell, t1$obs$cell)) > 1) {
        t1$obs <- dplyr::left_join(
          x = t1$obs,
          y = tcr %>% select(-channel),
          by = "cell"
        )
        tcr_cols <- c("TRBV", "TRAV", "TRB_cdr3", "TRA_cdr3")
        t1$obs$has_tcr <- with(t1$obs,
          !is.na(TRB_cdr3) & !is.na(TRA_cdr3)
        )
      }
    }
    if (!"has_bcr" %in% colnames(t1$obs)) {
      bcr_files <- Sys.glob(
        "analysis/terra/cellranger-vdj_2020-07-27/output/*_BCR*/all_contig_annotations.csv"
      )
      source("R/functions/load-bcr.R")
      bcr <- load_bcr(bcr_files)
      bcr$cell <- str_replace(bcr$cell, "_BCR", "_GEX")
      bcr$channel <- str_remove(bcr$channel, "_bcr")
      #
      if (length(intersect(bcr$cell, t1$obs$cell)) > 1) {
        t1$obs <- dplyr::left_join(
          x = t1$obs,
          y = bcr %>% select(-channel, -barcode),
          by = "cell"
        )
        t1$obs$has_bcr <- with(t1$obs,
          !is.na(IGH_cdr3) & (
            !is.na(IGL_cdr3) | !is.na(IGK_cdr3)
          )
        )
      }
    }
    stopifnot(all(t1$obs$cell == colnames(t1$counts)))
    qsave(t1, t1_file)
  }
}


# Add protein data to the existing objects
########################################################################
for (analysis_name in c("blood2_bcell4", "blood2_myeloid5", "blood2_tcell5")) {
  #
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    #
    out_dir <- glue("results/a20/{analysis_name}/figures")
    #
    if (!"adt_counts" %in% names(c1)) {
      antibodies <- read_excel(
        "data/biolegend/197 antibodies TOTAL-seqC_121819QG_99328_AntibodyList_Barcodes updated.xlsx"
      )
      antibodies$DNA_ID <- sprintf("ADT_C%04d", antibodies$DNA_ID)
      antibodies$id <- str_replace(str_replace(sprintf("%s %s",
        str_split_fixed(antibodies$Target, " ", 3)[,3],
        antibodies$DNA_ID
      ), "  ", " "), "^ ", "")
      ab_to_id <- unlist(split(antibodies$id, antibodies$DNA_ID))
      #
      adt_counts_file <-file.path(dirname(out_dir), "data", "adt_counts.qs")
      if (!file.exists(adt_counts_file)) {
        adt_files <- Sys.glob(
          glue(
            "/projects/irae_blood/cellranger_output/C*_fbc/*_fbc.csv"
          )
        )
        adt_batches <- basename(dirname(adt_files))
        source("R/functions/read-sparse-csv.R")
        read_adt <- function(adt_file, rna_cell_ids) {
          adt <- read_sparse_csv(adt_file)
          adt <- adt[!stringr::str_detect(rownames(adt), "^HT_"),]
          adt <- adt[!str_detect(rownames(adt), "^HT"),]
          adt <- adt[,Matrix::colSums(adt) > 0]
          batch <- str_replace(basename(dirname(adt_file)), "_fbc", "_gex")
          colnames(adt) <- sprintf("%s|%s", batch, colnames(adt))
          ix <- colnames(adt) %in% rna_cell_ids
          if (sum(ix) > 0) {
            adt <- adt[,ix]
            rownames(adt) <- ab_to_id[rownames(adt)]
            return(adt)
          }
          return(NA)
        }
        adt_list <- pblapply(adt_files, function(adt_file) {
          read_adt(adt_file, colnames(c1$counts))
        })
        adt_counts <- do.call(cbind, adt_list[!is.na(adt_list)])
        rm(adt_list)
        #
        qsave(adt_counts, adt_counts_file)
      } else {
        adt_counts <- qread(adt_counts_file)
      }
      #
      normalize_clr <- function(counts) {
        # cells are columns
        denom <- exp(colSums(log1p(counts)) / ncol(counts))
        counts@x <- log1p(counts@x / rep(denom, diff(counts@p)))
        return(counts)
      }
      #
      adt_clr <- normalize_clr(adt_counts)
      both_cell_ids <- intersect(colnames(adt_clr), colnames(c1$counts))
      length(both_cell_ids)
      length(colnames(c1$counts))
      length(colnames(adt_clr))
      #
      c1$adt_counts <- Matrix(
        nrow = nrow(adt_counts),
        ncol = ncol(c1$counts)
      )
      colnames(c1$adt_counts) <- colnames(c1$counts)
      rownames(c1$adt_counts) <- rownames(adt_counts)
      #
      ix <- which(colnames(c1$adt_counts) %in% colnames(adt_counts))
      for (i in seq_along(rownames(c1$adt_counts))) {
        c1$adt_counts[i,ix] <- adt_counts[i,]
      }
      c1$adt_counts <- as(c1$adt_counts, "dgCMatrix")
      c1$adt_clr <- normalize_clr(c1$adt_counts)
      #
      qsave(c1, c1_file)
    }
  }
}

# Plot proteins
########################################################################
for (analysis_name in c("blood2_bcell4", "blood2_myeloid5", "blood2_tcell5")) {
  #
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    #
    out_dir <- glue("results/a20/{analysis_name}/figures")
    #
    if ("adt_clr" %in% names(c1)) {
      #
      for (my_protein in rownames(c1$adt_clr)) {
        x <- left_join(
          x = c1$obs %>% select(cell, UMAP1, UMAP2),
          y = data.frame(cell = colnames(c1$adt_clr), protein = c1$adt_clr[my_protein,]),
          by = "cell"
        )
        x$protein[is.na(x$protein)] <- 0
        p <- plot_hexgene(
          x            = x$UMAP1,
          y            = x$UMAP2,
          z            = x$protein,
          bins         = 201,
          palette      = "davos",
          direction    = -1,
          use_quantile = FALSE,
          legend       = FALSE,
          italic       = FALSE
        ) +
        guides(
          fill = guide_colorbar(
            direction = "horizontal",
            title.position = "top",
            title = "CLR",
            barwidth = 15
          )
        ) +
        labs(subtitle = my_protein)
        my_ggsave(
          glue("umap-{safe(my_protein)}"),
          out_dir = glue("{out_dir}/adt-umap"),
          type = "pdf",
          plot = p,
          scale = 1, width = 4, height = 5, units = "in", dpi = 300
        )
      }
      c1$adt_clr@x[is.na(c1$adt_clr@x)] <- 0
      x <- cbind.data.frame(
        p1 = as.numeric(c1$adt_clr["CD8 ADT_C0046",]),
        p2 = as.numeric(c1$adt_clr["CD4 ADT_C0072",])
      )
      p <- ggplot(x) +
        aes(p1, p2) +
        stat_binhex(bins = 200) +
        scale_fill_gradientn(
          colors = scico::scico(20)[1:20],
          breaks = log_breaks(5), trans = "log10"
        ) +
        labs(
          x = "CD8", y = "CD4"
        ) +
        guides(fill = guide_colorbar(title = "Cells", barheight = 7))
      my_ggsave(
        glue("CD8-CD4"),
        out_dir = glue("{out_dir}/adt-pair"),
        type = "pdf",
        plot = p,
        scale = 1, width = 6, height = 5, units = "in", dpi = 300
      )
    }
  }
}



      # does not work on server
      # make_montage <- function(files, out, ncol = 10) {
      #   retval <- split(files, (seq_along(files) - 1) %/% ncol) %>%
      #     lapply(function(x) {
      #       lapply(x, magick::image_read) %>%
      #       {Reduce(c, .)} %>%
      #       magick::image_scale("x450") %>%
      #   #     image_annotate(
      #   #       text = x %>% basename,
      #   #       gravity = "north", location = "-20+40",
      #   #       size = 12, boxcolor = "white"
      #   #     ) %>%
      #       magick::image_append() #%>%
      #   #     image_trim()
      #     }) %>%
      #     {Reduce(c, .)} %>%
      #     magick::image_append(stack = TRUE) %>%
      #     magick::image_background("white", flatten = TRUE)
      #   magick::image_write(retval, path = out, format = "png", density = 300)
      # }
      # pdf_files <- Sys.glob(glue("{out_dir}/adt-umap/*.pdf"))
      # png_files <- str_replace(pdf_files, "pdf$", "png")
      # for (i in seq_along(pdf_files)) {
      #   animation::im.convert(pdf_files[i], output = png_files[i], extra.opts = "-density 150")
      # }
      # png_montage <- glue("{out_dir}/adt-umap.png")
      # make_montage(pdf_files, png_montage)


for (analysis_name in c("blood2_bcell4", "blood2_myeloid5", "blood2_tcell5")) {
  c1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  if (file.exists(c1_file)) {
    print_status(glue("Reading {c1_file}"))
    c1 <- qread(c1_file)
    print_status(glue("done"))
    my_leiden <- "leiden1.51"
    if (analysis_name == "blood2_myeloid5") {
      my_leiden <- "leiden0.644"
    } else if (analysis_name == "blood2_tcell5") {
      my_leiden <- "leiden1.08"
    }
    c1$obs$leiden <- c1$obs[[my_leiden]]
    assign(analysis_name, c1)
    source("R/functions/make-cellbrowser.R")
    make_cellbrowser(
      analysis_name     = analysis_name,
      out_dir           = glue("results/a20/{analysis_name}/figures_{my_leiden}"),
      cb_dir            = glue("cellbrowser/colitis/a20/{analysis_name}"),
      ensembl_to_symbol = ensembl_to_symbol
    )
  }
}

# }}}

# TCR sharing: Can we find the same TCR sequences in tissue and blood? {{{

t1_analysis_name <- "a12_4_4_t4_cd8_1_2"
t1_file <- glue("results/a20/{t1_analysis_name}/data/{t1_analysis_name}.qs")
file.exists(t1_file)
t1 <- qread(t1_file)
table(t1$obs$leiden)

b1_analysis_name <- "blood2_tcell5_cd8_5"
b1_file <- glue("results/a20/{b1_analysis_name}/data/{b1_analysis_name}.qs")
file.exists(b1_file)
b1 <- qread(b1_file)

"has_tcr" %in% colnames(t1$obs)
"has_tcr" %in% colnames(b1$obs)

t1$obs$tcr_seq <- with(t1$obs, sprintf("%s %s", TRA_cdr3_trim, TRB_cdr3_trim))
t1_tcrs <- unique(t1$obs$tcr_seq)
t1_tcrs <- t1_tcrs[t1_tcrs != "NA NA"]

b1$obs$tcr_seq <- with(b1$obs, sprintf("%s %s", TRA_cdr3_trim, TRB_cdr3_trim))
b1_tcrs <- unique(b1$obs$tcr_seq)
b1_tcrs <- b1_tcrs[b1_tcrs != "NA NA"]

length(t1_tcrs)
length(b1_tcrs)

both_tcrs <- intersect(t1_tcrs, b1_tcrs)
length(both_tcrs)

t1$obs$tcr_shared <- t1$obs$tcr_seq %in% both_tcrs
b1$obs$tcr_shared <- b1$obs$tcr_seq %in% both_tcrs

p <- plot_hexgene(
  x            = t1$obs$UMAP1,
  y            = t1$obs$UMAP2,
  z            = 100 * t1$obs$tcr_shared,
  # bins         = 201,
  # palette      = "grayC",
  bins         = 91,
  palette      = "bilbao",
  direction    = 1,
  use_quantile = FALSE,
  legend       = FALSE,
  italic       = FALSE
) +
scale_fill_gradientn(
  colors = scico::scico(palette = "grayC", direction = 1, n = 19)[2:19], trans = "log1p",
  breaks = c(0, 10, 50, 100)
) +
guides(
  fill = guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    title = "Percent of cells with TCR",
    barwidth = 15
  )
) +
labs(subtitle = "Tissue cells with TCRs seen in blood")
my_ggsave(
  "umap-tcr-tissue-blood",
  out_dir = file.path(dirname(dirname(b1_file)), "figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 4, height = 5, units = "in", dpi = 300
)

p <- plot_hexgene(
  x            = b1$obs$UMAP1,
  y            = b1$obs$UMAP2,
  z            = 100 * b1$obs$tcr_shared,
  # bins         = 201,
  # palette      = "grayC",
  bins         = 91,
  palette      = "bilbao",
  direction    = 1,
  use_quantile = FALSE,
  legend       = FALSE,
  italic       = FALSE
) +
scale_fill_gradientn(
  colors = scico::scico(palette = "grayC", direction = 1, n = 19)[2:19], trans = "log1p",
  breaks = c(0, 10, 50, 100)
) +
guides(
  fill = guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    title = "Percent of cells with TCR",
    barwidth = 15
  )
) +
labs(subtitle = "Blood cells with TCRs seen in tissue")
my_ggsave(
  "umap-tcr-blood-tissue",
  out_dir = file.path(dirname(dirname(b1_file)), "figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 4, height = 5, units = "in", dpi = 300
)

# per donor
for (donor in unique(t1$obs$donor)) {
  t_ix <- t1$obs$donor == donor
  t_z <- (100 * t1$obs$tcr_seq %in% both_tcrs)[t_ix]
  if (sum(t_z) > 0) {
    p <- plot_hexgene(
      x            = t1$obs$UMAP1[t_ix],
      y            = t1$obs$UMAP2[t_ix],
      z            = t_z,
      bins         = 201,
      palette      = "davos",
      direction    = -1,
      use_quantile = FALSE,
      legend       = FALSE,
      italic       = FALSE
    ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title = glue("{donor} Percent of cells with TCR"),
        barwidth = 15
      )
    ) +
    labs(subtitle = "TCRs seen in blood and tissue")
    my_ggsave(
      glue("umap-tcr-blood-tissue-{donor}"),
      out_dir = file.path(dirname(dirname(t1_file)), "figures/tcr/by-donor"),
      type = "pdf",
      plot = p,
      scale = 1, width = 4, height = 5, units = "in", dpi = 300
    )
  }
  b_ix <- b1$obs$donor == donor
  b_z <- (100 * b1$obs$tcr_seq %in% both_tcrs)[b_ix]
  if (sum(t_z) > 0) {
    p <- plot_hexgene(
      x            = b1$obs$UMAP1[b_ix],
      y            = b1$obs$UMAP2[b_ix],
      z            = b_z,
      bins         = 201,
      palette      = "davos",
      direction    = -1,
      use_quantile = FALSE,
      legend       = FALSE,
      italic       = FALSE
    ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title = glue("{donor} Percent of cells with TCR"),
        barwidth = 15
      )
    ) +
    labs(subtitle = "TCRs seen in blood and tissue")
    my_ggsave(
      glue("umap-tcr-blood-tissue-{donor}"),
      out_dir = file.path(dirname(dirname(b1_file)), "figures/tcr/by-donor"),
      type = "pdf",
      plot = p,
      scale = 1, width = 4, height = 5, units = "in", dpi = 300
    )
  }
}

t1$obs %>% filter(tcr_seq %in% both_tcrs) %>% count(case)
b1$obs %>% filter(tcr_seq %in% both_tcrs) %>% count(case)

x <- rbind(
  t1$obs %>% filter(tcr_seq %in% both_tcrs) %>%
    count(case, drug, donor, tcr_seq) %>%
    mutate(data = "tissue"),
  b1$obs %>% filter(tcr_seq %in% both_tcrs) %>%
    count(case, drug, donor, tcr_seq) %>%
    mutate(data = "blood")
)
x <- x %>% pivot_wider(names_from = "data", values_from = "n", values_fill = 0)
x %>% filter(tissue > 2 & blood > 2) %>% arrange(-tissue - blood) %>% head

x %>% filter(tcr_seq == "NGGYQKVT SLEAGGNEKLF")
x %>% filter(tcr_seq == "NRDDKII SPPVASATEAF")


# Percent of TCRs shared between blood and tissue
n_shared <- x %>% group_by(case, drug, donor) %>%
  summarize(
    n_shared_tcrs = sum(tissue >= 1 & blood >= 1),
    .groups = "drop"
  )
n_shared <- left_join(
  x = n_shared,
  y = t1$obs %>% group_by(case, drug, donor) %>% summarize(n_tissue_tcrs = length(unique(tcr_seq))),
  by = c("case", "drug", "donor")
)
n_shared <- left_join(
  x = n_shared,
  y = b1$obs %>% group_by(case, drug, donor) %>% summarize(n_blood_tcrs = length(unique(tcr_seq))),
  by = c("case", "drug", "donor")
)
n_shared <- n_shared %>% mutate(
  pct_shared_tcrs = 100 * n_shared_tcrs / (n_tissue_tcrs + n_blood_tcrs),
  total_tcrs = n_tissue_tcrs + n_blood_tcrs
)
n_shared_t <-t.test(formula = pct_shared_tcrs ~ case, data = n_shared)
n_shared_pval <- n_shared_t$p.value
p <- ggplot(n_shared) +
  geom_point(
    mapping = aes(y = case, x = pct_shared_tcrs, fill = case),
    shape = 21, size = 4, stroke = 0.1,
    position = position_quasirandom(width = 0.5, groupOnX = FALSE)
  ) +
  scale_fill_manual(values = pals::okabe()[2:1], guide = "none") +
  labs(
    y = NULL, x = "Percent",
    title = "Percent of TCRs shared\nbetween blood and tissue",
    subtitle = glue("P = {signif(n_shared_pval, 2)} (t-test)")
  )
my_ggsave(
  "percent-shared-blood-tissue",
  out_dir = glue("results/a20/{b1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 4, height = 2, units = "in", dpi = 300
)
n_shared_file <- glue("results/a20/{b1_analysis_name}/figures/tcr/n_shared_tcrs.tsv")
message(glue("Writing {n_shared_file}"))
fwrite(n_shared, n_shared_file, sep = "\t")


# per cluster
x <- rbind(
  t1$obs %>% filter(tcr_seq %in% both_tcrs) %>%
    count(leiden, case, donor, tcr_seq) %>%
    mutate(leiden = sprintf("t%s", leiden)),
  b1$obs %>% filter(tcr_seq %in% both_tcrs) %>%
    count(leiden, case, donor, tcr_seq) %>%
    mutate(leiden = sprintf("b%s", leiden))
)
x %>% filter(tcr_seq == "KQTDKLI SFGPGIQPQH")

x <- x %>% pivot_wider(names_from = "leiden", values_from = "n", values_fill = 0)
t_cols <- colnames(x)[4:15]
b_cols <- colnames(x)[16:31]

rowSums(x[,t_cols])
rowSums(x[,b_cols])

t1_hyper <- rbindlist(lapply(unique(t1$obs$donor), function(my_donor) {
  rbindlist(lapply(unique(t1$obs$leiden), function(my_cluster) {
    # my_donor <- "SIC_141"
    # my_cluster <- 1
    n_total <- t1$obs %>% filter(donor == my_donor) %>% nrow
    n_up <- t1$obs %>% filter(donor == my_donor, tcr_shared) %>% nrow
    n_drawn <- t1$obs %>% filter(donor == my_donor, leiden == my_cluster) %>% nrow
    n_overlap <- t1$obs %>% filter(donor == my_donor, leiden == my_cluster, tcr_shared) %>% nrow
    my_fold <- (n_overlap / (n_drawn + 1)) / (n_up / (n_total + 1))
    if (n_up == 0) {
      my_fold <- 1
    }
    my_pval <- phyper(n_overlap - 1, n_up, n_total - n_up, n_drawn, lower.tail = FALSE)
    list(
      donor = my_donor, case = t1$obs$case[t1$obs$donor == my_donor][1], leiden = my_cluster,
      n_total = n_total, n_up = n_up, n_drawn = n_drawn, n_overlap = n_overlap,
      fold = my_fold, pval = my_pval
    )
  }))
})) %>% filter(n_up > 0)

p <- ggplot(t1_hyper) +
  aes(x = log2(fold), y = case, group = donor, fill = case) +
  scale_x_continuous(labels = function(x) 2^x) +
  geom_vline(xintercept = 0, size = 0.3) +
  geom_point(
    shape = 21, size = 3, stroke = 0.1,
    # position = position_quasirandom(groupOnX = FALSE)
    position = position_jitter()
  ) +
  scale_fill_manual(values = pals::okabe()[2:1], guide = "none") +
  facet_wrap(vars(leiden)) +
  labs(x = "Fold", y = NULL) +
  theme(panel.spacing = unit(0.5, "lines"))
my_ggsave(
  "enrichment-of-blood-tcr-intissue",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 9, height = 4, units = "in", dpi = 300
)

# Enrichment of blood TCRs in tissue cells
########################################################################
{
  source("R/functions/composition.R")
  t1$obs$cluster <- t1$obs$leiden
  t1_masc_file <- glue("results/a20/{t1_analysis_name}/figures/tcr/masc.qs")
  if (!file.exists(t1_masc_file)) {
    t1_masc <- do_masc(
      t1$obs %>% filter(has_tcr),
      form1 = "is_cluster ~ 1 + (1|donor)",
      form2 = "is_cluster ~ 1 + tcr_shared + (1|donor)",
      mc.cores = 8
    )
    qsave(t1_masc, t1_masc_file)
  } else {
    t1_masc <- qread(t1_masc_file)
  }
  t1_coef <- summarize_masc(t1_masc)
  t1_coef <- t1_coef %>% group_by(value) %>% mutate(lrt_fdr = p.adjust(lrt_p, method = "fdr"))
  fwrite(t1_coef, glue("results/a20/{t1_analysis_name}/figures/tcr/masc-tcr-sharing.tsv"), sep = "\t")
  # conflict_prefer("geom_errorbarh", "ggplot2")
  #
  d_error <- t1_coef %>% filter(value == "tcr_sharedTRUE") %>% arrange(lrt_p)
  d_error <- left_join(
    x = d_error,
    y = t1$obs %>% filter(tcr_shared) %>% group_by(cluster) %>% summarize(n_donors = length(unique(donor))),
    by = "cluster"
  )
  # d_error %>% select(est, lrt_p, lrt_fdr, cluster, n_donors)
  ix_toosmall <- which(d_error$n_donors < 3)
  d_error$lrt_p[ix_toosmall] <- 1
  d_error$lrt_fdr[ix_toosmall] <- 1
  d_error$est[ix_toosmall] <- NA
  d_error$est_low[ix_toosmall] <- NA
  d_error$est_high[ix_toosmall] <- NA
  d_error <- d_error %>% arrange(lrt_p)
  #
  d_error$cluster <- factor(d_error$cluster, rev(d_error$cluster))
  d_error <- d_error %>%
    mutate(
      est      = log2(exp(est)),
      est_low  = log2(exp(est_low)),
      est_high = log2(exp(est_high))
    )
  vlines <- seq(-8, 8, by = 1)
  vlines <- vlines[vlines > min(d_error$est_low, na.rm = TRUE)]
  vlines <- vlines[vlines < max(d_error$est_high, na.rm = TRUE)]
  vlines <- vlines[vlines != 0]
  p_error <- ggplot() +
    ggforestplot::geom_stripes(data = d_error, aes(y = cluster)) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_vline(xintercept = vlines, size = 0.3, color = "white") +
    geom_errorbarh(
      data = d_error,
      mapping = aes(y = cluster, xmin = est_low, xmax = est_high, color = lrt_fdr < 0.05),
      height = 0
    ) +
    geom_point(
      data = d_error, aes(x = est, y = cluster, color = lrt_fdr < 0.05)
    ) +
    geom_text(
      data = d_error %>% filter(lrt_fdr < 0.05),
      mapping = aes(
        x = ifelse(est > 0, -Inf, Inf),
        hjust = ifelse(est > 0, 0, 1),
        y = cluster,
        label = sprintf(" %s ", signif(lrt_p, 1))
      ),
      size = 5, color = "grey30"
    ) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = seq(-8, 8, by = 1),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    labs(title = "TCR present in blood", y = NULL, x = "OR") +
    theme(plot.background = element_blank())
  my_ggsave(
    "errorbars",
    out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
    type = "pdf",
    plot = p_error,
    scale = 1, width = 5, height = 4, units = "in", dpi = 300
  )
  # d_error %>% select(cluster, est, lrt_p, lrt_fdr) %>%
  #   mutate(est = 2^est) %>% mutate_if(is.numeric, signif, 2)
  source("R/plot-composition.R")
  my_comp <- t1$obs %>%
    dplyr::filter(has_tcr) %>%
    dplyr::select(tcr_shared, case, donor, leiden) %>%
    dplyr::group_by(tcr_shared, case, donor, leiden) %>%
    dplyr::summarize(n = n(), .groups = "drop_last") %>%
    dplyr::mutate(freq = 100 * n / sum(n))
  # my_comp %>% group_by(tcr_shared, donor) %>%
  #   summarize(sum = sum(freq))
  my_comp$cluster <- factor(as.character(my_comp$leiden), levels(d_error$cluster))
  #
  p_boxplot <- plot_composition_h(
    d = t1$obs, fill = case, group = donor, x = cluster,
    comp = my_comp,
    legend.position = "right"
  ) +
  facet_wrap(~ tcr_shared) +
  scale_fill_manual(values = pals::okabe()[2:1]) +
  labs(title = glue("Composition of each donor (n = {length(unique(t1$obs$donor))})"))
  my_ggsave(
    "compositionh-case",
    out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
    # type = c("pdf", "png"),
    type = "pdf",
    plot = p_boxplot,
    scale = 1, height = 9, width = 8, units = "in", dpi = 300
  )
  #
  p_boxplot <- plot_composition_h(
    d = t1$obs, fill = tcr_shared, group = donor, x = cluster,
    comp = my_comp,
    legend.position = "right"
  ) +
  scale_fill_manual(
    values = pals::okabe()[c(6,8)],
    labels = c("FALSE" = "Tissue-only", "TRUE" = "Blood and Tissue"),
    name = "TCR"
  ) +
  theme(legend.position = "bottom") +
  labs(title = glue("Composition of each donor (n = {length(unique(t1$obs$donor))})"))
  my_ggsave(
    "compositionh",
    out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
    # type = c("pdf", "png"),
    type = "pdf",
    plot = p_boxplot,
    scale = 1, height = 8, width = 5, units = "in", dpi = 300
  )
  #
  p_both <- (
    p_error + 
      labs(title = "Tissue cells with TCR seen in blood", y = NULL, x = "OR")
  ) + (
    p_boxplot +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
  p_both <- p_both + plot_layout(widths = c(1, 1.5))
  my_ggsave(
    "lanes",
    out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
    type = "pdf",
    plot = p_both,
    scale = 1,
    height = length(unique(d_error$cluster)) * 0.3 + 1,
    width = 6,
    units = "in", dpi = 300
  )
}

# Enrichment of tissue TCRs in blood cells
########################################################################
{
  source("R/functions/composition.R")
  b1$obs$cluster <- b1$obs$leiden
  b1_masc_file <- glue("results/a20/{b1_analysis_name}/figures/tcr/masc.qs")
  if (!file.exists(b1_masc_file)) {
    b1_masc <- do_masc(
      b1$obs %>% filter(has_tcr),
      form1 = "is_cluster ~ 1 + (1|donor)",
      form2 = "is_cluster ~ 1 + tcr_shared + (1|donor)",
      mc.cores = 8
    )
    qsave(b1_masc, b1_masc_file)
  } else {
    b1_masc <- qread(b1_masc_file)
  }
  b1_coef <- summarize_masc(b1_masc)
  b1_coef <- b1_coef %>% group_by(value) %>% mutate(lrt_fdr = p.adjust(lrt_p, method = "fdr"))
  fwrite(b1_coef, glue("results/a20/{b1_analysis_name}/figures/tcr/masc-tcr-sharing.tsv"), sep = "\t")
  #
  d_error <- b1_coef %>% filter(value == "tcr_sharedTRUE") %>% arrange(lrt_p)
  d_error <- left_join(
    x = d_error,
    y = b1$obs %>% filter(tcr_shared) %>% group_by(cluster) %>% summarize(n_donors = length(unique(donor))),
    by = "cluster"
  )
  # d_error %>% select(est, lrt_p, lrt_fdr, cluster, n_donors)
  ix_toosmall <- which(d_error$n_donors < 3)
  d_error$lrt_p[ix_toosmall] <- 1
  d_error$lrt_fdr[ix_toosmall] <- 1
  d_error$est[ix_toosmall] <- NA
  d_error$est_low[ix_toosmall] <- NA
  d_error$est_high[ix_toosmall] <- NA
  d_error <- d_error %>% arrange(lrt_p)
  #
  d_error$cluster <- factor(d_error$cluster, rev(d_error$cluster))
  d_error <- d_error %>%
    mutate(
      est      = log2(exp(est)),
      est_low  = log2(exp(est_low)),
      est_high = log2(exp(est_high))
    )
  vlines <- seq(-8, 8, by = 1)
  vlines <- vlines[vlines > min(d_error$est_low)]
  vlines <- vlines[vlines < max(d_error$est_high)]
  vlines <- vlines[vlines != 0]
  p_error <- ggplot() +
    ggforestplot::geom_stripes(data = d_error, aes(y = cluster)) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_vline(xintercept = vlines, size = 0.3, color = "white") +
    ggplot2::geom_errorbarh(
      data = d_error,
      mapping = aes(y = cluster, xmin = est_low, xmax = est_high, color = lrt_fdr < 0.05),
      height = 0
    ) +
    geom_point(
      data = d_error, aes(x = est, y = cluster, color = lrt_fdr < 0.05)
    ) +
    geom_text(
      data = d_error %>% filter(lrt_fdr < 0.05),
      mapping = aes(
        x = ifelse(est > 0, -Inf, Inf),
        hjust = ifelse(est > 0, 0, 1),
        y = cluster,
        label = sprintf(" %s ", signif(lrt_p, 1))
      ),
      size = 5, color = "grey30"
    ) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = seq(-8, 8, by = 1),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    labs(title = "TCR present in blood", y = NULL, x = "OR") +
    theme(plot.background = element_blank())
  my_ggsave(
    "errorbars",
    out_dir = glue("results/a20/{b1_analysis_name}/figures/tcr"),
    type = "pdf",
    plot = p_error,
    scale = 1, width = 5, height = 4, units = "in", dpi = 300
  )
  # d_error %>% select(cluster, est, lrt_p, lrt_fdr) %>%
  #   mutate(est = 2^est) %>% mutate_if(is.numeric, signif, 2)
  source("R/plot-composition.R")
  my_comp <- b1$obs %>%
    dplyr::filter(has_tcr) %>%
    dplyr::select(tcr_shared, case, donor, leiden) %>%
    dplyr::group_by(tcr_shared, case, donor, leiden) %>%
    dplyr::summarize(n = n(), .groups = "drop_last") %>%
    dplyr::mutate(freq = 100 * n / sum(n))
  # my_comp %>% group_by(tcr_shared, donor) %>%
  #   summarize(sum = sum(freq))
  my_comp$cluster <- factor(as.character(my_comp$leiden), levels(d_error$cluster))
  #
  p_boxplot <- plot_composition_h(
    d = b1$obs, fill = case, group = donor, x = cluster,
    comp = my_comp,
    legend.position = "right"
  ) +
  facet_wrap(~ tcr_shared) +
  scale_fill_manual(values = pals::okabe()[2:1]) +
  labs(title = glue("Composition of each donor (n = {length(unique(b1$obs$donor))})"))
  my_ggsave(
    "compositionh-case",
    out_dir = glue("results/a20/{b1_analysis_name}/figures/tcr"),
    # type = c("pdf", "png"),
    type = "pdf",
    plot = p_boxplot,
    scale = 1, height = 9, width = 8, units = "in", dpi = 300
  )
  #
  p_boxplot <- plot_composition_h(
    d = b1$obs, fill = tcr_shared, group = donor, x = cluster,
    comp = my_comp,
    legend.position = "right"
  ) +
  scale_fill_manual(
    values = pals::okabe()[c(6,8)],
    labels = c("FALSE" = "Blood-only", "TRUE" = "Blood and Tissue"),
    name = "TCR"
  ) +
  theme(legend.position = "bottom") +
  labs(title = glue("Composition of each donor (n = {length(unique(b1$obs$donor))})"))
  my_ggsave(
    "compositionh",
    out_dir = glue("results/a20/{b1_analysis_name}/figures/tcr"),
    # type = c("pdf", "png"),
    type = "pdf",
    plot = p_boxplot,
    scale = 1, height = 8, width = 5, units = "in", dpi = 300
  )
  #
  p_both <- (
    p_error + 
      labs(title = "Blood cells with TCR seen in tissue", y = NULL, x = "OR")
  ) + (
    p_boxplot +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
  p_both <- p_both + plot_layout(widths = c(1, 1.5))
  my_ggsave(
    "lanes",
    out_dir = glue("results/a20/{b1_analysis_name}/figures/tcr"),
    type = "pdf",
    plot = p_both,
    scale = 1,
    height = length(unique(d_error$cluster)) * 0.3 + 1,
    width = 6,
    units = "in", dpi = 300
  )
}


# Do some TCRs appear in more clusters than we would expect by random chance?
# e.g. Do we see more TCRs in 5 clusters than expected?


# Do we see some cluster-networks more often than expected? (subgraph motifs in a graph)
# e.g. Do we see more TCRs with cluster-network (1,2,4,6) than expected?

t1$obs <- t1$obs %>%
  dplyr::mutate(tcr_id = paste(TRAV, TRBV, TRA_cdr3_trim, TRB_cdr3_trim)) %>%
  dplyr::mutate(tcr_id = ifelse(has_tcr, tcr_id, NA))

# Count the number of TCRs observed in each combination of clusters
# single clusters: 1, 10, 11, 12
# pairs of clusters: 1:2, 1:3
# three clusters: c(1, 2, 10)
# and so on...
ix <- t1$obs$has_tcr
my_clusters <- as.factor(t1$obs$leiden)
my_tcrs <- as.factor(t1$obs$tcr_id)
m <- sparseMatrix(
  i = as.integer(my_tcrs)[ix],
  j = as.integer(my_clusters)[ix]
)
retval <- data.frame(
  clusters = as.character(apply(m, 1, function(x) { which(x) }))
) %>% count(clusters) %>% mutate(pct = n / sum(n))

m_shuf <- rbindlist(lapply(1:1000, function(i) {
  shuf_clusters <- sample(my_clusters)
  m <- sparseMatrix(
    i = as.integer(my_tcrs)[ix],
    j = as.integer(shuf_clusters)[ix]
  )
  retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
  colnames(retval) <- c("c1", "c2")
  retval$shared <- apply(retval, 1, function(x) {
    sum(m[,x[1]] & m[,x[2]]) /
    sum(m[,x[1]] | m[,x[2]])
  })
  retval$iter <- i
  retval
}))

retval <- cbind(retval, rbindlist(apply(retval, 1, function(x) {
  val <- x[3]
  ix <- with(m_shuf, c1 == x[1] & c2 == x[2])
  num <- sum(m_shuf$shared[ix] > val)
  pval <- (num + 1) / (sum(ix) + 1)
  data.frame(greater = num, pval = pval, mean = mean(m_shuf$shared[ix]))
})))

retval %>% arrange(pval) %>% head


# Do some cluster pairs share more TCRs than we would expect by random chance?

t1$obs <- t1$obs %>%
  dplyr::mutate(tcr_id = paste(TRAV, TRBV, TRA_cdr3_trim, TRB_cdr3_trim)) %>%
  dplyr::mutate(tcr_id = ifelse(has_tcr, tcr_id, NA))

ix <- t1$obs$has_tcr
my_clusters <- as.factor(t1$obs$leiden)
my_tcrs <- as.factor(t1$obs$tcr_id)
m <- sparseMatrix(
  i = as.integer(my_tcrs)[ix],
  j = as.integer(my_clusters)[ix]
)
retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
colnames(retval) <- c("c1", "c2")
retval$shared <- apply(retval, 1, function(x) {
  sum(m[,x[1]] & m[,x[2]]) /
  sum(m[,x[1]] | m[,x[2]])
})

m_shuf <- rbindlist(lapply(1:1000, function(i) {
  shuf_clusters <- sample(my_clusters)
  m <- sparseMatrix(
    i = as.integer(my_tcrs)[ix],
    j = as.integer(shuf_clusters)[ix]
  )
  retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
  colnames(retval) <- c("c1", "c2")
  retval$shared <- apply(retval, 1, function(x) {
    sum(m[,x[1]] & m[,x[2]]) /
    sum(m[,x[1]] | m[,x[2]])
  })
  retval$iter <- i
  retval
}))

retval <- cbind(retval, rbindlist(apply(retval, 1, function(x) {
  val <- x[3]
  ix <- with(m_shuf, c1 == x[1] & c2 == x[2])
  num <- sum(m_shuf$shared[ix] > val)
  pval <- (num + 1) / (sum(ix) + 1)
  data.frame(greater = num, pval = pval, mean = mean(m_shuf$shared[ix]))
})))

retval %>% arrange(pval) %>% head

test_sharing <- function(cluster_ids, tcr_ids, n_iter = 1000) {
  stopifnot(
    length(cluster_ids) == length(tcr_ids)
  )
  #
  ix          <- !is.na(tcr_ids)
  my_clusters <- as.factor(cluster_ids)
  my_tcrs     <- as.factor(tcr_ids)
  #
  m <- Matrix::sparseMatrix(
    i = as.integer(my_tcrs)[ix],
    j = as.integer(my_clusters)[ix]
  )
  # All combinations of cluster ids
  retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
  colnames(retval) <- c("c1", "c2")
  # Number of shared TCRs for each pair of cluster ids
  # intersection / union
  retval$shared <- apply(retval, 1, function(x) {
    sum(m[,x[1]] & m[,x[2]]) /
    sum(m[,x[1]] | m[,x[2]])
  })
  # Shuffle the clusters and count again
  m_shuf <- data.table::rbindlist(pbapply::pblapply(seq(n_iter), function(i) {
    shuf_clusters <- sample(my_clusters)
    m <- Matrix::sparseMatrix(
      i = as.integer(my_tcrs)[ix],
      j = as.integer(shuf_clusters)[ix]
    )
    retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
    colnames(retval) <- c("c1", "c2")
    retval$shared <- apply(retval, 1, function(x) {
      sum(m[,x[1]] & m[,x[2]]) /
      sum(m[,x[1]] | m[,x[2]])
    })
    retval$iter <- i
    retval
  }))
  # Empirical p-values
  retval <- cbind(retval, rbindlist(apply(retval, 1, function(x) {
    val <- x[3]
    ix <- with(m_shuf, c1 == x[1] & c2 == x[2])
    gte <- sum(m_shuf$shared[ix] >= val)
    lte <- sum(m_shuf$shared[ix] <= val)
    data.frame(
      lte = lte,
      lte_pval = (lte + 1) / (sum(ix) + 1),
      gte = gte,
      gte_pval = (gte + 1) / (sum(ix) + 1),
      mean = mean(m_shuf$shared[ix])
    )
  })))
  list(
    shuf = m_shuf,
    pval = retval %>% arrange(gte_pval)
  )
}

t1_share <- test_sharing(
  cluster_ids = t1$obs$leiden,
  tcr_ids     = t1$obs$tcr_id,
  n_iter      = 1000
)

# Figure 2H {{{
t1_share2_file <- glue("results/a20/{t1_analysis_name}/figures/tcr/tcr-sharing.qs")
if (!file.exists(t1_share2_file)) {
  # elapsed 365 seconds
  system.time({
  t1_share2 <- test_sharing2(
    donor_ids   = t1$obs$donor,
    cluster_ids = t1$obs$leiden,
    tcr_ids     = t1$obs$tcr_id,
    n_iter      = 10000
  )
  })
  qsave(t1_share2, t1_share2_file)
} else {
  t1_share2 <- qread(t1_share2_file)
}

# Shuffle cluster labels within each donor
test_sharing2 <- function(donor_ids, cluster_ids, tcr_ids, n_iter = 1000) {
  stopifnot(
    length(donor_ids) == length(tcr_ids)
  )
  stopifnot(
    length(cluster_ids) == length(tcr_ids)
  )
  #
  ix          <- !is.na(tcr_ids)
  my_clusters <- as.factor(cluster_ids)
  my_tcrs     <- as.factor(tcr_ids)
  #
  m <- Matrix::sparseMatrix(
    i = as.integer(my_tcrs)[ix],
    j = as.integer(my_clusters)[ix]
  )
  # All combinations of cluster ids
  retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
  colnames(retval) <- c("c1", "c2")
  # Number of shared TCRs for each pair of cluster ids
  # intersection / union
  retval$shared <- apply(retval, 1, function(x) {
    sum(m[,x[1]] & m[,x[2]]) /
    sum(m[,x[1]] | m[,x[2]])
  })
  # Shuffle the clusters and count again
  m_shuf <- data.table::rbindlist(pbapply::pblapply(seq(n_iter), function(i) {
    shuf_clusters <- my_clusters
    for (donor_id in unique(donor_ids)) {
      ix_donor <- donor_id == donor_ids
      shuf_clusters[ix_donor] <- sample(my_clusters[ix_donor])
    }
    m <- Matrix::sparseMatrix(
      i = as.integer(my_tcrs)[ix],
      j = as.integer(shuf_clusters)[ix]
    )
    retval <- as.data.frame(t(combn(x = ncol(m), m = 2)))
    colnames(retval) <- c("c1", "c2")
    retval$shared <- apply(retval, 1, function(x) {
      sum(m[,x[1]] & m[,x[2]]) /
      sum(m[,x[1]] | m[,x[2]])
    })
    retval$iter <- i
    retval
  }))
  # Empirical p-values
  retval <- cbind(retval, rbindlist(apply(retval, 1, function(x) {
    val <- x[3]
    ix <- with(m_shuf, c1 == x[1] & c2 == x[2])
    gte <- sum(m_shuf$shared[ix] >= val)
    lte <- sum(m_shuf$shared[ix] <= val)
    data.frame(
      lte = lte,
      lte_pval = (lte + 1) / (sum(ix) + 1),
      gte = gte,
      gte_pval = (gte + 1) / (sum(ix) + 1),
      mean = mean(m_shuf$shared[ix])
    )
  })))
  list(
    shuf = m_shuf,
    pval = retval %>% arrange(gte_pval)
  )
}

  # cluster_pairs <- combn(unique(dd$cluster), 2)
  # tcr_pairs <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
  #   c1 <- cluster_pairs[1,i]
  #   c2 <- cluster_pairs[2,i]
  #   message(sprintf("%s %s", c1, c2))
  #   x <- dd %>%
  #     filter(has_tcr) %>%
  #     # filter(tcr_rank <= 10) %>%
  #     group_by(tcr_clone, case) %>%
  #     summarize(
  #       both_cluster = factor(all(c(c1, c2) %in% cluster)),
  #       expanded = factor(tcr_rank <= 10),
  #       .groups = "keep"
  #     ) %>% unique
  #   tab1 <- with(
  #     x %>% filter(case == "Case"),
  #     table(both_cluster, expanded)
  #   )
  #   tab2 <- with(
  #     x %>% filter(case == "Control"),
  #     table(both_cluster, expanded)
  #   )
  #   retval1 <- broom::tidy(fisher.test(tab1))
  #   retval1$case <- "Case"
  #   retval1$c1 <- c1
  #   retval1$c2 <- c2
  #   retval2 <- broom::tidy(fisher.test(with(
  #     x %>% filter(case == "Control"),
  #     table(both_cluster, expanded)
  #   )))
  #   retval2$case <- "Control"
  #   retval2$c1 <- c1
  #   retval2$c2 <- c2
  #   rbind(retval1, retval2)
  # }))

test_sharing_fisher <- function(case_labels, donor_ids, cluster_ids, tcr_ids, n_iter = 1000) {
  stopifnot(is.logical(case_labels))
  stopifnot(
    length(case_labels) == length(tcr_ids)
  )
  stopifnot(
    length(donor_ids) == length(tcr_ids)
  )
  stopifnot(
    length(cluster_ids) == length(tcr_ids)
  )
  #
  ix_case     <- case_labels
  ix          <- !is.na(tcr_ids)
  my_clusters <- as.factor(cluster_ids)
  my_tcrs     <- as.factor(tcr_ids)
  #
  m1 <- Matrix::sparseMatrix(
    i = as.integer(my_tcrs)[ix & ix_case],
    j = as.integer(my_clusters)[ix & ix_case]
  )
  m2 <- Matrix::sparseMatrix(
    i = as.integer(my_tcrs)[ix & !ix_case],
    j = as.integer(my_clusters)[ix & !ix_case]
  )
  # All combinations of cluster ids
  retval <- as.data.frame(t(combn(x = ncol(m1), m = 2)))
  colnames(retval) <- c("c1", "c2")
  # Number of shared TCRs for each pair of cluster ids
  # intersection / union
  retval$a1 <- apply(retval, 1, function(x) {
    sum(m1[,x[1]] & m1[,x[2]])
  })
  retval$b1 <- apply(retval, 1, function(x) {
    sum(m1[,x[1]] | m1[,x[2]])
  })
  retval$a2 <- apply(retval, 1, function(x) {
    sum(m2[,x[1]] & m2[,x[2]])
  })
  retval$b2 <- apply(retval, 1, function(x) {
    sum(m2[,x[1]] | m2[,x[2]])
  })
  htest <- rbindlist(apply(retval, 1, function(x) {
    mat <- matrix(
      c(x[[3]], x[[4]], x[[5]], x[[6]]),
      nrow = 2,
      dimnames = list(Shared = c("True", "False"), Case = c("True", "False"))
    )
    htest <- fisher.test(x = mat)
    broom::tidy(htest)
  }))
  retval <- cbind(retval, htest)
# TeaTasting <- matrix(
#   c(3, 1, 1, 3),
#   nrow = 2,
#   dimnames = list(Guess = c("Milk", "Tea"),
#   Truth = c("Milk", "Tea")))
# fisher.test(TeaTasting, alternative = "greater")
# ## => p = 0.2429, association could not be established
  # Shuffle the clusters and count again
  m_shuf <- data.table::rbindlist(pbapply::pblapply(seq(n_iter), function(i) {
    shuf_clusters <- my_clusters
    for (donor_id in unique(donor_ids)) {
      ix_donor <- donor_id == donor_ids
      shuf_clusters[ix_donor] <- sample(my_clusters[ix_donor])
    }
    m1 <- Matrix::sparseMatrix(
      i = as.integer(my_tcrs)[ix & ix_case],
      j = as.integer(shuf_clusters)[ix & ix_case]
    )
    m2 <- Matrix::sparseMatrix(
      i = as.integer(my_tcrs)[ix & !ix_case],
      j = as.integer(shuf_clusters)[ix & !ix_case]
    )
    # All combinations of cluster ids
    shuf <- as.data.frame(t(combn(x = ncol(m1), m = 2)))
    colnames(shuf) <- c("c1", "c2")
    # Number of shared TCRs for each pair of cluster ids
    # intersection / union
    shuf$a1 <- apply(shuf, 1, function(x) {
      sum(m1[,x[1]] & m1[,x[2]])
    })
    shuf$b1 <- apply(shuf, 1, function(x) {
      sum(m1[,x[1]] | m1[,x[2]])
    })
    shuf$a2 <- apply(shuf, 1, function(x) {
      sum(m2[,x[1]] & m2[,x[2]])
    })
    shuf$b2 <- apply(shuf, 1, function(x) {
      sum(m2[,x[1]] | m2[,x[2]])
    })
    htest <- rbindlist(apply(shuf, 1, function(x) {
      mat <- matrix(
        c(x[['a1']], x[['b1']], x[['a2']], x[['b2']]),
        nrow = 2,
        dimnames = list(Shared = c("True", "False"), Case = c("True", "False"))
      )
      htest <- fisher.test(x = mat)
      broom::tidy(htest)
    }))
    shuf <- cbind(shuf, htest)
    shuf$iter <- i
    shuf
  }))
  # # Empirical p-values
  # out <- cbind(retval, rbindlist(lapply(seq(nrow(retval)), function(i) {
  #   x <- retval[i,]
  #   val <- x[['p.value']]
  #   ix <- with(m_shuf, c1 == x[['c1']] & c2 == x[['c2']])
  #   # gte <- sum(m_shuf$p.value[ix] <= val)
  #   lte <- sum(m_shuf$p.value[ix] <= val)
  #   data.table(
  #     lte = lte,
  #     lte_pval = (lte + 1) / (sum(ix) + 1),
  #     # gte = gte,
  #     # gte_pval = (gte + 1) / (sum(ix) + 1),
  #     mean = mean(m_shuf$estimate[ix], na.rm = TRUE)
  #   )
  # })))
  # Empirical p-values
  out <- cbind(retval, rbindlist(lapply(seq(nrow(retval)), function(i) {
    x <- retval[i,]
    val <- x[['estimate']]
    ix <- with(m_shuf, c1 == x[['c1']] & c2 == x[['c2']])
    gte <- sum(m_shuf$estimate[ix] >= val)
    lte <- sum(m_shuf$estimate[ix] <= val)
    data.table(
      lte = lte,
      lte_pval = (lte + 1) / (sum(ix) + 1),
      gte = gte,
      gte_pval = (gte + 1) / (sum(ix) + 1),
      mean = mean(m_shuf$estimate[ix], na.rm = TRUE)
    )
  })))
  list(
    shuf = m_shuf,
    pval = out
  )
}

t1_fisher <- test_sharing_fisher(
  case_labels = t1$obs$case == "Case",
  donor_ids   = t1$obs$donor,
  cluster_ids = t1$obs$leiden,
  tcr_ids     = t1$obs$tcr_id
)

t1_fisher2 <- test_sharing_fisher(
  case_labels = t1$obs$case == "Case",
  donor_ids   = t1$obs$donor,
  cluster_ids = t1$obs$leiden,
  tcr_ids     = t1$obs$tcr_id
)

# Fisher p-values are inflated for this data
(
  t1_fisher2$shuf %>% group_by(c1, c2) %>%
  summarize(
    p001 = sum(p.value <= 0.001) / 1000,
    p01 = sum(p.value <= 0.01) / 1000,
    p05 = sum(p.value <= 0.05) / 1000 
  )
)

p <- ggplot(t1_fisher$shuf) +
  geom_histogram(aes(p.value), breaks = seq(0, 1, by = 1/50)) +
  labs(x = "p-value", y = "Count")
my_ggsave(
  "tcr-fisher-pvalues",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 4, height = 2, units = "in", dpi = 300
)

# Empirical p-values are well-calibrated
pvals <- rbindlist(lapply(seq(nrow(t1_share2$pval)), function(i) {
  x <- t1_share2$pval[i,]
  vals <- t1_share2$shuf[(c1 == x[['c1']] & c2 == x[['c2']]) | (c1 == x[['c2']] & c2 == x[['c1']])]$shared
  obs <- sample(vals, size = length(vals), replace = FALSE)
  print(i)
  data.frame(
    c1 = x[['c1']],
    c2 = x[['c2']],
    pval = sapply(obs, function(o) (1 + sum(vals >= o)) / (1 + length(vals)))
  )
}))

pvals %>% group_by(c1, c2) %>%
  summarize(
    p001 = sum(pval <= 0.001) / 1000,
    p01 = sum(pval <= 0.01) / 1000,
    p05 = sum(pval <= 0.05) / 1000
  )

p <- ggplot(pvals) +
  geom_histogram(aes(pval), breaks = seq(0, 1, by = 1/50)) +
  labs(x = "p-value", y = "Count")
my_ggsave(
  "tcr-empirical-pvalues",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 4, height = 2, units = "in", dpi = 300
)


t1_share2$pval %>% arrange(gte_pval) %>% head

t1_share2$pval %>% arrange(lte_pval) %>% head

# dir.create(glue("results/a20/{t1_analysis_name}/tcr_sharing"))

t1_share2$shuf <- t1_share2$shuf %>%
  mutate(pair = paste(c1, c2))
t1_share2$pval <- t1_share2$pval %>%
  mutate(pair = paste(c1, c2))
#
t1_share2$pval$pair <- factor(
  t1_share2$pval$pair,
  (t1_share2$pval %>% arrange(gte_pval))$pair
)
t1_share2$shuf$pair <- factor(
  t1_share2$shuf$pair,
  (t1_share2$pval %>% arrange(gte_pval))$pair
)
#
x <- t1_share2$pval %>%
  mutate(min_pval = gte_pval < lte_pval) %>%
  arrange(min_pval) %>%
  mutate(min_pval = ifelse(min_pval, gte_pval, lte_pval), direction = ifelse(shared > mean, 1, -1)) %>%
  mutate(sign_pval = direction * -log10 (min_pval)) %>%
  arrange(sign_pval)
#
t1_share2$shuf$pair <- factor(
  as.character(t1_share2$shuf$pair),
  as.character(x$pair)
)
#
p <- ggplot() +
  geom_histogram(
    data = t1_share2$shuf,
    mapping = aes(x = shared),
    bins = 30,
    fill = "grey80"
  ) +
  geom_vline(
    data = t1_share2$pval,
    mapping = aes(xintercept = shared),
    color = "red"
  ) +
  geom_text(
    data = x,
    mapping = aes(
      x = shared,
      hjust = ifelse(direction < 0, -0.1, 1.1),
      y = Inf,
      label = signif(shared, 2)
    ),
    vjust = 1.1,
    size = 6
  ) +
  geom_text(
    data = x,
    mapping = aes(
      x = ifelse(direction > 0, -Inf, Inf),
      hjust = ifelse(direction > 0, -0.1, 1.1),
      y = Inf,
      label = sprintf("p=%s", str_replace(signif(min_pval, 1), "-0", "-"))
    ),
    vjust = 1.1,
    size = 8
  ) +
  facet_wrap(vars(pair), scales = "free", ncol = 3, strip.position = "left") +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 24),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(x = NULL, y = NULL)
my_ggsave(
  "tcr-cluster-sharing",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1, width = 20, height = 20, units = "in", dpi = 300
)


t1_share2$pval <- t1_share2$pval %>%
  mutate(
    z = shared / mean,
    c1 = factor(c1),
    c2 = factor(c2)
  )


sorted_c1 <- umap_sorted_rows(
  data = rbind(
    t1_share2$pval,
    t1_share2$pval %>% mutate(c3 = c2) %>% mutate(c2 = c1, c1 = c3) %>% select(-c3)
  ),
  formula = "c1 ~ c2",
  value.var = "z"
)
# sorted_c2 <- umap_sorted_rows(data = t1_share2$pval, formula = "c2 ~ c1", value.var = "z")

# stopifnot(
#   all(levels(t1_share$pval$c1) %in% levels(t1_share$pval$c2))
# )

t1_share2$pval$c3 <- apply(t1_share2$pval[,c("c1","c2")], 1, function(x) min(as.numeric(x)))
t1_share2$pval$c4 <- apply(t1_share2$pval[,c("c1","c2")], 1, function(x) max(as.numeric(x)))
head(t1_share2$pval)
t1_share2$pval$c1 <- t1_share2$pval$c3
t1_share2$pval$c2 <- t1_share2$pval$c4

t1_share2$pval$c1 <- naturalfactor(t1_share2$pval$c1)
t1_share2$pval$c2 <- naturalfactor(t1_share2$pval$c2)

# t1_share2$pval <- rbind(
#   t1_share2$pval,
#   t1_share2$pval %>% mutate(c3 = c2) %>% mutate(c2 = c1, c1 = c3) %>% select(-c3)
# ) %>% filter(as.integer(c1) < as.integer(c2))

t1_share2$pval <- t1_share2$pval %>%
  # mutate(
  #   z = shared / mean,
  #   c1 = factor(as.character(c1), sorted_c1),
  #   c2 = factor(as.character(c2), sorted_c1)
  # ) %>%
  mutate(min_pval = gte_pval < lte_pval) %>%
  mutate(min_pval = ifelse(min_pval, gte_pval, lte_pval), direction = ifelse(shared > mean, 1, -1)) %>%
  mutate(min_fdr = p.adjust(min_pval, method = "fdr"))

x <- log2(t1_share2$pval$z)
col_fun <- circlize::colorRamp2(
  seq(-max(abs(x)), max(abs(x)), length.out = 11),
  # seq(min(abs(x)), max(abs(x)), length.out = 11),
  # quantile(x, probs = seq(0, 1, length.out = 11)),
  rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
)
# color_values <- scales::rescale(seq(min(x$logFC), max(x$logFC), length.out = 11))
color_values <- seq(min(x), max(x), length.out = 11)
# color_values <- quantile(x, probs = seq(0, 1, length.out = 11))
t1_share2$pval <- t1_share2$pval %>% mutate(bonf = min_pval < 0.05 / n())
# 0.05 / nrow(t1_share2$pval)
p <- ggplot(t1_share2$pval) +
  aes(x = c1, y = c2) +
  geom_tile(aes(fill = log2(z))) +
  geom_point(
    data = t1_share2$pval %>% filter(bonf),
    mapping = aes(x = c1, y = c2),
    shape = 21, size = 3, fill = "white", stroke = 0.3
  ) +
  # geom_text(
  #   data = t1_share2$pval %>% filter(min_fdr < 0.05),
  #   mapping = aes(label = str_replace(format.pval(min_pval, digits = 1), "-0", "-")),
  #   size = 3
  # ) +
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  scale_fill_gradientn(
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))[6:11],
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    # values = color_values,
    colors = col_fun(color_values),
    # limits = c(0, 1),
    # name = bquote("log"[2]~"FC"),
    name = "z",
    guide = guide_colorbar(barheight = 12),
    breaks = -8:8,
    labels = function(x) fractional::fractional(2^x)
  ) +
  labs(x = NULL, y = NULL)
n_clusters <- length(unique(t1_share2$pval$c1))
my_ggsave(
  "tcr-cluster-sharing-heatmap",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = 2 + n_clusters * 0.3,
  height = n_clusters * 0.3,
  units = "in",
  dpi = 300
)

x <- t1_share2$pval$shared
col_fun <- circlize::colorRamp2(
  seq(0, max(abs(x)), length.out = 9),
  # seq(min(abs(x)), max(abs(x)), length.out = 11),
  # quantile(x, probs = seq(0, 1, length.out = 11)),
  RColorBrewer::brewer.pal(name = "Greys", n = 9)
)
# color_values <- scales::rescale(seq(min(x$logFC), max(x$logFC), length.out = 11))
color_values <- seq(min(x), max(x), length.out = 9)
# color_values <- quantile(x, probs = seq(0, 1, length.out = 11))
p <- ggplot(t1_share2$pval) +
  aes(x = c1, y = c2) +
  geom_tile(aes(fill = shared)) +
  scale_fill_gradientn(
    colors = col_fun(color_values),
    # limits = c(0, 1),
    name = "Percent",
    guide = guide_colorbar(barheight = 12),
    # breaks = -8:8,
    # labels = function(x) fractional::fractional(2^x)
    labels = function(x) 100 * x
  ) +
  labs(x = NULL, y = NULL)
n_clusters <- length(unique(t1_share2$pval$c1))
my_ggsave(
  "tcr-cluster-sharing-heatmap-pct",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = 2 + n_clusters * 0.3,
  height = n_clusters * 0.3,
  units = "in",
  dpi = 300
)


t1_share$pval %>% head

p <- ggplot() +
  geom_vline(xintercept = 1, size = 0.3) +
  geom_point(
    data = t1_share$pval,
    mapping = aes(x = z, y = -log10(min_pval))
  ) +
  geom_point(
    data = t1_share$pval %>% filter(z > 1, min_fdr < 0.05),
    mapping = aes(x = z, y = -log10(min_pval)),
    color = "red"
  ) +
  geom_text_repel(
    data = t1_share$pval %>% filter(z > 1, min_fdr < 0.05),
    mapping = aes(x = z, y = -log10(min_pval), label = pair)
  ) +
  labs(
    title = "TCR sharing between pairs of clusters",
    x = "z", y = bquote("-Log"[10]~"p"))
my_ggsave(
  "tcr-cluster-sharing-volcano",
  out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = 4.5,
  height = 3.5,
  units = "in",
  dpi = 300
)

# show the cluster-pair sharing for each donor
#
# > t1_share2$pval %>% head
#   c1 c2     shared  lte  lte_pval  gte   gte_pval       mean
# 1  8 12 0.05555556 9998 0.9998000    2 0.00029997 0.02756358
# 2  1 12 0.02190722 9816 0.9816018  184 0.01849815 0.01566056
# 3  5  7 0.11646137 9542 0.9542046  458 0.04589541 0.10358374
# 4  4  7 0.11290323 9029 0.9029097  974 0.09749025 0.10342456
# 5  5 12 0.03030303 8818 0.8818118 1198 0.11988801 0.02410403
# 6  1  9 0.06115515 8354 0.8354165 1646 0.16468353 0.05625611
#
t1_analysis_name <- "a12_4_4_t4_cd8_1_2"
t1_file <- glue("results/a20/{t1_analysis_name}/data/{t1_analysis_name}.qs")
file.exists(t1_file)
t1 <- qread(t1_file)
table(t1$obs$leiden)
t1$obs <- t1$obs %>%
  dplyr::mutate(tcr_id = paste(TRAV, TRBV, TRA_cdr3_trim, TRB_cdr3_trim)) %>%
  dplyr::mutate(tcr_id = ifelse(has_tcr, tcr_id, NA))

c1 <- "8"
c2 <- "12"
for (c2 in c("12", "3")) {
  tcr_id1 <- t1$obs %>% filter(leiden == c1) %>% select(donor, case, tcr_id) %>% mutate(cluster = c1)
  tcr_id2 <- t1$obs %>% filter(leiden == c2) %>% select(donor, case, tcr_id) %>% mutate(cluster = c2)
  tcr_id <- rbind(tcr_id1, tcr_id2)
  pct_int <- function(x, y) length(intersect(x, y)) / length(union(x, y))
  my_pct <- tcr_id %>% 
    group_by(donor, case) %>%
    summarize(pct_shared = pct_int(tcr_id[cluster == c1], tcr_id[cluster == c2]))
  my_pct
  p <- ggplot(my_pct) +
    # aes(y = pct_shared, x = case, fill = case) +
    aes(y = pct_shared, x = "") +
    geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.3) +
    geom_point(
      size = 3, shape = 19, stroke = 0.3,
      position = position_quasirandom(width = 0.5)
    ) +
    scale_fill_manual(values = pals::okabe()[2:1], guide = "none") +
    scale_y_continuous(
      labels = \(x) 100 * x
    ) +
    labs(
      x = NULL, y = "Percent sharing",
      title = glue("Percent shared TCRs ({c1}, {c2})")
    ) +
    theme(
      plot.title.position = "plot"
    )
  my_ggsave(
    glue("boxploth-tcr-sharing-clusters-{c1}-{c2}"),
    out_dir = glue("results/a20/{t1_analysis_name}/figures/tcr"),
    type = "pdf",
    plot = p,
    scale = 1, width = 4, height = 5, units = "in", dpi = 300
  )
}

# }}}

# }}}


########################################################################
analysis_name <- "blood3"
a1_file <- as.character(glue("results/blood/{analysis_name}/{analysis_name}.qs"))
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
    equal_pcs     = FALSE,
    max_pcs       = 50,
    n_harmony     = 30,
    harmony_vars  = c("batch"),
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
    obs           = obs[ix_keep,],
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
  out_dir       = glue("results/blood/{analysis_name}/figures"),
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
  out_dir           = glue("results/blood/{analysis_name}/figures"),
  cb_dir            = glue("cellbrowser/colitis/blood/{analysis_name}"),
  ensembl_to_symbol = ensembl_to_symbol
)

# Analysis for each batch, separately
#########################################################################
#my_batches <- unique(obs$batch)
#for (my_batch in my_batches) {
#  analysis_name <- my_batch
#  a1_file <- as.character(glue("results/blood/{analysis_name}/{analysis_name}.qs"))
#  file.exists(a1_file)
#  dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#  #
#  ix_keep <- obs$batch == my_batch
#  stopifnot(sum(ix_keep) > 0)
#  message(glue("{analysis_name} {sum(ix_keep)} cells"))
#  #
#  source("R/functions/do-analysis.R")
#  if (!file.exists(a1_file)) {
#    #
#    params <- list(
#      analysis_name = analysis_name,
#      n_genes       = nrow(counts[,ix_keep]),
#      n_cells       = ncol(counts[,ix_keep]),
#      min_percent   = 100 * (50 / ncol(counts[,ix_keep])),
#      loess_span    = 0.01,
#      n_pcs         = 'mcv',
#      max_pcs       = 40,
#      n_harmony     = 0,
#      # harmony_vars  = c("channel"),
#      n_knn         = 50,
#      leiden_res    = seq(0.5, 1.8, length.out = 10),
#      leiden_iter   = 10,
#      umap_spread   = 1,
#      umap_min_dist = 0.25,
#      log_file      = as.character(glue("{dirname(a1_file)}/analysis.log"))
#    )
#    a1_params_file <- glue("{dirname(a1_file)}/params.json")
#    writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE), a1_params_file) 
#    #
#    a1 <- run_analysis(
#      obs           = obs[ix_keep,],
#      counts        = counts[,ix_keep],
#      exclude_genes = c(mito_genes, bcr_genes),
#      mito_genes    = mito_genes,
#      params        = params
#    )
#    a1$params <- params
#    print_status(glue("Writing {a1_file}"))
#    qsave(a1, a1_file)
#    print_status(glue("done"))
#  } else {
#    print_status(glue("Reading {a1_file}"))
#    a1 <- qread(a1_file)
#    print_status(glue("done"))
#  }
#  source("R/colors-int.R")
#  source("R/plot-analysis.R")
#  a1$obs$leiden <- a1$obs$leiden0.933
#  assign(analysis_name, a1)
#  try({
#  plot_analysis(
#    analysis_name = analysis_name,
#    out_dir       = glue("results/blood/{analysis_name}/figures"),
#    rowname_key   = ensembl_to_symbol,
#    exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes)),
#    do_pb         = FALSE
#  )
#  })
#  assign(analysis_name, a1)
#  # rm(list = c("a1", analysis_name))
#  source("R/functions/make-cellbrowser.R")
#  make_cellbrowser(
#    analysis_name     = analysis_name,
#    out_dir           = glue("results/blood/{analysis_name}/figures"),
#    cb_dir            = glue("cellbrowser/colitis/blood/{analysis_name}"),
#    ensembl_to_symbol = ensembl_to_symbol
#  )

obs <- qread(file.path(out_dir, "obs.qs"))
counts <- qread(file.path(out_dir, "counts.qs"))

rna_cell_ids <- colnames(counts)
rna_cell_ids[1:5]

normalize_clr <- function(counts) {
  # cells are columns
  denom <- exp(colSums(log1p(counts)) / ncol(counts))
  counts@x <- log1p(counts@x / rep(denom, diff(counts@p)))
  counts
}

antibodies <- read_excel(
  "data/biolegend/197 antibodies TOTAL-seqC_121819QG_99328_AntibodyList_Barcodes updated.xlsx"
)
antibodies$DNA_ID <- sprintf("ADT_C%04d", antibodies$DNA_ID)
antibodies$id <- str_replace(str_replace(sprintf("%s %s",
  str_split_fixed(antibodies$Target, " ", 3)[,3],
  antibodies$DNA_ID
), "  ", " "), "^ ", "")
ab_to_id <- unlist(split(antibodies$id, antibodies$DNA_ID))

adt_files <- Sys.glob(
  glue(
    "/projects/irae_blood/cellranger_output/C*_fbc/*_fbc.csv"
  )
)
adt_batches <- basename(dirname(adt_files))
read_adt <- function(adt_file) {
  adt <- read_sparse_csv(adt_file)
  adt <- adt[!stringr::str_detect(rownames(adt), "^HT_"),]
  adt <- adt[!str_detect(rownames(adt), "^HT"),]
  adt <- adt[,Matrix::colSums(adt) > 0]
  batch <- str_replace(basename(dirname(adt_file)), "_fbc", "_gex")
  colnames(adt) <- sprintf("%s|%s", batch, colnames(adt))
  ix <- colnames(adt) %in% rna_cell_ids
  stopifnot(sum(ix) > 0)
  adt <- adt[,ix]
  rownames(adt) <- ab_to_id[rownames(adt)]
  return(adt)
}
adt_counts <- do.call(cbind, pblapply(adt_files, read_adt))
adt_clr <- normalize_clr(adt_counts)

both_cell_ids <- intersect(colnames(adt_clr), colnames(counts))
length(both_cell_ids)
length(colnames(counts))
length(colnames(adt_clr))

for (my_protein in rownames(adt_clr)) {
  x <- left_join(
    x = a1$obs %>% select(cell, UMAP1, UMAP2),
    y = data.frame(cell = colnames(adt_clr), protein = adt_clr[my_protein,]),
    by = "cell"
  )
  x$protein[is.na(x$protein)] <- 0
  p <- plot_hexgene(
    x            = x$UMAP1,
    y            = x$UMAP2,
    z            = x$protein,
    bins         = 201,
    palette      = "davos",
    direction    = -1,
    use_quantile = FALSE,
    legend       = FALSE,
    italic       = FALSE
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title = "CLR",
      barwidth = 15
    )
  ) +
  labs(subtitle = my_protein)
  my_ggsave(
    glue("umap-{safe(my_protein)}"),
    out_dir = glue("{out_dir}/{analysis_name}/figures/adt-umap"),
    type = "pdf",
    plot = p,
    scale = 1, width = 4, height = 5, units = "in", dpi = 300
  )
}


adt_loess <- do_loess(adt_counts, NULL, loess_span = 0.1, min_percent = 1)

  p <- ggplot() +
    geom_point(
      data = adt_loess$counts_stats,
      mapping = aes(log10(mean), residuals)
    ) +
    geom_hline(yintercept = 0, size = 0.3, linetype = 2) +
    labs(
      x = bquote("Log"[10]~"Mean"),
      y = bquote("Log"[10]~"SD (resid.)"),
      title = glue::glue("Select {comma(sum(adt_loess$counts_stats$include))} of {comma(nrow(adt_loess$counts_stats))} total genes"),
      fill = "Genes"
    )
  p
my_ggsave(
  glue("mean-resid"),
  out_dir = glue("{out_dir}/{analysis_name}/figures/adt-loess"),
  type = "pdf",
  plot = p,
  scale = 1, width = 4, height = 3, units = "in", dpi = 300
)

exclude_proteins <- adt_loess$counts_stats$percent < 20 | str_detect(adt_loess$counts_stats$gene, "Ctrl")
exclude_proteins <- names(rownames(adt_counts))[exclude_proteins]

adt_counts2 <- adt_counts
rownames(adt_counts2) <- names(rownames(adt_counts))

id_to_feature <- c(ensembl_to_symbol, rownames(adt_counts))

########################################################################
analysis_name <- "blood2-rna-adt"
a1_file <- as.character(glue("results/blood/{analysis_name}/{analysis_name}.qs"))
file.exists(a1_file)
dir.create(dirname(a1_file), recursive = TRUE, showWarnings = FALSE)
#
ix_keep <- obs$cell %in% both_cell_ids
#
source("R/functions/do-analysis-citeseq.R")
if (!file.exists(a1_file)) {
  #
  params <- list(
    analysis_name = analysis_name,
    exclude_proteins = exclude_proteins,
    loess_span    = 0.01,
    n_pcs         = 'mcv',
    max_pcs       = 90,
    equal_pcs     = FALSE,
    n_harmony     = 30,
    harmony_vars  = c("batch"),
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
  a1 <- run_analysis_citeseq(
    obs              = obs[ix_keep,],
    rna_counts       = counts[,ix_keep],
    exclude_genes    = c(mito_genes, bcr_genes),
    mito_genes       = mito_genes,
    adt_count        = adt_counts2,
    exclude_proteins = exclude_proteins,
    params           = params
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
a1$counts <- rbind(a1$rna_counts, a1$adt_counts)
assign(analysis_name, a1)
try({
plot_analysis(
  analysis_name = analysis_name,
  out_dir       = glue("results/blood/{analysis_name}/figures"),
  rowname_key   = id_to_feature,
  exclude_genes = unique(c(tcr_genes, bcr_genes, mito_genes, exclude_proteins)),
  do_pb         = FALSE
)
})

