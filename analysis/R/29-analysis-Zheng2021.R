
# libraries {{{
pacman::p_load(
  dendsort,
  data.table,
  dplyr,
  forcats,
  ggplot2,
  ggforestplot,
  ggbeeswarm,
  ggrastr,
  ggrepel,
  ggtext,
  glue,
  janitor,
  limma,
  Matrix,
  magrittr,
  naturalsort,
  parameters,
  patchwork,
  pbapply,
  rhdf5,
  RColorBrewer,
  SummarizedExperiment,
  qs,
  readxl,
  scales,
  stringr,
  tibble,
  tidyr
)
source("R/functions/theme-kamil.R")
theme_set(theme_kamil)
source("R/functions/helpers.R")
source("R/functions/mpn65.R")
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(data.table::first)

# }}}

# globals {{{

hclust_dataframe <- function(d, x, y, z) {
  res_mat <- d %>%
    pivot_wider(id_cols = all_of(x), names_from = all_of(y), values_from = all_of(z)) %>%
    column_to_rownames(x) %>%
    as.matrix
  res_h1 <- dendsort(hclust(dist(res_mat)))
  res_h2 <- dendsort(hclust(dist(t(res_mat))))
  res_mat <- res_mat[res_h1$order, res_h2$order]
}

# Path with data
root_dir <- "/projects/external_data/Zheng2021"

# Path to output files (tsv, pdf)
out_dir <- "results/Zheng2021"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# }}}

# rds files {{{

cd8_rds_files <- Sys.glob(glue("{root_dir}/data/expression/CD8/byDataset/*.sce.rds"))
cd8_rds_names <- str_remove(basename(cd8_rds_files), "\\.sce\\.rds$")
cd8_rds_data <- pblapply(cd8_rds_files, readRDS)
names(cd8_rds_data) <- cd8_rds_names

cd4_rds_files <- Sys.glob(glue("{root_dir}/data/expression/CD4/byDataset/*.sce.rds"))
cd4_rds_names <- str_remove(basename(cd4_rds_files), "\\.sce\\.rds$")
cd4_rds_data <- pblapply(cd4_rds_files, readRDS)
names(cd4_rds_data) <- cd4_rds_names

# }}}

# CD8 dataset files {{{

gene_percents <- function(my_name, my_genes) {
  sce <- cd8_rds_data[[my_name]]
  sce_genes <- rowData(sce) %>% as.data.frame
  sce_meta <- colData(sce) %>% as.data.frame
  if ("counts" %in% names(assays(sce))) {
    sce_counts <- assay(sce, "counts")
  } else {
    sce_counts <- assay(sce, "norm_exprs")
  }
  if (!all(sce_genes$display.name == rownames(sce_counts))) {
    rownames(sce_counts) <- sce_genes$display.name
  }
  my_genes <- base::intersect(my_genes, rownames(sce_counts))
  for (my_gene in my_genes) {
    sce_meta[[glue("gene_{my_gene}")]] <- as.numeric(sce_counts[my_gene,])
  }
  sce_meta %>%
    select(cellID, patient, meta.cluster, starts_with("gene_")) %>%
    pivot_longer(starts_with("gene_"), names_to = "gene") %>%
    group_by(patient, meta.cluster, gene) %>%
    summarize(
      cells_total = n(),
      cells_pos = sum(value > 0),
      cells_freq = sum(value > 0) / n(),
      .groups = "keep"
    ) %>%
    mutate(dataset = my_name, gene = str_remove(gene, "gene_")) %>%
    relocate(dataset)
}
cd8_meta <- rbindlist(pblapply(cd8_rds_names, function(rds_name) {
  my_genes <- c("IL26", "IL17A", "CXCL13")
  gene_percents(rds_name, my_genes)
}))
cd8_meta %<>% mutate(cancer = str_split_fixed(dataset, "\\.", 2)[,1])

# bad boxplots {{{ 
my_d <- cd8_meta %>%
  filter(cells_total > 10) %>%
  group_by(dataset, patient, gene) %>%
  summarize(
    freq = sum(cells_pos) / sum(cells_total),
    cells = sum(cells_total),
    .groups = "keep"
  )
my_d
p <- ggplot(my_d) +
  aes(x = dataset, y = freq) +
  facet_grid(rows = vars(gene)) +
  geom_boxplot() +
  scale_x_discrete(position = "top") +
  scale_y_continuous(labels = \(x) signif(100 * x, 2)) +
  theme(
    axis.text.x.top = element_text(angle = 45, hjust = 0)
  )
my_ggsave(
  slug = glue("boxplot-IL17A-IL26-CXCL13"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = 15,
  height = 6,
  plot = p
)

# my_cancers <- c(
#   "HCC", "MM", "RC", "THCA", "BRCA", "CRC", "BCL", "UCEC", "NPC", "CHOL",
#   "LC", "OV", "PACA", "STAD", "SCC", "BCC", "ESCA"
# )
my_datasets <- c(
  "HCC.ChunhongZheng2017",
  "HCC.QimingZhang2019_10X",
  "HCC.QimingZhang2019_SS2",
  "MM.thisStudy",
  "RC.MatthewDYoung2018",
  "RC.thisStudy",
  "THCA.thisStudy",
  "BRCA.PeterSavas2018",
  "BRCA.thisStudy",
  "CRC.LeiZhang2018",
  "CRC.LeiZhang2020_10X",
  "BCL.thisStudy",
  "UCEC.thisStudy",
  "NPC.YangLiu2021",
  "CHOL.thisStudy",
  "LC.DietherLambrechts2018",
  "LC.QianqianSong2019",
  "LC.RapolasZilionis2019",
  "LC.XinyiGuo2018",
  "OV.thisStudy",
  "PACA.JunyaPeng2019",
  "PACA.thisStudy",
  "STAD.BoxiKang2019",
  "SCC.KathrynEYost2019",
  "BCC.KathrynEYost2019",
  "ESCA.thisStudy",
  "AML.PeterVanGalen2019",
  "FTC.thisStudy",
  "HNSCC.SidharthVPuram2017",
  "LIHC.LichunMa2019",
  "MELA.HanjieLi2019",
  "MELA.MosheSade-Feldman2018"
)
my_d <- cd8_meta %>%
  # filter(cells_total > 10) %>%
  filter(gene == "IL26") %>%
  # filter(str_detect(meta.cluster, "Tex")) %>%
  group_by(dataset, patient, gene) %>%
  summarize(
    freq = sum(cells_pos) / sum(cells_total),
    cells = sum(cells_total),
    .groups = "keep"
  ) %>%
  mutate(log10freq = log10(freq * 100 + 1))
my_d$dataset <- factor(my_d$dataset, my_datasets)
# my_d_inorder <- hclust_dataframe(
#   my_d %>% group_by(dataset, gene) %>% summarize(log10freq = median(log10freq)),
#   "dataset", "gene", "log10freq"
# )
# my_d$dataset <- factor(my_d$dataset, rownames(my_d_inorder))
p <- ggplot(my_d) +
  aes(y = dataset, x = freq) +
  facet_grid(cols = vars(gene)) +
  geom_stripes(odd = "#11111111", even = "#00000000") +
  # geom_vline(xintercept = c(0.1, 1, 10, 100), linewidth = 0.3) +
  geom_boxplot() +
  scale_y_discrete(position = "right") +
  # scale_x_continuous(trans = "log10", labels = \(x) signif(100 * x, 2)) +
  labs(
    title = "Percent of cells (per donor) positive for each gene",
    x = NULL,
    y = NULL,
    caption = "data from Zheng et al 2021"
  ) +
  theme(
    axis.text.y.right = element_text(hjust = 0)
  )
my_ggsave(
  slug = glue("boxploth-IL17A-IL26-CXCL13"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = 12,
  height = 9,
  plot = p
)
p <- ggplot(my_d) +
  aes(x = dataset, y = freq, fill = dataset) +
  geom_boxplot() +
  geom_hline(yintercept = c(0.05, 0.25), linewidth = 0.3, alpha = 0.3) +
  scale_fill_manual(values = mpn65, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    title = "Percent of cells (per donor) positive for each gene",
    x = NULL,
    y = NULL,
    caption = "data from Zheng et al 2021"
  ) +
  theme(
    axis.text.x.top = element_text(angle = 45, hjust = 0)
  )
my_ggsave(
  slug = glue("boxplot-IL17A-IL26-CXCL13"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = 12,
  height = 6,
  plot = p
)
# }}}

# boxplot-CD8-{my_gene} {{{
my_genes <- c("IL26", "IL17A", "CXCL13")
for (my_gene in my_genes) {
  # my_cancers <- c( "HCC", "MM", "RC", "THCA", "BRCA", "CRC", "BCL", "UCEC",
  #   "NPC", "CHOL", "LC", "OV", "PACA", "STAD", "SCC", "BCC", "ESCA", "AML",
  #   "FTC", "HNSCC", "LIHC", "MELA")
  my_datasets <- (cd8_meta %>% count(dataset, patient) %>% count(dataset) %>% filter(n > 1))$dataset
  my_d <- cd8_meta %>%
    filter(dataset %in% my_datasets) %>%
    filter(cells_total > 5) %>%
    filter(gene == my_gene) %>%
    # filter(str_detect(meta.cluster, "Tex")) %>%
    group_by(cancer, patient, gene) %>%
    summarize(freq = sum(cells_pos) / sum(cells_total), .groups = "keep") %>%
    mutate(log10freq = log10(freq * 100 + 1))
  # my_d$cancer <- factor(my_d$cancer, my_cancers)
  p <- ggplot(my_d) +
    aes(x = cancer, y = freq, fill = cancer) +
    geom_boxplot() +
    # geom_hline(yintercept = c(0.05, 0.25), linewidth = 0.3, alpha = 0.3) +
    scale_fill_manual(values = mpn65, guide = "none") +
    scale_x_discrete(position = "bottom") +
    scale_y_continuous(labels = function(d) signif(100 * d, 2)) +
    labs(
      title = glue("Percent of CD8 T cells positive for <i>{my_gene}</i>, per patient"),
      x = NULL,
      y = NULL,
      caption = "data from Zheng et al 2021"
    ) +
    theme(
      plot.title = element_markdown(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line()
    )
  my_ggsave(
    slug = glue("boxplot-CD8-{my_gene}"),
    out_dir = file.path(out_dir, "freq", "boxplot"),
    type = "pdf",
    width = 8,
    height = 3,
    plot = p
  )
}
# }}}

# boxplot-CD8-{paste(my_genes)} {{{
my_d <- cd8_meta %>%
  filter(cells_total > 5) %>%
  filter(gene %in% my_genes) %>%
  group_by(cancer, patient, gene) %>%
  summarize(freq = sum(cells_pos) / sum(cells_total), .groups = "keep") %>%
  mutate(log10freq = log10(freq * 100 + 1))
p <- ggplot(my_d) +
  aes(x = cancer, y = freq, fill = cancer) +
  facet_grid(rows = vars(gene), switch = "y") +
  geom_boxplot() +
  scale_fill_manual(values = mpn65, guide = "none") +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(labels = function(d) signif(100 * d, 2), position = "right") +
  labs(
    title = glue("Percent of CD8 T cells positive for each gene, per patient"),
    x = NULL,
    y = NULL,
    caption = "data from Zheng et al 2021"
  ) +
  theme(
    plot.title = element_markdown(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(),
    strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1)
  )
my_ggsave(
  slug = glue("boxplot-CD8-{paste(my_genes, collapse = '-')}"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = 8,
  height = 8,
  plot = p
)
# }}}

# boxploth-CD8-paste(my_genes) {{{
my_datasets <- (cd8_meta %>% count(dataset, patient) %>% count(dataset) %>% filter(n > 1))$dataset
my_d <- cd8_meta %>%
  filter(dataset %in% my_datasets) %>%
  # filter(cells_total > 5) %>%
  filter(gene %in% my_genes) %>%
  group_by(cancer, patient, gene) %>%
  summarize(
    freq = sum(cells_pos) / sum(cells_total),
    cells = sum(cells_total),
    .groups = "keep"
  ) %>%
  mutate(log10freq = log10(freq * 100 + 1)) %>%
  ungroup()
my_d <- left_join(
  my_d,
  my_d %>%
    select(cancer, patient) %>%
    unique %>%
    count(cancer) %>%
    mutate(cancer_n = glue("{cancer} (n = {n})")) %>%
    select(-n),
  by = "cancer"
)
#
my_d$cancer_n <- factor(my_d$cancer_n, (my_d %>%
  group_by(cancer_n, gene) %>%
  summarize(mean = mean(freq)) %>%
  group_by(cancer_n) %>%
  summarize(mean = mean(mean)) %>%
  arrange(mean)
)$cancer_n)
#
my_d_mean <- my_d %>%
  group_by(cancer_n, gene) %>%
  summarize(mean = mean(freq)) %>%
  mutate(label = glue("{signif(100 * mean, 2)}%"))
#
p <- ggplot(my_d) +
  aes(y = cancer_n, x = freq, fill = cancer_n) +
  facet_grid(cols = vars(gene), scales = "free_x") +
  geom_stripes(odd = "#00000000", even = "#00000011") +
  # geom_boxplot(outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = c(0), linewidth = 0.3) +
  geom_point(position = position_quasirandom(), size = 0.5, color = "grey40") +
  geom_crossbar(
    data = my_d_mean, mapping = aes(x = mean, xmin = mean, xmax = mean, y = cancer_n),
    linewidth = 0.3, fatten = 1.0, color = "red"
  ) +
  geom_text(
    data = my_d_mean,
    mapping = aes(x = 0, y = cancer_n, label = label, alpha = mean),
    hjust = 1.15, size = 5
  ) +
  scale_alpha_continuous(guide = "none", trans = "log10") +
  scale_fill_manual(values = mpn65, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0.8, 0.2)),
    breaks = function(d) { d[1] <- 0; pretty_breaks()(d, 3) },
    labels = function(d) ifelse(d >= 0, signif(100 * d, 2), "")
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = glue("Percent of CD8 T cells positive for each gene, per patient"),
    x = NULL,
    y = NULL,
    caption = "data from Zheng et al 2021"
  ) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.4, color = "#00000011"),
    panel.spacing = unit(0.5, "lines"),
    plot.title.position = "plot",
    plot.title = element_markdown(),
    strip.text.x = element_text(face = "italic")
    # strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1)
  )
my_height <- 1 + length(unique(my_d$cancer_n)) * 0.25
my_ggsave(
  slug = glue("boxploth-CD8-{paste(my_genes, collapse = '-')}"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = 8,
  height = my_height,
  plot = p
)

# Percent of CD8 T cells with IL17A in our data {{{

h5_file <- "paper/tissue-cd8.h5"
our_obs <- h5read(h5_file, "obs")
our_counts <- sparseMatrix(
  i = as.numeric(h5read(h5_file, "matrix/indices")),
  p = as.numeric(h5read(h5_file, "matrix/indptr")),
  x = as.numeric(h5read(h5_file, "matrix/data")),
  dims = as.numeric(h5read(h5_file, "matrix/shape")),
  index1 = FALSE
)
rownames(our_counts) <- h5read(h5_file, "matrix/features")
colnames(our_counts) <- h5read(h5_file, "matrix/barcodes")
our_counts[1:5,1:5]

genes <- fread("data/ensembl_id-symbol.tsv")

our_ens <- genes$ensembl_id[genes$symbol == "IL17A"]
our_obs$IL17A <- as.numeric(our_counts[our_ens,])

our_obs

for (my_gene in my_genes) {
  my_ens <- with(genes, ensembl_id[symbol == my_gene])
  our_obs[[glue("gene_{my_gene}")]] <- as.numeric(our_counts[my_ens,])
}

our_d <- our_obs %>%
  select(case, donor, starts_with("gene_")) %>%
  pivot_longer(starts_with("gene_"), names_to = "gene") %>%
  group_by(donor, case, gene) %>%
  summarize(
    freq = sum(value > 0) / n(),
    cells = n(),
    .groups = "keep"
  ) %>%
  mutate(log10freq = log10(freq * 100 + 1)) %>%
  ungroup() %>%
  mutate(
    gene = str_remove(gene, "gene_"),
    cancer = glue("irColitis {case}")
  ) %>%
  dplyr::rename(patient = donor) %>%
  relocate(cancer, patient, gene, freq) %>%
  select(-case)

our_d$cancer_n <- NULL
my_d$cancer_n <- NULL

our_d$data <- "Our data"
my_d$data <- "Zheng et al 2021"

head(our_d)
head(my_d)

comb_d <- rbind(our_d, my_d)
comb_d <- left_join(
  comb_d,
  comb_d %>%
    select(cancer, patient) %>%
    unique %>%
    count(cancer) %>%
    mutate(cancer_n = glue("{cancer} (n = {n})")) %>%
    select(-n),
  by = "cancer"
)
comb_d_mean <- comb_d %>%
  group_by(data, cancer_n, gene) %>%
  summarize(mean = mean(freq), .groups = "keep") %>%
  ungroup() %>%
  mutate(label = glue("{signif(100 * mean, 2)}%"))
#
p <- ggplot(comb_d) +
  aes(y = cancer_n, x = freq, fill = cancer_n) +
  facet_grid(
    rows = vars(data), cols = vars(gene), scales = "free", space = "free_y",
    switch = "y"
  ) +
  geom_stripes(odd = "#00000000", even = "#00000011") +
  # geom_boxplot(outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = c(0), linewidth = 0.3) +
  geom_point(position = position_quasirandom(), size = 0.5, color = "grey40") +
  geom_crossbar(
    data = comb_d_mean, mapping = aes(x = mean, xmin = mean, xmax = mean, y = cancer_n),
    linewidth = 0.3, fatten = 1.0, color = "red"
  ) +
  geom_text(
    data = comb_d_mean,
    mapping = aes(x = 0, y = cancer_n, label = label, alpha = mean),
    hjust = 1.15, size = 5
  ) +
  scale_alpha_continuous(guide = "none", trans = "log10") +
  scale_fill_manual(values = mpn65, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0.8, 0.2)),
    breaks = function(d) { d[1] <- 0; pretty_breaks()(d, 3) },
    labels = function(d) ifelse(d >= 0, signif(100 * d, 2), "")
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = glue("Percent of CD8 T cells positive for each gene, per patient"),
    x = NULL,
    y = NULL
  ) +
  theme(
    strip.text.y.left = element_text(angle = 0),
    panel.grid.major.x = element_line(linewidth = 0.4, color = "#00000011"),
    panel.spacing = unit(0.5, "lines"),
    plot.title.position = "plot",
    plot.title = element_markdown(),
    strip.text.x = element_text(face = "italic")
    # strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1)
  )
my_height <- 1 + length(unique(comb_d$cancer_n)) * 0.25
my_width <- 3 + 2 * length(unique(comb_d$gene))
my_ggsave(
  slug = glue("boxploth-CD8-{paste(unique(comb_d$gene), collapse = '-')}"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = my_width,
  height = my_height,
  plot = p
)

# }}}

# }}}

# }}}

# CD4 dataset files {{{

cd4_rds_files <- Sys.glob(glue("{root_dir}/data/expression/CD4/byDataset/*.sce.rds"))
cd4_rds_names <- str_remove(basename(cd4_rds_files), "\\.sce\\.rds$")
cd4_rds_data <- pblapply(cd4_rds_files, readRDS)
names(cd4_rds_data) <- cd4_rds_names

# What percent of Tregs have each of the genes we care about?
# Treg: IL10, TNFRSF9, CXCR3, TNFRSF4 (OX40)
treg_genes <- c(
  "TNFRSF4", "BATF", "TNFRSF18", "TNFRSF1B", "CTLA4", "IL2RA", "IL10", "FOXP3",
  "HIVEP1", "IRF4", "IKZF2", "CXCR3", "TNFRSF9", "LAG3"
)

gene_percents <- function(my_name, my_genes) {
  sce <- cd4_rds_data[[my_name]]
  sce_genes <- rowData(sce) %>% as.data.frame
  sce_meta <- colData(sce) %>% as.data.frame
  if ("counts" %in% names(assays(sce))) {
    sce_counts <- assay(sce, "counts")
  } else {
    sce_counts <- assay(sce, "norm_exprs")
  }
  if (!all(sce_genes$display.name == rownames(sce_counts))) {
    rownames(sce_counts) <- sce_genes$display.name
  }
  my_genes <- base::intersect(my_genes, rownames(sce_counts))
  for (my_gene in my_genes) {
    sce_meta[[glue("gene_{my_gene}")]] <- as.numeric(sce_counts[my_gene,])
  }
  sce_meta %>%
    select(cellID, patient, meta.cluster, starts_with("gene_")) %>%
    pivot_longer(starts_with("gene_"), names_to = "gene") %>%
    group_by(patient, meta.cluster, gene) %>%
    summarize(
      cells_total = n(),
      cells_pos = sum(value > 0),
      cells_freq = sum(value > 0) / n(),
      .groups = "keep"
    ) %>%
    mutate(dataset = my_name, gene = str_remove(gene, "gene_")) %>%
    relocate(dataset)
}
cd4_meta <- rbindlist(pblapply(cd4_rds_names, function(rds_name) {
  gene_percents(rds_name, treg_genes)
}))
cd4_meta %<>% mutate(cancer = str_split_fixed(dataset, "\\.", 2)[,1])

my_datasets <- (cd4_meta %>% count(dataset, patient) %>% count(dataset) %>% filter(n > 1))$dataset
my_d <- cd4_meta %>%
  filter(str_detect(meta.cluster, "Treg")) %>% # just the Tregs
  filter(dataset %in% my_datasets) %>%
  # filter(cells_total > 5) %>%
  group_by(cancer, patient, gene) %>%
  summarize(
    freq = sum(cells_pos) / sum(cells_total),
    cells = sum(cells_total),
    .groups = "keep"
  ) %>%
  mutate(log10freq = log10(freq * 100 + 1)) %>%
  ungroup()
#
my_d <- left_join(
  my_d,
  my_d %>%
    select(cancer, patient) %>%
    unique %>%
    count(cancer) %>%
    mutate(cancer_n = glue("{cancer} (n = {n})")) %>%
    select(-n),
  by = "cancer"
)
#
my_d$cancer_n <- factor(my_d$cancer_n, (my_d %>%
  group_by(cancer_n, gene) %>%
  summarize(mean = mean(freq)) %>%
  group_by(cancer_n) %>%
  summarize(mean = mean(mean)) %>%
  arrange(mean)
)$cancer_n)
#
my_d_mean <- my_d %>%
  group_by(cancer_n, gene) %>%
  summarize(mean = mean(freq)) %>%
  mutate(label = glue("{signif(100 * mean, 2)}%"))
#
p <- ggplot(my_d) +
  aes(y = cancer_n, x = freq, fill = cancer_n) +
  facet_grid(cols = vars(gene), scales = "free_x") +
  geom_stripes(odd = "#00000000", even = "#00000011") +
  # geom_boxplot(outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = c(0), linewidth = 0.3) +
  geom_point(position = position_quasirandom(), size = 0.5, color = "grey40") +
  geom_crossbar(
    data = my_d_mean, mapping = aes(x = mean, xmin = mean, xmax = mean, y = cancer_n),
    linewidth = 0.3, fatten = 1.0, color = "red"
  ) +
  geom_text(
    data = my_d_mean,
    mapping = aes(x = 0, y = cancer_n, label = label, alpha = mean),
    hjust = 1.15, size = 5
  ) +
  scale_alpha_continuous(guide = "none", trans = "log10") +
  scale_fill_manual(values = mpn65, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0.8, 0.2)),
    breaks = function(d) { d[1] <- 0; pretty_breaks()(d, 3) },
    labels = function(d) ifelse(d >= 0, signif(100 * d, 2), "")
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = glue("Percent of CD4 Tregs positive for each gene, per patient"),
    x = NULL,
    y = NULL,
    caption = "data from Zheng et al 2021"
  ) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.4, color = "#00000011"),
    panel.spacing = unit(0.5, "lines"),
    plot.title.position = "plot",
    plot.title = element_markdown(),
    strip.text.x = element_text(face = "italic")
    # strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1)
  )
my_height <- 1 + length(unique(my_d$cancer_n)) * 0.25
my_width <- 2 * length(treg_genes)
my_ggsave(
  slug = glue("boxploth-CD4-Treg-{paste(treg_genes, collapse = '-')}"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = my_width,
  height = my_height,
  plot = p
)

# Percent of CD4 T cells with Treg genes in our data {{{

h5_file <- "paper/tissue-cd4.h5"
our_obs <- h5read(h5_file, "obs")
our_counts <- sparseMatrix(
  i = as.numeric(h5read(h5_file, "matrix/indices")),
  p = as.numeric(h5read(h5_file, "matrix/indptr")),
  x = as.numeric(h5read(h5_file, "matrix/data")),
  dims = as.numeric(h5read(h5_file, "matrix/shape")),
  index1 = FALSE
)
rownames(our_counts) <- h5read(h5_file, "matrix/features")
colnames(our_counts) <- h5read(h5_file, "matrix/barcodes")
our_counts[1:5,1:5]

genes <- fread("data/ensembl_id-symbol.tsv")

for (my_sym in treg_genes) {
  our_ens <- genes$ensembl_id[genes$symbol == my_sym]
  our_obs[[glue("gene_{my_sym}")]] <- as.numeric(our_counts[our_ens,])
}

our_obs

our_d <- our_obs %>%
  filter(leiden %in% c(5, 8, 9)) %>%
  select(case, donor, starts_with("gene_")) %>%
  pivot_longer(starts_with("gene_"), names_to = "gene") %>%
  group_by(donor, case, gene) %>%
  summarize(
    freq = sum(value > 0) / n(),
    cells = n(),
    .groups = "keep"
  ) %>%
  mutate(log10freq = log10(freq * 100 + 1)) %>%
  ungroup() %>%
  mutate(
    gene = str_remove(gene, "gene_"),
    cancer = glue("irColitis {case}")
  ) %>%
  dplyr::rename(patient = donor) %>%
  relocate(cancer, patient, gene, freq) %>%
  select(-case)

our_d$cancer_n <- NULL
my_d$cancer_n <- NULL

our_d$data <- "Our data"
my_d$data <- "Zheng et al 2021"

head(our_d)
head(my_d)

comb_d <- rbind(our_d, my_d)
comb_d <- left_join(
  comb_d,
  comb_d %>%
    select(cancer, patient) %>%
    unique %>%
    count(cancer) %>%
    mutate(cancer_n = glue("{cancer} (n = {n})")) %>%
    select(-n),
  by = "cancer"
)
# comb_d <- comb_d %>%
#   filter(gene %in% c("IL10", "CXCR3", "TNFRSF9"))
comb_d_mean <- comb_d %>%
  group_by(data, cancer_n, gene) %>%
  summarize(mean = mean(freq), .groups = "keep") %>%
  ungroup() %>%
  mutate(label = glue("{signif(100 * mean, 2)}%"))
#
p <- ggplot(comb_d) +
  aes(y = cancer_n, x = freq, fill = cancer_n) +
  facet_grid(
    rows = vars(data), cols = vars(gene), scales = "free", space = "free_y",
    switch = "y"
  ) +
  geom_stripes(odd = "#00000000", even = "#00000011") +
  # geom_boxplot(outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = c(0), linewidth = 0.3) +
  geom_point(position = position_quasirandom(), size = 0.5, color = "grey40") +
  geom_crossbar(
    data = comb_d_mean, mapping = aes(x = mean, xmin = mean, xmax = mean, y = cancer_n),
    linewidth = 0.3, fatten = 1.0, color = "red"
  ) +
  geom_text(
    data = comb_d_mean,
    mapping = aes(x = 0, y = cancer_n, label = label, alpha = mean),
    hjust = 1.15, size = 5
  ) +
  scale_alpha_continuous(guide = "none", trans = "log10") +
  scale_fill_manual(values = mpn65, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0.8, 0.2)),
    breaks = function(d) { d[1] <- 0; pretty_breaks()(d, 3) },
    labels = function(d) ifelse(d >= 0, signif(100 * d, 2), "")
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = glue("Percent of CD4 Tregs positive for each gene, per patient"),
    x = NULL,
    y = NULL
  ) +
  theme(
    strip.text.y.left = element_text(angle = 0),
    panel.grid.major.x = element_line(linewidth = 0.4, color = "#00000011"),
    panel.spacing = unit(0.5, "lines"),
    plot.title.position = "plot",
    plot.title = element_markdown(),
    strip.text.x = element_text(face = "italic")
    # strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1)
  )
my_height <- 1 + length(unique(comb_d$cancer_n)) * 0.25
my_width <- 3 + 2 * length(unique(comb_d$gene))
my_ggsave(
  slug = glue("boxploth-CD4-Treg-{paste(unique(comb_d$gene), collapse = '-')}"),
  out_dir = file.path(out_dir, "freq", "boxplot"),
  type = "pdf",
  width = my_width,
  height = my_height,
  plot = p
)

# }}}


# }}}

