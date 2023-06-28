# libraries {{{
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(conflicted)
library(RColorBrewer)
library(dendsort)
library(data.table)
library(ggbeeswarm)
library(ggtext)
library(ggforestplot)
library(qs)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggstance)
library(glue)
library(harmony)
library(janitor)
library(limma)
library(magrittr)
library(naturalsort)
library(pals)
library(patchwork)
library(pbapply)
library(pheatmap)
library(purrr)
library(qs)
library(readxl)
library(rhdf5)
library(scales)
library(scattermore)
library(scico)
library(shadowtext)
library(sitools)
library(tidygraph)
library(tidyverse)
library(uwot)
library(fgsea)
#
# My functions
source("R/functions/helpers.R")
source("R/functions/do-analysis.R")
source("R/functions/theme-kamil.R")
#
source("R/plot-analysis.R")
theme_set(theme_kamil)
#
conflict_prefer("geom_errorbarh", "ggplot2")
conflict_prefer("discard", "purrr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("arrange", "dplyr")
#
source("R/functions/composition.R")
source("R/functions/do-pseudobulk-de.R")

# }}}

# Read data files {{{

data_dir <- "/projects/external_data/Smillie2019"
root_out_dir <- file.path(data_dir, "kamil")
dir.create(root_out_dir, showWarnings = FALSE, recursive = TRUE)

# Read the immune cell data (Imm)
my_type <- "Imm"
barcodes <- fread(glue("{data_dir}/{my_type}.barcodes2.tsv"), header = FALSE)$V1
genes <- fread(glue("{data_dir}/{my_type}.genes.tsv"), header = FALSE)$V1
m <- fread(glue("gzip -cd {data_dir}/gene_sorted-{my_type}.matrix.mtx.gz"),
  skip = 2, col.names = c("i", "j", "x"))
m_dims <- fread(
  glue("gzip -cd {data_dir}/gene_sorted-{my_type}.matrix.mtx.gz | head -n2"),
  skip = 1, col.names = c("rows", "columns", "values")
)
counts <- sparseMatrix(
  i = m$i, j = m$j, x = m$x,
  dims = c(m_dims$rows, m_dims$columns)
)
counts@Dimnames[[1]] <- genes
counts@Dimnames[[2]] <- barcodes
n_cells <- ncol(counts)

obs <- fread(file.path(data_dir, "all.meta2.txt.gz"))
obs <- obs[2:nrow(obs),]
obs$nGene <- as.numeric(obs$nGene)
obs$nUMI <- as.numeric(obs$nUMI)
colnames(obs) <- c(
  "cell", "cluster", "genes", "umis", "donor", "health", "location", "sample"
)
stopifnot(all(colnames(counts) %in% obs$cell))
obs <- obs[match(colnames(counts), obs$cell),]
stopifnot(all(obs$cell== colnames(counts)))

range(obs$genes)
range(obs$umis)

mito_genes <- genes[str_detect(genes, "^MT-")]
length(mito_genes)

counts[1:5,1:5]

# }}}

# Percent of CD8 T cells with IL17A in Smillie et al 2019 {{{

n_pos <- sum(counts["IL17A",] > 0) 
n_tot <- ncol(counts)

n_pos / n_tot

n_pos
n_tot

obs

obs$IL17A <- as.numeric(counts["IL17A",])

obs %>%
  count(health)

obs %>%
  filter(str_detect(cluster, "CD8")) %>%
  summarize(freq_IL17A = sum(IL17A > 0) / n())

obs %>%
  filter(str_detect(cluster, "CD8")) %>%
  group_by(health) %>%
  summarize(freq_IL17A = sum(IL17A > 0) / n())

out_dir <- "results/percent_IL17A"


obs_donor <- obs %>%
  filter(str_detect(cluster, "CD8")) %>%
  group_by(donor, health) %>%
  summarize(
    n_total = n(),
    n_IL17A = sum(IL17A > 0),
    freq_IL17A = sum(IL17A > 0) / n()
  )
obs_donor %<>% left_join(
  obs_donor %>%
  group_by(health) %>%
  summarize(median = median(freq_IL17A)) %>%
  mutate(health_freq = glue("{health} ({signif(100 * median, 2)}%)")),
  by = "health"
)

p_smillie <- ggplot(obs_donor) +
  aes(x = freq_IL17A, y = health_freq) +
  geom_stripes(odd = "#11111111", even = "#00000000") +
  geom_boxplot(linewidth = 0.3, outlier.shape = NA) +
  geom_point(position = position_quasirandom(), alpha = 0.5) +
  scale_y_discrete(
    position = "right"
  ) +
  scale_x_continuous(
    # trans = "log10",
    labels = \(x) signif(100 * x, 2)
  ) +
  theme(
    plot.title = element_markdown()
  ) +
  labs(
    title = "Percent of CD8 T cells with <i>IL17A</i>",
    x = "Percent",
    y = NULL,
    caption = "Smillie et al 2019"
  )
my_ggsave(
  slug = glue("Smillie2019"),
  out_dir = out_dir,
  type = "pdf",
  width = 6,
  height = 3,
  plot = p_smillie
)

# }}}

# Percent of CD8 T cells with IL17A in our data {{{

h5_file <- "paper/tissue-cd8.h5"
my_obs <- h5read(h5_file, "obs")
my_counts <- sparseMatrix(
  i = as.numeric(h5read(h5_file, "matrix/indices")),
  p = as.numeric(h5read(h5_file, "matrix/indptr")),
  x = as.numeric(h5read(h5_file, "matrix/data")),
  dims = as.numeric(h5read(h5_file, "matrix/shape")),
  index1 = FALSE
)
rownames(my_counts) <- h5read(h5_file, "matrix/features")
colnames(my_counts) <- h5read(h5_file, "matrix/barcodes")
my_counts[1:5,1:5]

genes <- fread("data/ensembl_id-symbol.tsv")

my_ens <- genes$ensembl_id[genes$symbol == "IL17A"]
my_obs$IL17A <- as.numeric(my_counts[my_ens,])

my_obs

my_obs_donor <- my_obs %>%
  group_by(donor, case) %>%
  summarize(freq_IL17A = sum(IL17A > 0) / n())
my_obs_donor %<>% left_join(
  my_obs_donor %>%
  group_by(case) %>%
  summarize(median = median(freq_IL17A)) %>%
  mutate(case_freq = glue("{case} ({signif(100 * median, 2)}%)")),
  by = "case"
)
p_ours <- ggplot(my_obs_donor) +
  aes(x = freq_IL17A, y = case_freq) +
  geom_stripes(odd = "#11111111", even = "#00000000") +
  geom_boxplot(linewidth = 0.3, outlier.shape = NA) +
  geom_point(position = position_quasirandom(), alpha = 0.5) +
  scale_y_discrete(
    position = "right"
  ) +
  scale_x_continuous(
    # trans = "log10",
    labels = \(x) signif(100 * x, 2)
  ) +
  theme(
    plot.title = element_markdown()
  ) +
  labs(
    title = "Percent of CD8 T cells with <i>IL17A</i>",
    x = "Percent",
    y = NULL,
    caption = "our data"
  )
my_ggsave(
  slug = glue("our-data"),
  out_dir = out_dir,
  type = "pdf",
  width = 6,
  height = 2.5,
  plot = p_ours
)

# }}}

# Percent of CD8 T cells with IL17A in Luoma et al 2021 {{{

h5_file <- "paper/luoma-cd8.h5"
lu_obs <- h5read(h5_file, "obs")
lu_counts <- sparseMatrix(
  i = as.numeric(h5read(h5_file, "matrix/indices")),
  p = as.numeric(h5read(h5_file, "matrix/indptr")),
  x = as.numeric(h5read(h5_file, "matrix/data")),
  dims = as.numeric(h5read(h5_file, "matrix/shape")),
  index1 = FALSE
)
rownames(lu_counts) <- h5read(h5_file, "matrix/features")
colnames(lu_counts) <- h5read(h5_file, "matrix/barcodes")
lu_counts[1:5,1:5]

genes <- fread("data/ensembl_id-symbol.tsv")

lu_ens <- genes$ensembl_id[genes$symbol == "IL17A"]
lu_obs$IL17A <- as.numeric(lu_counts[lu_ens,])

lu_obs

lu_obs_donor <- lu_obs %>%
  group_by(donor, case) %>%
  summarize(freq_IL17A = sum(IL17A > 0) / n())
lu_obs_donor %<>% left_join(
  lu_obs_donor %>%
  group_by(case) %>%
  summarize(median = median(freq_IL17A)) %>%
  mutate(case_freq = glue("{case} ({signif(100 * median, 2)}%)")),
  by = "case"
)
p_ours <- ggplot(lu_obs_donor) +
  aes(x = freq_IL17A, y = case_freq) +
  geom_stripes(odd = "#11111111", even = "#00000000") +
  geom_boxplot(linewidth = 0.3, outlier.shape = NA) +
  geom_point(position = position_quasirandom(), alpha = 0.5) +
  scale_y_discrete(
    position = "right"
  ) +
  scale_x_continuous(
    # trans = "log10",
    labels = \(x) signif(100 * x, 2)
  ) +
  theme(
    plot.title = element_markdown()
  ) +
  labs(
    title = "Percent of CD8 T cells with <i>IL17A</i>",
    x = "Percent",
    y = NULL,
    caption = "our data"
  )
my_ggsave(
  slug = glue("Luoma2020"),
  out_dir = out_dir,
  type = "pdf",
  width = 6,
  height = 2.5,
  plot = p_ours
)

# }}}


all_obs_donor <- rbind(
  rbind(
    obs_donor %>% mutate(data = "Smillie et al 2019") %>% rename(case = health),
    my_obs_donor %>% mutate(data = "Our Study")
  ),
  lu_obs_donor %>% mutate(data = "Luoma et al 2020")
) %>% ungroup()
#
all_obs_donor %<>% mutate(data_case = glue("{data}_{case}"))
data_case_map <- c(
  "Luoma et al 2020_Case" = "irColitis Case (n=6)",
  "Luoma et al 2020_Control" = "Immunotherapy Control (n=10)",
  "Our Study_Case" = "irColitis Case (n=13)",
  "Our Study_Control" = "Immunotherapy Control (n=14)",
  "Smillie et al 2019_Healthy" = "Healthy (n=12)",
  "Smillie et al 2019_Inflamed" = "Ulcerative Colitis Inflamed (n=18)",
  "Smillie et al 2019_Non-inflamed" = "Ulcerative Colitis Non-inflamed (n=18)"
)
#
all_obs_donor %<>% mutate(case = data_case_map[data_case])

my_median <- all_obs_donor %>%
  group_by(data, case) %>%
  summarize(mean = mean(freq_IL17A), median = median(freq_IL17A)) %>%
  mutate(label = glue("{signif(100 * mean, 2)}%"))
p_combined <- ggplot(all_obs_donor %>% mutate(freq_IL17A = ifelse(freq_IL17A > 0, freq_IL17A, 0.0001))) +
  aes(x = freq_IL17A, y = case) +
  geom_stripes(odd = "#11111111", even = "#00000000") +
  # geom_boxplot(linewidth = 0.3, outlier.shape = NA) +
  geom_crossbar(data = my_median, mapping = aes(x = mean, xmin = mean, xmax = mean, y = case)) +
  geom_point(
    mapping = aes(fill = freq_IL17A > 0.0001),
    shape = 21, size = 2, stroke = 0.3,
    position = position_quasirandom()
  ) +
  scale_fill_manual(guide = "none", values = c("white", "black")) +
  # facet_wrap(facets = vars(data), scales = "free_y", space = "free", ncol = 1) +
  facet_grid(rows = vars(data), scales = "free_y", space = "free", switch = "y") +
  geom_text(
    data = my_median,
    mapping = aes(x = 0.8, y = case, label = label),
    size = 5, hjust = 1
  ) +
  scale_y_discrete(
    position = "right"
  ) +
  scale_x_continuous(
    trans = "log10",
    labels = \(x) signif(100 * x, 2)
  ) +
  annotation_logticks(sides = "b", size = 0.3) +
  theme(
    plot.title.position = "plot",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    plot.title = element_markdown()
  ) +
  labs(
    title = "Percent of CD8 T cells with <i>IL17A</i>",
    x = "Percent",
    y = NULL,
    caption = "(for display, zeros were set to 0.01%)"
  )
my_ggsave(
  slug = glue("combined"),
  out_dir = out_dir,
  type = "pdf",
  width = 9,
  height = 5,
  plot = p_combined
)

