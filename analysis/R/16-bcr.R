library(Matrix)
library(RColorBrewer)
# library(cowplot)
library(dendsort)
library(data.table)
# library(doParallel)
# library(foreach)
library(ggbeeswarm)
library(qs)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggraph)
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
# library(princurve)
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
# library(waffle)
# Bioconductor
# library(monocle)
# library(SingleR)
# library(fgsea)
#
# My functions
source("R/functions/helpers.R")
source("R/functions/read-h5-files.R")
source("R/functions/write-mtx.R")
source("R/functions/do-analysis.R")
source("R/functions/write-de.R")
source("R/functions/theme-kamil.R")
#
source("R/plot-analysis.R")
source("R/colors-Luoma2020.R")
theme_set(theme_kamil)
#
discard <- purrr::discard
rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
#
source("R/functions/composition.R")
source("R/plot-composition.R")
source("R/functions/do-pseudobulk-de.R")
source("R/sample-ids.R")
#
rename_by_size <- function(x) {
  sizes <- sort(table(x), decreasing = TRUE)
  x_old_names <- names(sizes)
  x_new_names <- seq_along(x_old_names)
  names(x_new_names) <- x_old_names
  return(x_new_names[x])
}
#
recluster_cd8_leiden122 <- function(cd8_clusters) {
  # First round
  cd8_clusters[cd8_clusters %in% c("2", "6")] <- "2.6"
  cd8_clusters[cd8_clusters %in% c("1", "12")] <- "1.12"
  cd8_clusters[cd8_clusters %in% c("8", "9")] <- "8.9"
  #
  return(rename_by_size(cd8_clusters))
}
#
recluster_cd8_leiden151 <- function(cd8_clusters) {
  # First round
  cd8_clusters[cd8_clusters %in% c("1", "4", "13")] <- "1.4.13"
  cd8_clusters[cd8_clusters %in% c("10", "11", "15")] <- "10.11.15"
  cd8_clusters[cd8_clusters %in% c("2", "9")] <- "2.9"
  cd8_clusters[cd8_clusters %in% c("7", "16")] <- "7.16"
  #
  return(rename_by_size(cd8_clusters))
}

# Ensembl to gene symbol
########################################################################
my_h5_files <- Sys.glob(
  "analysis/terra/cellranger-per-channel/output/*/filtered_feature_bc_matrix.h5",
)
genes <- h5read(my_h5_files[1], "matrix/features")
genes <- tibble(ensembl_id = genes$id, symbol = genes$name)
ensembl_to_symbol <- unlist(with(genes, split(symbol, ensembl_id)))


# BCRdiv and alpha diversity
########################################################################

{

  analysis_name <- "a12_4_4_b5_1_3"
  params <- list(
    min_cells_in_cluster = 50,
    min_percent_of_cells_with_gene = 5
  )
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
  # my_leiden <- "leiden0.933"
  # my_leiden <- "leiden1.22"
  a1$obs$leiden <- a1$obs$leiden0.933
  a1$obs$cluster <- a1$obs$leiden
	#
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/composition-case-vs-control"))
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	#
	sample_info <- janitor::clean_names(read_excel(
			path = "data/luoma_villani_combined_colitis_samples2.xlsx",
			sheet = 3
	))
	donor_to_sex <- unlist(split(sample_info$sex, sample_info$donor))
  #
	a1$obs$sex <- donor_to_sex[a1$obs$donor]
	a1$obs$drug <- factor(a1$obs$drug, c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4"))
  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))

  vdj_files <- Sys.glob(
    "analysis/terra/cellranger-vdj_2020-07-27/output/*_BCR*/all_contig_annotations.csv"
  )
  length(vdj_files)
  source("R/functions/load-bcr.R")
  bcr <- load_bcr(vdj_files) # 30649
  bcr$cell <- str_replace(bcr$cell, "_BCR", "_GEX")
  n_bcr_cells <- length(intersect(bcr$cell, a1$obs$cell))
  stopifnot(n_bcr_cells > 0)

  a1$obs <- dplyr::left_join(
    x = a1$obs,
    y = bcr %>% select(-channel),
    by = "cell"
  )

  bcr_cols <- c("IGLC", "IGHC", "IGKC")
  if (all(bcr_cols %in% colnames(a1$obs))) {
    a1$obs$has_bcr <- with(a1$obs,
      !is.na(IGH_cdr3) & (
        !is.na(IGL_cdr3) | !is.na(IGK_cdr3)
      )
    )
  }
  table(a1$obs$has_bcr)

  # What percent of each channel's cells have a BCR?
  {
    x <- a1$obs %>%
      dplyr::count(has_bcr, channel) %>%
      dplyr::arrange(n)
    x$channel <- factor(x$channel, (
      x %>% dplyr::group_by(channel) %>%
        dplyr::summarize(n = sum(n)) %>%
        dplyr::arrange(n)
    )$channel)
    p1 <- ggplot(x) +
      aes(y = channel, x = n, group = has_bcr, fill = has_bcr) +
      geom_colh() +
      scale_fill_manual(
        name = "Has BCR:",
        values = c("grey80", "grey20"),
        labels = c("TRUE" = "T", "FALSE" = "F")
      ) +
      theme(
        legend.position = c(1, 0), legend.just = c(1, 0),
        # legend.position = "bottom",
        legend.background = element_blank()
      ) +
      scale_x_continuous(labels = scales::label_number_si()) +
      labs(
        x = "Cells", y = NULL
      )
    x2 <- a1$obs %>%
      group_by(channel) %>%
      summarize(
        has_bcr = sum(has_bcr),
        n = n(),
        pct = 100 * sum(has_bcr) / n()
      ) %>%
      ungroup() %>%
      mutate(channel = factor(channel))
    x2$channel <- factor(x2$channel, levels(x$channel))
    p2 <- ggplot(x2) +
      aes(y = channel, x = pct) +
      geom_colh(fill = "grey20") +
      labs(x = "Percent", y = NULL) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    # p1 + p2
    fig_height <- length(unique(x$channel)) * 0.28 + 0.5
    my_ggsave(
      "bcr-by-channel",
      out_dir = out_dir,
      type = "pdf",
      plot = p1 + p2 + plot_annotation(
        title = "Cells with BCR"
      ),
      scale = 1, width = 12, height = fig_height, units = "in", dpi = 300
    )
  }


  # What percent of each donor's cells have a BCR?
  {
    a1$obs$donor_pub <- pub_ids[a1$obs$donor]
    x <- a1$obs %>%
      dplyr::count(has_bcr, donor_pub) %>%
      dplyr::arrange(n)
    x$donor_pub <- factor(x$donor_pub, (
      x %>% dplyr::group_by(donor_pub) %>%
        dplyr::summarize(n = sum(n)) %>%
        dplyr::arrange(n)
    )$donor_pub)
    p1 <- ggplot(x) +
      aes(y = donor_pub, x = n, group = has_bcr, fill = has_bcr) +
      geom_colh() +
      scale_fill_manual(
        name = "Has BCR:",
        values = c("grey80", "grey20"),
        labels = c("TRUE" = "T", "FALSE" = "F")
      ) +
      theme(
        legend.position = c(1, 0), legend.just = c(1, 0),
        # legend.position = "bottom",
        legend.background = element_blank()
      ) +
      scale_x_continuous(labels = scales::label_number_si()) +
      labs(
        x = "Cells", y = NULL
      )
    x2 <- a1$obs %>%
      group_by(donor_pub) %>%
      summarize(
        has_bcr = sum(has_bcr),
        n = n(),
        pct = 100 * sum(has_bcr) / n()
      ) %>%
      ungroup() %>%
      mutate(donor_pub = factor(donor_pub))
    x2$donor_pub <- factor(x2$donor_pub, levels(x$donor_pub))
    p2 <- ggplot(x2) +
      aes(y = donor_pub, x = pct) +
      geom_colh(fill = "grey20") +
      labs(x = "Percent", y = NULL) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    # p1 + p2
    fig_height <- length(unique(x$donor_pub)) * 0.28 + 0.5
    my_ggsave(
      "bcr-by-donor_pub",
      out_dir = out_dir,
      type = "pdf",
      plot = p1 + p2 + plot_annotation(
        title = "Cells with BCR",
      ),
      scale = 1, width = 8, height = fig_height, units = "in", dpi = 300
    )
  }

  p <- plot_hexgene(
    x            = a1$obs$UMAP1,
    y            = a1$obs$UMAP2,
    z            = as.numeric(a1$obs$has_bcr),
    bins         = 91,
    palette      = "oslo",
    direction    = -1,
    italic       = FALSE,
    use_quantile = FALSE,
    text = TRUE
  ) +
  labs(title = "Cells with BCR")
  # theme(legend.position = "none")
  my_ggsave(
    "bcr-umap-has_bcr",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 4, height = 4.5, units = "in", dpi = 300
  )


  d <- a1$obs %>%
    dplyr::filter(has_bcr) %>%
    # dplyr::filter(donor %in% bcr_donors) %>%
    dplyr::group_by(donor) %>%
    dplyr::count(
      IGH_cdr3, IGL_cdr3, IGK_cdr3,
      case
    ) %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(
      rank = seq(length(n)),
      percent = 100 * n / sum(n),
      cumpercent = cumsum(100 * n / sum(n))
    )
  d_n <- d %>% select(donor, case) %>% unique %>%
    group_by(case) %>% count() %>%
    mutate(case_n = sprintf("%s (n = %s)", case, n)) %>%
    ungroup()
  d <- left_join(d, d_n %>% select(case, case_n), by = "case")
  # d <- inner_join(d, a1$obs %>% select(donor, case) %>% unique, by = "donor")
  #
  p1 <- ggplot(d) +
    aes(x = rank, y = cumpercent, group = donor, color = case_n) +
    geom_line(size = 1) +
    scale_color_manual(
      # values = my_colors$case
      # values = c("Case" = "#E69F00", "Control" = "#000000"),
      values = c("#E69F00", "#000000"),
      name = NULL
    ) +
    scale_x_log10() +
    # scale_y_log10() +
    annotation_logticks(sides = "b") +
    labs(x = "BCR Rank", y = "Percent", title = "Cum. percent of cells with BCR")
  my_ggsave(
    "bcr-rank-abundance",
    out_dir = out_dir,
    plot = p1,
    type = "pdf",
    scale = 0.75,
    width = 7,
    height = 4,
    units = "in", dpi = 300
  )
  #
  x <- a1$obs %>%
    filter(has_bcr) %>%
    # filter(donor %in% bcr_donors) %>%
    mutate(
      SAMPLE = donor,
      CLONE = paste(IGH_cdr3, IGL_cdr3, IGK_cdr3),
      donor = unname(pub_ids[donor])
    )
  # x <- bcr %>% filter(donor %in% unique(obs$donor)) %>%
  #   mutate(SAMPLE = donor, CLONE = TRB_cdr3_trim)

  # Very slow
  sample_curve_file <- glue("{out_dir}/alakazam-bcr-alpha-diversity.qs")
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
    qsave(sample_curve, sample_curve_file)
  } else {
    sample_curve <- qread(sample_curve_file)
  }

  #
  d_div <- sample_curve@diversity
  d_div <- left_join(d_div, a1$obs %>% select(donor, case) %>% unique, by = c("SAMPLE" = "donor"))
  #
  p2 <- ggplot(d_div) +
    aes(q, d, color = case, group = SAMPLE) +
    geom_line(size = 1) +
    scale_color_manual(
      # values = my_colors$case,
      values = c("Case" = "#E69F00", "Control" = "#000000"),
      name = NULL
    )
  # x <- left_join(x, a1$obs %>% select(donor, case) %>% unique, by = c("donor"))
  # x$donor_old <- x$donor
  # x$donor <- pub_ids[x$donor]
  d_div_counts <- x %>%
    group_by(donor, case) %>%
    summarise(n = length(CLONE), n_unique = length(unique(CLONE)))
  p3 <- ggplot(d_div_counts) +
    # geom_colh(
    #   mapping = aes(x = n, y = reorder(donor, n_unique), fill = case),
    #   color = "white", size = 0.5, alpha = 0.7
    # ) +
    geom_colh(
      mapping = aes(x = n_unique, y = reorder(donor, n_unique), fill = append_n(case)),
      color = "white", size = 0.5
    ) +
    scale_x_continuous(labels = comma) +
    # labs(title = "Number of BCRs per donor", x = "BCRs (Unique | Total)", y = NULL) +
    labs(title = "Number of BCRs per donor", x = "Unique BCRs", y = NULL) +
    scale_fill_manual(
      # values = my_colors$case,
      # values = c("Case" = "#E69F00", "Control" = "#000000"),
      values = c("#E69F00", "#000000"),
      name = NULL
    ) +
    theme(
      axis.text.y = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_blank()
    )
  # p <- (
  #     p1 + theme(legend.position = "none")
  #       # labs(title = "Rank abundance")
  #   ) + (
  #     p2 + labs(title = "Sample diversity") +
  #       theme(legend.position = "none")
  #   ) + p3
  p <- p3 + (
      p1 + theme(legend.position = "none")
        # labs(title = "Rank abundance")
    ) + (
      p2 + labs(title = "Diversity") +
        theme(legend.position = "none")
    )
  my_ggsave(
    "bcr-diversity",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.75, width = 14,
    # height = 6,
    height = 4,
    units = "in", dpi = 300
  )
