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
unname <- base::unname
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


# TCRdiv and alpha diversity
########################################################################

{

  # analysis_name <- "a12_4_4_t4_cd4_2_2"
  analysis_name <- "a12_4_4_t4_cd8_1_2"
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
  a1$obs$leiden <- recluster_cd8_leiden151(a1$obs$leiden1.51)
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
  ## Skip cell clusters that have too few cells.
  #exclude_clusters <- (
  #  a1$obs %>% count(cluster) %>% filter(n < params[["min_cells_in_cluster"]])
  #)$leiden
  ##
  #a1$obs <- a1$obs[!a1$obs$cluster %in% exclude_clusters,]
  #
  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))

  vdj_files <- Sys.glob(
    "analysis/terra/cellranger-vdj_2020-07-27/output/*_TCR*/all_contig_annotations.csv"
  )
  length(vdj_files)
  source("R/functions/load-tcr.R")
  tcr <- load_tcr(vdj_files) # 30649
  tcr$cell <- str_replace(tcr$cell, "_TCR", "_GEX")

  n_tcr_cells <- length(intersect(tcr$cell, a1$obs$cell))

  # a1$obs <- dplyr::left_join(
  #   x = a1$obs,
  #   y = tcr %>% select(-channel),
  #   by = "cell"
  # )

  # Write file for Neal.
  # xx <- a1$obs %>% select(cell, donor, case, n_features, mito_pct, UMAP1, UMAP2, cluster, starts_with("TR"))
  # fwrite(xx, file.path(out_dir, "obs_tcr.tsv.gz"), sep = "\t")

  #tcr_cols <- c("TRBV", "TRAV", "TRB_cdr3", "TRA_cdr3")
  #if (all(tcr_cols %in% colnames(a1$obs))) {
  #  a1$obs$has_tcr <- with(a1$obs,
  #    !is.na(TRB_cdr3) & !is.na(TRA_cdr3)
  #  )
  #}
  ##
  #bcr_cols <- c("IGLC", "IGHC", "IGKC")
  #if (all(bcr_cols %in% colnames(a1$obs))) {
  #  a1$obs$has_bcr <- with(a1$obs,
  #    !is.na(IGH_cdr3) & (
  #      !is.na(IGL_cdr3) | !is.na(IGK_cdr3)
  #    )
  #  )
  #}

  # What percent of each channel's cells have a TCR?
  {
    x <- a1$obs %>%
      dplyr::count(has_tcr, channel) %>%
      dplyr::arrange(n)
    x$channel <- factor(x$channel, (
      x %>% dplyr::group_by(channel) %>%
        dplyr::summarize(n = sum(n)) %>%
        dplyr::arrange(n)
    )$channel)
    p1 <- ggplot(x) +
      aes(y = channel, x = n, group = has_tcr, fill = has_tcr) +
      geom_colh() +
      scale_fill_manual(
        name = "Has TCR:",
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
        has_tcr = sum(has_tcr),
        n = n(),
        pct = 100 * sum(has_tcr) / n()
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
      "tcr-by-channel",
      out_dir = out_dir,
      type = "pdf",
      plot = p1 + p2 + plot_annotation(
        title = "Cells with TCR (alpha and beta CDR3)"
      ),
      scale = 1, width = 12, height = fig_height, units = "in", dpi = 300
    )
  }


  # What percent of each donor's cells have a TCR?
  {
    a1$obs$donor_pub <- pub_ids[a1$obs$donor]
    x <- a1$obs %>%
      dplyr::count(has_tcr, donor_pub) %>%
      dplyr::arrange(n)
    x$donor_pub <- factor(x$donor_pub, (
      x %>% dplyr::group_by(donor_pub) %>%
        dplyr::summarize(n = sum(n)) %>%
        dplyr::arrange(n)
    )$donor_pub)
    p1 <- ggplot(x) +
      aes(y = donor_pub, x = n, group = has_tcr, fill = has_tcr) +
      geom_colh() +
      scale_fill_manual(
        name = "Has TCR:",
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
        has_tcr = sum(has_tcr),
        n = n(),
        pct = 100 * sum(has_tcr) / n()
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
      "tcr-by-donor_pub",
      out_dir = out_dir,
      type = "pdf",
      plot = p1 + p2 + plot_annotation(
        title = "Cells with TCR (alpha and beta CDR3)"
      ),
      scale = 1, width = 8, height = fig_height, units = "in", dpi = 300
    )
  }

  p <- plot_hexgene(
    x            = a1$obs$UMAP1,
    y            = a1$obs$UMAP2,
    z            = as.numeric(a1$obs$has_tcr),
    bins         = 91,
    palette      = "oslo",
    direction    = -1,
    italic       = FALSE,
    use_quantile = FALSE,
    text = TRUE
  ) +
  labs(title = "Cells with TCR")
  # theme(legend.position = "none")
  my_ggsave(
    "tcr-umap-has_tcr",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 4, height = 4.5, units = "in", dpi = 300
  )


  d <- a1$obs %>%
    dplyr::filter(has_tcr) %>%
    # dplyr::filter(donor %in% tcr_donors) %>%
    dplyr::group_by(donor) %>%
    dplyr::count(
      TRBV,
      TRAV,
      TRB_cdr3_trim,
      TRA_cdr3_trim,
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
    labs(x = "TCR Rank", y = "Percent", title = "Cum. percent of cells with TCR")
  my_ggsave(
    "tcr-rank-abundance",
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
    filter(has_tcr) %>%
    # filter(donor %in% tcr_donors) %>%
    mutate(
      SAMPLE = donor,
      CLONE = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim),
      donor = unname(pub_ids[donor])
    )
  # x <- tcr %>% filter(donor %in% unique(obs$donor)) %>%
  #   mutate(SAMPLE = donor, CLONE = TRB_cdr3_trim)

  my_donors <- a1$obs$donor %>% unique
  tcr_metrics <- rbindlist(lapply(Sys.glob(
    "analysis/terra/cellranger-vdj_2020-07-27/output/*_TCR*/metrics_summary.csv"
  ), function(this_file) {
    this_donor <- this_file %>% dirname %>% basename
    this_donor <- str_split_fixed(this_donor, "_", 3)[1,1:2] %>% paste(collapse = "_")
    if (this_donor %in% my_donors) {
      retval <- fread(this_file)
      retval$donor <- this_donor
      retval <- janitor::clean_names(retval)
      return(retval)
    }
    return(NULL)
  }), fill = TRUE)
  tcr_metrics$read_pairs <- parse_number(tcr_metrics$number_of_read_pairs)
    
  x <- d %>%
    mutate(p = percent / 100) %>%
    group_by(donor, case) %>%
    summarize(
      entropy = -sum(p * log(p)),
      n = n(),
      gsi = 1 - sum(p * p),
      .groups = "drop"
    ) %>%
    arrange(-entropy) %>%
    mutate(clonality = 1 - entropy / log(n)) %>% arrange(clonality)
  # x <- left_join(x, tcr_metrics %>% select(donor, read_pairs), by = "donor")
  # p <- ggplot(x) +
  #   aes(x = read_pairs, y = n, fill = case) +
  #   geom_point(shape = 21, size = 3, stroke = 0.3) +
  #   # scale_x_log10() +
  #   scale_x_continuous(labels = scales::label_number_si()) +
  #   # annotation_logticks(sides = "b") +
  #   scale_fill_manual(
  #     values = c(Case = "#E69F00", Control = "#000000"),
  #     name = NULL
  #   )
  #   # labs(x = "Number of TCRs", y = "Entropy")
  # my_ggsave(
  #   "tcr-read_pairs-vs-n",
  #   out_dir = out_dir,
  #   plot = p,
  #   type = "pdf",
  #   scale = 0.75,
  #   width = 7,
  #   height = 4,
  #   units = "in", dpi = 300
  # )
  my_pval <- wilcox.test(clonality ~ case, x, exact = FALSE)$p.value
  p <- ggplot(x) +
    aes(y = case, x = clonality, fill = case) +
    geom_boxplot(outlier.shape = NA, size = 0.3, alpha = 0.3) +
    geom_point(shape = 21, size = 3, stroke = 0.3, position = position_quasirandom(groupOnX = FALSE)) +
    scale_fill_manual(
      values = c(Case = "#E69F00", Control = "#000000"),
      name = NULL
    ) +
    labs(
      x = NULL, y = NULL, title = "Clonality",
      subtitle = glue("Wilcoxon rank sum test P = {str_replace(format.pval(my_pval, 1), '-0', '-')}")
    ) +
    theme(legend.position = "none")
  my_ggsave(
    "tcr-clonality",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.75,
    width = 5,
    height = 2.5,
    units = "in", dpi = 300
  )

  # d_pie <- d %>% mutate(tcr = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim)) %>%
  #   select(donor, case, tcr, n, percent) %>%
  #   group_by(donor) %>%
  #   arrange(-n) %>%
  #   mutate(id = seq(n()))
  # d_pie$donor <- factor(d_pie$donor, (x %>% arrange(clonality))$donor)
  # p <- ggplot(d_pie) +
  #   aes(x = 0, y = percent, fill = id) +
  #   geom_bar(width = 1, stat = "identity") +
  #   scale_fill_gradientn(colors = scico::scico(100)) +
  #   coord_polar("y", start = 0) +
  #   facet_grid(rows = vars(donor), cols = vars(case)) +
  #   theme(
  #     panel.spacing = unit(0.5, "lines"),
  #     legend.position = "none",
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank()
  #   )
  # my_ggsave(
  #   "tcr-pies",
  #   out_dir = out_dir,
  #   plot = p,
  #   type = "pdf",
  #   scale = 1,
  #   width = 5,
  #   height = 15,
  #   units = "in", dpi = 300
  # )



  # Very slow
  # sample_curve_file <- glue("{out_dir}/alakazam-tcr-alpha-diversity.qs")
  sample_curve_file <- glue("{out_dir}/alakazam-tcr-alpha-diversity2.qs")
  if (!file.exists(sample_curve_file)) {
    x <- a1$obs %>%
      filter(has_tcr) %>%
      mutate(
        SAMPLE = case,
        CLONE = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim)
      )
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
    # labs(title = "Number of TCRs per donor", x = "TCRs (Unique | Total)", y = NULL) +
    labs(title = "Number of TCRs per donor", x = "Unique TCRs", y = NULL) +
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
      p2 + labs(title = "Sample diversity") +
        theme(legend.position = "none")
    )
  my_ggsave(
    "tcr-diversity",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.75, width = 14,
    # height = 6,
    height = 4,
    units = "in", dpi = 300
  )

  d <- a1$obs %>%
    dplyr::filter(has_tcr) %>%
    # dplyr::filter(donor %in% tcr_donors) %>%
    dplyr::group_by(donor) %>%
    dplyr::count(
      TRBV,
      TRAV,
      TRB_cdr3_trim,
      TRA_cdr3_trim,
      case
    ) %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(
      rank = seq(length(n)),
      percent = 100 * n / sum(n)
    ) %>%
    dplyr::mutate(
      # tcr_q = cut(rank, breaks = 5, labels = FALSE)
      tcr_q = cut(rank, breaks = 10, labels = FALSE)
    ) %>%
    dplyr::mutate(tcr_clone = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim)) %>%
    ungroup()
  #
  top_tcr_clones <- d$tcr_clone[d$rank <= 1]
  #
  a1$obs <- a1$obs %>%
    dplyr::mutate(tcr_clone = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim))
  a1$obs$tcr_top_clone <- a1$obs$tcr_clone %in% top_tcr_clones

  p <- a1$obs %>%
    group_by(leiden, case) %>%
    summarize(pct_tcr_top = 100 * sum(tcr_top_clone) / n()) %>%
    ungroup() %>%
    mutate(
      leiden = naturalfactor(as.character(leiden)),
      case = factor(case, c("Control", "Case"))
    ) %>%
  ggplot() +
  aes(x = pct_tcr_top, y = leiden, group = case, fill = case) +
  scale_fill_manual(
    values = c("Case" = "#E69F00", "Control" = "#000000"),
    name = NULL
  ) +
  geom_colh(size = 0.2, position = position_dodgev()) +
  theme(
    legend.position = "none"
  )
  my_ggsave(
    "tcr-diversity-bars",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.75, width = 4,
    # height = 6,
    height = 5,
    units = "in", dpi = 300
  )

  dd <- left_join(
    x = a1$obs,
    y = d %>% select(donor, tcr_rank = rank, tcr_clone, tcr_q, tcr_percent = percent, tcr_n = n),
    by = c("donor", "tcr_clone")
  )
    # p_tcr_quintiles <- plot_scattermore(
    #   x = dd$UMAP1,
    #   y = dd$UMAP2,
    #   group = dd$tcr_q,
    #   group_colors = cluster_colors,
    #   pixels = 1000,
    #   alpha = 0.35
    # )
  ix <- !is.na(dd$tcr_q)
    p_tcr_quintiles <- ggplot(dd[ix,]) +
      scattermore::geom_scattermore(
        mapping = aes(UMAP1, UMAP2, color = factor(tcr_q)),
        # color     = grDevices::adjustcolor(mpn65, alpha.f = 0.35)[dd$tcr_q[ix]],
        pointsize = 4,
        pixels    = c(1000, 1000)
      ) +
      facet_wrap(~ append_n(case)) +
      scale_color_manual(
        # values = scico::scico(n = 6, palette = "oslo")[1:5]
        # values = RColorBrewer::brewer.pal(n = 5, name = "Spectral"),
        values = RColorBrewer::brewer.pal(n = 10, name = "Spectral"),
        name = NULL
      ) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      theme(
        # plot.title = element_text(size = 32),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )
    my_ggsave(
      "umap-tcr-quintiles",
      out_dir = out_dir,
      plot = p_tcr_quintiles,
      type = "pdf",
      scale = 1, width = 9, height = 4, units = "in", dpi = 300
    )

    ix <- !is.na(dd$tcr_q)
    p <- plot_hexgene(
      x            = dd$UMAP1[ix],
      y            = dd$UMAP2[ix],
      # z            = as.numeric(dd$tcr_q == 1)[ix],
      z            = as.numeric(dd$tcr_rank <= 10)[ix],
      group        = append_n(dd$case[ix]),
      bins         = 91,
      palette      = "oslo",
      direction    = -1,
      italic       = FALSE,
      use_quantile = FALSE,
      text = TRUE
    ) +
    facet_wrap(~ group) +
    labs(title = "Cells with top 10 expanded TCRs per donor")
    # theme(legend.position = "none")
    my_ggsave(
      "umap-tcr-top-quintile",
      out_dir = out_dir,
      plot = p,
      type = "pdf",
      scale = 1, width = 6, height = 4.5, units = "in", dpi = 300
    )

  # What pairs of clusters share the same TCRs?

  cluster_pairs <- combn(unique(dd$cluster), 2)
  tcr_pairs <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
    c1 <- cluster_pairs[1,i]
    c2 <- cluster_pairs[2,i]
    message(sprintf("%s %s", c1, c2))
    x <- dd %>%
      filter(has_tcr) %>%
      # filter(tcr_rank <= 10) %>%
      group_by(tcr_clone, case) %>%
      summarize(
        both_cluster = factor(all(c(c1, c2) %in% cluster)),
        expanded = factor(tcr_rank <= 10),
        .groups = "keep"
      ) %>% unique
    tab1 <- with(
      x %>% filter(case == "Case"),
      table(both_cluster, expanded)
    )
    tab2 <- with(
      x %>% filter(case == "Control"),
      table(both_cluster, expanded)
    )
    retval1 <- broom::tidy(stats::fisher.test(tab1))
    retval1$case <- "Case"
    retval1$c1 <- c1
    retval1$c2 <- c2
    retval2 <- broom::tidy(stats::fisher.test(with(
      x %>% filter(case == "Control"),
      table(both_cluster, expanded)
    )))
    retval2$case <- "Control"
    retval2$c1 <- c1
    retval2$c2 <- c2
    rbind(retval1, retval2)
  }))

  x <- tcr_pairs %>%
    select(-method, -alternative) %>%
    mutate(logp = -log10(p.value)) %>%
    group_by(c1, c2) %>%
    mutate(max_p = max(logp), diffp = diff(logp)) %>%
    ungroup() %>%
    filter(max_p > 4) %>%
    # arrange(-max_p) %>%
    arrange(-diffp) %>%
    mutate(pair = paste(c1, c2))
  x$pair <- factor(x$pair)
  p <- ggplot(x) +
    aes(x = logp, y = reorder(pair, -diffp), fill = case, group = case) +
    geom_colh(position = position_dodgev()) +
    scale_fill_manual(values = pals::okabe()[2:1]) +
    labs(x = bquote("-Log"[10]~"P"), y = NULL)
  my_ggsave(
    "tcr-cluster-pairs-bars",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 4, height = 11, units = "in", dpi = 300
  )

  x <- inner_join(
    tcr_pairs %>% filter(case == "Control"),
    tcr_pairs %>% filter(case == "Case"),
    by = c("c1", "c2")
  ) %>%
  mutate(log10p.x = -log10(p.value.x), log10p.y = -log10(p.value.y))
  x$label <- paste(x$c1, x$c2)
  x$label[abs(-log10(x$p.value.x) + log10(x$p.value.y)) < 5] <- ""
  p <- ggplot(x) +
    aes(x = -log10(p.value.x), y = -log10(p.value.y), label = label) +
    geom_abline(size = 0.2) +
    geom_point() +
    geom_text_repel()
  my_ggsave(
    "tcr-cluster-pairs",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 6, height = 4.5, units = "in", dpi = 300
  )

  cluster_pairs <- combn(unique(dd$cluster), 2)
  tcr_pairs2 <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
    c1 <- cluster_pairs[1,i]
    c2 <- cluster_pairs[2,i]
    message(sprintf("%s %s", c1, c2))
    x <- dd %>%
      filter(has_tcr) %>%
      # filter(tcr_rank <= 10) %>%
      # group_by(tcr_clone, case) %>%
      group_by(tcr_clone) %>%
      summarize(
        both_cluster = factor(all(c(c1, c2) %in% cluster), c(FALSE, TRUE)),
        is_case = factor(case == "Case", c(FALSE, TRUE)),
        .groups = "keep"
      ) %>% unique
    tab <- with(x, table(both_cluster, is_case))
    retval <- broom::tidy(stats::fisher.test(tab))
    retval$c1 <- c1
    retval$c2 <- c2
    retval
  }))

  tcr_pairs2 %>% select(-method, -alternative) %>% arrange(p.value) %>% head(15)

  # Look for the 6 2 TCRs, show UMAP for cases and UMAP for controls
  my_clones <- (
    dd %>%
    filter(has_tcr) %>%
    filter(tcr_rank <= 10) %>%
    group_by(tcr_clone) %>%
    summarize(both_cluster = all(c(6, 2) %in% cluster)) %>%
    filter(both_cluster)
  )$tcr_clone %>% unique

  p <- plot_hexgene(
    x            = dd$UMAP1,
    y            = dd$UMAP2,
    z            = dd$tcr_clone %in% my_clones,
    group        = append_n(dd$case),
    bins         = 71,
    palette      = "bilbao",
    direction    = 1,
    italic       = FALSE,
    use_quantile = FALSE,
    text = TRUE
  ) +
  facet_wrap(~ group) +
  labs(title = "Cells with TCR") +
  theme(legend.position = "none")
  my_ggsave(
    slug = glue("umap-clusterpair-6-2"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.8,
    # width = 3.5,
    width = 7,
    height = 3.5, units = "in", dpi = 300
  )

  # What pairs of clusters share the same TCRs? Permutation test.
  cluster_pairs <- combn(unique(dd$cluster), 2)
  x <- dd %>%
    filter(has_tcr, tcr_percent >= 1, tcr_n >= 5) %>%
    select(donor, cluster, tcr_clone) %>%
    as.data.table
  tcr_stat_donor <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
    c1 <- cluster_pairs[1,i]
    c2 <- cluster_pairs[2,i]
    rbindlist(lapply(unique(x$donor), function(my_donor) {
      x3 <- x[donor == my_donor]
      x1 <- unique(x3[cluster == c1]$tcr_clone)
      x2 <- unique(x3[cluster == c2]$tcr_clone)
      s1 <- 0
      if (length(x1) || length(x2)) {
        s1 <- length(intersect(x1, x2)) / length(unique(x3$tcr_clone))
      }
      data.frame(
        c1 = c1, c2 = c2, donor = my_donor, stat = s1,
        n1 = length(x1), n2 = length(x2),
        n_intersect = length(intersect(x1, x2)) ,
        n_donor = length(unique(x3$tcr_clone))
      )
    }))
  }))
  tcr_stat_donor <- left_join(
    x = tcr_stat_donor,
    y = dd %>% select(donor, case) %>% unique,
    by = "donor"
  )
  fwrite(tcr_stat_donor, file.path(out_dir, "tcr_stat_donor.tsv"), sep = "\t")
  tcr_stat <- tcr_stat_donor %>%
    group_by(c1, c2) %>%
    summarize(stat = mean(stat), .groups = "drop")

  x_null <- x
  i <- 0
  tcr_stat_null <- replicate(n = 200, simplify = FALSE, expr = {
    i <<- i + 1
    message(i)
    rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
      c1 <- cluster_pairs[1,i]
      c2 <- cluster_pairs[2,i]
      rbindlist(lapply(unique(x_null$donor), function(my_donor) {
        x3 <- x_null[donor == my_donor]
        x3$cluster_null <- sample(x3$cluster)
        x1 <- unique(x3[cluster_null == c1]$tcr_clone)
        x2 <- unique(x3[cluster_null == c2]$tcr_clone)
        s1 <- 0
        if (length(x1) || length(x2)) {
          s1 <- length(intersect(x1, x2)) / length(unique(x3$tcr_clone))
        }
        data.frame(c1 = c1, c2 = c2, donor = my_donor, stat = s1)
      }))
    }))
  })
  for (i in seq_along(tcr_stat_null)) {
    tcr_stat_null[[i]]$rep <- i
    tcr_stat_null[[i]] <- tcr_stat_null[[i]] %>%
      group_by(c1, c2) %>% summarize(stat = mean(stat), .groups = "drop")
  }
  tcr_stat_null <- rbindlist(tcr_stat_null)

  tcr_stat$pval <- 1
  for (i in seq(nrow(tcr_stat))) {
    my_c1 <- tcr_stat$c1[i]
    my_c2 <- tcr_stat$c2[i]
    my_stat <- tcr_stat$stat[i]
    stat_null <- (tcr_stat_null %>% filter(c1 == my_c1, c2 == my_c2))$stat
    tcr_stat$pval[i] <- (1 + sum(stat_null >= my_stat)) / (1 + length(stat_null))
  }
  tcr_stat$fdr <- p.adjust(tcr_stat$pval, method = "fdr")

  tcr_stat %>% arrange(pval) %>% head(15)

  p <- gg_qqplot(tcr_stat$pval)
  my_ggsave(
    slug = glue("tcr_stat-qqplot2"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 4,
    height = 3.5, units = "in", dpi = 300
  )

  my_pval <- (
  tcr_stat_donor %>% group_by(c1, c2) %>% summarize(pval = wilcox.test(stat ~ case)$p.value) %>% arrange(pval)
  )$pval
  p <- gg_qqplot(my_pval[!is.na(my_pval)])
  my_ggsave(
    slug = glue("tcr_stat-qqplot3"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 4,
    height = 3.5, units = "in", dpi = 300
  )



  # What pairs of clusters share the same TCRs? Is it different for cases and controls?
  cluster_pairs <- combn(unique(dd$cluster), 2)
  tcr_stat <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
    c1 <- cluster_pairs[1,i]
    c2 <- cluster_pairs[2,i]
    rbindlist(lapply(c("Case", "Control"), function(my_case) {
      x1 <- (
        dd %>%
        # filter(has_tcr, tcr_rank <= 10, case == my_case, cluster == c1)
        filter(has_tcr, tcr_percent >= 1, tcr_n >= 5, case == my_case, cluster == c1)
      )$tcr_clone
      x2 <- (
        dd %>%
        # filter(has_tcr, tcr_rank <= 10, case == my_case, cluster == c2)
        filter(has_tcr, tcr_percent >= 1, tcr_n >= 5, case == my_case, cluster == c2)
      )$tcr_clone
      x3 <- (
        dd %>%
        # filter(has_tcr, tcr_rank <= 10, case == my_case)
        filter(has_tcr, tcr_percent >= 1, tcr_n >= 5, case == my_case)
      )$tcr_clone
      s1 <- length(intersect(x1, x2)) / length(x3)
      data.frame(c1 = c1, c2 = c2, case = my_case, stat = s1)
    }))
  }))
  tcr_stat <- tcr_stat %>% pivot_wider(names_from = "case", values_from = "stat")
  tcr_stat <- tcr_stat %>% mutate(diff = Case - Control)

  # tcr_stat_orig <- tcr_stat

  donor_null <- dd %>% select(donor, case_null = case) %>% unique
  i <- 0
  tcr_stat_null <- replicate(n = 500, expr = {
    i <<- i + 1
    message(i)
    donor_null$case_null <- sample(donor_null$case_null)
    dd$case_null <- NULL
    dd <- left_join(dd, donor_null, by = "donor")
    rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
      c1 <- cluster_pairs[1,i]
      c2 <- cluster_pairs[2,i]
      rbindlist(lapply(c("Case", "Control"), function(my_case) {
        x1 <- (
          dd %>%
          # filter(has_tcr, tcr_rank <= 10, case_null == my_case, cluster == c1)
          filter(has_tcr, tcr_percent >= 1, tcr_n >= 5, case_null == my_case, cluster == c1)
        )$tcr_clone
        x2 <- (
          dd %>%
          # filter(has_tcr, tcr_rank <= 10, case_null == my_case, cluster == c2)
          filter(has_tcr, tcr_percent >= 1, tcr_n >= 5, case_null == my_case, cluster == c2)
        )$tcr_clone
        x3 <- (
          dd %>%
          # filter(has_tcr, tcr_rank <= 10, case_null == my_case)
          filter(has_tcr, tcr_percent >= 1, tcr_n >= 5, case_null == my_case)
        )$tcr_clone
        s1 <- length(intersect(x1, x2)) / length(x3)
        data.frame(c1 = c1, c2 = c2, case = my_case, stat = s1)
      }))
    }))
  }, simplify = FALSE)
  for (i in seq_along(tcr_stat_null)) {
    tcr_stat_null[[i]]$rep <- i
    tcr_stat_null[[i]] <- tcr_stat_null[[i]] %>%
      pivot_wider(names_from = "case", values_from = "stat") %>%
      mutate(diff = Case - Control)
  }
  tcr_stat_null <- rbindlist(tcr_stat_null)

  tcr_stat$pval <- 1
  for (i in seq(nrow(tcr_stat))) {
    my_c1 <- tcr_stat$c1[i]
    my_c2 <- tcr_stat$c2[i]
    # my_case <- tcr_stat$case[i]
    # my_stat <- tcr_stat$stat[i]
    # tcr_stat$pval[i] <- (
    #   tcr_stat_null %>% filter(case == my_case, c1 == my_c1, c2 == my_c2) %>%
    #     summarize(n = (1 + sum(stat >= my_stat)) / (1 + n()))
    # )[1,1]
    my_diff <- tcr_stat$diff[i]
    diff_null <- (tcr_stat_null %>% filter(c1 == my_c1, c2 == my_c2))$diff
    if (my_diff > 0) {
      tcr_stat$pval[i] <- (1 + sum(diff_null >= my_diff)) / (1 + length(diff_null))
    } else {
      tcr_stat$pval[i] <- (1 + sum(diff_null <= my_diff)) / (1 + length(diff_null))
    }
  }
  tcr_stat$fdr <- p.adjust(tcr_stat$pval, method = "fdr")

  tcr_stat %>% arrange(pval) %>% head(15)

  p <- gg_qqplot(tcr_stat$pval)
  my_ggsave(
    slug = glue("heatmap-tcr_stat-qqplot"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 4,
    height = 3.5, units = "in", dpi = 300
  )

  p <- tcr_stat %>% pivot_longer(cols = c("Case", "Control")) %>%
    group_by(c1, c2) %>%
    mutate(c3 = min(c1, c2), c4 = max(c1, c2)) %>%
    ungroup() %>%
    ggplot() +
    aes(x = factor(c3), y = factor(c4), fill = value) +
    geom_tile() +
    facet_wrap(~ name, ncol = 2) +
    scale_fill_gradientn(colors = scico::scico(n = 20, pal = "bilbao")) +
    labs(x = NULL, y = NULL)
  my_ggsave(
    slug = glue("heatmap-tcr_stat"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 8,
    height = 3.5, units = "in", dpi = 300
  )

  tcr_stat %>% filter(c1 == 4, c2 == 7)

  # p <- tcr_stat_null %>% filter(c1 == 4, c2 == 2) %>%
  #   pivot_longer(cols = c("Case", "Control")) %>%
  #   ggplot() +
  #   aes(x = value, group = name, fill = name) +
  #   geom_histogram(bins = 30, position = "identity", alpha = 0.5) +
  #   labs(x = NULL, y = NULL)
  p <- tcr_stat_null %>% filter(c1 == 4, c2 == 7) %>%
    ggplot() +
    aes(x = diff) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    geom_vline(xintercept = (tcr_stat %>% filter(c1 == 4, c2 == 7))$diff) +
    labs(x = NULL, y = NULL, title = "4 7")
  my_ggsave(
    slug = glue("heatmap-tcr_stat-4-7"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 8,
    height = 3.5, units = "in", dpi = 300
  )


  # What pairs of clusters share the same TCRs?
  cluster_pairs <- combn(unique(dd$cluster), 2)
  tcr_stat2 <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
    c1 <- cluster_pairs[1,i]
    c2 <- cluster_pairs[2,i]
    rbindlist(lapply(unique(dd$donor), function(my_donor) {
      x1 <- (
        dd %>%
        filter(has_tcr, tcr_rank <= 10, donor == my_donor, cluster == c1)
      )$tcr_clone
      x2 <- (
        dd %>%
        filter(has_tcr, tcr_rank <= 10, donor == my_donor, cluster == c2)
      )$tcr_clone
      x3 <- (
        dd %>%
        filter(has_tcr, tcr_rank <= 10, donor == my_donor)
      )$tcr_clone
      s1 <- length(intersect(x1, x2)) / length(x3)
      data.frame(c1 = c1, c2 = c2, donor = my_donor, stat = s1)
    }))
  }))
  # tcr_stat2$stat[is.na(tcr_stat2$stat)] <- 0
  tcr_stat2 <- left_join(tcr_stat2, dd %>% select(donor, case) %>% unique, by = "donor")

  tcr_stat2_pval <- rbindlist(lapply(seq(ncol(cluster_pairs)), function(i) {
    my_c1 <- cluster_pairs[1,i]
    my_c2 <- cluster_pairs[2,i]
    wt <- wilcox.test(formula = stat ~ case, data = tcr_stat2[c1 == my_c1 & c2 == my_c2])
    data.frame(c1 = my_c1, c2 = my_c2, pval = wt$p.value)
  }))
  tcr_stat2_diff <- tcr_stat2 %>% group_by(c1, c2) %>%
    summarize(diff = mean(stat[case == "Case"], na.rm = TRUE) - mean(stat[case == "Control"], na.rm = TRUE))
  tcr_stat2_pval <- left_join(tcr_stat2_pval, tcr_stat2_diff, by = c("c1", "c2"))

  p <- tcr_stat2 %>% filter(c1 == 4, c2 == 7) %>%
    ggplot() +
    aes(x = case, y = stat, group = case) +
    geom_quasirandom(groupOnX = TRUE) +
    labs(x = NULL, y = NULL, title = "4 7",
      subtitle = sprintf("p = %s", signif(
          (tcr_stat2_pval %>% filter(c1 == 4, c2 == 7))$pval
      ))
    )
  my_ggsave(
    slug = glue("heatmap-tcr_stat2-4-7"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 3,
    height = 3.5, units = "in", dpi = 300
  )
  p <- tcr_stat2 %>% filter(c1 == 6, c2 == 2) %>%
    ggplot() +
    aes(x = case, y = stat, group = case) +
    geom_quasirandom(groupOnX = TRUE) +
    labs(x = NULL, y = NULL, title = "6 2",
      subtitle = sprintf("p = %s", signif(
          (tcr_stat2_pval %>% filter(c1 == 6, c2 == 2))$pval
      ))
    )
  my_ggsave(
    slug = glue("heatmap-tcr_stat2-6-2"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 3,
    height = 3.5, units = "in", dpi = 300
  )


  p <- plot_hexgene(
    x            = a1$obs$UMAP1,
    y            = a1$obs$UMAP2,
    z            = as.numeric(a1$obs$has_tcr),
    group        = append_n(a1$obs$case),
    bins         = 101,
    palette      = "davos",
    direction    = -1,
    italic       = FALSE,
    use_quantile = FALSE,
    text = TRUE
  ) +
  facet_wrap(~ group) +
  labs(title = "Cells with TCR") +
  theme(legend.position = "none")
  my_ggsave(
    slug = glue("umap-has_tcr"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.8,
    # width = 3.5,
    width = 7,
    height = 3.5, units = "in", dpi = 300
  )
  p <- plot_hexgene(
    x            = a1$obs$UMAP1,
    y            = a1$obs$UMAP2,
    z            = as.numeric(a1$obs$tcr_top_clone),
    group        = append_n(a1$obs$case),
    bins         = 101,
    palette      = "davos",
    direction    = -1,
    italic       = FALSE,
    use_quantile = FALSE,
    text = TRUE
  ) +
  facet_wrap(~ group) +
  labs(title = "Cells with top TCR") +
  theme(legend.position = "none")
  my_ggsave(
    slug = glue("umap-has_tcr_top_clones"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.8,
    # width = 3.5,
    width = 7,
    height = 3.5, units = "in", dpi = 300
  )
  
  x <- a1$obs %>%
    group_by(donor, leiden, case) %>%
    summarize(pct_tcr_top = 100 * sum(tcr_top_clone) / n()) %>%
    ungroup() %>%
    mutate(
      leiden = naturalfactor(as.character(leiden)),
      case = factor(case, c("Control", "Case"))
    )
  p <- ggplot(x) +
  aes(x = pct_tcr_top, y = leiden, group = donor, fill = case) +
  scale_fill_manual(
    values = c("Case" = "#E69F00", "Control" = "#000000"),
    name = NULL
  ) +
  geom_quasirandom(
    shape = 21, size = 3, groupOnX = FALSE, dodge.width = 0.5,
    stroke = 0.2
  ) +
  theme(
    legend.position = "none"
  )
  my_ggsave(
    "tcr-diversity-dots",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 0.75, width = 4,
    # height = 6,
    height = 5,
    units = "in", dpi = 300
  )

# How much gene expression variability do we see for each TCR?
########################################################################

  d <- a1$obs %>%
    dplyr::filter(has_tcr) %>%
    # dplyr::filter(donor %in% tcr_donors) %>%
    dplyr::group_by(donor) %>%
    dplyr::count(
      TRBV,
      TRAV,
      TRB_cdr3_trim,
      TRA_cdr3_trim,
      case
    ) %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(
      rank = seq(length(n)),
      percent = 100 * n / sum(n)
    ) %>%
    dplyr::mutate(
      tcr_q = cut(rank, breaks = 5, labels = FALSE)
    ) %>%
    dplyr::mutate(tcr_clone = paste(TRBV, TRAV, TRB_cdr3_trim, TRA_cdr3_trim)) %>%
    ungroup()
  #
  top_tcr_clones <- d$tcr_clone[d$rank <= 5]
  # d[d$n == 1 & d$rank == 2,]$donor
  # a1$obs %>% filter(donor == "SIC_134") %>% count(has_tcr)
  top_tcr_clones <- (
    a1$obs %>% 
    filter(tcr_clone %in% top_tcr_clones) %>%
    count(tcr_clone) %>%
    filter(n > 1)
  )$tcr_clone

  a1$logcpm <- do_log1p_cpm(a1$counts)

  tcr_sds <- sapply(top_tcr_clones, function(tcr_clone) {
    my_cells <- a1$obs$cell[a1$obs$tcr_clone == tcr_clone]
    my_sd <- mean(proxyC::rowSds(a1$logcpm[a1$pca$genes,my_cells]))
  })

  d$gexp_sd <- tcr_sds[d$tcr_clone]

  p <- d %>% filter(!is.na(gexp_sd)) %>%
    ggplot() +
    aes(x = percent, y = gexp_sd, fill = case) +
    geom_quasirandom(shape = 21, stroke = 0.2) +
    scale_fill_manual(values = pals::okabe()) +
    scale_x_log10()
  my_ggsave(
    "tcr-gexp_sd",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 5,
    height = 3,
    units = "in", dpi = 300
  )


# Do expanded TCRs have similar CDR3 sequences?
########################################################################
  donor_dists <- sapply(unique(d$donor), function(my_donor) {
    my_seqs <- (
      d %>% filter(donor == my_donor, rank <= 10) %>%
        mutate(seq = paste(TRB_cdr3_trim))
    )$seq
    retval <- Biostrings::stringDist(
      x                  = AAStringSet(my_seqs),
      method             = "substitutionMatrix",
      type               = "global",
      substitutionMatrix = "BLOSUM80",
      gapOpening         = 0,
      gapExtension       = 0
    ) %>% abs
    mean(retval)
    my_seqs <- (
      d %>% filter(donor == my_donor, rank > 10 & rank <= 20) %>%
        mutate(seq = paste(TRB_cdr3_trim))
    )$seq
    retval2 <- Biostrings::stringDist(
      x                  = AAStringSet(my_seqs),
      method             = "substitutionMatrix",
      type               = "global",
      substitutionMatrix = "BLOSUM80",
      gapOpening         = 0,
      gapExtension       = 0
    ) %>% abs
    mean(retval) / mean(retval2)
  })
  p <- d %>% select(donor, case) %>% unique %>% mutate(dist_top = donor_dists[donor]) %>%
    ggplot() +
    aes(x = case, y = dist_top, fill = case) +
    geom_quasirandom(shape = 21, stroke = 0.2, size = 3) +
    scale_fill_manual(values = pals::okabe())
  my_ggsave(
    "tcr-top_dist",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 5,
    height = 3,
    units = "in", dpi = 300
  )


# Heatmaps counting TCRs across GEX clusters
#######################################################################

n_clusters <- length(unique(a1$obs$leiden))
mats <- list(
  matrix(0, nrow = n_clusters, ncol = n_clusters),
  matrix(0, nrow = n_clusters, ncol = n_clusters)
)
case_levels <- c("Case", "Control")
mat_d <- a1$obs %>%
  filter(
    has_tcr#, tcr_top_clone
  )
for (my_case in seq_along(case_levels)) {
 for (i in seq(n_clusters)) {
   i_tcr_ids <- unique((mat_d %>% filter(case == case_levels[my_case], leiden == i))$tcr_clone)
   for (j in seq(n_clusters)) {
     if (i == j) {
       next 
     }
     j_tcr_ids <- unique((mat_d %>% filter(case == case_levels[my_case], leiden == j))$tcr_clone)
     n_ij <- length(intersect(i_tcr_ids, j_tcr_ids))
     mats[[my_case]][i, j] <- n_ij / length(i_tcr_ids)
   }
 }
}

 h_colors <- scico::scico(n = 100, palette = "grayC")
 # h_max <- max(c(mats[[1]], mats[[2]]))
 h_max <- 1
 h1 <- Heatmap(
   matrix = mats[[2]],
   col = colorRamp2(seq(0, h_max, length.out = length(h_colors)), h_colors),
   cluster_rows = FALSE,
   cluster_columns = FALSE,
   row_labels = seq(nrow(mats[[1]])),
   column_labels = seq(nrow(mats[[1]])),
   show_heatmap_legend = FALSE,
   column_title = "Control"
 )
 h2 <- Heatmap(
   matrix = mats[[1]],
   col = colorRamp2(seq(0, h_max, length.out = length(h_colors)), h_colors),
   cluster_rows = FALSE,
   cluster_columns = FALSE,
   row_labels = seq(nrow(mats[[1]])),
   column_labels = seq(nrow(mats[[1]])),
   # name = "TCRs",
   name = "Fraction",
   column_title = "Case"
 )
 h_list <- h1 + h2
 draw(h_list, 
   column_title = "Number of distinct TCRs shared between clusters",
   column_title_gp = gpar(fontsize = 16))


}

