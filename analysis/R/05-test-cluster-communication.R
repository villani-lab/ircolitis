#!/usr/bin/env Rscript

library(conflicted)
library(Matrix)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
# library(cowplot)
library(data.table)
# library(doParallel)
# library(foreach)
library(ggbeeswarm)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(ggstance)
library(ggtext)
library(glue)
library(harmony)
library(janitor)
library(magrittr)
library(naturalsort)
library(OmnipathR)
library(pals)
library(patchwork)
library(pbapply)
# library(princurve)
library(qs)
library(readxl)
library(rhdf5)
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
# library(fgsea)
# My functions
source("R/functions/helpers.R")
source("R/functions/read-h5-files.R")
source("R/functions/write-mtx.R")
source("R/functions/do-analysis.R")
source("R/plot-analysis.R")
source("R/functions/write-de.R")
source("R/colors-Luoma2020.R")
source("R/functions/theme-kamil.R")
source("R/load-sample-therapy.R")
source("R/sample-ids.R")
theme_set(theme_kamil)
#
conflict_prefer("filter", "dplyr")
conflict_prefer("discard", "purrr")
conflict_prefer("reduce", "purrr")
#
# Define important gene sets
########################################################################
source("R/mt-tcr-bcr-genes.R")
#
get_cluster_groups <- function(analysis_name) {
  cluster_groups <- list()
  if (analysis_name == "n3_2") {
    cluster_groups[["e1 Immature epithelial cells"]] <- c("8", "3", "14")
    cluster_groups[["e2 Absorptive epithelial cells"]] <- c("1", "5", "6", "16", "18", "9", "12")
    cluster_groups[["e3 Mature, absorptive epithelial cells"]] <- c("2", "11", "20")
    cluster_groups[["e4 Secretory cells"]] <- c("4", "10", "17", "15", "19")
    cluster_groups[["e5 Mesenchymal cells"]] <- c("7", "23", "21", "22", "13")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_t4_cd8_1_2") {
    cluster_groups[["t1 ITGB2"]] <- c("3", "11")
    # cluster_groups[["g2 ITGAE ID3 CD8 Trm"]] <- c("1", "4", "5", "9")
    # cluster_groups[["g3 ITGAE ZNF683 CD8 Trm"]] <- c("7", "6", "2")
    cluster_groups[["t2 ITGAE"]] <- c("1", "4", "5", "9", "7", "6", "2")
    cluster_groups[["t3 Other"]] <- c("8", "10", "12")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_t4_cd4_2_2") {
    cluster_groups[["cd4 Treg"]] <- c("8", "9", "5")
    cluster_groups[["cd4 other"]] <- c("1", "2", "3", "4")
    cluster_groups[["cd4 CXCL13"]] <- c("6", "7", "10")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "blood2_tcell5_cd8_5") {
    cluster_groups[["bt1 MHC-II"]] <- c("2", "12", "8")
    cluster_groups[["bt2 Innate"]] <- c("4", "14", "6", "10")
    cluster_groups[["bt3 CX3CR1"]] <- c("7", "13", "1", "9")
    cluster_groups[["bt4 IL7R"]] <- c("3", "11")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_m3_2") {
    cluster_groups[["m1"]] <- c("1", "6", "5", "2", "3")
    cluster_groups[["m2"]] <- c("7", "4")
    cluster_groups[["m3"]] <- c("8")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_b5_1_3") {
    cluster_groups[["b"]] <- c("2", "3", "5", "11", "12")
    cluster_groups[["plasma"]] <- c("1", "4", "6", "7", "8", "9", "10", "13")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  unlist(cluster_groups)
}


# Load data
########################################################################
file_ccc <- "cache/05-test-cluster-communication.qs"
if (file.exists(file_ccc)) {
  qload(file_ccc)
  log2cpm <- do_log2cpm(counts, median(colSums(counts)))
} else { 
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
  #
  a_slugs <- c(
    "a20/a12_4_4_b5_1_3",
    "a20/a12_4_4_m3_2",
    "a20/a12_4_4_t4_cd4_2_2",
    "a20/a12_4_4_t4_cd8_1_2",
    "a20/n3_2"
  )
  for (a in a_slugs) {
    a_dir <- dirname(a)
    a_name <- basename(a)
    a_file <- glue("results/{a_dir}/{a_name}/data/{a_name}.qs")
    #
    message(a_dir)
    message(a_name)
    message(a_file)
    message(file.exists(a_file))
    #
    stopifnot(file.exists(a_file))
    #
    print_status(glue("Reading {a_file}"))
    a1 <- qread(a_file)
    print_status(glue("done"))
    #
    if (a_name == "a12_4_4_t4_cd8_1_2") {
      a1$obs$leiden <- recluster_cd8_leiden151(a1$obs$leiden1.51)
    }
    if (a_name == "a12_4_4_t4_cd4_2_2") {
      a1$obs$leiden <- a1$obs$leiden0.933
    }
    if (a_name == "a12_4_4_m3_2") {
      a1$obs$leiden <- a1$obs$leiden0.933
    }
    if (a_name == "a12_4_4_b5_1_3") {
      a1$obs$leiden <- a1$obs$leiden0.933
    }
    if (a_name == "n3_2") {
      a1$obs$leiden <- a1$obs$leiden0.933
    }
    #
    if (!"drug" %in% colnames(a1$obs)) {
      a1$obs <- left_join(a1$obs, sample_therapy, by = "donor")
    }
    assign(a_name, a1)
  }
  a_labels <- c(
    "B cell subsets",
    "Myeloid cell subsets",
    "CD4 T cell subsets",
    "CD8 T cell subsets",
    "Epithelial cell subsets"
  )
  #
  counts <- lapply(a_slugs, function(a_slug) {
    a1 <- get(basename(a_slug))
    return(a1$counts)
  })
  x <- reduce(lapply(counts, rownames), intersect)
  stopifnot(all(rownames(counts[[1]]) == x))
  #
  counts[[5]] <- counts[[5]][rownames(counts[[1]]),]
  #
  counts <- reduce(counts, cbind)
  #
  log2cpms <- lapply(a_slugs, function(a_slug) {
    a1 <- get(basename(a_slug))
    log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
    return(log2cpm)
  })
  #
  lapply(log2cpms, nrow)
  lapply(log2cpms, dim)
  #
  x <- reduce(lapply(log2cpms, rownames), intersect)
  stopifnot(all(rownames(log2cpms[[1]]) == x))
  #
  log2cpms[[5]] <- log2cpms[[5]][rownames(log2cpms[[1]]),]
  #
  log2cpm <- reduce(log2cpms, cbind)
  rm(list = c("log2cpms", "a1", basename(a_slugs)))
  #
  a_prefixes <- c("B", "M", "T", "CT", "E")
  metas <- lapply(seq_along(a_slugs), function(i) {
    a_slug <- a_slugs[i]
    a1 <- get(basename(a_slug))
    meta <- a1$obs
    meta$leiden <- sprintf("%s%s", a_prefixes[i], meta$leiden)
    return(meta)
  })
  lapply(metas, dim)
  #
  meta_cols <- reduce(lapply(metas, colnames), intersect)
  # meta_cols <- meta_cols[!str_detect(meta_cols, "^[PC|UMAP]")]
  meta_cols <- meta_cols[!str_detect(meta_cols, "^PC")]
  #
  obs <- do.call(rbind, lapply(metas, function(this_meta) {
    this_meta[,meta_cols]
  }))
  rm(metas)
  #
  stopifnot(all(obs$cell == colnames(log2cpm)))
  #
  table(obs$leiden)
  obs$dataset <- str_remove(obs$leiden, "\\d+$")
  #
  qsavem(counts, obs, ensembl_to_symbol, file = file_ccc)
}


# Covariance of subcluster abundances
########################################################################
{

  cc0_file <- "results/a20/covarying-abundance/cc_level0.qs"
  cc0_case_file <- "results/a20/covarying-abundance/cc_level0_case.qs"
  if (file.exists(cc0_file) && file.exists(cc0_case_file)) {
    cc0 <- qread(cc0_file)
    cc0_case <- qread(cc0_case_file)
  } else {
    # Get all pairs of clusters.
    all_clusters <- sort(unique(obs$leiden))
    cc0 <- as.data.frame(t(combn(all_clusters, 2)))
    colnames(cc0) <- c("c1", "c2")
    comp <- rbindlist(lapply(unique(obs$dataset), function(this_dataset) {
      comp <- obs %>%
        dplyr::filter(dataset == this_dataset) %>%
        dplyr::group_by(donor, case, leiden) %>%
        dplyr::summarize(n = n(), .groups = "drop_last") %>%
        dplyr::mutate(percent = n / sum(n))
    }))
    # comp %>% group_by(donor) %>% summarize(total = sum(percent))
    cc0$estimate <- 0
    cc0$p.value <- 1
    cc0$n <- 0
    for (i in seq(nrow(cc0))) {
      x <- inner_join(
        x = comp %>% filter(leiden == cc0$c1[i]),
        y = comp %>% filter(leiden == cc0$c2[i]),
        by = "donor"
      )
      stat <- cor.test(x$percent.x, x$percent.y, method = "spearman", exact = FALSE)
      cc0$n[i] <- nrow(x)
      cc0$estimate[i] <- stat$estimate
      cc0$p.value[i] <- stat$p.value
    }
    cc0$fdr <- p.adjust(cc0$p.value, method = "fdr")
    cc0 <- as_tibble(cc0) %>% arrange(p.value)
    cc0 %<>% mutate(epi = substr(c1, 1, 1) == "E" | substr(c2, 1, 1) == "E")
    #
    cc0_case <- rbind(
      cc0 %>% select(c1, c2) %>% mutate(case = "Case"),
      cc0 %>% select(c1, c2) %>% mutate(case = "Control")
    )
    cc0_case$estimate <- 0
    cc0_case$p.value <- 1
    cc0_case$n <- 0
    for (i in seq(nrow(cc0_case))) {
      x <- inner_join(
        x = comp %>% filter(case == cc0_case$case[i], leiden == cc0_case$c1[i]),
        y = comp %>% filter(case == cc0_case$case[i], leiden == cc0_case$c2[i]),
        by = "donor"
      )
      stat <- cor.test(x$percent.x, x$percent.y, method = "spearman", exact = FALSE)
      cc0_case$n[i] <- nrow(x)
      cc0_case$estimate[i] <- stat$estimate
      cc0_case$p.value[i] <- stat$p.value
    }
    cc0_case <- as_tibble(cc0_case) %>% arrange(p.value)
    cc0_case <- cc0_case %>% group_by(case) %>% mutate(fdr = p.adjust(p.value, method = "fdr"))
    cc0_case %<>% mutate(epi = substr(c1, 1, 1) == "E" | substr(c2, 1, 1) == "E")
    qsave(cc0, cc0_file)
    qsave(cc0_case, cc0_case_file)
  }

  # Write tables with results
  cc0_both <- rbind(
    cc0 %>% mutate(case = "Both"),
    cc0_case
  ) %>%
  rename(group = case, n_donors = n, cluster1 = c1, cluster2 = c2) %>%
  select(-epi) %>%
  relocate(group, cluster1, cluster2, n_donors, estimate, p.value) %>%
  arrange(p.value)
  fix_names <- function(x) {
    x <- str_replace(x, "^CT", "CD8-")
    x <- str_replace(x, "^T", "CD4-")
    x <- str_replace(x, "^B", "B-")
    x <- str_replace(x, "^M", "MP-")
    x <- str_replace(x, "^E", "E-")
  }
  cc0_both$cluster1 <- fix_names(cc0_both$cluster1)
  cc0_both$cluster2 <- fix_names(cc0_both$cluster2)
  fwrite(cc0_both, "paper/spearman-cluster-abundance.tsv", sep = "\t")

  cc0_case %>% filter(n >= 10, epi)
  cc0_case %>% filter(c1 == "B11", c2 == "E11")

  pdf_file <- "results/a20/covarying-abundance/level0/heatmap-signif.pdf"
  cols <- naturalsort(unique(with(cc0, c(c1, c2))))
  top_res_mat <- matrix(0, ncol = length(cols), nrow = length(cols))
  rownames(top_res_mat) <- cols
  colnames(top_res_mat) <- cols
  for (i in seq(nrow(cc0))) {
    row <- cc0[i,]
    c_sort <- rev(naturalsort(c(row$c1, row$c2)))
    top_res_mat[c_sort[2], c_sort[1]] <- -log10(row$p.value)
  }
  top_res_mat <- top_res_mat + t(top_res_mat)
  #
  # legend_labs <- unique(round(seq(min(log10(top_res_mat + 1)), max(log10(top_res_mat + 1)), length.out = 5)))
  legend_labs <- unique(round(seq(min(top_res_mat), max(top_res_mat), length.out = 5)))
  ht <- Heatmap(
    matrix = log10(top_res_mat + 1),
    col = rev(scico(pal = "oslo", n = 20)),
    row_names_side = "left",
    row_split = str_remove(rownames(top_res_mat), "\\d+"),
    column_split = str_remove(colnames(top_res_mat), "\\d+"),
    row_gap = unit(5, "mm"),
    column_gap = unit(5, "mm"),
    border = TRUE,
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "-log10 P",
    column_names_rot = 45,
    heatmap_legend_param = list(
      at = legend_labs,
      labels = legend_labs
    )
  )
  message(pdf_file)
  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  unlink(pdf_file)
  pdf(pdf_file, width = 16, height = 13)
  draw(ht)
  dev.off()


  pdf_file <- "results/a20/covarying-abundance/level0/heatmap-cc0.pdf"
  # my_clusters <- unique(c(cc0$c1[cc0$fdr < 0.01], cc0$c2[cc0$fdr < 0.01]))
  my_clusters <- unique(c(
    # cc0_case$c1[cc0_case$fdr < 0.01 & cc0_case$case == "Case"],
    # cc0_case$c2[cc0_case$fdr < 0.01 & cc0_case$case == "Case"]
    cc0_case$c1[cc0_case$fdr < 0.01],
    cc0_case$c2[cc0_case$fdr < 0.01]
  ))
  my_donors <- (
    comp %>% filter(leiden %in% my_clusters) %>%
    select(donor, leiden) %>%
    mutate(leiden = str_remove(leiden, "\\d+")) %>%
    unique %>%
    count(donor) %>%
    filter(n > 1)
  )$donor
  mat <- dcast(
    data = as.data.table(comp %>% filter(leiden %in% my_clusters, donor %in% my_donors)),
    formula = leiden ~ donor,
    value.var = "percent",
    fill = 0
  )
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat * 100)
  rownames(mat) <- mat_rows
  #
  column_dat <- left_join(
    x = data.frame(donor = colnames(mat)),
    y = comp %>% select(donor, case) %>% unique,
    by = "donor"
  )
  stopifnot(all(column_dat$donor == colnames(mat)))
  column_ha <- HeatmapAnnotation(
    case = column_dat$case,
    col = my_colors
  )
  #
  ht <- Heatmap(
    matrix = mat,
    col = rev(scico(pal = "oslo", n = 20)),
    top_annotation = column_ha,
    row_names_side = "left",
    row_split = str_remove(rownames(mat), "\\d+"),
    row_gap = unit(5, "mm"),
    # column_split = str_remove(colnames(top_res_mat), "\\d+"),
    # column_gap = unit(5, "mm"),
    border = TRUE,
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    name = "Percent",
    column_names_rot = 45
  )
  message(pdf_file)
  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.3 + 1, height = nrow(mat) * 0.3 + 1)
  draw(ht)
  dev.off()

  tweak_name <- function(x) {
    x <- str_replace_all(x, "CT", "CD8-")
    x <- str_replace_all(x, "T", "CD4-")
    x <- str_replace_all(x, "E", "E-")
    x <- str_replace_all(x, "M", "M-")
    x <- str_replace_all(x, "B", "B-")
    x
  }

  my_cc0 <- cc0 %>%
    filter(substr(c1, 1, 1) != substr(c2, 1, 1), epi)
  for (i in 1:15) {
    my_c1 <- my_cc0$c1[i]
    my_c2 <- my_cc0$c2[i]
    x <- inner_join(
      x = comp %>% filter(leiden == my_c1),
      y = comp %>% filter(leiden == my_c2),
      by = "donor"
    )
    stat <- cor.test(x$percent.x, x$percent.y, method = "spearman", exact = FALSE)
    p <- ggplot(x) +
      aes(x = 100 * percent.x, y = 100 * percent.y, fill = case.x) +
      geom_point(shape = 21, size = 3, stroke = 0.2) +
      scale_fill_manual(values = rev(pals::okabe(2)), name = NULL) +
      labs(
        x = tweak_name(my_c1), y = tweak_name(my_c2),
        subtitle = glue("rho = {signif(stat$estimate, 2)}, P = {signif(stat$p.value, 2)}, FDR = {signif(my_cc0$fdr[i], 2)}")
      )
    my_ggsave(
      glue("spearman-{my_c1}-{my_c2}-p{signif(stat$p.value, 2)}"),
      out_dir = "results/a20/covarying-abundance/level0/spearman",
      types = "pdf",
      plot = p,
      scale = 1, width = 4.5, height = 3, units = "in", dpi = 300
    )
  }

  my_cc0 <- cc0 %>%
    filter(substr(c1, 1, 1) != substr(c2, 1, 1), epi) %>%
    mutate(cpair = paste(c1, c2)) %>%
    mutate(cpair = str_replace_all(cpair, "CT", "CD8-")) %>%
    mutate(cpair = str_replace_all(cpair, "T", "CD4-")) %>%
    mutate(cpair = str_replace_all(cpair, "E", "E-")) %>%
    mutate(cpair = str_replace_all(cpair, "M", "M-")) %>%
    mutate(cpair = str_replace_all(cpair, "B", "B-")) %>%
    mutate(label = ifelse(fdr < 0.05, cpair, ""))
  p <- ggplot(my_cc0) +
    aes(x = estimate, y = -log10(p.value), label = label) +
    geom_vline(xintercept = 0, size = 0.3, alpha = 0.3) +
    geom_point(size = 0.5) +
    geom_point(data = my_cc0 %>% filter(fdr < 0.05), color = "red") +
    geom_text_repel(
      min.segment.length = 0.2,
      segment.size = 0.3, size = 3
    ) +
    annotate(
      geom = "text",
      x = -Inf, y = 0.1,
      hjust = -0.5, vjust = 0,
      label = sprintf("%s pairs", nrow(my_cc0))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
    labs(
      x = "Spearman correlation",
      y = bquote("-Log"[10]~"P"),
      title = "Correlations of cell cluster abundances"
    )
  my_ggsave(
    "volcano-cc0",
    out_dir = "results/a20/covarying-abundance/level0",
    types = "pdf",
    plot = p,
    scale = 1, width = 5, height = 3.5, units = "in", dpi = 300
  )

  my_cc0 <- left_join(
    x = cc0_case %>% filter(case == "Case"),
    y = cc0_case %>% filter(case == "Control"),
    by = c("c1", "c2")
  ) %>%
    mutate(lab = glue("{c1} {c2}")) %>%
    filter(n.x >= 5 | n.y >= 5, substr(c1, 1, 1) != substr(c2, 1, 1), epi.x)
  my_cc0 %<>%
    group_by(c1, c2, case.x) %>%
    # mutate(fdr_min = min(fdr.x, fdr.y, na.rm = TRUE)) %>%
    mutate(signif = ifelse(fdr.x < 0.2 | fdr.y < 0.2, TRUE, FALSE)) %>%
    mutate(lab = ifelse(signif, lab, "")) %>%
    ungroup()
    # filter(str_detect(c1, "CT") | str_detect(c2, "CT"))
    # filter(n.x >= 5 | n.y >= 5, fdr.x < 0.05 | fdr.y < 0.05)
  # my_bonf <- 0.05 / length(unique(my_cc0$lab))
  p <- ggplot(my_cc0) +
    aes(x = sign(estimate.x) * -log10(p.value.x), y = sign(estimate.y) * -log10(p.value.y)) +
    # geom_hline(yintercept = -log10(my_bonf), size = 0.3, alpha = 0.3) + 
    # geom_vline(xintercept = -log10(my_bonf), size = 0.3, alpha = 0.3) + 
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.3) + 
    geom_vline(xintercept = 0, size = 0.3, alpha = 0.3) + 
    geom_point(aes(color = fdr.x < 0.2 | fdr.y < 0.2), size = 0.1) +
    geom_text_repel(
      mapping = aes(label = lab),
      size = 2, segment.size = 0.3, min.segment.length = 0.3
    ) +
    # scale_x_continuous(labels = function(x) ifelse(x == 0, "1", glue("10<sup>-{x}</sup>"))) +
    # scale_y_continuous(labels = function(x) ifelse(x == 0, "1", glue("10<sup>-{x}</sup>"))) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      name = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(
      x = "Case", y = "Control",
      title = "Signed Spearman p-values"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_markdown(),
      axis.text.y = element_markdown()
    )
  my_ggsave(
    glue("case-control-pvalues"),
    out_dir = "results/a20/covarying-abundance/level0",
    types = "pdf",
    plot = p,
    scale = 1, width = 3, height = 3.5, units = "in", dpi = 300
  )
  p <- ggplot(my_cc0) +
    aes(x = estimate.x, y = estimate.y) +
    geom_vline(xintercept = 0, size = 0.3, alpha = 0.3) +
    geom_hline(yintercept = 0, size = 0.3, alpha = 0.3) +
    geom_point(mapping = aes(color = signif), size = 0.5) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("FDR < 0.20", "FDR >= 0.20"),
      name = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    # geom_point(
    #   data = my_cc0 %>% filter(c1 == "CT5", c2 == "E18"),
    #   size = 1, color = "red"
    # ) +
    geom_text_repel(
      data = my_cc0 %>% filter(c1 == "CT5", c2 == "E18"),
      mapping = aes(x = estimate.x, y = estimate.y, label = lab),
      size = 5
    ) +
    annotate(
      geom = "text",
      x = -Inf, y = -0.8,
      hjust = -0.2, vjust = 0,
      label = sprintf("%s pairs", nrow(my_cc0))
    ) +
    theme(legend.position = "bottom") +
    coord_equal() +
    labs(
      x = "Within Cases", y = "Within Controls"
    )
  my_ggsave(
    glue("case-control-spearman"),
    out_dir = "results/a20/covarying-abundance/level0",
    types = "pdf",
    plot = p,
    scale = 1, width = 3, height = 3.5, units = "in", dpi = 300
  )

  for (i in which(my_cc0$fdr.x < 0.2 | my_cc0$fdr.y < 0.2)) {
    my_c1 <- my_cc0$c1[i]
    my_c2 <- my_cc0$c2[i]
    x <- inner_join(
      x = comp %>% filter(leiden == my_c1),
      y = comp %>% filter(leiden == my_c2),
      by = "donor"
    )
    ix <- x$case.x == "Case"
    stat <- cor.test(x$percent.x[ix], x$percent.y[ix], method = "spearman", exact = FALSE)
    sub1 <- glue(
      "rho = {signif(stat$estimate, 2)}, P = {signif(stat$p.value, 2)}"
    )
    stat <- cor.test(x$percent.x[!ix], x$percent.y[!ix], method = "spearman", exact = FALSE)
    sub2 <- glue(
      "rho = {signif(stat$estimate, 2)}, P = {signif(stat$p.value, 2)}"
    )
    p <- ggplot(x) +
      aes(100 * percent.x, 100 * percent.y, fill = case.x) +
      geom_point(shape = 21, size = 3, stroke = 0.2) +
      scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#333333"), name = NULL) +
      labs(
        x = my_c1, y = my_c2,
        subtitle = glue("Case: {sub1}\nControl: {sub2}")
      ) +
      facet_wrap(~ case.x)#, scales = "free")
    my_ggsave(
      glue("spearman-{my_c1}-{my_c2}"),
      out_dir = "results/a20/covarying-abundance/level0/spearman-case",
      types = "pdf",
      plot = p,
      scale = 1, width = 7, height = 3.5, units = "in", dpi = 300
    )
  }


  p <- ggplot(my_cc0) +
    aes(x = estimate.x, y = -log10(p.value.x)) +
    geom_point(size = 0.2)
  my_ggsave(
    glue("volcano-case"),
    out_dir = "results/a20/covarying-abundance/level0",
    types = "pdf",
    plot = p,
    scale = 1, width = 5, height = 3.5, units = "in", dpi = 300
  )

  my_cc0 %>% filter(
    n.x > 10 & n.y > 10,
    fdr.x < 0.05 | fdr.y < 0.05,
    sign(estimate.x) != sign(estimate.y)
  )

}

# Omnipath
########################################################################

{

my_genes <- c("IL26", "CD274", "CXCL11", "CXCL10", "CXCL1", "CXCL9", "CXCL13",
  "CXCR3", "CXCR5", "IL17A", "TNF", "IL10")

omni_full <- import_ligrecextra_interactions(organism = 9606)

# setdiff(my_genes, omni_full$source_genesymbol)
# setdiff(my_genes, omni_full$target_genesymbol)

"CD274" %in% omni_full$source_genesymbol

# icn <- OmnipathR::import_Omnipath_intercell()

icn_file <- "cache/omnipath-icn.rds"
if (!file.exists(icn_file)) {
  icn <- OmnipathR::import_intercell_network()
  saveRDS(icn, icn_file)
} else {
  icn <- readRDS(icn_file)
}

# This level of strictness seems to capture mostly receptor-ligands.
# For example, the non-pair IL10 - ICAM1 is discarded by this filter.
omni <- icn %>%
  filter(
    sources != "Wang",
    consensus_score_intercell_source >= 3,
    curation_effort >= 1,
    n_resources >= 4
  ) %>%
  select(
    target_genesymbol,
    source_genesymbol,
    # is_stimulation,
    # consensus_score_intercell_source,
    # sources,
    # n_references,
    # references
  ) %>%
  separate_rows(target_genesymbol, sep = "_") %>%
  unique
nrow(omni)
# setdiff(my_genes, omni$source_genesymbol)
# setdiff(my_genes, omni$target_genesymbol)
setdiff(my_genes, c(omni$source_genesymbol, omni$target_genesymbol))
#
symbol_to_ensembl <- genes$ensembl_id
names(symbol_to_ensembl) <- genes$symbol
omni$source_ens <- symbol_to_ensembl[omni$source_genesymbol]
omni$target_ens <- symbol_to_ensembl[omni$target_genesymbol]
sum(is.na(omni$source_ens))
sum(is.na(omni$target_ens))
omni <- omni %>% filter(!is.na(source_ens) & !is.na(target_ens))
nrow(omni)

# Try to get a better subset of the data
########################################################################

unique(icn %>% filter(target_genesymbol == "PDCD1") %>% select(sources))$sources

omni %>% filter(target_genesymbol == "PDCD1")

omni %>% filter(target_genesymbol == "CXCR3")

omni %>% filter(source_genesymbol == "IL26")

omni %>% filter(source_genesymbol == "IL10")

icn %>% filter(source_genesymbol == "IL26") %>% select(target_genesymbol, source_genesymbol, sources)

omni %>% filter(target_genesymbol == "NECTIN2")

icn %>% filter(target_genesymbol == "NECTIN2") %>% select(target_genesymbol, source_genesymbol, n_resources, sources)

icn %>% filter(target_genesymbol == "HAVCR2") %>% select(target_genesymbol, source_genesymbol, n_resources, sources)

# Get PMIDs
(
  # icn %>% filter(target_genesymbol == "KCNA3", source_genesymbol == "IL16") %>%
  icn %>% filter(target_genesymbol == "CXCR3") %>%
  # filter(category_intercell_source == "ligand") %>%
  # filter(category_intercell_target == "receptor") %>%
  filter(consensus_score_intercell_source >= 3, curation_effort >= 1) %>%
  select(references)
)$references %>% unlist %>% str_extract_all(":\\d+") %>% unlist %>% str_remove(":") %>%
unique

icn %>% filter(target_genesymbol == "CXCR3") %>%
  filter(consensus_score_intercell_source >= 3, curation_effort >= 1) %>%
  mutate(pmid = references %>%
  str_extract_all(":\\d+") %>%
  unlist %>%
  str_remove(":") %>%
  unique %>%
  str_flatten(collapse = ",")) %>%
  select(target_genesymbol, source_genesymbol, pmid)

my_genes <- strsplit("ITGA4 ITGB7 MADCAM1 VCAM1 FN1 CHD1 ITGAE IL12RB1 IL12RB2 IL23R TNFRSF1A TNFRSF1B JAK1 JAK2 TYK2 S1PR1 S1PR5 LTBR SLC5A11 TNFRSF21 TRAF2", " ")[[1]]
fwrite(
  file = "data/integrins.csv",
  x = omni %>% filter(target_genesymbol %in% my_genes | source_genesymbol %in% my_genes)
)

icn %>% filter(source_genesymbol == "CXCL10") %>% #, target_genesymbol == "CXCR3") %>%
  select(
    source_genesymbol, target_genesymbol,
    category_intercell_source, category_intercell_target
  )

icn %>% filter(source_genesymbol == "CXCL10", target_genesymbol == "DPP4") %>% t

  # select(
  #   source_genesymbol, target_genesymbol,
  #   category_intercell_source, category_intercell_target
  # )

# Cell communication
########################################################################

# # TODO Consider using the muti R package https://mdscheuerell.github.io/muti/

# ram <- read_excel(
#   "data/Ramilowski2015/SupplementaryData2.xlsx", sheet = "All.Pairs"
# )
# ram <- janitor::clean_names(ram)
# ix <- ram$ligand_approved_symbol %in% ensembl_to_symbol[rownames(log2cpm)]
# sum(ix) / nrow(ram)
# ram <- ram[ix,]
# ix <- ram$receptor_approved_symbol %in% ensembl_to_symbol[rownames(log2cpm)]
# sum(ix) / nrow(ram)
# ram <- ram[ix,]
# nrow(ram)

# # Get all pairs of clusters.
# obs$cluster <- obs$leiden
# all_clusters <- sort(unique(obs$cluster))
# cc <- expand.grid(all_clusters, all_clusters)
# cc <- cc[cc$Var1 != cc$Var2,]

# # Get Ramilowski genes that are present in our data.
# ram_genes <- unique(rbind(
#   genes %>% filter(symbol %in% ram$ligand_approved_symbol),
#   genes %>% filter(symbol %in% ram$receptor_approved_symbol)
# ))
# ram_genes <- ram_genes[!duplicated(ram_genes$symbol),]
# ram_genes <- ram_genes %>%
#   filter(ensembl_id %in% rownames(log2cpm))
# stopifnot(length(unique(ram_genes$symbol)) == nrow(ram_genes))
# stopifnot(all(ram_genes$ensembl_id %in% rownames(log2cpm)))

# my_genes <- c("CXCL11", "CXCL1", "CXCL9", "CXCL13", "IL6", "CXCR3")
# my_genes %in% ram_genes$symbol

# # log2cpm_means <- Matrix::rowMeans(log2cpm)
# # plot(hist(log2cpm_means[log2cpm_means > 0.1], breaks = 100))
# # sum(log2cpm_means > 0.1)

omni[1:5,]

omni_symbols <- unique(c(omni$source_genesymbol, omni$target_genesymbol))
length(omni_symbols)
setdiff(my_genes, omni_symbols)
# setdiff(my_genes, ram_genes$symbol)

omni_ens <- symbol_to_ensembl[omni_symbols] %>%
  discard(is.na) %>% unname %>% unique
omni_ens <- intersect(omni_ens, rownames(log2cpm))
setdiff(my_genes, ensembl_to_symbol[omni_ens])
stopifnot(all(omni_ens %in% rownames(log2cpm)))
omni_sym <- ensembl_to_symbol[omni_ens]
stopifnot(length(omni_sym) == length(unique(omni_sym)))

omni <- omni %>% filter(source_ens %in% rownames(counts), target_ens %in% rownames(counts))

# length(setdiff(omni_symbols, genes$symbol))
# head(setdiff(omni_symbols, genes$symbol))

}

# Product analysis
########################################################################

obs$cluster <- obs$dataset

gg1_file <- "results/a20/covarying-abundance/gg_level1_case.qs"
if (file.exists(gg1_file)) {
  message("Reading ", gg1_file)
  gg1 <- qread(gg1_file)
} else {
  source("R/covarying-composition.R")
  gg1 <- gene_gene_products(
    counts            = counts,
    cluster           = obs$dataset,
    donor             = obs$donor,
    case              = obs$case,
    source_ens        = omni$source_ens,
    target_ens        = omni$target_ens,
    min_donors        = 10,
    min_cells         = 25,
    ensembl_to_symbol = ensembl_to_symbol
  )
  gg1$limma %<>% mutate(epi = c1 == "E" | c2 == "E", inter = substr(c1, 1, 1) != substr(c2, 1, 1))
  qsave(gg1, gg1_file)
}
# gg1$limma %<>% group_by(c1, c2) %>% mutate(fdr = p.adjust(P.Value, method = "fdr"))

gg1$products[1:5,1:5]
# gg1_sum$products[1:5,1:5]

gg1$limma %>% select(
  c1, g1, c2, g2, logFC, P.Value, p1, p2, fdr
) %>% arrange(fdr)

# gg1_sum$limma %>% select(
#   c1, g1, c2, g2, logFC, P.Value, p1, p2, fdr
# ) %>% arrange(fdr)

gg1$limma %>%
  ungroup %>%
  filter(fdr < 0.05, p1 > 0.005, p2 > 0.005) %>%
  select(g1, g2) %>%
  unique %>%
  nrow

# all.equal(gg1$limma$P.Value, gg1_sum$limma$P.Value)

gg2_file <- "results/a20/covarying-abundance/gg_level2_case.qs"
if (file.exists(gg2_file)) {
  message("Reading ", gg2_file)
  gg2 <- qread(gg2_file)
} else {
  gg2 <- gene_gene_products(
    counts            = counts,
    cluster           = obs$leiden,
    donor             = obs$donor,
    case              = obs$case,
    source_ens        = omni$source_ens,
    target_ens        = omni$target_ens,
    min_donors        = 10,
    min_cells         = 25,
    ensembl_to_symbol = ensembl_to_symbol
  )
  gg2$limma %<>% mutate(epi = c1 == "E" | c2 == "E", inter = substr(c1, 1, 1) != substr(c2, 1, 1))
  gg2$limma %>% select(
    c1, g1, c2, g2, logFC, P.Value, p1, p2
  )
  qsave(gg2, gg2_file)
}

fwrite(
  x = gg1$limma %>% mutate_if(is.numeric, signif),
  file ="results/a20/covarying-abundance/ccc-sums-level1.tsv",
  sep = "\t"
)
fwrite(
  x = gg2$limma %>% mutate_if(is.numeric, signif),
  file ="results/a20/covarying-abundance/ccc-sums-level2.tsv",
  sep = "\t"
)

{
  wb <- openxlsx::createWorkbook()
  fname <- "results/a20/covarying-abundance/gg_level1_case.xlsx"
  unlink(fname)
  #
  this_sheet <- "Sheet 1"
  openxlsx::addWorksheet(wb, this_sheet)
  openxlsx::writeDataTable(
    wb,
    this_sheet,
    x = gg1$limma %>% mutate_if(is.numeric, signif),
    rowNames = FALSE, tableStyle = "TableStyleLight1"
  )
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}

{
  wb <- openxlsx::createWorkbook()
  fname <- "results/a20/covarying-abundance/gg_level2_case.xlsx"
  unlink(fname)
  #
  this_sheet <- "Sheet 1"
  openxlsx::addWorksheet(wb, this_sheet)
  openxlsx::writeDataTable(
    wb,
    this_sheet,
    x = gg2$limma %>% mutate_if(is.numeric, signif),
    rowNames = FALSE, tableStyle = "TableStyleLight1"
  )
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}


# Figure 7
########################################################################

# For epithelial cells, run another CCC with these clusters
cluster_groups <- list()
cluster_groups[["Immature epithelial cells"]] <- c("E8", "E3", "E14")
cluster_groups[["Absorptive epithelial cells"]] <- c("E1", "E5", "E6", "E16", "E18", "E9", "E12")
cluster_groups[["Mature absorptive epithelial cells"]] <- c("E2", "E11", "E20")
cluster_groups[["Secretory cells"]] <- c("E4", "E10", "E17", "E15", "E19")
cluster_groups[["Mesenchymal cells"]] <- c("E7", "E23", "E21", "E22")
cluster_groups[["Endothelial cells"]] <- c("E13")
cluster_groups <- unlist(split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups)))

ix <- obs$leiden %in% names(cluster_groups)
obs$cluster2 <- obs$leiden
obs$cluster2[ix] <- cluster_groups[obs$cluster2[ix]]
ix <- str_detect(obs$cluster2, "^(M|B|T)\\d")
obs$cluster2[ix] <- str_remove(obs$cluster2[ix], "\\d+$")
obs %>% count(cluster2) %>% arrange(-n) %>% head(10)

my_omni <- omni %>% filter(source_ens %in% rownames(counts), target_ens %in% rownames(counts))

if (!nrow(my_omni %>% filter(source_genesymbol == "CEACAM1", target_genesymbol == "HAVCR2"))) {
  my_omni <- rbind(omni,
    tibble(target_genesymbol = "HAVCR2", source_genesymbol = "CEACAM1",
    source_ens = "ENSG00000079385", target_ens = "ENSG00000135077")
  )
}
my_omni %>% filter(source_genesymbol == "CEACAM1")

# Make ICAM genes the target
ix <- which(with(my_omni, source_genesymbol %in% c("ICAM1", "ICAM2", "ICAM3", "CEACAM1")))
x <- my_omni[ix,]
my_omni[ix,] <- x[,c(2,1,4,3)]
ix <- which(with(my_omni, target_genesymbol %in% c("CXCR3", "PDCD1")))
x <- my_omni[ix,]
my_omni[ix,] <- x[,c(2,1,4,3)]

my_omni %>% filter(target_genesymbol == "PDCD1") # should be empty
my_omni %>% filter(source_genesymbol == "PDCD1")
my_omni %>% filter(source_genesymbol == "HAVCR2")

source("R/covarying-composition.R")
gg3_file <- "results/a20/covarying-abundance/gg_level3_case.qs"
if (file.exists(gg3_file)) {
  gg3 <- qread(gg3_file)
} else {
  gg3 <- gene_gene_products(
    counts            = counts,
    cluster           = obs$cluster2,
    donor             = obs$donor,
    case              = obs$case,
    source_ens        = my_omni$source_ens,
    target_ens        = my_omni$target_ens,
    min_donors        = 10,
    min_cells         = 25,
    ensembl_to_symbol = ensembl_to_symbol
  )
  gg3$limma %>% select(
    c1, g1, c2, g2, logFC, P.Value, p1, p2
  )
  qsave(gg3, gg3_file)
}

top_gg1_limma <- gg1$limma %>%
  arrange(P.Value) %>%
  # filter(inter, epi) %>%
  # filter(inter) %>%
  filter(p1 > 0.01, p2 > 0.01) %>%
  select(-t, -B, -CI.L, -CI.R, -donors1, -donors2, -id, -adj.P.Val, -epi, -inter) %>%
  # filter(g1 == "PDCD1") %>%
  head(100)
top_gg1_limma

gg1$limma %>%
  arrange(P.Value) %>%
  filter(inter) %>%
  filter(p1 > 0.01, p2 > 0.01) %>%
  select(-t, -B, -CI.L, -CI.R, -donors1, -donors2, -id, -adj.P.Val, -epi, -inter) %>%
  mutate(fdr = p.adjust(P.Value, method = "BH")) %>%
  group_by(c2, c1) %>%
  summarize(n = sum(fdr < 0.05)) %>%
  filter(c1 == "E" | c2 == "E")

# top_gg1_limma <- gg1$limma %>% filter(g1 == "PDCD1") %>% filter(fdr < 0.01)

# top_gg1_limma <- d1 %>% filter(g1 == "CXCL10")

to_name <- c("E", "B", "M", "CD8 T", "CD4 T")
names(to_name) <- c("E", "B", "M", "CT", "T")
for (i in 1:nrow(top_gg1_limma)) {
  # if (i > 5) {break}
  c1 <- top_gg1_limma$c1[i]
  c2 <- top_gg1_limma$c2[i]
  g1 <- top_gg1_limma$g1[i]
  g2 <- top_gg1_limma$g2[i]
  #
  e1 <- names(which(ensembl_to_symbol == g1))
  e2 <- names(which(ensembl_to_symbol == g2))
  gg1$pb_meta$g1 <- gg1$pb[e1,]
  gg1$pb_meta$g2 <- gg1$pb[e2,]
  x <- left_join(
    x = gg1$pb_meta %>% filter(cluster == c2) %>% select(-g1),
    y = gg1$pb_meta %>% filter(cluster == c1) %>% select(-g2),
    by = c("donor", "case")
  )
  x$score <- x$g1 + x$g2
  x$is_case <- x$case == "Case"
  # t.test(formula = score ~ is_case, data = x)
  subtitle <- glue("FC = {signif(2^abs(top_gg1_limma$logFC[i]), 2)}, P = {signif(top_gg1_limma$P.Value, 2)}")
  p <- ggplot() +
    geom_boxplot(
      data = x,
      mapping = aes(y = case, x = g1 + g2, fill = case),
      alpha = 0.3, size = 0.3,
      outlier.shape = NA, coef = 0
    ) +
    geom_quasirandom(
      data = x,
      mapping = aes(y = case, x = g1 + g2, fill = case),
      shape = 21, size = 3, stroke = 0.3, bandwidth = 0.5, nbins = 2,
      dodge.width = 1, width = 0.5,
      groupOnX = FALSE
    ) +
    # scale_x_log10() +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    # annotation_logticks(sides = "b", size = 0.3) +
    scale_fill_manual( values = okabe(8), guide = "none") +
    labs(
      x = NULL, y = NULL, title = glue("{to_name[c2]} *{g2}* and {to_name[c1]} *{g1}*"),
      subtitle = subtitle
    ) +
    theme(
      plot.title = element_markdown()
    )
  my_ggsave(
    glue("{c2}_{c1}_{g2}_{g1}_score"),
    out_dir = "results/a20/covarying-abundance/level1/product-dots",
    type = "pdf",
    plot = p,
    scale = 1,
    width = 5,
    height = 1.5,
    units = "in", dpi = 300
  )
  p <- ggplot() +
    geom_point(
      data = x,
      mapping = aes(x = g1, y = g2, fill = append_n(case)),
      shape = 21, size = 3, stroke = 0.3
    ) +
    scale_fill_manual(name = NULL, values = rev(okabe(2))) +#, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    scale_y_continuous(breaks = pretty_breaks(4)) +
    guides(
      fill = guide_legend(override.aes = list(size = 5))
    ) +
    labs(
      x = glue("{to_name[c1]} *{g1}*"),
      y = glue("{to_name[c2]} *{g2}*"),
      title = NULL
      # title = glue("{to_name[c2]} *{g2}* and {to_name[c1]} *{g1}*")
    ) +
    theme(
      plot.title = element_markdown(),
      axis.title.x = element_markdown(),
      axis.title.y = element_markdown()
    )
  my_ggsave(
    glue("{c2}_{c1}_{g2}_{g1}_scatter"),
    out_dir = "results/a20/covarying-abundance/level1/product-dots",
    type = "pdf",
    plot = p,
    scale = 1,
    width = 5.5,
    height = 3,
    units = "in", dpi = 300
  )
}

d1 <- gg1$limma %>%
  filter(p1 > 0.01, p2 > 0.01) %>%
	group_by(c1, c2) %>%
  mutate(
    fdr = p.adjust(P.Value, method = "BH")
	) %>%
  ungroup() %>%
  mutate(
    gene_pair = glue("{g2} - {g1}"),
    cell_pair = glue("{c2} - {c1}")
  )
my_gene_pairs <- unique((d1 %>% filter(fdr < 0.001, epi))$gene_pair)
length(my_gene_pairs)
d1 <- d1 %>% filter(gene_pair %in% my_gene_pairs)
# gene_pairs <- (
#   d1 %>% group_by(gene_pair) %>% summarize(min_pval = min(P.Value)) %>% arrange(-min_pval)
# )$gene_pair
# d1$gene_pair <- factor(as.character(d1$gene_pair), levels = gene_pairs)
set.seed(1)
sorted_gene_pairs <- umap_sorted_rows(d1, "gene_pair ~ cell_pair", "logFC")
d1$gene_pair <- factor(as.character(d1$gene_pair), levels = rev(sorted_gene_pairs))
#
col_fun <- circlize::colorRamp2(
  seq(-max(abs(d1$logFC)), max(abs(d1$logFC)), length.out = 11),
  rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
)
# color_values <- scales::rescale(seq(min(d1$signed_p), max(d1$signed_p), length.out = 11))
color_values <- seq(min(d1$logFC), max(d1$logFC), length.out = 11)
p <- ggplot() +
  geom_tile(
    data = d1,
    mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
  ) +
  geom_point(
    data = d1 %>% filter(fdr < 0.05),
    mapping = aes(x = cell_pair, y = gene_pair),
    shape = 19, color = "white"
  ) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_x_discrete(
    position = "t", name = NULL, expand = c(0, 0)
  ) +
  facet_wrap(~ c2, ncol = 5, scales = "free_x") +
  scale_fill_gradientn(
    colors = col_fun(color_values),
    # name = bquote("log"[2]~ "Fold-Change"),
    name = "Fold-Change",
    guide = guide_colorbar(barwidth = 15),
    breaks = pretty_breaks(4),
    labels = function(x) fractional::fractional(2^x)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0),
    axis.text.y = element_text(face = "italic"),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )
my_ggsave(
  "heatmap-top-product",
  out_dir = file.path("results/a20/covarying-abundance/level1"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = length(unique(d1$cell_pair)) * 0.4 + 1,
  height = length(unique(d1$gene_pair)) * 0.3 + 1.1,
  units = "in",
  dpi = 300
)


# Product analysis with curated pairs, major cell types (E, CT, T, M, B)
########################################################################

obs$cluster <- obs$dataset

{

  my_pairs <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  my_pairs$target_ens <- unname(symbol_to_ensembl[my_pairs$g_target])
  my_pairs$source_ens <- unname(symbol_to_ensembl[my_pairs$g_source])
  # my_pairs[which(is.na(my_pairs$source_ens)),]
  # my_pairs[which(is.na(my_pairs$target_ens)),]
  my_pairs <- my_pairs %>% filter(target_ens %in% rownames(counts), source_ens %in% rownames(counts))
  my_pairs <- my_pairs[!duplicated(with(my_pairs, paste(target_ens, source_ens))),]
  # my_pairs <- my_pairs %>% filter(
  #   !cat %in% c("Vedolizumab", "Vedolizumab, Ertolizumab", "Ustekinumab", "Anti-TNF")
  # )
  #
  source("R/covarying-composition.R")
  gg1_cur <- gene_gene_products(
    counts            = counts,
    cluster           = obs$dataset,
    donor             = obs$donor,
    case              = obs$case,
    source_ens        = my_pairs$target_ens,
    target_ens        = my_pairs$source_ens,
    min_donors        = 10,
    min_cells         = 25,
    ensembl_to_symbol = ensembl_to_symbol
  )
  gg1_cur$limma %<>% mutate(epi = c1 == "E" | c2 == "E", inter = substr(c1, 1, 1) != substr(c2, 1, 1))
  # gg1_cur$limma %<>% group_by(c1, c2) %>% mutate(fdr = p.adjust(P.Value, method = "fdr"))
  gg1_cur$limma %>% select(
    c1, g1, c2, g2, logFC, P.Value, p1, p2
  )
  stopifnot(sum(duplicated(gg1_cur$limmma$id)) == 0)

  # names(which(ensembl_to_symbol == "FN1")) # "ENSG00000115414"
  # names(which(ensembl_to_symbol == "ITGB7")) # "ENSG00000139626"

  e1 <- "ENSG00000115414"
  e2 <- "ENSG00000139626"

  sum(log2cpm[names(which(ensembl_to_symbol == "FN1")), obs$dataset == "E"] > 0)
  sum(log2cpm[names(which(ensembl_to_symbol == "ITGB7")), obs$dataset == "CT"] > 0)

  gg1_cur$limma %>%
    select(c1, g1, c2, g2, logFC, P.Value, p1, p2) %>%
    filter(g1 == "FN1", g2 == "ITGB7")

  #
  d1 <- gg1_cur$limma %>%
    filter(p1 > 0.005, p2 > 0.005) %>%
    # filter(p1 > 0.007, p2 > 0.007) %>%
    group_by(c1, c2) %>%
    mutate(
      fdr = p.adjust(P.Value, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
      gene_pair = glue("{g2} - {g1}"),
      cell_pair = glue("{c2} - {c1}")
    ) %>%
    left_join(
      my_pairs %>% select(g_target, g_source, cat), by = c("g2" = "g_target", "g1" = "g_source")
    )
  d1 %>% count(cat)

  d1 %>%
    select(c1, g1, c2, g2, logFC, P.Value, p1, p2) %>%
    filter(g1 == "FN1", g2 == "ITGB7")

  #
  stopifnot(sum(duplicated(d1$id)) == 0)
  # my_gene_pairs <- unique((d1 %>% filter(fdr < 0.1, epi))$gene_pair)
  # my_gene_pairs <- unique((d1 %>% filter(fdr < 0.05))$gene_pair)
  exclude_gene_pairs <- c(
    "CD48 - CD244",
    "COL3A1 - LAIR1",
    "TRAF2 - TNFRSF4",
    "TRAF2 - TNFRSF9",
    "CD58 - CD2",
    "TNFRSF4 - TNFRSF4",
    "SLAMF1 - SLAMF1",
    "TNFSF8 - TNFRSF8",
    "TNFRSF14 - TNFSF14",
    "LTBR - TNFSF14",
    "ITGAM - CD40LG",
    "TNFSF15 - TNFRSF25"
  )
  # length(my_gene_pairs)
  d1 <- d1 %>%
    # filter(gene_pair %in% my_gene_pairs) %>%
    filter(!gene_pair %in% exclude_gene_pairs)
  nrow(d1)
  #
  set.seed(1)
  sorted_gene_pairs <- naturalsort(unique(as.character(d1$gene_pair)))
  d1$gene_pair <- factor(as.character(d1$gene_pair), levels = rev(sorted_gene_pairs))
  d1$c2 <- factor(as.character(d1$c2), c("E", "CT", "T", "M", "B"))
  # 
  #
  d1$cell_pair <- str_replace_all(d1$cell_pair, "\\bCT\\b", "CD8")
  d1$cell_pair <- str_replace_all(d1$cell_pair, "\\bT\\b", "CD4")
  d1$cell_pair <- str_replace_all(d1$cell_pair, "\\bM\\b", "MP")
  d1$cell_pair <- str_replace_all(d1$cell_pair, "\\bE\\b", "E/M")
  #
  d1$cat <- factor(d1$cat,
    c(
      "Co-Inhibitory",
      "Co-Stimulatory",
      "Metabolism",
      "Cytokines",
      "Tertiary Lymphoid",
      "Homing",
      "Ustekinumab",
      "Vedolizumab",
      "Anti-TNF"
    )
  )
  #
  my_levels <- levels(d1$gene_pair)
  my_order <- naturalorder(my_levels)
  ix <- d1$cat %in% c("Co-Inhibitory", "Homing")
  d1_pairs <- unique(as.character(d1$gene_pair[ix]))
  io <- rev(naturalorder(str_split_fixed(d1_pairs, " - ", 2)[,2]))
  d1_pairs <- d1_pairs[io]
  my_levels <- c(d1_pairs, setdiff(my_levels, d1_pairs))
  d1$gene_pair <- factor(as.character(d1$gene_pair), my_levels)


  d2 <- d1 %>% filter(!cat %in% c("Ustekinumab", "Vedolizumab", "Anti-TNF"))
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(d2$logFC)), max(abs(d2$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(d2$signed_p), max(d2$signed_p), length.out = 11))
  color_values <- seq(min(d2$logFC), max(d2$logFC), length.out = 11)
  message(length(unique(d2$gene_pair)), " gene pairs")
  p <- ggplot() +
    geom_tile(
      data = d2,
      mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
    ) +
    geom_point(
      data = d2 %>% filter(fdr < 0.05),
      mapping = aes(x = cell_pair, y = gene_pair),
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    # facet_wrap(cat ~ c2, ncol = 5, scales = "free_x") +
    facet_grid(cat ~ c2, scales = "free", space = "free") + #, switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      # name = bquote("log"[2]~ "Fold-Change"),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 15),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      # strip.text.y.left = element_text(angle = 0, hjust = 1),
      strip.text.y = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  # longest <- unique(d2$cat[nchar(d2$cat) == max(nchar(d2$cat))])
  get_longest <- function(x) {
    x <- as.character(x)
    head(unique(x[nchar(x) == max(nchar(x))]), 1)
  }
  longest <- get_longest(d2$cat)
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-selected1",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(d2$cell_pair)) * 0.4 + 1 + text_width,
    height = length(unique(d2$gene_pair)) * 0.3 + 1.1,
    units = "in",
    dpi = 300
  )


  # Figure 7
  d3 <- d1 %>% filter(cat %in% c("Ustekinumab", "Vedolizumab", "Anti-TNF"))
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(d3$logFC)), max(abs(d3$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(d3$signed_p), max(d3$signed_p), length.out = 11))
  color_values <- seq(min(d3$logFC), max(d3$logFC), length.out = 11)
  d1 %>%
    select(c1, g1, p1, c2, g2, p2, logFC, P.Value, fdr) %>%
    # filter(cat == "Vedolizumab") %>%
    filter(g1 %in% c("FN1", "VCAM1")) %>% 
    arrange(c1, c2, g1, g2)

  p <- ggplot() +
    geom_tile(
      data = d3,
      mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
    ) +
    geom_point(
      data = d3 %>% filter(fdr < 0.05),
      mapping = aes(x = cell_pair, y = gene_pair),
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    # facet_wrap(cat ~ c2, ncol = 5, scales = "free_x") +
    facet_grid(cat ~ c2, scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      # name = bquote("log"[2]~ "Fold-Change"),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 15),
      breaks = c(0, 1, 2, 4, 8),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  # longest <- unique(d3$cat[nchar(d3$cat) == max(nchar(d3$cat))])
  get_longest <- function(x) {
    x <- as.character(x)
    head(unique(x[nchar(x) == max(nchar(x))]), 1)
  }
  longest <- get_longest(d3$cat)
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-selected2",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(d3$cell_pair)) * 0.4 + 1 + text_width,
    height = length(unique(d3$gene_pair)) * 0.4 + 1.1,
    units = "in",
    dpi = 300
  )

}

{
  d1_count <- d1 %>% filter(fdr < 0.05, p1 > 0.01, p2 > 0.01) %>% count(c1, c2)
  # set.seed(3)
  # d1_count$c1 <- factor(d1_count$c1, umap_sorted_rows(d1_count, "c1 ~ c2", "n", n_neighbors = 2))
  # d1_count$c2 <- factor(d1_count$c2, umap_sorted_rows(d1_count, "c2 ~ c1", "n", n_neighbors = 2))
  d1_count$c1 <- factor(d1_count$c1, c("M", "CT", "T", "E", "B"))
  d1_count$c2 <- factor(d1_count$c2, rev(c("M", "CT", "T", "E", "B")))
  p <- ggplot() +
    geom_tile(
      data = d1_count,
      mapping = aes(x = c1, y = c2, fill = n)
    ) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    scale_fill_gradientn(
      colors = scico::scico(palette = "grayC", n = 11),
      name = "N",
      guide = guide_colorbar(barwidth = 8),
      breaks = pretty_breaks(4)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      legend.position = "bottom"
    )
  my_ggsave(
    "heatmap-top-product-selected-count",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(d1_count$c1)) * 0.4 + 1,
    height = length(unique(d1_count$c2)) * 0.4 + 1 + 1,
    units = "in",
    dpi = 300
  )
}


# Differentially expressed genes in each cell type (E, CD8, CD4, M, B)
########################################################################

my_pair_genes <- unique(c(d1$g1, d1$g2))

analyses <- c(
  "a12_4_4_t4_cd8_1_2",
  "a12_4_4_t4_cd4_2_2",
  "a12_4_4_m3_2",
  "a12_4_4_b5_1_3",
  "n3_2"
)
de <- rbindlist(lapply(analyses, function(analysis_name) {
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
  de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
  de <- fread(de_donor_file)
  de$celltype <- analysis_name
  return(de)
}))
de$celltype <- str_remove(de$celltype, "a12_4_4_")
to_celltype <- c("B", "M", "E", "CD4", "CD8")
names(to_celltype) <- c("b5_1_3", "m3_2", "n3_2", "t4_cd4_2_2", "t4_cd8_1_2")
de$celltype <- to_celltype[de$celltype]

{

  de_select <- de %>% filter(Gene %in% my_pair_genes)
  #
  x1 <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  x1 <- pivot_longer(x1, cols = c("g_target", "g_source"))
  gene_to_cat <- with(
    x1 %>% select(cat, value) %>% unique(),
    split(cat, value)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  rm(x1)
  #
  de_select$cat <- gene_to_cat[de_select$Gene]

  set.seed(1)
  sorted_genes <- umap_sorted_rows(de_select, "Gene ~ celltype", "logFC")
  de_select$Gene <- factor(as.character(de_select$Gene), levels = rev(sorted_genes))

  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(de_select$logFC)), max(abs(de_select$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(de_select$signed_p), max(de_select$signed_p), length.out = 11))
  color_values <- seq(min(de_select$logFC), max(de_select$logFC), length.out = 11)
  p <- ggplot() +
    geom_tile(
      data = de_select,
      mapping = aes(x = celltype, y = Gene, fill = logFC)
    ) +
    geom_point(
      data = de_select %>% filter(adj.P.Val < 0.05),
      mapping = aes(x = celltype, y = Gene),
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    facet_grid(rows = vars(cat), scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 10),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  longest <- unique(de_select$cat[nchar(de_select$cat) == max(nchar(de_select$cat))])
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-top-product-selected-foldchange",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(de_select$celltype)) * 0.5 + 1 + text_width,
    height = length(unique(de_select$Gene)) * 0.3 + 1.1,
    units = "in",
    dpi = 300
  )

}

# Differentially expressed genes in each cell cluster
########################################################################

de <- rbindlist(lapply(analyses, function(analysis_name) {
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
  de_donor_file <- glue("{out_dir}/de_case-vs-control.tsv.gz")
  de <- fread(de_donor_file)
  de$celltype <- analysis_name
  return(de)
}))
de$celltype <- str_remove(de$celltype, "a12_4_4_")
to_celltype <- c("B", "M", "E", "CD4", "CD8")
names(to_celltype) <- c("b5_1_3", "m3_2", "n3_2", "t4_cd4_2_2", "t4_cd8_1_2")
de$celltype <- to_celltype[de$celltype]
#
de$cluster <- sprintf("%s-%s", de$celltype, de$cluster)

{

  de_select <- de %>% filter(Gene %in% my_pair_genes)
  #
  x1 <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  x1 <- pivot_longer(x1, cols = c("g_target", "g_source"))
  gene_to_cat <- with(
    x1 %>% select(cat, value) %>% unique(),
    split(cat, value)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  rm(x1)
  #
  de_select$cat <- gene_to_cat[de_select$Gene]

  set.seed(1)
  sorted_genes <- umap_sorted_rows(de_select, "Gene ~ cluster", "logFC")
  de_select$Gene <- factor(as.character(de_select$Gene), levels = rev(sorted_genes))
  sorted_clusters <- umap_sorted_rows(de_select, "cluster ~ Gene", "logFC")
  de_select$cluster <- factor(as.character(de_select$cluster), levels = rev(sorted_clusters))

  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(de_select$logFC)), max(abs(de_select$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(de_select$signed_p), max(de_select$signed_p), length.out = 11))
  color_values <- seq(min(de_select$logFC), max(de_select$logFC), length.out = 11)
  p <- ggplot() +
    geom_tile(
      data = de_select,
      mapping = aes(x = cluster, y = Gene, fill = logFC)
    ) +
    geom_point(
      data = de_select %>% filter(adj.P.Val < 0.05),
      mapping = aes(x = cluster, y = Gene),
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    facet_grid(rows = vars(cat), scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 10),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  longest <- unique(de_select$cat[nchar(de_select$cat) == max(nchar(de_select$cat))])
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-top-product-selected-foldchange-cluster",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(de_select$cluster)) * 0.3 + 1 + text_width,
    height = length(unique(de_select$Gene)) * 0.3 + 1.1,
    units = "in",
    dpi = 300
  )

}

{

  x1 <- readxl::read_excel("data/drug-targets.xlsx")
  gene_to_cat <- with(
    x1 %>% select(Category, Target) %>% unique(),
    split(Category, Target)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  rm(x1)
  de_select <- de %>% filter(Gene %in% names(gene_to_cat))
  #
  de_select$cat <- gene_to_cat[de_select$Gene]
  # filter_genes <- unique((de_select %>% filter(adj.P.Val < 0.05))$Gene)
  # de_select <- de_select %>% filter(Gene %in% filter_genes)
  set.seed(1)
  sorted_genes <- umap_sorted_rows(de_select, "Gene ~ cluster", "logFC")
  de_select$Gene <- factor(as.character(de_select$Gene), levels = rev(sorted_genes))
  sorted_clusters <- umap_sorted_rows(de_select, "cluster ~ Gene", "logFC")
  de_select$cluster <- factor(as.character(de_select$cluster), levels = rev(sorted_clusters))
  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(de_select$logFC)), max(abs(de_select$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(de_select$signed_p), max(de_select$signed_p), length.out = 11))
  color_values <- seq(min(de_select$logFC), max(de_select$logFC), length.out = 11)

  p_mat <- ggplot() +
    geom_tile(
      data = de_select,
      mapping = aes(x = cluster, y = Gene, fill = logFC)
    ) +
    geom_point(
      data = de_select %>% filter(adj.P.Val < 0.05),
      mapping = aes(x = cluster, y = Gene),
      # shape = 21, fill = "white", stroke = 0.2, size = 2
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    facet_grid(rows = vars(cat), scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 10),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  cluster_colors <- mpn65
  names(cluster_colors) <- as.character(seq_along(cluster_colors))
  p_bar <- ggplot(de_select %>% select(cluster, celltype) %>% unique) +
    geom_tile(
      aes(y = 1, x = cluster, fill = celltype)
    ) +
    scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    scale_fill_manual(values = mpn65, guide = "none") +
    theme(
      axis.text.x.top = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5)
    )
  layout <- "
A
B
"
  p_mat <- p_mat + theme(
    # axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(l = 0)
  )
  p <- p_bar + p_mat + plot_layout(
    design = layout,
    widths = c(1, 1),
    heights = c(1, length(unique(de_select$Gene)))
  )
  longest <- head(unique(de_select$cat[nchar(de_select$cat) == max(nchar(de_select$cat))]), 1)
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-top-product-drugtargets-foldchange-cluster",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    # scale = 1,
    width = length(unique(de_select$cluster)) * 0.17 + text_width,
    height = length(unique(de_select$Gene)) * 0.3 + 1,
    units = "in",
    dpi = 300
  )

}

{

  analyses <- c(
    CD8 = "a12_4_4_t4_cd8_1_2",
    CD4 = "a12_4_4_t4_cd4_2_2",
    M   = "a12_4_4_m3_2",
    B   = "a12_4_4_b5_1_3",
    E   = "n3_2"
  )
  cluster_to_group <- do.call(c, lapply(names(analyses), function(this_name) {
    retval <- get_cluster_groups(analyses[[this_name]])
    names(retval) <- sprintf("%s-%s", this_name, names(retval))
    retval
  }))

  de_select$cluster_group <- cluster_to_group[as.character(de_select$cluster)]

  x1 <- readxl::read_excel("data/drug-targets.xlsx")
  gene_to_cat <- with(
    x1 %>% select(Category, Target) %>% unique(),
    split(Category, Target)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  rm(x1)
  de_select <- de %>% filter(Gene %in% names(gene_to_cat))
  #
  de_select$cat <- gene_to_cat[de_select$Gene]
  # filter_genes <- unique((de_select %>% filter(adj.P.Val < 0.05))$Gene)
  # de_select <- de_select %>% filter(Gene %in% filter_genes)
  de_select <- de_select %>% group_by(Gene) %>% mutate(AveExpr = scale(AveExpr))
  #
  set.seed(1)
  sorted_genes <- umap_sorted_rows(de_select, "Gene ~ cluster", "AveExpr")
  de_select$Gene <- factor(as.character(de_select$Gene), levels = rev(sorted_genes))
  sorted_clusters <- umap_sorted_rows(de_select, "cluster ~ Gene", "AveExpr")
  de_select$cluster <- factor(as.character(de_select$cluster), levels = rev(sorted_clusters))
  de_select$cat <- factor(de_select$cat,
    c("S1P", "Fractalkine", "Cytokine", "Integrin", "JAK", "Other"))
  de_select %<>% ungroup

  de_select$celltype <- factor(de_select$celltype, c("CD8", "CD4", "E", "M", "B"))

  de_select$AveExpr_orig <- de_select$AveExpr

  de_select$AveExpr <- de_select$AveExpr_orig
  limit_max <- function(x, lim) {
    x[x > lim] <- lim
    return(x)
  }
  de_select$AveExpr <- limit_max(de_select$AveExpr, 3)

  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(de_select$AveExpr)), max(abs(de_select$AveExpr)), length.out = 11),
    # rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
    rev(RColorBrewer::brewer.pal(name = "PRGn", n = 11))
    # rev(RColorBrewer::brewer.pal(name = "PuOr", n = 11))
  )
  # color_values <- scales::rescale(seq(min(de_select$signed_p), max(de_select$signed_p), length.out = 11))
  color_values <- seq(min(de_select$AveExpr), max(de_select$AveExpr), length.out = 11)
  #
  p_mat <- ggplot() +
    geom_tile(
      data = de_select,
      mapping = aes(x = cluster, y = Gene, fill = AveExpr)
    ) +
    # geom_point(
    #   data = de_select %>% filter(adj.P.Val < 0.05),
    #   mapping = aes(x = cluster, y = Gene),
    #   # shape = 21, fill = "white", stroke = 0.2, size = 2
    #   shape = 19, color = "white"
    # ) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0),
      labels = function(x) str_split_fixed(x, "-", 2)[,2]
    ) +
    #facet_grid(
    #  cols = vars(celltype), rows = vars(cat), scales = "free", space = "free"
    #  #, switch = "y"
    #) +
    ggh4x::facet_nested(cols = vars(celltype, cluster_group), rows = vars(cat), scales = "free", space = "free") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      name = "Scaled expression across cell clusters",
      guide = guide_colorbar(barwidth = 10),
      breaks = pretty_breaks(4)
      # labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x.top = element_text(size = 12, angle = 90, hjust = 0, vjust = 0.5),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 0),
      # strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  cluster_colors <- mpn65
  names(cluster_colors) <- as.character(seq_along(cluster_colors))
  p_bar <- ggplot(de_select %>% select(cluster, celltype, cluster_group) %>% unique) +
    geom_tile(
      aes(y = 1, x = cluster, fill = celltype)
    ) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0),
      labels = function(x) str_split_fixed(x, "-", 2)[,2]
    ) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    scale_fill_manual(values = mpn65, guide = "none") +
    # facet_grid(cols = vars(celltype), scales = "free", space = "free") +
    ggh4x::facet_nested(cols = vars(celltype, cluster_group), scales = "free", space = "free") +
    theme(
      # axis.text.x.top = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
      # strip.text.x.top = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      axis.text.x.top = element_blank(),
      axis.ticks.x = element_blank()
    )
  layout <- "
A
B
"
  p_mat <- p_mat + theme(
    # axis.text.y = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(l = 0)
  )
  p <- p_bar + p_mat + plot_layout(
    design = layout,
    widths = c(1, 1),
    heights = c(1, length(unique(de_select$Gene)))
  )
  #
  get_longest <- function(x) {
    x <- as.character(x)
    head(unique(x[nchar(x) == max(nchar(x))]), 1)
  }
  #
  longest <- get_longest(de_select$cat)
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-top-product-drugtargets-zscore-cluster",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    # scale = 1,
    # width = length(unique(de_select$cluster)) * 0.17 + text_width,
    width = length(unique(de_select$cluster)) * 0.21 + text_width,
    height = length(unique(de_select$Gene)) * 0.33 + 1,
    units = "in",
    dpi = 300
  )

}


# Product analysis with curated pairs, cell clusters (B1, B2, etc.)
########################################################################

{

  my_pairs <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  my_pairs$target_ens <- unname(symbol_to_ensembl[my_pairs$g_target])
  my_pairs$source_ens <- unname(symbol_to_ensembl[my_pairs$g_source])
  # my_pairs[which(is.na(my_pairs$source_ens)),]
  # my_pairs[which(is.na(my_pairs$target_ens)),]
  my_pairs <- my_pairs %>% filter(target_ens %in% rownames(counts), source_ens %in% rownames(counts))
  my_pairs <- my_pairs[!duplicated(with(my_pairs, paste(target_ens, source_ens))),]

  #
  source("R/covarying-composition.R")
  gg2_cur_file <- "cache/gg2_cur.qs"
  if (file.exists(gg2_cur_file)) {
    gg2_cur <- qread(gg2_cur_file)
  } else {
    gg2_cur <- gene_gene_products(
      counts            = counts,
      cluster           = obs$leiden,
      donor             = obs$donor,
      case              = obs$case,
      source_ens        = my_pairs$target_ens,
      target_ens        = my_pairs$source_ens,
      min_donors        = 10,
      min_cells         = 25,
      ensembl_to_symbol = ensembl_to_symbol
    )
    gg2_cur$limma %<>% mutate(epi = c1 == "E" | c2 == "E", inter = substr(c1, 1, 1) != substr(c2, 1, 1))
    # gg2_cur$limma %<>% group_by(c1, c2) %>% mutate(fdr = p.adjust(P.Value, method = "fdr"))
    qsave(gg2_cur, gg2_cur_file)
  }
  gg2_cur$limma %>% select(
    c1, g1, c2, g2, logFC, P.Value, p1, p2
  ) %>% arrange(P.Value)
  stopifnot(sum(duplicated(gg2_cur$limmma$id)) == 0)

  #
  d2 <- gg2_cur$limma %>%
    filter(p1 > 0.01, p2 > 0.01) %>%
    group_by(c1, c2) %>%
    mutate(
      fdr = p.adjust(P.Value, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
      gene_pair = glue("{g2} - {g1}"),
      cell_pair = glue("{c2} - {c1}")
    )
  d2 <- left_join(
    x = d2,
    y = my_pairs %>% select(cat, source_ens, target_ens),
    by = c("e1" = "source_ens", "e2" = "target_ens")
  )
  stopifnot(sum(duplicated(d2$id)) == 0)

  d2 %>% group_by(cell_pair) %>% summarize(n = sum(fdr < 0.05)) %>% filter(n > 1) %>% arrange(-n)

  # my_gene_pairs <- unique((d2 %>% filter(fdr < 0.1, epi))$gene_pair)
  my_gene_pairs <- unique((d2 %>% filter(fdr < 0.05))$gene_pair)
  length(my_gene_pairs)
  d2 <- d2 %>% filter(gene_pair %in% my_gene_pairs)
  # gene_pairs <- (
  #   d2 %>% group_by(gene_pair) %>% summarize(min_pval = min(P.Value)) %>% arrange(-min_pval)
  # )$gene_pair
  # d2$gene_pair <- factor(as.character(d2$gene_pair), levels = gene_pairs)
  set.seed(1)
  sorted_gene_pairs <- umap_sorted_rows(d2, "gene_pair ~ cell_pair", "logFC")
  d2$gene_pair <- factor(as.character(d2$gene_pair), levels = rev(sorted_gene_pairs))
  # d2$c2 <- factor(as.character(d2$c2), c("E", "CT", "T", "M", "B"))

  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(d2$logFC)), max(abs(d2$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(d2$signed_p), max(d2$signed_p), length.out = 11))
  color_values <- seq(min(d2$logFC), max(d2$logFC), length.out = 11)
  p <- ggplot() +
    geom_tile(
      data = d2,
      mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
    ) +
    geom_point(
      data = d2 %>% filter(fdr < 0.05),
      mapping = aes(x = cell_pair, y = gene_pair),
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    # facet_wrap(cat ~ c2, ncol = 5, scales = "free_x") +
    facet_grid(cat ~ c2, scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      # name = bquote("log"[2]~ "Fold-Change"),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 15),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  longest <- unique(d2$cat[nchar(d2$cat) == max(nchar(d2$cat))])
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-top-product-selected",
    out_dir = file.path("results/a20/covarying-abundance/level2"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(d2$cell_pair)) * 0.4 + 1 + text_width,
    height = length(unique(d2$gene_pair)) * 0.3 + 1.1,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )

}

# Figure 7
{

  my_pairs <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  my_pairs$target_ens <- unname(symbol_to_ensembl[my_pairs$g_target])
  my_pairs$source_ens <- unname(symbol_to_ensembl[my_pairs$g_source])
  # my_pairs[which(is.na(my_pairs$source_ens)),]
  # my_pairs[which(is.na(my_pairs$target_ens)),]
  my_pairs <- my_pairs %>% filter(target_ens %in% rownames(counts), source_ens %in% rownames(counts))
  my_pairs <- my_pairs[!duplicated(with(my_pairs, paste(target_ens, source_ens))),]

  gg2 <- qread("results/a20/covarying-abundance/gg_level2_case.qs")
  stopifnot(sum(duplicated(gg2$limma$id)) == 0)

  ##
  #d2 <- gg2$limma %>%
  #  filter(p1 > 0.01, p2 > 0.01) %>%
  #  group_by(c1, c2) %>%
  #  mutate(
  #    fdr = p.adjust(P.Value, method = "BH")
  #  ) %>%
  #  ungroup() %>%
  #  mutate(
  #    gene_pair = glue("{g2} - {g1}"),
  #    cell_pair = glue("{c2} - {c1}")
  #  )
  #d2 <- left_join(
  #  x = d2,
  #  y = my_pairs %>% select(cat, source_ens, target_ens),
  #  by = c("e1" = "source_ens", "e2" = "target_ens")
  #)
  #stopifnot(sum(duplicated(d2$id)) == 0)
  ##
  #my_cluster <- "CT3"
  #d2 <- d2 %>% filter(c1 == my_cluster | c2 == my_cluster)
  ## d2 <- d2 %>% filter(!str_detect(cell_pair, "B"))
  #d2_left <- d2 %>% filter(c2 == my_cluster)
  #d2_right <- d2 %>% filter(c1 == my_cluster)

  # d2 %>% group_by(cell_pair) %>% summarize(n = sum(fdr < 0.05)) %>% filter(n > 1) %>% arrange(-n)

  #for (side in c("left", "right")) {
  #  this_d2 <- d2_right
  #  if (side == "left") {
  #    this_d2 <- d2_left
  #  }
  #  #
  #  # my_gene_pairs <- unique((d2 %>% filter(fdr < 0.1, epi))$gene_pair)
  #  my_gene_pairs <- unique((this_d2 %>% filter(fdr < 0.01))$gene_pair)
  #  length(my_gene_pairs)
  #  this_d2 <- this_d2 %>% filter(gene_pair %in% my_gene_pairs)
  #  # gene_pairs <- (
  #  #   this_d2 %>% group_by(gene_pair) %>% summarize(min_pval = min(P.Value)) %>% arrange(-min_pval)
  #  # )$gene_pair
  #  # this_d2$gene_pair <- factor(as.character(this_d2$gene_pair), levels = gene_pairs)
  #  set.seed(1)
  #  sorted_gene_pairs <- umap_sorted_rows(this_d2, "gene_pair ~ cell_pair", "logFC")
  #  this_d2$gene_pair <- factor(as.character(this_d2$gene_pair), levels = rev(sorted_gene_pairs))
  #  # this_d2$c2 <- factor(as.character(this_d2$c2), c("E", "CT", "T", "M", "B"))
  #  #
  #  col_fun <- circlize::colorRamp2(
  #    seq(-max(abs(this_d2$logFC)), max(abs(this_d2$logFC)), length.out = 11),
  #    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  #  )
  #  # color_values <- scales::rescale(seq(min(this_d2$signed_p), max(this_d2$signed_p), length.out = 11))
  #  color_values <- seq(min(this_d2$logFC), max(this_d2$logFC), length.out = 11)
  #  p <- ggplot() +
  #    geom_tile(
  #      data = this_d2,
  #      mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
  #    ) +
  #    geom_point(
  #      data = this_d2 %>% filter(fdr < 0.05),
  #      mapping = aes(x = cell_pair, y = gene_pair),
  #      shape = 19, color = "white"
  #    ) +
  #    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  #    scale_x_discrete(
  #      position = "t", name = NULL, expand = c(0, 0)
  #    ) +
  #    # facet_wrap(cat ~ c2, ncol = 5, scales = "free_x") +
  #    # facet_grid(cat ~ c2, scales = "free", space = "free", switch = "y") +
  #    facet_grid(cat ~ ., scales = "free", space = "free", switch = "y") +
  #    scale_fill_gradientn(
  #      colors = col_fun(color_values),
  #      # name = bquote("log"[2]~ "Fold-Change"),
  #      name = "Fold-Change",
  #      guide = guide_colorbar(barwidth = 15),
  #      breaks = pretty_breaks(4),
  #      labels = function(x) fractional::fractional(2^x)
  #    ) +
  #    theme(
  #      axis.text.x = element_text(angle = 60, hjust = 0),
  #      axis.text.y = element_text(face = "italic"),
  #      strip.text.x = element_blank(),
  #      strip.text.y.left = element_text(angle = 0, hjust = 1),
  #      panel.spacing = unit(0.5, "lines"),
  #      legend.position = "bottom"
  #    )
  #  longest <- unique(this_d2$cat[nchar(this_d2$cat) == max(nchar(this_d2$cat))])
  #  text_width <- strwidth(longest, font = 12, units = 'in')
  #  my_ggsave(
  #    glue("{my_cluster}-{side}-heatmap-top-product-selected"),
  #    out_dir = file.path("results/a20/covarying-abundance/level2/by-cluster"),
  #    type = "pdf",
  #    plot = p,
  #    scale = 1,
  #    width = length(unique(this_d2$cell_pair)) * 0.4 + 1 + text_width,
  #    height = length(unique(this_d2$gene_pair)) * 0.4 + 1.1,
  #    units = "in",
  #    dpi = 300,
  #    limitsize = FALSE
  #  )
  #}

  #
  d2 <- gg2$limma %>%
    filter(p1 > 0.01, p2 > 0.01) %>%
    group_by(c1, c2) %>%
    mutate(
      fdr = p.adjust(P.Value, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
      gene_pair = glue("{g2} - {g1}"),
      cell_pair = glue("{c2} - {c1}")
    )
  d2 <- left_join(
    x = d2,
    y = my_pairs %>% select(cat, source_ens, target_ens),
    by = c("e1" = "source_ens", "e2" = "target_ens")
  )
  stopifnot(sum(duplicated(d2$id)) == 0)

  my_clusters <- unique(c(d2$c1, d2$c2))

  if (FALSE) {
    # results/a20/covarying-abundance/level2/by-cluster
    for (my_cluster in my_clusters) {
  #     my_cluster <- "CT3"
      this_d2 <- d2 %>% filter(c1 == my_cluster | c2 == my_cluster)
      #
      # my_gene_pairs <- unique((d2 %>% filter(fdr < 0.1, epi))$gene_pair)
      # my_gene_pairs <- unique((this_d2 %>% filter(fdr < 0.01))$gene_pair)
      my_gene_pairs <- unique((this_d2 %>% group_by(cell_pair) %>% top_n(n = 2, wt = B))$gene_pair)
      length(my_gene_pairs)
      this_d2 <- this_d2 %>% filter(gene_pair %in% my_gene_pairs)
      # gene_pairs <- (
      #   this_d2 %>% group_by(gene_pair) %>% summarize(min_pval = min(P.Value)) %>% arrange(-min_pval)
      # )$gene_pair
      # this_d2$gene_pair <- factor(as.character(this_d2$gene_pair), levels = gene_pairs)
      set.seed(1)
      sorted_gene_pairs <- umap_sorted_rows(this_d2, "gene_pair ~ cell_pair", "logFC")
      this_d2$gene_pair <- factor(as.character(this_d2$gene_pair), levels = rev(sorted_gene_pairs))
      # this_d2$c2 <- factor(as.character(this_d2$c2), c("E", "CT", "T", "M", "B"))
      sorted_cell_pairs <- umap_sorted_rows(this_d2, "cell_pair ~ gene_pair", "logFC")
      this_d2$cell_pair <- factor(as.character(this_d2$cell_pair), levels = rev(sorted_cell_pairs))
      #
      col_fun <- circlize::colorRamp2(
        seq(-max(abs(this_d2$logFC)), max(abs(this_d2$logFC)), length.out = 11),
        rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
      )
      # color_values <- scales::rescale(seq(min(this_d2$signed_p), max(this_d2$signed_p), length.out = 11))
      color_values <- seq(min(this_d2$logFC), max(this_d2$logFC), length.out = 11)
      #
      this_d2$ct1 <- str_remove(this_d2$c1, "\\d+$")
      this_d2$ct2 <- str_remove(this_d2$c2, "\\d+$")
      this_d2$side <- this_d2$c2 != my_cluster
      this_d2$side_label <- with(this_d2, ifelse(side, ct2, ct1))
      #
      p_bar <- ggplot() +
        geom_tile(
          data = this_d2 %>% filter(fdr < 0.05),
          # data = this_d2,
          mapping = aes(x = cell_pair, y = 1, fill = side_label)
        ) +
        facet_grid(.~ side, scales = "free", space = "free", switch = "y") +
        scale_x_discrete(
          position = "t", name = NULL, expand = c(0, 0)
        ) +
        scale_fill_manual(values = mpn65, name = NULL) +
        theme(
          axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          strip.text.x = element_blank(),
          legend.position = "top"
        ) +
        labs(x = NULL, y = NULL, title = glue("{my_cluster}"))
      p <- ggplot() +
        geom_tile(
          # data = this_d2,
          data = this_d2 %>% filter(fdr < 0.05),
          mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
        ) +
        geom_point(
          data = this_d2 %>% filter(fdr < 0.05),
          mapping = aes(x = cell_pair, y = gene_pair),
          shape = 19, color = "white"
        ) +
        scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
        scale_x_discrete(
          position = "t", name = NULL, expand = c(0, 0)
        ) +
        # facet_wrap(cat ~ c2, ncol = 5, scales = "free_x") +
        # facet_grid(cat ~ c2, scales = "free", space = "free", switch = "y") +
        facet_grid(cat ~ side, scales = "free", space = "free", switch = "y") +
        scale_fill_gradientn(
          colors = col_fun(color_values),
          # name = bquote("log"[2]~ "Fold-Change"),
          name = "Fold-Change",
          guide = guide_colorbar(barwidth = 15),
          breaks = pretty_breaks(4),
          labels = function(x) fractional::fractional(2^x)
        ) +
        theme(
          # axis.text.x = element_text(angle = 60, hjust = 0),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(face = "italic"),
          strip.text.x = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          panel.spacing = unit(0.5, "lines"),
          legend.position = "bottom"
        )
      longest <- unique(this_d2$cat[nchar(this_d2$cat) == max(nchar(this_d2$cat))])
      text_width <- strwidth(longest, font = 12, units = 'in')
      my_ggsave(
        glue("{my_cluster}-both-heatmap-top-product-selected"),
        out_dir = file.path("results/a20/covarying-abundance/level2/by-cluster"),
        type = "pdf",
        plot = p_bar + p + plot_layout(ncol = 1, heights = c(1, length(unique(this_d2$gene_pair)))),
        scale = 1,
        width = length(unique(this_d2$cell_pair)) * 0.28 + 1 + text_width,
        height = length(unique(this_d2$gene_pair)) * 0.35 + 1.1,
        units = "in",
        dpi = 300,
        limitsize = FALSE
      )
    }
  }

  d3 <- gg3$limma %>% ungroup
  my_clusters <- c("CT3", "CT11", "CT7", "CT5", "CT2")
  my_gene_pairs <- c(
    "PDCD1 - CD274",
    "PDCD1 - PDCD1LG2",
    "HAVCR2 - CEACAM1",
    "CXCR3 - CXCL9",
    "CXCR3 - CXCL10",
    "CXCR3 - CXCL11",
    "ITGB2 - ICAM1",
    "ITGB2 - ICAM2",
    "ITGB2 - ICAM3",
    "IL17A - IL17RA",
    "IL17A - IL17RC",
    "IL26 - IL10RB",
    "IL26 - IL20RA",
    "CXCL13 - CXCR5"
  )
  this_d3 <- d3 %>%
    filter(
      c1 %in% my_clusters | c2 %in% my_clusters,
      c1 != c2,
      xor(str_detect(c1, "CT"), str_detect(c2, "CT"))
    ) %>%
    mutate(cell_pair = glue("{c2} - {c1}"), gene_pair = glue("{g2} - {g1}"))
  # this_d3 <- this_d3 %>% filter(g1 %in% my_genes | g2 %in% my_genes)
#   table((
#     this_d3 %>% select(c1, c2, g1, g2) %>% filter(g1 == "ICAM1")
#   )$g1)
  this_d3 <- this_d3 %>% filter(gene_pair %in% my_gene_pairs)
  # nrow(this_d3)
  # this_d3 %>% filter(g1 == "PDCD1" | g2 == "PDCD1") %>% count(gene_pair)
  # this_d3 %>% filter(g1 == "CEACAM1" | g2 == "CEACAM1") %>% count(gene_pair)

  #
  # my_cell_pairs <- unique((this_d3 %>% group_by(gene_pair) %>% top_n(n = 15, wt = B))$cell_pair)
  # this_d3 <- this_d3 %>% filter(cell_pair %in% my_cell_pairs)
  #
  gene_pair_cat <- c(
    "PDCD1 - CD274"    = "CD8 T Co-Inhibitory",
    "PDCD1 - PDCD1LG2" = "CD8 T Co-Inhibitory",
    "HAVCR2 - CEACAM1" = "CD8 T Co-Inhibitory",
    "CXCR3 - CXCL9"    = "CD8 T Homing Receptors",
    "CXCR3 - CXCL10"   = "CD8 T Homing Receptors",
    "CXCR3 - CXCL11"   = "CD8 T Homing Receptors",
    "ITGB2 - ICAM1"    = "CD8 T Homing Receptors",
    "ITGB2 - ICAM2"    = "CD8 T Homing Receptors",
    "ITGB2 - ICAM3"    = "CD8 T Homing Receptors",
    "IL17A - IL17RA"   = "CD8 T Secreted Factors",
    "IL17A - IL17RC"   = "CD8 T Secreted Factors",
    "IL26 - IL10RB"    = "CD8 T Secreted Factors",
    "IL26 - IL20RA"    = "CD8 T Secreted Factors",
    "CXCL13 - CXCR5"   = "CD8 T Secreted Factors"
  )
  this_d3$cat <- gene_pair_cat[this_d3$gene_pair]
  #
  this_d3$cell_pair <- as.character(this_d3$cell_pair) %>%
    str_replace_all("CT", "CD8-") %>%
    str_replace_all("T(\\d)", "CD4-\\1") %>%
    str_replace_all("M(\\d)", "MP-\\1") %>%
    str_replace_all("E(\\d)", "E-\\1")
  names(table(this_d3$cell_pair))

  #
  # my_gene_pairs <- unique((d3 %>% filter(fdr < 0.1, epi))$gene_pair)
  # my_gene_pairs <- unique((this_d3 %>% filter(fdr < 0.01))$gene_pair)
  # my_gene_pairs <- unique((this_d3 %>% group_by(cell_pair) %>% top_n(n = 2, wt = B))$gene_pair)
  # length(my_gene_pairs)
  # this_d3 <- this_d3 %>% filter(gene_pair %in% my_gene_pairs)
  # gene_pairs <- (
  #   this_d3 %>% group_by(gene_pair) %>% summarize(min_pval = min(P.Value)) %>% arrange(-min_pval)
  # )$gene_pair
  # this_d3$gene_pair <- factor(as.character(this_d3$gene_pair), levels = gene_pairs)
  set.seed(1)
  # sorted_gene_pairs <- umap_sorted_rows(this_d3, "gene_pair ~ cell_pair", "logFC")
  sorted_gene_pairs <- my_gene_pairs
  this_d3$gene_pair <- factor(as.character(this_d3$gene_pair), levels = rev(sorted_gene_pairs))
  # this_d3$c2 <- factor(as.character(this_d3$c2), c("E", "CT", "T", "M", "B"))
  # sorted_cell_pairs <- umap_sorted_rows(this_d3, "cell_pair ~ gene_pair", "logFC")
  # this_d3$cell_pair <- factor(as.character(this_d3$cell_pair), levels = rev(sorted_cell_pairs))
  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(this_d3$logFC)), max(abs(this_d3$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(this_d3$signed_p), max(this_d3$signed_p), length.out = 11))
  color_values <- seq(min(this_d3$logFC), max(this_d3$logFC), length.out = 11)
  this_d3$cell_pair <- str_replace_all(this_d3$cell_pair, "\\bM\\b", "Myeloid cells")
  this_d3$cell_pair <- str_replace_all(this_d3$cell_pair, "\\bB\\b", "B cells")
  this_d3$cell_pair <- str_replace_all(this_d3$cell_pair, "\\bT\\b", "CD4 T cells")
  this_d3$cell_pair <- str_replace_all(this_d3$cell_pair, "\\bCD8\\b", "CD8 T")
  this_d3$cell_pair <- str_replace_all(this_d3$cell_pair, "\\bSecretory cells\\b", "Secretory epithelial cells")
  this_d3$cell_pair <- str_replace_all(this_d3$cell_pair, "\\bMesenchymal cells\\b", "Fibroblasts")
  this_d3$cell_pair <- str_replace(this_d3$cell_pair, "epithelial cells", "epithelial")
  this_d3$cell_pair <- str_replace(this_d3$cell_pair, "absorptive", "abs.")
  this_d3$cell_pair
  #
  this_d3$ct1 <- str_remove(this_d3$c1, "\\d+$")
  this_d3$ct2 <- str_remove(this_d3$c2, "\\d+$")
  this_d3$side <- this_d3$c1 %in% my_clusters
  this_d3$side_label <- with(this_d3, ifelse(side, ct2, ct1))
  # table(this_d3$cell_pair)
  # table(this_d3$gene_pair)
  "CD274" %in% this_d3$g1
  "PDCD1LG2" %in% this_d3$g1

  this_d3$cd8_cell <- ""
  ix <- str_detect(this_d3$c1, "^CT")
  this_d3$cd8_cell[ix] <- this_d3$c1[ix]
  ix <- str_detect(this_d3$c2, "^CT")
  this_d3$cd8_cell[ix] <- this_d3$c2[ix]
  # table(this_d3$cd8_cell)
  my_levels <- c(
    naturalsort(unique(this_d3$cell_pair[this_d3$side])),
    naturalsort(unique(this_d3$cell_pair[!this_d3$side]))
  )
  this_d3$cell_pair <- factor(this_d3$cell_pair, my_levels)

  this_d3 <- this_d3[!this_d3$side,]
  this_d3$cell_pair1 <- str_split_fixed(this_d3$cell_pair, " - ", 2)[,1]
  this_d3$cell_pair2 <- str_split_fixed(this_d3$cell_pair, " - ", 2)[,2]
  this_d3$cell_pair1 <- naturalfactor(this_d3$cell_pair1)
  this_d3$cell_pair2 <- factor(this_d3$cell_pair2, rev(c(
    "CD4 T cells",
    "B cells",
    "Myeloid cells",
    "Immature epithelial",
    "Absorptive epithelial",
    "Mature abs. epithelial",
    "Secretory epithelial",
    "Endothelial cells",
    "Fibroblasts"
  )))
  this_d3 <- this_d3[with(this_d3, order(cell_pair1, cell_pair2)),]
  this_d3$cell_pair <- factor(
    as.character(this_d3$cell_pair),
    as.character(unique(this_d3$cell_pair))
  )

  this_d3$cd8_cell <- factor(
    as.character(this_d3$cd8_cell),
    c("CT11", "CT3", "CT7", "CT2", "CT5")
  )

  p_bar <- ggplot() +
    geom_tile(
      data = this_d3,
      mapping = aes(y = cell_pair, x = 1, fill = side_label)
    ) +
    facet_grid(rows = vars(cd8_cell), scales = "free", space = "free", switch = "y") +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      position = "r", name = NULL, expand = c(0, 0)
    ) +
    scale_fill_manual(values = mpn65, name = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "none",
      plot.background = element_blank()
    ) +
    labs(x = NULL, y = NULL)
  p <- ggplot() +
    geom_tile(
      data = this_d3,
      mapping = aes(y = cell_pair, x = gene_pair, fill = logFC)
    ) +
    geom_point(
      data = this_d3 %>% filter(fdr < 0.05),
      mapping = aes(y = cell_pair, x = gene_pair),
      shape = 19, color = "white"
    ) +
    scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
    scale_y_discrete(
      position = "l", name = NULL, expand = c(0, 0)
    ) +
    facet_grid(rows = vars(cd8_cell), cols = vars(cat), scales = "free", space = "free", switch = "x") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      # name = bquote("log"[2]~ "Fold-Change"),
      name = "FC",
      guide = guide_colorbar(barwidth = 15),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x.top = element_text(face = "italic", angle = 60, hjust = 0),
      # axis.text.x.top = element_text(face = "italic", angle = 90, hjust = 0, vjust = 0.5),
      strip.text.y = element_blank(),
      # strip.text.x.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_blank(),
      legend.position = "top"
    )
  longest <- head(unique(this_d3$cat[nchar(this_d3$cat) == max(nchar(this_d3$cat))]), 1)
  text_width <- strwidth(longest, font = 12, units = 'in')
  p_both <- p + p_bar + plot_layout(ncol = 2, widths = c(length(unique(this_d3$gene_pair)), 0.5))
  my_ggsave(
    glue("cd8-heatmap-top-product-selected-h"),
    out_dir = file.path("results/a20/covarying-abundance/level2"),
    type = "pdf",
    plot = p_both,
    height = length(unique(this_d3$cell_pair)) * 0.28 + 1 + text_width,
    width = length(unique(this_d3$gene_pair)) * 0.35 + 4,
    units = "in",
    limitsize = FALSE
  )


  # p_bar <- ggplot() +
  #   geom_tile(
  #     # data = this_d3 %>% filter(fdr < 0.05),
  #     data = this_d3,
  #     mapping = aes(x = cell_pair, y = 1, fill = side_label)
  #   ) +
  #   facet_grid(.~ side, scales = "free", space = "free", switch = "y") +
  #   scale_x_discrete(
  #     position = "t", name = NULL, expand = c(0, 0)
  #   ) +
  #   scale_fill_manual(values = mpn65, name = NULL) +
  #   theme(
  #     axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0),
  #     axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank(),
  #     panel.spacing = unit(0.5, "lines"),
  #     strip.text.x = element_blank(),
  #     legend.position = "top"
  #   ) +
  #   labs(x = NULL, y = NULL)
  # p <- ggplot() +
  #   geom_tile(
  #     data = this_d3,
  #     # data = this_d3 %>% filter(fdr < 0.05),
  #     mapping = aes(x = cell_pair, y = gene_pair, fill = logFC)
  #   ) +
  #   geom_point(
  #     data = this_d3 %>% filter(fdr < 0.05),
  #     mapping = aes(x = cell_pair, y = gene_pair),
  #     shape = 19, color = "white"
  #   ) +
  #   scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  #   scale_x_discrete(
  #     position = "t", name = NULL, expand = c(0, 0)
  #   ) +
  #   # facet_wrap(cat ~ c2, ncol = 5, scales = "free_x") +
  #   # facet_grid(cat ~ c2, scales = "free", space = "free", switch = "y") +
  #   facet_grid(cat ~ side, scales = "free", space = "free", switch = "y") +
  #   scale_fill_gradientn(
  #     colors = col_fun(color_values),
  #     # name = bquote("log"[2]~ "Fold-Change"),
  #     name = "Fold-Change",
  #     guide = guide_colorbar(barwidth = 15),
  #     breaks = pretty_breaks(4),
  #     labels = function(x) fractional::fractional(2^x)
  #   ) +
  #   theme(
  #     # axis.text.x = element_text(angle = 60, hjust = 0),
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     axis.text.y = element_text(face = "italic"),
  #     strip.text.x = element_blank(),
  #     strip.text.y.left = element_text(angle = 0, hjust = 1),
  #     panel.spacing = unit(0.5, "lines"),
  #     legend.position = "bottom"
  #   )
  # longest <- head(unique(this_d3$cat[nchar(this_d3$cat) == max(nchar(this_d3$cat))]), 1)
  # text_width <- strwidth(longest, font = 12, units = 'in')
  # p_both <- p_bar / p + plot_layout(ncol = 1, heights = c(1, length(unique(this_d3$gene_pair))))
  # my_ggsave(
  #   glue("cd8-heatmap-top-product-selected"),
  #   out_dir = file.path("results/a20/covarying-abundance/level2"),
  #   type = "pdf",
  #   plot = p_both,
  #   width = length(unique(this_d3$cell_pair)) * 0.28 + 1 + text_width,
  #   height = length(unique(this_d3$gene_pair)) * 0.35 + 2.1,
  #   units = "in",
  #   limitsize = FALSE
  # )
  
}

# Differentially expressed genes in each cell type (E, CD8, CD4, M, B)
########################################################################

d1 <- gg1$limma %>%
  filter(p1 > 0.01, p2 > 0.01) %>%
	group_by(c1, c2) %>%
  mutate(
    fdr = p.adjust(P.Value, method = "BH")
	) %>%
  ungroup() %>%
  mutate(
    gene_pair = glue("{g2} - {g1}"),
    cell_pair = glue("{c2} - {c1}")
  )
my_pair_genes <- unique(c(d1$g1, d1$g2))

analyses <- c(
  "a12_4_4_t4_cd8_1_2",
  "a12_4_4_t4_cd4_2_2",
  "a12_4_4_m3_2",
  "a12_4_4_b5_1_3",
  "n3_2"
)
de <- rbindlist(lapply(analyses, function(analysis_name) {
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
  de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
  de <- fread(de_donor_file)
  de$celltype <- analysis_name
  return(de)
}))
de$celltype <- str_remove(de$celltype, "a12_4_4_")
to_celltype <- c("B", "M", "E", "CD4", "CD8")
names(to_celltype) <- c("b5_1_3", "m3_2", "n3_2", "t4_cd4_2_2", "t4_cd8_1_2")
de$celltype <- to_celltype[de$celltype]

{

  de_select <- de %>% filter(Gene %in% my_pair_genes)
  #
  x1 <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  x1 <- pivot_longer(x1, cols = c("g_target", "g_source"))
  gene_to_cat <- with(
    x1 %>% select(cat, value) %>% unique(),
    split(cat, value)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  rm(x1)
  #
  de_select$cat <- gene_to_cat[de_select$Gene]

  set.seed(1)
  sorted_genes <- umap_sorted_rows(de_select, "Gene ~ celltype", "logFC")
  de_select$Gene <- factor(as.character(de_select$Gene), levels = rev(sorted_genes))

  #
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(de_select$logFC)), max(abs(de_select$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(de_select$signed_p), max(de_select$signed_p), length.out = 11))
  color_values <- seq(min(de_select$logFC), max(de_select$logFC), length.out = 11)
  p <- ggplot() +
    geom_tile(
      data = de_select,
      mapping = aes(x = celltype, y = Gene, fill = logFC)
    ) +
    geom_point(
      data = de_select %>% filter(adj.P.Val < 0.05),
      mapping = aes(x = celltype, y = Gene),
      shape = 19, color = "white"
    ) +
    scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
    scale_x_discrete(
      position = "t", name = NULL, expand = c(0, 0)
    ) +
    facet_grid(rows = vars(cat), scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      name = "Fold-Change",
      guide = guide_colorbar(barwidth = 10),
      breaks = pretty_breaks(4),
      labels = function(x) fractional::fractional(2^x)
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 0),
      axis.text.y = element_text(face = "italic"),
      strip.text.x = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
  longest <- unique(de_select$cat[nchar(de_select$cat) == max(nchar(de_select$cat))])
  text_width <- strwidth(longest, font = 12, units = 'in')
  my_ggsave(
    "heatmap-top-product-selected-foldchange",
    out_dir = file.path("results/a20/covarying-abundance/level1"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = length(unique(de_select$celltype)) * 1 + 1 + text_width,
    height = length(unique(de_select$Gene)) * 0.3 + 1.1,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )

}





































































#for (i in 1:nrow(top_gg1_limma)) {
#  c1 <- top_gg1_limma$c1[i]
#  c2 <- top_gg1_limma$c2[i]
#  g1 <- top_gg1_limma$g1[i]
#  g2 <- top_gg1_limma$g2[i]
#  #
#  e1 <- names(which(ensembl_to_symbol == g1))
#  e2 <- names(which(ensembl_to_symbol == g2))
#  #
#  x <- left_join(
#    x = pb_meta %>% filter(cluster == c1),
#    y = d[symbol == g1 & cluster == c1][,c("cluster", "donor","percent")],
#    by = c("cluster", "donor")
#  )
#  y <- left_join(
#    x = pb_meta %>% filter(cluster == c2),
#    y = d[symbol == g2 & cluster == c2][,c("cluster", "donor","percent")],
#    by = c("cluster", "donor")
#  )
#  #
#  x <- left_join(x = x, y = y, by = c("donor", "case"))
#  x$percent.x[is.na(x$percent.x)] <- 0
#  x$percent.y[is.na(x$percent.y)] <- 0
#  #
#  x$log.percent.x <- log(x$percent.x)
#  x$log.percent.y <- log(x$percent.y)
#  x$log.percent.x[!is.finite(x$log.percent.x)] <- 0
#  x$log.percent.y[!is.finite(x$log.percent.y)] <- 0
#  #
#  x <- x %>% mutate(score = log.percent.x + log.percent.y)
#  # x$score[!is.finite(x$score)] <- 0
#  # x$is_case <- x$case == "Case"
#  # t.test(formula = score ~ is_case, data = x)
#  # t.test(y = x$score[x$is_case], x = x$score[!x$is_case])
#  # summary(lm(formula = score ~ is_case, data = x))
#  # gg1$limma %>%
#  # arrange(P.Value) %>%
#  # filter(g1 == "PDCD1", g2 == "CD274", c1 == "CT", c2 == "E") %>% t
#  #
#  p <- ggplot() +
#    # geom_quasirandom(
#    #   data = x,
#    #   mapping = aes(y = case, x = exp(score), fill = case),
#    #   shape = 21, size = 5, stroke = 0.3, bandwidth = 0.5, nbins = 2,
#    #   dodge.width = 1, width = 0.3,
#    #   groupOnX = FALSE
#    # ) +
#    # scale_x_log10() +
#    # annotation_logticks(sides = "b", size = 0.3) +
#    geom_point(
#      data = x,
#      mapping = aes(x = percent.x, y = percent.y, fill = case),
#      shape = 21, size = 5, stroke = 0.3
#    ) +
#    # scale_fill_manual( values = okabe(8), guide = "none") +
#    scale_fill_manual( values = okabe(8)) +
#    labs(x = NULL, y = NULL, title = glue("{c1} {g1}, {c2} {g2}"))
#  my_ggsave(
#    glue("{c1}_{c2}_{g1}_{g2}"),
#    out_dir = "results/a20/covarying-abundance/level1/dots",
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 5,
#    height = 3,
#    units = "in", dpi = 300
#  )
#}


# Spearman correlation analysis
########################################################################

source("R/covarying-composition.R")

cc1_file <- "results/a20/covarying-abundance/cc_level1.qs"

if (file.exists(cc1_file)) {
  cc1 <- qread(cc1_file)
} else {
  ix_omni <- which(
    ensembl_to_symbol[rownames(log2cpm)] %in% c(omni$source_genesymbol, omni$target_genesymbol)
  )
  cc1 <- covarying_composition(
    log2cpm           = log2cpm[ix_omni,],
    cluster           = obs$dataset,
    donor             = obs$donor,
    min_donors        = 10,
    min_cells         = 25,
    source_genes      = omni$source_genesymbol,
    target_genes      = omni$target_genesymbol,
    ensembl_to_symbol = ensembl_to_symbol
  )
  cc1$correlations %<>% mutate(epi = c1 == "E" | c2 == "E", inter = substr(c1, 1, 1) != substr(c2, 1, 1))
  cc1$percents <- left_join(cc1$percents, obs %>% select(donor, case) %>% unique, by = "donor")
  cc1$correlations %<>% group_by(c1, c2) %>% mutate(fdr = p.adjust(p.value, method = "fdr"))
  # Drop duplicates
  cc1$correlations <- cc1$correlations %>%
    rowwise() %>%
    mutate(c_rank = diff(rank(c(c1, c2)))) %>%
    mutate(
      id = ifelse(c_rank == 1, glue("{c1} {g1} {c2} {g2}"), glue("{c2} {g2} {c1} {g1}"))
    ) %>%
    select(-c_rank)
  cc1$correlations <- cc1$correlations[!duplicated(cc1$correlations$id),]
  qsave(cc1, "results/a20/covarying-abundance/cc_level1.qs")
}

cc1_edges <- cc1$correlations %>% group_by(c1, c2) %>%
  filter(p1 > 0.01, p2 > 0.01) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  summarize(n_signif = sum(fdr < 0.05), .groups = "drop") %>% arrange(-n_signif)
cc1_edges <- cc1_edges %>%
  filter(n_signif > 0)
for (i in seq_len(nrow(cc1_edges))) {
  sorted <- naturalsort(c(cc1_edges$c1[i], cc1_edges$c2[i]))
  cc1_edges$c1[i] <- sorted[1]
  cc1_edges$c2[i] <- sorted[2]
}
cc1_edges %<>% group_by(c1, c2) %>% summarize(n_signif = sum(n_signif), .groups = "drop")


p <- ggplot(cc1$correlations) +
  aes(estimate, -log10(p.value), color = fdr < 0.05) +
  geom_scattermore() +
  geom_point(data = cc1$correlations %>% filter(fdr < 0.05), size = 0.5, color = "red") +
  scale_color_manual(
    name = "FDR < 0.05",
    values = c("TRUE" = "red", "FALSE" = "grey50")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  annotate(
    geom = "text",
    label = glue("{sum(cc1$correlations$fdr < 0.05)} / {comma(nrow(cc1$correlations))} gene-and-celltype pairs"),
    x = -Inf,
    y = Inf,
    hjust = -0.05, vjust = 1.5
  ) +
  theme(
    legend.position = "top",
    plot.title = element_textbox(width = unit(1, "npc"))
  ) +
  labs(
    x = "Spearman correlation",
    y = bquote("-log"[10]~"P"),
    title = "Correlations of percent of cells expressing pairs of genes"
  )
ggsave(
  filename = "results/a20/covarying-abundance/level1/spearman_volcano.pdf",
  plot = p,
  scale = 1, width = 4, height = 4, units = "in", dpi = 300
)

cc1_edges$c1[cc1_edges$c1 == "CT"] <- "CD8 T"
cc1_edges$c1[cc1_edges$c1 == "T"] <- "CD4 T"
cc1_edges$c1[cc1_edges$c1 == "E"] <- "Epi/Mes"
cc1_edges$c1[cc1_edges$c1 == "M"] <- "MP"
cc1_edges$c2[cc1_edges$c2 == "CT"] <- "CD8 T"
cc1_edges$c2[cc1_edges$c2 == "T"] <- "CD4 T"
cc1_edges$c2[cc1_edges$c2 == "E"] <- "Epi/Mes"
cc1_edges$c2[cc1_edges$c2 == "M"] <- "MP"

cc1_total <- length(unique(paste(cc1$correlations$g1, cc1$correlations$g2)))

cc1_graph <- graph_from_data_frame(
  cc1_edges,
  directed = FALSE
)
p <- ggraph(cc1_graph, layout = "kk") +
  geom_edge_link(aes(width = n_signif, color = n_signif), alpha = 1) +
  geom_node_point(size = 4) +
  scale_edge_width(range = c(0.3, 4)) +
  scale_edge_color_gradientn(colours = scico::scico(n = 11, palette = "lajolla", direction = -1)[1:9]) +
  geom_node_text(aes(label = name), nudge_x = c(0.2, 0.4, 0.5, 0.3, 0.4)) +#, repel = TRUE) +
  labs(
    title = glue("Correlated gene pairs"),
    subtitle = glue("of {comma(cc1_total)} total gene pairs"),
    edge_width = "FDR < 5%", edge_color = "FDR < 5%"
  ) +
  theme(
    plot.title = element_textbox(width = unit(1.1, "npc"))
  ) +
  # theme_graph() +
  coord_equal()
ggsave(
  filename = "results/a20/covarying-abundance/level1/spearman_graph.pdf",
  plot = p,
  scale = 1, width = 4, height = 3, units = "in", dpi = 300
)

cc1 <- qread("results/a20/covarying-abundance/cc_level1.qs")

fwrite(
  x = cc1$correlations %>% mutate_if(is.numeric, signif),
  file = "results/a20/covarying-abundance/ccc-spearman-level1.tsv",
  sep = "\t"
)

{
  wb <- openxlsx::createWorkbook()
  fname <- "results/a20/covarying-abundance/ccc-spearman-level1.xlsx"
  unlink(fname)
  #
  this_sheet <- "Sheet 1"
  openxlsx::addWorksheet(wb, this_sheet)
  openxlsx::writeDataTable(
    wb,
    this_sheet,
    x = cc1$correlations %>% mutate_if(is.numeric, signif),
    rowNames = FALSE, tableStyle = "TableStyleLight1"
  )
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}

# cc1$correlations %>% group_by(c1, c2) %>%
#   filter(p1 > 0.01, p2 > 0.01) %>%
#   mutate(fdr = p.adjust(p.value, method = "BH")) %>%
#   filter(fdr < 0.05, epi)

cc1_case <- rbindlist(lapply(c("Case", "Control"), function(this_case) {
  ix_case <- obs$case == this_case
  res <- covarying_composition(
    log2cpm           = log2cpm[ix_omni,ix_case],
    cluster           = obs$dataset[ix_case],
    donor             = obs$donor[ix_case],
    min_donors        = 10,
    min_cells         = 25,
    source_genes      = omni$source_genesymbol,
    target_genes      = omni$target_genesymbol,
    ensembl_to_symbol = ensembl_to_symbol
  )
  res$correlations$case <- this_case
  res$correlations %<>% mutate(epi = c1 == "E" | c2 == "E", inter = substr(c1, 1, 1) != substr(c2, 1, 1))
  # Drop duplicates
  res$correlations <- res$correlations %>%
    rowwise() %>%
    mutate(c_rank = diff(rank(c(c1, c2)))) %>%
    mutate(
      id = ifelse(c_rank == 1, glue("{c1} {g1} {c2} {g2}"), glue("{c2} {g2} {c1} {g1}"))
    ) %>%
    select(-c_rank)
  res$correlations <- res$correlations[!duplicated(res$correlations$id),]
  return(res$correlations)
}))
qsave(cc1_case, "results/a20/covarying-abundance/cc_level1_case.qs")

cc1_case %>% filter(epi, p1 >= 0.01, p2 >= 0.01) %>% arrange(p.value) %>%
  head(10)

cc2 <- covarying_composition(
  log2cpm           = log2cpm[ix_omni,],
  cluster           = obs$leiden,
  donor             = obs$donor,
  min_donors        = 10,
  min_cells         = 25,
  source_genes      = omni$source_genesymbol,
  target_genes      = omni$target_genesymbol,
  ensembl_to_symbol = ensembl_to_symbol
)
cc2$correlations %<>% mutate(
  epi = substr(c1, 1, 1) == "E" | substr(c2, 1, 1) == "E",
  inter = substr(c1, 1, 1) != substr(c2, 1, 1)
)
cc2$percents <- left_join(cc2$percents, obs %>% select(donor, case) %>% unique, by = "donor")
qsave(cc2, "results/a20/covarying-abundance/cc_level2.qs")

cc2$correlations %>% group_by(c1, c2) %>% summarize(n_signif = sum(fdr < 0.05)) %>% arrange(-n_signif)

cc2$correlations %>% filter(epi, inter, p1 >= 0.03, p2 >= 0.03) %>% head

cc2$correlations %>% filter(g1 == "ITGB1", g2 == "LAMB1", inter)

cc1$correlations %>% filter(g1 == "ITGB1", g2 == "LAMB1", inter) %>% select(-m1, -m2)

fwrite(
  x = cc2$correlations %>% mutate_if(is.numeric, signif),
  file = "results/a20/covarying-abundance/ccc-spearman-level2.tsv",
  sep = "\t"
)

cc2 <- qread("results/a20/covarying-abundance/cc_level2.qs")
for (my_cell in c("CT5", "CT3")) {
  f1 <- function(g1, g2, c1, c2) {
    plot_covarying_genes(
      cc2$percents, g1, g2, c1, c2,
      value.var = "mean", title = FALSE,
      out_dir = "results/a20/communication/cc"
    )
  }
  f1("ITGB2", "ICAM1", my_cell, "E2")
  f1("ITGAL", "ICAM1", my_cell, "E2")
  f1("IL17A", "IL17RA", my_cell, "E2")
  f1("IL17A", "IL17RC", my_cell, "E2")
  f1("IL26", "IL20RA", my_cell, "E2")
  f1("IL26", "IL10RB", my_cell, "E2")
  f1("CXCL13", "CXCR5", my_cell, "T6")
}


my_qs_files <- c(
"results/a20/a12_4_4_b5_1_3/figures/de-case-vs-control/pseudobulk_donor_cluster.qs",
"results/a20/a12_4_4_m3_2/figures/de-case-vs-control/pseudobulk_donor_cluster.qs",
"results/a20/a12_4_4_t4_cd4_2_2/figures/de-case-vs-control/pseudobulk_donor_cluster.qs",
"results/a20/a12_4_4_t4_cd8_1_2/figures/de-case-vs-control/pseudobulk_donor_cluster.qs",
"results/a20/n3_2/figures/de-case-vs-control/pseudobulk_donor_cluster.qs"
)
res <- lapply(my_qs_files, function(qs_file) {
  res <- qread(qs_file)
  my_slug <- qs_file %>% dirname %>% dirname %>% dirname %>% basename
  my_celltype <- switch(my_slug, "a12_4_4_b5_1_3" = "B", "a12_4_4_m3_2" = "M",
    "a12_4_4_t4_cd4_2_2" = "CD4T", "a12_4_4_t4_cd8_1_2" = "CD8T", "n3_2" = "E")
  res$meta$analysis <- my_slug
  res$meta$celltype <- my_celltype
  colnames(res$log2cpm) <- sprintf("%s:%s",
    my_celltype,
    str_remove_all(colnames(res$log2cpm), "factor\\([^)]+\\)")
  )
  x <- summary(res$log2cpm)
  x[[1]] <- rownames(res$log2cpm)[x[[1]]]
  x[[2]] <- colnames(res$log2cpm)[x[[2]]]
  res$x <- x
  res
})
names(res) <- my_qs_files %>% dirname %>% dirname %>% dirname %>% basename

pb <- list(
  meta = rbindlist(lapply(res, "[[", "meta")),
  x = rbindlist(lapply(res, "[[", "x"))
)
pb$x$symbol <- ensembl_to_symbol[as.character(pb$x$i)]
colnames(pb$x) <- c("ensembl_id", "analysis_cluster_donor", "log2cpm", "symbol")
x <- str_split_fixed(pb$x$analysis_cluster_donor, ":", 3)
pb$x$celltype <- x[,1]
pb$x$cluster <- x[,2]
pb$x$donor <- x[,3]
rm(x)
pb$x$mean <- pb$x$log2cpm
pb$x$cluster <- sprintf("%s-%s", pb$x$celltype, pb$x$cluster)
pb$meta$cluster <- sprintf("%s-%s", pb$meta$celltype, pb$meta$cluster)
pb$x <- pb$x[pb$meta[,c("cluster", "donor", "case")], on = .(cluster, donor)]

pb$x %>% filter(symbol == "IL26", cluster == "CD8T-5")
pb$x %>% filter(symbol == "IL10RB", cluster == "E-2")


gg2 <- qread(gg2_file)

# Figure 7
for (my_cell in c("CD8T-5", "CD8T-3")) {
  f1 <- function(g1, g2, c1, c2) {
    my_d <- left_join(pb$meta %>% filter(cluster %in% c(c1, c2)),
      pb$x[cluster == c1][symbol == g1][,c("symbol","cluster","donor","log2cpm")],
      by = c("cluster", "donor")
    )
    my_d <- left_join(my_d,
      pb$x[cluster == c2][symbol == g2][,c("symbol","cluster","donor","log2cpm")],
      by = c("cluster", "donor")
    )
    my_d$log2cpm.x[is.na(my_d$log2cpm.x)] <- 0
    my_d$log2cpm.y[is.na(my_d$log2cpm.y)] <- 0
    my_d$symbol <- ifelse(my_d$cluster == c1, g1, g2)
    my_d$log2cpm <- my_d$log2cpm.x + my_d$log2cpm.y
    my_d <- my_d %>% select(symbol, log2cpm, cluster, donor, case)
    my_d <- left_join(
      my_d %>% filter(symbol == g1),
      my_d %>% filter(symbol == g2),
      by = c("donor", "case")
    )
    my_d$log2cpm.x[is.na(my_d$log2cpm.x)] <- 0
    my_d$log2cpm.y[is.na(my_d$log2cpm.y)] <- 0
    x_lims <- c(0, max(my_d$log2cpm.x, na.rm = TRUE))
    y_lims <- c(0, max(my_d$log2cpm.y, na.rm = TRUE))
    x_lims[!is.finite(x_lims)] <- 0
    y_lims[!is.finite(y_lims)] <- 0
    # top_lim <- max(c(x_lims[2], y_lims[2]))
    # x_lims[2] <- top_lim
    # y_lims[2] <- top_lim
    if (x_lims[2] == 0) {
      x_lims[2] <- y_lims[2]
    }
    ##
    #my_g1 <- g2
    #my_g2 <- g1
    #my_c1 <- str_replace(str_replace(str_remove(c2, "-"), "CD8T", "CT"), "CD4T", "T")
    #my_c2 <- str_replace(str_replace(str_remove(c1, "-"), "CD8T", "CT"), "CD4T", "T")
    #message(glue("{my_c1} {my_g1} | {my_c2} {my_g2}"))
    #my_p <- gg2$limma %>% ungroup() %>%
    #  filter(g2 == my_g2, g1 == my_g1, c1 == my_c1, c2 == my_c2) %>%
    #  select(P.Value, fdr, adj.P.Val, c1, g1, c2, g2, logFC)
    #my_sub <- "P"
    ##
    #my_signif <- function(x, ...) {
    #  str_replace(str_replace(as.character(signif(x, ...)), "0.0$", "0"), "e-0", "e-")
    #}
    #if (nrow(my_p)) {
    #  # my_sub <- with(my_p, glue("P={my_signif(P.Value, 1)}, P<sub>adj</sub>={my_signif(adj.P.Val, 1)}"))
    #  my_sub <- with(my_p, glue("P = {my_signif(P.Value, 1)}"))
    #} else {
    #  my_d <- my_d %>% mutate(resp = (log2cpm.x + log2cpm.y) / 2)
    #  my_p <- infer::t_test(my_d, resp ~ case, order = c("Case", "Control"))
    #  my_sub <- with(my_p, glue("P = {my_signif(p_value, 1)}"))
    #}
      my_d <- my_d %>% mutate(resp = (log2cpm.x + log2cpm.y) / 2)
      my_sub <- "P"
      try({
      my_p <- infer::t_test(my_d, resp ~ case, order = c("Case", "Control"))
      my_sub <- with(my_p, glue("P = {my_signif(p_value, 1)}"))
      })
    #
    p <- ggplot(my_d) +
      geom_abline(slope = 1, intercept = 0, size = 0.3, alpha = 0.3) +
      # geom_segment(
      #   mapping = aes(
      #     x = log2cpm.x,
      #     y = log2cpm.y,
      #     xend = (log2cpm.x + log2cpm.y) / 2,
      #     yend = (log2cpm.x + log2cpm.y) / 2,
      #     group = "donor",
      #     color = append_n(case)
      #   )
      # ) +
      geom_point(
        mapping = aes(x = log2cpm.x, y = log2cpm.y, fill = append_n(case)),
        size = 5, shape = 21, stroke = 0.3
      ) +
      labs(
        x = glue("{c1} <i>{g1}</i>"),
        y = glue("{c2} <i>{g2}</i>"),
        subtitle = my_sub
      ) +
      scale_fill_manual(name = NULL, values = c("#E69F00", "#333333")) +
      scale_color_manual(name = NULL, values = c("#E69F00", "#333333")) +
      scale_x_continuous(breaks = x_lims, limits = x_lims, labels = my_signif(x_lims, 2)) +
      scale_y_continuous(breaks = y_lims, limits = y_lims, labels = my_signif(y_lims, 2)) +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      theme_kamil +
      theme(
        legend.position = "none",
        # axis.title.x = element_text(face = "italic"),
        # axis.title.y = element_text(face = "italic")
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        plot.subtitle = ggtext::element_markdown()
      )
    my_filename <- glue("{out_dir}/cc_log2cpm_{g1}_{g2}_{c1}_{c2}.pdf")
    message(my_filename)
    ggsave(
      filename = my_filename,
      plot = p,
      scale = scale, width = 4, height = 3, units = "in", dpi = 300
    )
    return(p)
  }
  f1("IL26", "IL20RA", my_cell, "E-2")
  f1("ITGB2", "ICAM1", my_cell, "E-2")
  f1("ITGAL", "ICAM1", my_cell, "E-2")
  f1("IL17A", "IL17RA", my_cell, "E-2")
  f1("IL17A", "IL17RC", my_cell, "E-2")
  f1("IL26", "IL10RB", my_cell, "E-2")
  f1("CXCL13", "CXCR5", my_cell, "CD4T-6")
  f1("HAVCR2", "CEACAM1", my_cell, "E-2")
}
#
p1 <- f1("ITGB2", "ICAM1", "CD8T-3", "E-2")
p2 <- f1("ITGAL", "ICAM1", "CD8T-3", "E-2")
p3 <- f1("ITGB2", "ICAM1", "CD8T-5", "E-2") + theme(axis.title.y = element_blank())
p4 <- f1("ITGAL", "ICAM1", "CD8T-5", "E-2") + theme(axis.title.y = element_blank())
p <- p1 + p2 + p3 + p4 + plot_layout(
  design = "ab
cd"
)
# + plot_annotation(title = "E-2 *ICAM1*", theme = theme(plot.title = element_markdown()))
my_filename <- glue("{out_dir}/cc_log2cpm_ICAM1.pdf")
message(my_filename)
ggsave(
  filename = my_filename,
  plot = p,
  scale = 0.8, width = 6, height = 6, units = "in", dpi = 300
)
#
p1 <- f1("IL26", "IL20RA", "CD8T-3", "E-2")
p2 <- f1("IL26", "IL10RB", "CD8T-3", "E-2")
p3 <- f1("IL26", "IL20RA", "CD8T-5", "E-2") + theme(axis.title.y = element_blank())
p4 <- f1("IL26", "IL10RB", "CD8T-5", "E-2") + theme(axis.title.y = element_blank())
p <- p1 + p2 + p3 + p4 + plot_layout(
  design = "ab
cd"
)
my_filename <- glue("{out_dir}/cc_log2cpm_IL26.pdf")
message(my_filename)
ggsave(
  filename = my_filename,
  plot = p,
  scale = 0.8, width = 6, height = 6, units = "in", dpi = 300
)
#
p1 <- f1("IL17A", "IL17RA", "CD8T-3", "E-2")
p2 <- f1("IL17A", "IL17RC", "CD8T-3", "E-2")
p3 <- f1("IL17A", "IL17RA", "CD8T-5", "E-2") + theme(axis.title.y = element_blank())
p4 <- f1("IL17A", "IL17RC", "CD8T-5", "E-2") + theme(axis.title.y = element_blank())
p <- p1 + p2 + p3 + p4 + plot_layout(
  design = "ab
cd"
)
my_filename <- glue("{out_dir}/cc_log2cpm_IL17A.pdf")
message(my_filename)
ggsave(
  filename = my_filename,
  plot = p,
  scale = 0.8, width = 6, height = 6, units = "in", dpi = 300
)
#
p1 <- f1("CXCL13", "CXCR5", "CD8T-3", "CD4T-6")
p2 <- f1("CXCL13", "CXCR5", "CD8T-5", "CD4T-6")
p <- p1 + p2 + plot_layout(
  design = "a
b"
)
my_filename <- glue("{out_dir}/cc_log2cpm_CXCL13.pdf")
message(my_filename)
ggsave(
  filename = my_filename,
  plot = p,
  scale = 0.8, width = 3, height = 6, units = "in", dpi = 300
)
#
p1 <- f1("HAVCR2", "CEACAM1", "CD8T-3", "E-2")
p2 <- f1("HAVCR2", "CEACAM1", "CD8T-5", "E-2")
p <- p1 + p2 + plot_layout(
  design = "a
b"
)
my_filename <- glue("{out_dir}/cc_log2cpm_CEACAM1.pdf")
message(my_filename)
ggsave(
  filename = my_filename,
  plot = p,
  scale = 0.8, width = 3, height = 6, units = "in", dpi = 300
)

cc2_case <- rbindlist(lapply(c("Case", "Control"), function(this_case) {
  ix_case <- obs$case == this_case
  res <- covarying_composition(
    log2cpm           = log2cpm[ix_omni,ix_case],
    cluster           = obs$leiden[ix_case],
    donor             = obs$donor[ix_case],
    min_donors        = 10,
    min_cells         = 25,
    source_genes      = omni$source_genesymbol,
    target_genes      = omni$target_genesymbol,
    ensembl_to_symbol = ensembl_to_symbol
  )
  res$correlations$case <- this_case
  res$correlations %<>% mutate(
    epi = substr(c1, 1, 1) == "E" | substr(c2, 1, 1) == "E",
    inter = substr(c1, 1, 1) != substr(c2, 1, 1)
  )
  return(res$correlations)
}))
qsave(cc2_case, "results/a20/covarying-abundance/cc_level2_case.qs")

# cd1 <- matrix_to_table(
#   log2cpm = log2cpm[ix_omni,],
#   cluster = obs$dataset,
#   donor = obs$donor,
#   min_donors = 10,
#   min_cells = 25,
#   ensembl_to_symbol = ensembl_to_symbol
# )
# cd1 <- left_join(cd1, obs %>% select(donor, case) %>% unique, by = "donor")

x <- cc1_case %>%
  filter(epi, p1 >= 0.01, p2 >= 0.01) %>%
  select(c1, g1, c2, g2, p.value, case) %>%
  pivot_wider(names_from = "case", values_from = "p.value", values_fill = 1)
  # select(-m1, -m2, -epi, -p1, -p2)
x %>% arrange(log10(Case) / log10(Control)) %>% head(10)
p <- ggplot(x %>% mutate(Case = -log10(Case), Control = -log10(Control))) +
  aes(x = Control, y = Case) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = 2) +
  geom_vline(xintercept = 0, size = 0.2, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.2, linetype = 2) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "1", glue("10<sup>-{x}</sup>"))) +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "1", glue("10<sup>-{x}</sup>"))) +
  labs(
    x = "Case", y = "Control",
    title = "Spearman p-values"
  ) +
  theme(
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown()
  )
my_ggsave(
  "spearman_case_control",
  out_dir = "results/a20/covarying-abundance/level1",
  types = "pdf",
  plot = p,
  scale = 1, width = 4, height = 3.5, units = "in", dpi = 300
)

# plot_covarying_genes(
#   d = cd, g1 = "CCR5", c1 = "M", g2 = "IL16", c2 = "CT",
#   out_dir = "results/a20/covarying-abundance/cc_major/separate"
# )

# plot_communication2(
#   percents = cc1$percents,
#   # correlations = cc1_case %>% filter(c1 == "M", g1 == "CX3CR1", c2 == "E", g2 == "CXCL2"),
#   correlations = cc1_case %>% filter(c1 == "E", g1 == "ICAM1", c2 == "CT", g2 == "IL10"),
#   out_dir = "results/a20/covarying-abundance/level1/spearman/case"
# )

plot_communication2(
  percents = cc1$percents,
  correlations = cc1_case %>% filter(g1 == "PDCD1", p1 > 0.01, p2 > 0.01),
  out_dir = "results/a20/covarying-abundance/level1/spearman/case"
)

plot_communication2(
  percents = cc1$percents,
  correlations = cc1_case %>%
    filter(p1 >= 0.01, p2 >= 0.01, inter) %>%
    group_by(case, c1, c2) %>%
    mutate(fdr = p.adjust(p.value, method = "BH")) %>%
    arrange(p.value) %>%
    filter(fdr < 0.2) %>%
    select(-m1, -m2, -donors1, -donors2) %>% head(20),
  out_dir = "results/a20/covarying-abundance/level1/spearman/case"
)

plot_communication2(
  percents = cc1$percents,
  correlations = cc1_case %>%
    filter(p1 >= 0.01, p2 >= 0.01, inter, epi) %>%
    group_by(case, c1, c2) %>%
    mutate(fdr = p.adjust(p.value, method = "BH")) %>%
    arrange(p.value) %>%
    filter(fdr < 0.2) %>%
    select(-m1, -m2, -donors1, -donors2) %>% head(20),
  out_dir = "results/a20/covarying-abundance/level1/spearman/case"
)

plot_communication2(
  percents = cc1$percents,
  correlations = cc1_case %>%
    filter(str_detect(g1, "ITG") | str_detect(g2, "ITG")) %>%
    arrange(p.value) %>% head(20),
  out_dir = "results/a20/covarying-abundance/level1/spearman/case"
)
plot_communication2(
  percents = cc1$percents,
  correlations = cc1_case %>%
    filter(str_detect(g1, "ITG") | str_detect(g2, "ITG"), epi) %>%
    arrange(p.value) %>% head(20),
  out_dir = "results/a20/covarying-abundance/level1/spearman/case"
)

# x <- inner_join(
#   x = cc1$percents %>% filter(symbol == "ICAM1", cluster == "E"),
#   y = cc1$percents %>% filter(symbol == "IFNG", cluster == "CT"),
#   by = c("donor", "case")
# )
# cor.test(log1p(x$percent.x), log1p(x$percent.y), method = "spearman", exact = FALSE)

plot_communication(
  percents = cc1$percents,
  correlations = cc1$correlations  %>%
      filter(p1 >= 0.01, p2 >= 0.01, inter, epi) %>%
      group_by(c1, c2) %>%
      mutate(fdr = p.adjust(p.value, method = "BH")) %>%
      arrange(p.value) %>%
      filter(fdr < 0.05) %>%
      select(-m1, -m2, -donors1, -donors2) %>% head(20),
  log_percent = TRUE,
  out_dir = "results/a20/covarying-abundance/level1/spearman"
)

plot_communication(
  percents = cc1$percents,
  correlations = cc1$correlations  %>%
      filter(p1 >= 0.01, p2 >= 0.01, inter) %>%
      group_by(c1, c2) %>%
      mutate(fdr = p.adjust(p.value, method = "BH")) %>%
      arrange(p.value) %>%
      filter(fdr < 0.05) %>%
      select(-m1, -m2, -donors1, -donors2) %>% head(20),
  log_percent = TRUE,
  out_dir = "results/a20/covarying-abundance/level1/spearman"
)

plot_communication(
  percents = cc2$percents,
  correlations = cc2$correlations  %>%
      filter(p1 >= 0.01, p2 >= 0.01, inter, epi) %>%
      group_by(c1, c2) %>%
      mutate(fdr = p.adjust(p.value, method = "BH")) %>%
      arrange(p.value) %>%
      filter(fdr < 0.05) %>%
      select(-m1, -m2, -donors1, -donors2) %>% head(20),
  log_percent = TRUE,
  out_dir = "results/a20/covarying-abundance/level2/spearman"
)

plot_communication(
  percents = cc2$percents,
  correlations = cc2$correlations  %>%
      filter(p1 >= 0.01, p2 >= 0.01, inter) %>%
      group_by(c1, c2) %>%
      mutate(fdr = p.adjust(p.value, method = "BH")) %>%
      arrange(p.value) %>%
      filter(fdr < 0.05) %>%
      select(-m1, -m2, -donors1, -donors2) %>% head(20),
  log_percent = TRUE,
  out_dir = "results/a20/covarying-abundance/level2/spearman"
)

# omni %>% filter(source_genesymbol == "IL26")

omni %>% filter(target_genesymbol == "IL15", source_genesymbol == "IFNG")

# res <- entrez_search(
#   db   = "pubmed",
#   term = "IL15RA AND (colitis[MeSH])"
# )

# res_case %>% filter(c1 == "T", g1 == "ITGB2", c2 == "B", g2 == "ICAM1")

# plot_communication2(
#   res_case %>%
#   filter(c1 == "E", g1 == "ICAM1", c2 == "CT", g2 == "IL10"),
#   out_dir = "results/a20/covarying-abundance/level1_case"
# )

# plot_communication2(
#   res_case %>% filter(g1 == "FCGR1A"),
#   out_dir = "results/a20/covarying-abundance/level1_case"
# )



# UMAP - turns out it is not very clear, better to use the scatter plots
#
# cc1_case %>% filter(g1 == "ITGB2", g2 == "ICAM1", c1 == "T", c2 == "B") 
# my_gs <- c("ITGB2", "ICAM1")
# my_cs <- c("T", "B")
# x <- inner_join(
#   x = cc1$percents %>% filter(symbol == my_gs[1], cluster == my_cs[1]),
#   y = cc1$percents %>% filter(symbol == my_gs[2], cluster == my_cs[2]),
#   by = c("donor", "case")
# )
# donor_order <- x$donor[order(x$percent.x + x$percent.y)]
# # cor.test(log1p(x$percent.x), log1p(x$percent.y), method = "spearman", exact = FALSE)
# for (i in 1:2) {
#   for (my_case in c("Case", "Control")) {
#     my_ix <- which(obs$dataset %in% c(my_cs[i]) & obs$donor %in% donor_order & obs$case == my_case)
#     my_id <- names(which(ensembl_to_symbol == my_gs[i]))
#     n_donors <- length(unique(obs$donor[my_ix]))
#     obs$my_donor <- factor(obs$donor, donor_order)
#     obs$my_donor <- factor(pub_ids[as.character(obs$my_donor)], pub_ids[levels(obs$my_donor)])
#     p <- plot_hexgene(
#       x = obs$UMAP1[my_ix],
#       y = obs$UMAP2[my_ix],
#       z = as.numeric(log2cpm[my_id,my_ix]),
#       group = obs$my_donor[my_ix],
#       text = FALSE,
#       bins = 15,
#       palette = "oslo"
#     ) +
#     facet_wrap(~ group, ncol = n_donors) +
#     labs(title = glue("<i>{ensembl_to_symbol[my_id]}</i>+ {my_cs[i]}")) +
#     theme(
#       panel.spacing = unit(0.5, "lines"),
#       plot.title = ggtext::element_markdown(face = "plain"),
#       legend.position = "none"
#     )
#     my_ggsave(
#       glue("umap-{safe(ensembl_to_symbol[my_id])}-{my_cs[i]}-{my_case}"),
#       out_dir = file.path("results/a20/covarying-abundance/level1/umap"),
#       type = "pdf",
#       plot = p,
#       scale = 0.7,
#       width = 1 + n_donors * 2,
#       height = 3,
#       units = "in",
#       dpi = 300
#     )
#   }
# }


d1 <- cc1$correlations %>%
  filter(p1 > 0.01, p2 > 0.01) %>%
  mutate(
    fdr = p.adjust(p.value, method = "BH"),
    gene_pair = glue("{g2} - {g1}"),
    cell_pair = glue("{c2} - {c1}"),
    signed_p = -log10(p.value) * sign(estimate)
  )
my_gene_pairs <- unique((d1 %>% filter(fdr < 0.01, epi))$gene_pair)
length(my_gene_pairs)
d1 <- d1 %>% filter(gene_pair %in% my_gene_pairs)
# gene_pairs <- (
#   d1 %>% group_by(gene_pair) %>% summarize(min_pval = min(p.value)) %>% arrange(-min_pval)
# )$gene_pair
# d1$gene_pair <- factor(as.character(d1$gene_pair), levels = gene_pairs)
set.seed(1)
sorted_gene_pairs <- umap_sorted_rows(d1, "gene_pair ~ cell_pair", "signed_p")
d1$gene_pair <- factor(as.character(d1$gene_pair), rev(sorted_gene_pairs))
d1$c2 <- factor(as.character(d1$c2), c("E", "CT", "T", "M", "B"))
#
col_fun <- circlize::colorRamp2(
  seq(-max(abs(d1$signed_p)), max(abs(d1$signed_p)), length.out = 11),
  rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
)
# color_values <- scales::rescale(seq(min(d1$signed_p), max(d1$signed_p), length.out = 11))
color_values <- seq(min(d1$signed_p), max(d1$signed_p), length.out = 11)
p <- ggplot() +
  geom_tile(
    data = d1,
    mapping = aes(x = cell_pair, y = gene_pair, fill = signed_p)
  ) +
  geom_point(
    data = d1 %>% filter(fdr < 0.05),
    mapping = aes(x = cell_pair, y = gene_pair),
    shape = 19, color = "white"
  ) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_x_discrete(
    position = "t", name = NULL, expand = c(0, 0)
  ) +
  facet_wrap(~ c2, ncol = 5, scales = "free_x") +
  scale_fill_gradientn(
    colors = col_fun(color_values),
    name = "Signed -log10 P",
    guide = guide_colorbar(barwidth = 10)
    # breaks = -8:8,
    # labels = function(x) fractional::fractional(10^x)
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 0),
    axis.text.y = element_text(face = "italic"),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )
my_ggsave(
  "spearman_heatmap-top",
  out_dir = file.path("results/a20/covarying-abundance/level1"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = length(unique(d1$cell_pair)) * 0.4 + 1,
  height = length(unique(d1$gene_pair)) * 0.3 + 1.1,
  units = "in",
  dpi = 300
)

plot_communication(
  percents = cc1$percents,
  correlations = d1 %>% filter(p1 > 0.01, p2 > 0.01, fdr < 0.05),
  log_percent = TRUE,
  out_dir = "results/a20/covarying-abundance/level1/spearman"
)


d1 <- cc1_case %>%
  filter(p1 >= 0.01, p2 >= 0.01, inter) %>%
  group_by(case, c1, c2) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value) %>%
  mutate(
    gene_pair = glue("{g1} - {g2}"),
    cell_pair = glue("{c1} - {c2}"),
    signed_p = -log10(p.value) * sign(estimate)
  ) %>%
  ungroup()
my_gene_pairs <- unique((d1 %>% filter(case == "Case", fdr < 0.20, epi))$gene_pair)
# my_gene_pairs <- unique(d1$gene_pair)
length(my_gene_pairs)
d1 <- d1 %>% filter(gene_pair %in% my_gene_pairs)
#
col_fun <- circlize::colorRamp2(
  seq(-max(abs(d1$signed_p)), max(abs(d1$signed_p)), length.out = 11),
  rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
)
# color_values <- scales::rescale(seq(min(d1$signed_p), max(d1$signed_p), length.out = 11))
color_values <- seq(min(d1$signed_p), max(d1$signed_p), length.out = 11)
p <- ggplot() +
  geom_tile(
    data = d1,
    mapping = aes(x = cell_pair, y = gene_pair, fill = signed_p)
  ) +
  geom_point(
    data = d1 %>% filter(fdr < 0.05),
    mapping = aes(x = cell_pair, y = gene_pair),
    shape = 19, color = "white"
  ) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_x_discrete(
    position = "b", name = NULL, expand = c(0, 0)
  ) +
  # facet_grid(c1 ~ case, scales = "free") +
  facet_wrap(case ~ c1, ncol = 5, scales = "free_x") +
  scale_fill_gradientn(
    colors = col_fun(color_values),
    name = "Signed -log10 P",
    guide = guide_colorbar(barwidth = 10)
    # breaks = -8:8,
    # labels = function(x) fractional::fractional(10^x)
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    # strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )
my_ggsave(
  "heatmap-top-case-epi",
  out_dir = file.path("results/a20/covarying-abundance/level1"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = length(unique(d1$cell_pair)) * 0.4 + 1,
  # height = length(unique(d1$gene_pair)) * 0.3 + 2.1,
  height = 2 * (length(unique(d1$gene_pair)) * 0.3 + 2.1),
  units = "in",
  dpi = 300
)

my_pairs <- readxl::read_excel("data/ligand-receptor.xlsx")
my_genes <- unique(c(my_pairs$g_target, my_pairs$g_source))
d1 <- cc1_case %>%
  filter(p1 >= 0.01, p2 >= 0.01, inter) %>%
  group_by(case, c1, c2) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value) %>%
  mutate(
    gene_pair = glue("{g1} - {g2}"),
    cell_pair = glue("{c1} - {c2}"),
    signed_p = -log10(p.value) * sign(estimate)
  ) %>%
  ungroup()
d1 <- d1 %>% filter(g1 %in% my_genes | g2 %in% my_genes, epi, fdr < 0.5)
col_fun <- circlize::colorRamp2(
  seq(-max(abs(d1$signed_p)), max(abs(d1$signed_p)), length.out = 11),
  rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
)
# color_values <- scales::rescale(seq(min(d1$signed_p), max(d1$signed_p), length.out = 11))
color_values <- seq(min(d1$signed_p), max(d1$signed_p), length.out = 11)
p <- ggplot() +
  geom_tile(
    data = d1,
    mapping = aes(x = cell_pair, y = gene_pair, fill = signed_p)
  ) +
  geom_point(
    data = d1 %>% filter(fdr < 0.05),
    mapping = aes(x = cell_pair, y = gene_pair),
    shape = 19, color = "white"
  ) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_x_discrete(
    position = "b", name = NULL, expand = c(0, 0)
  ) +
  # facet_grid(c1 ~ case, scales = "free") +
  facet_wrap(case ~ c1, ncol = 5, scales = "free_x") +
  scale_fill_gradientn(
    colors = col_fun(color_values),
    name = "Signed -log10 P",
    guide = guide_colorbar(barwidth = 10)
    # breaks = -8:8,
    # labels = function(x) fractional::fractional(10^x)
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    # strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )
my_ggsave(
  "heatmap-selected-case-epi",
  out_dir = file.path("results/a20/covarying-abundance/level1"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = length(unique(d1$cell_pair)) * 0.4 + 2.5,
  # height = length(unique(d1$gene_pair)) * 0.3 + 2.1,
  height = 2 * (length(unique(d1$gene_pair)) * 0.3 + 2.1),
  units = "in",
  dpi = 300
)



d1 <- cc1_case %>%
  filter(p1 >= 0.01, p2 >= 0.01, inter) %>%
  group_by(case, c1, c2) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value) %>%
  mutate(
    gene_pair = sprintf("%s - %s", g1, g2),
    cell_pair = sprintf("%s - %s", c1, c2),
    signed_p = -log10(p.value) * sign(estimate)
  ) %>%
  ungroup()
my_gene_pairs <- unique((d1 %>% filter(case == "Case", fdr < 0.1))$gene_pair)
# my_gene_pairs <- unique(d1$gene_pair)
length(my_gene_pairs)
d1 <- d1 %>% filter(gene_pair %in% my_gene_pairs) %>%
  mutate(cell_pair_case = sprintf("%s %s", cell_pair, case))
  # mutate(c1 = sprintf("%s %s", c1, case))
#
col_fun <- circlize::colorRamp2(
  seq(-max(abs(d1$signed_p)), max(abs(d1$signed_p)), length.out = 11),
  rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
)
# color_values <- scales::rescale(seq(min(d1$signed_p), max(d1$signed_p), length.out = 11))
color_values <- seq(min(d1$signed_p), max(d1$signed_p), length.out = 11)
p <- ggplot() +
  geom_tile(
    data = d1,
    mapping = aes(x = cell_pair_case, y = gene_pair, fill = signed_p)
  ) +
  geom_point(
    data = d1 %>% filter(fdr < 0.05),
    mapping = aes(x = cell_pair_case, y = gene_pair),
    shape = 19, color = "white"
  ) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_x_discrete(
    position = "t", name = NULL, expand = c(0, 0)
  ) +
  # facet_grid(c1 ~ case, scales = "free") +
  facet_grid(~ cell_pair, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colors = col_fun(color_values),
    name = "Signed -log10 P",
    guide = guide_colorbar(barwidth = 10)
    # breaks = -8:8,
    # labels = function(x) fractional::fractional(10^x)
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 0),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )
my_ggsave(
  "heatmap-top-case",
  out_dir = file.path("results/a20/covarying-abundance/level1"),
  type = "pdf",
  plot = p,
  scale = 1,
  width = length(unique(d1$cell_pair)) * 0.8 + 1,
  height = length(unique(d1$gene_pair)) * 0.3 + 2.1,
  units = "in",
  dpi = 300
)




















# Results at the level of major cell type
########################################################################
obs$cluster <- obs$leiden
d <- as.data.table(summary(log2cpm[omni_ens,]))
d$symbol <- omni_sym[d$i]
d$cluster <- str_remove(obs$cluster[d$j], "\\d+$")
d$donor <- obs$donor[d$j]
y <- obs %>% group_by(cluster = str_remove(cluster, "\\d+$"), donor) %>% count() %>% as.data.table
d <- y[d, on = c("cluster", "donor")]
d <- d[
  ,
  .(
    mean = mean(c(x, rep(0, n[1] - length(x)))),
    percent = sum(x > 0) / n[1]
  ),
  by = .(cluster, donor, symbol)
]
setkey(d, symbol, cluster, donor)
# d$percent_logit <- car::logit(d$percent)
# The gene is robustly detected in some cluster for 20 donors.
min_donors <- 10
d_nonzero <- d[,.(nonzero = sum(mean > 0)),by = .(cluster, symbol)]
robust_genes <- d_nonzero[nonzero >= min_donors]$symbol %>% unique
intersect(my_genes, robust_genes)
#
d <- d[symbol %in% robust_genes]
# Add zeros when we don't have anything.
d <- tidyr::complete(
  d, cluster, donor, symbol, fill = list(mean = 0, percent = 0)
) %>% as.data.table
x <- obs %>% select(donor, case) %>% unique
donor_to_case <- x$case
names(donor_to_case) <- x$donor
d$case <- donor_to_case[d$donor]
#
# ix <- ram$ligand_approved_symbol %in% robust_genes &
#       ram$receptor_approved_symbol %in% robust_genes
# ram <- ram[ix,]
ix <- omni$source_genesymbol %in% robust_genes &
  omni$target_genesymbol %in% robust_genes
omni <- omni[ix,]
# Add the number of cells each donor has in each cluster
d <- left_join(d, obs %>% count(donor, cluster = str_remove(cluster, "\\d+$")), by = c("cluster", "donor"))
d <- as.data.table(d)
d <- d[!is.na(d$n)]
d <- d[d$n >= 25]

d2 <- as.data.table(summary(log2cpm[omni_ens,]))
d2$symbol <- omni_sym[d2$i]
d2$cluster <- str_remove(obs$cluster[d2$j], "\\d+$")
y <- obs %>% group_by(cluster = str_remove(cluster, "\\d+$")) %>% count() %>% as.data.table
d2 <- y[d2, on = c("cluster")]
d2 <- d2[
  ,
  .(
    mean = mean(c(x, rep(0, n[1] - length(x)))),
    percent = sum(x > 0) / n[1]
  ),
  by = .(cluster, symbol)
]

# Exploring Simpson paradox
########################################################################
# by_case <- d %>%
#   group_by(cluster, symbol, case) %>%
#   summarize(percent = 100 * sum(percent * n) / sum(n), .groups = "drop") %>%
#   pivot_wider(names_from = "case", values_from = "percent") %>%
#   mutate(diff = Case - Control)
# by_donor <- d %>%
#   mutate(percent = 100 * percent) %>%
#   group_by(cluster, symbol) %>%
#   summarize(
#     pval = t.test(log1p(percent) ~ case)$p.value,
#     diffmeans = mean(percent[case == "Case"]) - mean(percent[case == "Control"]),
#     .groups = "drop"
#   )
# by_case <- left_join(by_case, by_donor, by = c("cluster", "symbol"))
# by_case %<>% arrange(-abs(diff))

# by_case  %>% arrange(pval)

# by_case %>% filter(pval > 0.05)

# d %>% filter(cluster == "T", symbol == "TNFRSF18")

# my_pairs <- list(
#   c(cluster = "B", symbol = "TNFRSF17"),
#   c(cluster = "B", symbol = "IL21R"),
#   c(cluster = "B", symbol = "EZR"),
#   c(cluster = "B", symbol = "LTB"),
#   c(cluster = "T", symbol = "IL7"),
#   c(cluster = "T", symbol = "EBI3"),
#   c(cluster = "CT", symbol = "IL17RA")
# )
# for (my_pair in my_pairs) {
#   my_cluster <- my_pair['cluster']
#   my_symbol <- my_pair['symbol']
#   my_pval <- (by_case %>% filter(cluster == my_cluster, symbol == my_symbol))$pval
#   d1 <- d %>% filter(cluster == my_cluster, symbol == my_symbol)
#   d1_mean <- d1 %>% group_by(case) %>% summarize(percent = 100 * mean(percent), .groups = "drop")
#   d1_median <- d1 %>% group_by(case) %>% summarize(percent = 100 * median(percent), .groups = "drop")
#   d2 <- d %>% filter(cluster == my_cluster, symbol == my_symbol) %>% group_by(case) %>%
#     summarize(percent = 100 * sum(percent * n) / sum(n), .groups = "drop")
#   p1 <- ggplot() +
#     geom_errorbar(
#       aes(x = case, ymin = percent, ymax = percent),
#       data = d1_mean, width = 0.8, size = 0.3
#     ) +
#     geom_quasirandom(
#       data = d1,
#       aes(x = case, y = 100 * percent, fill = case),
#       shape = 21, size = 3, stroke = 0.2,
#       groupOnX = TRUE,
#       width = 0.4
#     ) +
#     scale_y_log10() +
#     annotation_logticks(sides = "l", size = 0.3) +
#     annotate(
#       geom = "text",
#       x = -Inf, y = Inf,
#       hjust = -0.05, vjust = 1.2,
#       label = glue("P = {signif(my_pval, 2)}")
#     ) +
#     scale_fill_manual(values = pals::okabe()[2:1]) +
#     theme(legend.position = "none") +
#     labs(
#       x = NULL, y = "Percent", title = "Percent of cells by donor",
#       subtitle = glue(
#         "{signif(d1_mean$percent[1], 2)}% of cells in a Case\n{signif(d1_mean$percent[2], 2)}% of cells in a Control"
#       )
#     )
#   p2 <- ggplot(d2) +
#     aes(x = case, y = percent, fill = case) +
#     geom_col() +
#     scale_fill_manual(values = pals::okabe()[2:1]) +
#     theme(legend.position = "none") +
#     labs(
#       x = NULL, y = "Percent", title = "Percent of cells by group",
#       subtitle = glue(
#         "{signif(d2$percent[1], 2)}% of Case cells\n{signif(d2$percent[2], 2)}% of Control cells"
#       )
#     )
#   my_ggsave(
#     glue("{my_cluster}_{my_symbol}"),
#     out_dir = "results/a20/covarying-abundance/denominator",
#     types = "png",
#     plot = p1 + p2 + plot_annotation(title = glue("{my_cluster} {my_symbol}")),
#     scale = 1, width = 8, height = 4, units = "in", dpi = 300
#   )
# }

# by_case %>% filter(diff < 0, diffmeans > 0) %>% arrange(-pval)

# by_case %>% arrange(-abs(diff - diffmeans))

# p <- ggplot(by_case) +
#   aes(x = diff, y = diffmeans) +
#   geom_abline(yintercept = 0, slope = 1, type = 2, size = 0.2) +
#   geom_vline(xintercept = 0, size = 0.2) +
#   geom_hline(yintercept = 0, size = 0.2) +
#   geom_point(size = 0.5) +
#   labs(
#     x = "Difference by group",
#     y = "Difference by donor"
#   )
# my_ggsave(
#   glue("diff_diffmeans"),
#   out_dir = "results/a20/covarying-abundance/denominator",
#   types = "png",
#   plot = p,
#   scale = 1, width = 5, height = 4, units = "in", dpi = 300
# )

# res <- pblapply(seq(nrow(ram)), function(ram_i) {
#   # g1 <- "ADAM15"
#   # g2 <- "ITGAV"
#   g1 <- ram$ligand_approved_symbol[ram_i]
#   g2 <- ram$receptor_approved_symbol[ram_i]
res <- pblapply(seq(nrow(omni)), function(omni_i) {
  # g1 <- "ADAM15"
  # g2 <- "ITGAV"
  g1 <- omni$source_genesymbol[omni_i]
  g2 <- omni$target_genesymbol[omni_i]
  x <- d[symbol == g1]
  y <- d[symbol == g2]
  stopifnot(nrow(x) == nrow(y))
  xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
  # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
  xx_donor <- xx$donor
  xx <- as.matrix(xx[,2:ncol(xx)])
  rownames(xx) <- xx_donor
  yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
  # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
  yy_donor <- yy$donor
  yy <- as.matrix(yy[,2:ncol(yy)])
  rownames(yy) <- yy_donor
  # Remove zeros
  xx[xx == 0] <- NA
  yy[yy == 0] <- NA
  sp <- Hmisc::rcorr(xx, yy, "spearman")
  sp_r <- sp$r[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
  sp_p <- sp$P[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
  ut <- upper.tri(sp_r, diag = FALSE)
  res <- data.table(
    g1 = g1,
    g2 = g2,
    c1 = rownames(sp_r)[row(sp_r)[ut]],
    c2 = rownames(sp_r)[col(sp_r)[ut]],
    estimate = sp_r[ut],
    p.value = sp_p[ut]
  )
})
res <- rbindlist(res)
res <- res[!is.na(res$estimate)]
res$p.value[res$p.value == 0] <- 1
d_n <- d[, .(donors = sum(mean > 0)), by = .(cluster, symbol)]
d_n$cluster <- as.character(d_n$cluster)
res2 <- d_n[res, on = c("symbol" = "g1", "cluster" = "c1")]
res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "g2")] 
colnames(res3) <- c("c1", "g1", "donors1", "c2", "g2", "donors2", "estimate", "p.value")
res <- res3
rm(res2)
rm(res3)
res %<>% arrange(p.value)
res <- unique(res)
res <- res[donors1 >= min_donors & donors2 >= min_donors,]

top_res <- res %>%
  mutate(inter = substr(c1, 1, 1) != substr(c2, 1, 1)) %>%
  mutate(epi = xor(c1 == "E", c2 == "E")) %>%
  filter(donors1 >= min_donors & donors2 >= min_donors) %>%
  # filter(inter) %>%
  arrange(p.value)
top_res %<>%
  group_by(c1, c2) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  mutate(signif = fdr < 0.05) %>%
  ungroup()
  # mutate(signif = p.value < 0.05 / nrow(top_res))
# top_res %>% filter(signif) %>% count(c1, c2)
top_res <- left_join(
  top_res, d2 %>% select(c1 = cluster, g1 = symbol, p1 = percent, m1 = mean),
  by = c("c1", "g1")
)
top_res <- left_join(
  top_res, d2 %>% select(c2 = cluster, g2 = symbol, p2 = percent, m2 = mean),
  by = c("c2", "g2")
)
top_res

# top_res <- left_join(
#   x = top_res,
#   y = omni_full %>% select(source_genesymbol, target_genesymbol, is_stimulation, is_inhibition),
#   by = c("g1" = "source_genesymbol", "g2" = "target_genesymbol")
# )

qsave(top_res, "results/a20/covarying-abundance/res_major.qs")

plot_communication(
  top_res %>% filter(epi, fdr < 0.05) %>% select(-inter, -epi, -fdr, -signif, -donors1, -donors2) %>%
  filter(p1 >= 0.01 & p2 >= 0.01),
  out_dir = "results/a20/covarying-abundance/cc_major"
)

plot_communication(
  top_res %>% filter(fdr < 0.05) %>% select(-inter, -epi, -fdr, -signif, -donors1, -donors2) %>%
  filter(p1 >= 0.01 & p2 >= 0.01),
  out_dir = "results/a20/covarying-abundance/cc_major"
)

  # top_res %>% filter(fdr < 0.05) %>% select(-inter, -epi, -fdr, -signif, -donors1, -donors2) %>%
  # filter(p1 >= 0.01 & p2 >= 0.01) %>%
  # filter(g1 == "CXCR3")
  # # filter(g1 == "FCGR1A")
  # # filter(c2 == "B")




# Results at the level of cell clusters within each major cell type
########################################################################

# d <- as.data.table(summary(log2cpm[ram_genes$ensembl_id,]))
# d$symbol <- ram_genes$symbol[d$i]
obs$cluster <- obs$leiden
# obs$cluster <- ifelse(str_detect(obs$leiden, "^E"), "epithelial", "immune")
# obs$cluster <- str_remove(obs$leiden, "[0-9]+")
table(obs$cluster)
d <- as.data.table(summary(log2cpm[omni_ens,]))
d$symbol <- omni_sym[d$i]
d$cluster <- obs$cluster[d$j]
d$donor <- obs$donor[d$j]
y <- obs %>% group_by(cluster, donor) %>% count() %>% as.data.table
d <- y[d, on = c("cluster", "donor")]
d <- d[
  ,
  .(
    mean = mean(c(x, rep(0, n[1] - length(x)))),
    percent = sum(x > 0) / n[1]
  ),
  by = .(cluster, donor, symbol)
]
setkey(d, symbol, cluster, donor)
# d$percent_logit <- car::logit(d$percent)
# The gene is robustly detected in some cluster for 20 donors.
min_donors <- 10
d_nonzero <- d[,.(nonzero = sum(mean > 0)),by = .(cluster, symbol)]
robust_genes <- d_nonzero[nonzero >= min_donors]$symbol %>% unique
intersect(my_genes, robust_genes)
#
d <- d[symbol %in% robust_genes]
# Add zeros when we don't have anything.
d <- tidyr::complete(
  d, cluster, donor, symbol, fill = list(mean = 0, percent = 0)
) %>% as.data.table
x <- obs %>% select(donor, case) %>% unique
donor_to_case <- x$case
names(donor_to_case) <- x$donor
d$case <- donor_to_case[d$donor]
#
# ix <- ram$ligand_approved_symbol %in% robust_genes &
#       ram$receptor_approved_symbol %in% robust_genes
# ram <- ram[ix,]
ix <- omni$source_genesymbol %in% robust_genes &
  omni$target_genesymbol %in% robust_genes
omni <- omni[ix,]

setdiff(my_genes, d$symbol)

# Add the number of cells each donor has in each cluster
d <- left_join(d, obs %>% count(donor, leiden), by = "donor")
d <- as.data.table(d)

# res <- pblapply(seq(nrow(ram)), function(ram_i) {
#   # g1 <- "ADAM15"
#   # g2 <- "ITGAV"
#   g1 <- ram$ligand_approved_symbol[ram_i]
#   g2 <- ram$receptor_approved_symbol[ram_i]
res <- pblapply(seq(nrow(omni)), function(omni_i) {
  # g1 <- "ADAM15"
  # g2 <- "ITGAV"
  g1 <- omni$source_genesymbol[omni_i]
  g2 <- omni$target_genesymbol[omni_i]
  x <- d[symbol == g1]
  y <- d[symbol == g2]
  xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
  # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
  xx_donor <- xx$donor
  xx <- as.matrix(xx[,2:ncol(xx)])
  rownames(xx) <- xx_donor
  yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
  # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
  yy_donor <- yy$donor
  yy <- as.matrix(yy[,2:ncol(yy)])
  rownames(yy) <- yy_donor
  # Remove zeros
  xx[xx == 0] <- NA
  yy[yy == 0] <- NA
  sp <- Hmisc::rcorr(xx, yy, "spearman")
  sp_r <- sp$r[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
  sp_p <- sp$P[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
  ut <- upper.tri(sp_r, diag = FALSE)
  res <- data.table(
    g1 = g1,
    g2 = g2,
    c1 = rownames(sp_r)[row(sp_r)[ut]],
    c2 = rownames(sp_r)[col(sp_r)[ut]],
    estimate = sp_r[ut],
    p.value = sp_p[ut]
  )
})
res <- rbindlist(res)
res <- res[!is.na(res$estimate)]
res$p.value[res$p.value == 0] <- 1
d_n <- d[, .(donors = sum(mean > 0)), by = .(cluster, symbol)]
d_n$cluster <- as.character(d_n$cluster)
res2 <- d_n[res, on = c("symbol" = "g1", "cluster" = "c1")]
res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "g2")] 
colnames(res3) <- c("c1", "g1", "donors1", "c2", "g2", "donors2", "estimate", "p.value")
res <- res3
rm(res2)
rm(res3)
res %<>% arrange(p.value)
res <- unique(res)
res <- res[donors1 >= min_donors & donors2 >= min_donors,]

top_res <- res %>%
  mutate(inter = substr(c1, 1, 1) != substr(c2, 1, 1)) %>%
  filter(donors1 >= min_donors & donors2 >= min_donors) %>%
  # filter(inter) %>%
  arrange(p.value)
top_res %<>%
  group_by(c1, c2) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  mutate(signif = fdr < 0.05) %>%
  ungroup()
  # mutate(signif = p.value < 0.05 / nrow(top_res))
# top_res %>% filter(signif) %>% count(c1, c2)
top_res

qsave(top_res, "results/a20/covarying-abundance/res.qs")

pdf_file <- "results/a20/covarying-abundance/heatmap-signif.pdf"
top_res_dat <- top_res %>%
  group_by(c1, c2) %>%
  summarize(n = sum(signif), .groups = "drop") %>%
  as.data.table
cols <- naturalsort(unique(with(top_res_dat, c(c1, c2))))
top_res_mat <- matrix(0, ncol = length(cols), nrow = length(cols))
rownames(top_res_mat) <- cols
colnames(top_res_mat) <- cols
for (i in seq(nrow(top_res_dat))) {
  row <- top_res_dat[i,]
  c_sort <- rev(naturalsort(c(row$c1, row$c2)))
  top_res_mat[c_sort[2], c_sort[1]] <- row$n
}
top_res_mat <- top_res_mat + t(top_res_mat)
#
legend_labs <- unique(round(seq(min(log10(top_res_mat + 1)), max(log10(top_res_mat + 1)), length.out = 5)))
ht <- Heatmap(
  matrix = log10(top_res_mat + 1),
  col = rev(scico(pal = "oslo", n = 20)),
  row_names_side = "left",
  row_split = str_remove(rownames(top_res_mat), "\\d+"),
  column_split = str_remove(colnames(top_res_mat), "\\d+"),
  row_gap = unit(5, "mm"),
  column_gap = unit(5, "mm"),
  border = TRUE,
  row_title_rot = 0,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  name = "Gene Pairs",
  column_names_rot = 45,
  heatmap_legend_param = list(
    at = legend_labs,
    labels = 10 ^ legend_labs
  )
)
message(pdf_file)
unlink(pdf_file)
pdf(pdf_file, width = 16, height = 13)
draw(ht)
dev.off()

pdf_file <- "results/a20/covarying-abundance/heatmap-signif-major.pdf"
top_res_mat <- top_res %>%
  mutate(c1 = str_remove(c1, "\\d+"), c2 = str_remove(c2, "\\d+")) %>%
  group_by(c1, c2) %>%
  summarize(n = sum(signif), .groups = "drop") %>%
  mutate(c1 = naturalfactor(c1), c2 = naturalfactor(c2)) %>%
  as.data.table
top_res_mat <- dcast(top_res_mat, c1 ~ c2, value.var = "n")
top_res_mat_rownames <- top_res_mat[[1]]
top_res_mat[[1]] <- NULL
top_res_mat <- as.matrix(top_res_mat)
rownames(top_res_mat) <- top_res_mat_rownames
top_res_mat[is.na(top_res_mat)] <- 0
top_res_mat[1:5,1:5]
#
legend_labs <- unique(round(seq(min(log10(top_res_mat + 1)), max(log10(top_res_mat + 1)), length.out = 5)))
ht <- Heatmap(
  matrix = log10(top_res_mat + 1),
  col = rev(scico(pal = "oslo", n = 20)),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  name = "Gene Pairs",
  row_names_side = "left",
  column_names_rot = 0,
  column_names_centered = TRUE,
  heatmap_legend_param = list(
    at = legend_labs,
    labels = 10 ^ legend_labs
  )
)
pdf(pdf_file, width = 3, height = 2)
draw(ht)
dev.off()

top_res %>% filter(g1 %in% my_genes | g2 %in% my_genes) %>%
  head(15)

my_ggsave(
  "qqplot",
  out_dir = "results/a20/covarying-abundance",
  types = "png",
  plot = gg_qqplot(top_res$p.value),
  scale = 1, width = 4, height = 3, units = "in", dpi = 300
)

table(top_res$signif)

# my_res <- top_res %>%
#   # filter(xor(str_detect(c1, "^N"), str_detect(c2, "^N"))) %>%
#   filter(c1 %in% c("N7") | c2 %in% c("N7")) %>%
#   filter(c1 %in% c("C11", "C13") | c2 %in% c("C11", "C13"))
# my_res %<>% mutate(signif = p.value < 0.05 / nrow(my_res))
# head(my_res)


# write_tsv(
#   top_res[1:100,] %>% mutate_if(is.numeric, signif, 2),
#   "analysis/kamil/colitis_immune_cell_communication.tsv.new"
# )
write_tsv(
  top_res %>% mutate_if(is.numeric, signif, 2) %>% head(1e4),
  "results/a20/covarying-abundance/communication.tsv"
)

# Case-only
########################################################################
res_case <- rbindlist(lapply(c("Case", "Control"), function(this_case) {
  res <- pblapply(seq(nrow(omni)), function(omni_i) {
    # lig1 <- "ADAM15"
    # rec2 <- "ITGAV"
    lig1 <- omni$source_genesymbol[omni_i]
    rec2 <- omni$target_genesymbol[omni_i]
    x <- d[symbol == lig1][case == this_case]
    y <- d[symbol == rec2][case == this_case]
    xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
    # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
    xx_donor <- xx$donor
    xx <- as.matrix(xx[,2:ncol(xx)])
    rownames(xx) <- xx_donor
    yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
    # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
    yy_donor <- yy$donor
    yy <- as.matrix(yy[,2:ncol(yy)])
    rownames(yy) <- yy_donor
    # Remove zeros
    xx[xx == 0] <- NA
    yy[yy == 0] <- NA
    sp <- Hmisc::rcorr(xx, yy, "spearman")
    sp_r <- sp$r[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
    sp_p <- sp$P[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
    ut <- upper.tri(sp_r, diag = FALSE)
    res <- data.table(
      lig1 = lig1,
      rec2 = rec2,
      c1 = rownames(sp_r)[row(sp_r)[ut]],
      c2 = rownames(sp_r)[col(sp_r)[ut]],
      estimate = sp_r[ut],
      p.value = sp_p[ut]
    )
  })
  res <- rbindlist(res)
  res <- res[!is.na(res$estimate)]
  res$p.value[res$p.value == 0] <- 1
  d_n <- d[case == this_case, .(donors = sum(mean > 0)), by = .(cluster, symbol)]
  d_n$cluster <- as.character(d_n$cluster)
  res2 <- d_n[res, on = c("symbol" = "lig1", "cluster" = "c1")]
  res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "rec2")] 
  colnames(res3) <- c("c1", "g1", "donors1", "c2", "g2", "donors2", "estimate", "p.value")
  res <- res3
  rm(res2)
  rm(res3)
  res$case <- this_case
  return(res)
}))
res_case <- unique(res_case)
res_case$name <- with(res_case, glue("{c1}|{g1}|{c2}|{g2}"))
shared_names <- intersect(
  res_case[case == "Case"]$name,
  res_case[case != "Case"]$name
)
res_case$shared <- res_case$name %in% shared_names

res_case <- res_case %>%
  mutate(inter = substr(c1, 1, 1) != substr(c2, 1, 1)) %>%
  # filter(donors1 >= min_donors & donors2 >= min_donors) %>%
  # filter(inter) %>%
  arrange(p.value)
res_case %<>%
  group_by(c1, c2, case) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  mutate(signif = fdr < 0.05) %>%
  ungroup()
  # mutate(signif = p.value < 0.05 / nrow(res_case))
res_case <- as.data.table(res_case)
qsave(res_case, "results/a20/covarying-abundance/res_case.qs")

res_case %>% filter(signif) %>% count(c1, c2)
res_case

plot_communication(res_case %>% filter(g2 == "IL26", c1 == "E25", c2 == "CT6"))

plot_communication(res_case %>% filter(g1 == "PTPRS", c1 == "E5", c2 == "CT9"))

plot_communication(res_case %>% filter(signif) %>% head)

plot_communication(res_case %>% filter(g1 == "IL7R", g2 == "IL7") %>% head(2))

plot_communication(res %>% filter(g1 == "IL7R", g2 == "IL7") %>% head(10))




# res_case <- res_case[with(res_case, order(case, name)),]

# stopifnot(all(
#   res_case[case == "Case"]$name == res_case[case == "Control"]$name
# ))

# x <- data.frame(
#   case_rho = res_case[case == "Case"]$estimate,
#   control_rho = res_case[case != "Case"]$estimate,
#   case_p = res_case[case == "Case"]$p.value,
#   control_p = res_case[case != "Case"]$p.value
# )
# ggplot(x) +
#   aes(-log10(case_p), -log10(control_p)) +
#   stat_binhex(bins = 100) +
#   scale_fill_gradientn(
#     trans = "log10",
#     colors = scico::scico(
#       n = 20, palette = "bilbao", direction = 1
#     )[4:20]
#   )

plot_communication(
  res_case[which(x$case_p < 10^-4.5 & x$control_p > 1e-2),]
)

# g1 <- "PTPRS"
# g2 <- "TNF"
# c1 <- "E5"
# c2 <- "CT9"
#     x <- d[cluster == c1][symbol == g1]# %>% filter(cluster == c1, symbol == g1)
#     y <- d[cluster == c2][symbol == g2]#%>% filter(cluster == c2, symbol == g2)
#     x$percent[x$percent == 0] <- NA
#     y$percent[y$percent == 0] <- NA
#     x_lm <- data.frame(x = x$percent, y = y$percent, case = x$case)
# anova(lm(y ~ x, x_lm), lm(y ~ 0 + case + x:case, x_lm))

res_case_wide <- dcast(
  res_case %>% as.data.table,
  g1 + g2 + c1 + c2 + donors1 + donors2 + inter ~ case,
  value.var = "p.value",
  fill = 1
) %>% arrange(Case)
res_case_wide <- res_case_wide %>% filter(Case < 1 | Control < 1) %>% mutate(Case = -log10(Case), Control = -log10(Control))
p <- ggplot(res_case_wide) +
  geom_scattermost(res_case_wide %>% select(Case, Control))
my_ggsave(
  "point-pvalue-case-vs-control",
  out_dir = "results/a20/covarying-abundance/case",
  types = "pdf",
  plot = p,
  scale = 1, width = 8, height = 6, units = "in", dpi = 300
)

res_case_wide <- dcast(
  res_case %>% as.data.table,
  g1 + g2 + c1 + c2 + donors1 + donors2 + inter ~ case,
  value.var = "estimate",
  fill = 0
) %>% arrange(Case)
p <- ggplot(res_case_wide) +
  geom_scattermost(res_case_wide %>% select(Case, Control)) +
  labs(x = "Case", y = "Control", title = "Spearman's rho for all pairs of genes and clusters")
my_ggsave(
  "point-spearman-case-vs-control",
  out_dir = "results/a20/covarying-abundance/case",
  types = "pdf",
  plot = p,
  scale = 1, width = 8, height = 6, units = "in", dpi = 300
)


plot_communication2(res_case %>% head(10))

plot_communication2(res_case %>% filter(g1 == "PTPRS") %>% head(2))

d %>% filter(symbol == "TNF", cluster == "CT9") %>% arrange(percent)

obs %>% count(donor, leiden) %>% filter(leiden == "CT9")

# plot_communication2(
#   res_case[which(x$case_p < 10^-5 & x$control_p > 1e-2),]
# )

# Plots
########################################################################


# plot_communication(top_res %>% filter(signif))
# plot_communication(top_res %>% head(300))
plot_communication(top_res %>% filter(signif, g1 %in% my_genes | g2 %in% my_genes))

plot_communication(top_res %>% head(15))

plot_communication(top_res %>% filter(g1 %in% c("IL26") | g2 %in% "IL26") %>% head(5))

plot_communication(top_res %>% filter(g1 %in% c("CXCR6") | g2 %in% "CXCR6") %>% head(5))

top_res %>% dplyr::filter(lig1 == "LILRB1", p.value < 1e-3) %>% plot_communication

cc_graph <- top_res %>% filter(signif)
cc_graph <- cc_graph %>%
  mutate(from = c1, to = c2) %>%
  count(c1, c2) %>%
  arrange(n)
cc_graph <- tbl_graph(edges = cc_graph)
p <- ggraph(cc_graph, layout = "linear", circular = TRUE) + 
  geom_edge_arc(
    aes(
      start_cap = label_rect(node1.name, padding = margin(3, 3, 3, 3, "mm")),
      end_cap = label_rect(node2.name, padding = margin(3, 3, 3, 3, "mm")),
      # colour = stat(index),
      color = n,
      edge_width = stat(index),
      # alpha = stat(index)
      alpha = n
    )
  ) +
  # scale_edge_colour_viridis(option = "Grays", direction = 1) +
  scale_edge_color_gradientn(colors = scico::scico(10)[c(2:10)]) + #, guide = FALSE) +
  # scale_edge_color_manual(values = scico::scico(10)[c(2:10)]) +
  scale_edge_width_continuous(range = c(0.1, 2), guide = FALSE) +
  scale_edge_alpha_continuous(range = c(0.5, 0.8), guide = FALSE) +
  geom_node_text(aes(label = name), size = 8) +
  guides(edge_color = guide_edge_colorbar(title = "Edges", barheight = 10)) +
  coord_fixed(clip = "off") +
  theme(panel.border = element_blank())
my_ggsave(
  glue("cc_web"),
  out_dir = "figures/case_control/communication",
  types = "pdf",
  plot = p,
  scale = 1, width = 4, height = 3, units = "in", dpi = 300
)

set.seed(42)
p <- ggraph(cc_graph, layout = "fr") +
  geom_edge_arc(
    aes(
      start_cap = label_rect(node1.name, padding = margin(3, 3, 3, 3, "mm")),
      end_cap = label_rect(node2.name, padding = margin(3, 3, 3, 3, "mm")),
      # colour = stat(index),
      color = n,
      edge_width = stat(index),
      # alpha = stat(index)
      alpha = n
    ),
    strength = 0.33
  ) +
  # scale_edge_colour_viridis(option = "Grays", direction = 1) +
  scale_edge_color_gradientn(colors = scico::scico(10)[c(2:10)]) + #, guide = FALSE) +
  # scale_edge_color_manual(values = scico::scico(10)[c(2:10)]) +
  scale_edge_width_continuous(range = c(0.1, 2), guide = FALSE) +
  scale_edge_alpha_continuous(range = c(0.6, 0.8), guide = FALSE) +
  geom_node_text(aes(label = name), size = 8) +
  guides(edge_color = guide_edge_colorbar(title = "Edges", barheight = 10)) +
  coord_fixed(clip = "off") +
  theme(panel.border = element_blank())
#p
my_ggsave(
  glue("cc_graph_fr"),
  out_dir = "figures/case_control/communication",
  types = "pdf",
  plot = p,
  scale = 1, width = 5, height = 4, units = "in", dpi = 300
)


# Cell communication for all pairs of genes
########################################################################

# Get all pairs of clusters.
obs$cluster <- obs$leiden
all_clusters <- sort(unique(obs$cluster))
cc <- expand.grid(all_clusters, all_clusters)
cc <- cc[cc$Var1 != cc$Var2,]

log2cpm_means <- Matrix::rowMeans(log2cpm)

ix <- log2cpm_means > 0.5
plot(hist(log2cpm_means[ix], breaks = 100))
sum(ix)
my_ensembl_ids <- rownames(log2cpm)[ix]
#
d <- as.data.table(summary(log2cpm[my_ensembl_ids,]))
d$symbol <- my_ensembl_ids[d$i]
d$cluster <- obs$cluster[d$j]
d$donor <- obs$donor[d$j]
y <- obs %>% group_by(cluster, donor) %>% count() %>% as.data.table
d <- y[d, on = c("cluster", "donor")]
d <- d[
  ,
  .(
    mean = mean(c(x, rep(0, n[1] - length(x)))),
    percent = sum(x > 0) / n[1]
  ),
  by = .(cluster, donor, symbol)
]
setkey(d, symbol, cluster, donor)
# d$percent_logit <- car::logit(d$percent)
# The gene is robustly detected in some cluster for 20 donors.
d_nonzero <- d[,.(nonzero = sum(mean > 0)),by = .(cluster, symbol)]
robust_genes <- d_nonzero[nonzero > 20]$symbol %>% unique
d <- d[symbol %in% robust_genes]
# Add zeros when we don't have anything.
d <- tidyr::complete(
  d, cluster, donor, symbol, fill = list(mean = 0, percent = 0)
) %>% as.data.table
x <- obs %>% select(donor, case) %>% unique
donor_to_case <- x$case
names(donor_to_case) <- x$donor
d$case <- donor_to_case[d$donor]
stopifnot(all(d$symbol %in% my_ensembl_ids))
my_ensembl_ids <- intersect(my_ensembl_ids, d$symbol)
#
my_ensembl_id_pairs <- expand.grid(my_ensembl_ids, my_ensembl_ids)
my_ensembl_id_pairs <- my_ensembl_id_pairs[
  my_ensembl_id_pairs$Var1 != my_ensembl_id_pairs$Var2,
]
my_ensembl_id_pairs$Var1 <- as.character(my_ensembl_id_pairs$Var1)
my_ensembl_id_pairs$Var2 <- as.character(my_ensembl_id_pairs$Var2)
stopifnot(all( my_ensembl_ids %in% d$symbol))


res_file <- "cache/immune_nuclei_communication_225_genes.rds"
if (!file.exists(res_file)) {
  res <- pblapply(seq(nrow(my_ensembl_id_pairs)), function(ens_i) {
    # lig1 <- "ADAM15"
    # rec2 <- "ITGAV"
    # lig1 <- ram$ligand_approved_symbol[ram_i]
    # rec2 <- ram$receptor_approved_symbol[ram_i]
    lig1 <- my_ensembl_id_pairs[ens_i,1]
    rec2 <- my_ensembl_id_pairs[ens_i,2]
    x <- d[symbol == lig1]    
    y <- d[symbol == rec2]
    xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
    # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
    xx_donor <- xx$donor
    xx <- as.matrix(xx[,2:ncol(xx)])
    rownames(xx) <- xx_donor
    yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
    # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
    yy_donor <- yy$donor
    yy <- as.matrix(yy[,2:ncol(yy)])
    rownames(yy) <- yy_donor
    # Remove zeros
    xx[xx == 0] <- NA
    yy[yy == 0] <- NA
    sp <- Hmisc::rcorr(xx, yy, "spearman")
    sp_r <- sp$r[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
    sp_p <- sp$P[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
    ut <- upper.tri(sp_r, diag = FALSE)
    res <- data.table(
      lig1 = lig1,
      rec2 = rec2,
      c1 = rownames(sp_r)[row(sp_r)[ut]],
      c2 = rownames(sp_r)[col(sp_r)[ut]],
      estimate = sp_r[ut],
      p.value = sp_p[ut]
    )
  })
  res <- rbindlist(res)
  res <- res[!is.na(res$estimate)]
  d_n <- d[, .(donors = sum(mean > 0)), by = .(cluster, symbol)]
  d_n$cluster <- as.character(d_n$cluster)
  res2 <- d_n[res, on = c("symbol" = "lig1", "cluster" = "c1")]
  res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "rec2")] 
  colnames(res3) <- c("c1", "lig1", "donors1", "c2", "rec2", "donors2", "estimate", "p.value")
  res <- res3
  rm(res2)
  rm(res3)
  saveRDS(res, res_file)
} else {
  res <- readRDS(res_file)
}

top_res <- res %>% filter(donors1 >= 15, donors2 >= 15) %>% arrange(p.value)
top_res <- top_res %>% mutate(signif = p.value < 0.05 / nrow(top_res))
top_res$inter <- substr(top_res$c1, 1, 1) != substr(top_res$c2, 1, 1)
top_res %>% filter(signif) %>% count(c1, c2)

my_res <- top_res %>%
  # filter(xor(str_detect(c1, "^N"), str_detect(c2, "^N"))) %>%
  filter(c1 %in% c("N7") | c2 %in% c("N7")) %>%
  filter(c1 %in% c("C11", "C13") | c2 %in% c("C11", "C13"))
my_res %<>% mutate(signif = p.value < 0.05 / nrow(my_res))

my_res$lig1_symbol <- ensembl_to_symbol[my_res$lig1]
my_res$rec2_symbol <- ensembl_to_symbol[my_res$rec2]
ix <- (
  !str_detect(my_res$lig1_symbol, "^(RP|MT)") &
  !str_detect(my_res$rec2_symbol, "^(RP|MT)")
)
head(my_res[ix,])

table(top_res$signif)

top_res$lig1_symbol <- ensembl_to_symbol[top_res$lig1]
top_res$rec2_symbol <- ensembl_to_symbol[top_res$rec2]
ix <- (
  !str_detect(top_res$lig1_symbol, "^(RP|MT)") &
  !str_detect(top_res$rec2_symbol, "^(RP|MT)")
)

top_res[ix,] %>%
  filter(inter) %>%
  head(20) %>%
  select(c1, c2, lig1_symbol, rec2_symbol, p.value, signif)

# write_tsv(
#   top_res[1:100,] %>% mutate_if(is.numeric, signif, 2),
#   "analysis/kamil/colitis_immune_cell_communication.tsv.new"
# )
write_tsv(
  top_res %>% mutate_if(is.numeric, signif, 2),
  "results/allgenes_immune_nuclei_cell_communication.tsv"
)

#left_join(d, obs %>% select(donor, case), by = "donor")

# gg_qqplot(res$p.value)


## Linear model for cases and controls
########################################################################

# Lots of bugs here... something is wrong

x <- d[cluster == c1][symbol == lig1]
y <- d[cluster == c2][symbol == rec2]
x$g1_pct <- x$percent
x$g2_pct <- y$percent
plot(x$g1_pct, x$g2_pct)

res %>% filter(g1 == "IL26" | g2 == "IL26") %>%
  arrange(p.value)

  select(g1, g2) %>% unique

res_lm <- pblapply(seq(nrow(omni)), function(omni_i) {

  # lig1 <- "IL10RB"
  lig1 <- "IL10RA"
  rec2 <- "IL26"

  # lig1 <- omni$source_genesymbol[omni_i]
  # rec2 <- omni$target_genesymbol[omni_i]

  # https://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
  # This works, about 30 seconds for one pair of genes
  fit <- rbindlist(pblapply(seq(nrow(cc)), function(i) {
    # i <- 1129
    c1 <- cc[i,1]
    c2 <- cc[i,2]
    # c1 <- "E1"
    g1 <- "IL20RA"
    # c2 <- "CT4"
    g2 <- "IL26"
    # x <- d[cluster == "E1"][symbol == "IL20RA"]
    # y <- d[cluster == "CT4"][symbol == "IL26"]
    x <- d[cluster == c1][symbol == g1]
    y <- d[cluster == c2][symbol == g2]
    x$g1 <- x$percent
    x$g2 <- y$percent
    x$case <- x$case == "Case"
    if (sum(x$g1 > 0) > 5 && sum(x$g2 > 0) > 5) {
      # plot(x$g1, x$g2)
      fit1 <- lm(g2 ~ g1, x, x = TRUE)
      fit2 <- lm(g2 ~ g1 * case, x, x = TRUE)
      r1 <- broom::tidy(fit1)
      r1$test <- "fit1"
      r2 <- broom::tidy(fit2)
      r2$test <- "fit2"
      r3 <- broom::tidy(anova(fit1, fit2))
      r3$test <- "anova"
      retval <- rbindlist(list(r1, r2, r3), fill = TRUE)
      retval$c1 <- c1
      retval$g1 <- g1
      retval$c2 <- c2
      retval$g2 <- g2
      return(retval)
    }
    return(NULL)
  }))

  a <- cbind(
    intercept = rep(1, length(x$g1)),
    g1        = x$g1,
    case      = 1 * x$case,
    g1_case   = (1 * x$case) * x$g1
  )
  b <- x$g2
  solve(crossprod(a), crossprod(a, b))

  n <- length(x$g1)
  a <- cbind(
    intercept = rep(rep(1, n), 2),
    g1        = c(x$g1, rep(0, n)),
    case      = c(1 * x$case, rep(0, n)),
    g1_case   = c((1 * x$case) * x$g1, rep(0, n)),
    g3        = c(rep(0, n), x$g1),
    case      = c(rep(0, n), 1 * x$case),
    g3_case   = c(rep(0, n), (1 * x$case) * x$g1)
  )
  b <- c(x$g2, x$g2)
  solve(crossprod(a), crossprod(a, b))

  # solve(crossprod(fit2$x), crossprod(fit2$x,x$g2))

  fit %>% arrange(p.value) %>% head

  x <- d[symbol == lig1][cluster == c1]
  y <- d[symbol == rec2][cluster == c2]
  # stopifnot(all(x$cluster == y$cluster))
  stopifnot(all(x$donor == y$donor))

  plot(x$percent, y$percent)

  x$g1 <- x$percent
  x$g2 <- y$percent

  all_clusters <- sort(unique(obs$cluster))
  cc <- expand.grid(all_clusters, all_clusters)
  cc <- cc[cc$Var1 != cc$Var2,]
  cc[,1] <- as.character(cc[,1])
  cc[,2] <- as.character(cc[,2])

  # cc_IL20RA_IL26_E1_CT4

  fit <- rbindlist(lapply(seq(nrow(cc)), function(i) {
    # c1 <- cc[i,1]
    # c2 <- cc[i,2]

    g1 <- x$percent[x$cluster == c1]
    # g1[g1 == 0] <- NA
    g2 <- y$percent[y$cluster == c2]
    # g2[g2 == 0] <- NA
    case <- y$case[y$cluster == c2] == "Case"
    fit <- broom::tidy(lm(g1 ~ g2))
    fit$g1 <- lig1
    fit$c1 <- c1
    fit$g2 <- rec2
    fit$c2 <- c2
    fit

    plot(g1, g2)

    cor(g1, g2, na.rm = TRUE)

top_res %>% filter(g1 %in% c("IL26") | g2 %in% "IL26") %>% head(5)

  }))

  fit %>%
    arrange(p.value) %>%
    filter(term == "caseTRUE:g2") %>%
    head

  plot(g1, g2, col = ifelse(case, "red", "blue"))

  xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
  # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
  xx_donor <- xx$donor
  xx <- as.matrix(xx[,2:ncol(xx)])
  rownames(xx) <- xx_donor
  yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
  # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
  yy_donor <- yy$donor
  yy <- as.matrix(yy[,2:ncol(yy)])
  rownames(yy) <- yy_donor
  # Remove zeros
  xx[xx == 0] <- NA
  yy[yy == 0] <- NA

  sp <- Hmisc::rcorr(xx, yy, "spearman")
  sp_r <- sp$r[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
  sp_p <- sp$P[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
  ut <- upper.tri(sp_r, diag = FALSE)
  res <- data.table(
    lig1 = lig1,
    rec2 = rec2,
    c1 = rownames(sp_r)[row(sp_r)[ut]],
    c2 = rownames(sp_r)[col(sp_r)[ut]],
    estimate = sp_r[ut],
    p.value = sp_p[ut]
  )
})





#source("R/recluster.R")
#a_slugs <- c(
#  "subset-bcell/a12_4_4_b5_1_min_genes500_n_pcs16",
#  "subset-myeloid/a12_4_4_m3_min_genes500_n_pcs8",
#  "subset-tcell/a12_4_4_t4_cd4_2_min_genes500_n_pcs12",
#  "subset-tcell/a12_4_4_t4_cd8_1_min_genes500_n_pcs12",
#  "nuclei/n3_min_genes500_n_pcs30"
#)
#for (a in a_slugs) {
#  a_dir <- dirname(a)
#  a_name <- basename(a)
#  a_file <- glue("results/{a_dir}/{a_name}/data/{a_name}.rds")
#  #
#  message(a_dir)
#  message(a_name)
#  message(a_file)
#  message(file.exists(a_file))
#  #
#  stopifnot(file.exists(a_file))
#  #
#  print_status(glue("Reading {a_file}"))
#  a1 <- readRDS(a_file)
#  print_status(glue("done"))
#  #
#  a1 <- recluster(a1, a)
#  #
#  if (!"drug" %in% colnames(a1$obs)) {
#    a1$obs <- left_join(a1$obs, sample_therapy, by = "donor")
#  }
#  assign(a_name, a1)
#}
#a_labels <- c(
#  "B cell subsets",
#  "Myeloid cell subsets",
#  "CD4 T cell subsets",
#  "CD8 T cell subsets",
#  "Epithelial cell subsets"
#)
## for (i in seq_along(a_labels)) {
##   de_files <- Sys.glob(
##     glue("results/{a_slugs[i]}/figures/pb_cluster/case-control/de_case-vs-control_cluster-*.tsv.gz")
##   )
##   dat <- rbindlist(lapply(de_files, function(de_file) {
##     x <- fread(de_file)
##     x$cluster <- str_remove(str_split_fixed(basename(de_file), "-", 4)[,4], "\\..+")
##     x
##   }))
## }

# # Nuclei data
# ########################################################################
# min_genes <- 200
# n_pcs <- 25
# analysis_name <- as.character(glue("n5_min_genes{min_genes}_n_pcs{n_pcs}"))
# nuclei_file <- as.character(glue("results/nuclei/{analysis_name}/data/{analysis_name}.rds"))
# file.exists(nuclei_file)
# nuclei <- readRDS(nuclei_file)
# stopifnot(all(nuclei$obs$cell == colnames(nuclei$counts)))
# 
# 
# # Immune cell data
# ########################################################################
# min_genes <- 500
# n_pcs <- 20
# analysis_name <- as.character(glue("a12_4_4_min_genes{min_genes}_n_pcs{n_pcs}"))
# a12_4_4_file <- as.character(glue("results/{analysis_name}/data/{analysis_name}.rds"))
# file.exists(a12_4_4_file)
# cells <- readRDS(a12_4_4_file)
# stopifnot(all(cells$obs$cell == colnames(cells$counts)))
# 
# 
# # Integrated analysis
# ########################################################################
# 
# meta_cols <- intersect(colnames(nuclei$obs), colnames(cells$obs))
# meta_cols <- meta_cols[!str_detect(meta_cols, "^PC")]
# # meta_cols <- meta_cols[!str_detect(meta_cols, "^UMAP")]
# # meta_cols <- meta_cols[!str_detect(meta_cols, "^leiden")]
# meta_cols
# 
# counts_rows <- intersect(rownames(nuclei$counts), rownames(cells$counts))
# length(counts_rows) / nrow(nuclei$counts)
# length(counts_rows) / nrow(cells$counts)
# 
# nuclei$obs$leiden <- sprintf("N%s", nuclei$obs$leiden)
# cells$obs$leiden <- sprintf("C%s", cells$obs$leiden)
# 
# obs <- rbind(cells$obs[,meta_cols], nuclei$obs[,meta_cols])
# #
# counts <- cbind(cells$counts[counts_rows,], nuclei$counts[counts_rows,])
# log2cpm <- do_log2cpm(counts, median(colSums(counts)))
# rm(counts)
# #
# stopifnot(all(obs$cell == colnames(log2cpm)))
# 
# shared_donors <- intersect(
#   nuclei$obs$donor,
#   cells$obs$donor
# )
# 
# length(shared_donors)
# length(unique(nuclei$obs$donor))
# length(unique(cells$obs$donor))
# 
# ix <- obs$donor %in% shared_donors
# obs <- obs[ix,]
# log2cpm <- log2cpm[,ix]
# stopifnot(all(obs$cell == colnames(log2cpm)))
# 
# genes <- genes[genes$ensembl_id %in% rownames(log2cpm),]
# 
# 
# # Analysis
# ########################################################################
# 
# # Free memory
# rm(nuclei)
# rm(cells)
# 
# # This takes about 3 hours on my macbook
# keep_channels <- (
#   obs %>%
#     dplyr::count(channel) %>%
#     dplyr::filter(n > 100)
# )$channel
# keep_cells <- which(obs$channel %in% keep_channels)
# length(keep_cells)
# length(keep_cells) / nrow(obs)
# obs <- obs[keep_cells,]
# log2cpm <- log2cpm[,keep_cells]
# stopifnot(all(obs$cell == colnames(log2cpm)))
# #

