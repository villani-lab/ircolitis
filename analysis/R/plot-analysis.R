source("R/functions/mpn65.R")
# source("R/colors.R")

#' Plot all results from do_analysis()
#' @param analysis_name A string with the name of the analysis variable.
#' @param rowname_key A named character vector to translate rownames to meaningful names.
plot_analysis <- function(analysis_name, out_dir, rowname_key, exclude_genes, do_pb = FALSE) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  #
  # Get the analysis variable from the parent environment
  a2 <- get(analysis_name)
  # Log2CPM
  a2$log2cpm <- do_log2cpm(a2$counts, total = median(colSums(a2$counts)))
  x1 <- sort(unique(a2$de$group))
  x2 <- sort(unique(a2$obs$leiden))
  if (length(x1) != length(x2) || !all(x1 == x2)) {
    a2$de <- presto::wilcoxauc(a2$log2cpm, a2$obs$leiden)
  }
  #
  stopifnot(all(colnames(a2$log2cpm) == a2$obs$cell))
  n_clusters <- length(unique(a2$de$group))
  # cluster_colors <- S4Vectors::unname(pals::alphabet())
  cluster_colors <- mpn65
  cluster_colors <- rep(cluster_colors, length.out = n_clusters)
  # names(cluster_colors) <- seq_along(cluster_colors)
  names(cluster_colors) <- unique(a2$de$group)
  #
  # a2$obs$chemistry <- "3p"
  # a2$obs$chemistry[str_detect(a2$obs$channel, "_5p")] <- "5p"
  #
  my_cols <- c("TRBV", "TRAV", "TRB_cdr3", "TRA_cdr3")
  if (all(my_cols %in% colnames(a2$obs))) {
    a2$obs$has_tcr <- !is.na(a2$obs$TRB_cdr3) & !is.na(a2$obs$TRA_cdr3)
  }
  #
  my_cols <- c("IGLC", "IGHC", "IGKC")
  if (all(my_cols %in% colnames(a2$obs))) {
    a2$obs$has_bcr <- !is.na(a2$obs$IGH_cdr3) & (
      !is.na(a2$obs$IGL_cdr3) | !is.na(a2$obs$IGK_cdr3)
    )
  }
  #
  n_donors <- length(unique(a2$obs$donor))
  n_channels <- length(unique(a2$obs$channel))
  n_samples <- length(unique(a2$obs$sample))
  a2$obs$chemistry <- "3p"
  a2$obs$chemistry[str_detect(a2$obs$channel, "_5p")] <- "5p"
  a2$de$symbol <- rowname_key[a2$de$feature]

  stopifnot(nrow(a2$obs) == ncol(a2$counts))
  y <- with(a2$obs, model.matrix(~ 0 + factor(leiden)))
  y <- as(y, "dgCMatrix")
  leiden_percents <- as((a2$counts > 0) %*% y, "dgCMatrix")
  colnames(leiden_percents) <- str_replace(colnames(leiden_percents), "factor\\(leiden\\)", "")
  n_cells_per_leiden <- a2$obs %>% count(leiden)
  stopifnot(n_cells_per_leiden$leiden == colnames(leiden_percents))
  leiden_percents <- sweep(leiden_percents, 2, as.numeric(n_cells_per_leiden$n), "/")

  if ("mcv" %in% names(a2)) {
    d_mcv <- a2$mcv %>% pivot_longer( cols = c('mcv_loss', 'rec_loss'))
    p <- ggplot() +
      stat_summary(
        data = d_mcv,
        mapping = aes(x = as.integer(k), y = value, group = name, fill = name),
        fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95),
        geom = "ribbon",
        alpha = 0.3
      ) +
      stat_summary(
        data = d_mcv,
        mapping = aes(x = as.integer(k), y = value, group = name, color = name),
        geom = "line",
        alpha = 0.8,
        size = 0.5
      ) +
      geom_vline(data = NULL, xintercept = as.integer(a2$mcv$optimal_k[1]), size = 0.2) +
      geom_point(
        data = d_mcv,
        mapping = aes(x = as.integer(k), y = value, group = name, color = name),
        size = 0.1, position = position_jitter(width = 0.1)
      ) +
      scale_color_manual(
        values = c("red", "grey50"), name = NULL
      ) +
      scale_fill_manual(
        values = c("red", "grey50"), name = NULL
      ) +
      labs(
        title = glue("MCV loss for top k PCs ({max(d_mcv$rep)} reps, best is {a2$mcv$optimal_k[1]})"),
        x = "Top k PCs",
        y = "MCV loss"
      ) +
      theme(
        legend.position = c(1, 1),
        legend.just = c(1, 1),
        legend.background = element_blank()
      )
    my_ggsave(
      "mcv-loss",
      out_dir = out_dir,
      plot = p,
      type = "pdf",
      scale = 1.8, width = 3, height = 2, units = "in", dpi = 300
    )
  }

  #if ("pca_h" %in% names(a2)) {
  #  n_pcs <- ncol(a2$pca_h)
  #  for (i in seq(n_pcs - 1)) {
  #    i <- 1
  #    pc_i <- sprintf("PC%s", i)
  #    pc_j <- sprintf("PC%s", i + 1)
  #    #
  #    p <- ggplot(a2$obs %>% mutate(leiden = naturalfactor(leiden))) +
  #      aes_string(x = pc_i, y = pc_j, color = "leiden") +
  #      scale_color_manual(values = mpn65, guide = "none") +
  #      scattermore::geom_scattermore(pixels = c(1024, 1024), pointsize = 2)
  #    my_ggsave(
  #      glue("pca-{pc_i}-{pc_j}"),
  #      out_dir = file.path(out_dir, "pca"),
  #      plot = p,
  #      type = "pdf",
  #      scale = 1.8, width = 3, height = 3, units = "in", dpi = 300
  #    )
  #  }
  #}

  if ("pca_h" %in% names(a2)) {
    n_pcs <- ncol(a2$pca_h)
    # for (i in seq(n_pcs)) {

      x <- a2$obs %>% select(leiden, starts_with("PC")) %>% pivot_longer(-leiden) %>%
        mutate(
          leiden = naturalfactor(leiden),
          name = naturalfactor(name)
        )
      p <- ggplot(x) +
        aes(y = leiden, x = value, fill = leiden) +
        scale_fill_manual(values = mpn65, guide = "none") +
        stat_summary(geom = "bar", fun.data = mean_se) +
        scale_y_discrete(limits = rev(levels(x$leiden))) +
        facet_row(vars(name), scales = "free_x") +
        labs(x = "PC", y = NULL) +
        theme(
          panel.spacing = unit(0.5, "lines"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
      my_ggsave(
        "pca-bars",
        out_dir = file.path(out_dir, "pca"),
        plot = p,
        type = "pdf",
        scale = 1, width = 20, height = 10, units = "in", dpi = 300
      )

    # }
  }

  if (all(c("donor", "case") %in% colnames(a2$obs)) && n_donors > 1) {
    x <- a2$obs %>%
      dplyr::group_by(donor, case) %>%
      dplyr::summarize(n_cells = n())
    p <- ggplot(x) +
      aes(x = n_cells, y = reorder(donor, n_cells), fill = case) +
      geom_colh() +
      scale_fill_manual(values = pals::okabe()[2:3]) +
      # scale_x_continuous(labels = scales::comma) +
      scale_x_continuous(trans = "log10") +
      # scale_y_discrete(expand = c(0.035, 0)) +
      annotation_logticks(side = "b") +
      labs(y = NULL, x = "Cells", fill = NULL) +
      theme(panel.grid.major.x = element_line(size = 1))
    fig_height <- 0.2 * nrow(x)
    my_ggsave(
      "bars-donor",
      out_dir = out_dir,
      plot = p,
      type = "pdf",
      scale = 1.8, width = 3.5, height = fig_height, units = "in", dpi = 300
    )
  }

  if (all(c("leiden") %in% colnames(a2$obs))) {
    x <- a2$obs %>%
      dplyr::mutate(leiden = append_n(naturalfactor(leiden))) %>%
      dplyr::group_by(leiden) %>%
      dplyr::summarize(n_cells = n())
    p <- ggplot(x) +
      aes(x = n_cells, y = reorder(leiden, n_cells), fill = naturalfactor(leiden)) +
      geom_colh() +
      scale_fill_manual(values = mpn65, guide = "none") +
      scale_x_continuous(labels = scales::comma) +
      scale_y_discrete(position = "r") +
      # scale_x_continuous(trans = "log10", labels = function(x) scales::comma(x, accuracy = 1)) +
      # annotation_logticks(side = "b") +
      # scale_y_discrete(expand = c(0.035, 0)) +
      labs(y = NULL, x = "Cells", fill = NULL) +
      theme(panel.grid.major.x = element_line(size = 1))
    fig_height <- 0.2 * nrow(x)
    my_ggsave(
      "bars-cluster",
      out_dir = out_dir,
      plot = p,
      type = "pdf",
      scale = 1.8, width = 3.5, height = fig_height, units = "in", dpi = 300
    )
  }

  if (all(c("channel", "case") %in% colnames(a2$obs)) && n_channels > 1) {
    x <- a2$obs %>%
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
      labs(y = NULL, x = "Cells", fill = NULL) +
      theme(panel.grid.major.x = element_line(size = 1))
    fig_height <- 0.2 * nrow(x)
    my_ggsave(
      "bars-channel",
      out_dir = out_dir,
      plot = p,
      type = "pdf",
      scale = 1.8, width = 6.5, height = fig_height, units = "in", dpi = 300
    )
  }

# Abundance of each cluster per donor heatmap
########################################################################
my_cols <- c("leiden", "donor", "case", "class_short") #, "chemistry")
if (all(my_cols %in% colnames(a2$obs)) && n_donors > 1) {

  x <- a2$obs %>%
    # dplyr::group_by(leiden, donor, case, class_short, chemistry) %>%
    dplyr::group_by(leiden, donor, case, class_short) %>%
    dplyr::count() %>%
    # dplyr::group_by(leiden) %>%
    # dplyr::mutate(percent_cluster = 100 * n / sum(n)) %>%
    dplyr::group_by(donor) %>%
    dplyr::mutate(
      percent_donor = 100 * n / sum(n),
      total_cells   = sum(n)
    )
  #
  set.seed(2)
  # x <- seriate_dataframe(x, "leiden", "donor", "percent_donor")
  #
  # x$percent_donor <- log10(x$percent_donor)
  # x$percent_donor[!is.finite(x$percent_donor)] <- 0
  #
  sorted_donor <- umap_sorted_rows(x, "donor ~ leiden", "percent_donor")
  sorted_leiden <- umap_sorted_rows(x, "leiden ~ donor", "percent_donor")
  x$donor <- factor(x$donor, sorted_donor)
  x$leiden <- factor(x$leiden, sorted_leiden)
  #
  x <- as.data.table(x)
  #
  xmat <- dcast(x, donor ~ leiden, value.var = "percent_donor", fill = 0)
  xmat_rownames <- xmat$donor
  xmat <- as.matrix(xmat[,2:ncol(xmat)])
  rownames(xmat) <- xmat_rownames
  #
  mat_col <- data.frame(cluster = colnames(xmat), stringsAsFactors = FALSE)
  rownames(mat_col) <- colnames(xmat)
  mat_colors <- list(
    cluster = cluster_colors,
    # chemistry = my_colors$chemistry[c("3p", "5p")],
    # class_short = my_colors$class_short,
    # case = my_colors$case[c("Case", "Control")]
    case = c("Case" = cbPalette[2], "Control" = cbPalette[3])
    # case = my_colors$case[c("Case", "Control")]
  )
  #
  mat_row <- x %>%
    # group_by(donor, case, class_short, chemistry) %>%
    # group_by(donor, case, class_short) %>%
    group_by(donor, case) %>%
    summarize(cells = sum(n)) %>%
    as.data.frame
  rownames(mat_row) <- mat_row$donor
  mat_row$donor <- NULL
  #
  pheatmap::pheatmap(
    mat            = xmat,
    color          = c("white", scico(palette = "batlow", n = 101, direction = -1)),
    border_color   = "grey90",
    cluster_cols   = FALSE,
    cluster_rows   = FALSE,
    annotation_col = mat_col,
    annotation_row = mat_row,
    annotation_colors = mat_colors,
    drop_levels    = TRUE,
    legend_labels  = 10 ^ c(-1, -0.5, 0, 0.5, 1, 1.5),
    labels_row     = str_remove(rownames(xmat), "_[35]p_GEX.*$"),
    fontsize       = 14,
    main           = "Percent of donor's cells across clusters",
    filename       = glue("{out_dir}/heatmap-percent_donor.pdf"),
    width          = ncol(xmat) * 0.3 + 2,
    height         = nrow(xmat) * 0.18 + 2
  )

}

# # Test if the order of the channels in the heatmap is associated with some variable
# sample_test <- sample_info %>%
# filter(channel %in% a2$obs$channel) %>%
# mutate(
#   channel = fct_relevel(channel, levels(xo$channel)),
#   channel_i = as.numeric(channel)
# )
# nrow(sample_test)
# 
# # Exclude columns with too many NA values
# my_cols <- colnames(sample_test)[which(lapply(colnames(sample_test), function(my_col) {
#   sum(is.na(sample_test[[my_col]])) / nrow(sample_test)
# }) < 0.5)]
# 
# # Exclude colinear columns
# my_cols <- setdiff(my_cols, c("sample", "donor", "sample_i"))
# 
# res <- list()
# for (my_col in my_cols) {
#   try({
#     my_formula <- sprintf("channel_i ~ %s", my_col)
#     x <- anova(lm(as.formula(my_formula), sample_test))
#     res[[my_col]] <- broom::tidy(x)[1,]
#   })
# }
# 
# # star_reads_mapped_confidently_to_genome
# do.call(rbind, res) %>% as_tibble %>% arrange(p.value) %>%
# filter(p.value < 0.005)
# 
# # sample_test %>% dplyr::count(dead_cell_depletion_method)

# Abundance of cells per cluster heatmap
########################################################################
my_cols <- c("leiden", "channel", "case", "class_short", "chemistry")
if (all(my_cols %in% colnames(a2$obs)) && n_channels > 1) {
  x <- a2$obs %>%
    dplyr::group_by(leiden, channel, case, class_short, chemistry) %>%
    dplyr::count() %>%
    # dplyr::group_by(leiden) %>%
    # dplyr::mutate(percent_cluster = 100 * n / sum(n)) %>%
    dplyr::group_by(leiden) %>%
    dplyr::mutate(
      percent_cluster = 100 * n / sum(n),
      total_cells   = sum(n)
    )
  #
  set.seed(2)
  x <- seriate_dataframe(x, "leiden", "channel", "percent_cluster")
  x <- as.data.table(x)
  #
  xmat <- dcast(x, channel ~ leiden, value.var = "percent_cluster", fill = 0)
  xmat_rownames <- xmat$channel
  xmat <- as.matrix(xmat[,2:ncol(xmat)])
  rownames(xmat) <- xmat_rownames
  #
  mat_col <- data.frame(cluster = colnames(xmat), stringsAsFactors = FALSE)
  rownames(mat_col) <- colnames(xmat)
  mat_colors <- list(
    cluster = cluster_colors,
    chemistry = my_colors$chemistry[c("3p", "5p")],
    class_short = my_colors$class_short,
    case = my_colors$case[c("Case", "Control")]
  )
  #
  mat_row <- x %>%
    group_by(channel, case, class_short, chemistry) %>%
    summarize(cells = sum(n)) %>%
    as.data.frame
  rownames(mat_row) <- mat_row$channel
  mat_row$channel <- NULL
  #
  pheatmap::pheatmap(
    mat            = xmat,
    color          = scico(palette = "davos", n = 101, direction = -1),
    border_color   = "grey90",
    cluster_cols   = FALSE,
    cluster_rows   = FALSE,
    annotation_col = mat_col,
    annotation_row = mat_row,
    annotation_colors = mat_colors,
    drop_levels    = TRUE,
    labels_row     = str_remove(rownames(xmat), "_[35]p_GEX.*$"),
    fontsize       = 14,
    main           = "Percent of each cluster's cells in each cluster",
    filename       = glue("{out_dir}/heatmap-percent_cluster.pdf"),
    width          = ncol(xmat) * 0.3 + 8,
    height         = nrow(xmat) * 0.3 + 2
  )
}

# Abundance of cells heatmap
########################################################################
my_cols <- c("leiden", "channel", "case", "class_short", "chemistry")
if (all(my_cols %in% colnames(a2$obs)) && n_channels > 1) {
  x <- a2$obs %>%
    dplyr::group_by(leiden, channel, case, class_short, chemistry) %>%
    dplyr::count()
    # dplyr::mutate(n = log10(n + 1))
  #
  set.seed(2)
  x <- seriate_dataframe(x, "leiden", "channel", "n")
  x <- as.data.table(x)
  #
  xmat <- dcast(x, channel ~ leiden, value.var = "n", fill = 0)
  xmat_rownames <- xmat$channel
  xmat <- as.matrix(xmat[,2:ncol(xmat)])
  rownames(xmat) <- xmat_rownames
  #
  mat_col <- data.frame(cluster = colnames(xmat), stringsAsFactors = FALSE)
  rownames(mat_col) <- colnames(xmat)
  mat_colors <- list(
    cluster = cluster_colors,
    chemistry = my_colors$chemistry[c("3p", "5p")],
    class_short = my_colors$class_short,
    case = my_colors$case[c("Case", "Control")]
  )
  #
  mat_row <- x %>%
    group_by(channel, case, class_short, chemistry) %>%
    summarize(cells = sum(n)) %>%
    as.data.frame
  rownames(mat_row) <- mat_row$channel
  mat_row$channel <- NULL
  #
  pheatmap::pheatmap(
    mat            = xmat,
    color          = scico(palette = "davos", n = 101, direction = -1),
    # breaks         = quantile_breaks(xmat, n = 11),
    border_color   = "grey90",
    cluster_cols   = FALSE,
    cluster_rows   = FALSE,
    annotation_col = mat_col,
    annotation_row = mat_row,
    annotation_colors = mat_colors,
    drop_levels    = TRUE,
    labels_row     = str_remove(rownames(xmat), "_[35]p_GEX.*$"),
    fontsize       = 14,
    main           = "Percent of each cluster's cells in each cluster",
    filename       = glue("{out_dir}/heatmap-n.pdf"),
    width          = ncol(xmat) * 0.3 + 8,
    height         = nrow(xmat) * 0.3 + 2
  )
}


# Highly variable gene selection
########################################################################

if ("log2cpm_stats" %in% names(a2)) {
p <- ggplot() +
  stat_binhex(data = a2$log2cpm_stats, aes(mean, sd), bins = 40) +
  geom_line(data = a2$log2cpm_stats %>% dplyr::filter(include), aes(mean, fitted)) +
  geom_point(
    data = a2$log2cpm_stats %>% dplyr::filter(include, rank <= 2000),
    mapping = aes(mean, sd),
    size = 1
  ) +
  scale_fill_gradientn(
    colors = scico::scico(100)[10:100],
    trans = "log10",
    breaks = scales::log_breaks(7)
  ) +
  guides(
    fill = guide_colorbar(barheight = 10)
  ) +
  labs(
    x = bquote("Mean Log"[2]~"CPM"),
    y = bquote("SD Log"[2]~"CPM"),
    title = glue::glue("Select {comma(nrow(a2$log2cpm_stats %>% dplyr::filter(include, rank <= 2000)))} of {comma(nrow(a2$log2cpm_stats))} total genes"),
    fill = "Genes"
  )
my_ggsave(
  "log2cpm_stats",
  # out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1.8, width = 3.5, height = 2, units = "in", dpi = 300
)
}

# Highly variable gene selection
########################################################################

if ("counts_stats" %in% names(a2)) {
  if (nrow(a2$counts_stats) < 1000) {
    p <- ggplot() +
      geom_point(
        data = a2$counts_stats,
        mapping = aes(log10(mean), residuals),
        size = 1
      ) +
      geom_hline(yintercept = 0, size = 0.3) +
      geom_point(
        data = a2$counts_stats[a2$ix_genes,],
        mapping = aes(log10(mean), residuals),
        size = 1.5, color = "red"
      ) +
      # scale_fill_gradientn(
      #   colors = scico::scico(100)[10:100],
      #   trans = "log10",
      #   breaks = scales::log_breaks(7)
      # ) +
      guides(
        fill = guide_colorbar(barheight = 10)
      ) +
      labs(
        x = bquote("Log"[10]~"Mean"),
        y = bquote("Log"[10]~"SD (resid.)"),
        title = glue::glue("Select {comma(length(a2$ix_genes))} of {comma(nrow(a2$counts_stats))} total features"),
        fill = "Features"
      )
  } else {
    p <- ggplot() +
      stat_binhex(
        data = a2$counts_stats,
        mapping = aes(log10(mean), residuals),
        size = 1, bins = 131
      ) +
      geom_hline(yintercept = 0, size = 0.3, linetype = 2) +
      # geom_point(
      #   data = a2$counts_stats[a2$ix_genes,],
      #   mapping = aes(log10(mean), residuals),
      #   size = 0.1, color = "red", alpha = 0.1
      # ) +
      scale_fill_gradientn(
        colors = scico::scico(20)[4:20],
        trans = "log10",
        breaks = scales::log_breaks(7)
      ) +
      guides(
        fill = guide_colorbar(barheight = 10)
      ) +
      labs(
        x = bquote("Log"[10]~"Mean"),
        y = bquote("Log"[10]~"SD (resid.)"),
        title = glue::glue("Select {comma(length(a2$ix_genes))} of {comma(nrow(a2$counts_stats))} total genes"),
        fill = "Genes"
      )
  }
  my_ggsave(
    "counts_stats",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1.8, width = 3.5, height = 2, units = "in", dpi = 300
  )
}

if (n_donors > 1) {
  p1 <- plot_hexmix(
    x = a2$obs$UMAP1,
    y = a2$obs$UMAP2,
    group = a2$obs$donor,
    group_colors = cluster_colors,
    group_labels = FALSE,
    group_legend = TRUE,
    bins = 301
  ) +
  labs(
    title = glue(
      "{comma(nrow(a2$obs))} cells from {length(unique(a2$obs$donor))} donors"
    )
  )
  my_ggsave(
    "umap-donors",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    plot = p1,
    type = "pdf",
    scale = 1, width = 7.5, height = 5, units = "in", dpi = 300
  )
}

leiden_cols <- colnames(a2$obs)[str_detect(colnames(a2$obs), "leiden")]
for (my_leiden in leiden_cols) {
  p1 <- plot_scattermore(
    x = a2$obs$UMAP1,
    y = a2$obs$UMAP2,
    group = a2$obs[[my_leiden]],
    group_colors = cluster_colors,
    pixels = 1000,
    alpha = 0.35
  ) +
  labs(
    title = glue(
      "{length(unique(a2$obs[[my_leiden]]))} clusters of {comma(nrow(a2$obs))} cells from {n_donors} donor{ifelse(n_donors > 1, 's', '')}"
    )
  )
  my_ggsave(
    glue("umap-{my_leiden}"),
    out_dir = out_dir,
    plot = p1,
    type = "pdf",
    scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
  )
}

p1 <- plot_hexmix(
  x = a2$obs$UMAP1,
  y = a2$obs$UMAP2,
  group = a2$obs$leiden,
  group_colors = cluster_colors,
  bins = 301
) +
labs(
  title = glue(
    "{length(unique(a2$obs$leiden))} clusters of {comma(nrow(a2$obs))} cells from {n_donors} donor{ifelse(n_donors > 1, 's', '')}"
  )
)
my_ggsave(
  "umap-clusters-hex",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p1,
  type = "pdf",
  scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
)

plots <- lapply(naturalsort::naturalsort(unique(a2$obs$leiden)), function(this_leiden) {
  ggplot(a2$obs) +
    # stat_summary_hex(bins = 71, color = "grey50", size = 0.1) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = 1,
        group = -2
      ),
      bins = 71, fill = NA, color = "grey50"
    ) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = leiden == this_leiden,
        group = -1
      ),
      bins = 71, color = NA
    ) +
    # guides(fill = guide_colorbar(barheight = 10)) +
    scale_fill_gradientn(
      guide = "none",
      name = NULL,
      limits = c(0.01, 1),
      colors = scico::scico(20),
      na.value = "white",
      trans = "log10"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      title = glue::glue(
        "{this_leiden} (n = {comma(sum(a2$obs$leiden == this_leiden))})"
      )
    )
})
plots[[length(plots)]] <- plots[[length(plots)]] + 
  guides(fill = guide_colorbar(barwidth = 10)) +
  theme(legend.position = "bottom")
ncol <- which.min(abs(sapply(1:16, function(i) {
  i - (16 / 9) * (length(plots) / i)
})))
p <- wrap_plots(plots, ncol = ncol) +
  plot_annotation(
    title = "Proportion of each hexagon occupied by each cluster"
  )
my_ggsave(
  slug = "umap-by-cluster",
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1 + 0.25 * log10(length(unique(a2$obs$leiden))),
  limitsize = FALSE,
  width = 18, height = 9, units = "in", dpi = 300
)

p <- ggplot(a2$obs) +
  aes(UMAP1, UMAP2, z = n_features) +
  stat_summary_hex(bins = 301) +
  guides(fill = guide_colorbar(barheight = 20)) +
  scale_fill_gradientn(
    name = "Features",
    # colors = scico::scico(100)[5:100],
    colors = scico::scico(n = 100, palette = "batlow", direction = -1),
    trans = "log10"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
my_ggsave(
  slug = "umap-n_features",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1, width = 6.5, height = 5, units = "in", dpi = 300
)

if ("mito_pct" %in% colnames(a2$obs)) {
  p <- ggplot(a2$obs) +
    aes(UMAP1, UMAP2, z = mito_pct) +
    stat_summary_hex(bins = 301) +
    guides(fill = guide_colorbar(barheight = 20)) +
    scale_fill_gradientn(
      name = "% Mito.",
      # colors = scico::scico(100)[5:100],
      colors = scico::scico(n = 100, palette = "batlow", direction = -1),
      trans = "log10"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  my_ggsave(
    slug = "umap-percent_mito",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1, width = 6.5, height = 5, units = "in", dpi = 300
  )
}

x0 <- a2$obs %>%
  select(leiden, starts_with("PC")) %>%
  pivot_longer(starts_with("PC"))
x0 <- x0 %>%
  group_by(leiden, name) %>%
  summarize(value = mean(value), .groups = "drop") %>%
  mutate(pc = str_remove(name, "^PC"))
set.seed(42)
x0 <- seriate_dataframe(x0, "pc", "leiden", "value", method = "BEA_TSP")
x0$leiden <- factor(
  x0$leiden,
  (
    x0 %>%
      group_by(leiden) %>%
      summarize(max = max(abs(value))) %>%
      arrange(-max)
  )$leiden
)
color_limit <- max(abs(x0$value)) * c(-1, 1)
p <- ggplot(x0) +
  geom_tile(
    mapping = aes(x = pc, y = leiden, fill = value)
  ) +
  scale_fill_gradientn(
    colors = scico::scico(n = 11, palette = "vik", direction = 1),
    limit = color_limit
  ) +
  guides(fill = guide_colorbar(barheight = 10)) +
  labs(x = "PC", y = "Cluster", fill = "Mean")
my_ggsave(
  slug = "pca-cluster-means",
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1,
  width = length(unique(x0$pc)) * 0.3 + 1.5,
  height = length(unique(x0$leiden)) * 0.25,
  units = "in", dpi = 300
)


cluster_names <- read_excel("data/cluster_names.xlsx")
ix <- cluster_names$Analysis == analysis_name
cluster_names <- cluster_names[ix,]
cluster_names %<>% mutate(
  label = glue("{`Cluster Number`}. {`Short Cluster Name`}")
)
cluster_names <- cluster_names$label
names(cluster_names) <- as.character(1:length(cluster_names))

# my_gene <- "IFNA1"
# my_gene <- "PIGR"
# my_gene <- "ISG15"
my_genes <- c(
  "PTPRC", "PIGR", "EPCAM", "KRT18", "JCHAIN", "IGKC", "FOXP3", "CD14", "CD16",
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD79A", "CD79B", "LYZ", "CD1C", "CLU",
  "MKI67", "CD274", "PDCD1LG2"
)
my_genes <- intersect(my_genes, rowname_key[rownames(a2$counts)])
for (my_gene in my_genes) {
  this_ens <- names(which(rowname_key == my_gene))
  if (my_gene %in% rowname_key && this_ens %in% rownames(a2$counts)) {
    a2$obs$gene <- a2$counts[this_ens,]
    x <- a2$obs %>%
      group_by(leiden) %>%
      summarize(gene_pct = 100 * sum(gene > 0) / length(gene), .groups = "drop")
    x$leiden <- factor(
      x$leiden,
      (
        x %>% arrange(gene_pct)
      )$leiden
    )
    p <- ggplot(x) +
      geom_colh(
        mapping = aes(x = gene_pct, y = leiden)
      ) +
      scale_y_discrete(labels = cluster_names, position = "right") +
      labs(
        x = "Percent",
        y = NULL,
        title = glue("Percent of cells with {my_gene}")
      )
    my_ggsave(
      slug = glue(
        "bars-cluster-percent-{janitor::make_clean_names(my_gene, case = 'none')}"
      ),
      out_dir = glue("{out_dir}/bars"),
      plot = p,
      type = "pdf",
      scale = 1,
      width = 5,
      height = length(x$leiden) * 0.25 + 1,
      units = "in", dpi = 300
    )
  }
}

# my_gene <- "IFNA1"
# my_gene <- "PIGR"
# my_gene <- "ISG15"
my_genes <- c(
  "PTPRC", "PIGR", "EPCAM", "KRT18", "JCHAIN", "IGKC", "FOXP3", "CD14", "CD16",
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD79A", "CD79B", "LYZ", "CD1C", "CLU",
  "MKI67", "CD274", "PDCD1LG2", "IFNA1", "ISG15", "ISG20"
)
my_genes <- intersect(my_genes, rowname_key[rownames(a2$counts)])
for (my_gene in my_genes) {
  this_ens <- names(which(rowname_key == my_gene))
  if (my_gene %in% rowname_key && this_ens %in% rownames(a2$counts)) {
    a2$obs$gene <- a2$counts[this_ens,]
    x <- a2$obs %>%
      group_by(donor, case, leiden) %>%
      summarize(gene_pct = 100 * sum(gene > 0) / length(gene), .groups = "drop")
    x2 <- a2$obs %>%
      group_by(leiden) %>%
      summarize(gene_pct = 100 * sum(gene > 0) / length(gene), .groups = "drop")
    x$leiden <- factor(x$leiden, ( x2 %>% arrange(gene_pct) )$leiden)
    p <- ggplot(x) +
      aes(x = gene_pct, y = leiden, fill = case) +
      ggforestplot::geom_stripes(even = "#ffffff", odd = "#eeeeee") +
      geom_boxplot(linewidth = 0.3, outlier.size = 0.3) +
      geom_point(position = position_quasirandom(dodge.width = 1), size = 0.3) +
      scale_fill_manual(name = NULL, values = pals::okabe(3)[c(2,3)], guide = guide_legend(reverse = TRUE)) +
      scale_y_discrete(labels = cluster_names, position = "right") +
      labs(
        x = "Percent",
        y = NULL,
        title = glue("Percent of cells with <i>{my_gene}</i>")
      ) +
      theme(
        legend.position = "top",
        plot.title = ggtext::element_markdown()
      )
    my_ggsave(
      slug = glue(
        "boxplot-cluster-percent-{janitor::make_clean_names(my_gene, case = 'none')}"
      ),
      out_dir = glue("{out_dir}/boxplot"),
      plot = p,
      type = "pdf",
      scale = 1,
      width = 5,
      height = length(unique(x$leiden)) * 0.25 + 1.5,
      units = "in", dpi = 300
    )
  }
}

# log2cpm boxplot {{{
y <- with(a2$obs, model.matrix(~ 0 + factor(leiden):factor(donor)))
y <- as(y, "dgCMatrix")
# y <- sweep(y, 2, colSums(y), "/")
pb <- as(a2$counts %*% y, "dgCMatrix")
pb <- do_log2cpm(pb, median(colSums(pb)))
#
library(limma)
pb_meta <- str_split_fixed(colnames(pb), ":", 2)
colnames(pb_meta) <- c("leiden", "donor")
pb_meta <- as_tibble(pb_meta)
pb_meta %<>%
  mutate(
    leiden = str_replace(leiden, "factor\\(leiden\\)", ""),
    donor = str_replace(donor, "factor\\(donor\\)", "")
  )
pb_meta <- left_join(
  pb_meta,
  a2$obs %>%
    select(
      donor, case
    ) %>%
    group_by(donor, case) %>%
    summarize_if(is.numeric, mean),
  by = "donor"
)
pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
stopifnot(nrow(pb_meta) == ncol(pb))
colnames(pb) <- str_replace(str_replace(colnames(pb), "factor\\(leiden\\)", ""), "factor\\(donor\\)", "")
stopifnot(all(with(pb_meta, glue("{leiden}:{donor}")) == colnames(pb)))

# pb_de <- fread(glue("results/a20/{analysis_name}/figures/de-case-vs-control/de_case.tsv.gz"))
pb_de <- fread(glue("results/a20/{analysis_name}/figures/de-case-vs-control/de_case-vs-control.tsv.gz"))
my_de <- pb_de %>% filter(Gene %in% c("CD274", "PDCD1LG2")) %>%
  rename(gene = Gene) %>%
  mutate(
    cluster = as.character(cluster),
    signif = ifelse(adj.P.Val < 0.05, "T", "F")
  )
for (my_gene in c("CD274", "PDCD1LG2")) {
  my_ens <- names(which(ensembl_to_symbol == my_gene))
  pb_meta$gene <- pb[my_ens,]
  leiden_inorder <- (
    pb_meta %>% group_by(leiden) %>% summarize(mean = mean(gene)) %>% arrange(mean)
  )$leiden %>% as.character
  pb_meta$leiden <- factor(as.character(pb_meta$leiden), leiden_inorder)
  my_de$cluster <- factor(as.character(my_de$cluster), leiden_inorder)
  p_error <- ggplot(my_de %>% filter(gene == my_gene)) +
    aes(x = logFC, y = cluster) +
    ggforestplot::geom_stripes(even = "#ffffff", odd = "#eeeeee") +
    geom_vline(xintercept = 0, size = 0.3) +
    ggplot2::geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), linewidth = 0.3, height = 0) +
    geom_point(aes(size = signif, fill = signif), stroke = 0.3, shape = 21) +
    scale_size_manual(values = c(2, 3), guide = "none") +
    scale_fill_manual(
      values = c("white", "black"), name = "FDR < 0.05",
      guide = guide_legend(override.aes = list(size = 3))
    ) +
    scale_x_continuous(
      labels = \(x) signif(2^x, 2),
      name = "Fold Change",
    ) +
    labs(
      title = glue("<i>{my_gene}</i>")
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = ggtext::element_markdown(),
      legend.position = "top"
    )
  p_boxplot <- ggplot(pb_meta) +
    aes(x = gene, y = leiden, fill = case) +
    ggforestplot::geom_stripes(even = "#ffffff", odd = "#eeeeee") +
    geom_boxplot(linewidth = 0.3, outlier.size = 0.3) +
    geom_point(position = position_quasirandom(dodge.width = 1), size = 0.3) +
    scale_fill_manual(
      name = NULL,
      values = c("gray50", pals::okabe(3)[2]),
      guide = guide_legend(reverse = TRUE)
    ) +
    scale_y_discrete(labels = cluster_names, position = "right") +
    labs(
      x = bquote("Log"[2]~"CPM"),
      y = NULL
    ) +
    theme(
      legend.position = "bottom"
    )
  my_ggsave(
    glue("boxplot-log2cpm-{my_gene}"),
    #out_dir = glue("figures/{analysis_name}/cluster_top"),
    out_dir = glue("{out_dir}/boxplot"),
    plot = p_error + p_boxplot,
    # plot = p_boxplot,
    type = "pdf",
    scale = 1, units = "in", dpi = 300,
    width = 7,
    height = length(unique(x$leiden)) * 0.25 + 2.5
  )
}

my_gene <- "CD274"
my_ens <- names(which(ensembl_to_symbol == my_gene))
pb_meta$gene <- pb[my_ens,]
pb_meta %>% filter(leiden == 4) %>% arrange(-gene)


# }}}

my_genepairs <- as.list(as.data.frame(combn(c("CD8A", "CD4", "CD14", "CD79A"), 2)))
for (i in seq_along(my_genepairs)) {
  my_gene1 <- my_genepairs[[i]][1]
  my_gene2 <- my_genepairs[[i]][2]
  this_ens1 <- names(which(rowname_key == my_gene1))
  this_ens2 <- names(which(rowname_key == my_gene2))
  if (
    my_gene1 %in% rowname_key && this_ens1 %in% rownames(a2$counts) &&
    my_gene2 %in% rowname_key && this_ens2 %in% rownames(a2$counts)
  ) {
    a2$obs$gene1 <- a2$counts[this_ens1,]
    a2$obs$gene2 <- a2$counts[this_ens2,]
    x <- a2$obs %>%
      group_by(leiden) %>%
      summarize(
        gene1_pct = 100 * sum(gene1 > 0) / length(gene1),
        gene2_pct = 100 * sum(gene2 > 0) / length(gene2),
        .groups = "drop"
      )
    # x$leiden <- factor(
    #   x$leiden,
    #   (
    #     x %>% arrange(gene1_pct)
    #   )$leiden
    # )
    x$leiden <- naturalfactor(x$leiden)
    p <- ggplot(x) +
      geom_point(
        mapping = aes(x = gene1_pct, y = gene2_pct, fill = leiden),
        shape = 21, size = 5, stroke = 0.1
      ) +
      scale_fill_manual(values = mpn65, name = "Cluster") +
      ggrepel::geom_text_repel(
        mapping = aes(x = gene1_pct, y = gene2_pct, label = leiden),
        size = 3
      ) +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      labs(
        x = glue("Percent of cells with *{my_gene1}*"), y = glue("*{my_gene2}*"),
        title = glue("Percent of cells with *{my_gene1}* or <i>{my_gene2}</i>")
      ) +
      theme(
        plot.title = ggtext::element_markdown(),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown()
      )
    my_ggsave(
      slug = glue(
        "{janitor::make_clean_names(my_gene1, case = 'none')}_{janitor::make_clean_names(my_gene2, case = 'none')}"
      ),
      out_dir = glue("{out_dir}/gene-pair"),
      plot = p,
      type = "pdf",
      scale = 1,
      width = 6,
      height = 5,
      units = "in", dpi = 300
    )
  }
}

if ("mito_pct" %in% colnames(a2$obs)) {
# Genes vs percent mitochondrial reads
########################################################################
p <- ggplot(a2$obs) +
  aes(n_features, mito_pct) +
  stat_binhex(bins = 50) +
  # geom_hline(yintercept = 20, size = 0.3) +
  # geom_vline(xintercept = 500, size = 0.3) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_continuous(
    breaks = pretty_breaks(5),
    expan = expansion(mult = c(0.1, 0.1))
  ) +
  # scale_fill_scico(breaks = log_breaks(5), trans = "log10") +
  scale_fill_gradientn(
    colors = scico::scico(20)[3:20],
    breaks = log_breaks(5), trans = "log10"
  ) +
  annotation_logticks(sides = "b", base = 10, color = "grey20", size = 0.4) +
  labs(
    title = "Genes and mito. reads for each cell",
    x = "Number of genes", y = "Percent mito. reads"
  ) +
  guides(fill = guide_colorbar(title = "Cells", barheight = 7))
my_ggsave(
  "genes-vs-mito",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.33, width = 16, height = 9, units = "in", dpi = 300
)
if ("case" %in% colnames(a2$obs)) {
  my_ggsave(
    "genes-vs-mito-by-case",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p + facet_wrap(~ case),
    scale = 0.5, width = 16, height = 9, units = "in", dpi = 300
  )
}
my_ggsave(
  "genes-vs-mito-by-leiden",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p + facet_wrap(~ leiden),
  scale = 0.7, width = 16, height = 14, units = "in", dpi = 300
)
}

# UMIs vs genes
########################################################################
p <- ggplot(a2$obs) +
  aes(x = n_counts, y = n_features) +
  stat_binhex(bins = 50) +
  # geom_hline(yintercept = 500, size = 0.3) +
  # scale_fill_scico(breaks = log_breaks(5), trans = "log10") +
  scale_fill_gradientn(
    colors = scico::scico(20)[3:20],
    breaks = log_breaks(5), trans = "log10"
  ) +
  scale_y_log10(
    labels = label_number_si(), expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_x_log10(
    labels = label_number_si(), expand = expansion(mult = c(0.1, 0.1))
  ) +
  labs(
    title = "UMIs and features for each cell",
    x = "UMIs", y = "Features"
  ) +
  annotation_logticks(base = 10, color = "grey20", size = 0.4) +
  guides(fill = guide_colorbar(title = "Cells", barheight = 7))
my_ggsave(
  "umis-vs-features",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.33, width = 16, height = 9, units = "in", dpi = 300
)
if ("class" %in% colnames(a2$obs)) {
  my_ggsave(
    "umis-vs-features-by-case",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p + facet_wrap(~ class),
    scale = 0.5, width = 16, height = 9, units = "in", dpi = 300
  )
}
my_ggsave(
  "umis-vs-features-by-leiden",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p + facet_wrap(~ leiden),
  scale = 0.7, width = 22, height = 14, units = "in", dpi = 300
)


# Bars by case and class
########################################################################

if (n_donors > 1) {
  my_cols <- c("class", "class2", "donor")
  if (all(my_cols %in% colnames(a2$obs))) {
  # if (all(intersect(my_cols, colnames(a2$obs)) == my_cols)) {
    x <- a2$obs %>%
      dplyr::group_by(class, class2, donor) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(donor = fct_reorder(donor, n))
    p <- ggplot(x) +
      aes(y = donor, x = n, fill = class2) +
      geom_colh() +
      labs(
        x = "Cells", y = "Donor",
        title = glue::glue("Cells per donor (n = {length(unique(a2$obs$donor))})")
      ) +
      scale_x_continuous(labels = scales::label_number_si())
    # if ("class2" %in% names(my_colors)) {
    #   p <- p + scale_fill_manual(name = NULL, values = my_colors$class2)
    # } else {
    #   p <- p + scale_fill_manual(name = NULL, values = pals::okabe())
    # }
    fig_height <- length(unique(a2$obs$donor)) * 0.25
    my_ggsave(
      "n_cells-by-donor",
      #out_dir = glue("figures/{analysis_name}"),
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1, width = 5, height = fig_height, units = "in", dpi = 300
    )
  }
  #
  if (all(c("facs_sorting", "donor", "leiden") %in% colnames(a2$obs))) {
  d <- a2$obs %>%
    dplyr::select(facs_sorting, donor, leiden) %>%
    # dplyr::mutate(facs_sorting = class2) %>%
    dplyr::group_by(facs_sorting, donor, leiden) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(freq = 100 * n / sum(n))
  d_t <- rbindlist(lapply(unique(d$leiden), function(this_leiden) {
    x <- d %>% dplyr::filter(leiden == this_leiden)
    if (length(table(x$facs_sorting)) != 2) {
      return(
        tibble(
          statistic = NA, p.value = 1, method = NA,
          alternative = NA, leiden = this_leiden
        )
      )
    }
    # x <- broom::tidy(t.test(n ~ facs_sorting, x))
    # x <- broom::tidy(t.test(freq ~ facs_sorting, x))
    x <- broom::tidy(wilcox.test(log(freq) ~ facs_sorting, x))
    x$leiden <- this_leiden
    return(x)
  }))
  #
  d_label <- d %>%
    dplyr::group_by(facs_sorting) %>%
    summarise(n = sum(n)) %>%
    dplyr::mutate(label = glue("{facs_sorting} (n = {comma(n)})"))
  d_t <- inner_join(
    d_t,
    d %>%
      dplyr::group_by(leiden) %>%
      dplyr::mutate(y = mean(freq) + 4) %>%
      dplyr::select(leiden, y) %>% unique,
    by = "leiden"
  )
  #
  p <- ggplot() +
    # geom_col(
    #   data = d %>%
    #     dplyr::group_by(leiden, facs_sorting) %>%
    #     dplyr::summarize(freq = median(freq)),
    #   mapping = aes(x = factor(leiden), y = freq, fill = facs_sorting),
    #   position = position_dodge()
    # ) +
    # geom_errorbar(
    geom_segment(
      data = data.frame(
        y = c(10^(-1:2), 0.5 * 10^(-1:2))
      ),
      mapping = aes(
        y = y, yend = y, x = -Inf, xend = Inf
      ),
      size = 0.2, alpha = 0.3
    ) +
    # geom_segment(
    #   data = data.frame(
    #     x = factor(seq_along(unique(d$leiden))),
    #     color = rep(letters[1:2], length.out = length(unique(d$leiden)))
    #   ),
    #   mapping = aes(
    #     y = 0, yend = Inf, x = x, xend = x, color = color
    #   ),
    #   size = length(unique(d$leiden)) * 0.65, alpha = 0.1
    # ) +
    scale_color_manual(
      values = c(NA, "black"), guide = "none"
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by(leiden, facs_sorting) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        x = factor(leiden), ymin = ymin, y = y, ymax = ymax,
        group = facs_sorting, fill = facs_sorting
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    annotate(
      geom = "rect",
      ymin = 0, ymax = Inf,
      xmin = as.integer(sort(unique(d$leiden))) - 0.5,
      xmax = as.integer(sort(unique(d$leiden))) + 0.5,
      fill = rep(
        c(NA, "black"),
        length.out = length(unique(d$leiden))
      ),
      alpha = 0.1
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by(leiden, facs_sorting) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        x = factor(leiden), ymin = ymin, y = y, ymax = ymax,
        group = facs_sorting, fill = facs_sorting
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    geom_point(
      data = d,
      mapping = aes(x = factor(leiden), y = freq, fill = facs_sorting),
      position = position_quasirandom(groupOnX = TRUE, dodge.width = 0.9),
      alpha = 0.4, shape = 21
    ) +
    geom_text(
      data = d_t,
      mapping = aes(
        x = factor(leiden),
        y = 100,
        # y = y + 2,
        label = ifelse(
          # p.value < 0.05 / 12,
          p.value < 0.5 / length(unique(a2$obs$leiden)),
          scientific_10(p.value),
          ""
        )
      ),
      size = 5, parse = TRUE
    ) +
    scale_fill_manual(
      # name = NULL, values = pals::okabe(3)[c(2,3)],
      name = NULL,
      values = pals::okabe(length(unique(a2$obs$facs_sorting))),
      labels = d_label$label
    ) +
    # scale_y_continuous(name = "Percent") +
    scale_y_log10(
      name = "Percent",
      #labels = scales::label_number()
      labels = function(x) signif(x, 3)
    ) +
    annotation_logticks(side = "l") +
    labs(x = NULL) +
    theme(
      # legend.position = c(1, 1),
      # legend.justification = c(1, 1),
      legend.position = "bottom",
      legend.background = element_blank()
      # panel.grid.major.y = element_line(size = 0.1, color = "#00000055"),
      # panel.grid.minor.y = element_line(size = 0.1, color = "#00000055")
    )
  fig_width <- 2 + length(unique(a2$obs$leiden)) * 0.7
  my_ggsave(
    "composition-leiden-sort",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p +
      labs(
        title = "Per-donor percent of cells in each cluster (Wilcoxon Rank Sum P-value)",
        y = NULL
      ),
    scale = 1, width = fig_width, height = 4, units = "in", dpi = 300
  )
  }
}

# Composition analysis
########################################################################

if (n_donors > 1 && all(c("case", "donor", "leiden") %in% colnames(a2$obs))) {
  d <- a2$obs %>%
    dplyr::select(case, donor, leiden) %>%
    # dplyr::mutate(case = class2) %>%
    dplyr::group_by(case, donor, leiden) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(freq = 100 * n / sum(n))
  d_t <- rbindlist(lapply(unique(d$leiden), function(this_leiden) {
    x <- d %>% dplyr::filter(leiden == this_leiden)
    if (length(table(x$case)) != 2) {
      return(
        tibble(
          statistic = NA, p.value = 1, method = NA,
          alternative = NA, leiden = this_leiden
        )
      )
    }
    # x <- broom::tidy(t.test(n ~ case, x))
    # x <- broom::tidy(t.test(freq ~ case, x))
    x <- broom::tidy(wilcox.test(log(freq) ~ case, x))
    x$leiden <- this_leiden
    return(x)
  }))
  d_label <- d %>%
    dplyr::group_by(case) %>%
    summarise(n = sum(n)) %>%
    dplyr::mutate(label = glue("{case} (n = {comma(n)})"))
  d_t <- inner_join(
    d_t,
    d %>% dplyr::group_by(leiden) %>% dplyr::mutate(y = mean(freq) + 4) %>% select(leiden, y) %>% unique,
    by = "leiden"
  )
  scientific_10 <- function(x) {
    ifelse(x < 0.01,
      gsub("e", "%*%10^", scales::scientific_format(digits = 1)(x)),
      signif(x, 1)
    )
  }
  p <- ggplot() +
    # geom_col(
    #   data = d %>%
    #     dplyr::group_by(leiden, case) %>%
    #     dplyr::summarize(freq = median(freq)),
    #   mapping = aes(x = factor(leiden), y = freq, fill = case),
    #   position = position_dodge()
    # ) +
    # geom_errorbar(
    geom_segment(
      data = data.frame(
        y = c(10^(-1:2), 0.5 * 10^(-1:2))
      ),
      mapping = aes(
        y = y, yend = y, x = -Inf, xend = Inf
      ),
      size = 0.2, alpha = 0.3
    ) +
    # geom_segment(
    #   data = data.frame(
    #     x = factor(seq_along(unique(d$leiden))),
    #     color = rep(letters[1:2], length.out = length(unique(d$leiden)))
    #   ),
    #   mapping = aes(
    #     y = 0, yend = Inf, x = x, xend = x, color = color
    #   ),
    #   size = length(unique(d$leiden)) * 0.65, alpha = 0.1
    # ) +
    scale_color_manual(
      values = c(NA, "black"), guide = "none"
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by(leiden, case) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        x = factor(leiden), ymin = ymin, y = y, ymax = ymax,
        group = case, fill = case
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    annotate(
      geom = "rect",
      ymin = 0, ymax = Inf,
      xmin = as.integer(sort(unique(d$leiden))) - 0.5,
      xmax = as.integer(sort(unique(d$leiden))) + 0.5,
      fill = rep(
        c(NA, "black"),
        length.out = length(unique(d$leiden))
      ),
      alpha = 0.1
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by(leiden, case) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        x = factor(leiden), ymin = ymin, y = y, ymax = ymax,
        group = case, fill = case
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    geom_point(
      data = d,
      mapping = aes(x = factor(leiden), y = freq, fill = case),
      position = position_quasirandom(groupOnX = TRUE, dodge.width = 0.9),
      alpha = 0.4, shape = 21
    ) +
    geom_text(
      data = d_t,
      mapping = aes(
        x = factor(leiden),
        y = 100,
        # y = y + 2,
        label = ifelse(
          # p.value < 0.05 / 12,
          p.value < 0.05 / length(unique(a2$obs$leiden)),
          scientific_10(p.value),
          ""
        )
      ),
      size = 5, parse = TRUE
    ) +
    scale_fill_manual(
      name = NULL, values = pals::okabe(3)[c(2,3)],
      labels = d_label$label
    ) +
    # scale_y_continuous(name = "Percent") +
    scale_y_log10(
      name = "Percent",
      #labels = scales::label_number()
      labels = function(x) signif(x, 3)
    ) +
    annotation_logticks(side = "l") +
    labs(x = NULL) +
    theme(
      # legend.position = c(1, 1),
      # legend.justification = c(1, 1),
      legend.position = "bottom",
      legend.background = element_blank()
      # panel.grid.major.y = element_line(size = 0.1, color = "#00000055"),
      # panel.grid.minor.y = element_line(size = 0.1, color = "#00000055")
    )
  fig_width <- 2 + length(unique(a2$obs$leiden)) * 0.7
  my_ggsave(
    "composition-leiden",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p +
      labs(
        title = "Per-donor percent of cells in each cluster (Wilcoxon Rank Sum P-value)",
        y = NULL
      ),
    scale = 1, width = fig_width, height = 4, units = "in", dpi = 300
  )
  dd <- a2$obs %>%
    dplyr::group_by(donor) %>%
    dplyr::summarize(
      total_n = n(),
      n_1 = sum(leiden == 1)
    )
  d <- left_join(d, dd, by = "donor")
  # ddd <- fread("results/a12_4_4_min_genes500_n_pcs20/data/cells.tsv.gz")
  # d <- left_join(
  #   d,
  #   ddd %>% group_by(donor) %>% summarize(n_plasma = sum(leiden == 1) / n()),
  #   by = "donor"
  # )
  p <- ggplot(d) +
    annotation_logticks() +
    aes(x = total_n, y = freq, fill = case) +
    geom_point(shape = 21, size = 3) +
    facet_wrap(~ leiden, ncol = 5) +
    scale_fill_manual(
      name = NULL, values = pals::okabe(3)[c(2,3)]
    ) +
    scale_y_log10() +
    scale_x_log10(labels = scales::label_number_si())
  fig_height <- 2 + ceiling(n_clusters / 5) * 0.7
  my_ggsave(
    "composition-scatter",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p +
      labs(
        title = "Per-donor percent of cells in each cluster",
        y = "Percent",
        x = "Total cells from each donor"
      ),
    scale = 2, width = 6, height = fig_height, units = "in", dpi = 300
  )
  p <- ggplot(d) +
    annotation_logticks() +
    aes(x = n_1, y = freq, fill = case) +
    geom_point(shape = 21, size = 3) +
    facet_wrap(~ leiden, ncol = 5) +
    scale_fill_manual(
      name = NULL, values = pals::okabe(3)[c(2,3)]
    ) +
    scale_y_log10() +
    scale_x_log10(labels = scales::label_number_si())
  fig_height <- 2 + ceiling(n_clusters / 5) * 0.7
  my_ggsave(
    "composition-scatter-c1",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p +
      labs(
        title = "Per-donor percent of cells in each cluster",
        y = "Percent",
        x = "Total cells from cluster 1"
      ),
    scale = 2, width = 6, height = fig_height, units = "in", dpi = 300
  )
}



# this_ens <- names(which(rowname_key == "PTPRC"))
# a2$obs$gene1 <- a2$counts[this_ens,]
# x1 <- a2$obs %>%
#   group_by(leiden) %>%
#   summarize(gene1_pct = 100 * sum(gene1 > 0) / length(gene1), .groups = "drop")
# this_ens <- names(which(rowname_key == "EPCAM"))
# a2$obs$gene2 <- a2$counts[this_ens,]
# x2 <- a2$obs %>%
#   group_by(leiden) %>%
#   summarize(gene2_pct = 100 * sum(gene2 > 0) / length(gene2), .groups = "drop")
# x <- left_join(x1, x2, by = "leiden")
# ggplot(x) +
#   aes(gene1_pct, gene2_pct, label = leiden) +
#   geom_point() +
#   geom_text_repel() +
#   labs(x = "PTPRC", y = "EPCAM")

# ggplot(a2$obs) +
#   aes(PC11, PC27, fill = as.character(leiden)) +
#   geom_point(size = 3, shape = 21, stroke = 0.2) +
#   scale_fill_manual(values = rep(S4Vectors::unname(pals::alphabet()), 2))


# ggplot(x0) +
#   geom_col(
#     mapping = aes(x = name, y = value, fill = as.character(name)),
#     color = "black", size = 0.2
#   ) +
#   facet_wrap(~ leiden, scales = "free_y") +
#   guides(fill = "none") +
#   scale_fill_manual(values = rep(S4Vectors::unname(pals::alphabet(26)), 2)) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )

if (n_donors > 1) {
plots <- lapply(sort(unique(a2$obs$donor)), function(this_donor) {
  ggplot(a2$obs) +
    # stat_summary_hex(bins = 71, color = "grey50", size = 0.1) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = 1,
        group = -2
      ),
      bins = 71, fill = NA, color = "grey50"
    ) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = donor == this_donor,
        group = -1
      ),
      bins = 71, color = NA
    ) +
    # guides(fill = guide_colorbar(barheight = 10)) +
    scale_fill_gradientn(
      guide = "none",
      name = NULL,
      limits = c(0.01, 1),
      colors = scico::scico(20),
      na.value = "white",
      trans = "log10"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      title = glue::glue(
        "{this_donor} (n = {comma(sum(a2$obs$donor == this_donor))})"
      )
    )
})
plots[[length(plots)]] <- plots[[length(plots)]] + 
  guides(fill = guide_colorbar(barwidth = 10)) +
  theme(legend.position = "bottom")
ncol <- which.min(abs(sapply(1:16, function(i) {
  i - (16 / 9) * (length(plots) / i)
})))
p <- wrap_plots(plots, ncol = ncol) +
  plot_annotation(
    title = "Proportion of each hexagon occupied by each donor"
  )
my_ggsave(
  slug = "umap-by-donor",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1 + 0.25 * log10(length(unique(a2$obs$donor))),
  limitsize = FALSE,
  width = 18, height = 9, units = "in", dpi = 300
)
}

if (n_donors > 1 && "class2" %in% colnames(a2$obs)) {
plots <- lapply(sort(unique(a2$obs$class2)), function(this_class2) {
  ggplot(a2$obs) +
    # aes(UMAP1, UMAP2, z = class2 == this_class2, group = -1) +
    # stat_summary_hex(bins = 71, color = "grey50", size = 0.1) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = 1,
        group = -2
      ),
      bins = 71, fill = NA, color = "grey50"
    ) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = class2 == this_class2,
        group = -1
      ),
      bins = 71, color = NA
    ) +
    # guides(fill = guide_colorbar(barheight = 10)) +
    scale_fill_gradientn(
      guide = "none",
      name = NULL,
      limits = c(0.01, 1),
      colors = scico::scico(100),
      na.value = "white",
      trans = "log10"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      title = glue::glue(
        "{this_class2} (n = {comma(sum(a2$obs$class2 == this_class2))})"
      )
    )
})
# plots[[length(plots)]] <- plots[[length(plots)]] + 
#   guides(fill = guide_colorbar(barheight = 10))
ncol <- which.min(abs(sapply(1:16, function(i) {
  i - (16 / 9) * (length(plots) / i)
})))
plots[[ncol]] <- plots[[ncol]] +
  guides(fill = guide_colorbar(barheight = 10))
p <- wrap_plots(plots, ncol = ncol)
my_ggsave(
  slug = "umap-by-class2",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 0.8, width = 16 + 2, height = 9, units = "in", dpi = 300
)
}

if (n_donors > 1 && "batch" %in% colnames(a2$obs)) {
plots <- lapply(sort(unique(a2$obs$batch)), function(this_batch) {
  ggplot(a2$obs) +
    # aes(UMAP1, UMAP2, z = batch == this_batch, group = -1) +
    # stat_summary_hex(bins = 71, color = "grey50", size = 0.1) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = 1,
        group = -2
      ),
      bins = 71, fill = NA, color = "grey50"
    ) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = batch == this_batch,
        group = -1
      ),
      bins = 71, color = NA
    ) +
    # guides(fill = guide_colorbar(barheight = 10)) +
    scale_fill_gradientn(
      guide = "none",
      name = NULL,
      limits = c(0.01, 1),
      colors = scico::scico(100),
      na.value = "white",
      trans = "log10"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      title = glue::glue(
        "{this_batch} (n = {comma(sum(a2$obs$batch == this_batch))})"
      )
    )
})
# plots[[length(plots)]] <- plots[[length(plots)]] + 
#   guides(fill = guide_colorbar(barheight = 10))
ncol <- which.min(abs(sapply(1:16, function(i) {
  i - (16 / 9) * (length(plots) / i)
})))
plots[[ncol]] <- plots[[ncol]] +
  guides(fill = guide_colorbar(barheight = 10))
p <- wrap_plots(plots, ncol = ncol)
my_ggsave(
  slug = "umap-by-batch",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 0.8, width = 16 + 2, height = 9, units = "in", dpi = 300
)
}

if (n_donors > 1 && "facs_sorting" %in% colnames(a2$obs)) {
plots <- lapply(sort(unique(a2$obs$facs_sorting)), function(this_facs_sorting) {
  ggplot(a2$obs) +
    # aes(UMAP1, UMAP2, z = facs_sorting == this_facs_sorting, group = -1) +
    # stat_summary_hex(bins = 71, color = "grey50", size = 0.1) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = 1,
        group = -2
      ),
      bins = 71, fill = NA, color = "grey50"
    ) +
    stat_summary_hex(
      mapping = aes(
        x = UMAP1,
        y = UMAP2,
        z = facs_sorting == this_facs_sorting,
        group = -1
      ),
      bins = 71, color = NA
    ) +
    # guides(fill = guide_colorbar(barheight = 10)) +
    scale_fill_gradientn(
      guide = "none",
      name = NULL,
      limits = c(0.01, 1),
      colors = scico::scico(100),
      na.value = "white",
      trans = "log10"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      title = glue::glue(
        "{this_facs_sorting} (n = {comma(sum(a2$obs$facs_sorting == this_facs_sorting))})"
      )
    )
})
plots[[length(plots)]] <- plots[[length(plots)]] + 
  guides(fill = guide_colorbar(barheight = 10))
ncol <- which.min(abs(sapply(1:16, function(i) {
  i - (16 / 9) * (length(plots) / i)
})))
p <- wrap_plots(plots, ncol = ncol)
my_ggsave(
  slug = "umap-by-facs_sorting",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 0.4, width = 16 + 3, height = 9, units = "in", dpi = 300
)
}


# Number of marker genes per cluster
########################################################################
min_auc <- 0.6
x <- a2$de %>%
  dplyr::filter(auc >= min_auc) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    bin = cut(auc, breaks = seq(min_auc, 1, by = 0.1))
  ) %>%
  select(auc, bin) %>%
  dplyr::count(bin)
# x$bin <- fct_rev(x$bin)
x$group <- factor(x$group, rev(sort(as.integer(unique(x$group)))))
p <- ggplot() +
  annotate(
    geom = "rect",
    xmin = 0, xmax = Inf,
    ymin = as.integer(levels(x$group)) - 0.5,
    ymax = as.integer(levels(x$group)) + 0.5,
    fill = rep(c(NA, "black"), length.out = length(levels(x$group))),
    alpha = 0.1
  ) +
  scale_y_continuous(
    breaks = as.integer(levels(x$group)),
    expand = expansion(c(0.01, 0.0025))
  ) +
  # geom_colh(
  #   data = x,
  #   mapping = aes(
  #     y = as.integer(as.character(group)),
  #     x = log10(n + 1), fill = bin
  #   ),
  #   position = position_dodge2(preserve = "single")
  # ) +
  geom_point(
    data = x,
    mapping = aes(
      y = as.integer(as.character(group)),
      x = n, fill = bin
    ),
    shape = 21, size = 5,
    position = position_dodge2(width = 0.5, preserve = "single")
  ) +
  # geom_point(shape = 21, size = 5) +
  scale_x_log10() +
  # scale_x_continuous(
  #   breaks = log10(c(1, 2, 11, 101)),
  #   labels = c(0, 1, 10, 100)
  # ) +
  annotation_logticks(sides = "b") +
  scale_fill_manual(
    guide = guide_legend(reverse = TRUE),
    name = "AUC",
    # values = rev(scico::scico(5)[1:5])
    values = scico::scico(5)[2:5]
  ) +
  labs(y = NULL, x = "Genes",
    title = "Number of marker genes per cluster")
fig_height <- 1 + n_clusters * 0.3
my_ggsave(
  slug = "marker-genes-per-cluster",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1, width = 5, height = fig_height, units = "in", dpi = 300
)

# x[x$bin == levels(x$bin)[1],]

# my_genes <- unique(
#   (
#     a2$de %>%
#       dplyr::group_by(group) %>%
#       top_n(n = 3, wt = auc * (pct_in - pct_out))
#   )$symbol
# )
# plots <- lapply(my_genes, function(this_gene) {
#   plot_hexgene(
#     x = a2$obs$UMAP1,
#     y = a2$obs$UMAP2,
#     z = as.numeric(log2cpm[names(which(rowname_key == this_gene)),]),
#   ) + labs(title = this_gene) +
#   theme(legend.position = "none")
# })
# p <- wrap_plots(c(list(p1), plots), ncol = 8)
# my_ggsave(
#   slug = "umap-top_marker",
#   #out_dir = glue("figures/{analysis_name}"),
#   out_dir = out_dir,
#   plot = p,
#   type = "pdf",
#   scale = 1.8, width = 16, height = 9, units = "in", dpi = 300
# )


a2$log2cpm_means <- cluster_means(a2$log2cpm, a2$obs$leiden)
x <- reshape2::melt(
  cor(x = as.matrix(a2$log2cpm_means[a2$ix_genes,]), method = "spearman")
)
set.seed(1)
seriate_method <- "BEA_TSP"
if (!all(x$value >= 0)) {
  seriate_method <- "PCA"
}
x <- seriate_dataframe(x, "Var1", "Var2", "value", method = seriate_method)

p <- ggplot(x) +
  aes(Var1, Var2, fill = value) +
  geom_tile() +
  scale_fill_scico(name = "Spearman\nCorrelation") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(barheight = 10)) +
  labs(x = NULL, y = NULL)
fig_width <- 1.5 + n_clusters / 2
fig_height <- n_clusters / 2
my_ggsave(
  slug = "cluster_correlation",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 0.8, width = fig_width, height = fig_height, units = "in", dpi = 300
)

# TCR
########################################################################

if ("has_tcr" %in% colnames(a2$obs)) {
# What percent of each cluster's cells have a TCR?
x <- a2$obs %>% dplyr::count(has_tcr, leiden) %>%
  dplyr::arrange(n)
x$leiden <- factor(x$leiden, (
  x %>% dplyr::group_by(leiden) %>%
    dplyr::summarize(n = sum(n)) %>%
    dplyr::arrange(n)
)$leiden)
p1 <- ggplot(x) +
  aes(y = leiden, x = n, group = has_tcr, fill = has_tcr) +
  geom_colh() +
  scale_fill_manual(
    name = "Has TCR:",
    values = c("grey80", "grey20"),
    labels = c("TRUE" = "T", "FALSE" = "F")
  ) +
  theme(
    # legend.position = c(1, 1), legend.just = c(1, 1),
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number_si()) +
  labs(
    title = "Cells with TCR",
    x = "Cells", y = "Cluster"
  )
x2 <- a2$obs %>%
  dplyr::group_by(leiden) %>%
  dplyr::summarize(
    has_tcr = sum(has_tcr),
    n = n(),
    pct = 100 * sum(has_tcr) / n()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(leiden = factor(leiden))
x2$leiden <- factor(x2$leiden, levels(x$leiden))
p2 <- ggplot(x2) +
  aes(y = leiden, x = pct) +
  geom_colh(fill = "grey20") +
  labs(x = "Percent", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
fig_height <- length(unique(x$leiden)) * 0.28 + 0.5
my_ggsave(
  "tcr-per-cluster",
  out_dir = out_dir,
  type = "pdf",
  plot = p1 + p2,
  scale = 1, width = 6, height = fig_height, units = "in", dpi = 300
)

# What percent of each donor's cells have a TCR?
if (n_channels > 1 && "has_tcr" %in% colnames(a2$obs)) {
  x <- a2$obs %>% dplyr::count(has_tcr, channel) %>%
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
      # legend.position = c(1, 1), legend.just = c(1, 1),
      legend.position = "bottom",
      legend.background = element_blank()
    ) +
    scale_x_continuous(labels = scales::label_number_si()) +
    labs(
      title = "Cells with TCR",
      x = "Cells", y = "Cluster"
    )
  x2 <- a2$obs %>%
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
    "channel-tcr",
    out_dir = file.path(out_dir, "tcr"),
    type = "pdf",
    plot = p1 + p2,
    scale = 1, width = 12, height = fig_height, units = "in", dpi = 300
  )
}

# What percent of each donor's cells have a BCR?
if (n_channels > 1 && "has_bcr" %in% colnames(a2$obs)) {
  x <- a2$obs %>% dplyr::count(has_bcr, channel) %>%
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
      # legend.position = c(1, 1), legend.just = c(1, 1),
      legend.position = "bottom",
      legend.background = element_blank()
    ) +
    scale_x_continuous(labels = scales::label_number_si()) +
    labs(
      title = "Cells with BCR",
      x = "Cells", y = "Cluster"
    )
  x2 <- a2$obs %>%
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
    "channel-bcr",
    out_dir = file.path(out_dir, "bcr"),
    type = "pdf",
    plot = p1 + p2,
    scale = 1, width = 12, height = fig_height, units = "in", dpi = 300
  )
}

p <- plot_hexgene(
  x            = a2$obs$UMAP1,
  y            = a2$obs$UMAP2,
  z            = as.numeric(a2$obs$has_tcr),
  bins         = 301,
  palette      = "davos",
  direction    = -1,
  use_quantile = FALSE,
  legend       = FALSE
) +
guides(
  fill = guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    title = "Proportion of cells with TCR",
    barwidth = 15
  )
)
my_ggsave(
  "umap-tcr",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1, width = 6, height = 6, units = "in", dpi = 300
)

}

#if (n_channels > 1) {
#  my_cols <- c("has_tcr", "leiden", "channel", "case", "class_short", "chemistry")
#  if (all(my_cols %in% colnames(a2$obs))) {
#    x <- a2$obs %>%
#      dplyr::group_by(leiden, channel, case, class_short, chemistry) %>%
#      dplyr::count()
#    x_tcr <- a2$obs %>%
#      dplyr::group_by(leiden, channel, has_tcr) %>%
#      dplyr::summarize(n_tcr = n()) %>%
#      dplyr::filter(has_tcr)
#    x <- left_join(x, x_tcr %>% dplyr::select(-has_tcr), by = c("leiden", "channel"))
#    x$n_tcr[is.na(x$n_tcr)] <- 0
#    x$percent_donor <- 100 * x$n_tcr / x$n
#    #
#    set.seed(2)
#    x <- seriate_dataframe(x, "leiden", "channel", "percent_donor")
#    x <- as.data.table(x)
#    #
#    xmat <- dcast(x, channel ~ leiden, value.var = "percent_donor", fill = 0)
#    xmat_rownames <- xmat$channel
#    xmat <- as.matrix(xmat[,2:ncol(xmat)])
#    rownames(xmat) <- xmat_rownames
#    #
#    mat_col <- data.frame(cluster = colnames(xmat), stringsAsFactors = FALSE)
#    rownames(mat_col) <- colnames(xmat)
#    mat_colors <- list(
#      cluster = cluster_colors,
#      chemistry = my_colors$chemistry[c("3p", "5p")],
#      class_short = my_colors$class_short,
#      case = my_colors$case[c("Case", "Control")]
#    )
#    #
#    mat_row <- x %>%
#      group_by(channel, case, class_short, chemistry) %>%
#      summarize(cells = sum(n)) %>%
#      as.data.frame
#    rownames(mat_row) <- mat_row$channel
#    mat_row$channel <- NULL
#    if (var(as.numeric(xmat)) > 0) {
#      pheatmap::pheatmap(
#        mat            = xmat,
#        color          = scico(palette = "davos", n = 21, direction = -1),
#        border_color   = "grey90",
#        cluster_cols   = FALSE,
#        cluster_rows   = FALSE,
#        annotation_col = mat_col,
#        annotation_row = mat_row,
#        annotation_colors = mat_colors,
#        drop_levels    = TRUE,
#        labels_row     = str_remove(rownames(xmat), "_[35]p_GEX.*$"),
#        fontsize       = 14,
#        main           = "Percent of cells with TCR, grouped by channel and cluster",
#        filename       = glue("{out_dir}/tcr/heatmap-percent.pdf"),
#        width          = ncol(xmat) * 0.3 + 8,
#        height         = nrow(xmat) * 0.3 + 2
#      )
#    }
#  }
#}

# BCR
########################################################################

if ("has_bcr" %in% colnames(a2$obs)) {
x <- data.frame(with(a2$obs, table(has_bcr, leiden)))
x$leiden <- fct_rev(x$leiden)
p1 <- ggplot(x) +
  aes(y = leiden, x = Freq, group = has_bcr, fill = has_bcr) +
  geom_colh() +
  scale_fill_manual(
    name = "Has BCR:",
    values = c("grey80", "grey20"),
    labels = c("TRUE" = "T", "FALSE" = "F")
  ) +
  theme(
    # legend.position = c(1, 1), legend.just = c(1, 1),
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number_si()) +
  labs(
    title = "Cells with BCR",
    x = "Cells", y = "Cluster"
  )
x <- a2$obs %>%
  group_by(leiden) %>%
  summarize(
    has_bcr = sum(has_bcr),
    n = n(),
    pct = 100 * sum(has_bcr) / n()
  ) %>%
  ungroup() %>%
  mutate(leiden = factor(leiden))
x$leiden <- fct_rev(x$leiden)
p2 <- ggplot(x) +
  aes(y = leiden, x = pct) +
  geom_colh(fill = "grey20") +
  labs(x = "Percent", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
fig_height <- length(unique(x$leiden)) * 0.28 + 0.5
my_ggsave(
  "bcr-per-cluster",
  out_dir = out_dir,
  type = "pdf",
  plot = p1 + p2,
  scale = 1, width = 6, height = fig_height, units = "in", dpi = 300
)
#
p <- plot_hexgene(
  x            = a2$obs$UMAP1,
  y            = a2$obs$UMAP2,
  z            = as.numeric(a2$obs$has_bcr),
  bins         = 201,
  palette      = "davos",
  direction    = -1,
  use_quantile = FALSE,
  legend       = FALSE
) +
guides(
  fill = guide_colorbar(
    direction = "horizontal",
    title.position = "top",
    title = "Proportion of cells with BCR",
    barwidth = 15
  )
)
my_ggsave(
  "umap-bcr",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1, width = 6, height = 6, units = "in", dpi = 300
)
}

if ("has_tcr" %in% colnames(a2$obs) && sum(a2$obs$has_tcr) / nrow(a2$obs) > 0.25) {
tcr_clusters <- (
  a2$obs %>%
  dplyr::group_by(leiden) %>%
  dplyr::summarize(
    has_tcr = sum(has_tcr),
    n = n(),
    pct = 100 * sum(has_tcr) / n()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(pct > 10)
)$leiden
if (length(tcr_clusters)) {
x <- a2$obs %>%
  dplyr::filter(!is.na(TRA_cdr3)) %>%
  dplyr::filter(leiden %in% tcr_clusters) %>%
  dplyr::group_by(leiden) %>%
  dplyr::count(
    TRBV,
    TRAV,
    TRB_cdr3_trim,
    TRA_cdr3_trim
  ) %>%
  dplyr::arrange(-n) %>%
  dplyr::mutate(rank = seq(length(n)), percent = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(leiden = factor(leiden)) %>%
  dplyr::mutate(
    tcr = paste(sep = "|", TRAV, TRBV, TRA_cdr3_trim, TRB_cdr3_trim)
  )
  top_tcrs <- unique((
    x %>%
      dplyr::filter(rank < 10)
  )$tcr)
  if (length(top_tcrs) > 1) {
    x <- x %>%
      dplyr::filter(tcr %in% top_tcrs) %>%
      dplyr::select(leiden, tcr, percent)
    xw <- pivot_wider(x, names_from = "tcr", values_from = "percent")
    mat <- as.matrix(xw[,2:ncol(xw)])
    rownames(mat) <- xw$leiden
    mat[is.na(mat)] <- 0
    mat <- mat[,apply(mat, 2, max) >= 1]
    # Order the heatmap
    set.seed(1)
    xo <- seriation::seriate(mat, method = "BEA_TSP")
    x$leiden <- factor(x$leiden, rownames(mat)[xo[[1]]])
    x$tcr <- factor(x$tcr, colnames(mat)[xo[[2]]])
    # Plot
    p2 <- ggplot(x) +
    geom_tile(
      aes(y = leiden, x = tcr, fill = percent)
    ) +
    scale_x_discrete(position = "b", name = "TCR", expand = c(0, 0)) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    # scale_fill_viridis_c(
    #   name = "AUC", guide = guide_colorbar(barheight = 40), breaks = pretty_breaks(10),
    #   limits = c(0, 1)
    # ) +
    scale_fill_scico(
      palette = "davos",
      direction = -1,
      name = "Percent",
      guide = guide_colorbar(barheight = length(tcr_clusters)),
      breaks = pretty_breaks(5)
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  p1 <- ggplot(x) +
  geom_tile(
    aes(x = 1, y = leiden, fill = leiden)
  ) +
  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
  scale_fill_manual(
    values = cluster_colors, guide = "none"
  )
  p <- (
    p1 + theme(plot.margin = margin(b = 0))
  ) + (
    p2 + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0)
    )
  ) + plot_layout(widths = c(1, ncol(mat))) +
  plot_annotation(title = "Abundance of top 10 TCRs for each cell cluster")
  fig_width <- ncol(mat) * 0.25
  fig_height <- nrow(mat) * 0.35
  my_ggsave(
    slug = "heatmap-cluster-top10-tcr",
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )
}
}

# bcr_clusters <- (
#   a2$obs %>%
#   dplyr::group_by(leiden) %>%
#   dplyr::summarize(
#     has_bcr = sum(has_bcr),
#     n = n(),
#     pct = 100 * sum(has_bcr) / n()
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(pct > 10)
# )$leiden
# if (length(bcr_clusters)) {
# x <- a2$obs %>%
#   dplyr::filter(!is.na(TRA_cdr3)) %>%
#   dplyr::filter(leiden %in% bcr_clusters) %>%
#   dplyr::group_by(leiden) %>%
#   dplyr::count(
#     IGLC,
#     IGLV
#   ) %>%
#   dplyr::arrange(-n) %>%
#   dplyr::mutate(rank = seq(length(n)), percent = 100 * n / sum(n)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(leiden = factor(leiden)) %>%
#   dplyr::mutate(
#     bcr = paste(sep = "|", IGLC, IGLV)
#   )
# top_bcrs <- unique((
#   x %>%
#     dplyr::filter(rank < 10)
# )$bcr)
# if (length(top_bcrs) > 1) {
#   x <- x %>%
#     dplyr::filter(bcr %in% top_bcrs) %>%
#     dplyr::select(leiden, bcr, percent)
#   xw <- pivot_wider(x, names_from = "bcr", values_from = "percent")
#   mat <- as.matrix(xw[,2:ncol(xw)])
#   rownames(mat) <- xw$leiden
#   mat[is.na(mat)] <- 0
#   mat <- mat[,apply(mat, 2, max) >= 1]
#   # Order the heatmap
#   set.seed(1)
#   xo <- seriation::seriate(mat, method = "BEA_TSP")
#   x$leiden <- factor(x$leiden, rownames(mat)[xo[[1]]])
#   x$bcr <- factor(x$bcr, colnames(mat)[xo[[2]]])
#   # Plot
#   p2 <- ggplot(x) +
#   geom_tile(
#     aes(y = leiden, x = bcr, fill = percent)
#   ) +
#   scale_x_discrete(position = "b", name = "BCR", expand = c(0, 0)) +
#   scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
#   # scale_fill_viridis_c(
#   #   name = "AUC", guide = guide_colorbar(barheight = 40), breaks = pretty_breaks(10),
#   #   limits = c(0, 1)
#   # ) +
#   scale_fill_scico(
#     palette = "davos",
#     direction = -1,
#     name = "Percent",
#     guide = guide_colorbar(barheight = length(bcr_clusters)),
#     breaks = pretty_breaks(5)
#   ) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )
# }
# p1 <- ggplot(x) +
# geom_tile(
#   aes(x = 1, y = leiden, fill = leiden)
# ) +
# scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
# scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
# scale_fill_manual(
#   values = cluster_colors, guide = "none"
# )
# p <- (
#   p1 + theme(plot.margin = margin(b = 0))
# ) + (
#   p2 + theme(
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin = margin(t = 0)
#   )
# ) + plot_layout(widths = c(1, ncol(mat))) +
# plot_annotation(title = "Abundance of top 10 BCRs for each cell cluster")
# fig_width <- ncol(mat) * 0.25 + 1
# fig_height <- nrow(mat) * 0.35 + 1
# my_ggsave(
#   slug = "heatmap-cluster-top10-bcr",
#   out_dir = out_dir,
#   type = "pdf",
#   plot = p,
#   scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
# )
# }

# TODO

if ("TRAV" %in% colnames(a2$obs)) {
a2$obs$TRAV <- str_replace(a2$obs$TRAV, "/", "")
this_gene <- "TRAV29DV5"
this_ens  <- names(which(ensembl_to_symbol == this_gene))
ensembl_to_symbol[which(str_detect(ensembl_to_symbol, "TRAV29"))]
a2$obs$x <- as.numeric(a2$log2cpm[this_ens,])

my_genes <- ensembl_to_symbol[str_detect(ensembl_to_symbol, "^TRAV")]
my_ens <- intersect(names(my_genes), rownames(a2$log2cpm))
my_genes <- S4Vectors::unname(ensembl_to_symbol[my_ens])
ix <- apply(a2$counts[my_ens,], 2, which.max)
a2$obs$TRAV_gex <- my_genes[ix]
a2$obs$TRAV_gex[colSums(a2$counts[my_ens,]) == 0] <- ""
# table(a2$obs$TRAV_gex)["TRAV2"]
# Count cells and how many times gene expression matches TCR
x <- a2$obs %>%
  dplyr::filter(TRAV_gex != "") %>%
  dplyr::group_by(TRAV_gex, TRAV == TRAV_gex & !is.na(TRAV)) %>%
  dplyr::summarize(cells = n())  %>%
  dplyr::ungroup()
colnames(x) <- c("gene", "match", "cells")
# Sort genes
x$gene <- factor(x$gene, (
  x %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(n = sum(cells)) %>%
  dplyr::arrange(n)
)$gene)
# Bar plot
p1 <- ggplot(x) +
  aes(y = gene, x = cells, fill = match) +
  geom_colh() +
  scale_fill_manual(
    name = "TCR GEX matches SEQ",
    values = c("grey80", "grey20"),
    labels = c("TRUE" = "T", "FALSE" = "F")
  ) +
  theme(
    # legend.position = c(1, 1), legend.just = c(1, 1),
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number_si()) +
  labs(
    title = "Cells with TRAV gene expression",
    x = "Cells", y = NULL
  )
# Percent matched
x2 <- x %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(total = sum(cells)) %>%
  dplyr::filter(match) %>%
  dplyr::mutate(pct = 100 * cells / total)
x2$gene <- factor(x2$gene, levels(x$gene))
x2 <- left_join(x %>% select(gene), x2) %>% unique
# Bar plot
p2 <- ggplot(x2) +
  aes(y = gene, x = pct) +
  geom_colh(fill = "grey20") +
  theme(
    # legend.position = c(1, 1), legend.just = c(1, 1),
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number_si()) +
  labs(
    title = "Percent matching TCR seq",
    x = "Percent", y = NULL
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
fig_height <- 0.5 * length(unique(x$gene)) * 0.3
my_ggsave(
  slug = "tcr-TRAV-gex-vs-seq",
  out_dir = out_dir,
  plot = p1 + p2,
  type = "pdf",
  scale = 1.8, width = 6, height = fig_height, units = "in", dpi = 300
)
}


if ("IGLV" %in% colnames(a2$obs)) {
my_genes <- ensembl_to_symbol[str_detect(ensembl_to_symbol, "^IGLV")]
my_ens <- intersect(names(my_genes), rownames(a2$log2cpm))
my_genes <- S4Vectors::unname(ensembl_to_symbol[my_ens])
ix <- apply(a2$counts[my_ens,], 2, which.max)
a2$obs$IGLV_gex <- my_genes[ix]
a2$obs$IGLV_gex[colSums(a2$counts[my_ens,]) == 0] <- ""
# table(a2$obs$IGLV_gex)["IGLV2"]
# Count cells and how many times gene expression matches BCR
x <- a2$obs %>%
  dplyr::filter(IGLV_gex != "") %>%
  dplyr::group_by(IGLV_gex, IGLV == IGLV_gex & !is.na(IGLV)) %>%
  dplyr::summarize(cells = n())  %>%
  dplyr::ungroup()
colnames(x) <- c("gene", "match", "cells")
# Sort genes
x$gene <- factor(x$gene, (
  x %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(n = sum(cells)) %>%
  dplyr::arrange(n)
)$gene)
# Bar plot
p1 <- ggplot(x) +
  aes(y = gene, x = cells, fill = match) +
  geom_colh() +
  scale_fill_manual(
    name = "BCR GEX matches SEQ",
    values = c("grey80", "grey20"),
    labels = c("TRUE" = "T", "FALSE" = "F")
  ) +
  theme(
    # legend.position = c(1, 1), legend.just = c(1, 1),
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number_si()) +
  labs(
    title = "Cells with IGLV gene expression",
    x = "Cells", y = NULL
  )
# Percent matched
x2 <- x %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(total = sum(cells)) %>%
  dplyr::filter(match) %>%
  dplyr::mutate(pct = 100 * cells / total)
x2$gene <- factor(x2$gene, levels(x$gene))
x2 <- left_join(x %>% select(gene), x2) %>% unique
# Bar plot
p2 <- ggplot(x2) +
  aes(y = gene, x = pct) +
  geom_colh(fill = "grey20") +
  theme(
    # legend.position = c(1, 1), legend.just = c(1, 1),
    legend.position = "bottom",
    legend.background = element_blank()
  ) +
  scale_x_continuous(labels = scales::label_number_si()) +
  labs(
    title = "Percent matching BCR seq",
    x = "Percent", y = NULL
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
fig_height <- 0.5 * length(unique(x$gene)) * 0.3
my_ggsave(
  slug = "bcr-IGLV-gex-vs-seq",
  out_dir = out_dir,
  plot = p1 + p2,
  type = "pdf",
  scale = 1.8, width = 6, height = fig_height, units = "in", dpi = 300
)
}

# x <- a2$obs %>%
#   dplyr::filter(!is.na(TRA_cdr3)) %>%
#   dplyr::group_by(leiden, case, donor) %>%
#   dplyr::count(
#     TRBV,
#     TRAV,
#     TRB_cdr3_trim,
#     TRA_cdr3_trim
#   ) %>%
#   dplyr::arrange(-n) %>%
#   dplyr::mutate(rank = seq(length(n)), percent = 100 * n / sum(n))
# p1 <- ggplot(x) +
#   aes(x = rank, y = percent, group = donor, color = case) +
#   geom_line(size = 1) +
#   scale_color_manual(values = my_colors$case) +
#   scale_x_log10() +
#   scale_y_log10() +
#   facet_wrap(~ leiden) +
#   annotation_logticks(sides = "bl") +
#   labs(x = "Rank", y = "Percent")
# p1

# x <- a2$obs %>%
#   group_by(chemistry) %>%
#   summarize(
#     has_tcr = sum(has_tcr),
#     n = n(),
#     pct = 100 * sum(has_tcr) / n()
#   ) %>%
#   ungroup()
# x
# x <- a2$obs %>%
#   group_by(chemistry) %>%
#   summarize(
#     has_bcr = sum(has_bcr),
#     n = n(),
#     pct = 100 * sum(has_bcr) / n()
#   ) %>%
#   ungroup()
# x



# Hierarchical dendrogram with clusters
########################################################################

# a2$cluster_dist <- as.dist(
#   1 - cor(
#     as.matrix(a2$log2cpm_means[a2$ix_genes,]),
#     method = "spearman"
#   )
# )
# hc <- dendsort::dendsort(hclust(a2$cluster_dist, method = "complete"))
# # plot(as.dendrogram(hc), horiz = TRUE)
# d <- tidygraph::as_tbl_graph(hc)
# ggraph(d, layout = 'dendrogram', circular = FALSE) + 
#   geom_edge_elbow() +
#   geom_node_point(
#     mapping = aes(
#       filter = leaf,
#       x      = x,
#       y      = y - 0.2,
#       colour = label
#     ),
#     size = 8
#   ) +
#   geom_node_text(
#     mapping = aes(
#       x = x - 0.01, y = y - 0.2, filter = leaf, label = label
#     ),
#     size = 5, fontface = "bold", alpha = 0.5
#   ) +
#   geom_node_text(
#     mapping = aes(
#       x = x + 0.01, y = y - 0.2, filter = leaf, label = label
#     ),
#     size = 5, fontface = "bold", alpha = 0.5
#   ) +
#   geom_node_text(
#     mapping = aes(
#       x = x, y = y - 0.2 - 0.02, filter = leaf, label = label
#     ),
#     size = 5, fontface = "bold", alpha = 0.5
#   ) +
#   geom_node_text(
#     mapping = aes(
#       x = x, y = y - 0.2 + 0.02, filter = leaf, label = label
#     ),
#     size = 5, fontface = "bold", alpha = 0.5
#   ) +
#   geom_node_text(
#     mapping = aes(
#       x = x, y = y - 0.2, filter = leaf, label = label
#     ),
#     size = 5, color = "white"
#   ) +
#   scale_colour_manual(values = cluster_colors) +
#   theme_void() +
#   theme(
#     legend.position = "none",
#     plot.margin     = unit(c(0,0,0,0),"cm"),
#   ) +
#   expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))


# UMAP with top markers
########################################################################
p1 <- plot_hexmix(
  x = a2$obs$UMAP1,
  y = a2$obs$UMAP2,
  group = a2$obs$leiden,
  group_colors = cluster_colors,
  bins = 301
) +
labs(
  title = glue(
    "{length(unique(a2$obs$leiden))} clusters of {comma(nrow(a2$obs))} cells from {length(unique(a2$obs$donor))} donors"
  )
)

##
#a2$de_top <- a2$de %>%
#  dplyr::group_by(group) %>%
#  dplyr::filter(pct_in > 10) %>%
#  dplyr::mutate(
#    rank1 = rank(100 * auc - pct_out),
#    rank2 = rank(logFC * (pct_in - pct_out))
#  ) %>%
#  # top_n(n = 3, wt = (rank1 + rank2))
#  # dplyr::top_n(n = 3, wt = 100 * auc - pct_out) %>%
#  dplyr::top_n(n = 10, wt = auc) %>%
#  dplyr::arrange(as.integer(group))
#  # top_n(n = 2, wt = logFC * (pct_in - pct_out))
#a2$de_top %>%
#  dplyr::select(
#    feature, symbol, group, rank1, rank2, logFC, auc, pct_in, pct_out
#  ) %>%
#  arrange(-auc)
#  # dplyr::filter(str_detect(symbol, "^AL"))
#a2$de_top <- a2$de_top[!duplicated(a2$de_top$feature),]
#a2$de_top <- a2$de_top %>%
#  dplyr::group_by(group) %>%
#  dplyr::top_n(n = 3, wt = auc) %>%
#  dplyr::arrange(as.integer(group), -auc)
#table(a2$de_top$group)
##

# Similar to the heatmap
a2$de_top <- a2$de %>%
  dplyr::group_by(group) %>%
  dplyr::filter(pct_in > 10) %>%
  dplyr::mutate(
    rank1 = rank(100 * auc - pct_out),
    rank2 = rank(logFC * (pct_in - pct_out))
  ) %>%
  dplyr::top_n(n = 3, wt = (rank1 + rank2))
these_genes <- unique(a2$de_top$feature)
plots <- lapply(seq(nrow(a2$de_top)), function(i) {
  this_gene <- a2$de_top$feature[i]
  this_cluster <- a2$de_top$group[i]
  plot_hexgene(
    x            = a2$obs$UMAP1,
    y            = a2$obs$UMAP2,
    z            = as.numeric(a2$log2cpm[this_gene,]),
    bins         = 101,
    palette      = "davos",
    direction    = -1,
    use_quantile = FALSE
  ) +
  labs(
    title = glue("{rowname_key[this_gene]} ({this_cluster})")
  ) +
  theme(legend.position = "none")
})
ncol <- which.min(abs(sapply(1:16, function(i) {
  i - (16 / 9) * (length(plots) / i)
})))
p <- wrap_plots(c(list(p1), plots), ncol = ncol)
my_ggsave(
  slug = "umap-top_marker-auc",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  type = "pdf",
  scale = 1.8, width = 16, height = 9, units = "in", dpi = 300
)


# Heatmap of known a priori markers for each cluster
########################################################################
#
# monaco <- janitor::clean_names(readxl::read_excel(
#   path = "data/Monaco2019/Monaco2019_TableS2.xlsx",
#   sheet = "FDR TPM",
#   skip = 2
# ))
# monaco_types <- colnames(monaco)[4:length(colnames(monaco))]
# monaco_markers <- list()
# for (x in monaco_types) {
#   r <- rank(monaco[[x]])
#   monaco_markers[[x]] <- stringr::str_split_fixed(
#     monaco$ensembl_id[which(r <= 3)], "\\.", 2
#   )[,1]
# }
# these_genes <- unique(unlist(monaco_markers))
# 
known_markers <- list(
  "Immune"        = c("PTPRC"),
  "Pre-T"         = c("CD3E", "SH2D1A"),
  "CD4 Mem"       = c("CD4", "CCR7"),
  "CD4 Naive"     = c("SELL"),
  "TCR"           = c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2"),
  "Treg"          = c("IL2RA"),
  "HLA"           = c("HLA-DRA"),
  #
  "Cycling"       = c("STMN1", "TYMS", "MKI67"),
  #
  "CD8 Naive"     = c("CD8A"),
  "CD8 Effector"  = c("CCL5"),
  "NK"            = c("NKG7", "KLRB1"),
  "MAIT"          = c("TRAV1-2"),
  #
  "Naive B"       = c("CD79A", "CD79B", "MS4A1", "CD19"),
  "pre-B"         = c("VPREB3", "MME"),
  "Memory B"      = c("CD27"),
  "Plasma"        = c("XBP1", "JCHAIN"),
  #
  "DC"            = c("FCER1A"),
  "pDC"           = c("DERL3"),
  "GMP"           = c("MPO"),
  "CD14 Monocyte" = c("CD14", "LYZ"),
  "CD16 Monocyte" = c("FCG3RA"),
  "PPBP Monocyte" = c("PPBP"),
  "Mk"            = c("PF4"),
  "Erythr"        = c("HBB"),
  "HSC/MPP"       = c("CD34", "AVP")
)
invert_list <- function(lov) split(rep(names(lov), lengths(lov)), unlist(lov))
gene_to_name <- unlist(invert_list(known_markers))

# known_markers <- c(known_markers, list(
#   "treg1" = c("IL12RB2", "TNFRSF4", "CXCR3", "TNFRSF18", "STAT1"),
#   "treg2" = c("KLF2", "LEF1", "ITGB1", "ZFP36L2"),
#   "treg3" = c("VIM", "LGALS1", "GPR15", "TXNIP", "FOSB"),
#   "treg4" = c("IL10", "KLRB1", "LAG3", "IKZF3", "GPR25"),
#   "treg5" = c("SLC2A3", "CCR4", "NEAT1", "CD38", "CCR7")
# ))

# known_markers <- c(known_markers, list(
#   "custom" = strsplit(
#   "CD3E TRAC TRBC1 TRBC2 TRDC TRGC1 TRGC2 CD4 CD8A CD8B CD79A CD79B NCAM1 HLA-DRA LYZ CST3", " ")[[1]]
# ))

these_genes <- sort(unique(unlist(known_markers)))
these_genes <- names(rowname_key[rowname_key %in% these_genes])
these_genes <- intersect(these_genes, a2$de$feature)
these_genes <- these_genes[order(ensembl_to_symbol[these_genes])]
fig_height <- length(these_genes) * 0.3 + 1
fig_width <- length(unique(a2$obs$leiden)) * 0.4 + 2
#
# Subset the table
x <- a2$de %>% dplyr::filter(feature %in% these_genes)
x <- as.data.frame(dcast.data.table(
  data = as.data.table(x), formula = symbol ~ group, value.var = "auc"
))
rownames(x) <- x$symbol
x$symbol <- NULL
# Order the heatmap
set.seed(1)
xo <- seriation::seriate(as.matrix(x), method = "BEA_TSP")
xd <- a2$de %>% dplyr::filter(feature %in% these_genes)
xd$symbol <- factor(xd$symbol, rownames(x)[xo[[1]]])
xd$group <- factor(xd$group, colnames(x)[xo[[2]]])
xd$symbol_group <- factor(gene_to_name[as.character(xd$symbol)], names(known_markers))
# Plot
p2 <- ggplot(xd) +
geom_tile(
  aes(y = symbol, x = group, fill = auc)
) +
facet_grid(rows = "symbol_group", scales = "free_y", space = "free_y") +
# facet_wrap(~ symbol_group, scales = "free_y", ncol = 1) +
scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
# scale_fill_viridis_c(
#   name = "AUC", guide = guide_colorbar(barheight = 40), breaks = pretty_breaks(10),
#   limits = c(0, 1)
# ) +
# scale_fill_scico(
#   palette = "vik", limits = c(0, 1), direction = 1,
#   name = "AUC", guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(9)
# ) +
scale_fill_gradientn(
  colors = rev(RColorBrewer::brewer.pal(n = 11, "RdBu")),
  limits = c(0, 1),
  name = "AUC", guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(9)
) +
theme(
  strip.text.y = element_text(angle = 0, hjust = 0, size = 12),
  panel.spacing.y = unit(0.33, "lines"),
  axis.text.y = element_text(size = 12, face = "italic"),
#   axis.text.y = element_blank(), axis.ticks.y = element_blank()
) +
labs(x = NULL, y = NULL)
p1 <- ggplot(xd) +
geom_tile(
  aes(y = 1, x = group, fill = group)
) +
scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
scale_fill_manual(
  values = cluster_colors, guide = "none"
)
p <- (
  p1 + theme(plot.margin = margin(b = 0), axis.text.x = element_text(size = 12))
) / (
  p2 + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
) + plot_layout(
  heights = c(0.8, length(these_genes))
) + plot_annotation(
  title = glue("{length(unique(xd$symbol))} selected genes")
)
#p
my_ggsave(
  slug = "heatmap-known-markers",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
)


# my_symbol <- "TRAV1-2"
# this_gene <- names(ensembl_to_symbol)[ensembl_to_symbol == my_symbol]


# UMAP with known markers
########################################################################
plots <- lapply(seq_along(these_genes), function(i) {
  this_gene    <- these_genes[i]
  plot_hexgene(
    x            = a2$obs$UMAP1,
    y            = a2$obs$UMAP2,
    z            = as.numeric(a2$log2cpm[this_gene,]),
    bins         = 101,
    palette      = "davos",
    direction    = -1,
    use_quantile = FALSE
  ) +
  labs(
    title = rowname_key[this_gene]
  ) +
  theme(legend.position = "none")
})
#
x <- wrap_plots(plots, ncol = 5)
my_ggsave(
  slug = "umap-known-markers",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = x,
  # plot = wrap_plots(p[1:9], ncol = 1),
  type = "pdf",
  limitsize = FALSE,
  scale = 1.8,
  width = 5 * 1.8,
  height = 1.5 * ceiling(length(these_genes) / 5),
  units = "in", dpi = 300
)

# this_gene <- names(rowname_key[rowname_key == "FCER1A"])
#   plot_hexgene(
#     x            = a2$obs$UMAP1,
#     y            = a2$obs$UMAP2,
#     z            = as.numeric(a2$log2cpm[this_gene,]),
#     bins         = 101,
#     palette      = "davos",
#     direction    = -1,
#     use_quantile = TRUE
#   )

# Heatmap of best markers for each cluster
########################################################################
# Top markers for each cluster
a2$de_top <- a2$de %>%
  dplyr::group_by(group) %>%
  dplyr::filter(pct_in > 10) %>%
  # dplyr::filter(group %in% c(12, 14, 16, 17, 19, 25, 29)) %>%
  dplyr::mutate(
    rank1 = rank(100 * auc - pct_out),
    rank2 = rank(logFC * (pct_in - pct_out))
  ) %>%
  dplyr::top_n(n = 3, wt = (rank1 + rank2))
these_genes <- unique(a2$de_top$feature)
#
# Subset the table
x <- a2$de %>% dplyr::filter(feature %in% these_genes)
x <- as.data.frame(dcast.data.table(
  data = as.data.table(x), formula = symbol ~ group, value.var = "auc"
))
rownames(x) <- x$symbol
x$symbol <- NULL
# Order the heatmap
set.seed(1)
xo <- seriation::seriate(as.matrix(x), method = "BEA_TSP")
xd <- a2$de %>% dplyr::filter(feature %in% these_genes)
xd$symbol <- factor(xd$symbol, rownames(x)[xo[[1]]])
xd$group <- factor(xd$group, colnames(x)[xo[[2]]])
# Plot
p2 <- ggplot(xd) +
geom_tile(
  aes(y = symbol, x = group, fill = auc)
) +
scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
# scale_fill_viridis_c(
#   name = "AUC", guide = guide_colorbar(barheight = 40), breaks = pretty_breaks(10),
#   limits = c(0, 1)
# ) +
scale_fill_scico(
  palette = "vik", limits = c(0, 1), direction = 1,
  name = "AUC", guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(9)
) +
theme(
  axis.text.y = element_text(size = 12, face = "italic")
#   axis.text.y = element_blank(), axis.ticks.y = element_blank()
)
p1 <- ggplot(xd) +
geom_tile(
  aes(y = 1, x = group, fill = group)
) +
scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
scale_fill_manual(
  values = cluster_colors, guide = "none"
)
p <- (
  p1 + theme(plot.margin = margin(b = 0))
) / (
  p2 + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 0)
  )
) + plot_layout(heights = c(1, length(these_genes)))
#p
fig_height <- length(these_genes) * 0.2
fig_width <- length(unique(a2$de_top$group)) * 0.4 + 3
my_ggsave(
  slug = "heatmap-top-markers",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
)


# Pseudobulk
########################################################################

if ("case" %in% colnames(a2$obs) && n_donors > 1) {

  y <- with(a2$obs, model.matrix(~ 0 + factor(leiden):factor(donor)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/")
  pb <- as(a2$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(colSums(pb)))
  #
  library(limma)
  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
  colnames(pb_meta) <- c("leiden", "donor")
  pb_meta <- as_tibble(pb_meta)
  pb_meta %<>%
    mutate(
      leiden = str_replace(leiden, "factor\\(leiden\\)", ""),
      donor = str_replace(donor, "factor\\(donor\\)", "")
    )
  pb_meta <- left_join(
    pb_meta,
    a2$obs %>%
      select(
        donor, case
      ) %>%
      group_by(donor, case) %>%
      summarize_if(is.numeric, mean),
    by = "donor"
  )
  # pb_meta <- left_join(
  #   pb_meta,
  #   sample_info %>%
  #     select(
  #       # FIXME: this needs to be an rgument to a function
  #       donor, case#, chemistry, qubit_library_quantification_ng_ul, facs_sorting
  #     ) %>%
  #     # FIXME: this needs to be an rgument to a function
  #     # group_by(donor, case, chemistry) %>%
  #     group_by(donor, case) %>%
  #     summarize_if(is.numeric, mean),
  #   by = "donor"
  # )
  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
  stopifnot(nrow(pb_meta) == ncol(pb))

  # this_ens <- names(which(ensembl_to_symbol == "PTPRC"))
  # pb_meta$gene <- pb[this_ens,]
  # pb_meta$leiden <- fct_reorder(
  #   pb_meta$leiden, pb_meta$gene,
  #   .fun = mean
  # )
  # ggplot(pb_meta) +
  #   aes(y = leiden, x = gene) +
  #   geom_quasirandom(, groupOnX = FALSE)
  # What percent of cells have this gene?
  # a2$obs$gene <- a2$log2cpm[this_ens,]
  # a2$obs %>%
  #   group_by(leiden) %>%
  #   summarize(pct = sum(gene > 0) / length(gene)) %>%
  # ggplot() +
  # aes(y = as.character(leiden), x = pct) +
  # geom_colh()

  x <- table(a2$obs$leiden)
  min_percent <- 100 * min(x) * 0.25 / sum(x)
  keep_ens <- a2$counts_stats$gene[a2$counts_stats$percent >= min_percent]

  # All versus All (AVA)
  # Test all pairs of clusters
  print_status("Finding marker genes with pseudobulk")
  pb_meta$x <- factor(pb_meta$leiden)
  des1 <- with(pb_meta, model.matrix(
    # FIXME: this needs to be an argument to a function
    ~ 0 + x #+ chemistry + qubit_library_quantification_ng_ul
  ))
  ob <- as.matrix(pb[keep_ens,])
  # ob <- ob[rowMeans(ob) > 0.5,]
  fit1 <- lmFit(object = ob, design = des1)
  fit1 <- eBayes(fit1)
  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  cluster_pairs <- t(combn(levels(pb_meta$x), 2))
  cont <- makeContrasts(contrasts = lapply(seq(nrow(cluster_pairs)), function(i) {
    glue("x{cluster_pairs[i,1]} - x{cluster_pairs[i,2]}")
  }), levels = des1)
  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
  fit2 <- contrasts.fit(fit1, cont)
  fit2 <- eBayes(fit2)
  de_ava <- rbindlist(lapply(colnames(cont), function(this_coef) {
    x <- topTable(fit2, coef = this_coef, number = nrow(fit1$coefficients))
    this_coef <- str_replace_all(this_coef, "x", "")
    this_coef <- str_replace(this_coef, "vs", " vs ")
    x$coef <- this_coef
    x$ensembl_id <- rownames(x)
    x
  }))
  # FIXME: this needs to be an argument to a function
  de_ava_file <- glue("{out_dir}/pseudobulk_de_ava.tsv.gz")
  print_status(glue("Writing {de_ava_file}"))
  data.table::fwrite(de_ava, de_ava_file, sep = "\t")

  # One versus All (OVA)
  print_status("Finding marker genes with pseudobulk")
  de_ova <- rbindlist(
    pblapply(sort(unique(pb_meta$leiden)), function(this_cluster) {
      pb_meta$x <- pb_meta$leiden == this_cluster
      des1 <- with(pb_meta, model.matrix(
        # FIXME: this needs to be an argument to a function
        ~ x #+ chemistry + qubit_library_quantification_ng_ul
      ))
      ob <- as.matrix(pb[keep_ens,])
      # ob <- ob[rowMeans(ob) > 0.5,]
      fit1 <- lmFit(object = ob, design = des1)
      fit1 <- eBayes(fit1)
      fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
      res <- topTable(fit1, coef = 2, number = 1e6)
      res$coef <- sprintf("%s vs all", this_cluster)
      res$ensembl_id <- rownames(res)
      return(res)
    })
  )
  de_ova <- left_join(
    x = de_ova %>% mutate(group = str_remove(coef, ' vs all')),
    y = a2$de %>% select(group, feature, auc, pct_in, pct_out),
    by = c("ensembl_id" = "feature", "group" = "group")
  )
  de_ova_file <- glue("{out_dir}/pseudobulk_de_ova.tsv.gz")
  print_status(glue("Writing {de_ova_file}"))
  data.table::fwrite(de_ova, de_ova_file, sep = "\t")

  # head(de_ova, 3)
  # head(a2$de, 3)

  min_percent <- 1
  max_pval <- 0.05 / nrow(de_ova)
  min_logfc <- log2(1.5)
  lab_pval <- scientific_10(max_pval)
  # x <- qvalue::qvalue(de_ova$P.Value)
  # de_ova$fdr <- x$qvalues
  x <- de_ova %>%
    group_by(coef) %>%
    summarize(n = sum(pct_in > min_percent & P.Value < max_pval & logFC > min_logfc)) %>%
    mutate(cluster = str_split_fixed(coef, " ", 2)[,1])
  p <- ggplot(x) +
    aes(x = n, y = reorder(coef, n), fill = cluster) +
    geom_colh() +
    geom_text(aes(label = n), hjust = -0.1, size = 5) +
    scale_x_continuous(expand = expansion(c(0, 0.1))) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    labs(
      x = "Differentially expressed genes",
      y = "Contrast",
      title = "Differentially expressed genes for each cluster",
      subtitle = parse(text = glue(
        'P < {scientific_10(max_pval)} ~ "and Fold-change >" ~ {signif(2^min_logfc, 2)} ~ "and Percent >" ~ {min_percent}'
      ))
      # subtitle = parse(text = bquote("Fold-change > " ~ .(2^min_logfc) ~ " and P < "~lab_pval))
    )
  my_ggsave(
    slug = "cluster-de-genes-pseudobulk-ova",
    out_dir = out_dir,
    plot = p,
    type = "pdf",
    scale = 1,
    width = 6,
    height = 0.5 + 0.25 * length(unique(de_ova$coef)),
    units = "in", dpi = 300
  )

  # TODO: Turn this into volcano for each ova cluster
  for (this_cluster in unique(pb_meta$leiden)) {
    top1 <- de_ova %>%
      mutate(cluster = str_split_fixed(coef, " ", 2)[,1]) %>%
      filter(cluster == this_cluster, pct_in > min_percent) %>%
      mutate(Gene = ID)
    p <- plot_limma_volcano(top1) +
      labs(title = glue("Cluster {this_cluster} vs All"))
    my_ggsave(
      glue("volcano-{this_cluster}"),
      #out_dir = glue("figures/{analysis_name}/cluster_volcano"),
      out_dir = glue("{out_dir}/cluster_volcano"),
      plot = p,
      type = "pdf",
      scale = 1, width = 8, height = 6, units = "in", dpi = 300
    )
  }

  print_status("skipping individual gene plots in cluster_volcano/genes")
  if (FALSE) {
    for (this_cluster in unique(pb_meta$leiden)) {
      top1 <- de_ova %>%
        mutate(cluster = str_split_fixed(coef, " ", 2)[,1]) %>%
        filter(cluster == this_cluster, pct_in > min_percent) %>%
        mutate(Gene = ID)
      my_genes <- rev(top1$ensembl_id[1:10])
      for (my_gene in my_genes) {
        # my_gene <- names(ensembl_to_symbol[ensembl_to_symbol == "LINC01781"])
        a2$obs$value <- as.numeric(a2$log2cpm[my_gene,])
        x <- a2$obs %>%
          dplyr::mutate(leiden = as.character(leiden)) %>%
          dplyr::group_by(donor, leiden) %>%
          dplyr::summarize(
            mean = mean(value),
            median = median(value),
            pct  = 100 * sum(value > 0) / n(),
            xmin = quantile(value, 0.75),
            x    = quantile(value, 0.5),
            xmax = quantile(value, 0.25)
          )
        p <- ggplot(x) +
          aes(x = mean, y = reorder(leiden, mean), group = donor) +
          geom_point(alpha = 0) +
          annotate(
            geom = "rect",
            xmin = -Inf, xmax = Inf,
            ymin = seq(1, length(unique(x$leiden))) - 0.5,
            ymax = seq(1, length(unique(x$leiden))) + 0.5,
            fill = rep(c("white", "grey90"), length.out = length(unique(x$leiden)))
          ) +
            ggstance::geom_crossbarh(
              data = x %>%
                dplyr::group_by(leiden) %>%
                dplyr::summarize(
                  xmin = quantile(mean, 0.75),
                  x    = quantile(mean, 0.5),
                  xmax = quantile(mean, 0.25)
                ),
              mapping = aes(
                y = leiden, xmin = xmin, x = x, xmax = xmax,
                group = leiden, fill = leiden
              ),
              color = "grey20",
              width = 0.8,
              alpha = 0.6,
              position = position_dodge(width = 0.9)
            ) +
          scale_fill_manual(values = cluster_colors, guide = "none") +
          geom_quasirandom(
            groupOnX = FALSE,
            width = 0.2
          ) +
          labs(
            x = bquote("Mean Log"[2]~"CPM"),
            y = NULL,
            title = ensembl_to_symbol[my_gene]
          ) +
          theme(plot.title = element_text(face = "italic"))
        my_ggsave(
          glue("{ensembl_to_symbol[my_gene]}-{my_gene}"),
          out_dir = glue("{out_dir}/cluster_volcano/genes"),
          plot = p,
          type = "pdf",
          scale = 1,
          width = 4,
          height = 0.5 + 0.3 * length(unique(x$leiden)),
          units = "in", dpi = 300
        )
      }
    }
  }

  de_case <- rbindlist(lapply(
    unique(pb_meta$leiden),
    function(this_cluster) {
      # this_cluster <- 13
      ix <- which(pb_meta$leiden == this_cluster)
      des1 <- with(
        pb_meta[ix,],
        model.matrix(~ case)
      )
      fit1 <- lmFit(object = as.matrix(pb[rowMeans(pb) > 0.5, ix]), design = des1)
      fit1 <- eBayes(fit1)
      fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
      top1 <- topTable(
        fit1, coef = "caseCase", number = 1e6, confint = FALSE)
      top1$coef <- glue("case vs control c{this_cluster}")
      top1$ensembl_id <- rownames(top1)
      return(top1)
    }
  ))

  for (this_cluster in unique(pb_meta$leiden)) {
    # this_cluster <- 13
    ix <- which(pb_meta$leiden == this_cluster)
    des1 <- with(
      pb_meta[ix,],
      # FIXME: This needs to be an argument to a function
      # model.matrix(~ case + chemistry + qubit_library_quantification_ng_ul)
      model.matrix(~ case)
    )
    fit1 <- lmFit(object = as.matrix(pb[rowMeans(pb) > 0.5, ix]), design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
    top1 <- topTable(fit1, coef = "caseCase", number = 1e6)
    # top1 <- topTable(fit1, coef = "chemistry5pV2", number = 1e6)
    # if ("Gene" %in% colnames(top1)) {
    #   top1 %<>% rename(ID = "Gene")
    # }
    # top1 %<>% dplyr::rename(Gene = ID)
    colnames(top1)[1] <- "Gene"
    p <- plot_limma_volcano(top1) +
      labs(title = glue("Case vs Control (Cluster {this_cluster})"))
    my_ggsave(
      glue("volcano-{this_cluster}"),
      #out_dir = glue("figures/{analysis_name}/cluster_volcano"),
      out_dir = glue("{out_dir}/case_control/cluster_volcano"),
      plot = p,
      type = "pdf",
      scale = 1, width = 8, height = 6, units = "in", dpi = 300
    )
    my_genes <- rev(rownames(top1)[1:10])
    x <- cbind.data.frame(pb_meta[ix,], t(as.matrix(pb[my_genes,ix])))
    x <- reshape2::melt(x, id.vars = c("donor", "case"), measure.vars = my_genes)
    p <- ggplot(x) +
      annotate(
        geom = "rect",
        xmin = -Inf, xmax = Inf,
        ymin = seq(1, length(my_genes)) - 0.5,
        ymax = seq(1, length(my_genes)) + 0.5,
        fill = rep(c("white", "grey90"), length.out = length(my_genes))
      ) +
      ggstance::geom_crossbarh(
        data = x %>% dplyr::group_by(variable, case) %>%
          summarize(
            xmin = quantile(value, 0.75),
            x    = quantile(value, 0.5),
            xmax = quantile(value, 0.25)
          ),
        mapping = aes(
          y = variable, xmin = xmin, x = x, xmax = xmax,
          group = case, fill = case
        ),
        color = "grey20",
        width = 0.8,
        alpha = 0.6,
        position = position_dodge(width = 0.9)
      ) +
      geom_quasirandom(
        data = x,
        mapping = aes(x = value, y = variable, fill = case, group = factor(case)),
        dodge.width = 0.8,
        groupOnX = FALSE,
        width = 0.2,
        shape = 21,
        size = 3,
        stroke = 0.2
      ) +
      scale_fill_manual(values = my_colors$case) +
      guides(fill = guide_legend(title = NULL, reverse = TRUE)) +
      scale_y_discrete(labels = ensembl_to_symbol, expand = c(0, 0)) +
      labs(
        x = bquote("Log"[2]~"CPM"),
        y = NULL,
        title = glue("Case vs Control (Cluster {this_cluster})")
      ) +
      theme(
        axis.ticks.y = element_blank(),
        # legend.position = c(1, 0),
        # legend.justification = c(1, 0),
        # legend.background = element_rect(color = "grey10", fill = NA, size = 0.2),
        axis.text.y = element_text(face = "italic")
      )
    my_ggsave(
      glue("top-genes-{this_cluster}"),
      #out_dir = glue("figures/{analysis_name}/cluster_top"),
      out_dir = glue("{out_dir}/case_control/cluster_top"),
      plot = p,
      type = "pdf",
      scale = 1, width = 8, height = 6, units = "in", dpi = 300
    )
    print_status("skipping individual gene plots in cluster_top/genes")
    if (FALSE) {
      for (my_gene in my_genes) {
        # my_gene <- names(ensembl_to_symbol[ensembl_to_symbol == "LINC01781"])
        a2$obs$value <- as.numeric(a2$log2cpm[my_gene,])
        x <- a2$obs %>%
          dplyr::mutate(leiden = as.character(leiden)) %>%
          dplyr::group_by(donor, leiden) %>%
          dplyr::summarize(
            mean = mean(value),
            median = median(value),
            pct  = 100 * sum(value > 0) / n(),
            xmin = quantile(value, 0.75),
            x    = quantile(value, 0.5),
            xmax = quantile(value, 0.25)
          )
        p <- ggplot(x) +
          aes(x = mean, y = reorder(leiden, mean), group = donor) +
          geom_point(alpha = 0) +
          annotate(
            geom = "rect",
            xmin = -Inf, xmax = Inf,
            ymin = seq(1, length(unique(x$leiden))) - 0.5,
            ymax = seq(1, length(unique(x$leiden))) + 0.5,
            fill = rep(c("white", "grey90"), length.out = length(unique(x$leiden)))
          ) +
            ggstance::geom_crossbarh(
              data = x %>%
                dplyr::group_by(leiden) %>%
                dplyr::summarize(
                  xmin = quantile(mean, 0.75),
                  x    = quantile(mean, 0.5),
                  xmax = quantile(mean, 0.25)
                ),
              mapping = aes(
                y = leiden, xmin = xmin, x = x, xmax = xmax,
                group = leiden, fill = leiden
              ),
              color = "grey20",
              width = 0.8,
              alpha = 0.6,
              position = position_dodge(width = 0.9)
            ) +
          scale_fill_manual(values = cluster_colors, guide = "none") +
          geom_quasirandom(
            groupOnX = FALSE,
            width = 0.2
          ) +
          labs(
            x = bquote("Mean Log"[2]~"CPM"),
            y = NULL,
            title = ensembl_to_symbol[my_gene]
          ) +
          theme(plot.title = element_text(face = "italic"))
        my_ggsave(
          glue("{ensembl_to_symbol[my_gene]}-{my_gene}"),
          out_dir = glue("{out_dir}/case_control/cluster_top/genes"),
          plot = p,
          type = "pdf",
          scale = 1,
          width = 4,
          height = 0.5 + 0.3 * length(unique(x$leiden)),
          units = "in", dpi = 300
        )
      }
    }
  }

}

# Pseudobulk by donor

if (do_pb && n_donors > 1) {

  y <- with(a2$obs, model.matrix(~ 0 + donor))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/")
  pb <- as(a2$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "^donor", "")
  #
  pb_meta <- a2$obs %>% dplyr::select(
    case, class, class2, class_short,
    donor, drug
  ) %>% unique
  pb_meta <- inner_join(
    data.frame(donor = colnames(pb)),
    pb_meta,
    by = "donor"
  )
  pb_meta <- as_tibble(pb_meta)
  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
  all(pb_meta$donor == colnames(pb))
  # if (!"drug" %in% colnames(pb_meta)) {
  #   pb_meta <- left_join(
  #     pb_meta, sample_therapy, by = "donor"
  #   )
  # }
  #

  x <- table(a2$obs$leiden)
  min_percent <- 100 * min(x) * 0.25 / sum(x)
  keep_ens <- a2$counts_stats$gene[a2$counts_stats$percent >= min_percent]
  #
  library(limma)
  des1 <- with(pb_meta, model.matrix(
    ~ case #+ chemistry + qubit_library_quantification_ng_ul
  ))
  ob <- as.matrix(pb[keep_ens,])
  # ob <- ob[rowMeans(ob) > 0.5,]
  fit1 <- lmFit(object = ob, design = des1)
  fit1 <- eBayes(fit1)
  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  res <- topTable(fit1, coef = 2, number = 1e6, confint = TRUE)
  res %<>% dplyr::rename(Gene = ID)
  res$ensembl_id <- rownames(res)

  dir.create(glue("{out_dir}/pb_donor"), showWarnings = FALSE)
  fwrite(
    res, glue("{out_dir}/pb_donor/de_case-vs-control.tsv.gz"),
    sep = "\t"
  )
  #
  p <- plot_limma_volcano(res, fdr = 0.05, n_text = 20) +
    labs(title = glue("Case vs Control"))
  p
  # Molly grant
  # text_ids <- res$ensembl_id[
  #   # which(res$Gene %in% c("IFNG", "IL26", "TNF", "IL17A", "GZMB"))
  #   # which(res$Gene %in% strsplit("CD274, ISG20, STAT1, GBP5, CXCL9, CXCL10, CXCL11", ", ")[[1]])
  #   which(res$Gene %in% strsplit("CXCL9, CXCL10, CXCL11, HLA-DRA, CD274, GBP4, IDO1, ISG15, ISG20", ", ")[[1]])
  # ]
  # p <- plot_limma_volcano(res, fdr = 0.05, n_text = 0, text_ids = text_ids) +
  #   labs(title = glue("Case vs Control"))
  # p
  my_ggsave(
    glue("volcano-case"),
    out_dir = glue("{out_dir}/pb_donor"),
    plot = p,
    type = "pdf",
    scale = 1, width = 8, height = 6, units = "in", dpi = 300
  )

  # my_ens <- res$ensembl_id[res$adj.P.Val < 0.05]
  # my_ens <- head(my_ens, 40)
  # length(my_ens)
  my_ens <- head(res$ensembl_id, 40)
  mat <- as.matrix(pb[my_ens,])
  all(colnames(mat) == pb_meta$donor)
  a_col <- as.data.frame(pb_meta[,c("class", "case", "drug")])
  # a_col$facs_sorting <- as.character(a_col$facs_sorting)
  a_colors <- list(
    facs_sorting = c(
      "TRUE"  = "grey60",
      "FALSE" = "white"
    ),
    case = c(
      "Case"    = Okabe_Ito[1],
      "Control" = Okabe_Ito[2]
    ),
    class = c(
      "irColitis"             = Okabe_Ito[1],
      "On ICI therapy"        = Okabe_Ito[4],
      "Chronic diarrhea"      = Okabe_Ito[7],
      "Screening colonoscopy" = Okabe_Ito[3]
    ),
    drug = c(
      "None"        = Okabe_Ito[8],
      "CTLA-4"      = Okabe_Ito[4],
      "PD-1"        = Okabe_Ito[7],
      "PD-1/CTLA-4" = Okabe_Ito[6]
    )
  )
  rownames(a_col) <- pb_meta$donor
  labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
    sprintf("%s%s", x[1], x[2])
  })
  # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
  # mat_o <- seriation::seriate(mat, method = "PCA")
  # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
  heatmap_file <- glue("{out_dir}/pb_donor/heatmap-case-top40.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.15 + 5,
    height = nrow(mat) * 0.15 + 2,
    mat = mat,
    hclust_method = "complete",
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    # show_colnames = FALSE,
    # color = rev(scico::scico(palette = "davos", n = 20)),
    scale = "row",
    color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )

  # plots <- lapply(head(my_ens, 10), function(this_ens) {
  #   pb_meta$gene <- ob[this_ens,]
  #   ggplot() +
  #     geom_crossbar(
  #       data = pb_meta %>%
  #         group_by(case) %>%
  #         summarize(
  #           y = quantile(gene, 0.5),
  #           ymin = quantile(gene, 0.25),
  #           ymax = quantile(gene, 0.75),
  #           .groups = "drop"
  #         ),
  #       mapping = aes(x = case, y = y, ymin = ymin, ymax = ymax, color = case)
  #     ) +
  #     geom_quasirandom(
  #       data = pb_meta,
  #       mapping = aes(x = case, y = gene, group = case, color = case),
  #       width = 0.2,
  #       size = 3
  #     ) +
  #     guides(color = "none") +
  #     labs(
  #       x = NULL, y = bquote("Log"[2]~"CPM"),
  #       title = rowname_key[this_ens]
  #     ) +
  #     theme(title = element_text(face = "italic"))
  # })
  # wrap_plots(plots, ncol = 5)

  plot_ens <- head(my_ens, 20)
  x <- rbindlist(lapply(plot_ens, function(this_ens) {
    cbind(pb_meta,
      ens = this_ens,
      gene = ob[this_ens,],
      symbol = rowname_key[this_ens]
    )
  }))
  x$symbol <- factor(x$symbol, rowname_key[plot_ens])
  p <- ggplot() +
    geom_crossbar(
      data = x %>%
        group_by(symbol, case) %>%
        summarize(
          y = quantile(gene, 0.5),
          ymin = quantile(gene, 0.25),
          ymax = quantile(gene, 0.75),
          .groups = "drop"
        ),
      mapping = aes(x = case, y = y, ymin = ymin, ymax = ymax, color = case)
    ) +
    geom_quasirandom(
      data = x,
      mapping = aes(x = case, y = gene, group = case, color = case),
      width = 0.4,
      size = 3
    ) +
    facet_wrap(~ symbol, scales = "free_y", ncol = 5) +
    guides(color = "none") +
    labs(
      x = NULL, y = bquote("Log"[2]~"CPM")
    ) +
    theme(strip.text = element_text(face = "italic"))
  my_ggsave(
    glue("point-case-top20"),
    out_dir = glue("{out_dir}/pb_donor"),
    plot = p,
    type = "pdf",
    scale = 1.4, width = 9, height = 6, units = "in", dpi = 300
  )

  # pb_fit <- do_loess(pb, exclude_genes, loess_span = 0.3, min_percent = 0.05)
  # ix <- which(pb_fit$rank < 5000)
  # pb_pca <- RSpectra::svds(
  #   A    = t(pb[ix,]),
  #   k    = 10,
  #   opts = list(
  #     center = TRUE, scale = TRUE, maxitr = 2000, tol = 1e-10
  #   )
  # )
  # pb_pca$genes <- rownames(pb)[ix]
  # pb_meta$PC1 <- pb_pca$u[,1]
  # pb_meta$PC2 <- pb_pca$u[,2]
  # pb_meta$PC3 <- pb_pca$u[,3]
  # pb_meta$PC4 <- pb_pca$u[,4]
  # pb_meta$PC5 <- pb_pca$u[,5]
  # pb_meta$PC6 <- pb_pca$u[,6]
  # ggplot(pb_meta) +
  #   aes(PC1, PC2, fill = class) +
  #   geom_point(shape = 21, size = 5)
  # ggplot(pb_meta) +
  #   aes(PC3, PC4, fill = class) +
  #   geom_point(shape = 21, size = 5)
  # ggplot(pb_meta) +
  #   aes(PC5, PC6, fill = class) +
  #   geom_point(shape = 21, size = 5)

  library(limma)
  pb_meta$classf <- factor(
    pb_meta$class,
    c(
      "Screening colonoscopy",
      "Chronic diarrhea",
      "On ICI therapy",
      "irColitis"
    )
    # pb_meta$class_short,
    # c("SC", "CD", "ICI", "IRC")
  )
  des2 <- with(pb_meta, model.matrix(
    ~ classf # case # + chemistry#  + facs_sorting
  ))
  ob <- as.matrix(pb[keep_ens,])
  # ob <- ob[rowMeans(ob) > 0.5,]
  fit2 <- lmFit(object = ob, design = des2)
  fit2 <- eBayes(fit2)
  fit2$genes <- ensembl_to_symbol[rownames(fit2$coefficients)]
  # x <- lmCompare(
  #   object = as.matrix(pb[rowMeans(pb) > 0.5,]),
  #   fit1, fit2
  # ) %>% arrange(pval)
  # gg_qqplot(x$pval)
  res <- topTable(fit2, number = 1e6, confint = TRUE)
  colnames(res)[1] <- "Gene"
  res$ensembl_id <- rownames(res)
  head(res)
  dir.create(glue("{out_dir}/pb_donor"), showWarnings = FALSE)
  fwrite(
    res, glue("{out_dir}/pb_donor/de_class-anova.tsv.gz"),
    sep = "\t"
  )
  # my_ens <- head(res[order(res$P.Value),]$ensembl_id, 40)
  # my_ens <- res$ensembl_id[res$adj.P.Val < 0.01]
  # my_ens <- head(res$ensembl_id[res$adj.P.Val < 0.2], 40)
  my_ens <- head(res$ensembl_id, 40)
  if (length(my_ens) > 1) {
  mat <- as.matrix(pb[my_ens,])
  all(colnames(mat) == pb_meta$donor)
  a_col <- as.data.frame(pb_meta[,c("class", "case", "drug")])
  # a_col$facs_sorting <- as.character(a_col$facs_sorting)
  a_colors <- list(
    facs_sorting = c(
      "TRUE"  = "grey60",
      "FALSE" = "white"
    ),
    case = c(
      "Case"    = Okabe_Ito[1],
      "Control" = Okabe_Ito[2]
    ),
    class = c(
      "irColitis"             = Okabe_Ito[1],
      "On ICI therapy"        = Okabe_Ito[4],
      "Chronic diarrhea"      = Okabe_Ito[7],
      "Screening colonoscopy" = Okabe_Ito[3]
    ),
    drug = c(
      "None"        = Okabe_Ito[8],
      "CTLA-4"      = Okabe_Ito[4],
      "PD-1"        = Okabe_Ito[7],
      "PD-1/CTLA-4" = Okabe_Ito[6]
    )
  )
  rownames(a_col) <- pb_meta$donor
  labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
    sprintf("%s%s", x[1], x[2])
  })
  # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
  # mat_o <- seriation::seriate(mat, method = "PCA")
  # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
  heatmap_file <- glue("{out_dir}/pb_donor/heatmap-class-anova-top40.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.15 + 4,
    height = nrow(mat) * 0.15 + 1,
    mat = mat,
    hclust_method = "complete",
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    # show_colnames = FALSE,
    # color = rev(scico::scico(palette = "davos", n = 20)),
    scale = "row",
    color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )
  }

  # Drug
  ########################################################################
  pb_meta$drug <- factor(pb_meta$drug,
    c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4")
  )
  # pb_meta$drug2 <- as.character(pb_meta$drug)
  # pb_meta$drug2[pb_meta$drug2 == "PD-1/CTLA-4"] <- "PD-1"
  # pb_meta$drug2 <- factor(pb_meta$drug2,
  #   c("None", "CTLA-4", "PD-1")
  # )
  ix <- pb_meta$case == "Case" & pb_meta$drug != "CTLA-4"
  pb_meta$drug <- as.character(pb_meta$drug)
  des1 <- with(pb_meta[ix,], model.matrix(
    ~ drug
  ))
  ob <- as.matrix(pb[keep_ens,ix])
  # ob <- ob[rowMeans(ob) > 0.5,]
  fit1 <- lmFit(object = ob, design = des1)
  fit1 <- eBayes(fit1)
  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  res <- topTable(fit1, coef = "drugPD-1/CTLA-4", number = 1e6)
  # res <- topTable(fit1, number = 1e6)
  # res %<>% dplyr::rename(Gene = ID)
  colnames(res)[1] <- "Gene"
  res$ensembl_id <- rownames(res)
  p <- plot_limma_volcano(res, fdr = 0.1, n_text = 20) +
    labs(title = glue("CTLA-4 vs None"))
  my_ggsave(
    glue("volcano-CTLA-4"),
    out_dir = glue("{out_dir}/pb_donor"),
    plot = p,
    type = "pdf",
    scale = 1, width = 8, height = 6, units = "in", dpi = 300
  )

  # my_ens <- res$ensembl_id[res$adj.P.Val < 0.1]
  my_ens <- head(res$ensembl_id, 40)
  if (length(my_ens) > 1) {
  plot_ens <- head(my_ens, 20)
  x <- rbindlist(lapply(plot_ens, function(this_ens) {
    cbind(pb_meta[ix,],
      ens = this_ens,
      gene = ob[this_ens,],
      symbol = head(rowname_key[this_ens], 1)
    )
  }))
  x$symbol <- factor(x$symbol, rowname_key[plot_ens])
  p <- ggplot() +
    geom_crossbar(
      data = x %>%
        group_by(symbol, drug) %>%
        summarize(
          y = quantile(gene, 0.5),
          ymin = quantile(gene, 0.25),
          ymax = quantile(gene, 0.75),
          .groups = "drop"
        ),
      mapping = aes(x = drug, y = y, ymin = ymin, ymax = ymax, color = drug)
    ) +
    geom_quasirandom(
      data = x,
      mapping = aes(x = drug, y = gene, group = drug, color = drug),
      width = 0.4,
      size = 3
    ) +
    facet_wrap(~ symbol, scales = "free_y", ncol = 5) +
    guides(color = "none") +
    labs(
      x = NULL, y = bquote("Log"[2]~"CPM")
    ) +
    theme(
      strip.text = element_text(face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  my_ggsave(
    glue("point-CTLA-4-top20"),
    out_dir = glue("{out_dir}/pb_donor"),
    plot = p,
    type = "pdf",
    scale = 1.4, width = 9, height = 6, units = "in", dpi = 300
  )
  }
  
  # res <- topTable(fit1, coef = "drugPD-1", number = 1e6)
  # # res <- topTable(fit1, number = 1e6)
  # # res %<>% dplyr::rename(Gene = ID)
  # colnames(res)[1] <- "Gene"
  # res$ensembl_id <- rownames(res)
  # p <- plot_limma_volcano(res, fdr = 0.2, n_text = 20) +
  #   labs(title = glue("PD-1 vs None"))
  # my_ggsave(
  #   glue("volcano-PD-1"),
  #   out_dir = glue("{out_dir}/pb_donor"),
  #   plot = p,
  #   type = "pdf",
  #   scale = 1, width = 8, height = 6, units = "in", dpi = 300
  # )
  # res <- topTable(fit1, coef = "drugPD-1/CTLA-4", number = 1e6)
  # # res <- topTable(fit1, number = 1e6)
  # # res %<>% dplyr::rename(Gene = ID)
  # colnames(res)[1] <- "Gene"
  # res$ensembl_id <- rownames(res)
  # p <- plot_limma_volcano(res, fdr = 0.2, n_text = 20) +
  #   labs(title = glue("PD-1 and CTLA-4 vs None"))
  # my_ggsave(
  #   glue("volcano-PD-1-and-CTLA-4"),
  #   out_dir = glue("{out_dir}/pb_donor"),
  #   plot = p,
  #   type = "pdf",
  #   scale = 1, width = 10, height = 6, units = "in", dpi = 300
  # )

  res <- topTable(fit1, number = 1e6)
  colnames(res)[1] <- "Gene"
  res$ensembl_id <- rownames(res)
  # my_ens <- res$ensembl_id[res$adj.P.Val < 0.2]
  my_ens <- head(res$ensembl_id, 40)
  length(my_ens)
  mat <- as.matrix(pb[my_ens,])
  all(colnames(mat) == pb_meta$donor)
  a_col <- as.data.frame(pb_meta[,c("class", "case", "drug")])
  # a_col$facs_sorting <- as.character(a_col$facs_sorting)
  a_colors <- list(
    facs_sorting = c(
      "TRUE"  = "grey60",
      "FALSE" = "white"
    ),
    case = c(
      "Case"    = Okabe_Ito[1],
      "Control" = Okabe_Ito[2]
    ),
    class = c(
      "irColitis"             = Okabe_Ito[1],
      "On ICI therapy"        = Okabe_Ito[4],
      "Chronic diarrhea"      = Okabe_Ito[7],
      "Screening colonoscopy" = Okabe_Ito[3]
    ),
    drug = c(
      "None"        = Okabe_Ito[8],
      "CTLA-4"      = Okabe_Ito[4],
      "PD-1"        = Okabe_Ito[7],
      "PD-1/CTLA-4" = Okabe_Ito[6]
    )
  )
  rownames(a_col) <- pb_meta$donor
  labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
    sprintf("%s%s", x[1], x[2])
  })
  # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
  # mat_o <- seriation::seriate(mat, method = "PCA")
  # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
  heatmap_file <- glue("{out_dir}/pb_donor/heatmap-drug-anova-top40.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.15 + 4,
    height = nrow(mat) * 0.15 + 2,
    mat = mat,
    hclust_method = "complete",
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    # show_colnames = FALSE,
    # color = rev(scico::scico(palette = "davos", n = 20)),
    scale = "row",
    color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )

  # Only Controls
  ix <- pb_meta$case == "Control"
  des1 <- with(pb_meta[ix,], model.matrix(
    ~ drug != "None"
  ))
  colnames(des1)[2] <- "drug"
  ob <- as.matrix(pb[keep_ens,ix])
  # ob <- ob[rowMeans(ob) > 0.5,]
  fit1 <- lmFit(object = ob, design = des1)
  fit1 <- eBayes(fit1)
  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  res <- topTable(fit1, coef = "drug", number = 1e6)
  # res <- topTable(fit1, number = 1e6)
  # res %<>% dplyr::rename(Gene = ID)
  colnames(res)[1] <- "Gene"
  res$ensembl_id <- rownames(res)
  p <- plot_limma_volcano(res, fdr = 0.2) +
    labs(title = glue("Drug vs None"))
  p
  my_ggsave(
    glue("volcano_control-only_drug"),
    out_dir = glue("{out_dir}/pb_donor"),
    plot = p,
    type = "pdf",
    scale = 1, width = 8, height = 6, units = "in", dpi = 300
  )

}

# Heatmap with clusters, averaged over donors

if (do_pb && n_donors > 1) {

  y <- with(a2$obs, model.matrix(~ 0 + class_short:factor(leiden)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/")
  pb <- as(a2$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "^class_short", "")
  colnames(pb) <- str_replace(colnames(pb), "factor\\(leiden\\)", "")
  colnames(pb)
  #
  pb_meta <- a2$obs %>% dplyr::select(
    class_short, leiden, case, class
  ) %>% unique
  stopifnot(nrow(pb_meta) == ncol(pb))
  pb_meta$id <- glue("{pb_meta$class_short}:{pb_meta$leiden}")
  pb_meta <- pb_meta[match(colnames(pb), pb_meta$id),]
  stopifnot(all(pb_meta$id == colnames(pb)))
  pb_meta <- as_tibble(pb_meta)
  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))

  x <- table(a2$obs$leiden)
  min_percent <- 100 * min(x) * 0.25 / sum(x)
  keep_ens <- a2$counts_stats$gene[a2$counts_stats$percent >= min_percent]

    my_ens <- res$ensembl_id[1:40]
    mat <- as.matrix(pb[my_ens,ix])
    all(colnames(mat) == glue("{pb_meta$donor[ix]}:{this_cluster}"))
    colnames(mat) <- str_split_fixed(colnames(mat), ":", 2)[,1]
    a_col <- as.data.frame(pb_meta[ix,c("class", "case", "drug")])
    # a_col$facs_sorting <- as.character(a_col$facs_sorting)
    a_colors <- list(
      facs_sorting = c(
        "TRUE"  = "grey60",
        "FALSE" = "white"
      ),
      case = c(
        "Case"    = Okabe_Ito[1],
        "Control" = Okabe_Ito[2]
      ),
      class = c(
        "irColitis"             = Okabe_Ito[1],
        "On ICI therapy"        = Okabe_Ito[4],
        "Chronic diarrhea"      = Okabe_Ito[7],
        "Screening colonoscopy" = Okabe_Ito[3]
      ),
      drug = c(
        "None"        = Okabe_Ito[8],
        "CTLA-4"      = Okabe_Ito[4],
        "PD-1"        = Okabe_Ito[7],
        "PD-1/CTLA-4" = Okabe_Ito[6]
      )
    )
    rownames(a_col) <- pb_meta$donor[ix]
    labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
      sprintf("%s%s", x[1], x[2])
    })
    # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
    # mat_o <- seriation::seriate(mat, method = "PCA")
    # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
    heatmap_file <- glue("{my_out_dir}/heatmap_case-30_cluster-{this_cluster}_rowscale.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      # color = rev(scico::scico(palette = "davos", n = 20)),
      scale = "row",
      color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
    heatmap_file <- glue("{my_out_dir}/heatmap_case-30_cluster-{this_cluster}.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      color = rev(scico::scico(palette = "davos", n = 20)),
      # scale = "row",
      # color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )

  # Molly grant
  my_ens <- names(rowname_key[
    # which(res$Gene %in% c("IFNG", "IL26", "TNF", "IL17A", "GZMB"))
    # which(res$Gene %in% strsplit("CD274, ISG20, STAT1, GBP5, CXCL9, CXCL10, CXCL11", ", ")[[1]])
    which(rowname_key %in% strsplit("CXCL9, CXCL10, CXCL11, HLA-DRA, CD274, GBP4, IDO1, ISG15, ISG20", ", ")[[1]])
  ])
  mat <- as.matrix(pb[my_ens,])
  stopifnot(all(colnames(mat) == pb_meta$id))
  a_col <- as.data.frame(pb_meta[,c("leiden", "class", "case")])
  # a_col$facs_sorting <- as.character(a_col$facs_sorting)
  a_colors <- list(
    leiden = mpn65,
    # donor = tail(mpn65, length(unique(pb_meta$donor))),
    facs_sorting = c(
      "TRUE"  = "grey60",
      "FALSE" = "white"
    ),
    case = c(
      "Case"    = Okabe_Ito[1],
      "Control" = Okabe_Ito[2]
    ),
    class = c(
      "irColitis"             = Okabe_Ito[1],
      "On ICI therapy"        = Okabe_Ito[4],
      "Chronic diarrhea"      = Okabe_Ito[7],
      "Screening colonoscopy" = Okabe_Ito[3]
    ),
    drug = c(
      "None"        = Okabe_Ito[8],
      "CTLA-4"      = Okabe_Ito[4],
      "PD-1"        = Okabe_Ito[7],
      "PD-1/CTLA-4" = Okabe_Ito[6]
    )
  )
  # names(a_colors$donor) <- naturalfactor(unique(pb_meta$donor))
  names(a_colors$leiden) <- naturalfactor(unique(pb_meta$leiden))
  rownames(a_col) <- pb_meta$id
  labels_col <- colnames(mat)
  # labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
  #   sprintf("%s%s", x[1], x[2])
  # })
  # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
  # mat_o <- seriation::seriate(mat, method = "PCA")
  # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
  heatmap_file <- glue("{out_dir}/pb_cluster/heatmap-class_cluster-rowscale.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.1 + 4,
    height = nrow(mat) * 0.1 + 2,
    mat = mat,
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    hclust_method = "average",
    show_colnames = FALSE,
    # color = rev(scico::scico(palette = "davos", n = 20)),
    scale = "row",
    color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    # labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )
  heatmap_file <- glue("{out_dir}/pb_cluster/heatmap-class_cluster.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.1 + 4,
    height = nrow(mat) * 0.1 + 2,
    mat = mat,
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    hclust_method = "average",
    show_colnames = FALSE,
    color = rev(scico::scico(palette = "davos", n = 20)),
    # scale = "row",
    # color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    # labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )

}

# Case, class, drug pseudobulk by cluster

if (do_pb && n_donors > 1) {

  y <- with(a2$obs, model.matrix(~ 0 + donor:factor(leiden)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/")
  pb <- as(a2$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "^donor", "")
  colnames(pb) <- str_replace(colnames(pb), "factor\\(leiden\\)", "")
  colnames(pb)
  #
  pb_meta <- a2$obs %>% dplyr::select(
    case, class, class2, class_short,
    donor, drug
  ) %>% unique
  pb_meta <- inner_join(
    data.frame(
      donor = str_split_fixed(colnames(pb), ":", 2)[,1],
      leiden = str_split_fixed(colnames(pb), ":", 2)[,2]
    ),
    pb_meta,
    by = "donor"
  )
  pb_meta <- as_tibble(pb_meta)
  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
  all(glue("{pb_meta$donor}:{pb_meta$leiden}") == colnames(pb))

  x <- table(a2$obs$leiden)
  min_percent <- 100 * min(x) * 0.25 / sum(x)
  keep_ens <- a2$counts_stats$gene[a2$counts_stats$percent >= min_percent]

  for (this_cluster in unique(pb_meta$leiden)) {
    # this_cluster <- 1
    ix <- pb_meta$leiden == this_cluster
    library(limma)
    des1 <- with(pb_meta[ix, ], model.matrix(
      ~ case #+ chemistry + qubit_library_quantification_ng_ul
    ))
    object <- pb[keep_ens,ix]
    # object <- as.matrix(object[rowMeans(object) > 0.5,])
    fit1 <- lmFit(object = object, design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
    res <- topTable(fit1, coef = 2, number = 1e6, confint = TRUE)
    colnames(res)[1] <- "Gene"
    res$ensembl_id <- rownames(res)
    my_out_dir <- glue("{out_dir}/pb_cluster/case-control")
    dir.create(my_out_dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(
      res, glue("{my_out_dir}/de_case-vs-control_cluster-{this_cluster}.tsv.gz"),
      sep = "\t"
    )
    #
    p <- plot_limma_volcano(res) +
      labs(title = glue("Case vs Control - Cluster {this_cluster}"))
    p
    my_ggsave(
      glue("volcano-case_cluster-{this_cluster}"),
      out_dir = my_out_dir,
      plot = p,
      type = "pdf",
      scale = 1, width = 8, height = 6, units = "in", dpi = 300
    )
    # my_ens <- res$ensembl_id[res$adj.P.Val < 0.05]
    # length(my_ens)
    my_ens <- res$ensembl_id[1:40]
    mat <- as.matrix(pb[my_ens,ix])
    all(colnames(mat) == glue("{pb_meta$donor[ix]}:{this_cluster}"))
    colnames(mat) <- str_split_fixed(colnames(mat), ":", 2)[,1]
    a_col <- as.data.frame(pb_meta[ix,c("class", "case", "drug")])
    # a_col$facs_sorting <- as.character(a_col$facs_sorting)
    a_colors <- list(
      facs_sorting = c(
        "TRUE"  = "grey60",
        "FALSE" = "white"
      ),
      case = c(
        "Case"    = Okabe_Ito[1],
        "Control" = Okabe_Ito[2]
      ),
      class = c(
        "irColitis"             = Okabe_Ito[1],
        "On ICI therapy"        = Okabe_Ito[4],
        "Chronic diarrhea"      = Okabe_Ito[7],
        "Screening colonoscopy" = Okabe_Ito[3]
      ),
      drug = c(
        "None"        = Okabe_Ito[8],
        "CTLA-4"      = Okabe_Ito[4],
        "PD-1"        = Okabe_Ito[7],
        "PD-1/CTLA-4" = Okabe_Ito[6]
      )
    )
    rownames(a_col) <- pb_meta$donor[ix]
    labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
      sprintf("%s%s", x[1], x[2])
    })
    # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
    # mat_o <- seriation::seriate(mat, method = "PCA")
    # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
    heatmap_file <- glue("{my_out_dir}/heatmap_case-30_cluster-{this_cluster}_rowscale.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      # color = rev(scico::scico(palette = "davos", n = 20)),
      scale = "row",
      color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
    heatmap_file <- glue("{my_out_dir}/heatmap_case-30_cluster-{this_cluster}.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      color = rev(scico::scico(palette = "davos", n = 20)),
      # scale = "row",
      # color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
  }

  # Molly grant
  my_ens <- names(rowname_key[
    # which(res$Gene %in% c("IFNG", "IL26", "TNF", "IL17A", "GZMB"))
    # which(res$Gene %in% strsplit("CD274, ISG20, STAT1, GBP5, CXCL9, CXCL10, CXCL11", ", ")[[1]])
    which(rowname_key %in% strsplit("CXCL9, CXCL10, CXCL11, HLA-DRA, CD274, GBP4, IDO1, ISG15, ISG20", ", ")[[1]])
  ])
  mat <- as.matrix(pb[my_ens,])
  stopifnot(all(colnames(mat) == glue("{pb_meta$donor}:{pb_meta$leiden}")))
  a_col <- as.data.frame(pb_meta[,c("leiden", "class", "case", "drug")])
  # a_col$facs_sorting <- as.character(a_col$facs_sorting)
  a_colors <- list(
    leiden = mpn65,
    donor = tail(mpn65, length(unique(pb_meta$donor))),
    facs_sorting = c(
      "TRUE"  = "grey60",
      "FALSE" = "white"
    ),
    case = c(
      "Case"    = Okabe_Ito[1],
      "Control" = Okabe_Ito[2]
    ),
    class = c(
      "irColitis"             = Okabe_Ito[1],
      "On ICI therapy"        = Okabe_Ito[4],
      "Chronic diarrhea"      = Okabe_Ito[7],
      "Screening colonoscopy" = Okabe_Ito[3]
    ),
    drug = c(
      "None"        = Okabe_Ito[8],
      "CTLA-4"      = Okabe_Ito[4],
      "PD-1"        = Okabe_Ito[7],
      "PD-1/CTLA-4" = Okabe_Ito[6]
    )
  )
  names(a_colors$donor) <- naturalfactor(unique(pb_meta$donor))
  names(a_colors$leiden) <- naturalfactor(unique(pb_meta$leiden))
  rownames(a_col) <- glue("{pb_meta$donor}:{pb_meta$leiden}")
  labels_col <- colnames(mat)
  # labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
  #   sprintf("%s%s", x[1], x[2])
  # })
  # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
  # mat_o <- seriation::seriate(mat, method = "PCA")
  # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
  heatmap_file <- glue("{out_dir}/pb_cluster/heatmap-donor_cluster-rowscale.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.01 + 4,
    height = nrow(mat) * 0.1 + 2,
    mat = mat,
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    hclust_method = "average",
    show_colnames = FALSE,
    # color = rev(scico::scico(palette = "davos", n = 20)),
    scale = "row",
    color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    # labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )
  heatmap_file <- glue("{out_dir}/pb_cluster/heatmap-donor_cluster.png")
  message(glue("Writing {heatmap_file}"))
  pheatmap::pheatmap(
    filename = heatmap_file,
    width = ncol(mat) * 0.01 + 4,
    height = nrow(mat) * 0.1 + 2,
    mat = mat,
    # cluster_col = FALSE,
    # cluster_row = FALSE,
    hclust_method = "average",
    show_colnames = FALSE,
    color = rev(scico::scico(palette = "davos", n = 20)),
    # scale = "row",
    # color = rev(scico::scico(palette = "roma", n = 20)),
    border_color = NA,
    labels_row = ensembl_to_symbol[rownames(mat)],
    # labels_col = labels_col,
    annotation_col = a_col,
    annotation_colors = a_colors
  )

  for (this_cluster in unique(pb_meta$leiden)) {
    # this_cluster <- 1
    ix <- pb_meta$leiden == this_cluster
    pb_meta$classf <- factor(
      pb_meta$class,
      c(
        "Screening colonoscopy",
        "Chronic diarrhea",
        "On ICI therapy",
        "irColitis"
      )
      # pb_meta$class_short,
      # c("SC", "CD", "ICI", "IRC")
    )
    pb_meta$class_shortf <- factor(
      pb_meta$class_short,
      c("SC", "CD", "ICI", "IRC")
    )
    des1 <- with(pb_meta[ix,], model.matrix(
      ~ classf # case # + chemistry#  + facs_sorting
    ))
    object <- pb[keep_ens,ix]
    object <- as.matrix(object[rowMeans(object) > 0.5,])
    fit1 <- lmFit(object = object, design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
    res <- topTable(fit1, number = 1e6, confint = TRUE)
    colnames(res)[1] <- "Gene"
    res$ensembl_id <- rownames(res)
    my_out_dir <- glue("{out_dir}/pb_cluster/class")
    dir.create(my_out_dir, showWarnings = FALSE)
    fwrite(
      res, glue("{my_out_dir}/de_class_cluster-{this_cluster}.tsv.gz"),
      sep = "\t"
    )
    #
    for (i in 2:4) {
      res <- topTable(fit1, coef = i, number = 1e6, confint = TRUE)
      colnames(res)[1] <- "Gene"
      res$ensembl_id <- rownames(res)
      p <- plot_limma_volcano(res, fdr = 0.1) +
        labs(title = glue("{levels(pb_meta$classf)[i]} vs Screening colonoscopy - Cluster {this_cluster}"))
      my_ggsave(
        glue("volcano_class-{levels(pb_meta$class_shortf)[i]}_cluster-{this_cluster}"),
        out_dir = my_out_dir,
        plot = p,
        type = "pdf",
        scale = 1, width = 8, height = 6, units = "in", dpi = 300
      )
    }
    # my_ens <- res$ensembl_id[res$adj.P.Val < 0.05]
    # length(my_ens)
    res <- topTable(fit1, number = 1e6, confint = TRUE)
    colnames(res)[1] <- "Gene"
    res$ensembl_id <- rownames(res)
    my_ens <- res$ensembl_id[1:40]
    mat <- as.matrix(pb[my_ens,ix])
    all(colnames(mat) == glue("{pb_meta$donor[ix]}:{this_cluster}"))
    colnames(mat) <- str_split_fixed(colnames(mat), ":", 2)[,1]
    a_col <- as.data.frame(pb_meta[ix,c("class", "case", "drug")])
    # a_col$facs_sorting <- as.character(a_col$facs_sorting)
    a_colors <- list(
      facs_sorting = c(
        "TRUE"  = "grey60",
        "FALSE" = "white"
      ),
      case = c(
        "Case"    = Okabe_Ito[1],
        "Control" = Okabe_Ito[2]
      ),
      class = c(
        "irColitis"             = Okabe_Ito[1],
        "On ICI therapy"        = Okabe_Ito[4],
        "Chronic diarrhea"      = Okabe_Ito[7],
        "Screening colonoscopy" = Okabe_Ito[3]
      ),
      drug = c(
        "None"        = Okabe_Ito[8],
        "CTLA-4"      = Okabe_Ito[4],
        "PD-1"        = Okabe_Ito[7],
        "PD-1/CTLA-4" = Okabe_Ito[6]
      )
    )
    rownames(a_col) <- pb_meta$donor[ix]
    labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
      sprintf("%s%s", x[1], x[2])
    })
    # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
    # mat_o <- seriation::seriate(mat, method = "PCA")
    # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
    heatmap_file <- glue("{my_out_dir}/heatmap_class_cluster-{this_cluster}_rowscale.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      # color = rev(scico::scico(palette = "davos", n = 20)),
      scale = "row",
      color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
    heatmap_file <- glue("{my_out_dir}/heatmap_class_cluster-{this_cluster}.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      color = rev(scico::scico(palette = "davos", n = 20)),
      # scale = "row",
      # color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
  }
  
  for (this_cluster in unique(pb_meta$leiden)) {
    # this_cluster <- 1
    ix <- pb_meta$leiden == this_cluster
    pb_meta$drug <- factor(pb_meta$drug,
      c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4")
    )
    # pb_meta$drug2 <- as.character(pb_meta$drug)
    # pb_meta$drug2[pb_meta$drug2 == "PD-1/CTLA-4"] <- "PD-1"
    # pb_meta$drug2 <- factor(pb_meta$drug2,
    #   c("None", "CTLA-4", "PD-1")
    # )
    des1 <- with(pb_meta[ix,], model.matrix(
      ~ drug
    ))
    object <- pb[keep_ens,ix]
    object <- as.matrix(object[rowMeans(object) > 0.5,])
    fit1 <- lmFit(object = object, design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
    my_out_dir <- glue("{out_dir}/pb_cluster/drug")
    dir.create(my_out_dir, showWarnings = FALSE)
    #
    for (i in 2:4) {
      res <- topTable(fit1, coef = i, number = 1e6, confint = TRUE)
      colnames(res)[1] <- "Gene"
      res$ensembl_id <- rownames(res)
      p <- plot_limma_volcano(res, fdr = 0.1) +
        labs(title = glue("{levels(pb_meta$drug)[i]} vs None - Cluster {this_cluster}"))
      my_drug <- str_replace_all(
        janitor::make_clean_names(levels(pb_meta$drug)[i], case = "none"),
        "_", ""
      )
      fwrite(
        res, glue("{my_out_dir}/de_drug-{my_drug}_cluster-{this_cluster}.tsv.gz"),
        sep = "\t"
      )
      my_ggsave(
        glue("volcano_drug-{my_drug}_cluster-{this_cluster}"),
        out_dir = my_out_dir,
        plot = p,
        type = "pdf",
        scale = 1, width = 8, height = 6, units = "in", dpi = 300
      )
    }
    # my_ens <- res$ensembl_id[res$adj.P.Val < 0.05]
    # length(my_ens)
    res <- topTable(fit1, number = 1e6, confint = TRUE)
    colnames(res)[1] <- "Gene"
    res$ensembl_id <- rownames(res)
    my_ens <- res$ensembl_id[1:40]
    mat <- as.matrix(pb[my_ens,ix])
    all(colnames(mat) == glue("{pb_meta$donor[ix]}:{this_cluster}"))
    colnames(mat) <- str_split_fixed(colnames(mat), ":", 2)[,1]
    a_col <- as.data.frame(pb_meta[ix,c("class", "case", "drug")])
    # a_col$facs_sorting <- as.character(a_col$facs_sorting)
    a_colors <- list(
      facs_sorting = c(
        "TRUE"  = "grey60",
        "FALSE" = "white"
      ),
      case = c(
        "Case"    = Okabe_Ito[1],
        "Control" = Okabe_Ito[2]
      ),
      class = c(
        "irColitis"             = Okabe_Ito[1],
        "On ICI therapy"        = Okabe_Ito[4],
        "Chronic diarrhea"      = Okabe_Ito[7],
        "Screening colonoscopy" = Okabe_Ito[3]
      ),
      drug = c(
        "None"        = Okabe_Ito[8],
        "CTLA-4"      = Okabe_Ito[4],
        "PD-1"        = Okabe_Ito[7],
        "PD-1/CTLA-4" = Okabe_Ito[6]
      )
    )
    rownames(a_col) <- pb_meta$donor[ix]
    labels_col <- sapply(str_split(colnames(mat), "_"), function(x) {
      sprintf("%s%s", x[1], x[2])
    })
    # mat_o <- seriation::seriate(mat, method = "BEA_TSP")
    # mat_o <- seriation::seriate(mat, method = "PCA")
    # mat <- mat[rev(mat_o[[1]]), mat_o[[2]]]
    heatmap_file <- glue("{my_out_dir}/heatmap_drug_cluster-{this_cluster}_rowscale.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      # color = rev(scico::scico(palette = "davos", n = 20)),
      scale = "row",
      color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
    heatmap_file <- glue("{my_out_dir}/heatmap_drug_cluster-{this_cluster}.png")
    message(glue("Writing {heatmap_file}"))
    pheatmap::pheatmap(
      filename = heatmap_file,
      width = ncol(mat) * 0.15 + 4,
      height = nrow(mat) * 0.15 + 1,
      mat = mat,
      # cluster_col = FALSE,
      # cluster_row = FALSE,
      hclust_method = "average",
      # show_colnames = FALSE,
      color = rev(scico::scico(palette = "davos", n = 20)),
      # scale = "row",
      # color = rev(scico::scico(palette = "roma", n = 20)),
      border_color = NA,
      labels_row = ensembl_to_symbol[rownames(mat)],
      labels_col = labels_col,
      annotation_col = a_col,
      annotation_colors = a_colors
    )
  }

  # # pb_fit <- do_loess(pb[,ix], exclude_genes, loess_span = 0.3, min_percent = 0.05)
  # pb_fit <- do_loess(pb[,], NULL, loess_span = 0.3, min_percent = 0.05)
  # keep_genes <- which(pb_fit$counts_stats$rank < 5000)
  # pb_pca <- RSpectra::svds(
  #   A    = t(pb[keep_genes,]),
  #   k    = 10,
  #   opts = list(
  #     center = TRUE, scale = TRUE, maxitr = 2000, tol = 1e-10
  #   )
  # )
  # pb_pca$genes <- rownames(pb)[keep_genes]
  # pb_meta$PC1 <- pb_pca$u[,1]
  # pb_meta$PC2 <- pb_pca$u[,2]
  # pb_meta$PC3 <- pb_pca$u[,3]
  # pb_meta$PC4 <- pb_pca$u[,4]
  # pb_meta$PC5 <- pb_pca$u[,5]
  # pb_meta$PC6 <- pb_pca$u[,6]
  # ggplot(pb_meta) +
  #   aes(PC1, PC2, fill = class) +
  #   geom_point(shape = 21, size = 5)
  # ggplot(pb_meta) +
  #   aes(PC3, PC4, fill = class) +
  #   geom_point(shape = 21, size = 5)
  # ggplot(pb_meta) +
  #   aes(PC5, PC6, fill = class) +
  #   geom_point(shape = 21, size = 5)

}

# All vs All (ava)
########################################################################

# x <- de_ava %>%
#   filter(P.Value < 0.05 / nrow(de_ava)) %>%
#   filter(logFC > log2(2)) %>%
#   group_by(coef) %>%
#   summarize(n = n())
# x$group1 <- as.integer(str_split_fixed(x$coef, " vs ", 2)[,1])
# x$group2 <- as.integer(str_split_fixed(x$coef, " vs ", 2)[,2])
# x$coef <- NULL
# y <- de_ava %>%
#   filter(P.Value < 0.05 / nrow(de_ava)) %>%
#   filter(logFC < log2(1/2)) %>%
#   group_by(coef) %>%
#   summarize(n = n())
# y$group2 <- as.integer(str_split_fixed(y$coef, " vs ", 2)[,1])
# y$group1 <- as.integer(str_split_fixed(y$coef, " vs ", 2)[,2])
# y$coef <- NULL
# x <- rbind(
#   x[,c("n", "group1", "group2")],
#   y[,c("n", "group1", "group2")]
# )
# x <- reshape2::dcast(x, group1 ~ group2, value.var = "n")
# rownames(x) <- x$group1
# x$group1 <- NULL

if ("de_ava" %in% ls() && n_donors > 1) {

x <- de_ava %>%
  # filter(P.Value < 0.05 / nrow(de_ava)) %>%
  filter(P.Value < 0.05 / S4Vectors::unname(table(de_ava$coef)[1])) %>%
  filter(abs(logFC) > log2(2)) %>%
  group_by(coef) %>%
  summarize(n = n())
x$group1 <- as.character(str_split_fixed(x$coef, " vs ", 2)[,1])
x$group2 <- as.character(str_split_fixed(x$coef, " vs ", 2)[,2])
x$coef <- NULL
y <- x
colnames(y) <- c("n", "group2", "group1")
x <- rbind(
  x[,c("n", "group1", "group2")],
  y[,c("n", "group1", "group2")]
)
x <- reshape2::dcast(x, group1 ~ group2, value.var = "n")
rownames(x) <- x$group1
x$group1 <- NULL

x[is.na(x)] <- 0
x <- as.dist(x)

hc <- dendsort::dendsort(hclust(x, method = "complete"))

d <- tidygraph::as_tbl_graph(hc)
p <- ggraph(d, layout = 'dendrogram', height = height, circular = FALSE) + 
  geom_edge_elbow() +
  geom_node_point(
    mapping = aes(
      filter = leaf,
      x      = x,
      y      = y - max(hc$height) * 0.05,
      colour = label
    ),
    size = 8
  ) +
  geom_node_text(
    mapping = aes(
      x      = x,
      y      = y - max(hc$height) * 0.015,
      filter = leaf,
      label  = label
    ),
    size = 5, color = "black"
  ) +
  scale_y_reverse() +
  # scale_y_continuous(trans = trans_reverser(sqrt_trans())) +
  scale_colour_manual(values = cluster_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin     = unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) +
  coord_flip() 
my_ggsave(
  slug = "cluster-hclust-pseudobulk-de",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  plot = p,
  # plot = wrap_plots(p[1:9], ncol = 1),
  type = "pdf",
  limitsize = FALSE,
  scale = 1,
  width = 10,
  height = 0.5 * ncol(as.matrix(x)),
  units = "in", dpi = 300
)
}

if ("de_ava" %in% ls() && n_donors > 1) {
  for (this_coef in unique(de_ava$coef)) {
    top1 <- de_ava %>%
      filter(coef == this_coef) %>%
      # mutate(g1 = str_split_fixed(coef, " vs ", 2)[,1]) %>%
      # mutate(g2 = str_split_fixed(coef, " vs ", 2)[,2]) %>%
      # filter(g1 == 1 | g2 == 1) %>%
      # filter(g1 == 7 | g2 == 7) %>%
      mutate(Gene = ID)
    p <- plot_limma_volcano(top1) +
      labs(title = this_coef)
    my_ggsave(
      glue("volcano-{str_replace_all(this_coef, ' ', '_')}"),
      #out_dir = glue("figures/{analysis_name}/cluster_volcano"),
      out_dir = glue("{out_dir}/cluster_volcano/ava"),
      plot = p,
      type = "pdf",
      scale = 1, width = 8, height = 6, units = "in", dpi = 300
    )
  }
}

gene_percents <- rowSums(a2$counts > 0) / ncol(a2$counts)
gene_percents <- data.frame(
  ensembl_id = names(gene_percents),
  percent = gene_percents * 100
)
de_ova <- left_join(de_ova, gene_percents, by = "ensembl_id")

# UMAP with top markers from pseudobulk
########################################################################
print_status("skipping huge umap figure: umap-top_marker-pseudobulk.pdf")
if (TRUE) {
  if (n_donors > 1) {
  print_status("UMAP with top markers from pseudobulk")
  # p1 <- plot_hexmix(
  #   x = a2$obs$UMAP1,
  #   y = a2$obs$UMAP2,
  #   group = a2$obs$leiden,
  #   group_colors = cluster_colors,
  #   bins = 301
  # ) +
  # labs(
  #   title = glue(
  #     "{length(unique(a2$obs$leiden))} clusters of {comma(nrow(a2$obs))} cells from {length(unique(a2$obs$donor))} donors"
  #   )
  # )
  de_ova_top <- de_ova %>%
    group_by(coef) %>%
    filter(percent > 1) %>%
    filter(logFC > 0) %>%
    filter(!ensembl_id %in% exclude_genes) %>%
    # top_n(n = 200, wt = -log10(P.Value)) %>%
    mutate(cluster = as.integer(str_split_fixed(coef, " ", 2)[,1])) %>%
    arrange(P.Value, cluster)
  # de_ova_top <- de_ova_top[!duplicated(de_ova_top$ensembl_id),]
  table(de_ova_top$cluster)
  de_ova_top <- de_ova_top %>%
    group_by(coef) %>%
    top_n(n = 5, wt = -log10(P.Value))
  table(de_ova_top$cluster)
  #
  plots <- lapply(seq(nrow(de_ova_top)), function(i) {
    this_gene    <- de_ova_top$ensembl_id[i]
    this_cluster <- de_ova_top$cluster[i]
    this_logfc   <- de_ova_top$logFC[i]
    this_pval    <- de_ova_top$P.Value[i]
    plot_hexgene(
      x            = a2$obs$UMAP1,
      y            = a2$obs$UMAP2,
      z            = as.numeric(a2$log2cpm[this_gene,]),
      # bins         = 101,
      bins         = 51,
      palette      = "davos",
      direction    = -1,
      use_quantile = TRUE
    ) +
    labs(
      # title = glue("{rowname_key[this_gene]} ({this_cluster})")
      title = rowname_key[this_gene],
      subtitle = parse(text = glue(
        '"P =" ~ {scientific_10(this_pval)} ~ ", FC =" ~ {signif(2^this_logfc, 2)}'
      ))
    ) +
    theme(legend.position = "none")
  })
  plots <- split(plots, de_ova_top$cluster)
  #
  p <- lapply(seq_along(plots), function(i) {
    my_group <- names(plots)[i]
    z <- as.integer(a2$obs$leiden == my_group)
    p <- suppressMessages(
      plot_hexgene(
        x = a2$obs$UMAP1,
        y = a2$obs$UMAP2,
        z = z,
        # bins = 101
        bins = 51
      ) +
      scale_fill_gradientn(
        colors = scico::scico(
          n = 20, direction = 1, palette = "grayC"
        )[3:20],
        labels = NULL
      ) +
      labs(title = glue("Cluster {my_group} (n={comma(sum(z))})")) +
      theme(plot.title = element_text(face = "italic")) +
      guides(fill = "none")
    )
    res <- wrap_plots(
      c(list(p), plots[[i]]),
      ncol = max(1, 1 + length(plots[[i]]))
    )
  })
  #
  x <- wrap_plots(p, ncol = 1)
  my_ggsave(
    slug = "umap-top_marker-pseudobulk",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    plot = x,
    # plot = wrap_plots(p[1:9], ncol = 1),
    type = "pdf",
    limitsize = FALSE,
    scale = 1.8,
    width = 6 * 1.8,
    height = 1.5 * length(unique(de_ova_top$cluster)),
    units = "in", dpi = 300
  )
  print_status("done writing.")
  }
}


#ncol <- which.min(abs(sapply(1:16, function(i) {
#  i - (16 / 9) * (length(plots) / i)
#})))
#p <- wrap_plots(c(list(p1), plots), ncol = ncol)
#my_ggsave(
#  slug = "umap-top_marker-pseudobulk",
#  #out_dir = glue("figures/{analysis_name}"),
#  out_dir = out_dir,
#  plot = p,
#  type = "pdf",
#  scale = 1.8, width = 16, height = 9, units = "in", dpi = 300
#)

# my_gene <- names(ensembl_to_symbol[ensembl_to_symbol == "IGKV1D-12"])
#   plot_hexgene(
#     x = a2$obs$UMAP1,
#     y = a2$obs$UMAP2,
#     z = as.numeric(a2$log2cpm[my_gene,]),
#     bins = 101,
#     palette = "davos",
#     direction = -1,
#     use_quantile = TRUE
#   )


# Heatmap of best markers for each cluster
########################################################################
# Top markers for each cluster
if ("de_ova" %in% ls() && n_donors > 1) {
  x <- de_ova %>%
    dplyr::mutate(group = str_split_fixed(coef, " vs ", 2)[,1]) %>%
    dplyr::group_by(group) %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::top_n(n = 200, wt = -log10(P.Value))
  x <- x[!duplicated(x$ensembl_id),]
  x <- x %>%
    group_by(group) %>%
    top_n(n = 5, wt = -log10(P.Value))
  table(x$group)
    # dplyr::top_n(n = 5, wt = -log10(P.Value))
  these_genes <- x$ensembl_id
  fig_height <- length(these_genes) * 0.2
  fig_width <- length(unique(de_ova$coef)) * 0.4 + 3
  #
  # Subset the table
  x <- de_ova %>% dplyr::filter(ensembl_id %in% these_genes) %>%
    dplyr::mutate(group = str_split_fixed(coef, " vs ", 2)[,1])
  x <- as.data.frame(dcast.data.table(
    data = as.data.table(x), formula = ID ~ group, value.var = "logFC"
  ))
  rownames(x) <- x$ID
  x$ID <- NULL
  # Order the heatmap
  set.seed(1)
  xo <- seriation::seriate(as.matrix(x) + abs(min(x)), method = "BEA_TSP")
  xd <- de_ova %>% dplyr::filter(ensembl_id %in% these_genes) %>%
    dplyr::mutate(group = str_split_fixed(coef, " vs ", 2)[,1])
  xd$ID <- factor(xd$ID, rownames(x)[xo[[1]]])
  xd$group <- factor(xd$group, colnames(x)[xo[[2]]])
  # Plot
  logfc_limits <- range(xd$logFC)
  logfc_limits[1] <- min(-logfc_limits[2], logfc_limits[1])
  logfc_limits[2] <- max(-logfc_limits[1], logfc_limits[2])
  p2 <- ggplot(xd) +
  geom_tile(
    aes(y = ID, x = group, fill = logFC)
  ) +
  scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  # scale_fill_viridis_c(
  #   name = "AUC", guide = guide_colorbar(barheight = 40), breaks = pretty_breaks(10),
  #   limits = c(0, 1)
  # ) +
  scale_fill_scico(
    # palette = "vik", limits = c(0, 1), direction = 1,
    palette = "vik", limits = logfc_limits, direction = 1,
    name = "Log2FC", guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(9)
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "italic")
  #   axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )
  p1 <- ggplot(xd) +
  geom_tile(
    aes(y = 1, x = group, fill = group)
  ) +
  scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_fill_manual(
    values = cluster_colors, guide = "none"
  )
  p <- (
    p1 + theme(plot.margin = margin(b = 0))
  ) / (
    p2 + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0)
    )
  ) + plot_layout(heights = c(1, length(these_genes)))
  #p
  my_ggsave(
    slug = "heatmap-top-markers-pseudobulk",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    limitsize = FALSE,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )
}


# Top markers for each group of cells
########################################################################
if ("de_ova" %in% ls() && n_donors > 1) {
  for (my_group in sort(unique(a2$obs$leiden))) {
    # pdf_file <- glue("figures/{analysis_name}/cluster_umap/umap-cluster-{my_group}.pdf")
    # if (file.exists(pdf_file)) {
    #   print_status(glue("File exists '{pdf_file}'"))
    #   next
    # }
    # Select top genes for this group
    x1 <- a2$de %>%
      dplyr::filter(group == my_group) %>%
      dplyr::filter(logFC > 0) %>%
      dplyr::mutate(
        rank1 = rank(100 * auc - pct_out),
        rank2 = rank(logFC * (pct_in - pct_out))
      ) %>%
      dplyr::arrange(-auc) %>%
      head(14)
      # top_n(n = 31, wt = (rank1 + rank2))
      # dplyr::top_n(n = 31, wt = auc)
    my_genes1 <- x1$feature
    x2 <- a2$de %>%
      dplyr::filter(group == my_group) %>%
      dplyr::filter(logFC > 0) %>%
      dplyr::mutate(
        rank1 = rank(100 * auc - pct_out),
        rank2 = rank(logFC * (pct_in - pct_out))
      ) %>%
      dplyr::top_n(n = 14, wt = (rank1 + rank2))
    my_genes2 <- x2$feature
    de_ova$exclude <- de_ova$ensembl_id %in% exclude_genes
    x3 <- de_ova %>%
      dplyr::filter(str_detect(coef, glue("^{my_group} vs all$"))) %>%
      dplyr::filter(!exclude) %>%
      dplyr::filter(logFC > 0) %>%
      dplyr::top_n(n = 14, wt = logFC * -log10(P.Value))
      # dplyr::top_n(n = 14, wt = -log10(P.Value))
    my_genes3 <- x3$ensembl_id
    gene_lists <- list(
      # "auc"        = my_genes1,
      # "rank"       = my_genes2,
      "pseudobulk" = my_genes3
    )
    z <- as.integer(a2$obs$leiden == my_group)
    p1 <- suppressMessages(
        plot_hexgene(
        x = a2$obs$UMAP1,
        y = a2$obs$UMAP2,
        z = z,
        bins = 91
      ) +
      scale_fill_gradientn(
        colors = scico::scico(
          n = 20, direction = 1, palette = "grayC"
        )[3:20],
        labels = NULL
      ) +
      labs(title = glue("Cluster {my_group} (n={comma(sum(z))})")) +
      theme(plot.title = element_text(face = "italic")) +
      guides(fill = FALSE)
    )
    # guides(fill = guide_colorbar(title = NULL, barwidth = 5, ticks = FALSE))
    for (i in names(gene_lists)) {
      plots <- lapply(gene_lists[[i]], function(this_gene) {
        p <- plot_hexgene(
          x = a2$obs$UMAP1,
          y = a2$obs$UMAP2,
          z = as.numeric(a2$log2cpm[this_gene,]),
          bins = 91,
          palette = "davos",
          direction = -1
        ) + labs(title = rowname_key[this_gene]) +
        theme(legend.position = "none")
        if (i %in% c("auc", "rank")) {
          this_auc <- a2$de$auc[a2$de$group == my_group & a2$de$feature == this_gene]
          # p <- p + annotate(
          #   geom = "text",
          #   x = Inf, y = Inf,
          #   hjust = 1.05, vjust = 1.5,
          #   label = glue("{round(this_auc, 2)} AUC"),
          #   size = 6
          # )
          p <- p + labs(subtitle = glue("AUROC = {signif(this_auc, 2)}"))
        } else {
          this_pval  <- x3$P.Value[x3$ensembl_id == this_gene]
          this_logfc <- x3$logFC[x3$ensembl_id == this_gene]
          # p <- p + annotate(
          #   geom = "text",
          #   x = Inf, y = Inf,
          #   hjust = 1.05, vjust = 1.5,
          #   label = glue("P = {signif(this_pval, 2)}"),
          #   size = 6
          # )
          # p <- p + labs(subtitle = glue("P = {signif(this_pval, 2)}"))
          p <- p + labs(
            subtitle = parse(text = glue(
              '"P =" ~ {scientific_10(this_pval)} ~ ", FC =" ~ {signif(2^this_logfc, 2)}'
            ))
          ) 
        }
        return(p)
      })
      ncol <- which.min(abs(sapply(1:16, function(i) {
        i - (16 / 9) * (length(plots) / i)
      })))
      p <- wrap_plots(c(list(p1), plots), ncol = ncol)
      my_ggsave(
        slug = glue("umap-cluster-{i}-{my_group}"),
        #out_dir = glue("figures/{analysis_name}/cluster_umap"),
        out_dir = glue("{out_dir}/cluster_umap"),
        plot = p,
        type = "pdf",
        scale = 1.5, width = 16, height = 9, units = "in", dpi = 300
      )
    }
  }
}
# UMAP of best AUC markers for each cluster

## PCA coordinates displayed on UMAP
#########################################################################
#for (i in seq_along(colnames(a2$pca_h))) {

#  my_pc <- colnames(a2$pca_h)[i]

#  p1 <- suppressMessages(
#    plot_hexgene(
#      x = a2$obs$UMAP1,
#      y = a2$obs$UMAP2,
#      z = a2$obs[[my_pc]],
#      bins = 201
#    ) +
#    scale_fill_gradientn(
#      colors = scico::scico(
#        n = 20, direction = -1, palette = "roma"
#      )[1:20],
#      labels = NULL,
#      name = my_pc
#    ) +
#    # guides(fill = guide_colorbar(barwidth = 7))
#    labs(title = my_pc) +
#    guides(fill = "none") +
#    theme(plot.title = element_text(face = "plain"))
#  )
#  p1

#  ix <- which(rank(-a2$pca$v[,2]) < 10)
#  my_ens <- rownames(a2$counts)[a2$ix_genes[ix]]
#  p1 <- suppressMessages(
#    plot_hexgene(
#      x = a2$obs$UMAP1,
#      y = a2$obs$UMAP2,
#      z = a2$log2cpm[my_ens[1],],
#      bins = 101
#    ) +
#    scale_fill_gradientn(
#      colors = scico::scico(
#        n = 20, direction = -1, palette = "davos"
#      )[3:20],
#      labels = NULL,
#      name = my_pc
#    ) +
#    # guides(fill = guide_colorbar(barwidth = 7))
#    # labs(title = rowname_key[my_ens[1]]) +
#    guides(fill = "none")
#  )
#  p1

#  # guides(fill = guide_colorbar(title = NULL, barwidth = 5, ticks = FALSE))
#  for (i in names(gene_lists)) {
#    plots <- lapply(gene_lists[[i]], function(this_gene) {
#      p <- plot_hexgene(
#        x = a2$obs$UMAP1,
#        y = a2$obs$UMAP2,
#        z = as.numeric(a2$log2cpm[this_gene,]),
#        bins = 91,
#        palette = "davos",
#        direction = -1
#      ) + labs(title = rowname_key[this_gene]) +
#      theme(legend.position = "none")
#      if (i %in% c("auc", "rank")) {
#        this_auc <- a2$de$auc[a2$de$group == my_group & a2$de$feature == this_gene]
#        # p <- p + annotate(
#        #   geom = "text",
#        #   x = Inf, y = Inf,
#        #   hjust = 1.05, vjust = 1.5,
#        #   label = glue("{round(this_auc, 2)} AUC"),
#        #   size = 6
#        # )
#        p <- p + labs(subtitle = glue("AUROC = {signif(this_auc, 2)}"))
#      } else {
#        this_pval  <- x3$P.Value[x3$ensembl_id == this_gene]
#        this_logfc <- x3$logFC[x3$ensembl_id == this_gene]
#        # p <- p + annotate(
#        #   geom = "text",
#        #   x = Inf, y = Inf,
#        #   hjust = 1.05, vjust = 1.5,
#        #   label = glue("P = {signif(this_pval, 2)}"),
#        #   size = 6
#        # )
#        # p <- p + labs(subtitle = glue("P = {signif(this_pval, 2)}"))
#        p <- p + labs(
#          subtitle = parse(text = glue(
#            '"P =" ~ {scientific_10(this_pval)} ~ ", FC =" ~ {signif(2^this_logfc, 2)}'
#          ))
#        ) 
#      }
#      return(p)
#    })
#    ncol <- which.min(abs(sapply(1:16, function(i) {
#      i - (16 / 9) * (length(plots) / i)
#    })))
#    p <- wrap_plots(c(list(p1), plots), ncol = ncol)
#    my_ggsave(
#      slug = glue("umap-cluster-{i}-{my_group}"),
#      #out_dir = glue("figures/{analysis_name}/cluster_umap"),
#      out_dir = glue("{out_dir}/cluster_umap"),
#      plot = p,
#      type = "pdf",
#      scale = 1.5, width = 16, height = 9, units = "in", dpi = 300
#    )
#  }
#}

# # SingleR
# ########################################################################
# 
# dice <- SingleR::DatabaseImmuneCellExpressionData()
# 
# y <- model.matrix(~ 0 + factor(leiden), data = a2$obs)
# y <- sweep(y, 2, colSums(y), "/")
# z <- a2$log2cpm %*% y
# rownames(z) <- rowname_key[rownames(z)]
# colnames(z) <- 1:ncol(z)
# res <- SingleR::SingleR(
#   test      = z,
#   ref       = dice,
#   labels    = dice$label.fine,
#   de.method = "wilcox"
# )
# 
# SingleR::plotScoreHeatmap(
#   results       = res,
#   show_colnames = TRUE,
#   fontsize      = 20,
#   filename      = glue("figures/{analysis_name}/dice.pdf"),
#   width         = 6 + length(unique(a2$obs$leiden)) * 0.6,
#   height        = 6
# )


# Composition analysis of TCR
########################################################################

if (n_donors > 1) {
  if (
    all(c("case", "donor", "leiden", "TRAV") %in% colnames(a2$obs)) &&
    sum(is.na(a2$obs$TRAV)) < 0.8 * nrow(a2$obs)
  ) {
    d <- a2$obs %>%
      dplyr::filter(!is.na(TRAV)) %>%
      dplyr::select(case, donor, TRAV) %>%
      # dplyr::mutate(case = class2) %>%
      dplyr::group_by(case, donor, TRAV) %>%
      dplyr::summarize(n = n()) %>%
      dplyr::mutate(freq = 100 * n / sum(n))
    d_t <- rbindlist(lapply(unique(d$TRAV), function(this_TRAV) {
      x <- d %>% dplyr::filter(TRAV == this_TRAV)
      if (length(table(x$case)) != 2) {
        return(
          tibble(
            statistic = NA, p.value = 1, method = NA,
            alternative = NA, TRAV = this_TRAV
          )
        )
      }
      # x <- broom::tidy(t.test(n ~ case, x))
      # x <- broom::tidy(t.test(freq ~ case, x))
      x <- broom::tidy(wilcox.test(log(freq) ~ case, x))
      x$TRAV <- this_TRAV
      return(x)
    }))
    d_label <- d %>%
      dplyr::group_by(case) %>%
      summarise(n = sum(n)) %>%
      dplyr::mutate(label = glue("{case} (n = {comma(n)})"))
    d_t <- inner_join(
      d_t,
      d %>% dplyr::group_by(TRAV) %>% dplyr::mutate(y = mean(freq) + 4) %>% select(TRAV, y) %>% unique,
      by = "TRAV"
    )
    scientific_10 <- function(x) {
      ifelse(x < 0.01,
        gsub("e", "%*%10^", scales::scientific_format(digits = 1)(x)),
        signif(x, 1)
      )
    }
    d_t$TRAV <- factor(d_t$TRAV, d_t$TRAV[order(-d_t$p.value)])
    d$TRAV <- factor(d$TRAV, levels(d_t$TRAV))
    p <- ggplot() +
      geom_segment(
        data = data.frame(
          y = c(10^(-1:2), 0.5 * 10^(-1:2))
        ),
        mapping = aes(
          x = y, xend = y, y = -Inf, yend = Inf
        ),
        size = 0.2, alpha = 0.3
      ) +
      scale_color_manual(
        values = c(NA, "black"), guide = "none"
      ) +
      geom_crossbar(
        data = d %>% dplyr::group_by(TRAV, case) %>%
          dplyr::summarize(
            ymin = quantile(freq, 0.75),
            y = quantile(freq, 0.5),
            ymax = quantile(freq, 0.25)
          ),
        mapping = aes(
          y = factor(TRAV), xmin = ymin, x = y, xmax = ymax,
          group = case, fill = case
        ),
        color = "grey20",
        width = 0.5,
        position = position_dodge(width = 0.9)
      ) +
      annotate(
        geom = "rect",
        xmin = 0, xmax = Inf,
        ymin = as.integer(sort(unique(d$TRAV))) - 0.5,
        ymax = as.integer(sort(unique(d$TRAV))) + 0.5,
        fill = rep(
          c(NA, "black"),
          length.out = length(unique(d$TRAV))
        ),
        alpha = 0.1
      ) +
      geom_crossbar(
        data = d %>% dplyr::group_by(TRAV, case) %>%
          dplyr::summarize(
            ymin = quantile(freq, 0.75),
            y = quantile(freq, 0.5),
            ymax = quantile(freq, 0.25)
          ),
        mapping = aes(
          y = factor(TRAV), xmin = ymin, x = y, xmax = ymax,
          group = case, fill = case
        ),
        color = "grey20",
        width = 0.5,
        position = position_dodge(width = 0.9)
      ) +
      geom_point(
        data = d,
        mapping = aes(y = TRAV, x = freq, fill = case),
        position = position_quasirandom(groupOnX = FALSE, dodge.width = 0.9),
        alpha = 0.4, shape = 21
      ) +
      geom_text(
        data = d_t,
        mapping = aes(
          y = TRAV,
          x = 85,
          # y = y + 2,
          label = ifelse(
            # p.value < 0.05 / 12,
            p.value < 0.5 / length(unique(a2$obs$TRAV)),
            scientific_10(p.value),
            ""
          )
        ),
        size = 5, parse = TRUE
      ) +
      scale_fill_manual(
        name = NULL, values = pals::okabe(3)[c(2,3)],
        labels = d_label$label
      ) +
      # scale_y_continuous(name = "Percent") +
      scale_x_log10(
        name = "Percent",
        #labels = scales::label_number()
        labels = function(x) signif(x, 3)
      ) +
      annotation_logticks(side = "b") +
      labs(y = NULL) +
      theme(
        # legend.position = c(1, 1),
        # legend.justification = c(1, 1),
        legend.position = "bottom",
        legend.background = element_blank()
        # panel.grid.major.y = element_line(size = 0.1, color = "#00000055"),
        # panel.grid.minor.y = element_line(size = 0.1, color = "#00000055")
      )
    fig_height <- 2 + length(unique(a2$obs$TRAV)) * 0.35
    my_ggsave(
      "composition-TRAV",
      #out_dir = glue("figures/{analysis_name}"),
      out_dir = file.path(out_dir, "tcr"),
      type = "pdf",
      plot = p +
        labs(
          title = "Per-donor percent of cells with each TRAV (Wilcoxon P-value)",
          y = NULL
        ),
      scale = 1, width = 8, height = fig_height, units = "in", dpi = 300
    )
  }
}

if (
  n_donors > 1 &&
  all(c("case", "donor", "leiden", "TRBV") %in% colnames(a2$obs)) &&
  sum(is.na(a2$obs$TRBV)) < 0.8 * nrow(a2$obs)
) {
  d <- a2$obs %>%
    dplyr::filter(!is.na(TRBV)) %>%
    dplyr::select(case, donor, TRBV) %>%
    # dplyr::mutate(case = class2) %>%
    dplyr::group_by(case, donor, TRBV) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(freq = 100 * n / sum(n))
  d_t <- rbindlist(lapply(unique(d$TRBV), function(this_TRBV) {
    x <- d %>% dplyr::filter(TRBV == this_TRBV)
    if (length(table(x$case)) != 2) {
      return(
        tibble(
          statistic = NA, p.value = 1, method = NA,
          alternative = NA, TRBV = this_TRBV
        )
      )
    }
    # x <- broom::tidy(t.test(n ~ case, x))
    # x <- broom::tidy(t.test(freq ~ case, x))
    x <- broom::tidy(wilcox.test(log(freq) ~ case, x))
    x$TRBV <- this_TRBV
    return(x)
  }))
  d_label <- d %>%
    dplyr::group_by(case) %>%
    summarise(n = sum(n)) %>%
    dplyr::mutate(label = glue("{case} (n = {comma(n)})"))
  d_t <- inner_join(
    d_t,
    d %>% dplyr::group_by(TRBV) %>% dplyr::mutate(y = mean(freq) + 4) %>% select(TRBV, y) %>% unique,
    by = "TRBV"
  )
  scientific_10 <- function(x) {
    ifelse(x < 0.01,
      gsub("e", "%*%10^", scales::scientific_format(digits = 1)(x)),
      signif(x, 1)
    )
  }
  d_t$TRBV <- factor(d_t$TRBV, d_t$TRBV[order(-d_t$p.value)])
  d$TRBV <- factor(d$TRBV, levels(d_t$TRBV))
  p <- ggplot() +
    geom_segment(
      data = data.frame(
        y = c(10^(-1:2), 0.5 * 10^(-1:2))
      ),
      mapping = aes(
        x = y, xend = y, y = -Inf, yend = Inf
      ),
      size = 0.2, alpha = 0.3
    ) +
    scale_color_manual(
      values = c(NA, "black"), guide = "none"
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by(TRBV, case) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        y = factor(TRBV), xmin = ymin, x = y, xmax = ymax,
        group = case, fill = case
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    annotate(
      geom = "rect",
      xmin = 0, xmax = Inf,
      ymin = as.integer(sort(unique(d$TRBV))) - 0.5,
      ymax = as.integer(sort(unique(d$TRBV))) + 0.5,
      fill = rep(
        c(NA, "black"),
        length.out = length(unique(d$TRBV))
      ),
      alpha = 0.1
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by(TRBV, case) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        y = factor(TRBV), xmin = ymin, x = y, xmax = ymax,
        group = case, fill = case
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    geom_point(
      data = d,
      mapping = aes(y = TRBV, x = freq, fill = case),
      position = position_quasirandom(groupOnX = FALSE, dodge.width = 0.9),
      alpha = 0.4, shape = 21
    ) +
    geom_text(
      data = d_t,
      mapping = aes(
        y = TRBV,
        x = 85,
        # y = y + 2,
        label = ifelse(
          # p.value < 0.05 / 12,
          p.value < 0.5 / length(unique(a2$obs$TRBV)),
          scientific_10(p.value),
          ""
        )
      ),
      size = 5, parse = TRUE
    ) +
    scale_fill_manual(
      name = NULL, values = pals::okabe(3)[c(2,3)],
      labels = d_label$label
    ) +
    # scale_y_continuous(name = "Percent") +
    scale_x_log10(
      name = "Percent",
      #labels = scales::label_number()
      labels = function(x) signif(x, 3)
    ) +
    annotation_logticks(side = "b") +
    labs(y = NULL) +
    theme(
      # legend.position = c(1, 1),
      # legend.justification = c(1, 1),
      legend.position = "bottom",
      legend.background = element_blank()
      # panel.grid.major.y = element_line(size = 0.1, color = "#00000055"),
      # panel.grid.minor.y = element_line(size = 0.1, color = "#00000055")
    )
  fig_height <- 2 + length(unique(a2$obs$TRBV)) * 0.35
  my_ggsave(
    "composition-TRBV",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = file.path(out_dir, "tcr"),
    type = "pdf",
    plot = p +
      labs(
        title = "Per-donor percent of cells with each TRBV (Wilcoxon P-value)",
        y = NULL
      ),
    scale = 1, width = 8, height = fig_height, units = "in", dpi = 300
  )
}


# Number of markers vs number of cells
########################################################################

d <- left_join(
  a2$de,
  a2$obs %>%
    dplyr::count(leiden, name = "n_cells") %>%
    dplyr::mutate(group = as.character(leiden)),
  by = "group"
)

x <- d %>%
  dplyr::group_by(n_cells, group) %>%
  dplyr::summarize(n_markers = sum(auc >= 0.7))

p <- ggplot(x) +
  aes(n_markers, n_cells, fill = group) +
  geom_point(size = 10, shape = 21) +
  geom_shadowtext(aes(label = group), size = 5) +
  scale_fill_manual(values = cluster_colors, guide = "none") + 
  scale_y_log10(labels = comma) +
  scale_x_log10(labels = comma) +
  annotation_logticks(side = "bl") +
  labs(
    x = "Markers with AUC >= 0.7", y = "Cells",
    title = glue("Markers and cells per cluster (n = {length(unique(d$leiden))})")
  )
my_ggsave(
  "markers-per-cluster",
  #out_dir = glue("figures/{analysis_name}"),
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1, width = 6, height = 4, units = "in", dpi = 300
)


# Heatmap of all features for each cluster
########################################################################
these_genes <- unique(rownames(a2$log2cpm))
if (length(these_genes) < 300) {
  fig_height <- length(these_genes) * 0.2
  fig_width <- length(unique(a2$de$group)) * 0.4 + 3
  #
  # Subset the table
  x <- a2$de %>% dplyr::filter(feature %in% these_genes)
  x <- as.data.frame(dcast.data.table(
    data = as.data.table(x), formula = symbol ~ group, value.var = "auc"
  ))
  rownames(x) <- x$symbol
  x$symbol <- NULL
  # Order the heatmap
  set.seed(1)
  xo <- seriation::seriate(as.matrix(x), method = "BEA_TSP")
  xd <- a2$de %>% dplyr::filter(feature %in% these_genes)
  xd$symbol <- factor(xd$symbol, rownames(x)[xo[[1]]])
  xd$group <- factor(xd$group, colnames(x)[xo[[2]]])
  # Plot
  p2 <- ggplot(xd) +
  geom_tile(
    aes(y = symbol, x = group, fill = auc)
  ) +
  scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  # scale_fill_viridis_c(
  #   name = "AUC", guide = guide_colorbar(barheight = 40), breaks = pretty_breaks(10),
  #   limits = c(0, 1)
  # ) +
  scale_fill_scico(
    palette = "vik", limits = c(0, 1), direction = 1,
    name = "AUC", guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(9)
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "italic"),
  #   axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )
  p1 <- ggplot(xd) +
  geom_tile(
    aes(y = 1, x = group, fill = group)
  ) +
  scale_x_discrete(position = "t", name = "Cluster", expand = c(0, 0)) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_fill_manual(
    values = cluster_colors, guide = "none"
  )
  p <- (
    p1 + theme(plot.margin = margin(b = 0))
  ) / (
    p2 + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0)
    )
  ) + plot_layout(heights = c(1, length(these_genes)))
  p
  my_ggsave(
    slug = "heatmap-all-markers",
    #out_dir = glue("figures/{analysis_name}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )
}

# stop here
return()

# Plot PC density before and after Harmony
########################################################################

plot_pc_density <- function(pca_dat, i = 4, this_color = "chemistry") {
  my_colors <- pals::okabe()
  if (length(unique(pca_dat[[this_color]])) > length(my_colors)) {
    my_colors <- rep(pals::alphabet(), 10)
  }
  p <- ggplot(pca_dat) +
  # scale_color_manual(name = this_color, values = my_colors[[this_color]]) +
  scale_color_manual(
    name = this_color,
    values = S4Vectors::unname(my_colors)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
  p1 <- p +
  geom_line(
    mapping = aes_string(x = sprintf("PC%s", i), color = this_color),
    stat = "density", adjust = 2, size = 1
  ) +
  theme(legend.position = "none") +
  labs(
    x = glue::glue("PC{i}"),
    # x = sprintf(
    #   "PC%s (%s%%)", i,
    #   signif(100 * uns$pca$variance_ratio[[i]], 2)
    # ),
    y = "Density"
  )
  p1
}

if ("pca" %in% names(a2)) {
  my_pc_groups <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))
  pca_dat <- as.data.frame(a2$pca$u)
  colnames(pca_dat) <- sprintf("PC%s", 1:ncol(pca_dat))
  pca_dat_channel <- as.data.frame(a2$pca_h)
  for (this_color in c("case", "donor")) {
    if (this_color %in% colnames(a2$obs)) {
      for (my_pcs in my_pc_groups) {
        i <- my_pcs[1]
        j <- my_pcs[2]
        pca_dat[[this_color]] <- a2$obs[[this_color]]
        pca_dat_channel[[this_color]] <- a2$obs[[this_color]]
        p1 <- plot_pc_density(pca_dat, i, this_color)
        p2 <- plot_pc_density(pca_dat, j, this_color) +
          theme(legend.position = "right")
        my_ggsave(
          slug = glue("pca-{this_color}-{i}-{j}"),
          #out_dir = glue("figures/{analysis_name}/pca_density"),
          out_dir = glue("{out_dir}/pca_density"),
          type = "pdf",
          plot = p1 + p2 + plot_layout(ncol = 2) + plot_annotation(),
          scale = 1, width = 7, height = 2.5, units = "in", dpi = 300
        )
        p1 <- plot_pc_density(pca_dat_channel, i, this_color)
        p2 <- plot_pc_density(pca_dat_channel, j, this_color) +
          theme(legend.position = "right")
        my_ggsave(
          glue("pca-harmonized_channel-{this_color}-{i}-{j}"),
          #out_dir = glue("figures/{analysis_name}/pca_density"),
          out_dir = glue("{out_dir}/pca_density"),
          type = "pdf",
          plot = p1 + p2 + plot_layout(ncol = 2) + plot_annotation(),
          scale = 1, width = 7, height = 2.5, units = "in", dpi = 300
        )
      }
    }
  }
  rm(list = c("pca_dat", "pca_dat_channel"))
}


# UMAP for each PC, including genes with large loading values
########################################################################

extreme_n <- function(xs, nlow = 5, nhigh = nlow) {
  xr <- rank(xs)
  c(which(xr <= nlow), which(xr >= max(xr) - nhigh + 1))
}

if (FALSE) {
for (my_pc in seq(ncol(a2$pca$v))) {
  pdf_file <- glue("figures/{analysis_name}/pca_umap/umap-pc-{my_pc}.pdf")
  if (file.exists(pdf_file)) {
    print_status(glue("File exists '{pdf_file}'"))
    next
  }
  # my_cor <- proxyC::simil(
  #   x = Matrix(a2$pca$u, sparse = TRUE),
  #   y = t(a2$log2cpm[a2$pca$genes,]),
  #   margin = 2,
  #   rank = 14,
  #   method = "correlation"
  # )
  # my_genes <- a2$pca$genes[extreme_n(my_cor[my_pc,], nlow = 7)]
  # my_pc <- 1
  # Select genes with most extreme loading values
  # my_genes <- a2$pca$genes[extreme_n(a2$pca$v[,my_pc], nlow = 7)]
  my_genes <- rownames(a2$log2cpm)[
    a2$ix_genes[extreme_n(a2$pca$v[,my_pc], nlow = 7)]
  ]
  # my_genes <- rowname_key[my_genes]
  p1 <- plot_hexgene(
    x = a2$obs$UMAP1,
    y = a2$obs$UMAP2,
    z = a2$pca$u[,my_pc],
    bins = 51
  ) +
  scale_fill_gradientn(
    colors = scico::scico(n = 100, direction = -1, palette = "roma"),
    labels = NULL
  ) +
  labs(title = glue("PC{my_pc}")) +
  guides(fill = guide_colorbar(title = NULL, barwidth = 5, ticks = FALSE))
  plots <- lapply(my_genes, function(this_gene) {
    plot_hexgene(
      x = a2$obs$UMAP1,
      y = a2$obs$UMAP2,
      z = as.numeric(a2$log2cpm[this_gene,]),
      bins = 51
    ) + labs(title = rowname_key[this_gene]) +
    theme(legend.position = "none")
  })
  ncol <- which.min(abs(sapply(1:16, function(i) {
    i - (16 / 9) * (length(plots) / i)
  })))
  p <- wrap_plots(c(list(p1), plots), ncol = ncol)
  my_ggsave(
    slug = glue("umap-pc-{my_pc}"),
    #out_dir = glue("figures/{analysis_name}/pca_umap"),
    out_dir = glue("{out_dir}/pca_umap"),
    plot = p,
    type = "pdf",
    scale = 1, width = 16, height = 9, units = "in", dpi = 300
  )
}
}



# All proteins ordered by BEA_TSP
########################################################################

  #entropy <- function(target) {
  #  freq <- table(target)/length(target)
  #  # vectorize
  #  vec <- as.data.frame(freq)[,2]
  #  #drop 0 to avoid NaN resulting from log2
  #  vec <- vec[vec > 0]
  #  #compute entropy
  #  -sum(vec * log2(vec))
  #}

  # my_genes_auc <- (
  #   a2$de %>%
  #     dplyr::group_by(feature) %>%
  #     # dplyr::summarize(score = max(auc / sum(auc))) %>%
  #     dplyr::summarize(score = entropy(auc)) %>%
  #     dplyr::arrange(-score)
  # )$feature

  # # AUC compared to next best AUC
  # my_genes_auc <- (
  #   a2$de %>%
  #     dplyr::group_by(feature) %>%
  #     dplyr::summarize(
  #       score = max(auc) / sort(auc, decreasing = TRUE)[2]
  #     ) %>%
  #     dplyr::arrange(-score)
  # )$feature

  # BEA_TSP
  x <- seriate_dataframe(a2$de, "feature", "group", "auc")
  my_genes_auc <- levels(x$feature)

  ## F statistic
  #print_status("Sort features by F statistic across groups")
  #y <- with(a2$obs, model.matrix(~ 0 + factor(leiden):factor(donor)))
  #a2$counts <- counts[,colnames(a2$log2cpm)]
  #pb <- a2$counts %*% y
  #pb <- as(pb, "dgCMatrix")
  #pb <- do_log2cpm(pb, median(colSums(pb)))
  #library(limma)
  ##
  #pb_meta <- str_split_fixed(colnames(pb), ":", 2)
  #colnames(pb_meta) <- c("leiden", "donor")
  #pb_meta <- as_tibble(pb_meta)
  #pb_meta %<>%
  #  mutate(
  #    leiden = str_replace(leiden, "factor\\(leiden\\)", ""),
  #    donor = str_replace(donor, "factor\\(donor\\)", "")
  #  )
  #pb_meta <- left_join(
  #  pb_meta, samples, by = "donor"
  #)
  #pb_meta$case <- factor(pb_meta$case, c("healthy", "mild", "severe"))
  #stopifnot(nrow(pb_meta) == ncol(pb))
  #pb_meta
  ##
  #des1 <- with(
  #  pb_meta,
  #  model.matrix(~ 0 + leiden)
  #)
  #colnames(des1) <- str_replace(colnames(des1), "^leiden", "")
  #fit1 <- lmFit(object = as.matrix(pb[rowMeans(pb) > 0.5,]), design = des1)
  #fit1 <- eBayes(fit1)
  #fit1$genes <- rowname_key[rownames(fit1$coefficients)]
  #top1 <- topTable(fit1, number = 1e6)
  #top1$feature <- rownames(top1)
  #top1 %<>% rename(ProbeID = "Gene") %>% as_tibble
  #my_genes_auc <- top1$feature

  if (length(my_genes_auc) < 300) {

    # # x <- a2$log2cpm[1,]
    # x <- sort(apply(a2$log2cpm, 1, function(x) {
    #   sum(x[rank(x) > 0.8 * length(x)]) / sum(x)
    # }), decreasing = TRUE)
    # my_genes_auc <- names(x)

    x <- ceiling(seq(0, length(my_genes_auc), length.out = 20))
    for (k in 2:length(x)) {
      i <- x[k - 1] + 1
      j <- x[k]
      pdf_file <- glue("figures/{analysis_name}/all_{nrow(a2$log2cpm)}_features/features_{i}-{j}.pdf")
      # if (file.exists(pdf_file)) {
      #   print_status(glue("File exists '{pdf_file}'"))
      #   next
      # }
      my_genes <- my_genes_auc[c(i:j)]
      plots <- lapply(my_genes, function(this_gene) {
        plot_hexgene(
          x = a2$obs$UMAP1,
          y = a2$obs$UMAP2,
          z = as.numeric(a2$log2cpm[this_gene,]),
          bins = 91,
          palette = "davos",
          direction = -1
        ) + labs(
          title = glue::glue("{this_gene}, {rowname_key[this_gene]}")
        ) +
        theme(legend.position = "none")
      })
      ncol <- which.min(abs(sapply(1:16, function(i) {
        i - (16 / 9) * (length(plots) / i)
      })))
      # p <- wrap_plots(c(list(p1), plots), ncol = ncol)
      p <- wrap_plots(plots, ncol = ncol)
      my_ggsave(
        slug = glue("features_{i}-{j}"),
        # out_dir = glue("figures/{analysis_name}/all_{nrow(a2$log2cpm)}_features"),
        out_dir = glue("{out_dir}/all_{nrow(a2$log2cpm)}_features"),
        plot = p,
        type = "pdf",
        scale = 1.5, width = 16, height = 9, units = "in", dpi = 300
      )
    }

  }
  
  # Not worth pursuing...
  #
  # a2$obs$cluster <- factor(a2$obs$leiden)
  # pc_lm <- rbindlist(lapply(seq(ncol(a2$pca$u)), function(my_pc) {
  #   my_form <- glue::glue("PC{my_pc} ~ 0 + cluster")
  #   res <- broom::tidy(summary(lm(as.formula(my_form), a2$obs)))
  #   # res <- broom::tidy(anova(aov(as.formula(my_form), a2$obs)))[1,]
  #   res$pc <- my_pc
  #   res
  # }))
  # pc_lm <- seriate_dataframe(pc_lm, "pc", "term", "statistic")
  # ggplot(pc_lm) +
  #   aes(x = pc, y = term, fill = p.value < 0.005 / nrow(pc_lm)) +
  #   geom_tile()
  #   # scale_fill_scico(palette = "nuuk", direction = -1)

## Gene set enrichment on PC loadings
#########################################################################
#
## MSigDB from Broad Institute
#msig <- readRDS("data/gene_sets.rds")
#
#do_gsea <- function(
#  dat, col_group, col_value, col_name, pathways
#) {
#  res <- list()
#  for (this_group in unique(dat[[col_group]])) {
#    ix <- dat[[col_group]] == this_group
#    x <- structure(
#      .Data = dat[[col_value]][ix],
#      .Names = dat[[col_name]][ix]
#    )
#    # msig_genes <- msig[[this_set]] %>% unlist %>% unname %>% unique
#    # x <- x[intersect(msig_genes, names(x))]
#    d <- fgsea::fgseaMultilevel(
#      pathways = pathways,
#      stats = x,
#      minSize = 15,
#      maxSize = 500
#    )
#    d$group <- this_group
#    res[[this_group]] <- d %>%
#      dplyr::arrange(padj) %>%
#      dplyr::filter(!is.na(pval))
#  }
#  do.call(rbind, res)
#}
#
#msig_names <- c(
#  # "h_all", "c1_all", "c2_cp_kegg", "c3_all", "c4_all", "c5_mf", "c6_all", "c7_all"
#  "h_all", "c2_cp_kegg", "c5_mf", "c5_bp", "c5_cc", "c7_all"
#)
##
#for (this_set in msig_names) {
#  print_status(
#    glue("Gene set enrichment for MSigDB gene set {this_set}")
#  )
#  x <- as.data.frame(a2$pca$v)
#  colnames(x) <- sprintf("%s", seq(ncol(x)))
#  x <- cbind(feature = a2$pca$genes, x)
#  x <- reshape2::melt(x, id.vars = "feature")
#  colnames(x) <- c("feature", "group", "value")
#  x$feature <- as.character(x$feature)
#  res <- do_gsea(
#    x, "group", "value", "feature", msig[[this_set]]
#  )
#  #
#  my_terms <- unique((
#    res %>%
#    filter(padj < 0.05) %>%
#    group_by(group) %>%
#    top_n(wt = abs(NES), n = 10)
#  )$pathway)
#  #
#  x <- res %>% filter(pathway %in% my_terms)
#  # x$sscore2 <- x$NES + abs(min(x$NES))
#  x$sscore2 <- -log10(x$pval) * sign(x$NES)
#  x$term <- factor(as.character(x$pathway))
#  x$cluster <- factor(as.character(x$group))
#  #
#  offset <- abs(min(x$sscore2))
#  x$sscore2 <- x$sscore2 + offset
#  set.seed(2)
#  xo <- order_dataframe(x, "cluster ~ term", "sscore2", method = "BEA_TSP")
#  x$sscore2 <- x$sscore2 - offset
#  #
#  x$cluster <- fct_relevel(x$cluster, levels(x$cluster)[xo[[1]]])
#  x$term <- fct_relevel(x$term, levels(x$term)[xo[[2]]])
#  #
#  n_clusters <- length(unique(x$cluster))
#  n_chars <- max(sapply(my_terms, nchar))
#  fig_width <- 1.5 + 0.25 * n_clusters + 0.18 * n_chars
#  fig_height <- 1 + 0.25 * length(my_terms)
#  p <- ggplot(x) +
#  aes(x = cluster, fill = sscore2, y = term) +
#  # aes(x = cluster, fill = -log10(padj), y = term) +
#  # aes(x = cluster, fill = NES, y = term) +
#  geom_tile() +
#  scale_y_discrete(
#    expand = c(0, 0),
#    labels = function(x) {
#      str_replace_all(str_replace(x, "HALLMARK_", ""), "_", " ")
#    }
#  ) +
#  scale_x_discrete(
#    expand = c(0, 0)
#  ) +
#  # scale_fill_viridis_c(
#  #   # name = bquote("log"[10]~"P"),
#  #   guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(10),
#  # ) +
#  # scale_fill_scico(
#  #   palette = "vik", direction = 1,
#  #   guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(5)
#  # ) +
#  scale_fill_scico(
#    name = "-log10p",
#    palette = "vik", direction = 1,
#    limits = c(-1, 1) * max(abs(x$sscore2)),
#    guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(5)
#  ) +
#  labs(
#    x = "PC", y = NULL,
#    title = glue("MSigDB {this_set} enrichment results")
#  ) +
#  theme(axis.ticks.x = element_blank())
#  my_ggsave(
#    slug = glue("fgsea-pca-{this_set}"),
#    # out_dir = glue("figures/{analysis_name}/gsea"),
#    out_dir = glue("{out_dir}/gsea"),
#    type = "pdf",
#    plot = p,
#    scale = 1, width = fig_width, height = fig_height,
#    units = "in", dpi = 300, limitsize = FALSE
#  )
#}
#
#
## Gene set enrichment on DE genes
#########################################################################
#
## msig_names <- c(
##   # "h_all", "c1_all", "c2_cp_kegg", "c3_all", "c4_all", "c5_mf", "c6_all", "c7_all"
##   "h_all", "c2_cp_kegg", "c5_mf", "c7_all"
## )
##
#for (this_set in msig_names) {
#  print_status(
#    glue("Gene set enrichment for MSigDB gene set {this_set}")
#  )
#  res <- do_gsea(
#    a2$de, "group", "logFC", "feature", msig[[this_set]]
#  )
#  #
#  my_terms <- unique((
#    res %>%
#    filter(padj < 0.05) %>%
#    group_by(group) %>%
#    top_n(wt = abs(NES), n = 10)
#  )$pathway)
#  #
#  x <- res %>% filter(pathway %in% my_terms)
#  # x$sscore2 <- x$NES + abs(min(x$NES))
#  x$sscore2 <- -log10(x$pval) * sign(x$NES)
#  x$term <- factor(as.character(x$pathway))
#  x$cluster <- factor(as.character(x$group))
#  #
#  offset <- abs(min(x$sscore2))
#  x$sscore2 <- x$sscore2 + offset
#  set.seed(2)
#  xo <- order_dataframe(x, "cluster ~ term", "sscore2", method = "BEA_TSP")
#  x$sscore2 <- x$sscore2 - offset
#  #
#  x$cluster <- fct_relevel(x$cluster, levels(x$cluster)[xo[[1]]])
#  x$term <- fct_relevel(x$term, levels(x$term)[xo[[2]]])
#  #
#  n_clusters <- length(unique(x$cluster))
#  n_chars <- max(sapply(my_terms, nchar))
#  fig_width <- 1.5 + 0.25 * n_clusters + 0.18 * n_chars
#  fig_height <- 1 + 0.25 * length(my_terms)
#  p <- ggplot(x) +
#  aes(x = cluster, fill = sscore2, y = term) +
#  # aes(x = cluster, fill = -log10(padj), y = term) +
#  # aes(x = cluster, fill = NES, y = term) +
#  geom_tile() +
#  scale_y_discrete(
#    expand = c(0, 0),
#    labels = function(x) {
#      str_replace_all(str_replace(x, "HALLMARK_", ""), "_", " ")
#    }
#  ) +
#  scale_x_discrete(
#    expand = c(0, 0)
#  ) +
#  # scale_fill_viridis_c(
#  #   # name = bquote("log"[10]~"P"),
#  #   guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(10),
#  # ) +
#  # scale_fill_scico(
#  #   palette = "vik", direction = 1,
#  #   guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(5)
#  # ) +
#  # scale_fill_scico(
#  #   guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(5)
#  # ) +
#  scale_fill_scico(
#    name = "-log10p",
#    palette = "vik", direction = 1,
#    limits = c(-1, 1) * max(abs(x$sscore2)),
#    guide = guide_colorbar(barheight = 20), breaks = pretty_breaks(5)
#  ) +
#  labs(
#    x = "Cluster", y = NULL,
#    title = glue("MSigDB {this_set} enrichment results")
#  ) +
#  theme(axis.ticks.x = element_blank())
#  my_ggsave(
#    slug = glue("fgsea-logfc-{this_set}"),
#    # out_dir = glue("figures/{analysis_name}/gsea"),
#    out_dir = glue("{out_dir}/gsea"),
#    type = "pdf",
#    plot = p,
#    scale = 1, width = fig_width, height = fig_height,
#    units = "in", dpi = 300
#  )
#}

}

