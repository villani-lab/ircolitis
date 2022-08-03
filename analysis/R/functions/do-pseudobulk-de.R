
library(data.table)
library(limma)
library(pbapply)

do_pseudobulk_de <- function(a2, out_dir) {

  stopifnot("cluster" %in% colnames(a2$obs))
  stopifnot("donor" %in% colnames(a2$obs))

  y <- with(a2$obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/") # means
  pb <- as(a2$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
  #
  library(limma)
  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
  colnames(pb_meta) <- c("cluster", "donor")
  pb_meta <- as_tibble(pb_meta)
  pb_meta %<>%
    mutate(
      cluster = str_replace(cluster, "factor\\(cluster\\)", ""),
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
  # pb_meta$cluster <- fct_reorder(
  #   pb_meta$cluster, pb_meta$gene,
  #   .fun = mean
  # )
  # ggplot(pb_meta) +
  #   aes(y = cluster, x = gene) +
  #   geom_quasirandom(, groupOnX = FALSE)
  # What percent of cells have this gene?
  # a2$obs$gene <- a2$log2cpm[this_ens,]
  # a2$obs %>%
  #   group_by(cluster) %>%
  #   summarize(pct = sum(gene > 0) / length(gene)) %>%
  # ggplot() +
  # aes(y = as.character(cluster), x = pct) +
  # geom_colh()

  x <- table(a2$obs$cluster)
  min_percent <- 100 * min(x) * 0.25 / sum(x)
  keep_ens <- a2$counts_stats$gene[a2$counts_stats$percent >= min_percent]

  # All versus All (AVA)
  # Test all pairs of clusters
  print_status("Finding marker genes with pseudobulk")
  pb_meta$x <- factor(pb_meta$cluster)
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
    pblapply(sort(unique(pb_meta$cluster)), function(this_cluster) {
      pb_meta$x <- pb_meta$cluster == this_cluster
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
  de_ova_file <- glue("{out_dir}/pseudobulk_de_ova.tsv.gz")
  print_status(glue("Writing {de_ova_file}"))
  data.table::fwrite(de_ova, de_ova_file, sep = "\t")

}

do_log1p_cpm <- function(counts, total = NULL) {
  stopifnot(is(counts, "dgCMatrix"))
  if (is.null(total)) {
    total <- median(Matrix::colSums(counts))
  }
  counts@x <- counts@x / rep.int(Matrix::colSums(counts), diff(counts@p))
  counts@x <- total * counts@x
  counts@x <- log1p(counts@x)
  return(counts)
}

make_pseudobulk <- function(counts, factor1, factor2) {
  stopifnot(length(factor1) == ncol(counts))
  stopifnot(length(factor2) == ncol(counts))
  y <- model.matrix(~ 0 + factor(factor1):factor(factor2))
  y <- as(y, "dgCMatrix")
  pb <- as(counts %*% y, "dgCMatrix")
  pb_sums <- Matrix::colSums(pb)
  pb <- pb[,pb_sums > 0]
  # pb <- do_log1p_cpm(pb, median(pb_sums[pb_sums > 0]))
  pb <- do_log2cpm(pb, median(pb_sums[pb_sums > 0]))
  #
  pb_meta <- as.data.frame(str_split_fixed(colnames(pb), ":", 2))
  colnames(pb_meta) <- c("factor1", "factor2")
  pb_meta$factor1 <- str_replace(pb_meta$factor1, "factor\\(factor1\\)", "")
  pb_meta$factor2 <- str_replace(pb_meta$factor2, "factor\\(factor2\\)", "")
  colnames(pb) <- str_replace(colnames(pb), "factor\\(factor1\\)", "")
  colnames(pb) <- str_replace(colnames(pb), "factor\\(factor2\\)", "")
  list(logcpm = pb, obs = pb_meta)
}

do_de <- function(counts, cluster, donor, min_mean = 0.5) {
  pb <- make_pseudobulk(counts, cluster, donor)
  # Fit a linear model with one coefficient for each cluster.
  pb$obs$xx <- factor(pb$obs$factor1)
  # print(levels(pb$obs$xx))
  # print(mean(rowMeans(pb$logcpm)))
  des1 <- with(pb$obs, model.matrix(~ 0 + xx))
  fit1 <- lmFit(object = pb$logcpm[rowMeans(pb$logcpm) > min_mean,], design = des1)
  fit1 <- eBayes(fit1)
  # All versus All (AVA)
  # Test all pairs of clusters
  cluster_pairs <- t(combn(levels(pb$obs$xx), 2))
  cont <- makeContrasts(contrasts = lapply(seq(nrow(cluster_pairs)), function(i) {
    glue("xx{cluster_pairs[i,1]} - xx{cluster_pairs[i,2]}")
  }), levels = des1)
  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
  fit2 <- contrasts.fit(fit1, cont)
  fit2 <- eBayes(fit2)
  de_ava <- rbindlist(pblapply(colnames(cont), function(this_contrast) {
    res <- topTable(fit2, coef = this_contrast, number = nrow(fit1$coefficients), confint = TRUE)
    this_contrast <- str_replace_all(this_contrast, "xx", "")
    this_contrast <- str_replace(this_contrast, "vs", " vs ")
    res$contrast <- this_contrast
    res$feature <- rownames(res)
    return(res)
  }))
  # One versus All (OVA)
  de_ova <- rbindlist(
    pblapply(sort(unique(pb$obs$factor1)), function(this_cluster) {
      pb$obs$x <- pb$obs$factor1 == this_cluster
      des1 <- with(pb$obs, model.matrix(~ x))
      fit1 <- lmFit(object = pb$logcpm[rowMeans(pb$logcpm) > min_mean,], design = des1)
      fit1 <- eBayes(fit1)
      res <- topTable(fit1, coef = 2, number = 1e6, confint = TRUE)
      res$contrast <- sprintf("%s vs all", this_cluster)
      res$feature <- rownames(res)
      return(res)
    })
  )
  pb$obs$xx <- NULL
  list(pb = pb, ava = de_ava, ova = de_ova)
}

#' Write differential expression results to an Excel file.
#' @param d Dataframe with differential expression results from presto::wilcoxauc()
#' @param fname Filename of the output Excel file.
write_de <- function(d, fname, n = 500) {
  wb <- openxlsx::createWorkbook()
  #fname <- "analysis/nuclei/louvain_de.xlsx"
  unlink(fname)
  for (this_group in sort(unique(d$group))) {
    openxlsx::addWorksheet(wb, as.character(this_group))
    x <- d %>%
      dplyr::filter(group == this_group) %>%
      dplyr::top_n(n = n, wt = abs(0.5 - auc)) %>%
      dplyr::arrange(-auc) %>%
      # top_n(n = 300, wt = -logFC * (pct_in - pct_out)) %>%
      as.data.frame
    openxlsx::writeDataTable(
      wb,
      as.character(this_group),
      x = x, rowNames = FALSE, tableStyle = "TableStyleLight1"
    )
  }
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}

write_de_xlsx <- function(d, fname, col = NULL) {
  wb <- openxlsx::createWorkbook()
  #fname <- "analysis/nuclei/louvain_de.xlsx"
  unlink(fname)
  if (!is.null(col)) {
    for (this_sheet in naturalsort::naturalsort(unique(d[[col]]))) {
      openxlsx::addWorksheet(wb, as.character(this_sheet))
      x <- d[d[[col]] == this_sheet,,drop=FALSE] %>%
        # head(n) %>%
        # dplyr::top_n(n = n, wt = -log10(P.Value)) %>%
        # arrange(P.Value) %>%
        as.data.frame
      openxlsx::writeDataTable(
        wb,
        as.character(this_sheet),
        x = x, rowNames = FALSE, tableStyle = "TableStyleLight1"
      )
    }
  } else {
    openxlsx::addWorksheet(wb, "Sheet1")
    openxlsx::writeDataTable(
      wb,
      "Sheet1",
      x = as.data.frame(d),
      rowNames = FALSE, tableStyle = "TableStyleLight1"
    )
  }
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}

#make_de_files <- function(counts, clusters, samples) {
#  logcpm <- do_log1p_cpm(counts)
#  de_auc <- presto::wilcoxauc(logcpm, clusters)
#  #
#  de <- do_de(counts, clusters, samples)
#  #
#  ova <- de$ova %>%
#    mutate(contrast = str_split_fixed(contrast, " vs ", 2)[,1]) %>%
#    select(-B, -t) %>%
#    rename(cluster = contrast, pval = P.Value, fdr = adj.P.Val, mean = AveExpr, log2fc = logFC) %>%
#    mutate(feature = glue("{feature}|{ensembl_to_symbol[feature]}"))
#  #
#  ova_auc <- a1$de %>% select(group, feature, auc, pct_in, pct_out) %>%
#    mutate(pct_ratio = pct_in / pct_out)
#  #
#  ova <- left_join(ova, ova_auc, by = c("cluster" = "group", "feature"))
#  #
#  ova <- ova %>%
#    mutate(gene = str_split_fixed(feature, "\\|", 2)[,2]) %>%
#    relocate(cluster, gene, auc) %>%
#    select(-feature) %>%
#    group_by(cluster) %>%
#    top_n(n = 3000, wt = auc) %>%
#    ungroup() %>%
#    arrange(cluster, -auc) %>%
#    mutate_if(is.numeric, signif, 4) %>%
#    mutate_if(is.numeric, as.character)

  #ava <- de$ava %>%
  #  mutate(gene = ensembl_to_symbol[feature]) %>%
  #  select(contrast, gene, logFC, AveExpr, P.Value, adj.P.Val, feature) %>%
  #  mutate(feature = glue("{feature}|{gene}")) %>%
  #  rename(cluster = contrast, pval = P.Value, fdr = adj.P.Val, mean = AveExpr, log2fc = logFC)
  #ava$c1 <- str_split_fixed(ava$cluster, " vs ", 2)[,1]
  #ava$c2 <- str_split_fixed(ava$cluster, " vs ", 2)[,2]
  #ava <- left_join(
  #  x = ava, 
  #  y = a1$de[,c("feature", "group", "pct_in")],
  #  by = c("feature" = "feature", "c1" = "group")
  #)
  #ava <- left_join(
  #  x = ava, 
  #  y = a1$de[,c("feature", "group", "pct_in")],
  #  by = c("feature" = "feature", "c2" = "group")
  #)
  #ava <- ava %>%
  #  rename(pct_in = pct_in.x, pct_out = pct_in.y) %>%
  #  mutate(pct_ratio = pct_in / pct_out) %>%
  #  mutate(cluster = str_replace_all(cluster, ' ', '-')) %>%
  #  # select(cluster, gene, pval, pct_in, pct_out) %>%
  #  select(-c1, -c2, -feature) %>%
  #  group_by(cluster) %>%
  #  top_n(n = 3000, wt = -log10(pval)) %>%
  #  ungroup() %>%
  #  arrange(cluster, pval) %>%
  #  mutate_if(is.numeric, signif, 2) %>%
  #  mutate_if(is.numeric, as.character) %>%
  #  as_tibble
  ##
  #ova_markers_file <- file.path(cb_dir, glue("{leiden_col}-ova-markers.tsv"))
  #print_status(glue("Writing {length(unique(ova$gene))} genes to {ova_markers_file}"))
  #write_tsv(ova, ova_markers_file)
  ##
  #ava_markers_file <- file.path(cb_dir, glue("{leiden_col}-ava-markers.tsv"))
  #print_status(glue("Writing {length(unique(ava$gene))} genes to {ava_markers_file}"))
  #write_tsv(ava, ava_markers_file)
#}

