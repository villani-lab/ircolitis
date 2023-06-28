#!/usr/bin/env Rscript

# libraries {{{
pacman::p_load(
  dendsort,
  data.table,
  dplyr,
  forcats,
  ggplot2,
  ggrastr,
  ggrepel,
  ggtext,
  glue,
  janitor,
  limma,
  magrittr,
  naturalsort,
  parameters,
  pbapply,
  RColorBrewer,
  qs,
  readxl,
  scales,
  stringr,
  tibble,
  tidyr
)
source("R/functions/theme-kamil.R")
theme_set(theme_kamil)

palettes <- list(
  "Giles2022" = c(
    "CD4_BulkNaive"    = "#5e215e",
    "CD4_BulkNonNaive" = "#8c658e",
    "CD4_Tfh"          = "#354c96",
    "CD4_Treg"         = "#c0adce",
    "CD8_BulkNaive"    = "#576cae",
    "CD8_BulkNonNaive" = "#7c95cd",
    "CD8_CM"           = "#e8ea52",
    "CD8_EM1"          = "#d7b03e",
    "CD8_EM2"          = "#c78130",
    "CD8_EMRA"         = "#a31518",
    "CD8_Naive"        = "#629c55",
    "CD8_PD1+CD39+"    = "#c6548b",
    "CD8_SCM-R3-"      = "#bfd399",
    "CD8_SCM-R3+"      = "#8db972"
  )
)

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
root_dir <- "/projects/external_data/Giles2022"

# Path to output files (tsv, pdf)
out_dir <- "results/Giles2022"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "published-volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# }}}

# read expression data csv file {{{

d <- fread(glue("{root_dir}/GSE179609_HD_RNAseq_counts_normalizedLog2BatchCorrected.csv.gz"))
mat <- as.matrix(d[,3:ncol(d)])
rownames(mat) <- with(d, glue("{id}-{geneSymbol}"))
mat[1:5,1:5]

meta <- read_excel(glue("{root_dir}/metadata.xlsx"))
meta <- clean_names(meta)
meta$cell_subset <- factor(meta$cell_subset, c("CD4_BulkNaive",
    "CD4_BulkNonNaive", "CD4_Treg", "CD4_Tfh", "CD8_BulkNaive",
    "CD8_BulkNonNaive", "CD8_Naive", "CD8_SCM-R3+", "CD8_SCM-R3-", "CD8_CM",
    "CD8_EM1", "CD8_EM2", "CD8_EMRA", "CD8_PD1+CD39+"))

colnames(meta)

meta %>% count(cell_subset) # 14 subsets

celltypes <- c(
  "CD4_BulkNaive", "CD4_BulkNonNaive", "CD4_Tfh", "CD4_Treg",
  "CD8_BulkNaive", "CD8_BulkNonNaive", "CD8_CM", "CD8_EM1", "CD8_EM2",
  "CD8_EMRA", "CD8_Naive", "CD8_PD1+CD39+", "CD8_SCM-R3-", "CD8_SCM-R3+"
)

meta %>% count(tissue_type) # PBMC
meta %>% count(cell_type) # CD4, CD8

meta %>% count(sample_acquisiton_location) # NIA, UPenn

setdiff(meta$rna_id, colnames(mat)) # "NA"
setdiff(colnames(mat), meta$rna_id) # character(0)

meta <- meta[match(colnames(mat), meta$rna_id),]
stopifnot(all(meta$rna_id ==  colnames(mat)))

meta$x <- factor(make_clean_names(meta$cell_subset, allow_dupes = TRUE, replace = c("\\+"="p","-"="n")))
# meta %>% count(x, cell_subset)

# }}}

# volcanos for published DE results in supplemental table 3 {{{

de_file <- "/projects/external_data/Giles2022/HD-CD8T-DGE.xlsx"
my_sheets <- excel_sheets(de_file)
de <- rbindlist(lapply(my_sheets, function(sheet) {
  read_excel(de_file, sheet = sheet)
}))

head(de)

de$gene <- with(de, glue("{id}-{geneSymbol}"))

de$c1 <- with(de, str_split_fixed(subset_comparison, "_vs_", 2)[,1])
de$c2 <- with(de, str_split_fixed(subset_comparison, "_vs_", 2)[,2])

de %>% count(subset_comparison)

de %>% count(c2)

my_comp <- "naive2_vs_centralMemory"

for (my_comp in unique(de$subset_comparison)) {
  message(my_comp)
  my_de <- de %>% filter(subset_comparison == my_comp)
  my_labels <- my_de %>% arrange(padj) %>% head(20)
  p <- ggplot(my_de) +
    aes(x = log2FoldChange, y = -log10(padj)) +
    geom_vline(xintercept = 0, alpha = 0.3, linewidth = 0.3) +
    rasterise(geom_point(size = 0.1), dpi = 300) +
    rasterise(
      geom_point(data = my_labels, size = 0.3, color = "red"),
      dpi = 300
    ) +
    geom_text_repel(
      data = my_labels,
      mapping = aes(label = geneSymbol),
      size = 2, fontface = "italic",
      segment.size = 0.3
    ) +
    labs(title = my_comp)
  p_file <- glue("{out_dir}/published-volcano/volcano-{my_comp}.pdf")
  message(p_file)
  ggsave(p_file, p, width = 5, height = 5)
}


de %>% group_by(c1) %>%
  filter(log2FoldChange > 0) %>%
  arrange(c1, -log2FoldChange)

# }}}

# Giles2022 limma OVA and AVA {{{

# All versus All (AVA)
# Test all pairs of clusters
des1 <- with(meta, model.matrix(~ 0 + x))
fit1 <- lmFit(object = mat, design = des1)
fit1 <- eBayes(fit1)
cluster_pairs <- t(combn(levels(meta$x), 2))
cont <- makeContrasts(contrasts = lapply(seq(nrow(cluster_pairs)), function(i) {
  glue("x{cluster_pairs[i,1]} - x{cluster_pairs[i,2]}")
}), levels = des1)
colnames(cont) <- str_replace(colnames(cont), " - ", "_vs_")
fit2 <- contrasts.fit(fit1, cont)
fit2 <- eBayes(fit2)
de_ava <- rbindlist(lapply(colnames(cont), function(this_coef) {
  x <- topTable(fit2, coef = this_coef, number = nrow(fit1$coefficients))
  this_coef <- str_replace_all(this_coef, "x", "")
  this_coef <- str_replace(this_coef, "_vs_", " vs ")
  x$coef <- this_coef
  x$ensembl_id <- rownames(x)
  x
}))
de_ava_file <- glue("{out_dir}/tables/limma-ava.tsv.gz")
message(glue("Writing {de_ava_file}"))
data.table::fwrite(de_ava, de_ava_file, sep = "\t")

# One versus All (OVA)
de_ova <- rbindlist(
  pblapply(sort(unique(meta$cell_subset)), function(this_cluster) {
    meta$is_cluster <- meta$cell_subset == this_cluster
    des1 <- with(meta, model.matrix(~ is_cluster))
    fit1 <- lmFit(object = mat, design = des1)
    fit1 <- eBayes(fit1)
    res <- topTable(fit1, coef = 2, number = nrow(mat))
    res$coef <- sprintf("%s vs all", this_cluster)
    res$gene <- rownames(res)
    return(res)
  })
)
de_ova_file <- glue("{out_dir}/tables/limma-ova.tsv.gz")
message(glue("Writing {de_ova_file}"))
data.table::fwrite(de_ova, de_ova_file, sep = "\t")

de_ova$ensembl_id <- str_split_fixed(de_ova$gene, "-", 2)[,1]

# }}}

# Giles2022 limma OVA and AVA (matched CD8 and CD4) {{{

# Test all pairs of clusters
meta$x <- factor(make_clean_names(meta$cell_subset, allow_dupes = TRUE, replace = c("\\+"="p","-"="n")))

# CD8
ix <- which(with(meta, str_detect(x, "cd8_[^b]")))
my_meta <- meta[ix,]
my_mat <- mat[,ix]
# One versus All (OVA)
de_ova <- rbindlist(
  pblapply(sort(unique(my_meta$cell_subset)), function(this_cluster) {
    my_meta$is_cluster <- my_meta$cell_subset == this_cluster
    des1 <- with(my_meta, model.matrix(~ is_cluster))
    fit1 <- lmFit(object = my_mat, design = des1)
    fit1 <- eBayes(fit1)
    res <- topTable(fit1, coef = 2, number = nrow(my_mat))
    res$coef <- sprintf("%s vs all", this_cluster)
    res$gene <- rownames(res)
    return(res)
  })
)
de_ova$type <- "CD8"
de_ova$ensembl_id <- str_split_fixed(de_ova$gene, "-", 2)[,1]
de_ova_file <- glue("{out_dir}/tables/limma-ova-CD8.tsv.gz")
message(glue("Writing {de_ova_file}"))
data.table::fwrite(de_ova, de_ova_file, sep = "\t")

# CD4
ix <- which(with(meta, str_detect(x, "cd4_")))
meta$x[ix]
my_meta <- meta[ix,]
my_mat <- mat[,ix]
# One versus All (OVA)
de_ova <- rbindlist(
  pblapply(sort(unique(my_meta$cell_subset)), function(this_cluster) {
    my_meta$is_cluster <- my_meta$cell_subset == this_cluster
    des1 <- with(my_meta, model.matrix(~ is_cluster))
    fit1 <- lmFit(object = my_mat, design = des1)
    fit1 <- eBayes(fit1)
    res <- topTable(fit1, coef = 2, number = nrow(my_mat))
    res$coef <- sprintf("%s vs all", this_cluster)
    res$gene <- rownames(res)
    return(res)
  })
)
de_ova$type <- "CD4"
de_ova$ensembl_id <- str_split_fixed(de_ova$gene, "-", 2)[,1]
de_ova_file <- glue("{out_dir}/tables/limma-ova-CD4.tsv.gz")
message(glue("Writing {de_ova_file}"))
data.table::fwrite(de_ova, de_ova_file, sep = "\t")

de_ova <- rbindlist(lapply(Sys.glob(glue("{out_dir}/tables/limma-ova-CD*.tsv.gz")), function(file) {
  fread(file)
}))

# }}}

# compare Giles2022 OVA to ircolitis OVA within CD4 and CD8 {{{

ova <- fread("paper/de-ova.tsv.gz") %>%
  mutate(contrast2 = glue("{analysis} {contrast}")) %>%
  select(-c("analysis", "contrast"))
ova %>% count(contrast2)

ova %>%
  filter(contrast2 == "Blood CD8 T cells 11 vs all") %>%
  arrange(-auc) %>%
  filter(!str_detect(gene, "^RP")) %>%
  # filter(gene == "ITGB2") %>%
  head(10)

res_file <- glue("{out_dir}/tables/limma-cor-t.qs")
if (file.exists(res_file)) {
  res <- qread(res_file)
} else {
  res <- rbindlist(pblapply(unique(ova$contrast2), function(my_contrast1) {
    rbindlist(lapply(unique(de_ova$coef), function(my_contrast2) {
      # ircolitis
      # my_contrast1 <- "Blood B cells 1 vs all"
      my_ova1 <- ova %>% filter(contrast2 == my_contrast1)
      # Giles2022
      # my_contrast2 <- "CD4_BulkNaive vs all"
      my_ova2 <- clean_names(de_ova %>% filter(coef == my_contrast2))
      # head(my_ova1)
      # head(my_ova2)
      my_ova <- full_join(my_ova1, my_ova2, by = "ensembl_id")
      res <- suppressWarnings({
        # with(my_ova, cor.test(log_fc.x, log_fc.y, method = "spearman")) %>%
        with(my_ova, cor.test(t.x, t.y, method = "spearman", na.action = "na.omit")) %>%
        parameters %>%
        as.data.frame
      })
      res$contrast.x <- my_contrast1
      res$contrast.y <- my_contrast2
      return(res)
    }))
  }))
  qsave(res, res_file)
}
res %<>%
  mutate(
    analysis = str_extract(contrast.x, "^.+cells"),
    con.x = str_extract(contrast.x, " \\d+ "),
    con.y = str_extract(contrast.y, "[^ ]+")
  )
res_mat <- hclust_dataframe(res, "contrast.x", "contrast.y", "rho")
res$contrast.x <- factor(res$contrast.x, rownames(res_mat))
res$contrast.y <- factor(res$contrast.y, colnames(res_mat))

head(res)


# my_res <- res %>%
#   filter(str_detect(analysis, "(CD8|CD4)"))
# p <- ggplot(my_res) +
#   aes(y = contrast.x, x = contrast.y, fill = rho) +
#   facet_wrap(~ analysis, scales = "free") +
#   scale_fill_gradientn(
#     colors = rev(brewer.pal("RdBu", n = 11)),
#     rescaler = ~ rescale_mid(.x, mid = 0)
#   ) +
#   guides(fill = guide_colorbar(barwidth = 20)) +
#   geom_tile() +
#   geom_point(
#     data = ~ .x %>%
#       group_by(contrast.x) %>%
#       top_n(n = 1, wt = sign(rho) * -log10(p)),
#     size = 5, shape = 21, fill = "white", stroke = 0.3
#   ) +
#   scale_x_discrete(labels = ~ str_remove(.x, " vs all")) +
#   scale_y_discrete(labels = ~ str_remove(.x, " vs all")) +
#   theme(
#     legend.position = "top",
#     axis.text.x = element_text(angle = 30, hjust = 1)
#   )
# p_file <- glue("{out_dir}/limma-cor.pdf")
# my_width <- length(unique(my_res$contrast.y)) * 1.2 + 5
# my_height <- length(unique(my_res$contrast.x)) * 0.3 + 2
# message(p_file)
# ggsave(p_file, p, width = my_width, height = my_height, limitsize = FALSE)

filters <- list(
  "CD4" = c("cell_subset" = "CD4", "cell_cluster" = "Blood CD4"),
  "CD8" = c("cell_subset" = "CD8", "cell_cluster" = "Blood CD8")
)
for (my_type in names(filters)) {
  my_res2 <- res %>%
    filter(str_detect(analysis, "(CD8|CD4)")) %>%
    filter(
      str_detect(contrast.y, filters[[my_type]]["cell_subset"]),
      str_detect(contrast.x, filters[[my_type]]["cell_cluster"])
    )
  my_res2_mat <- hclust_dataframe(my_res2, "contrast.x", "contrast.y", "rho")
  my_res2$contrast.x <- factor(my_res2$contrast.x, rownames(my_res2_mat))
  my_res2$contrast.y <- factor(my_res2$contrast.y, colnames(my_res2_mat))
  p <- ggplot(my_res2) +
    aes(y = contrast.x, x = contrast.y, fill = rho) +
    scale_fill_gradientn(
      colors = rev(brewer.pal("RdBu", n = 11)),
      rescaler = ~ rescale_mid(.x, mid = 0)
    ) +
    guides(fill = guide_colorbar(barwidth = 15, title = "Spearman\n")) +
    geom_tile() +
    geom_point(
      data = ~ .x %>%
        group_by(contrast.x) %>%
        # top_n(n = 1, wt = sign(rho) * -log10(p)),
        top_n(n = 1, wt = rho),
      size = 5, shape = 21, fill = "white", stroke = 0.3
    ) +
    scale_x_discrete(
      expand = expansion(0, 0),
      labels = ~ str_replace_all(str_remove(.x, " vs all"), "_", " "),
      position = "top"
    ) +
    scale_y_discrete(
      expand = expansion(0, 0),
      labels = ~ str_remove(.x, " vs all"),
      position = "right"
    ) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 30, hjust = 0)
    ) +
    labs(x = NULL, y = NULL)
  p_file <- glue("{out_dir}/ova-cor/ova-cor-{my_type}.pdf")
  dir.create(dirname(p_file), recursive = TRUE, showWarnings = FALSE)
  my_width <- length(unique(my_res2$contrast.y)) * 0.4 + 4
  my_height <- length(unique(my_res2$contrast.x)) * 0.4 + 2
  message(p_file)
  ggsave(p_file, p, width = my_width, height = my_height, limitsize = FALSE)
}

# my_type <- "CD8"
for (my_type in c("CD4", "CD8")) {
  my_res2 <- res %>%
    filter(str_detect(analysis, my_type)) %>%
    filter(
      str_detect(contrast.y, filters[[my_type]]["cell_subset"]),
      str_detect(contrast.x, filters[[my_type]]["cell_cluster"])
    )
  my_res2_mat <- hclust_dataframe(my_res2, "contrast.x", "contrast.y", "rho")
  my_res2$contrast.x <- factor(my_res2$contrast.x, rownames(my_res2_mat))
  my_res2$contrast.y <- factor(my_res2$contrast.y, colnames(my_res2_mat))
  #
  my_res3 <- my_res2 %>%
    group_by(contrast.x) %>%
    top_n(n = 1, wt = rho) %>%
    select(contrast.x, contrast.y)
  my_res3
  #
  for (i in seq(nrow(my_res3))) {
    #
    my_contrast1 <- my_res3$contrast.x[i]
    my_contrast2 <- my_res3$contrast.y[i]
    #
    my_ova1 <- ova %>% filter(contrast2 == my_contrast1)
    my_ova2 <- clean_names(de_ova %>% filter(coef == my_contrast2))
    my_ova <- full_join(my_ova1, my_ova2, by = "ensembl_id")
    my_ova$t.x[is.na(my_ova$t.x)] <- 0
    my_ova$t.y[is.na(my_ova$t.y)] <- 0
    #
    my_cor <- suppressWarnings({
      with(my_ova, cor.test(t.x, t.y, method = "spearman", na.action = "na.omit")) %>%
      parameters %>%
      as.data.frame
    })
    #
    p <- ggplot(my_ova) +
      aes(x = t.x, y = t.y) +
      geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.3) +
      geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.3) +
      geom_abline(intercept = 0, slope = c(1, -1), linewidth = 0.3, alpha = 0.3) +
      rasterise(geom_point(size = 0.1, alpha = 0.3), dpi = 300) +
      labs(
        x = my_contrast1, y = my_contrast2,
        title = with(my_cor, glue("rho = {signif(rho, 2)}, P = {signif(p, 1)}"))
      )
    p_file <- glue("{out_dir}/ova-cor/scatter/scatter - {my_contrast1} - {my_contrast2}.pdf")
    dir.create(dirname(p_file), recursive = TRUE, showWarnings = FALSE)
    message(p_file)
    ggsave(p_file, p, width = 5, height = 4, limitsize = FALSE)
  }
} 

# TODO: redo the limma-cor.pdf heatmap, but:
# - do not show the Tissue cells
# - only show CD4 for our Blood CD4, and CD8 for our Blood CD8

# boxplots {{{

# my_contrast1 <- "Tissue CD8 T cells 11 vs all"

my_contrasts1 <- unique(ova$contrast2)
my_contrasts1 <- my_contrasts1[str_detect(my_contrasts1, "Blood CD")]
for (my_contrast1 in my_contrasts1) {
  my_ova1 <- ova %>% filter(contrast2 == my_contrast1) %>% filter(ave_expr > 0.5)
  my_genes <- (
    my_ova1 %>%
      group_by(contrast2) %>%
      # arrange(-auc) %>%
      arrange(p_value) %>%
      top_n(n = 5, wt = auc) %>%
      ungroup()
  )$gene
  dir.create(glue("{out_dir}/boxplot/{my_contrast1}"), showWarnings = FALSE)
  # my_genes <- c("ICOS", "IL2RA", "FOXP3", "CTLA4", "TYROBP", "TRDC", "NCAM1", "KLRC1", "KIR2DL4", "EOMES", "CX3CR1", "ZNF683", "ITGB2", "ITGAE", "HAVCR2", "NR4A1", "NR4A2", "CD69", "CCL4L2", "CCL4")
  my_ix <- unlist(sapply(my_genes, \(x) which(str_detect(rownames(mat), glue("\\b{x}\\b")))))
  for (my_i in my_ix ) {
    my_gene <- rownames(mat)[my_i]
    meta$my_gene <- mat[my_gene,]
    p <- ggplot(meta) +
      aes(y = fct_rev(cell_subset), x = my_gene, fill = cell_subset) +
      geom_boxplot(alpha = 0.8, linewidth = 0.3) +
      facet_grid(rows = vars(cell_type), space = "free", scales = "free", switch = "y") +
      scale_fill_manual(values = palettes$Giles2022, guide = "none") +
      scale_x_continuous(breaks = pretty_breaks(3)) +
      scale_y_discrete(
        labels = ~ str_split_fixed(.x, "_", 2)[,2],
        position = "right"
      ) +
      labs(x = NULL, y = NULL, title = str_split_fixed(my_gene, "-", 2)[,2]) +
      theme(
        strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(face = "italic"),
        panel.spacing = unit(0.5, "lines")
      )
    p_file <- glue("{out_dir}/boxplot/{my_contrast1}/boxplot-{my_gene}.pdf")
    my_width <- 4
    my_height <- 4
    message(p_file)
    ggsave(p_file, p, width = my_width, height = my_height, limitsize = FALSE)
  }
}

# res %>%
#   filter(rho > 0, p < 0.05) %>%
#   group_by(contrast.x) %>%
#   top_n(n = 1, wt = rho) %>%
#   select(rho, p, contrast.x, contrast.y)

# }}}

# }}}

# Simple Spearman {{{

pb <- qread("paper/pseudobulk_donor_cluster.qs")

mat[1:5,1:5]
meta[1:5,]

mat_genes <- as.data.frame(str_split_fixed(rownames(mat), "-", 2))
colnames(mat_genes) <- c("ensembl_id", "gene")

# Average of each of the 14 cell types
y <- with(meta, model.matrix(~ 0 + x))
y <- t(t(y) / colSums(y))
my_mat1 <- mat %*% y
colnames(my_mat1) <- str_remove(colnames(my_mat1), "^x")
rownames(my_mat1) <- str_split_fixed(rownames(my_mat1), "-", 2)[,1]
my_mat1[1:5,]

# Average of each cell cluster log2cpm
pb$meta$dataset_cluster <- str_replace_all(with(pb$meta, glue("{dataset}_{cluster}")), " ", "_")
pb$meta %>% count(dataset_cluster)
y <- with(pb$meta, model.matrix(~ 0 + dataset_cluster))
y <- t(t(y) / colSums(y))
my_mat2 <- as.matrix(pb$log2cpm) %*% y
colnames(my_mat2) <- str_remove(colnames(my_mat2), "^dataset_cluster")

my_mat1[1:5,1:5]
my_mat2[1:5,1:5]

my_mat1 <- my_mat1[rowMeans(my_mat1) > 0.5,]
my_mat2 <- my_mat2[rowMeans(my_mat2) > 0.5,]

my_genes <- intersect(rownames(my_mat1), rownames(my_mat2))
length(my_genes)

my_mat1 <- my_mat1[my_genes,]
my_mat2 <- my_mat2[my_genes,]

cor2df <- function(m) {
  data.frame(
    row  = rownames(m)[row(m)[upper.tri(m)]],
    col  = colnames(m)[col(m)[upper.tri(m)]],
    corr = m[upper.tri(m)]
  )
}

my_cor <- cor2df(cor(my_mat1, my_mat2, method = "spearman"))
colnames(my_cor) <- c("cell_subset", "cell_cluster", "spearman")
#
res_mat <- my_cor %>%
  pivot_wider(id_cols = "cell_subset", names_from = "cell_cluster", values_from = "spearman") %>%
  column_to_rownames("cell_subset") %>%
  as.matrix
res_h1 <- dendsort(hclust(dist(res_mat)))
res_h2 <- dendsort(hclust(dist(t(res_mat))))
res_mat <- res_mat[res_h1$order, res_h2$order]
# res_mat[1:5,1:5]
#
my_cor$cell_subset <- factor(my_cor$cell_subset, rownames(res_mat))
my_cor$cell_cluster <- factor(my_cor$cell_cluster, colnames(res_mat))

# my_cor %>%
#   group_by(cell_cluster) %>%
#   top_n(n = 1, wt = spearman) %>%
#   filter(str_detect(cell_cluster, "Blood_CD")) %>%
#   arrange(cell_subset, cell_cluster) %>%
#   print(n = 100)

cd8_names <- c(
  "Blood_CD8_T_cells_1"  = "1\\. <i>CX3CR1 FGFBP2 ZNF683</i>",
  "Blood_CD8_T_cells_2"  = "2\\. <i>HLA-DRA GZMK</i>",
  "Blood_CD8_T_cells_3"  = "3\\. <i>IL7R SELL</i>",
  "Blood_CD8_T_cells_4"  = "4\\. GDT <i>KLRC3</i>",
  "Blood_CD8_T_cells_5"  = "5\\. <i>CX3CR1 FGFBP2</i>",
  "Blood_CD8_T_cells_6"  = "6\\. NK <i>FCER1G</i>",
  "Blood_CD8_T_cells_7"  = "7\\. <i>CX3CR1 FGFBP2 TPM4</i>",
  "Blood_CD8_T_cells_8"  = "8\\. <i>HLA-DRA GZMK ZNF683</i>",
  "Blood_CD8_T_cells_9"  = "9\\. CD8 T <i>CX3CR1 FGFBP2 NEAT1</i>",
  "Blood_CD8_T_cells_10" = "10\\. MAIT",
  "Blood_CD8_T_cells_11" = "11\\. <i>IL7R SELL CCR7</i>",
  "Blood_CD8_T_cells_12" = "12\\. Cycling",
  "Blood_CD8_T_cells_13" = "13\\. <i>CX3CR1 FGFBP2 TRBV12-4</i>",
  "Blood_CD8_T_cells_14" = "14\\. NK <i>NCAM1</i>, GDT, NKT",
  #
  "Blood_CD4_T_cells_1"  = "1\\. <i>CD45RA CCR7</i>",
  "Blood_CD4_T_cells_2"  = "2\\. <i>SELL IL7R</i>",
  "Blood_CD4_T_cells_3"  = "3\\. <i>HLA-DRA PDCD1</i>",
  "Blood_CD4_T_cells_4"  = "4\\. <i>LGALS3</i>",
  "Blood_CD4_T_cells_5"  = "5\\. <i>GZMK CD69</i>",
  "Blood_CD4_T_cells_6"  = "6\\. <i>AKT ADAM19</i>",
  "Blood_CD4_T_cells_7"  = "7\\. <i>HLA-DRA FOXP3</i>"
)
filters <- list(
  "CD4" = c("cell_subset" = "cd4", "cell_cluster" = "Blood_CD4"),
  "CD8" = c("cell_subset" = "cd8", "cell_cluster" = "Blood_CD8")
)
for (my_type in names(filters)) {
  my_cor2 <- my_cor %>%
    filter(
      str_detect(cell_subset, filters[[my_type]]["cell_subset"]),
      str_detect(cell_cluster, filters[[my_type]]["cell_cluster"])
    )
  my_cor2_mat <- hclust_dataframe(my_cor2, "cell_subset", "cell_cluster", "spearman")
  my_cor2$cell_subset <- factor(my_cor2$cell_subset, rownames(my_cor2_mat))
  my_cor2$cell_cluster <- factor(my_cor2$cell_cluster, colnames(my_cor2_mat))
  p <- ggplot(my_cor2) +
    aes(x = cell_subset, y = fct_rev(cell_cluster), fill = spearman) +
    scale_fill_gradientn(
      # colors = rev(brewer.pal("RdBu", n = 11)),
      # rescaler = ~ rescale_mid(.x, mid = 0)
      colors = brewer.pal("Reds", n = 9),
      breaks = pretty_breaks(3)
    ) +
    guides(fill = guide_colorbar(barwidth = 13, title = "Spearman\n")) +
    geom_tile() +
    geom_point(
      data = ~ .x %>%
        group_by(cell_cluster) %>%
        top_n(n = 1, wt = spearman),
      size = 5, shape = 21, fill = "white", stroke = 0.3
    ) +
    scale_x_discrete(
      expand = expansion(0, 0),
      labels = ~ toupper(str_replace_all(str_replace_all(.x, "_", " "), "cd[84] ", "")),
      position = "top"
    ) +
    scale_y_discrete(
      expand = expansion(0, 0),
      # labels = ~ str_replace_all(.x, "_", " "),
      labels = cd8_names,
      position = "right"
    ) +
    labs(x = NULL, y = NULL, title = "Correlation with Giles et al 2022") +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 0),
      axis.text.y.right = element_markdown()
    )
  p_file <- glue("{out_dir}/spearman-cor-{my_type}.pdf")
  my_width <- length(unique(my_cor2$cell_subset)) * 0.4 + 4
  my_height <- length(unique(my_cor2$cell_cluster)) * 0.4 + 2
  message(p_file)
  ggsave(p_file, p, width = my_width, height = my_height, limitsize = FALSE)
}

# }}}

# Spearman on cells, then visualize on UMAP {{{

library(rhdf5)
library(Matrix)
source("R/functions/do-analysis.R")
source("R/functions/helpers.R")

h5_file <- "paper/blood-cd8.h5"
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

my_log2cpm <- do_log2cpm(my_counts)
my_log2cpm <- my_log2cpm[rowMeans(my_log2cpm) > 0.05,]
dim(my_log2cpm)

ix <- which(with(meta, str_detect(x, "cd8_[^b]")))
my_meta <- meta[ix,]
my_mat <- mat[,ix]
my_mat_genes <- str_split_fixed(rownames(my_mat), "-", 2)

my_ens <- intersect(rownames(my_log2cpm), my_mat_genes[,1])
length(my_ens)

my_log2cpm <- as.matrix(my_log2cpm[my_ens,])
my_mat <- my_mat[match(my_ens, my_mat_genes[,1]),]

my_meta$x <- factor(my_meta$x)
y <- with(my_meta, model.matrix(~ 0 + x))
y <- t(t(y) / colSums(y))
colnames(y) <- str_remove(colnames(y), "^x")
my_mat <- my_mat %*% y

range(my_mat)

dim(my_log2cpm)
dim(my_mat)

system.time({
my_cor <- cor(my_log2cpm, my_mat, method = "spearman")
})
# centered and scaled correlation values for each cell
my_cor2 <- t(scale(t(my_cor)))

dim(my_cor)

my_top <- colnames(my_cor)[apply(my_cor, 1, which.max)]
sort(table(my_top))

all(my_obs$cell == rownames(my_cor))
my_obs$top_cor <- my_top

my_obs %>%
  left_join(my_meta %>% select(rna_id, cell_subset), by = c("top_cor" = "rna_id")) %>%
  count(leiden, cell_subset) %>%
  group_by(leiden) %>%
  top_n(n = 1, wt = n)

# my_ref <- colnames(my_cor)[1]
for (my_ref in colnames(my_cor)) {
  message(my_ref)
  p <- plot_hexgene(
      x = my_obs$UMAP1,
      y = my_obs$UMAP2,
      z = my_cor[,my_ref],
      italic = FALSE,
      text = FALSE
    ) +
    scale_fill_gradientn(
      colors = scico::scico(palette = "roma", direction = -1, n = 11),
      breaks = scales::pretty_breaks(3),
      limits = range(my_cor),
      rescaler = function(...) rescale_mid(..., mid = 0)
    ) + 
    guides(fill = guide_colorbar(title = "Spearman", barwidth = 10)) +
    labs(title = glue("Correlation with {my_ref}"))
  p_file <- glue("{out_dir}/umap-cor/blood-CD8-umap-spearman-{my_ref}.pdf")
  dir.create(dirname(p_file), showWarnings = FALSE, recursive = TRUE)
  message(p_file)
  ggsave(p_file, p, width = 5, height = 5)
}

# }}}

# Linear Discriminant Analysis (LDA) on CD8 subsets {{{
# 
# # Test all pairs of clusters
# meta$x <- factor(make_clean_names(meta$cell_subset, allow_dupes = TRUE, replace = c("\\+"="p","-"="n")))
# 
# # CD8
# ix <- which(with(meta, str_detect(x, "cd8_[^b]")))
# my_meta <- meta[ix,]
# my_meta$x <- factor(my_meta$x)
# my_mat <- mat[,ix]
# #
# des1 <- with(my_meta, model.matrix(~ x))
# fit1 <- lmFit(object = my_mat, design = des1)
# fit1 <- eBayes(fit1)
# res <- topTable(fit1, number = nrow(my_mat))
# res$signif <- res$P.Value < 0.05 / nrow(res)
# my_genes <- rownames(res)[res$signif]
# length(my_genes)
# my_mat <- my_mat[my_genes,]
# #
# my_ens <- str_split_fixed(my_genes, "-", 2)[,1]
# ix <- which(my_ens %in% rownames(pb$log2cpm))
# my_genes <- my_genes[ix]
# my_ens <- my_ens[ix]
# length(my_genes)
# my_mat <- my_mat[my_genes,]
# 
# my_means <- rowMeans(my_mat)
# my_sds <- matrixStats::rowSds(my_mat)
# 
# # my_pca1 <- summary(prcomp(t(my_mat[my_genes,]), center = TRUE, scale = TRUE))
# my_pca <- summary(prcomp(t(my_mat[my_genes,]), center = my_means, scale = my_sds))
# 
# # res <- MASS::lda(
# #   x        = my_pca$x[,1:10],
# #   # x = t(my_mat[my_genes,]),
# #   grouping = my_meta$x,
# #   CV       = TRUE
# # )
# # # o <- seriation::seriate(res$posterior)
# # # image(res$posterior[o[[1]], o[[2]]])
# # res_dat <- cbind.data.frame(
# #  cell_subset = as.character(res$class),
# #  sample_id = rownames(res$posterior),
# #  res$posterior,
# #  stringsAsFactors = FALSE
# # ) %>% reshape2::melt(id.vars = c("cell_subset", "sample_id"))
# 
# # Average of each cell cluster log2cpm
# pb$meta$dataset_cluster <- str_replace_all(with(pb$meta, glue("{dataset}_{cluster}")), " ", "_")
# pb$meta %>% count(dataset_cluster)
# y <- with(pb$meta, model.matrix(~ 0 + dataset_cluster))
# y <- t(t(y) / colSums(y))
# my_mat2 <- as.matrix(pb$log2cpm %*% y)
# colnames(my_mat2) <- str_remove(colnames(my_mat2), "^dataset_cluster")
# my_mat2 <- as.matrix(my_mat2[,str_detect(colnames(my_mat2), "Blood_CD8")])
# my_mat2 <- my_mat2[my_ens,]
# 
# my_pca2 <- summary(prcomp(t(my_mat2), center = my_means, scale = my_sds))
# 
# my_lda <- MASS::lda(
#   formula = x ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC8 + PC9 + PC10,
#   data = cbind.data.frame(x = my_meta$x, my_pca$x[,1:10])
# )
# 
# my_lda2 <- predict(my_lda, as.data.frame(my_pca2$x[,1:10]))
# 
# res_dat <- cbind.data.frame(
#  cell_subset = as.character(my_lda2$class),
#  cell_cluster = rownames(my_lda2$posterior),
#  my_lda2$posterior,
#  stringsAsFactors = FALSE
# ) %>% reshape2::melt(id.vars = c("cell_subset", "cell_cluster"))
# 
# res_dat %>%
#   group_by(cell_cluster) %>%
#   top_n(n = 2, wt = value) %>%
#   arrange(cell_cluster)
# 
# # x <- reshape2::melt(x, id.vars = c("cell", "class")) %>% as_tibble
# # x$class <- as.character(x$class)
# # x$variable <- as.character(x$variable)
# # ggplot(x %>% filter(class == variable)) +
# #  aes(x = variable, y = value) +
# #  geom_boxplot()
# # x <- seriate_dataframe(x, "cell", "variable", "value")
# # ggplot(x) +
# #  aes(x = cell, y = variable, fill = value) +
# #  geom_tile() +
# #  scale_fill_scico()

# }}}

