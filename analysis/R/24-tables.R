# libraries, scripts {{{
library(anndata)
library(circlize)
library(ComplexHeatmap)
library(conflicted)
library(data.table)
library(dendsort)
library(devtools)
# library(doParallel)
library(dplyr)
# library(foreach)
library(ggbeeswarm)
library(ggforce)
library(ggplot2)
library(ggraph)
library(ggrepel)
library(ggstance)
library(glue)
library(glue)
library(harmony)
library(janitor)
library(janitor)
library(limma)
library(magrittr)
library(Matrix)
library(naturalsort)
library(openxlsx)
library(pals)
library(patchwork)
library(pbapply)
library(pheatmap)
# library(princurve)
library(purrr)
library(qs)
library(qs)
library(RColorBrewer)
library(readxl)
library(rhdf5)
library(scales)
library(scattermore)
library(scico)
library(shadowtext)
library(sitools)
library(stringr)
library(tidygraph)
library(tidyr)
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
count <- dplyr::count
arrange <- dplyr::arrange
conflict_prefer("geom_errorbarh", "ggplot2")
conflict_prefer("read_excel", "readxl")
#
source("R/functions/composition.R")
source("R/plot-composition.R")
source("R/functions/do-pseudobulk-de.R")
source("R/sample-ids.R")
# }}}

# Ensembl to gene symbol {{{
########################################################################
my_h5_files <- Sys.glob(
  "analysis/terra/cellranger-per-channel/output/*/filtered_feature_bc_matrix.h5",
)
genes <- h5read(my_h5_files[1], "matrix/features")
genes <- tibble(ensembl_id = genes$id, symbol = genes$name)
ensembl_to_symbol <- unlist(with(genes, split(symbol, ensembl_id)))
# }}}

# globals {{{
analysis_name_pub <- c(
  "a12_4_4_t4_cd8_1_2"  = "Tissue CD8 T cells",
  "a12_4_4_t4_cd4_2_2"  = "Tissue CD4 T cells",
  "a12_4_4_m3_2"        = "Tissue Myeloid cells",
  "a12_4_4_b5_1_3"      = "Tissue B cells",
  "n3_2"                = "Tissue epithelial cells",
  "blood2_myeloid5"     = "Blood Myeloid cells",
  "blood2_bcell5"       = "Blood B cells",
  "blood2_tcell5_cd4_5" = "Blood CD4 T cells",
  "blood2_tcell5_cd8_5" = "Blood CD8 T cells",
  "luoma_cd45_a5_tcell2_cd8_3" = "Luoma Tissue CD8 T cells",
  "luoma_cd45_a5_tcell2_cd4_3" = "Luoma Tissue CD4 T cells"
)
analyses <- names(analysis_name_pub)
# }}}

# Write expression files

# .h5 files {{{
analysis_name_h5 <- c(
  "a12_4_4_t4_cd8_1_2"  = "paper/tissue-cd8.h5",
  "a12_4_4_t4_cd4_2_2"  = "paper/tissue-cd4.h5",
  "a12_4_4_m3_2"        = "paper/tissue-myeloid.h5",
  "a12_4_4_b5_1_3"      = "paper/tissue-b.h5",
  "n3_2"                = "paper/tissue-epithelial.h5",
  "blood2_myeloid5"     = "paper/blood-myeloid.h5",
  "blood2_bcell5"       = "paper/blood-b.h5",
  "blood2_tcell5_cd4_5" = "paper/blood-cd4.h5",
  "blood2_tcell5_cd8_5" = "paper/blood-cd8.h5",
  "luoma_cd45_a5_tcell2_cd8_3" = "paper/luoma-cd8.h5",
  "luoma_cd45_a5_tcell2_cd4_3" = "paper/luoma-cd4.h5"
)
# TODO: write features like "ENGS1234|SYMBOL"
for (analysis_name in analyses) {
  print_status(analysis_name)
  # analysis_name <- analyses[1]
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  stopifnot(file.exists(a1_file))
  h5_file <- str_replace(a1_file, ".qs$", ".h5")
  if (file.exists(h5_file)) {
    print_status(glue("Found {h5_file}"))
    h5_file_paper <- glue("{analysis_name_h5[analysis_name]}")
    file.copy(h5_file, h5_file_paper, overwrite = TRUE)
  } else {
    a1 <- qread(a1_file)
    a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
    print_status(glue("Writing {h5_file}"))
    unlink(h5_file)
    h5createFile(h5_file)
    print_status(glue("Writing obs, de, mcv, pca, pca_h, knn"))
    h5write(a1$obs %>% mutate(across(where(is.factor), as.character)), h5_file, "obs")
    h5write(a1$counts_stats, h5_file, "genes")
    h5write(a1$de, h5_file, "de")
    h5write(a1$mcv, h5_file, "mcv")
    h5write(a1$pca, h5_file, "pca")
    h5write(a1$pca_h, h5_file, "pca_h")
    h5createGroup(h5_file, "knn")
    h5createGroup(h5_file, "knn/simil_cells")
    h5write(a1$knn$simil_cells@i, h5_file, "knn/simil_cells/indices")
    h5write(a1$knn$simil_cells@p, h5_file, "knn/simil_cells/indptr")
    h5write(a1$knn$simil_cells@x, h5_file, "knn/simil_cells/data")
    h5createGroup(h5_file, "knn/nn_method")
    h5write(a1$knn$nn_method$idx, h5_file, "knn/nn_method/idx")
    h5write(a1$knn$nn_method$dist, h5_file, "knn/nn_method/dist")
    print_status(glue("Writing counts matrix"))
    h5createGroup(h5_file, "matrix")
    h5write(dim(a1$counts), h5_file, "matrix/shape")
    h5write(a1$counts@i, h5_file, "matrix/indices")
    h5write(a1$counts@p, h5_file, "matrix/indptr")
    h5write(a1$counts@x, h5_file, "matrix/data")
    h5write(rownames(a1$counts), h5_file, "matrix/features")
    h5write(colnames(a1$counts), h5_file, "matrix/barcodes")
    # h5createGroup(h5_file, "log2cpm")
    # h5write(dim(a1$log2cpm), h5_file, "log2cpm/shape")
    # h5write(a1$log2cpm@i, h5_file, "log2cpm/indices")
    # h5write(a1$log2cpm@p, h5_file, "log2cpm/indptr")
    # h5write(a1$log2cpm@x, h5_file, "log2cpm/data")
    # h5write(rownames(a1$log2cpm), h5_file, "log2cpm/features")
    # h5write(colnames(a1$log2cpm), h5_file, "log2cpm/barcodes")
    h5closeAll()
  }
}
# }}}

# h5ad files {{{
analysis_name_h5ad <- c(
  "a12_4_4_t4_cd8_1_2"  = "paper/tissue-cd8.h5ad",
  "a12_4_4_t4_cd4_2_2"  = "paper/tissue-cd4.h5ad",
  "a12_4_4_m3_2"        = "paper/tissue-myeloid.h5ad",
  "a12_4_4_b5_1_3"      = "paper/tissue-b.h5ad",
  "n3_2"                = "paper/tissue-epithelial.h5ad",
  "blood2_myeloid5"     = "paper/blood-myeloid.h5ad",
  "blood2_bcell5"       = "paper/blood-b.h5ad",
  "blood2_tcell5_cd4_5" = "paper/blood-cd4.h5ad",
  "blood2_tcell5_cd8_5" = "paper/blood-cd8.h5ad",
  "luoma_cd45_a5_tcell2_cd8_3" = "paper/luoma-cd8.h5ad",
  "luoma_cd45_a5_tcell2_cd4_3" = "paper/luoma-cd4.h5ad"
)
analyses <- names(analysis_name_h5ad)
#
for (analysis_name in analyses) {
  # analysis_name <- analyses[1]
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  h5ad_file <- analysis_name_h5ad[analysis_name]
  #
  stopifnot(file.exists(a1_file))
  # unlink(h5ad_file)
  if (file.exists(h5ad_file)) {
    print_status(glue("Found {h5ad_file}"))
  } else {
    print_status(glue("Reading {a1_file}"))
    a1 <- qread(a1_file)
    # a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
    #
    exclude_cols <- c("unwanted", "donor2", "is_cell", "fdr_not_empty", "pass_qc")
    include_cols <- setdiff(colnames(a1$obs), exclude_cols)
    a1$obs <- as.data.frame(a1$obs)[,include_cols]
    stopifnot("cell" %in% colnames(a1$obs))
    var_names <- sprintf("%s|%s", rownames(a1$counts), ensembl_to_symbol[rownames(a1$counts)])
    #
    print_status(glue("Writing {h5ad_file}"))
    anndata::write_h5ad(
      anndata = anndata::AnnData(
        # X = t(a1$log2cpm),
        # raw = anndata::AnnData(
        #   X = t(a1$counts)
        # ),
        X = t(a1$counts),
        obs = data.frame(a1$obs, row.names = a1$obs$cell),
        var = data.frame(a1$counts_stats, row.names = var_names),
        uns = list(
          "mcv"   = a1$mcv,
          "pca"   = a1$pca,
          "pca_h" = a1$pca_h,
          "knn"   = a1$knn,
          "de"    = a1$de
        )
      ),
      filename = h5ad_file,
      compression = "gzip"
    )
    print_status("done")
  }
}
# }}}

# Write analysis result files

# Cell cluster abundance analysis {{{
########################################################################

ix_luoma <- str_detect(analyses, "luoma")

retval_case <- rbindlist(lapply(analyses, function(analysis_name) {
  if (str_detect(analysis_name, "luoma")) {
    out_dir <- as.character(glue("results/Luoma2020/{analysis_name}/figures/composition-case-vs-control"))
  } else {
    out_dir <- as.character(glue("results/a20/{analysis_name}/figures/composition-case-vs-control"))
  }
  tsv_file <- glue("{out_dir}/masc_1-case-complete.tsv")
  retval <- fread(tsv_file)
  # retval <- fread(file.path(out_dir, "masc_1-case-Case.tsv"))
  retval$analysis <- analysis_name_pub[analysis_name]
  retval
}), fill = TRUE) %>% relocate(analysis, cluster, value, lrt_p, lrt_fdr)
head(retval_case)
fwrite(retval_case, "paper/cluster-abundance_case-vs-control.tsv", sep = "\t")
#
retval_drug <- rbindlist(lapply(analyses[!ix_luoma], function(analysis_name) {
  retval <- fread(glue("results/a20/{analysis_name}/figures/drug/masc_1-case-drug-complete.tsv"))
  # retval <- fread(glue("results/a20/{analysis_name}/figures/drug/masc_drugPD_1_CTLA_4.tsv"))
  retval$analysis <- analysis_name_pub[analysis_name]
  retval
}), fill = TRUE) %>% relocate(analysis, cluster, value, lrt_p, lrt_fdr) %>%
mutate(OR = exp(est), `OR_2.5` = exp(est_low), `OR_97.5` = exp(est_high))
head(retval_drug)
fwrite(retval_drug, "paper/cluster-abundance_drug-within-cases.tsv", sep = "\t")

# retval_control <- rbindlist(lapply(analyses[!ix_luoma], function(analysis_name) {
#   slug <- glue("results/a20/{analysis_name}")
#   retval <- fread(glue("{slug}/figures/drug/masc-control-drugPD_1.tsv"))
#   retval$analysis <- analysis_name_pub[analysis_name]
#   retval
# })) %>% relocate(analysis) %>% select(-p)
# fwrite(retval_control, "paper/cluster-abundance_control.tsv", sep = "\t")

# abundance <- rbind(
#   retval_case %>% mutate(contrast = "Case vs Control"),
#   retval_drug %>% mutate(contrast = "Case PD-1/CTLA-4 vs Case PD-1"),
#   retval_control %>% mutate(contrast = "Control PD-1 vs Control None")
# ) %>% relocate(contrast)
# fwrite(abundance, "paper/cluster-abundance.tsv", sep = "\t")

## Get the labels to paste into the text
#retval_case_labels <- rbindlist(lapply(analyses, function(analysis_name) {
#	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/composition-case-vs-control"))
#  retval <- fread(file.path(out_dir, "masc_1-case-Case.tsv"))
#  retval$analysis <- analysis_name_pub[analysis_name]
#  retval
#}), fill = TRUE) %>% select(analysis, label)
#fwrite(retval_case_labels, "paper/cluster-abundance_case-vs-control_labels.tsv", sep = "\t")
##
#retval_drug <- rbindlist(lapply(analyses[!ix_luoma], function(analysis_name) {
#  retval <- fread(glue("results/a20/{analysis_name}/figures/drug/masc_drugPD_1_CTLA_4.tsv"))
#  retval$analysis <- analysis_name_pub[analysis_name]
#  retval
#}), fill = TRUE) %>% select(analysis, label)
#fwrite(retval_drug, "paper/cluster-abundance_drug-within-cases_labels.tsv", sep = "\t")


#{
#  wb <- openxlsx::createWorkbook()
#  fname <- "paper/abundance.xlsx"
#  unlink(fname)
#  #
#  this_sheet <- "Sheet 1"
#  openxlsx::addWorksheet(wb, this_sheet)
#  openxlsx::writeDataTable(
#    wb,
#    this_sheet,
#    x = abundance %>% mutate_if(is.numeric, signif),
#    rowNames = FALSE, tableStyle = "TableStyleLight1"
#  )
#  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
#}


# }}}

# Summary of the number of differentially expressed genes for each dataset {{{
retval <- rbindlist(lapply(analyses, function(analysis_name) {
  slug <- "a20"
  if (str_detect(analysis_name, "luoma")) {
    slug <- "Luoma2020"
  }
	retval <- fread(as.character(glue(
    "results/{slug}/{analysis_name}/figures/de-case-vs-control/de_summary_case-vs-control.tsv"
  )))
  retval$analysis <- analysis_name_pub[analysis_name]
  retval
})) %>% relocate(analysis)
fwrite(retval, "paper/genes_case-vs-control.tsv", sep = "\t")
# }}}

# de-case-vs-control.tsv.gz - This is updated to use the model with sex {{{
de_case <- rbindlist(lapply(analyses, function(analysis_name) {
  slug <- "a20"
  if (str_detect(analysis_name, "luoma")) {
    slug <- "Luoma2020"
  }
  x_file <- glue(
    "results/{slug}/{analysis_name}/figures/de-case-vs-control/de_donor-sex-case.tsv.gz"
  )
  if (!file.exists(x_file)) {
    x_file <- glue(
      "results/{slug}/{analysis_name}/figures/de-case-vs-control/de_donor_case-vs-control.tsv.gz"
    )
  }
  x <- fread(x_file)
  x$analysis <- analysis_name_pub[analysis_name]
  x$cluster <- "all cells"
  #
  y_file <- glue(
    "results/{slug}/{analysis_name}/figures/de-case-vs-control/de_sex-case.tsv.gz"
  )
  if (!file.exists(y_file)) {
    y_file <- glue(
      "results/{slug}/{analysis_name}/figures/de-case-vs-control/de_case-vs-control.tsv.gz"
    )
  }
  y <- fread(y_file)
  y$analysis <- analysis_name_pub[analysis_name]
  y$GeneName <- NULL
  x <- rbind(x, y)
})) %>% relocate(analysis, cluster, Gene, percent)
de_case <- clean_names(de_case)
de_case <- de_case %>% mutate_if(is.numeric, signif, 4) # saves 50% file size
fwrite(de_case, "paper/de-case-vs-control.tsv.gz", sep = "\t")

#{
#  wb <- openxlsx::createWorkbook()
#  fname <- "paper/de-case-vs-control.xlsx"
#  unlink(fname)
#  #
#  for (this_sheet in sort(unique(de_case$analysis))) {
#    openxlsx::addWorksheet(wb, this_sheet)
#    openxlsx::writeDataTable(
#      wb,
#      this_sheet,
#      x = de_case %>% filter(analysis == this_sheet) %>% mutate_if(is.numeric, signif, 3),
#      rowNames = FALSE, tableStyle = "TableStyleLight1"
#    )
#  }
#  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
#}

# }}}

# de-drug.tsv.gz {{{
de_drug <- rbindlist(lapply(analyses[!str_detect(analyses, "luoma")], function(analysis_name) {
  x <- fread(glue(
    "results/a20/{analysis_name}/figures/drug/de_donor.tsv.gz"
  ))
  x$analysis <- analysis_name_pub[analysis_name]
  x$cluster <- "all cells"
  #
  y <- fread(glue(
    "results/a20/{analysis_name}/figures/drug/de_cluster.tsv.gz"
  ))
  y$analysis <- analysis_name_pub[analysis_name]
  if ("GeneName" %in% colnames(y)) {
    y$GeneName <- NULL
  }
  x <- rbind(x, y, fill = TRUE)
})) %>% relocate(analysis, cluster, Gene, percent)
de_drug <- clean_names(de_drug)
de_drug <- de_drug %>% mutate_if(is.numeric, signif, 3) # saves 50% file size
fwrite(de_drug, "paper/de-drug.tsv.gz", sep = "\t")

de_drug_bx <- fread("paper/de-drug.tsv.gz")
de_drug_bx <- de_drug_bx %>%
  filter(!(
      str_detect(de_drug_bx$analysis, "Blood") &
      contrast != "ControlPD1vsControlNone"
  )) %>%
  filter(
    contrast != "CasePD1vsControlPD1"
  )
de_drug_bx %>% count(analysis, contrast)
fwrite(de_drug_bx, "paper/de-drug-biorxiv.tsv.gz")

#{
#  wb <- openxlsx::createWorkbook()
#  fname <- "paper/de-drug.xlsx"
#  unlink(fname)
#  #
#  for (this_sheet in sort(unique(de_drug$analysis))) {
#    openxlsx::addWorksheet(wb, this_sheet)
#    openxlsx::writeDataTable(
#      wb,
#      this_sheet,
#      x = de_drug %>% filter(analysis == this_sheet) %>% mutate_if(is.numeric, signif, 3),
#      rowNames = FALSE, tableStyle = "TableStyleLight1"
#    )
#  }
#  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
#}

# }}}

# de-contrasts.tsv.gz {{{
x <- fread("paper/de-case-vs-control.tsv.gz")
x$contrast <- "CasevsControl"
y <- fread("paper/de-drug.tsv.gz")
x <- rbind(x, y)
x$contrast <- str_replace(x$contrast, "vs", " vs ")
x <- x[x$contrast != "CasePD1 vs ControlPD1"]
with(x, table(analysis, contrast))
x$cluster <- sprintf("%s-%s",
  x$analysis %>%
  str_replace("Tissue ", "T-") %>%
  str_replace("Blood ", "B-") %>%
  str_replace(" T", "") %>%
  str_replace(" cells", "") %>%
  str_replace(" T", "") %>%
  str_replace("epithelial", "E") %>%
  str_replace("Myeloid", "MP"),
  x$cluster %>% str_replace("all cells", "all")
)
table(x$cluster)
fwrite(x, "paper/de-contrasts.tsv.gz", sep = "\t")
# }}}

# ova.tsv.gz {{{
analyses <- names(analysis_name_pub)
de_ova <- rbindlist(lapply(analyses, function(analysis_name) {
  x <- fread(glue(
    "results/a20/{analysis_name}/figures/pseudobulk_de_ova.tsv.gz"
  ))
  x$analysis <- analysis_name_pub[analysis_name]
  x
})) %>% relocate(analysis, contrast, Gene, auc)
de_ova <- clean_names(de_ova)
de_ova <- de_ova %>% mutate_if(is.numeric, signif, 3) # saves 50% file size
fwrite(de_ova, "paper/ova.tsv.gz", sep = "\t")

#{
#  wb <- openxlsx::createWorkbook()
#  fname <- "paper/de-ova.xlsx"
#  unlink(fname)
#  #
#  for (this_sheet in sort(unique(de_ova$analysis))) {
#    openxlsx::addWorksheet(wb, this_sheet)
#    openxlsx::writeDataTable(
#      wb,
#      this_sheet,
#      x = de_ova %>% filter(analysis == this_sheet) %>% mutate_if(is.numeric, signif, 3),
#      rowNames = FALSE, tableStyle = "TableStyleLight1"
#    )
#  }
#  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
#}

# # Do we see TFs with differential expression for a cluster?
# jaspar <- fread("/projects/external_data/jaspar.genereg.net/jaspar-ensembl.tsv")
# de_ova %>%
#   filter(
#     analysis == "Tissue CD8 T cells",
#     ensembl_id %in% jaspar$ensembl_id,
#     # contrast %in% c("5 vs all", "1 vs all"),
#     contrast %in% c("11 vs all"),
#     adj_p_val < 0.05
#   ) %>% arrange(-auc) %>%
#   head(12)

# }}}

# ccc-spearman-cell-lineages.tsv.gz {{{
cc1 <- qread("results/a20/covarying-abundance/cc_level1.qs")$correlations
cell_key <- list("CT" = "CD8", "T" = "CD4")
cc1 <- cc1 %>% dplyr::mutate(c1 = recode(c1, !!!cell_key), c2 = recode(c2, !!!cell_key))
cc1 <- cc1 %>% dplyr::rename(
  percent1 = p1, percent2 = p2, mean1 = m1, mean2 = m2
)
cc1 <- cc1 %>% dplyr::mutate(
  id = glue("{c2} {g2} {c1} {g1}")
)
fwrite(cc1, "paper/ccc-spearman-cell-lineages.tsv.gz", sep = "\t")
# }}}

# ccc-limma-cell-lineages.tsv.gz {{{

# # TODO: Why is Myeloid absent from these results?
# cc2 <- qread("results/a20/covarying-abundance/cc_level2.qs")$correlations
# cc2$c1 <- str_replace(cc2$c1, "^CT", "CD8-")
# cc2$c1 <- str_replace(cc2$c1, "^T", "CD4-")
# cc2$c1 <- str_replace(cc2$c1, "^B", "B-")
# cc2$c1 <- str_replace(cc2$c1, "^E", "E-")
# cc2 <- cc2 %>% dplyr::rename(
#   percent1 = p1, percent2 = p2, mean1 = m1, mean2 = m2
# )
# cc2 <- cc2 %>% dplyr::mutate(
#   id = glue("{c2} {g2} {c1} {g1}")
# )
# cc2 <- cc2 %>% mutate_if(is.numeric, function(x) sprintf("%.5g", x))
# fwrite(cc2, "paper/ccc-spearman-cell-clusters.tsv.gz", sep = "\t")

gg1 <- fread("results/a20/covarying-abundance/ccc-sums-level1.tsv")
cell_key <- list("CT" = "CD8", "T" = "CD4")
gg1 <- gg1 %>% dplyr::mutate(c1 = recode(c1, !!!cell_key), c2 = recode(c2, !!!cell_key))
gg1 <- gg1 %>% dplyr::rename(
  percent1 = p1, percent2 = p2
)
gg1 <- gg1 %>% dplyr::mutate(
  id = glue("{c2} {g2} {c1} {g1}")
)
fwrite(gg1, "paper/ccc-limma-cell-lineages.tsv.gz", sep = "\t")

# gg1 <- clean_names(qread("results/a20/covarying-abundance/gg_level1_case.qs")$limma)
# cell_key <- list("CT" = "CD8", "T" = "CD4")
# gg1 <- gg1 %>% dplyr::mutate(c1 = recode(c1, !!!cell_key), c2 = recode(c2, !!!cell_key))
# gg1 <- gg1 %>% dplyr::rename(
#   percent1 = p1, percent2 = p2
# )
# fwrite(gg1, "paper/ccc-sum.tsv.gz", sep = "\t")

# {
#   wb <- openxlsx::createWorkbook()
#   fname <- "paper/ccc.xlsx"
#   unlink(fname)
#   {
#     this_sheet <- "Sums of gene pairs"
#     openxlsx::addWorksheet(wb, this_sheet)
#     openxlsx::writeDataTable(
#       wb,
#       this_sheet,
#       x = gg1 %>% mutate_if(is.numeric, signif),
#       rowNames = FALSE, tableStyle = "TableStyleLight1"
#     )
#   }
#   {
#     this_sheet <- "Spearman of gene pairs"
#     openxlsx::addWorksheet(wb, this_sheet)
#     openxlsx::writeDataTable(
#       wb,
#       this_sheet,
#       x = cc1 %>% mutate_if(is.numeric, signif),
#       rowNames = FALSE, tableStyle = "TableStyleLight1"
#     )
#   }
#   openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
# }

# }}}

# packages.tsv {{{
x <- devtools::session_info()$packages
x <- x %>% select(package, loadedversion, date, source)
fwrite(x, "paper/packages.tsv", sep = "\t")
# }}}

# summary.tsv, summary-dge.tsv, summary-abundance.tsv {{{
retval <- rbindlist(lapply(analyses[!ix_luoma], function(analysis_name) {
	x <- read_excel(as.character(glue(
    "results/a20/{analysis_name}/figures/de-case-vs-control/de_donor_case-vs-control.xlsx"
  )))
  n_case_up <- with(x, sum(adj.P.Val < 0.05 & logFC > log2(1.5)))
  n_case_dn <- with(x, sum(adj.P.Val < 0.05 & logFC < log2(1 / 1.5)))
  #
  y <- read_excel(glue("results/a20/{analysis_name}/figures/drug/de_donor.xlsx"))
  n_drug_up <- with(y, sum(adj.P.Val < 0.05 & logFC > log2(1.5)))
  n_drug_dn <- with(y, sum(adj.P.Val < 0.05 & logFC < log2(1 / 1.5)))
  #
	z <- fread(glue("results/a20/{analysis_name}/figures/composition-case-vs-control/masc_1-case-Case.tsv"))
  n_cluster_case_up <- with(z, sum(OR > 1 & lrt_fdr < 0.05))
  n_cluster_case_dn <- with(z, sum(OR < 1 & lrt_fdr < 0.05))
  #
  w <- fread(glue("results/a20/{analysis_name}/figures/drug/masc_drugPD_1_CTLA_4.tsv"))
  n_cluster_drug_up <- with(w, sum(OR > 1 & lrt_fdr < 0.05))
  n_cluster_drug_dn <- with(w, sum(OR < 1 & lrt_fdr < 0.05))
  #
  # data.frame(list(
  #   analysis = analysis_name_pub[analysis_name],
  #   n_genes_case_up = n_case_up,
  #   n_genes_case_dn = n_case_dn,
  #   n_cluster_case_up = n_cluster_case_up,
  #   n_cluster_case_dn = n_cluster_case_dn,
  #   n_genes_drug_up = n_drug_up,
  #   n_genes_drug_dn = n_drug_dn,
  #   n_cluster_drug_up = n_cluster_drug_up,
  #   n_cluster_drug_dn = n_cluster_drug_dn
  # ))
  data.frame(list(
    analysis      = analysis_name_pub[analysis_name],
    case_genes    = glue("{n_case_up}, {n_case_dn}"),
    case_clusters = glue("{n_cluster_case_up}, {n_cluster_case_dn}"),
    drug_genes    = glue("{n_drug_up}, {n_drug_dn}"),
    drug_clusters = glue("{n_cluster_drug_up}, {n_cluster_drug_dn}")
  ))
}))
fwrite(retval, "paper/summary.tsv", sep = "\t")

{
  de_case <- fread("paper/de-case-vs-control.tsv.gz")
  de_drug <- fread("paper/de-drug.tsv.gz")
  n_case <- de_case %>% filter(cluster == "all cells") %>%
    group_by(analysis) %>%
    summarize(n = sum(abs(log_fc) > log2(1.5) & adj_p_val < 0.05 & percent > 1), .groups = "drop") %>%
    mutate(contrast = "CasevsControl")
  n_drug <- de_drug %>% filter(cluster == "all cells") %>%
    group_by(analysis, contrast) %>%
    summarize(n = sum(abs(log_fc) > log2(1.5) & adj_p_val < 0.05 & percent > 1), .groups = "drop") %>%
    filter(contrast != "CasePD1vsControlPD1")
  d <- rbind(n_case, n_drug)
  d <- d %>% pivot_wider(names_from = contrast, values_from = n)
}
fwrite(d, "paper/summary-dge.tsv", sep = "\t")

de_drug %>%
  filter(
    contrast == "CasePD1CTLA4vsCasePD1",
    analysis == "Blood CD8 T cells",
    cluster == "all cells",
    log_fc > 0,
    percent > 1
  ) %>%
  arrange(p_value) %>% head

de_drug %>%
  filter(
    contrast == "CasePD1CTLA4vsCasePD1",
    analysis == "Tissue epithelial cells",
    cluster == "all cells",
    log_fc > 0,
    percent > 1
  ) %>%
  arrange(p_value) %>% head


{
  retval_case <- fread("paper/cluster-abundance_case-vs-control.tsv")
  retval_case$contrast <- "Case vs Control"
  #
  retval_drug <-  fread("paper/cluster-abundance_drug-within-cases.tsv")
  retval_drug$contrast <- "Case PD-1/CTLA-4 vs Case PD-1"
  #
  # retval_control <- fread("paper/cluster-abundance_control.tsv")
  # retval_control$contrast <- "Control PD-1 vs Control None"
  #
  # d <- rbind(retval_case, retval_drug, retval_control)
  d <- rbind(retval_case, retval_drug, fill = TRUE)
  #
  d <- d %>%
    filter(value %in% c("drugPD-1/CTLA-4", "Case")) %>%
    group_by(analysis, contrast) %>%
    summarize(n = sum(lrt_fdr < 0.05)) %>%
    pivot_wider(names_from = contrast, values_from = n)
}
fwrite(d, "paper/summary-abundance.tsv", sep = "\t")

# read_excel("paper/abundance.xlsx", sheet = 1) %>%
#   # filter(analysis == "Blood CD8 T cells") %>%
#   # filter(analysis == "Blood CD4 T cells") %>%
#   # filter(analysis == "Tissue CD4 T cells") %>%
#   # filter(analysis == "Tissue Myeloid cells") %>%
#   filter(analysis == "Blood B cells") %>%
#   # filter(analysis == "Tissue epithelial cells") %>%
#   # mutate(label = glue("cluster {cluster}, OR = {signif(OR,2)}, P = {signif(lrt_p, 1)}")) %>%
#   mutate(label = glue("cluster {cluster}, OR = {signif(OR,2)}, 95% CI {signif(OR_2.5,2)} to {signif(OR_97.5,2)}")) %>%
#   select(analysis, label, lrt_fdr)

# }}}

# How to read the .h5 files {{{
if (FALSE) {

  library(rhdf5)

  a1 <- qread("paper/pseudobulk_donor.qs")
  h5_file <- "paper/pseudobulk_donor.h5"
  unlink(h5_file)
  h5createFile(h5_file)
  h5createGroup(h5_file, "matrix")
  h5write(a1$meta %>% mutate(across(where(is.factor), as.character)), h5_file, "obs")
  h5write(dim(a1$log2cpm), h5_file, "matrix/shape")
  h5write(a1$log2cpm@i, h5_file, "matrix/indices")
  h5write(a1$log2cpm@p, h5_file, "matrix/indptr")
  h5write(a1$log2cpm@x, h5_file, "matrix/data")
  h5write(rownames(a1$log2cpm), h5_file, "matrix/features")
  h5write(colnames(a1$log2cpm), h5_file, "matrix/barcodes")
  h5closeAll()

  # Read it like this:
  h5 <- rhdf5::h5read(h5_file, "matrix")
  counts <- Matrix::sparseMatrix(
    dims   = h5$shape,
    i      = as.numeric(h5$indices),
    p      = as.numeric(h5$indptr),
    x      = as.numeric(h5$data),
    index1 = FALSE
  )


  library(rhdf5)

  a1 <- qread("paper/pseudobulk_donor_cluster.qs")
  h5_file <- "paper/pseudobulk_donor_cluster.h5"
  unlink(h5_file)
  h5createFile(h5_file)
  h5createGroup(h5_file, "matrix")
  h5write(a1$meta %>% mutate(across(where(is.factor), as.character)), h5_file, "obs")
  h5write(dim(a1$log2cpm), h5_file, "matrix/shape")
  h5write(a1$log2cpm@i, h5_file, "matrix/indices")
  h5write(a1$log2cpm@p, h5_file, "matrix/indptr")
  h5write(a1$log2cpm@x, h5_file, "matrix/data")
  h5write(rownames(a1$log2cpm), h5_file, "matrix/features")
  h5write(colnames(a1$log2cpm), h5_file, "matrix/barcodes")
  h5closeAll()

  h5read(h5_file, "obs") %>% as.data.table

}
# }}}

