
# libraries {{{
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(conflicted)
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
count <- dplyr::count
arrange <- dplyr::arrange
conflict_prefer("geom_errorbarh", "ggplot2")
conflict_prefer("select", "dplyr")
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
fwrite(genes %>% arrange(symbol), "data/ensembl_id-symbol.tsv", sep = "\t")
ensembl_to_symbol <- unlist(with(genes, split(symbol, ensembl_id)))
#
rename_by_size <- function(x) {
  sizes <- sort(table(x), decreasing = TRUE)
  x_old_names <- names(sizes)
  x_new_names <- seq_along(x_old_names)
  names(x_new_names) <- x_old_names
  return(x_new_names[x])
}
# tissue cd8
recluster_cd8_leiden122 <- function(cd8_clusters) {
  # First round
  cd8_clusters[cd8_clusters %in% c("2", "6")] <- "2.6"
  cd8_clusters[cd8_clusters %in% c("1", "12")] <- "1.12"
  cd8_clusters[cd8_clusters %in% c("8", "9")] <- "8.9"
  #
  return(rename_by_size(cd8_clusters))
}
# tissue cd8
recluster_cd8_leiden151 <- function(cd8_clusters) {
  # First round
  cd8_clusters[cd8_clusters %in% c("1", "4", "13")] <- "1.4.13"
  cd8_clusters[cd8_clusters %in% c("10", "11", "15")] <- "10.11.15"
  cd8_clusters[cd8_clusters %in% c("2", "9")] <- "2.9"
  cd8_clusters[cd8_clusters %in% c("7", "16")] <- "7.16"
  #
  return(rename_by_size(cd8_clusters))
}
recluster_blood2_tcell5_cd4_5_leiden18 <- function(xs) {
  stopifnot(all(1:16 %in% xs))
  # First round
  xs[xs %in% c("1", "2", "4", "11", "12")] <- "1.2.4.11.12"
  xs[xs %in% c("3", "7", "14", "9")] <- "3.7.14.9"
  xs[xs %in% c("6", "16")] <- "6.16"
  xs[xs %in% c("8", "15")] <- "8.15"
  #
  return(rename_by_size(xs))
}
recluster_blood2_myeloid5_leiden108 <- function(xs) {
  stopifnot(all(1:15 %in% xs))
  # First round
  xs[xs %in% c("2", "11")] <- "2.11"
  #
  return(rename_by_size(xs))
}
recluster_blood2_bcell5_leiden151 <- function(xs) {
  stopifnot(all(1:12 %in% xs))
  # First round
  xs[xs %in% c("9", "10", "11", "4")] <- "9.10.11.4"
  xs[xs %in% c("1", "2", "5", "6", "7", "8")] <- "1.2.5.6.7.8"
  #
  return(rename_by_size(xs))
}
get_cluster_groups <- function(analysis_name) {
  cluster_groups <- list()
  if (analysis_name == "n3_2") {
    cluster_groups[["g1 Immature epithelial cells"]] <- c("8", "3", "14")
    cluster_groups[["g2 Absorptive epithelial cells"]] <- c("1", "5", "6", "16", "18", "9", "12")
    cluster_groups[["g3 Mature, absorptive epithelial cells"]] <- c("2", "11", "20")
    cluster_groups[["g4 Secretory cells"]] <- c("4", "10", "17", "15", "19")
    cluster_groups[["g5 Mesenchymal cells"]] <- c("7", "23", "21", "22", "13")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_t4_cd8_1_2") {
    cluster_groups[["g1 ITGB2"]] <- c("3", "11")
    # cluster_groups[["g2 ITGAE ID3 CD8 Trm"]] <- c("1", "4", "5", "9")
    # cluster_groups[["g3 ITGAE ZNF683 CD8 Trm"]] <- c("7", "6", "2")
    cluster_groups[["g2 ITGAE"]] <- c("1", "4", "5", "9", "7", "6", "2")
    cluster_groups[["g3 Other"]] <- c("8", "10", "12")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_t4_cd4_2_2") {
    cluster_groups[["g1 Treg"]] <- c("8", "9", "5")
    cluster_groups[["g2"]] <- c("1", "2", "3", "4")
    cluster_groups[["g4 CXCL13"]] <- c("6", "7", "10")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "blood2_tcell5_cd8_5") {
    cluster_groups[["g1 MHC-II"]] <- c("2", "12", "8")
    cluster_groups[["g2 Innate"]] <- c("4", "14", "6", "10")
    cluster_groups[["g3 CX3CR1"]] <- c("7", "13", "1", "9")
    cluster_groups[["g4 IL7R"]] <- c("3", "11")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_m3_2") {
    cluster_groups[["g1"]] <- c("1", "6", "5", "2", "3")
    cluster_groups[["g2"]] <- c("7", "4")
    cluster_groups[["g3"]] <- c("8")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  else if (analysis_name == "a12_4_4_b5_1_3") {
    cluster_groups[["b"]] <- c("2", "3", "5", "11", "12")
    cluster_groups[["plasma"]] <- c("1", "4", "6", "7", "8", "9", "10", "13")
    cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  }
  unlist(cluster_groups)
}

# }}}

# Summary of assays {{{
########################################################################

# retval <- rbindlist(lapply(analyses, function(analysis_name) {
#   a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
#   print_status(glue("Reading {a1_file}"))
#   a1 <- qread(a1_file)
#   my_cols <- intersect(colnames(a1$obs), c("donor", "case", "has_tcr", "has_bcr"))
#   retval <- a1$obs %>% group_by(across(all_of(my_cols))) %>%
#     summarise(n = n(), .groups = 'drop_last')
#   retval$analysis <- analysis_name
#   retval
# }), fill = TRUE)
# retval %>% count(analysis)

# setdiff(t2$patient, retval$donor)
# setdiff(retval$donor, t2$patient)

# naturalsort(unique((retval %>% filter(analysis == "a12_4_4_m3_2"))$donor)) %>% length

# naturalsort(unique((t2 %>% filter(name == "scRNA-seq"))$patient))

# t2$patient[
#   t2$tissue_immune_cell_sc_rna_seq_performed == "Yes" &
#   t2$sc_nuc_seq_performed == "No"
# ]
# retval %>% filter(donor %in% c("SIC_71", "SIC_19", "SIC_94"))

t2 <- clean_names(read_excel("data/Table_S2_scRNA-seq_snRNA-seq_summary.xlsx"))
my_cols <- c(
"tissue_immune_cell_sc_rna_seq_performed",
"tissue_sc_tcr_seq_performed",
"tissue_sc_bcr_seq_performed",
"sc_nuc_seq_performed",
"blood_sc_rna_seq_performed",
"blood_sc_tcr_seq_performed",
"blood_sc_bcr_seq_performed",
"blood_cite_seq_performed",
"tissue_microscopy",
"blood_luminex")
t2 <- t2 %>% select(
  -biopsy_location_of_immune_cell_gex_libraries,
  -location_of_biopsy_for_sc_nuc_seq
) %>% pivot_longer(all_of(my_cols))
t2$name <- t2$name %>% str_replace("_seq_performed", "")
table(t2$name)
t2$name[t2$name == "sc_nuc"] <- "Tissue Non-immune snRNA-seq"
t2$name[t2$name == "tissue_immune_cell_sc_rna"] <- "Tissue Immune scRNA-seq"
t2$name[t2$name == "tissue_sc_bcr"] <- "Tissue Immune BCR"
t2$name[t2$name == "tissue_sc_tcr"] <- "Tissue Immune TCR"
t2$name[t2$name == "blood_sc_tcr"] <- "Blood Immune TCR"
t2$name[t2$name == "blood_sc_bcr"] <- "Blood Immune BCR"
t2$name[t2$name == "blood_sc_rna"] <- "Blood Immune scRNA-seq"
t2$name[t2$name == "blood_cite"] <- "Blood Immune CITE-seq"
t2$name[t2$name == "blood_luminex"] <- "Blood Whole Luminex"
t2$name[t2$name == "tissue_microscopy"] <- "Tissue Whole Microscopy"
# t2 <- t2[t2$name != "blood_cite",]
t2$name <- naturalfactor(t2$name)
t2$colitis_status <- str_replace(t2$colitis_status, " - ", "\n")
ix <- t2$value == "Yes" & t2$colitis_status == "irColitis Case"
t2$value[ix] <- "irColitis Case"
ix <- t2$value == "Yes" & t2$colitis_status == "Control\nHealthy"
t2$value[ix] <- "Control\nHealthy"
ix <- t2$value == "Yes" & t2$colitis_status == "Control\nOn ICI Therapy"
t2$value[ix] <- "Control\nOn ICI"
ix <- t2$colitis_status == "Control\nOn ICI Therapy"
t2$colitis_status[ix] <- "Control\nOn ICI"
t2$name <- factor(t2$name, c(
  "Blood Whole Luminex",
  "Blood Immune BCR",
  "Blood Immune TCR",
  "Blood Immune CITE-seq",
  "Blood Immune scRNA-seq",
  "Tissue Immune BCR",
  "Tissue Immune TCR",
  "Tissue Immune scRNA-seq",
  "Tissue Non-immune snRNA-seq",
  "Tissue Whole Microscopy"
))
t2 <- t2 %>% left_join(
    t2 %>% select(patient, colitis_status) %>% unique %>%
    count(colitis_status) %>%
    mutate(colitis_status_n = as.character(glue("{colitis_status}\n(n = {n})"))),
    by = "colitis_status"
  )
x <- str_split_fixed(t2$name, " ", 3)
t2$tissue <- str_remove(paste(x[,1], x[,2]), " $")
t2$name <- x[,3]
t2$name <- factor(
  as.character(t2$name),
  rev(c("Microscopy", "CITE-seq", "scRNA-seq", "snRNA-seq", "TCR", "BCR", "Luminex"))
)
t2$pub <- pub_ids[t2$patient]
stopifnot(!any(is.na(t2$pub)))
blood_donors <- c("MC_2", "SIC_100", "SIC_109", "SIC_126", "SIC_13", "SIC_134",
  "SIC_140", "SIC_141", "SIC_172", "SIC_186", "SIC_188", "SIC_19",
  "SIC_31", "SIC_32", "SIC_36", "SIC_40", "SIC_43", "SIC_71", "SIC_89",
  "SIC_94", "SIC_97")
ix_remove <- !(t2$patient %in% blood_donors) & t2$tissue == "Blood Immune"
t2 <- t2[!ix_remove,]
t2$tissue <- factor(
  t2$tissue,
  rev(c("Blood Whole", "Blood Immune", "Tissue Immune", "Tissue Non-immune", "Tissue Whole"))
)
t2$patient <- factor(
  t2$patient,
  (t2 %>% filter(value != "No") %>% count(patient) %>% arrange(-n))$patient
)
p <- ggplot(t2) +
  aes(x = patient, y = name, fill = value) +
  geom_tile(color = "white", size = 0.15) +
  scale_fill_manual(
    values = c(
      "No"               = "white",
      "irColitis Case"   = okabe(8)[2],
      # "Control\nHealthy" = okabe(8)[3],
      # "Control\nOn ICI"  = okabe(8)[4]
      "Control\nHealthy" = "grey50",
      "Control\nOn ICI"  = "grey20"
    ),
    guide = "none"
  ) +
  facet_grid(
    tissue ~ colitis_status_n, scales = "free", space = "free"
  ) +
  scale_y_discrete(position = "l", expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle = 90, size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.text.x = element_text(hjust = 0.5, vjust = 0)
  ) +
  labs(
    y = NULL,
    x = glue("{length(unique(t2$patient))} Patients")
  )
my_ggsave(
  "assay-table",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 8,
  height = 3.5,
  units = "in", dpi = 300
)

library(tidyr)
library(purrr)
conflict_prefer("extract", "tidyr")
t1 <- clean_names(read_excel("data/Table_S1_patient_profiles.xlsx"))
t1 <- t1[!is.na(t1$patient_id),]
colnames(t1)
t1$ctla4 <- str_detect(replace_na(t1$anti_ctla_4_dose, ""), "mg")
t1$pd1 <- (
  str_detect(replace_na(t1$anti_pd_1_dose, ""), "mg") |
  str_detect(replace_na(t1$anti_pd_l1_dose, ""), "mg") 
)
t1$pd1_only <- t1$pd1 & !t1$ctla4
t1$ctla4_only <- !t1$pd1 & t1$ctla4
t1$dual <- t1$ctla4 & t1$pd1
# t1 <- t1[c("anti_ctla_4_dose","anti_pd_1_dose","anti_pd_l1_dose", "ctla4", "pd1", "dual")]
t1$mel <- str_detect(t1$indication_for_immunotherapy, "Squamous cell cancer") |
 str_detect(t1$indication_for_immunotherapy, "Melanoma")
t1$doses <- rowSums(t1 %>%
  mutate(z = extract_numeric(number_checkpoint_cycles_prior_to_symptom_onset)) %>%
  extract(
    number_checkpoint_cycles_prior_to_symptom_onset,
    c("x", "y"), "([\\d.]+) dose.+ ([\\d.]+) dose", convert = TRUE) %>%
  select(x, y, z) %>%
  as.matrix,
  na.rm = TRUE
)
t1$grade <- readr::parse_number(t1$colitis_grade)
t1$weeks <- readr::parse_number(t1$duration_from_last_cycle_to_symptom_onset_weeks, na = c("-", "NA"))
my_cols <- c(
  "case_vs_control",
  "patient_age",
  "sex",
  "mel",
  "pd1_only",
  "ctla4_only",
  "dual",
  "doses",
  "weeks",
  "grade",
  "mayo_colitis_score"
)
t1 <- t1[,my_cols]
my_render_cont <- function(x) {
    with(stats.apply.rounding(stats.default(x), digits=2), c("",
        "Mean (SD)"=sprintf("%s (±%s)", MEAN, SD)))
}
library(table1)
t1_out <- knitr::kable(
  table1::table1(
    ~ . | case_vs_control,
    data = t1,
    overall = FALSE,
    render.continuous = my_render_cont
  )
)
write_file(paste(as.character(t1_out), collapse = "\n"), "paper/table1.txt")

# }}}

# Table S4. Quantitation of CD45+ immune cells and CD45+ CD3+ T cells from colon mucosal tissue {{{
facs <- clean_names(fread("data/facs.tsv"))
#
# x <- facs$case_control == 'irColitis Case'
# y <- facs$avg_number_cells_biopy_x10_5 * 1e5
# my_t <- t.test(y ~ x)
# my_ci <- (my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
# my_est <- my_t$estimate[2] / my_t$estimate[1]
# glue("({signif(my_est, 2)}-fold, 95% CI {signif(my_ci[2], 2)} to {signif(my_ci[1], 2)}, p = {signif(my_t$p.value, 1)})")
#
cells <- data.frame(
  # "total" = readr::parse_number(facs$total_live_trypan_cells_from_pooled_biopsies_hemocytometer_count),
  "case" = facs$case_control,
  "total" = facs$avg_number_cells_biopy_x10_5 * 1e5,
  "cd45" = facs$avg_cd45_cells_biopsy_x10_5 * 1e5,
  "cd3" = facs$avg_cd45_cd3_cells_biopsy_x10_5 * 1e5
)
#
cor1 <- cor.test(cells$total, cells$cd45)
cor2 <- cor.test(cells$total, cells$cd45 / cells$total)
#
p1 <- ggplot(cells) +
  aes(x = total, y = cd45, fill = case) +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  annotate(
    geom = "text",
    x = 0,
    y = Inf,
    label = glue("R^2 = {signif(cor1$estimate ^ 2, 2)} p = {signif(cor1$p.value, 1)}"),
    hjust = -0.1,
    vjust = 2
  ) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_log10(labels = label_number_si()) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
p2 <- ggplot(cells) +
  aes(x = total, y = cd45 / total, fill = case) +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  annotate(
    geom = "text",
    x = 0,
    y = Inf,
    label = glue("R^2 = {signif(cor2$estimate ^ 2, 2)} p = {signif(cor2$p.value, 1)}"),
    hjust = -0.1,
    vjust = 2
  ) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_continuous(labels = label_percent()) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    )
  ) +
  annotation_logticks(sides = "b")
p <- p1 / p2
my_ggsave(
  "facs-cells",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 7,
  height = 5,
  units = "in", dpi = 300
)

cells$case2 <- ifelse(str_detect(cells$case, "Case"), "Case", "Control")
# my_t <- t.test(y = log10(cells$total), x = cells$case2 == "Case")
my_t <- t.test(log2(cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
p_total <- ggplot(cells) +
  aes(x = total, y = case2, fill = case, group = paste(case2, total)) +
  geom_point(
    shape = 21, size = 3, stroke = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  annotate(
    geom = "text",
    x = 1e6,
    y = 1.5,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_log10(labels = label_number_si()) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(x = NULL, y = NULL)
p_scatter <- ggplot(cells) +
  aes(x = total, y = cd45 / total, fill = case) +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_continuous(labels = label_percent()) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  labs(x = "Mean cells per biopsy", y = "% CD45+ EPCAM-")
my_t <- t.test(log2(cells$cd45 / cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
p_cd45 <- ggplot(cells) +
  aes(x = case2, y = cd45 / total, fill = case, group = paste(case2, total)) +
  geom_point(
    shape = 21, size = 3, stroke = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  annotate(
    geom = "text",
    x = 1.5,
    y = 0.75,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_y_continuous(labels = label_percent()) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    )
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = NULL, y = NULL)
p <- p_scatter + p_total + p_cd45 + plot_layout(design =
"B#
AC", widths = c(1, 0.3, 0.3), heights = c(0.3, 1, 1))
my_ggsave(
  "facs-cells-cd45",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 10,
  height = 4,
  units = "in", dpi = 300
)

cells$case2 <- ifelse(str_detect(cells$case, "Case"), "Case", "Control")
# my_t <- t.test(y = log10(cells$total), x = cells$case2 == "Case")
my_t <- t.test(log2(cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
p_total <- ggplot(cells) +
  aes(x = total, y = case2, fill = case, group = paste(case2, total)) +
  geom_point(
    shape = 21, size = 3, stroke = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  annotate(
    geom = "text",
    x = 1e6,
    y = 1.5,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_log10(labels = label_number_si()) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(x = NULL, y = NULL)
p_scatter <- ggplot(cells) +
  aes(x = total, y = cd3 / total, fill = case) +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  scale_x_log10(labels = label_number_si()) +
  scale_y_continuous(labels = label_percent()) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  labs(x = "Mean cells per biopsy", y = "% CD45+ CD3+ EPCAM-")
my_t <- t.test(log2(cells$cd3 / cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
p_cd3 <- ggplot(cells) +
  aes(x = case2, y = cd3 / total, fill = case, group = paste(case2, total)) +
  geom_point(
    shape = 21, size = 3, stroke = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  annotate(
    geom = "text",
    x = 1.5,
    y = 0.55,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_y_continuous(labels = label_percent()) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control on ICI Therapy"  = "grey20"
    )
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = NULL, y = NULL)
p <- p_scatter + p_total + p_cd3 + plot_layout(design =
"B#
AC", widths = c(1, 0.3, 0.3), heights = c(0.3, 1, 1))
my_ggsave(
  "facs-cells-cd3",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 10,
  height = 4,
  units = "in", dpi = 300
)

cells$id <- 1:nrow(cells)
cells$case2 <- ifelse(str_detect(cells$case, "Case"), "Case", "Control")
# my_t <- t.test(y = log10(cells$total), x = cells$case2 == "Case")
my_t <- t.test(log2(cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
signif(2 ^ as.numeric(my_t[c("estimate", "conf.low", "conf.high")]), 2)
p_total <- ggplot(cells) +
  aes(x = total, y = case2, fill = case, group = naturalfactor(paste(case2, id))) +
  geom_point(
    shape = 21, size = 4, stroke = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(aes(group = case2), fun = "mean", shape = 3, size = 1) +
  annotate(
    geom = "text",
    x = 1e6,
    y = 1.5,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_log10(labels = label_number_si(), limits = c(5e3, 2e6)) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      "Control Healthy" = "grey70",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  theme(
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "Mean cells per biopsy")
my_t <- t.test(log2(cells$cd45) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
signif(2 ^ as.numeric(my_t[c("estimate", "conf.low", "conf.high")]), 2)
p_total_cd45 <- ggplot(cells) +
  aes(x = cd45, y = case2, fill = case, group = naturalfactor(paste(case2, id))) +
  geom_point(
    shape = 21, size = 4, stroke = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(aes(group = case2), fun = "mean", shape = 3, size = 1) +
  annotate(
    geom = "text",
    x = 1e6,
    y = 1.5,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_log10(labels = label_number_si(), limits = c(5e3, 2e6)) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      "Control Healthy" = "grey70",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  theme(
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "CD45+ EPCAM- cells per biopsy")
my_t <- t.test(log2(cells$cd3) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
signif(2 ^ as.numeric(my_t[c("estimate", "conf.low", "conf.high")]), 2)
p_total_cd3 <- ggplot(cells) +
  aes(x = cd3, y = case2, fill = case, group = naturalfactor(paste(case2, id))) +
  geom_point(
    shape = 21, size = 4, stroke = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(aes(group = case2), fun = "mean", shape = 3, size = 1) +
  annotate(
    geom = "text",
    x = 1e6,
    y = 1.5,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_log10(labels = label_number_si(), limits = c(5e3, 2e6)) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      "Control Healthy" = "grey70",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  annotation_logticks(sides = "b") +
  theme(
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "CD45+ CD3+ EPCAM- cells per biopsy")
my_t <- t.test(log2(cells$cd45 / cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
signif(2 ^ as.numeric(my_t[c("estimate", "conf.low", "conf.high")]), 2)
p_cd45 <- ggplot(cells) +
  aes(y = case2, x = cd45 / total, fill = case, group = naturalfactor(paste(case2, id))) +
  geom_point(
    shape = 21, size = 4, stroke = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(aes(group = case2), fun = "mean", shape = 3, size = 1) +
  annotate(
    geom = "text",
    y = 1.5,
    x = 0.75,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1)) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      "Control Healthy" = "grey70",
      "Control on ICI Therapy"  = "grey20"
    ),
    guide = "none"
  ) +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "% CD45+ EPCAM- cells")
my_t <- t.test(log2(cells$cd3 / cells$total) ~ as.integer(cells$case2 != "Case"))
my_t <- as.list(broomExtra::tidy(my_t))
signif(2 ^ as.numeric(my_t[c("estimate", "conf.low", "conf.high")]), 2)
p_cd3 <- ggplot(cells) +
  aes(y = case2, x = cd3 / total, fill = case, group = naturalfactor(paste(case2, id))) +
  geom_point(
    shape = 21, size = 4, stroke = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(aes(group = case2), fun = "mean", shape = 3, size = 1) +
  annotate(
    geom = "text",
    y = 1.5,
    x = 0.75,
    label = glue("{signif(2 ^ my_t$estimate, 2)}-fold\np = {signif(my_t$p.value, 1)}"),
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1)) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "irColitis Case"  = okabe(8)[2],
      "Control Healthy" = "grey70",
      "Control on ICI Therapy"  = "grey20"
    )
  ) +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "% CD45+ CD3+ EPCAM- cells")
p <- p_total + p_total_cd45 + p_total_cd3 + p_cd45 + p_cd3 + plot_layout(
  design = "A
B
C
D
E"
)
my_ggsave(
  "facs-cells3",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 9,
  height = 6,
  units = "in", dpi = 300
)

#  [1] "patient_id"
#  [2] "simple_pt_id"
#  [3] "case_control"
#  [4] "total_number_of_biopsy_pieces_enzymatically_digested"
#  [5] "total_live_trypan_cells_from_pooled_biopsies_hemocytometer_count"
#  [6] "avg_number_cells_biopy_x10_5"
#  [7] "percent_dapi_cells_by_facs"
#  [8] "percent_cd45_epcam_cells_by_facs_of_total_dapi_cd235a"
#  [9] "percent_cd45_cd3_cells_by_facs_of_total_dapi_cd235a"
# [10] "avg_cd45_cells_biopsy_x10_5"
# [11] "avg_cd45_cd3_cells_biopsy_x10_5"
facs <- facs %>%
  select(c("case_control", "patient_id", contains("total_dapi"))) %>%
  pivot_longer(contains("percent_"))
colnames(facs) <- c("sample_type", "patient_id", "name", "value")
facs$sample_type[facs$sample_type == "Control on ICI Therapy"] <- "Control ICI"
facs$sample_type[facs$sample_type == "Healthy Control"] <- "Control Healthy"
facs$name[facs$name == "percent_cd45_epcam_cells_by_facs_of_total_dapi_cd235a"] <- "CD45+ cells"
facs$name[facs$name == "percent_cd45_cd3_cells_by_facs_of_total_dapi_cd235a"] <- "CD45+ CD3+ cells"
facs$sample_type <- factor(facs$sample_type)
facs$value <- facs$value / 100
#
facs_text <- data.frame(
  name = c(
    # "CD45+ cells",
    "CD45+ cells",
    # "CD45+ CD3+ cells",
    "CD45+ CD3+ cells"
  ),
  pval = c(
    # {
    #   ix <- facs$sample_type %in% c("Control ICI", "irColitis Case") & facs$name == "CD45+ cells"
    #   wilcox.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    # },
    {
      ix <- facs$sample_type %in% c("Control Healthy", "Control ICI", "irColitis Case") & facs$name == "CD45+ cells"
      t.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    },
    # {
    #   ix <- facs$sample_type %in% c("Control ICI", "irColitis Case") & facs$name == "CD45+ CD3+ cells"
    #   wilcox.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    # },
    {
      ix <- facs$sample_type %in% c("Control Healthy", "Control ICI", "irColitis Case") & facs$name == "CD45+ CD3+ cells"
      t.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    }
  ),
  nudge_x = c(
    # 0.33 / 2,
    0
  ),
  nudge_y = c(
    # 0.4,
    0.8
  )
)
#
#
ix <- facs$sample_type != "irColitis Case" & facs$name == "CD45+ cells"
wilcox.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "Control ICI"))
my_t <- t.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "Control ICI"))
my_ci <- (my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
my_est <- my_t$estimate[1] / my_t$estimate[2]
glue("Control ICI vs Control Healthy: lower abundance of CD45+ cells ({signif(my_est, 2)}-fold, 95% CI {signif(my_ci[2], 2)} to {signif(my_ci[1], 2)}, p = {signif(my_t$p.value, 1)})")
#
ix <- facs$name == "CD45+ cells"
wilcox.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
my_t <- t.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
my_ci <- (my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
my_est <- my_t$estimate[1] / my_t$estimate[2]
glue("irColitis Case vs Control: greater abundance of CD45+ cells ({signif(my_est, 2)}-fold, 95% CI {signif(my_ci[2], 2)} to {signif(my_ci[1], 2)}, p = {signif(my_t$p.value, 1)})")
#
ix <- facs$sample_type != "irColitis Case" & facs$name == "CD45+ CD3+ cells"
wilcox.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "Control ICI"))
my_t <- t.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "Control ICI"))
my_ci <- (my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
my_est <- my_t$estimate[1] / my_t$estimate[2]
glue("Control ICI vs Control Healthy: lower abundance of CD45+ CD3+ cells ({signif(my_est, 2)}-fold, 95% CI {signif(my_ci[2], 2)} to {signif(my_ci[1], 2)}, p = {signif(my_t$p.value, 1)})")
#
ix <- facs$name == "CD45+ CD3+ cells"
wilcox.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
my_t <- t.test(log2(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
my_ci <- (my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
my_est <- my_t$estimate[1] / my_t$estimate[2]
glue("irColitis Case vs Control: greater abundance of CD45+ CD3+ cells ({signif(my_est, 2)}-fold, 95% CI {signif(my_ci[2], 2)} to {signif(my_ci[1], 2)}, p = {signif(my_t$p.value, 1)})")
#
#
facs_mean <- facs %>% group_by(sample_type, name) %>% summarize(mean = mean(value))
#
facs$name <- factor(facs$name, c("CD45+ cells", "CD45+ CD3+ cells"))
#
p <- ggplot() +
  geom_quasirandom(
    data = facs,
    mapping = aes(x = name, y = value, fill = sample_type, group = sample_type),
    shape = 21, size = 5, stroke = 0.3,
    dodge.width = 1, width = 0.1
  ) +
	stat_summary(
		data = facs,
    mapping = aes(x = name, y = value, fill = sample_type, group = sample_type),
    shape = 3,
    fun = 'mean',
    position = position_dodge(width = 1)
	) +
	scale_y_continuous(
    limits = c(-0.1, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = label_percent()
  ) +
  # scale_y_log10(labels = label_percent()) +
  # annotation_logticks(sides = "l", size = 0.3) +
  geom_text(
    data = facs_text,
    mapping = aes(x = name, y = 0.1, label = signif(pval, 2)),
    nudge_x = facs_text$nudge_x,
    nudge_y = facs_text$nudge_y,
    size = 5
  ) +
  annotate(
    geom = "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 0, ymax = -Inf,
    fill = "grey90"
  ) +
  geom_text(
    data = facs_mean,
    mapping = aes(x = name, y = -0.08, label = sprintf("%s%%", signif(100 * mean, 2)), group = sample_type),
    position = position_dodge(width = 1),
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control ICI"  = "grey20"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Percent of cells")
my_ggsave(
  "facs",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 5,
  height = 3,
  units = "in", dpi = 300
)



facs <- clean_names(fread("data/simplified_total_cell_counts_per_biopsy.csv"))
facs <- facs %>% pivot_longer(contains("biopsy"))
facs$sample_type[facs$sample_type == "Control on ICI"] <- "Control ICI"
facs$sample_type[facs$sample_type == "Healthy Control"] <- "Control Healthy"
facs$name[facs$name == "immune_cells_x10_5_biopsy"] <- "CD45+ cells"
facs$name[facs$name == "t_cells_x10_5_biopsy"] <- "CD45+ CD3+ cells"
facs$sample_type <- factor(facs$sample_type)
#
facs_text <- data.frame(
  name = c(
    "CD45+ cells",
    "CD45+ cells",
    "CD45+ CD3+ cells",
    "CD45+ CD3+ cells"
  ),
  pval = c(
    {
      ix <- facs$sample_type %in% c("Control ICI", "irColitis Case") & facs$name == "CD45+ cells"
      wilcox.test(facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    },
    {
      ix <- facs$sample_type %in% c("Control Healthy", "irColitis Case") & facs$name == "CD45+ cells"
      wilcox.test(facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    },
    {
      ix <- facs$sample_type %in% c("Control ICI", "irColitis Case") & facs$name == "CD45+ CD3+ cells"
      wilcox.test(facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    },
    {
      ix <- facs$sample_type %in% c("Control Healthy", "irColitis Case") & facs$name == "CD45+ CD3+ cells"
      wilcox.test(facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))$p.value
    }
  ),
  nudge_x = c(
    0.33 / 2,
    0
  ),
  nudge_y = c(
    0.4,
    0.8
  )
)
#
ix <- facs$name == "CD45+ cells"
wilcox.test(1e5 * facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
my_t <- t.test(1e5 * facs[ix,]$value ~ as.integer(facs[ix,]$sample_type != "irColitis Case"))
(my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
my_t$estimate[1] / my_t$estimate[2]
#
ix <- facs$name == "CD45+ CD3+ cells"
wilcox.test(1e5 * facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
my_t <- t.test(1e5 * facs[ix,]$value ~ as.integer(facs[ix,]$sample_type != "irColitis Case"))
(my_t$conf.int + my_t$estimate[2]) / my_t$estimate[2]
my_t$estimate[1] / my_t$estimate[2]
#
facs$name <- factor(facs$name, c("CD45+ cells", "CD45+ CD3+ cells"))
#
p <- ggplot() +
  geom_quasirandom(
    data = facs,
    mapping = aes(x = name, y = value * 1e5, fill = sample_type, group = sample_type),
    shape = 21, size = 5, stroke = 0.3,
    dodge.width = 1, width = 0.1
  ) +
  scale_y_log10(labels = label_number_si()) +
  annotation_logticks(sides = "l", size = 0.3) +
  geom_text(
    data = facs_text,
    mapping = aes(x = name, y = 1e6, label = signif(pval, 2)),
    nudge_x = facs_text$nudge_x,
    nudge_y = facs_text$nudge_y,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control ICI"  = "grey20"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Cells per biopsy")
my_ggsave(
  "facs",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 5,
  height = 3,
  units = "in", dpi = 300
)

# ix <- facs$name == "CD45+ cells"
# t.test(1e5 * facs[ix,]$value ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
# f1 <- lm(log10(facs[ix,]$value) ~ as.integer(facs[ix,]$sample_type == "irColitis Case"))
# f1
# xs <- signif(10 ^ c(f1$coefficients[2], as.vector(confint(f1)[2,])), 3)
# glue("{xs[1]} (CI {xs[2]}-{xs[3]})")

facs$case <- str_split_fixed(facs$sample_type, " ", 2)[,1]
facs$logvalue <- log10(facs$value * 1e5)
my_t <- list()
ix <- facs$name == "CD45+ cells"
my_t[["CD45+ cells"]] <- t.test(logvalue ~ case, data = facs[ix,])
ix <- facs$name == "CD45+ CD3+ cells"
my_t[["CD45+ CD3+ cells"]] <- t.test(logvalue ~ case, data = facs[ix,])
my_t
#
facs_text <- data.frame(
  name = c(
    "CD45+ cells",
    "CD45+ CD3+ cells"
  ),
  pval = c(
    my_t[["CD45+ cells"]]$p.value,
    my_t[["CD45+ CD3+ cells"]]$p.value
  ),
  nudge_x = c(
    0,
    0
  ),
  nudge_y = c(
    0.4,
    0.4
  )
)
#
est1 <- 10 ^ (my_t[[1]]$estimate[2] - my_t[[1]]$estimate[1])
est_low1 <- 10 ^ (-my_t[[1]]$conf.int[2])
est_high1 <- 10 ^ (-my_t[[1]]$conf.int[1])
sprintf("%.0f (%.0f to %.0f)", est1, est_low1, est_high1)
#
est1 <- 10 ^ (my_t[[2]]$estimate[2] - my_t[[2]]$estimate[1])
est_low1 <- 10 ^ (-my_t[[2]]$conf.int[2])
est_high1 <- 10 ^ (-my_t[[2]]$conf.int[1])
sprintf("%.0f (%.0f to %.0f)", est1, est_low1, est_high1)
#
facs$name <- factor(facs$name, c("CD45+ cells", "CD45+ CD3+ cells"))
#
p <- ggplot() +
  geom_quasirandom(
    data = facs,
    mapping = aes(x = name, y = value * 1e5, fill = sample_type, group = sample_type),
    shape = 21, size = 5, stroke = 0.3,
    dodge.width = 1, width = 0.1
  ) +
  scale_y_log10(labels = label_number_si()) +
  annotation_logticks(sides = "l", size = 0.3) +
  geom_text(
    data = facs_text,
    mapping = aes(x = name, y = 1e6, label = signif(pval, 1)),
    hjust = 0.5,
    nudge_x = facs_text$nudge_x,
    nudge_y = facs_text$nudge_y,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "irColitis Case"  = okabe(8)[2],
      # "Control Healthy" = okabe(8)[3],
      # "Control ICI"     = okabe(8)[4]
      "Control Healthy" = "grey50",
      "Control ICI"  = "grey20"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Cells per biopsy")
my_ggsave(
  "facs",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 5,
  height = 3,
  units = "in", dpi = 300
)

# }}}

# Luminex {{{
########################################################################
lum <- fread("data/luminex_replicates_clean_v1.csv")
my_cols <- c("Patient ID", "Collection Date", "Condition")
lum <- lum %>% pivot_longer(cols = !my_cols)
colnames(lum) <- c("patient", "date", "condition", "name", "value")
lum$condition[lum$condition == "Case irColitis - Follow Up"] <- "irColitis Case Followup"
# Average replicates
lum <- lum %>%
  group_by(patient, date, condition, name) %>%
  summarize(value = mean(value), .groups = "drop")
lum$name <- factor(
  lum$name,
  (lum %>% group_by(name) %>% summarize(mean = sum(value > 0)) %>% arrange(-mean))$name
)
lum_label <- lum %>% select(patient, condition) %>% unique
lum_label$label <- append_n(lum_label$condition)
lum <- left_join(lum, lum_label %>% select(condition, label) %>% unique, by = "condition")
p <- ggplot(lum) +
  aes(x = value, y = label, group = patient, fill = label) +
  geom_point(size = 3, shape = 21, stroke = 0.2, position = position_dodge(width = 0.5)) +
  scale_x_log10(breaks = log_breaks(5)) +
  scale_fill_manual(
    values = c("grey20", "grey60", pals::okabe(8)[2:8])
  ) +
  guides(
    fill = guide_legend(reverse = TRUE, title = "", override.aes = list(size = 5))
  ) +
  facet_wrap(vars(name), ncol = 5, scales = "free_x") +
  labs(
    x = "Analyte concentration (pg/mL)",
    title = glue("{length(unique(lum$name))} factors from serum measured with Luminex")
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines")
  )
my_ggsave(
  "luminex",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 17,
  height = 16,
  units = "in", dpi = 300
)

lum$log_value <- log10(lum$value)
lum$log_value[!is.finite(lum$log_value)] <- NA
lum_test <- lum %>%
  filter(condition != "irColitis Case Followup") %>%
  mutate(case = str_detect(condition, "Case"))
lum_exclude <- as.character(unique((c(
  as.character((
    lum_test %>% filter(!is.na(log_value)) %>%
      count(name, case) %>% arrange(n) %>% filter(n < 5)
  )$name),
  as.character((
    lum_test %>% filter(!is.na(log_value)) %>%
      group_by(name) %>%
      summarize(n_cases = length(unique(case))) %>% filter(n_cases == 1)
  )$name),
  "IL-31"
))))
lum_include <- setdiff(unique(lum$name), lum_exclude)
lum_t <- rbindlist(lapply(lum_include, function(this_name) {
  retval <- broom::tidy(with(
    lum_test %>% filter(name == this_name),
    t.test(log_value ~ case)
  ))
  retval$name <- this_name
  retval
})) %>% relocate(name, p.value)
lum_t <- lum_t %>% arrange(p.value)
lum_t$fdr <- p.adjust(lum_t$p.value, method = "fdr")

lum_select <- lum %>%
  filter(
    name %in% c("MIG (CXCL9)", "IP-10 (CXCL10)", "BLC (CXCL13)", "IL-17A (CTLA-8)"),
    condition %in% c("Control  On ICI Therapy", "Healthy Control", "irColitis Case")
  ) %>%
  left_join(lum_t %>% select(name, p.value), by = "name") %>%
  mutate(name = glue("{name}\np={signif(p.value, 1)}"))
p <- ggplot() +
  geom_point(
    data = lum_select,
    mapping = aes(x = value, y = label, group = patient, fill = label),
    size = 3, shape = 21, stroke = 0.2, position = position_dodge(width = 0.5)
  ) +
  scale_x_log10(breaks = log_breaks(5)) +
  annotation_logticks(sides = "b") +
  scale_fill_manual(
    values = c("grey20", "grey60", pals::okabe(8)[2:8])
  ) +
  guides(
    fill = guide_legend(reverse = TRUE, title = "", override.aes = list(size = 5))
  ) +
  facet_wrap(vars(name), ncol = 5, scales = "free_x") +
  labs(
    x = "Analyte concentration (pg/mL)",
    title = glue("{length(unique(lum_select$name))} factors from serum measured with Luminex"),
    subtitle = "Two-sided t-test p-values"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines")
  )
my_ggsave(
  "luminex-select",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 18,
  height = 3,
  units = "in", dpi = 300
)

# }}}

# Merge clusters {{{
########################################################################
analyses <- c(
  "a12_4_4_t4_cd8_1_2",
  "a12_4_4_t4_cd4_2_2",
  "a12_4_4_m3_2",
  "a12_4_4_b5_1_3",
  "n3_2",
  "blood2_myeloid5",
  "blood2_bcell5",
  "blood2_tcell5_cd4_5",
  "blood2_tcell5_cd8_5"
)
  # "blood2",
  # "blood2_tcell5",
  # "blood2_bcell4",
for (analysis_name in analyses) {
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  stopifnot(file.exists(a1_file))
  a1 <- qread(a1_file)
  switch(analysis_name,
    "a12_4_4_t4_cd8_1_2" = {
      a1$obs[["leiden"]] <- recluster_cd8_leiden151(a1$obs[["leiden1.51"]])
    },
    "blood2_myeloid5" = {
      a1$obs[["leiden"]] <- recluster_blood2_myeloid5_leiden108(a1$obs[["leiden1.08"]])
    },
    "blood2_tcell5" = {
      a1$obs[["leiden"]] <- a1$obs[["leiden1.08"]]
    },
    "blood2_tcell5_cd8_5" = {
      a1$obs[["leiden"]] <- a1$obs[["leiden1.51"]]
    },
    "blood2_tcell5_cd4_5" = {
      a1$obs[["leiden"]] <- recluster_blood2_tcell5_cd4_5_leiden18(a1$obs[["leiden1.8"]])
    },
    "blood2_bcell4" = {
      a1$obs[["leiden"]] <- a1$obs[["leiden1.51"]]
    },
    "blood2_bcell5" = {
      a1$obs[["leiden"]] <- recluster_blood2_bcell5_leiden151(a1$obs[["leiden1.51"]])
    },
    #
    "luoma_cd3_a5" = {
        a1$obs[["leiden"]] <- a1$obs[["leiden0.933"]]
    },
    "luoma_cd45_a5_tcell2_cd8_3" = {
        a1$obs[["leiden"]] <- a1$obs[["leiden0.933"]]
    },
    "luoma_cd45_a5_tcell2_cd8_4" = {
        a1$obs[["leiden"]] <- a1$obs[["leiden0.933"]]
    },
    "luoma_cd45_a5_tcell2_cd4_3" = {
        a1$obs[["leiden"]] <- a1$obs[["leiden0.933"]]
    },
    #
    {
      a1$obs[["leiden"]] <- a1$obs[["leiden0.933"]]
    }
  )
  message(analysis_name)
  print(table(a1$obs[["leiden"]]))
  qsave(a1, a1_file)
}

# }}}

# OVA and AVA volcanos {{{
########################################################################

# analysis_name <- "a12_4_4_t4_cd8_1_2"
# analysis_name <- "a12_4_4_t4_cd4_2_2"
# analysis_name <- "a12_4_4_m3_2"
# analysis_name <- "a12_4_4_b5_1_3"
# analysis_name <- "n3_2"

for (analysis_name in analyses) {
  message(analysis_name)
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  stopifnot(file.exists(a1_file))
  a1 <- qread(a1_file)
  out_dir <- glue("results/a20/{analysis_name}/figures")
  a1$logcpm <- do_log1p_cpm(a1$counts)
  #
  # Compute pseudobulk OVA and AVA
  de <- do_de(a1$counts, a1$obs[["leiden"]], a1$obs$donor, min_mean = 0)
  #
  # de$pb$logcpm[1:5,1:5]
  # de$pb$obs[1:5,]
  #
  de$ava <- de$ava %>%
    dplyr::mutate(Gene = ensembl_to_symbol[feature]) %>%
    dplyr::rename(ensembl_id = feature) %>%
    dplyr::relocate(Gene, ensembl_id, contrast)
  #
  # Compute AVA auc
  # ix_genes <- rownames(a1$logcpm) %in% de$ava$ensembl_id
  # x <- rbindlist(pbapply::pblapply(unique(de$ava$contrast), function(this_contrast) {
  #   groups <- as.vector(str_split_fixed(this_contrast, " vs ", 2))
  #   ix <- a1$obs[["leiden"]] %in% groups
  #   x <- presto::wilcoxauc(X = a1$logcpm[ix_genes,ix], y = a1$obs[["leiden"]][ix] == groups[1])
  #   x$contrast <- this_contrast
  #   x
  # }))
  # x <- x[x$group == "TRUE",]
  # de$ava <- left_join(
  #   de$ava, x[,c("feature", "contrast", "auc")],
  #   by = c("ensembl_id" = "feature", "contrast")
  # )
  # de$ava <- de$ava %>% dplyr::relocate(Gene, ensembl_id, contrast, auc)
  write_de_xlsx(
    d = de$ava %>%
      # group_by(contrast) %>% top_n(n = 2000, wt = abs(auc - 0.5)) %>% ungroup %>%
      group_by(contrast) %>% top_n(n = 2000, wt = -log10(P.Value)) %>% arrange(contrast, P.Value) %>%
      ungroup %>% mutate_if(is.numeric, signif, 5),
    fname = glue("{out_dir}/pseudobulk_de_ava.xlsx"),
    col = "contrast"
  )
  fwrite(de$ava, glue("{out_dir}/pseudobulk_de_ava.tsv.gz"), sep = "\t")
  #
  de$ova <- de$ova %>%
    dplyr::mutate(Gene = ensembl_to_symbol[feature]) %>%
    dplyr::rename(ensembl_id = feature) %>%
    dplyr::relocate(Gene, ensembl_id, contrast)
  # Compute OVA auc
  de$ova_auc <- presto::wilcoxauc(X = a1$logcpm, y = a1$obs[["leiden"]])
  de$ova_auc[["contrast"]] <- sprintf("%s vs all", de$ova_auc$group)
  de$ova <- left_join(
    de$ova, de$ova_auc[,c("feature", "contrast", "auc")],
    by = c("ensembl_id" = "feature", "contrast" = "contrast")
  )
  de$ova <- de$ova %>% dplyr::relocate(Gene, ensembl_id, contrast, auc)
  write_de_xlsx(
    d = de$ova %>%
      group_by(contrast) %>% top_n(n = 2000, wt = auc) %>% arrange(contrast, -auc) %>%
      ungroup %>% mutate_if(is.numeric, signif, 5),
    fname = glue("{out_dir}/pseudobulk_de_ova.xlsx"),
    col = "contrast"
  )
  fwrite(de$ova, glue("{out_dir}/pseudobulk_de_ova.tsv.gz"), sep = "\t")
  #
  # Plot a volcano for each AVA
  # for (this_contrast in unique(de$ava$contrast)) {
  #   top1 <- de$ava %>% filter(contrast == this_contrast)
  #   p <- plot_limma_volcano(top1) +
  #     labs(title = this_contrast)
  #   my_ggsave(
  #     glue("volcano-{str_replace_all(this_contrast, ' ', '_')}"),
  #     out_dir = glue("{out_dir}/cluster_volcano/ava"),
  #     plot = p,
  #     type = "pdf",
  #     scale = 1, width = 5, height = 4, units = "in", dpi = 300
  #   )
  # }
  #
  #de_ava <- fread(glue("{out_dir}/pseudobulk_de_ava.tsv.gz"))
  #for (this_coef in unique(de_ava$coef)) {
  #  top1 <- de_ava %>%
  #    filter(coef == this_coef) %>%
  #    # mutate(g1 = str_split_fixed(coef, " vs ", 2)[,1]) %>%
  #    # mutate(g2 = str_split_fixed(coef, " vs ", 2)[,2]) %>%
  #    # filter(g1 == 1 | g2 == 1) %>%
  #    # filter(g1 == 7 | g2 == 7) %>%
  #    mutate(Gene = ID)
  #  p <- plot_limma_volcano(top1) +
  #    labs(title = this_coef)
  #  my_ggsave(
  #    glue("volcano-{str_replace_all(this_coef, ' ', '_')}"),
  #    #out_dir = glue("figures/{analysis_name}/cluster_volcano"),
  #    out_dir = glue("{out_dir}/cluster_volcano/ava"),
  #    plot = p,
  #    type = "pdf",
  #    scale = 1, width = 8, height = 6, units = "in", dpi = 300
  #  )
  #}
}

# }}}

# Full figures with gene AUC, logistic OR, and per-donor composition {{{
# Figure 2
########################################################################

# analysis_name <- "a12_4_4_t4_cd8_1_2"
# analysis_name <- "a12_4_4_t4_cd4_2_2"
# analysis_name <- "a12_4_4_m3_2"
# analysis_name <- "a12_4_4_b5_1_3"
# analysis_name <- "n3_2"
# analyses <- c(
#   "a12_4_4_t4_cd8_1_2",
#   "a12_4_4_t4_cd4_2_2",
#   "a12_4_4_m3_2",
#   "a12_4_4_b5_1_3",
#   "n3_2",
#   "blood2",
#   "blood2_bcell4",
#   "blood2_myeloid5",
#   "blood2_tcell5"
# )

analyses <- c(
  # "a12_4_4_t4_cd8_1_2",
  # "a12_4_4_t4_cd4_2_2",
  # "a12_4_4_m3_2",
  # "a12_4_4_b5_1_3",
  # "n3_2",
  # "blood2_myeloid5",
  "blood2_bcell5",
  "blood2_tcell5_cd4_5"
  # "blood2_tcell5_cd8_5"
)

cluster_names <- clean_names(read_excel("data/cluster_names.xlsx"))
colnames(cluster_names) <- c("analysis", "number", "name", "top3", "markers", "extra")
# cluster_names$name <- str_replace(cluster_names$name, " - ", ". ")
# cluster_names$name <- str_replace_all(cluster_names$name, ", ", " ")
cluster_names$name <- with(cluster_names, glue("{number}. {name}"))
# cluster_names$name %>% str_replace("CD8 T ", "")

# for (analysis_name in analyses) {
#   out_dir <- glue("results/a20/{analysis_name}/figures")
#   res3_file <- glue("{out_dir}/composition.qs")
#   unlink(res3_file)
# }

# analysis_name <- "luoma_cd3_a5"
# analyses <- c("luoma_cd45_a5_tcell2_cd8_3", "luoma_cd45_a5_tcell2_cd4_3")
# analyses <- c("luoma_cd45_a5_tcell2_cd8_4")
# analysis_name <- analyses[2]
# for (analysis in analyses) {
#   system(glue("ln -s ../Luoma2020/{analysis} results/a20/{analysis}"))
#   system(glue("mkdir results/a20/{analysis}/data"))
#   system(glue("ln -s ../{analysis}.qs results/a20/{analysis}/data/{analysis}.qs"))
# }

print(analyses)

for (analysis_name in analyses) {

  params <- list(
    min_cells_in_cluster = 50,
    min_percent_of_cells_with_gene = 5
  )
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  if ("drug" %in% colnames(a1$obs)) {
    a1$obs$drug <- factor(a1$obs$drug, c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4"))
  }
  a1$analysis_name <- analysis_name
  print_status(glue("done"))
  #
  out_dir <- glue("results/a20/{analysis_name}/figures")
  source("R/functions/composition.R")
  source("R/plot-composition.R")
	#
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/composition-case-vs-control"))
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  #
  if (analysis_name == "luoma_cd3_a5") {
    a1$obs$leiden <- a1$obs$leiden0.933
    a1$obs$cluster <- a1$obs$leiden0.933
  }
	#
	sample_info <- janitor::clean_names(read_excel(
			path = "data/luoma_villani_combined_colitis_samples2.xlsx",
			sheet = 3
	))
	donor_to_sex <- unlist(split(sample_info$sex, sample_info$donor))
  a1$obs$sex <- donor_to_sex[a1$obs$donor]
  a1$obs$sex[a1$obs$donor == "SIC_33"] <- "M"
  stopifnot(all(!is.na(a1$obs$sex)))
  my_obs <- a1$obs
  #
  # my_drug <- c("None", "PD-1", "PD-1/CTLA-4")
  # my_obs <- my_obs %>% filter(drug %in% my_drug)
  # my_obs$drug <- factor(my_obs$drug, my_drug)
  # ix <- with(
  #   my_obs, case == "Case" | (case == "Control" & drug == "PD-1") |
  #     (case == "Control" & drug == "None")
  # )
  # my_obs <- my_obs[ix,]
  # Skip cell clusters that have too few cells.
  exclude_clusters <- (
    my_obs %>% dplyr::count(leiden) %>% filter(n < params[["min_cells_in_cluster"]])
  )$leiden
  #
  my_obs <- my_obs[!my_obs$leiden %in% exclude_clusters,]
  #
  my_obs <- as_tibble(my_obs)
  my_obs$case <- factor(my_obs$case, c("Control", "Case"))
  # with(my_obs, table(sex, case))
  if (analysis_name %in% cluster_names$analysis) {
    num_to_name <- unlist(with(
      cluster_names %>% dplyr::filter(analysis == analysis_name),
      split(name, as.character(number))
    ))
  } else {
    num_to_name <- levels(xd$group)
    names(num_to_name) <- as.character(levels(xd$group))
  }
  a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))

  # We can guess the sex if needed.
  # sapply(unique(a1$obs$donor), function(this_donor) {
  #   ix <- a1$obs$donor == this_donor
  #   x <- rowSums(a1$log2cpm[,ix])
  #   m <- x[names(which(ensembl_to_symbol == "RPS4Y1"))]
  #   f <- x[names(which(ensembl_to_symbol == "XIST"))]
  #   m / f
  # })
  # setdiff(a1$obs$donor, names(donor_to_sex))

  # Gating on single-cell gene expression does not work
  # if (analysis_name == "a12_4_4_t4_cd8_1_2") {
  #   my_symbols <- c("IL17A", "IL26", "CXCL13")
  #   # my_symbols <- c("IL26", "CXCL13")
  #   my_ids <- names(ensembl_to_symbol[which(ensembl_to_symbol %in% my_symbols)])
  #   a1$obs$triple <- colSums(a1$log2cpm[my_ids,] > 0)
  #   a1$obs %>% count(cluster, triple) %>% filter(triple == 3)
  # }

  p <- plot_umap_by_factor(a1$obs, "leiden", palette = "grayC")
  my_ggsave(
    "umap-facet-by-cluster",
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1 + 0.05 * length(unique(a1$obs$leiden)),
    width = 10, height = 6,
    units = "in", dpi = 300
  )

  p <- plot_relative_hex(
    x = a1$obs$UMAP1,
    y = a1$obs$UMAP2,
    group = a1$obs$case == "Case",
    bins = 91,
    # palette = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
    # palette = rev(RColorBrewer::brewer.pal(name = "PiYG", n = 11))
    # palette = rev(RColorBrewer::brewer.pal(name = "Spectral", n = 11))
    palette = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 11))
    # palette = "roma",
    # palette = "berlin",
    # direction = 1
  ) +
  labs(title = "Enrichment of cells from Cases")
  my_ggsave(
    glue("umap-case-enrichment"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1,
    width = 4.2,
    height = 3,
    units = "in",
    dpi = 300
  )

  source("R/functions/composition.R")
  print_status("do_masc()")
  my_obs$cluster <- my_obs$leiden
  res5_file <- glue("{out_dir}/composition.qs")
  if (!file.exists(res5_file)) {
    res5 <- do_masc(
      my_obs,
      form1 = "is_cluster ~ 1 + (1|donor)",
      form2 = "is_cluster ~ 1 + case + (1|donor)",
      mc.cores = 8
    )
    qsave(res5, res5_file)
  } else {
    res5 <- qread(res5_file)
  }

  source("R/functions/composition.R")
  print_status("do_masc()")
  my_obs$cluster <- my_obs$leiden
  res3_file <- glue("{out_dir}/composition-sex-case.qs")
  if (!file.exists(res3_file)) {
    res3 <- do_masc(
      my_obs,
      form1 = "is_cluster ~ 1 + sex + (1|donor)",
      form2 = "is_cluster ~ 1 + sex + case + (1|donor)",
      mc.cores = 8
    )
    qsave(res3, res3_file)
  } else {
    res3 <- qread(res3_file)
  }

  comp <- get_composition(d = my_obs, fill = case, group = donor, x = cluster)
  comp$freq <- comp$freq / 100
  if ("sex" %in% colnames(my_obs)) {
    comp <- left_join(comp, my_obs %>% select(donor, sex) %>% unique, by = "donor")
  } else {
    comp <- left_join(comp, my_obs %>% select(donor) %>% unique, by = "donor")
  }
  comp_file <- glue("{out_dir}/composition.tsv")
  fwrite(comp, comp_file, sep = "\t")
  print_status("done")

  res3_coef <- summarize_masc(res3)
  res3_coef$value <- str_replace(res3_coef$value, "case", "")
  res3_coef <- res3_coef %>%
    dplyr::group_by(value) %>%
    dplyr::mutate(
      OR      = exp(est),
      OR_2.5  = exp(est_low),
      OR_97.5 = exp(est_high),
      lrt_fdr = p.adjust(lrt_p, method = "fdr")
    ) %>% ungroup
  #
  tsv_file1 <- glue("{out_dir}/masc_1-case-complete.tsv")
  fwrite(res3_coef, tsv_file1, sep = "\t")

  #
  for (key in safe(unique(res3_coef$value))) {
    tsv_file2 <- glue("{out_dir}/masc_1-case-{safe(key)}.tsv")
    xx <- res3_coef %>% filter(value == key)
    xx %>%
      dplyr::select(cluster, OR, OR_2.5, OR_97.5, p, lrt_p, lrt_fdr) %>%
      dplyr::mutate(label = glue("{safe(key)} cluster {cluster} {signif(OR,2)}-fold (95% CI {signif(OR_2.5,2)} to {signif(OR_97.5,2)}, P = {str_replace(signif(lrt_p, 2), '-0', '-')})")) %>%
      dplyr::mutate_if(is.numeric, signif, 5) %>%
      dplyr::arrange(lrt_p) %>%
      write_tsv(tsv_file2)
  }
  #
  res3_coef$cluster <- factor(
    as.character(res3_coef$cluster),
    (
      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
    )$cluster
  )

  cluster_groups <- get_cluster_groups(analysis_name)
  #if (analysis_name == "n3_2") {
  #  cluster_groups[["Immature epithelial cells"]] <- c("8", "3", "14")
  #  cluster_groups[["Absorptive epithelial cells"]] <- c("1", "5", "6", "16", "18", "9", "12")
  #  cluster_groups[["Mature, absorptive epithelial cells"]] <- c("2", "11", "20")
  #  cluster_groups[["Secretory cells"]] <- c("4", "10", "17", "15", "19")
  #  cluster_groups[["Mesenchymal cells"]] <- c("7", "23", "21", "22", "13")
  #  stopifnot(sort(as.integer(unlist(cluster_groups))) == sort(unique(a1$obs$leiden)))
  #  cluster_groups <- split(rep(names(cluster_groups), lengths(cluster_groups)), unlist(cluster_groups))
  #  #
  #  a1$obs$cluster_group <- unname(unlist(cluster_groups[as.character(a1$obs$leiden)]))
  #  my_obs$cluster_group <- unname(unlist(cluster_groups[as.character(my_obs$leiden)]))
  #  res3_coef$cluster_group <- unname(unlist(cluster_groups[as.character(res3_coef$cluster)]))
  #  #
  #  cluster_group_levels <- (
  #    res3_coef %>% filter(value == "Case") %>% group_by(cluster_group) %>%
  #      summarize(min_p = min(lrt_p)) %>% arrange(min_p)
  #  )$cluster_group
  #  res3_coef$cluster_group <- factor(res3_coef$cluster_group, rev(cluster_group_levels))
  #  res3_coef$cluster <- factor(
  #    as.character(res3_coef$cluster),
  #    as.character((res3_coef %>% filter(value == "Case") %>% arrange(cluster_group, -lrt_p))$cluster)
  #  )
  #}
  saveRDS(levels(res3_coef$cluster), file.path(out_dir, "cluster_order.rds"))

  if (length(cluster_groups) > 0) {
    res3_coef$cluster_group <- naturalfactor(cluster_groups[as.character(res3_coef$cluster)])
  }
  n_clusters <- length(unique(my_obs$leiden))
  d_error <- res3_coef %>% filter(value == "Case") %>%
    mutate(
      est      = log2(exp(est)),
      est_low  = log2(exp(est_low)),
      est_high = log2(exp(est_high))
    )
  # if ("cluster_group" %in% colnames(res3_coef)) {
  #   d_error$cluster_group <- factor(d_error$cluster_group, rev(levels(d_error$cluster_group)))
  # }
  vlines <- seq(-8, 8, by = 1)
  vlines <- vlines[vlines > min(d_error$est_low)]
  vlines <- vlines[vlines < max(d_error$est_high)]
  vlines <- vlines[vlines != 0]
  p_error <- ggplot(d_error) +
    aes(x = est, y = cluster) +
    ggforestplot::geom_stripes() +
		geom_vline(
			xintercept = vlines,
			size = 0.3, color = "white"
		) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_point(aes(color = lrt_fdr < 0.05)) +
    geom_errorbarh(
      mapping = aes(xmin = est_low, xmax = est_high, color = lrt_fdr < 0.05),
      height = 0
    ) +
    geom_text(
      data = d_error %>% filter(lrt_fdr < 0.05),
      mapping = aes(
        x = ifelse(est > 0, -Inf, Inf),
        hjust = ifelse(est > 0, 0, 1),
        y = cluster,
        label = str_replace(sprintf(" %s ", signif(lrt_p, 1)), "-0", "-")
      ),
      size = 5, color = "grey30"
    ) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = seq(-8, 8, by = 1),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    # annotation_logticks(sides = "b", size = 0.3) +
    # expand_limits(y = c(0.5, max(as.integer(res3_coef$cluster)) + 0.5)) +
    # facet_grid(~ value) +#, scales = "free_x") +
    labs(title = "Case vs Control", y = NULL, x = "OR") +
    theme(plot.background = element_blank())
  if ("cluster_group" %in% colnames(res3_coef)) {
    p_error <- p_error +
      facet_grid(rows = vars(cluster_group), scales = "free_y", space = "free") +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
  }
  my_ggsave(
    "lanes_1-case",
    out_dir = out_dir,
    # type = c("pdf", "png"),
    type = "pdf",
    plot = p_error,
    scale = 1,
    width = 4,
    height = length(unique(res3_coef$cluster)) * 0.25 + 1,
    units = "in", dpi = 300
  )
  #
  source("R/functions/okabe-ito.R")
  my_obs$cluster <- factor(
    as.character(my_obs$leiden),
    (
      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
    )$cluster
  )
  if ("cluster_group" %in% colnames(res3_coef)) {
    my_obs$cluster_group <- factor(
      cluster_groups[as.character(my_obs$leiden)], levels(res3_coef$cluster_group)
    )
  }
  if ("cluster_group" %in% colnames(my_obs)) {
    my_obs$cluster <- factor(my_obs$cluster, levels(d_error$cluster))
    my_comp <- my_obs %>%
      dplyr::select(case, donor, cluster, cluster_group) %>%
      dplyr::count(donor, cluster, case, cluster_group)
    my_comp <- my_comp %>% dplyr::select(case, donor, cluster, cluster_group, n) %>%
      dplyr::group_by(donor) %>%
      dplyr::mutate(freq = 100 * n / sum(n))
    p_boxplot <- plot_composition_h(
      # d = my_obs, fill = case, group = donor, x = cluster,
      fill = case, group = donor, x = cluster,
      comp = my_comp,
      legend.position = "right"
    ) +
    labs(title = glue("Composition of each donor (n = {length(unique(my_obs$donor))})"))
    p_boxplot <- p_boxplot +
      facet_grid(rows = vars(cluster_group), scales = "free_y", space = "free") +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
  } else {
    p_boxplot <- plot_composition_h(
      d = my_obs, fill = case, group = donor, x = cluster,
      legend.position = "right"
    ) +
    labs(title = glue("Composition of each donor (n = {length(unique(my_obs$donor))})"))
  }
	my_ggsave(
		"compositionh",
		out_dir = out_dir,
		# type = c("pdf", "png"),
		type = "pdf",
		plot = p_boxplot,
		scale = 1, height = n_clusters * 0.5 + 1, width = 6, units = "in", dpi = 300
	)
  #
  p_both <- (
    p_error + 
      labs(title = "Case vs Control", y = NULL, x = "OR")
  ) + (
    p_boxplot +
      theme(
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
  p_both <- p_both + plot_layout(widths = c(1, 3))
	my_ggsave(
		"lanes_1-case-compositionh",
		out_dir = out_dir,
		# type = c("pdf", "png"),
		type = "pdf",
		plot = p_both,
		scale = 1, height = n_clusters * 0.5 + 1, width = 10, units = "in", dpi = 300
	)

  if ("sex" %in% colnames(my_obs)) {
    my_d <- left_join(
      x = get_composition(my_obs, case, donor, cluster),
      y = my_obs %>% select(donor, sex) %>% unique,
      by = "donor"
    )
    p <- ggplot(my_d) +
      geom_point(
        mapping = aes(y = cluster, x = freq, fill = case, group = case),
        position = position_quasirandom(groupOnX = FALSE, dodge.width = 0.9),
        shape = 21, size = 3, stroke = 0.2
      ) +
      facet_grid(cols = vars(sex)) +
      scale_fill_manual(
        name = NULL, values = pals::okabe(8)[1:8]
      ) +
      scale_x_log10(
        name = "Percent",
        #labels = scales::label_number()
        labels = function(x) signif(x, 3)
      ) +
      annotation_logticks(side = "b") +
      labs(x = NULL)
    my_ggsave(
      "case-sex",
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1, height = n_clusters * 0.5 + 1, width = 10, units = "in", dpi = 300
    )
  }

  # Percent of cells from each leiden
  cluster_colors <- mpn65
  names(cluster_colors) <- seq_along(cluster_colors)
  donor_cells <- a1$obs %>%
    group_by(donor) %>%
    count(leiden)
  donor_cells <- donor_cells %>% group_by(donor) %>% mutate(pct = 100 * n / sum(n))
  # donor_cells$donor <- pub_ids[donor_cells$donor]
  if (all(donor_cells$donor %in% names(pub_ids))) {
    donor_cells$pub <- pub_ids[donor_cells$donor]
    donor_cells$pub <- naturalfactor(donor_cells$pub)
    donor_cells$pub <- factor(donor_cells$pub, rev(levels(donor_cells$pub)))
  } else {
    donor_cells$pub <- factor(donor_cells$donor)
  }
  donor_cells$leiden <- naturalfactor(donor_cells$leiden)
  p_percents <- ggplot(donor_cells) +
    aes(x = pct, y = pub, fill = leiden) +
    geom_colh() +
    scale_fill_manual(name = NULL, values = cluster_colors) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      legend.position = "none"
    ) +
    labs(x = "Percent", y = NULL)
  donor_cells_pub <- donor_cells %>% group_by(pub) %>% summarize(n = sum(n))
  p_counts <- ggplot(donor_cells) +
    geom_colh(
      data = donor_cells,
      mapping = aes(x = n, y = pub, fill = leiden)
    ) +
    geom_text(
      data = donor_cells_pub,
      mapping = aes(x = n, y = pub, label = comma(n, accuracy = 1)),
      hjust = 0, nudge_x = 100
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.35)), labels = label_number_si()
    ) +
    scale_fill_manual(name = NULL, values = mpn65) +
    labs(x = "Cells", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  my_ggsave(
    "bars-donor_cells",
    out_dir = out_dir,
    type = "pdf",
    plot = p_percents + p_counts + plot_annotation(title = "Composition of each donor's cells"),
    scale = 0.75,
    width = 11, height = 10,
    units = "in", dpi = 300
  )
  donor_cells_cluster <- donor_cells %>% group_by(leiden) %>% summarize(n = sum(n))
  n_pub <- length(unique(donor_cells$pub))
  p_counts2 <- ggplot(donor_cells) +
    # geom_colh(
    #   data = donor_cells_cluster,
    #   # mapping = aes(x = n, y = leiden, fill = pub)
    #   mapping = aes(x = -n, y = rev(leiden), fill = leiden)
    # ) +
    geom_text(
      data = donor_cells_cluster,
      # mapping = aes(x = -n, y = rev(leiden), label = comma(n, accuracy = 1)),
      # hjust = 1, nudge_x = -100
      mapping = aes(x = 0, y = rev(leiden), label = comma(n, accuracy = 1)),
      hjust = 1, size = 5
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(1, 0)), labels = label_number_si()
    ) +
    # scale_fill_manual(
    #   name = NULL,
    #   values = scico::scico(direction = -1, n = n_pub + 1, palette = "grayC")[1:n_pub]
    # ) +
    scale_fill_manual(
      name = NULL,
      values = cluster_colors
    ) +
    labs(x = NULL, y = NULL) +
    theme_void()
    # theme(
    #   legend.position = "none",
    #   # legend.position = c(1, 0),
    #   # legend.just = c(1, 0),
    #   # legend.background = element_blank()
    #   axis.text.y = element_blank(),
    #   axis.ticks.y = element_blank(),
    #   axis.text.x = element_blank(),
    #   axis.ticks.x = element_blank()
    # )
  donor_cells$leiden <- factor(
    x = as.character(donor_cells$leiden),
    levels = rev(unique(naturalsort(donor_cells$leiden)))
  )
  donor_cells$pub <- factor(
    x = as.character(donor_cells$pub),
    levels = unique(naturalsort(donor_cells$pub))
  )
  p_heatmap <- ggplot(donor_cells) +
  geom_tile(
    aes_string(x = "pub", y = "leiden", fill = "pct")
  ) +
  scale_x_discrete(position = "b", name = NULL, expand = c(0, 0)) +
  scale_y_discrete(position = "r", name = NULL, expand = c(0, 0)) +
  scale_fill_gradientn(
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    colors = scico::scico(n = 21, pal = "batlow"),
    trans = "log10",
    # breaks = log_breaks,
    name = "Percent",
    guide = guide_colorbar(barwidth = 20)
  ) +
  # scale_fill_gradientn(
  #   # colors = RColorBrewer::brewer.pal(name = "Greys", n = 9),
  #   colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
  #   guide = guide_colorbar(barwidth = 15),
  #   breaks = pretty_breaks(5)
  # ) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0,0,0,0, "lines"),
    legend.margin = margin(0, 0, 0, 0, "lines"),
    legend.background = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )
  p_clusterbar <- ggplot(donor_cells) +
  geom_tile(
    aes(x = 1, y = leiden, fill = leiden)
  ) +
  scale_x_discrete(position = "b", name = NULL, expand = c(0, 0)) +
  scale_y_discrete(
    position = "r", name = NULL, expand = c(0, 0),
    labels = function(x) num_to_name[x]
  ) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  theme(axis.ticks.y = element_blank())
  p <- (
    p_heatmap
  ) + (
    p_clusterbar + theme(plot.margin = margin(r = 0))
  ) + plot_layout(widths = c(length(unique(donor_cells$pub)), 0.5))
  fig_width  <- 2 + length(unique(donor_cells$pub)) * 0.4
  fig_height <- length(unique(donor_cells$leiden)) * 0.3 + 2
  my_ggsave(
    slug = glue("heatmap-donor_cells"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )
  layout <- "
A##
B##
CDE
"
  p_percents <- ggplot(donor_cells) +
    aes(x = pct, y = pub, fill = leiden) +
    geom_colh() +
    scale_fill_manual(name = NULL, values = cluster_colors) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    labs(x = "Percent", y = NULL) +
    coord_flip()
  p_counts <- ggplot(donor_cells) +
    geom_colh(
      data = donor_cells,
      mapping = aes(x = n, y = pub, fill = leiden)
    ) +
    geom_text(
      data = donor_cells_pub,
      mapping = aes(x = n, y = pub, label = comma(n, accuracy = 1)),
      angle = 90, hjust = 0, nudge_x = 100
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.45)), labels = label_number_si()
    ) +
    scale_fill_manual(name = NULL, values = cluster_colors) +
    labs(x = "Cells", y = NULL) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    coord_flip()
  p <- p_counts + p_percents + p_heatmap + p_counts2 + p_clusterbar + plot_layout(
    design = layout,
    widths = c(length(unique(donor_cells$pub)), 2.2, 0.5),
    heights = c(1, 1, 3),
  )
  fig_width  <- 2 + length(unique(donor_cells$pub)) * 0.4
  fig_height <- length(unique(donor_cells$leiden)) * 0.4 + 3
  my_ggsave(
    slug = glue("heatmap-donor_cells_full"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )

  a1$de <- presto::wilcoxauc(a1$log2cpm, a1$obs$leiden)
  #
  de_file <- glue("results/a20/{analysis_name}/figures/pseudobulk_de_ova.tsv.gz")
  stopifnot(file.exists(de_file))
  de <- fread(de_file)
  de$group <- str_split_fixed(de$contrast, " ", 2)[,1]
  a1$de <- inner_join(
    x = a1$de %>% select(group, feature, auc, pct_in, pct_out, avgExpr),
    y = de %>% select(-auc),
    by = c("group", "feature" = "ensembl_id")
  )
  if ("cluster_group" %in% colnames(res3_coef)) {
    #
    a1$de$cluster_group <- cluster_groups[as.character(a1$de$group)]
    a1$de$cluster_group <- naturalfactor(a1$de$cluster_group)
    a1$de$group <- factor(as.character(a1$de$group), as.character(levels(res3_coef$cluster)))
  }

  n_donors <- length(unique(a1$obs$donor))
  umap_alpha <- 0.4
  umap_pointsize <- 1
  if (nrow(a1$obs) < 5e4) {
    umap_pointsize <- 1.5
    umap_alpha <- 0.6
  } 
  if (nrow(a1$obs) < 3e4) {
    umap_pointsize <- 2
    umap_alpha <- 0.6
  } 
  p_umap <- plot_scattermore(
    x = a1$obs$UMAP1,
    y = a1$obs$UMAP2,
    group = a1$obs$leiden,
    group_colors = cluster_colors,
    pixels = 1000,
    pointsize = umap_pointsize,
    alpha = umap_alpha
  ) +
  labs(
    title = glue(
      "{length(unique(a1$obs$leiden))} clusters of {comma(nrow(a1$obs))} cells from {n_donors} donor{ifelse(n_donors > 1, 's', '')}"
    )
  )
  my_ggsave(
    "umap-clusters",
    out_dir = out_dir,
    plot = p_umap,
    type = "pdf",
    scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
  )

  # Heatmap of best markers for each cluster
  ########################################################################
  cluster_colors <- mpn65
  names(cluster_colors) <- as.character(seq_along(cluster_colors))
  n_top <- 3
  # if (analysis_name == "n3_2") {
  #   n_top <- 4
  # }
  # Top markers for each cluster
  a1$de <- a1$de %>%
    mutate(symbol = ensembl_to_symbol[feature])
  a1$de_top <- a1$de %>%
    dplyr::filter(!group %in% exclude_clusters) %>%
    dplyr::group_by(group) %>%
    dplyr::filter(pct_in > 1) %>%
    # dplyr::filter(group %in% c(12, 14, 16, 17, 19, 25, 29)) %>%
    dplyr::mutate(
      rank1 = rank(100 * auc - pct_out),
      rank2 = rank(logFC * (pct_in - pct_out)),
    ) %>%
    dplyr::top_n(n = n_top, wt = (rank1 + rank2))
  #
  these_genes <- unique(a1$de_top$feature)
  #
  if (analysis_name %in% cluster_names$analysis) {
    selected_genes <- (
      cluster_names %>%
      filter(analysis == analysis_name)
    )$markers
    selected_genes <- unique(unlist(str_split(selected_genes, ", ")))
    selected_genes <- selected_genes[selected_genes %in% a1$de$symbol]
    if (length(selected_genes) > 1) {
      these_genes <- unique(
        (
          a1$de %>% select(symbol, feature) %>% filter(symbol %in% selected_genes) %>% unique
        )$feature
      )
    }
  }
  these_genes <- these_genes[!duplicated(ensembl_to_symbol[these_genes])]

  if (analysis_name == "n3_2") {
    #which(ensembl_to_symbol == "SELENOP")
    # which(ensembl_to_symbol == "CX3CL1")
    these_genes <- c(these_genes, "ENSG00000250722", "ENSG00000006210")
  }

  #
  # Subset the table
  # my_value <- "auc"
  # my_value <- "z"
  for (my_value in c("auc", "z")) {

    x <- a1$de %>% dplyr::filter(feature %in% these_genes)
    x <- x %>% group_by(feature) %>% mutate(z = scale(avgExpr)) %>% ungroup
    x <- as.data.frame(dcast.data.table(
      data = as.data.table(x), formula = symbol ~ group, value.var = my_value
    ))
    rownames(x) <- x$symbol
    x$symbol <- NULL
    if ("cluster_group" %in% colnames(a1$de)) {
      x <- x[,levels(a1$de$group)]
    }
    # Order the heatmap
    set.seed(1)
    # xo <- seriation::seriate(as.matrix(x), method = "BEA_TSP")
    my_cols <- levels(res3_coef$cluster)
    x <- x[,my_cols]
    xd <- a1$de %>% dplyr::filter(feature %in% these_genes)
    xd <- xd %>% group_by(feature) %>% mutate(z = scale(avgExpr)) %>% ungroup
    # xd$symbol <- factor(xd$symbol, rownames(x)[xo[[1]]])
    xd$symbol <- factor(as.character(xd$symbol), names(sort(apply(x, 1, which.max))) %>% rev)
    # xd$group <- factor(xd$group, colnames(x)[xo[[2]]])
    xd$group <- factor(as.character(xd$group), levels(res3_coef$cluster))
    #
    if (analysis_name %in% cluster_names$analysis) {
      num_to_name <- unlist(with(
        cluster_names %>% filter(analysis == analysis_name),
        split(name, as.character(number))
      ))
    } else {
      num_to_name <- levels(xd$group)
      names(num_to_name) <- as.character(levels(xd$group))
    }
    # Plot
    my_limits <- c(0, 1)
    my_name <- "AUC\n\n"
    col_fun <- circlize::colorRamp2(
      seq(my_limits[1], my_limits[2], length.out = 11),
      # rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
      scico::scico(palette = "bam", n = 11, direction = -1)
    )
    if (my_value == "z") {
      my_limits <- range(xd[[my_value]])
      my_name <- "Scaled expression across cell clusters\n\n"
      col_fun <- circlize::colorRamp2(
        seq(-max(abs(my_limits)), max(abs(my_limits)), length.out = 11),
        # rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
        scico::scico(palette = "bam", n = 11, direction = -1)
      )
    }
    color_values <- seq(my_limits[1], my_limits[2], length.out = 11)
    if (length(cluster_groups) > 0) {
      xd$cluster_group <- cluster_groups[as.character(xd$group)]
      xd <- xd %>% group_by(feature) %>% mutate(gene_group = cluster_group[which.max(auc)])
    }
    p_heatmap <- ggplot(xd) +
    geom_tile(
      aes_string(x = "symbol", y = "group", fill = my_value)
    ) +
    scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
    scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
    scale_fill_gradientn(
      colors = col_fun(color_values),
      limits = my_limits,
      name = my_name,
      guide = guide_colorbar(barwidth = 15), breaks = pretty_breaks(5)
    ) +
    theme(
      axis.text.x = element_text(size = 12, face = "italic", angle = 60, hjust = 0),
      legend.position = "bottom",
      plot.margin = margin(0,0,0,0, "lines"),
      legend.margin = margin(0, 0, 0, 0, "lines"),
      legend.background = element_blank()
    )
    #
    p_clusterbar <- ggplot(xd) +
    geom_tile(
      aes(x = 1, y = group, fill = group)
    ) +
    scale_x_discrete(position = "b", name = NULL, expand = c(0, 0)) +
    scale_y_discrete(
      position = "r", name = NULL, expand = c(0, 0),
      labels = function(x) num_to_name[x]
    ) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = 18))
    if (length(cluster_groups) > 0)  {
      p_clusterbar <- p_clusterbar +
        facet_grid(rows = vars(cluster_group), scales = "free", space = "free") +
        theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
      p_heatmap <- p_heatmap +
        facet_grid(cols = vars(gene_group), rows = vars(cluster_group), scales = "free", space = "free") +
        theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
    }
    p <- (
      p_clusterbar + theme(plot.margin = margin(r = 0))
    ) + (
      p_heatmap + theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(l = 0)
      )
    ) + plot_layout(widths = c(1, length(these_genes)))
    fig_width  <- 2 + length(these_genes) * 0.3
    fig_height <- length(unique(a1$de_top$group)) * 0.4 + 2
    my_ggsave(
      slug = glue("heatmap-top-markers_h-{my_value}"),
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
    )
    #
    my_height <- n_clusters * 0.6 + 1
    if (n_clusters < 9) {
      my_height <- n_clusters * 0.8 + 1
    }
    #
    my_widths <- c(4, 0.15, length(these_genes) * 0.2, 1.5, 3)
    if (analysis_name == "a12_4_4_t4_cd8_1_2") {
      my_widths <- c(4, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "a12_4_4_t4_cd4_2_2") {
      my_widths <- c(3, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "a12_4_4_m3_2") {
      my_widths <- c(2.5, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "a12_4_4_b5_1_3") {
      my_widths <- c(3, 0.15, length(these_genes) * 0.2, 1.5, 3)
      my_height <- n_clusters * 0.5 + 1
    }
    if (analysis_name == "n3_2") {
      my_widths <- c(6, 0.15, length(these_genes) * 0.2, 1.5, 3)
      my_height <- n_clusters * 0.5 + 1
    }
    if (analysis_name == "blood2") {
      my_widths <- c(4, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "blood2_tcell5") {
      my_widths <- c(3, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "blood2_myeloid5") {
      my_widths <- c(4, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "blood2_bcell5") {
      my_widths <- c(1.5, 0.15, length(these_genes) * 0.2, 1.5, 3)
      my_height <- n_clusters * 1 + 1
    }
    if (analysis_name == "blood2_tcell5_cd4_5") {
      my_widths <- c(2, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    if (analysis_name == "blood2_tcell5_cd8_5") {
      my_widths <- c(4, 0.15, length(these_genes) * 0.2, 1.5, 3)
    }
    #
    p <- (
      p_umap
    ) + (
      p_clusterbar + theme(
        axis.text.y = element_text(hjust = 0),
        plot.margin = margin(r = 0, b = 0),
        legend.margin = margin(t = 0, b = 0)
      )
    ) + (
      p_heatmap + theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(l = 0)
      )
    ) + (
      p_error + labs(title = "Case vs Control", y = NULL, x = "OR")
    ) + (
      p_boxplot +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom"
        )
    ) + plot_layout(widths = my_widths)
    my_width <- 13 + length(these_genes) * 0.4
    if (analysis_name %in% c("n3_2", "blood2_tcell5_cd8_5")) {
      my_width <- 15 + length(these_genes) * 0.4
    }
    my_ggsave(
      glue("full-{my_value}"),
      out_dir = out_dir,
      # type = c("pdf", "png"),
      type = "pdf",
      plot = p,
      scale = 0.8,
      height = my_height,
      width = my_width,
      units = "in", dpi = 300,
      limitsize = FALSE
    )

  }

  umap_n_bins <- 37
  if (nrow(a1$obs) > 5e4) {
    umap_n_bins <- 57
  }
  if (nrow(a1$obs) > 1e5) {
    umap_n_bins <- 97
  }

  if (analysis_name == "luoma_cd45_a5_tcell2_cd8_3") {
    my_ids <- names(
      ensembl_to_symbol[ensembl_to_symbol %in% c("ITGAE", "ITGB2", "CX3CR1", "GZMK")]
    )
    for (my_id in my_ids) {
      my_gene <- ensembl_to_symbol[my_id]
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

  if (analysis_name == "luoma_cd45_a5_tcell2_cd4_3") {
    my_ids <- names(
      ensembl_to_symbol[ensembl_to_symbol %in% c("CCL5", "TCF7", "FOXP3", "PDCD1")]
    )
    for (my_id in my_ids) {
      my_gene <- ensembl_to_symbol[my_id]
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

  for (i in seq(1:nrow(a1$de_top))) {
    my_id <- a1$de_top$feature[i]
    my_cluster <- a1$de_top$group[i]
    p <- plot_hexgene(
      x = a1$obs$UMAP1,
      y = a1$obs$UMAP2,
      z = as.numeric(a1$log2cpm[my_id,]),
      bins = umap_n_bins,
      palette = "oslo"
    ) +
    labs(title = ensembl_to_symbol[my_id]) +
    theme(
      legend.position = "none",
      strip.text = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(face = "italic")
    )
    my_ggsave(
      glue("umap-cluster-{str_pad(my_cluster, 2, pad = '0')}-{safe(ensembl_to_symbol[my_id])}"),
      out_dir = file.path(out_dir, "umap/slim"),
      type = "pdf",
      plot = p,
      scale = 0.8,
      width = 3.5,
      height = 3,
      units = "in",
      dpi = 300
    )
  }

  if (analysis_name == "a12_4_4_t4_cd8_1_2") {
    cd8_genes <- list(
      "Naive/Effector" = c("IL7R", "TCF7", "GZMA", "GZMK", "GZMB"),
      "Circulation" = c("KLF2", "SELL", "S1PR1", "S1PR5"),
      "Exhaustion" = c("HAVCR2", "KIR2DL4", "ENTPD1", "TIGIT"),
      "Homing" = c("ITGAE", "KLRG1", "ITGB2", "CXCR3", "CX3CR1"),
      "TF" = c("ZNF683", "EOMES", "IKZF2", "ID3")
    )
    for (my_name in unlist(cd8_genes)) {
      my_id <- names(which(ensembl_to_symbol == my_name))
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-cluster--{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

  if (analysis_name == "a12_4_4_t4_cd4_2_2") {
    cd4_genes = list(
      "Naive/Tcm" = c("ITGA1", "IL23R", "TCF7", "CCL5"),
      "Effector" = c("IL17A", "IFNG", "TNFRSF18", "CXCR3"),
      "Tfh" = c("BCL6", "CXCR5", "CXCL13", "PDCD1"),
      "Treg" = c("FOXP3", "TNFRSF4", "SELL", "LAG3")
    )
    for (my_name in unlist(cd4_genes)) {
      my_id <- names(which(ensembl_to_symbol == my_name))
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-cluster--{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

  if (analysis_name == "blood2_tcell5_cd8_5") {
    blood_cd8_genes <- list(
      "Lineage" = c("TRAC", "TRDC", "NCAM1", "CD8B"),
      "Innate" = c("TYROBP", "FCER1G", "KLRC3", "TRAV1-2"),
      "Intravascular" = c("CX3CR1", "FGFBP2", "GZMH", "S1PR5"),
      "Effector" = c("GZMK", "HLA-DRA", "MKI67", "CD69"),
      "Naive" = c("SELL", "CCR7", "IL7R", "LTB")
    )
    for (my_name in unlist(blood_cd8_genes)) {
      my_id <- names(which(ensembl_to_symbol == my_name))
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-cluster--{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

  if (analysis_name == "n3_2") {
    my_names <- c("CHGA", "EPCAM", "CCL21", "SH2D6")
    for (my_name in my_names) {
      my_id <- names(which(ensembl_to_symbol == my_name))
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-cluster--{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

}

# }}}

# Case vs Control DE tables {{{
# Figure 3
########################################################################

plot_de_summary <- function(
  pb_donor_meta, pb_donor, de_donor, de, ensembl_ids, out_file, cluster_order = NA, title = NA
) {
  d <- as_tibble(reshape2::melt(as.matrix(pb_donor)[ensembl_ids,]))
  colnames(d) <- c("ensembl_id", "donor", "value")
  d$GeneName <- factor(
    as.character(unname(ensembl_to_symbol[as.character(d$ensembl_id)])),
    rev(as.character(unname(ensembl_to_symbol[as.character(ensembl_ids)])))
  )
  d$Gene <- as.integer(d$GeneName)
  d <- left_join(d, pb_donor_meta, by = "donor")
  #
  d_error <- left_join(d, de_donor, by = c("ensembl_id"))
  print(head(d_error))
  p_error <- ggplot(d_error) +
    # annotate(
    #   geom = "rect",
    #   xmin = -Inf,
    #   xmax = Inf,
    #   ymin = seq(from = 1, to = length(ensembl_ids), by = 2) - 0.5,
    #   ymax = seq(from = 1, to = length(ensembl_ids), by = 2) + 0.5,
    #   alpha = 0.2
    # ) +
    # aes(x = 2^logFC, xmin = 2^CI.L, xmax = 2^CI.R, y = Gene.x) +
    aes(x = logFC, xmin = CI.L, xmax = CI.R, y = GeneName) +
    ggforestplot::geom_stripes(aes(xmin = -Inf, xmax = Inf)) +
    scale_y_discrete(expand = c(0, 0.5)) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_point() +
    geom_errorbarh(height = 0) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(3),
      # labels = function(x) 2^x
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    # annotation_logticks(side = "b") +
    # scale_y_continuous(
    #   expand = c(0, 0)
    #   # breaks = seq(1, max(d_error$Gene.x)),
    #   # labels = levels(d_error$GeneName)
    # ) +
    # expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
    labs(x = "FC", y = NULL) +
    theme(
      axis.text.y = element_text(face = "italic"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  p_donor <- ggplot(d_error) +
      # annotate(
      #   geom = "rect",
      #   xmin = -Inf,
      #   xmax = Inf,
      #   ymin = seq(from = 1, to = length(ensembl_ids), by = 2) - 0.5,
      #   ymax = seq(from = 1, to = length(ensembl_ids), by = 2) + 0.5,
      #   alpha = 0.2
      # ) +
      # geom_vline(xintercept = 0, size = 0.3) +
      ggforestplot::geom_stripes(aes(y = GeneName)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      geom_boxplot(
        data = d_error, mapping = aes(x = value, y = GeneName, fill = case), size = 0.3,
        outlier.size = 0.5
      ) +
      # geom_point(
      #   data = d_error,
      #   mapping = aes(x = value, y = GeneName, group = case, fill = case),
      #   position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE),
      #   shape = 21, size = 2.5, stroke = 0.2
      # ) +
      scale_fill_manual(values = pals::okabe(8)) +
      scale_x_continuous(
        breaks = scales::pretty_breaks(3)
      ) +
      # scale_x_continuous(
      #   breaks = seq(-10, 10, by = 1),
      #   labels = function(x) fractional::fractional(2 ^ x)
      # ) +
      # scale_y_continuous(
      #   expand = c(0, 0)
      #   # breaks = seq(1, max(d_error$Gene.x)),
      #   # labels = levels(d_error$GeneName)
      # ) +
      guides(fill = guide_legend(reverse = TRUE, title = NULL, override.aes = list(size = 5))) +
      # annotation_logticks(sides = "b", size = 0.3) +
      # expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
      labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
      theme(
        axis.text.y = element_text(face = "italic"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)
      )
  if ("cat" %in% colnames(d_error)) {
    p_error <- p_error + facet_grid(rows = vars(cat), scales = "free_y", space = "free_y") +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
    p_donor <- p_donor + facet_grid(rows = vars(cat), scales = "free_y", space = "free_y") +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
  }
  # my_ggsave(
  #   "de_donor",
  #   out_dir = out_dir,
  #   type = "pdf",
  #   plot = p_error + (p_donor + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())),
  #   scale = 1,
  #   width = 5,
  #   height = length(ensembl_ids) * 0.3 + 1,
  #   units = "in",
  #   dpi = 300
  # )
  cluster_colors <- mpn65
  names(cluster_colors) <- as.character(seq_along(cluster_colors))
  x <- as_tibble(de %>% filter(ensembl_id %in% ensembl_ids))
  if (!all(ensembl_ids %in% x$ensembl_id)) {
    missing_ids <- setdiff(ensembl_ids, x$ensembl_id)
    stop(sprintf("de at cluster level is missing some genes: %s", missing_ids))
  }
  x$ensembl_id = factor(x$ensembl_id, rev(ensembl_ids))
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(x$logFC)), max(abs(x$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
    # colorRampPalette(c("black", "white", pals::okabe()[2]))(11)
  )
  # color_values <- scales::rescale(seq(min(x$logFC), max(x$logFC), length.out = 11))
  color_values <- seq(min(x$logFC), max(x$logFC), length.out = 11)
  if (!any(is.na(cluster_order))) {
    stopifnot(length(cluster_order) == length(unique(x$cluster)))
    x$cluster <- factor(x$cluster, cluster_order)
  } else {
    # x$cluster <- factor(x$cluster, (naturalsort(unique(x$cluster))))
    x_mat <- as.matrix(
      x %>%
        select(ensembl_id, cluster, logFC) %>%
        pivot_wider(names_from = "cluster", values_from = "logFC") %>%
        select(-ensembl_id)
    )
    set.seed(1)
    n_neighbors <- ifelse(ncol(x_mat) >= 15, 15, ncol(x_mat) - 1)
    order_umap <- order(
      uwot::umap(t(x_mat), n_components = 1, n_neighbors = n_neighbors)
    )
    clusters_inorder <- colnames(x_mat)[order_umap]
    x$cluster <- factor(as.character(x$cluster), clusters_inorder)
  }
  # p_mat_barheight <- max(10, length(unique(x$cluster)))
  p_mat_barheight <- 7
  if (length(unique(x$ensembl_id)) > 10) {
    p_mat_barheight <- 10
  }
  if (length(unique(x$ensembl_id)) > 20) {
    p_mat_barheight <- 15
  }
  p_mat <- ggplot(x) +
  geom_tile(
    aes(x = cluster, y = ensembl_id, fill = logFC)
  ) +
  geom_point(
    aes(x = cluster, y = ensembl_id, alpha = 1 * adj.P.Val < 0.05),
    color = "white"
  ) +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
  scale_y_discrete(
    position = "l", name = NULL, expand = c(0, 0),
    labels = function(x) ensembl_to_symbol[x]
  ) +
  scale_fill_gradientn(
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))[6:11],
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    # values = color_values,
    colors = col_fun(color_values),
    # limits = c(0, 1),
    # name = bquote("log"[2]~"FC"),
    name = "FC",
    guide = guide_colorbar(barheight = p_mat_barheight),
    # breaks = pretty_breaks(5),
    # breaks = c(1 / 2^(1:8), 2^(1:8)),
    breaks = -8:8,
    labels = function(x) fractional::fractional(2^x)
  ) +
  theme(
    axis.text.y = element_text(size = 14, face = "italic"),
    legend.position = "right"
  )
  de$cluster <- factor(de$cluster, (naturalsort(unique(de$cluster))))
  if ("cluster_group" %in% colnames(x)) {
    p_bar <- ggplot(x %>% select(cluster, cluster_group) %>% unique)
  } else {
    p_bar <- ggplot(x %>% select(cluster) %>% unique)
  }
  p_bar <- p_bar + geom_tile(
    aes(y = 1, x = cluster, fill = cluster)
  ) +
  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
  scale_fill_manual(values = cluster_colors, guide = "none")
  if ("cat" %in% colnames(x) && "cluster_group" %in% colnames(x)) {
    p_bar <- p_bar + facet_grid(cols = vars(cluster_group), scales = "free", space = "free") +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
    p_mat <- p_mat + facet_grid(
      rows = vars(cat),
      cols = vars(cluster_group),
      scales = "free", space = "free"
    ) +
    theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
  }
  else if ("cat" %in% colnames(x)) {
    p_mat <- p_mat + facet_grid(rows = vars(cat), scales = "free_y", space = "free_y") +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
  }
  # p <- (
  #     plot_spacer() + p_bar + theme(plot.margin = margin(l = 2))
  # ) / (
  #   (
  #     p_donor + theme(legend.position = "none")
  #   ) + (
  #     p_mat + theme(
  #       axis.text.y = element_blank(),
  #       axis.text.x = element_blank(),
  #       axis.ticks.x = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       axis.title.y = element_blank(),
  #       plot.margin = margin(l = 0)
  #     )
  #   )
  # ) + plot_layout(heights = c(1, length(ensembl_ids)))
  layout <- "
##C
ABD
"
  p_donor <- p_donor + theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
  p_mat <- p_mat + theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(l = 0)
  )
  p_mat_width <- 2.8
  if (length(unique(de$cluster)) > 15) {
    p_mat_width <- 6
  }
  p <- p_error + p_donor + p_bar + p_mat + plot_layout(
    design = layout,
    widths = c(1, 1, p_mat_width),
    heights = c(1, length(ensembl_ids))
  ) + plot_annotation(title = title)
  # fig_height <- length(ensembl_ids) * 0.3 + 2
  fig_height <- length(ensembl_ids) * 0.275 + 2
  fig_width <- length(unique(de$cluster)) * 0.3 + 5
  my_ggsave(
    slug = str_remove(basename(out_file), "\\.[^.]+$"),
    out_dir = dirname(out_file),
    type = "pdf",
    plot = p,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )
}

plot_de_summary_h <- function(
  pb_donor_meta, pb_donor, de_donor, de, ensembl_ids, out_file, cluster_order = NA, title = NA
) {
  d <- as_tibble(reshape2::melt(as.matrix(pb_donor)[ensembl_ids,]))
  colnames(d) <- c("ensembl_id", "donor", "value")
  d$GeneName <- factor(
    as.character(unname(ensembl_to_symbol[as.character(d$ensembl_id)])),
    rev(as.character(unname(ensembl_to_symbol[as.character(ensembl_ids)])))
  )
  d$Gene <- as.integer(d$GeneName)
  d <- left_join(d, pb_donor_meta, by = "donor")
  #
  cluster_colors <- mpn65
  names(cluster_colors) <- as.character(seq_along(cluster_colors))
  x <- de %>% filter(ensembl_id %in% ensembl_ids)
  x$ensembl_id = factor(x$ensembl_id, ensembl_ids)
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(x$logFC)), max(abs(x$logFC)), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  # color_values <- scales::rescale(seq(min(x$logFC), max(x$logFC), length.out = 11))
  color_values <- seq(min(x$logFC), max(x$logFC), length.out = 11)
  if (!any(is.na(cluster_order))) {
    x$cluster <- factor(x$cluster, cluster_order)
  } else {
    x$cluster <- factor(x$cluster, (naturalsort(unique(x$cluster))))
  }
  p_mat_barwidth <- 5
  if (length(unique(x$ensembl_id)) > 10) {
    p_mat_barwidth <- 15
  }
  if (length(unique(x$ensembl_id)) > 20) {
    p_mat_barwidth <- 30
  }
  p_mat <- ggplot(x) +
  geom_tile(
    aes(y = cluster, x = ensembl_id, fill = logFC)
  ) +
  geom_point(
    aes(y = cluster, x = ensembl_id, alpha = 1 * adj.P.Val < 0.05),
    color = "white"
  ) +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
  scale_x_discrete(
    position = "t", name = NULL, expand = c(0, 0),
    labels = function(x) ensembl_to_symbol[x]
  ) +
  scale_fill_gradientn(
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))[6:11],
    # colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    # values = color_values,
    colors = col_fun(color_values),
    # limits = c(0, 1),
    # name = bquote("log"[2]~"FC"),
    name = "FC",
    guide = guide_colorbar(barwidth = p_mat_barwidth),
    # breaks = pretty_breaks(5),
    # breaks = c(1 / 2^(1:8), 2^(1:8)),
    breaks = -8:8,
    labels = function(x) fractional::fractional(2^x)
  ) +
  theme(
    axis.text.x = element_text(size = 14, face = "italic", angle = 60, hjust = 0),
    legend.position = "bottom"
  )
  if ("cat" %in% colnames(x)) {
    p_mat <- p_mat + facet_grid(~ cat, scales = "free_x", space = "free_x", switch = "x") +
      theme(panel.spacing = unit(0.5, "lines"))
  }
  de$cluster <- factor(de$cluster, (naturalsort(unique(de$cluster))))
  # p_bar <- ggplot(de %>% count(cluster)) +
  p_bar <- ggplot(x %>% select(cluster) %>% unique) +
  geom_tile(
    aes(x = 1, y = cluster, fill = cluster)
  ) +
  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
  scale_fill_manual(values = cluster_colors, guide = "none")
  layout <- "AB"
  p_mat <- p_mat + theme(
    axis.text.y = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(l = 0)
  )
  p <- p_bar + p_mat + plot_layout(
    design = layout,
    widths = c(1, length(unique(x$ensembl_id))),
  ) + plot_annotation(title = title)
  fig_width <- length(ensembl_ids) * 0.3 + 1
  fig_height <- length(unique(de$cluster)) * 0.25 + 2
  my_ggsave(
    slug = str_remove(basename(out_file), "\\.[^.]+$"),
    out_dir = dirname(out_file),
    type = "pdf",
    plot = p,
    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
  )
}


# analysis_name <- "a12_4_4_t4_cd8_1_2"
# analysis_name <- "a12_4_4_t4_cd4_2_2"
# # analysis_name <- "a12_4_4_m3_2"
# # analysis_name <- "a12_4_4_b5_1_3"
# # analysis_name <- "n3_2"

for (analysis_name in analyses) {

  params <- list(
    min_cells_in_cluster = 50,
    min_percent_of_cells_with_gene = 5
  )
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
  if (analysis_name == "luoma_cd3_a5") {
    a1$obs$leiden <- a1$obs$leiden0.933
  }
  a1$obs$cluster <- a1$obs$leiden
	#
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	#
	sample_info <- janitor::clean_names(read_excel(
			path = "data/luoma_villani_combined_colitis_samples2.xlsx",
			sheet = 3
	))
	donor_to_sex <- unlist(split(sample_info$sex, sample_info$donor))
  #
  a1$obs$sex <- NA
  a1$obs$sex <- donor_to_sex[a1$obs$donor]
  #
  # a1$obs %>% select(donor, sex) %>% unique
  #
  if ("drug" %in% colnames(a1$obs)) {
    drug_levels <- c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4")
    a1$obs$drug <- factor(a1$obs$drug, drug_levels)
  } else {
    a1$obs$drug <- NA
  }
  ## Skip cell clusters that have too few cells.
  #exclude_clusters <- (
  #  a1$obs %>% count(cluster) %>% filter(n < params[["min_cells_in_cluster"]])
  #)$leiden
  ##
  #a1$obs <- a1$obs[!a1$obs$cluster %in% exclude_clusters,]
  #
  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))
  a1$obs$sex[a1$obs$donor == "SIC_33"] <- "M"
  gene_pct <- 100 * rowSums(a1$counts > 0) / ncol(a1$counts)
  cluster_colors <- head(mpn65, length(unique(a1$obs$leiden)))
  names(cluster_colors) <- as.character(seq(length(cluster_colors)))
  a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))

  # Pseudobulk at the donor level
  ########################################################################
  y <- with(a1$obs, model.matrix(~ 0 + factor(donor)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/") # means
  pb_donor <- as(a1$counts %*% y, "dgCMatrix")
  pb_donor <- do_log2cpm(pb_donor, median(Matrix::colSums(pb_donor)))
  #
  pb_donor_meta <- data.frame(donor = colnames(pb_donor))
  pb_donor_meta <- as_tibble(pb_donor_meta)
  pb_donor_meta %<>%
    mutate(
      donor = str_replace(donor, "factor\\(donor\\)", "")
    )
  colnames(pb_donor) <- str_replace(colnames(pb_donor), "factor\\(donor\\)", "")
  if ("drug" %in% colnames(a1$obs)) {
    pb_donor_meta <- left_join(
      pb_donor_meta,
      a1$obs %>%
      select(
        donor, case, drug, sex
      ) %>% unique,
      by = "donor"
    )
  } else {
    pb_donor_meta <- left_join(
      pb_donor_meta,
      a1$obs %>%
      select(
        donor, case, sex
      ) %>% unique,
      by = "donor"
    )
  }
  pb_donor_meta$case <- factor(pb_donor_meta$case, c("Control", "Case"))
  stopifnot(nrow(pb_donor_meta) == ncol(pb_donor))
  qsave(list(meta = pb_donor_meta, log2cpm = pb_donor), file.path(out_dir, "pseudobulk_donor.qs"))
  # x <- table(a1$obs$cluster)
  # min_percent <- 100 * min(x) * 0.25 / sum(x)
  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
  keep_ens <- rownames(pb_donor)[rowMeans(pb_donor) > 0.5]
  #

  models <- list(
    "case" = "~ case",
    "sex-case" = "~ sex + case"
  )
  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    model_form <- models[[i]]
    des1 <- with(pb_donor_meta, model.matrix(
      as.formula(model_form)
    ))
    ob <- as.matrix(pb_donor[keep_ens,])
    fit1 <- lmFit(object = ob, design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
    de_donor <- topTable(
      fit = fit1, coef = "caseCase", number = 1e6, confint = TRUE
    ) %>% rename(Gene = ID)
    de_donor$ensembl_id <- rownames(de_donor)
    de_donor <- as_tibble(de_donor)
    de_donor$percent <- gene_pct[de_donor$ensembl_id]
    # de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
    de_donor_file <- file.path(out_dir, glue("de_donor-{model_name}.tsv.gz"))
    fwrite(x = de_donor, file = de_donor_file, sep = "\t")
    xlsx_file <- file.path(out_dir, glue("de_donor-{model_name}.xlsx"))
    print_status(glue("Writing {xlsx_file}"))
    write_de_xlsx(de_donor %>% mutate_if(is.numeric, signif), xlsx_file)
  }

  extra_genes <- list(
    "n3_2" = c("STAT1", "CD274", "AQP7", "GBP4"),
    "a12_4_4_t4_cd8_1_2" = c("IFNG", "CTLA4"),
    "a12_4_4_t4_cd4_2_2" = c("CXCL13")
  )
  my_de <- de_donor %>%
    filter(percent > 1, abs(logFC) >= log2(1.5), adj.P.Val < 0.05) %>%
    mutate(score = -log10(P.Value) * logFC) %>%
    group_by(logFC > 0) %>%
    top_n(n = 10, wt = abs(score)) %>%
    ungroup %>%
    select(Gene, logFC, P.Value, score)
  my_genes <- my_de$Gene
  if (analysis_name %in% names(extra_genes)) {
    my_genes <- c(my_genes, extra_genes[[analysis_name]])
  }

  n_donors <- table(pb_donor_meta$case)
  p <- plot_limma_volcano(
    de_donor %>% filter(percent > 1),
    # n_text = 20,
    expand_x = expansion(mult = 0.15),
    n_text = 0,
    text_ids = which((de_donor %>% filter(percent > 1))$Gene %in% my_genes),
    text_size = 5
  ) +
    labs(title = glue("Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
  my_ggsave(
    glue("volcano-de_donor"),
    out_dir = glue("{out_dir}"),
    plot = p,
    type = "pdf",
    scale = 1.1, width = 6, height = 5, units = "in", dpi = 300
  )
  my_ggsave(
    glue("volcano-de_donor-tall"),
    out_dir = glue("{out_dir}"),
    plot = p,
    type = "pdf",
    scale = 1.1, width = 6.5, height = 6, units = "in", dpi = 300
  )

  # n_donors <- table(pb_donor_meta$case)
  # p <- plot_limma_volcano(de_donor %>% filter(percent > 1), n_text = 20) +
  #   labs(title = glue("Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
  # my_ggsave(
  #   glue("volcano-de_donor"),
  #   out_dir = glue("{out_dir}"),
  #   plot = p,
  #   type = "pdf",
  #   scale = 1.1, width = 6, height = 5, units = "in", dpi = 300
  # )
  # my_ggsave(
  #   glue("volcano-de_donor-tall"),
  #   out_dir = glue("{out_dir}"),
  #   plot = p,
  #   type = "pdf",
  #   scale = 1.1, width = 6, height = 6, units = "in", dpi = 300
  # )

  # Crypt-axis score
  # The crypt-axis score was assigned for each cell and was based on the
  # expression of a previously defined set of genes: SELENOP, CEACAM7, PLAC8,
  # CEACAM1, TSPAN1, CEACAM5, CEACAM6, IFI27, DHRS9, KRT20, RHOC, CD177, PKIB,
  # HPGD, LYPD8 (Parikh et al., 2019). For each gene, expression within a cell
  # was divided by the max expression across all cells to mitigate the weight
  # of highly expressed genes. The max expression normalized values places the
  # gene’s expression on a 0 to 1 scale. The crypt-axis score is the summation
  # across all genes of the max-normalized expression values.
  if (analysis_name == "n3_2") {
    crypt_genes <- c(
      "SELENOP", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6",
      "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8"
    )
    crypt_genes <- ensembl_to_symbol[ensembl_to_symbol %in% crypt_genes]
    crypt_ens <- names(crypt_genes)
    #
    m <- a1$log2cpm[which(rownames(a1$log2cpm) %in% crypt_ens),]
    gene_maxes <- apply(m, 1, max)
    for (i in seq(nrow(m))) {
      m[i,] <- m[i,] / gene_maxes[i]
    }
    a1$obs$crypt <- colSums(m)
    cluster_groups <- get_cluster_groups(analysis_name)
    a1$obs$cluster_group <- cluster_groups[as.character(a1$obs$leiden)]
    a1$obs$leiden <- factor(as.character(a1$obs$leiden), cluster_order)
    #
    crypt_colors <- sapply(unique(a1$obs$leiden), function(my_leiden) {
      median(a1$obs$crypt[a1$obs$leiden == my_leiden])
    })
    crypt_colors <- colorRamp2(
      # breaks = quantile(crypt_colors, probs = ppoints(11)),
      breaks = seq(0, max(crypt_colors), length.out = 11),
      colors = scico(13)[1:11]
    )(crypt_colors)
    names(crypt_colors) <- unique(a1$obs$leiden)
    #
    p <- ggplot(a1$obs %>% filter(cluster_group != "g5 Mesenchymal cells")) +
      aes(y = leiden, x = crypt, fill = leiden) +
      # aes(y = reorder(leiden, crypt), x = crypt) +
      ggforestplot::geom_stripes() +
      # geom_boxplot(aes(fill = naturalfactor(leiden)), size = 0.3, outlier.size = 0.1) +
      geom_boxplot(size = 0.3, outlier.size = 0.1) +
      # scale_fill_manual(values = cluster_colors, guide = "none") +
      scale_fill_manual(values = crypt_colors, guide = "none") +
      # scale_fill_manual(values = scico(11), guide = "none") +
      facet_grid(rows = vars(cluster_group), scales = "free", space = "free") +
      labs( x = NULL, y = NULL, title = "Crypt-axis Score") +
      theme(
        plot.subtitle = ggtext::element_markdown(),
        panel.spacing = unit(0.5, "lines"),
        # strip.text.y = element_text(angle = 0, hjust = 0)
        strip.text.y = element_blank()
      )
    my_ggsave(
      glue("crypt-score-small"),
      out_dir = glue("{out_dir}"),
      plot = p,
      type = "pdf",
      scale = 1.1, width = 1, height = 6.5, units = "in", dpi = 300
    )
    p <- plot_hexgene(
      x = a1$obs$UMAP1,
      y = a1$obs$UMAP2,
      z = a1$obs$crypt,
      text = FALSE,
      bins = 57,
      palette = "bilbao", direction = 1
    ) +
    labs(title = "Crypt-axis score") +
    guides(fill = guide_colorbar(title = NULL, barwidth = 10)) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "plain")
    )
    my_ggsave(
      glue("crypt-score-umap2"),
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 0.8,
      width = 3.5,
      height = 3 + 1,
      units = "in",
      dpi = 300
    )
  }


  ## NRG1 fusion partners
  #nrg1_fusion <- strsplit("CD74 ATP1B1 SDC4 ADAM9 CDH1 COX10-AS1 DIP2B DPYSL2 GDF15 HMBOX1 MDK-NRG1 MRPL13 NOTCH2 PARP8 POMK ROCK1 SETD4 SLC3A2 TNC TSHZ2 VTCN1 WHSC1L1 ZMYM2", " ")[[1]]
  #de_donor %>% filter(Gene %in% nrg1_fusion)
  ##
  #n_overlap <- nrow(de_donor %>% filter(Gene %in% nrg1_fusion, logFC > log2(1.5), adj.P.Val < 0.05))
  #n_drawn <- length(nrg1_fusion)
  #n_up <- nrow(de_donor %>% filter(logFC > 0, adj.P.Val < 0.05))
  #n_total <- length(unique(c(nrg1_fusion, de_donor$Gene)))
  #length(unique(de_donor$Gene))
  #phyper(n_overlap - 1, n_up, n_total - n_up, n_drawn, lower.tail = FALSE)

  # top_de_donor_ids <- (
  #   de_donor %>%
  #   filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
  #   arrange(P.Value)
  #   # arrange(-logFC)
  # )$ensembl_id %>% head(20)

  top_de_donor_ids <- (
    de_donor %>%
    filter(
      abs(logFC) > log2(1.5) & adj.P.Val < 0.05 & percent > 0.5
    ) %>%
    # top_n(wt = abs(logFC) * -log10(P.Value), n = 60) %>%
    arrange(-logFC * -log10(P.Value))
  )$ensembl_id

  #d <- as_tibble(reshape2::melt(as.matrix(pb_donor)[my_ids,]))
  #colnames(d) <- c("ensembl_id", "donor", "value")
  #d$GeneName <- factor(
  #  as.character(unname(ensembl_to_symbol[as.character(d$ensembl_id)])),
  #  rev(as.character(unname(ensembl_to_symbol[as.character(my_ids)])))
  #)
  #d$Gene <- as.integer(d$GeneName)
  #d <- left_join(d, pb_donor_meta, by = "donor")
  ##
  #p <- ggplot(d) +
  #    annotate(
  #      geom = "rect",
  #      xmin = -Inf,
  #      xmax = Inf,
  #      ymin = seq(from = 1, to = length(my_ids), by = 2) - 0.5,
  #      ymax = seq(from = 1, to = length(my_ids), by = 2) + 0.5,
  #      alpha = 0.2
  #    ) +
  #    # geom_vline(xintercept = 0, size = 0.3) +
  #    geom_point(
  #      data = d,
  #      mapping = aes(x = value, y = Gene, group = case, fill = case),
  #      position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE),
  #      shape = 21, size = 2, stroke = 0.2
  #    ) +
  #    scale_fill_manual(values = okabe(8)) +
  #    # scale_x_continuous(
  #    #   breaks = seq(-10, 10, by = 1),
  #    #   labels = function(x) fractional::fractional(2 ^ x)
  #    # ) +
  #    scale_y_continuous(
  #      expand = c(0, 0),
  #      breaks = seq(1, max(d$Gene)),
  #      labels = levels(d$GeneName)
  #    ) +
  #    guides(fill = guide_legend(reverse = TRUE, title = NULL, override.aes = list(size = 5))) +
  #    # annotation_logticks(sides = "b", size = 0.3) +
  #    expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
  #    labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
  #    theme(
  #      axis.text.y = element_text(face = "italic"),
  #      legend.margin = margin(0, 0, 0, 0),
  #      legend.box.margin = margin(0, 0, 0, 0)
  #    )
  #my_ggsave(
  #  "de_donor",
  #  out_dir = out_dir,
  #  type = "pdf",
  #  plot = p,
  #  scale = 1,
  #  width = 5,
  #  height = length(my_ids) * 0.3 + 1,
  #  units = "in",
  #  dpi = 300
  #)

  umap_n_bins <- 37
  if (nrow(a1$obs) > 5e4) {
    umap_n_bins <- 57
  }
  if (nrow(a1$obs) > 1e5) {
    umap_n_bins <- 97
  }

  plot_umap_figures <- FALSE

  if (plot_umap_figures) {

    for (my_id in head(top_de_donor_ids, 50)) {
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        group = append_n(a1$obs$case),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      facet_wrap(~ group) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/de_donor"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 7,
        height = 5,
        units = "in",
        dpi = 300
      )
    }

    extra_genes <- c("PDCD1", "CD274", "PDCD1LG2","CXCL13", "IFNG", "IL17A", "IL26",
      "HIF1A", "CX3CL1", "S1PR1")
    my_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in% extra_genes])
    # for (my_id in c(head(top_de_donor_ids, 50), my_ids)) {
    for (my_id in my_ids) {
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        group = append_n(a1$obs$case),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      facet_wrap(~ group) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/de_donor/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 6.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }

  }

  # Pseudobulk at the cluster level
  ########################################################################
  y <- with(a1$obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/") # means
  pb <- as(a1$counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
  #
  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
  colnames(pb_meta) <- c("cluster", "donor")
  pb_meta <- as_tibble(pb_meta)
  pb_meta %<>%
    mutate(
      cluster = str_replace(cluster, "factor\\(cluster\\)", ""),
      donor = str_replace(donor, "factor\\(donor\\)", "")
    )
  if ("drug" %in% colnames(a1$obs)) {
    pb_meta <- left_join(
      pb_meta,
      a1$obs %>%
      select(
        donor, case, drug, sex
      ) %>% unique,
      by = "donor"
    )
  } else {
    pb_meta <- left_join(
      pb_meta,
      a1$obs %>%
      select(
        donor, case, sex
      ) %>% unique,
      by = "donor"
    )
  }
  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
  stopifnot(nrow(pb_meta) == ncol(pb))
  #
  y <- with(a1$obs, model.matrix(~ 0 + factor(cluster)))
  y <- as(y, "dgCMatrix")
  pct <- as((a1$counts > 0) %*% y, "dgCMatrix")
  pct <- sweep(pct, 2, colSums(y), "/")
  pct <- as(pct, "dgCMatrix")
  colnames(pct) <- str_replace(colnames(pct), "factor\\(cluster\\)", "")
  pct_rows <- rownames(pct)
  pct_cols <- colnames(pct)
  pct <- summary(pct)
  pct[[1]] <- pct_rows[pct[[1]]]
  pct[[2]] <- pct_cols[pct[[2]]]
  colnames(pct) <- c("ensembl_id", "cluster", "percent")
  pct$percent <- pct$percent * 100
  pct <- as.data.table(pct)
  # Keep the gene if it is expressed in >1% of cells within at least one cluster.
  keep_ens <- pct[,.(max=max(percent)),by=ensembl_id][max >= 1]$ensembl_id
  #
  qsave(list(meta = pb_meta, log2cpm = pb), file.path(out_dir, "pseudobulk_donor_cluster.qs"))

  models <- list(
    "case" = "~ case",
    "sex-case" = "~ sex + case"
  )
  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    model_form <- models[[i]]
    de <- rbindlist(
      pblapply(sort(unique(pb_meta$cluster)), function(this_cluster) {
        # my_donors <- comp$donor[comp$cluster == this_cluster & comp$n >= 10]
        ix_cluster <- pb_meta$cluster == this_cluster #& pb_meta$donor %in% my_donors
        des1 <- with(pb_meta[ix_cluster,], model.matrix(
          as.formula(model_form)
        ))
        # keep_ens <- pct[cluster == this_cluster][percent >= 1]$ensembl_id
        ob <- as.matrix(pb[keep_ens,ix_cluster])
        # ob <- ob[rowMeans(ob) > 0.5,] # we lose IL17A in CD8 T cells if we use this filter
        fit1 <- lmFit(object = ob, design = des1)
        fit1 <- eBayes(fit1)
        fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
        res <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE)
        res$cluster <- this_cluster
        res <- res %>% rename(Gene = ID)
        res$ensembl_id <- rownames(res)
        return(res)
      })
    )
    de$GeneName <- de$Gene
    de$cluster <- naturalfactor(de$cluster)
    # de_file <- glue("{out_dir}/de_case-vs-control.tsv.gz")
    de <- left_join(de, pct, by = c("cluster", "ensembl_id"))
    de$percent[is.na(de$percent)] <- 0
    fwrite(de, glue("{out_dir}/de_{model_name}.tsv.gz"), sep = "\t")
    xlsx_file <- glue("{out_dir}/de_{model_name}.xlsx")
    print_status(glue("Writing {xlsx_file}"))
    write_de_xlsx(
      de %>% select(-GeneName) %>% mutate_if(is.numeric, signif),
      xlsx_file,
      col = "cluster"
    )
  }

  de_summary <- de %>% group_by(cluster) %>%
    filter(percent > 1) %>%
    summarize(
      n_up = sum(logFC > log2(1.5) & adj.P.Val < 0.05 & AveExpr > 0.5),
      n_down = sum(-logFC > log2(1.5) & adj.P.Val < 0.05 & AveExpr > 0.5)
    )
  de_summary_file <- glue("{out_dir}/de_summary_case-vs-control.tsv")
  fwrite(de_summary, de_summary_file, sep = "\t")

  file_cluster_order <- file.path(out_dir, "../composition-case-vs-control/cluster_order.rds")
  if (file.exists(file_cluster_order)) {
    cluster_order <- readRDS(file_cluster_order)
  } else {
    d_masc <- fread(glue("{out_dir}/../composition-case-vs-control/masc_1-case-Case.tsv"))
    cluster_order <- rev(factor(d_masc$cluster))
  }
  #
  cluster_groups <- get_cluster_groups(analysis_name)
  if (!is.null(cluster_groups)) {
    de_summary$cluster_group <- naturalfactor(cluster_groups[as.character(de_summary$cluster)])
    de_summary$cluster <- factor(as.character(de_summary$cluster), cluster_order)
    p <- ggplot(
      de_summary %>% pivot_longer(cols = c("n_up", "n_down")) %>%
        mutate(value = ifelse(name == "n_down", -value, value))
    ) +
      aes(x = value, y = cluster, fill = name) +
      ggforestplot::geom_stripes() +
      geom_vline(xintercept = 0, size = 0.3) +
      geom_colh() +
      # scale_y_discrete(limits = rev(levels(naturalfactor(de_summary$cluster)))) +
      facet_grid(rows = vars(cluster_group), scales = "free", space = "free") +
      scale_fill_manual(
        # values = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[c(9,3)],
        values = okabe(2),
        guide = "none"
      ) +
      scale_x_continuous(labels = abs) +
      labs(x = NULL, y = NULL) +
      theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())
    my_ggsave(
      "bars-summary",
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1,
      width = 3,
      height = nrow(de_summary) * 0.2 + 0.5,
      units = "in",
      dpi = 300
    )
    my_ggsave(
      "bars-summary-narrow",
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1,
      width = 1.5,
      height = nrow(de_summary) * 0.2 + 1.5,
      units = "in",
      dpi = 300
    )
    # my_ggsave(
    #   "bars-summary-order",
    #   out_dir = out_dir,
    #   type = "pdf",
    #   plot = p + scale_y_discrete(limits = cluster_order),
    #   scale = 1,
    #   width = 2,
    #   height = nrow(de_summary) * 0.2 + 0.5,
    #   units = "in",
    #   dpi = 300
    # )
  }

  for (this_cluster in unique(pb_meta$cluster)) {
    ix <- pb_meta$cluster == this_cluster
    n_donors <- table(pb_meta$case[ix])
    #
    p <- plot_limma_volcano(de[de$cluster == this_cluster,]) +
      labs(title = glue("{this_cluster}: Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
    my_ggsave(
      glue("volcano-cluster-{this_cluster}"),
      out_dir = glue("{out_dir}/volcano-cluster"),
      plot = p,
      type = "pdf",
      scale = 1, width = 5, height = 4, units = "in", dpi = 300
    )
    #
    my_ids <- c(
      (
        de %>%
        filter(cluster == this_cluster) %>%
        filter(logFC > 0) %>%
        top_n(n = 20, wt = -log10(P.Value) * logFC)
      )$ensembl_id,
      (
        de %>%
        filter(cluster == this_cluster) %>%
        filter(logFC < 0) %>%
        top_n(n = 20, wt = -log10(P.Value) * -logFC)
      )$ensembl_id
    )
    # my_ids <- head((
    #   de %>%
    #   filter(cluster == this_cluster) %>%
    #   filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
    #   arrange(abs(logFC) + -log10(P.Value))
    # )$ensembl_id, 50)
    mat <- as.matrix(pb[my_ids,ix])
    a_col <- as.data.frame(pb_meta)[ix,]
    if (all(c("donor", "drug") %in% colnames(a1$obs))) {
      # a_col <- left_join(a_col, unique(a1$obs[,c("donor","drug")]), by = "donor")
      rownames(a_col) <- colnames(mat)
      a_col <- a_col %>% select(-cluster)
      h1 <- dendsort(hclust(dist(mat), method = "complete"))
      h2 <- dendsort(hclust(dist(t(mat)), method = "complete"))
      a_col <- a_col[h2$order,]
      mat <- mat[h1$order,h2$order]
      #
      a_colors <- list()
      a_colors[["drug"]] <- okabe(8)[4:7]
      names(a_colors[["drug"]]) <- c("None", "PD-1", "CTLA-4", "PD-1/CTLA-4")
      a_colors[["case"]] <- okabe(2)
      names(a_colors[["case"]]) <- c("Control", "Case")
      #
      heatmap_file <- glue("{out_dir}/heatmap-cluster/heatmap-cluster-{this_cluster}.pdf")
      dir.create(dirname(heatmap_file), recursive = TRUE, showWarnings = FALSE)
      message(heatmap_file)
      try({
      pheatmap::pheatmap(
        filename = heatmap_file,
        width = ncol(mat) * 0.2 + 3,
        height = nrow(mat) * 0.1 + 3,
        mat = mat,
        cluster_col = FALSE,
        cluster_row = FALSE,
        hclust_method = "average",
        show_colnames = FALSE,
        show_rownames = TRUE,
        scale = "row",
        # color = rev(scico::scico(palette = "roma", n = 20)),
        color = rev(RColorBrewer::brewer.pal("RdBu", n = 11)),
        border_color = NA,
        labels_row = ensembl_to_symbol[rownames(mat)],
        # labels_col = labels_col,
        annotation_col = a_col[,c("drug","case")],
        annotation_colors = a_colors,
        main = glue("Cluster {this_cluster}")
      )
      })
    }
  }

  de_top_genes <- unique(c(
    (
      de %>% group_by(cluster) %>%
      filter(
        abs(logFC) > log2(1.5) & adj.P.Val < 0.01
      ) %>% ungroup
    )$ensembl_id
    # (
    #   de_donor %>%
    #   filter(logFC > log2(2), adj.P.Val < 0.05)
    # )$ensembl_id
  ))
  if (length(de_top_genes)) {
    de_top <- de %>% filter(ensembl_id %in% de_top_genes)
    # "IL17A" %in% de_top$Gene
    #
    mat <- dcast(data = de_top, formula = ensembl_id ~ cluster, value.var = "logFC")
    # mat <- dcast(data = de_top, formula = ensembl_id ~ cluster, value.var = "AveExpr")
    mat_rows <- mat[[1]]
    mat[[1]] <- NULL
    mat <- as.matrix(mat)
    mat[is.na(mat)] <- 0
    rownames(mat) <- ensembl_to_symbol[mat_rows]
    #
    # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
    my_genes <- (de %>% arrange(P.Value) %>% head(30))$Gene
    row_ha <- rowAnnotation(
      gene = anno_mark(
        at = match(my_genes, rownames(mat)),
        labels = my_genes,
        labels_gp = gpar(fontface = "italic")
      )
    )
    cluster_colors <- mpn65[1:length(unique(de$cluster))]
    names(cluster_colors) <- as.character(seq_along(cluster_colors))
    column_ha <- HeatmapAnnotation(
      cluster = as.character(colnames(mat)),
      col = list(cluster = cluster_colors)
    )
    col_fun <- circlize::colorRamp2(
      seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
      rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
    )
    ht <- Heatmap(
      matrix = mat,
      name = "Log2FC",
      top_annotation = column_ha,
      right_annotation = row_ha,
      # col = rev(scico(pal = "", n = 20)),
      # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
      col = col_fun,
      row_split = 5,
      show_row_names = FALSE,
      # show_row_dend = FALSE,
      show_column_names = TRUE,
      column_names_side = "top",
      column_names_rot = 90,
      # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
      show_column_dend = FALSE
    )
    pdf_file <- file.path(out_dir, "heatmap-de-top.pdf")
    message(pdf_file)
    unlink(pdf_file)
    pdf(pdf_file)
    draw(ht)
    dev.off()
  }

  if (length(top_de_donor_ids) > 1) {
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor,
      de,
      head(top_de_donor_ids, 20),
      glue("{out_dir}/heatmap_top_de_donor.pdf")
    )
  }

  top_de_donor_ids2 <- c(
    (
      de_donor %>%
      filter(
        logFC > log2(1.5) & adj.P.Val < 0.05 & Gene %in% de$Gene
      ) %>%
      top_n(wt = abs(logFC) * -log10(P.Value), n = 10) %>%
      arrange(-logFC)
    )$ensembl_id,
    (
      de_donor %>%
      filter(
        logFC < log2(1.5) & adj.P.Val < 0.05 & Gene %in% de$Gene
      ) %>%
      top_n(wt = abs(logFC) * -log10(P.Value), n = 10) %>%
      arrange(logFC)
    )$ensembl_id
  )
  if (length(top_de_donor_ids2) > 1) {
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor,
      de,
      top_de_donor_ids2,
      glue("{out_dir}/heatmap_top_de_donor2.pdf"),
      cluster_order = rev(cluster_order)
    )
  }

  max_fdr <- 0.05
  # if (analysis_name == "a12_4_4_b5_1_3") {
  #   max_fdr <- 0.05
  # }
  top_de_ids <- (
    de %>% group_by(cluster) %>%
    filter(
      abs(logFC) > log2(1.5) & adj.P.Val < max_fdr
    ) %>%
    top_n(wt = abs(logFC) * -log10(P.Value), n = 5) %>%
    arrange(-logFC)
  )$ensembl_id %>% unique
  if (length(top_de_ids) > 1) {
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor,
      de,
      top_de_ids,
      glue("{out_dir}/heatmap_top_de.pdf"),
      cluster_order = rev(cluster_order)
    )
  }

  if (FALSE) {
    umap_top_de_ids <- (
      de %>% group_by(cluster) %>%
      filter(
        abs(logFC) > log2(1.5) & adj.P.Val < max_fdr
      ) %>%
      top_n(wt = abs(logFC) * -log10(P.Value), n = 20) %>%
      arrange(-logFC)
    )$ensembl_id %>% unique
    if (length(umap_top_de_ids) > 1) {
      for (my_id in umap_top_de_ids) {
        p <- plot_hexgene(
          x = a1$obs$UMAP1,
          y = a1$obs$UMAP2,
          z = as.numeric(a1$log2cpm[my_id,]),
          group = append_n(a1$obs$case),
          bins = umap_n_bins,
          palette = "oslo"
        ) +
        facet_wrap(~ group) +
        labs(title = ensembl_to_symbol[my_id]) +
        theme(
          panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(face = "italic")
        )
        my_ggsave(
          glue("umap-{safe(ensembl_to_symbol[my_id])}"),
          out_dir = file.path(out_dir, "umap/de_cluster"),
          type = "pdf",
          plot = p,
          scale = 0.8,
          width = 7,
          height = 5,
          units = "in",
          dpi = 300
        )
      }
      for (my_id in umap_top_de_ids) {
        p <- plot_hexgene(
          x = a1$obs$UMAP1,
          y = a1$obs$UMAP2,
          z = as.numeric(a1$log2cpm[my_id,]),
          group = append_n(a1$obs$case),
          bins = umap_n_bins,
          palette = "oslo"
        ) +
        facet_wrap(~ group) +
        labs(title = ensembl_to_symbol[my_id]) +
        theme(
          legend.position = "none",
          strip.text = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(face = "italic")
        )
        my_ggsave(
          glue("umap-{safe(ensembl_to_symbol[my_id])}"),
          out_dir = file.path(out_dir, "umap/de_cluster/slim"),
          type = "pdf",
          plot = p,
          scale = 0.8,
          width = 6.5,
          height = 3,
          units = "in",
          dpi = 300
        )
      }
    }
  }

  if (analysis_name == "a12_4_4_t4_cd4_2_2") {

    filename <- "data/colitis/CD4_DGE_lists.xlsx"
    sheets <- readxl::excel_sheets(filename)
    for (sheet in sheets) {
      message(sheet)
      # x <- rbindlist(lapply(sheets, function(sheet) {
      x <- readxl::read_excel(filename, sheet = sheet)[,1:2]
      colnames(x) <- c("gene", "cat")
      # })) %>% unique()
      x <- left_join(
        x = x,
        y = de_donor %>% select(Gene, ensembl_id, P.Value, logFC, CI.L, CI.R),
        by = c("gene" = "Gene")
      )
      x$cat <- str_replace(x$cat, "&", "and")
      x$cat <- str_split_fixed(x$cat, ",", 2)[,1]
      x <- x %>% group_by(cat) %>% mutate(wt = abs(logFC) * -log10(P.Value)) %>% arrange(-wt) %>% ungroup()
      x <- x[!duplicated(x$ensembl_id),]
      # cat_ids <- split(x$ensembl_id, x$cat)
      #
      # for (cat_name in names(cat_ids)) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          # cat_ids[[cat_name]],
          x$ensembl_id,
          glue("{out_dir}/cat/heatmap_{safe(sheet)}.pdf"),
          # title = cat_name
          title = sheet
        )
      # }
    }

    my_cats <- list(
      "coinhibitory" = str_split("MYO7A CTLA4 MT1E LAIR2 HAVCR2 PDCD1 MT2A TOX2 CD160 CD38 BTLA LAG3 TOX LAIR1 ENTPD1 CD244 TIGIT", " ")[[1]],
      "costimulatory" = str_split("TNFRSF8 TNFRSF4 CD27 TNFRSF25 TNFRSF14 CD28 SLAMF1 CD40LG TNFRSF9 TNFSF14 CD2 CD226 ICOS TNFRSF18", " ")[[1]]
    )
    for (cat in names(my_cats)) {
      x <- left_join(x = data.frame(Gene = my_cats[[cat]]), y = de_donor, by = "Gene")
      plot_de_summary(
        pb_donor_meta,
        pb_donor,
        de_donor,
        de,
        x$ensembl_id,
        glue("{out_dir}/cat/heatmap_{safe(cat)}.pdf"),
        title = cat
      )
    }

    # for (my_cluster in unique(de$cluster)) {
    #   n <- 25
    #   my_ens <- (
    #     de %>% filter(cluster == my_cluster) %>%
    #       filter(adj.P.Val < 0.05) %>%
    #       mutate(wt = abs(logFC) * -log10(P.Value)) %>%
    #       top_n(wt = wt, n = n)  %>%
    #       arrange(-logFC)
    #   )$ensembl_id
    #   if (length(my_ens) > 1) {
    #     plot_de_summary(
    #       pb_donor_meta,
    #       pb_donor,
    #       de_donor,
    #       de,
    #       my_ens,
    #       glue("{out_dir}/cat/heatmap_top_genes_cluster{my_cluster}.pdf"),
    #       title = glue("Top {n} genes in cluster {my_cluster}")
    #     )
    #   }
    # }

    my_cats <- list(
      # "cytokine" = str_split("IL10 IL26 IL17A IL17F IL22 EBI3 IL23R IL12RB1 IL1R2 STAT1 JAK3 IFNG GZMB GZMH CXCL13 CCL3 CCL4 CXCR6 ITGB2 ITGAX ITGA2 CADM1", " ")[[1]],
      "upin9" = c("IL1R2", "CD74", "IL21R", "IL12RB2", "CXCR6", "CXCR3", "CD7",
        "IL12RB1", "IL18R1", "CCL20", "ILF2", "CD320", "CD58", "ICAM2",
        "CD226", "ITGAL", "LGALS1", "CD81", "TBX21", "LGALS3"),
      "upin7" = c("CXCL13", "MAF", "IL17A", "IFNG", "RBPJ", "IL21", "CCL4", "SLAMF1",
        "JAK3", "ITGB1", "NOTCH1", "IL4R", "CD151"),
      "coinhibitory" = str_split("MYO7A CTLA4 MT1E LAIR2 HAVCR2 PDCD1 MT2A TOX2 CD160 CD38 BTLA LAG3 TOX LAIR1 ENTPD1 CD244 TIGIT", " ")[[1]],
      "costimulatory" = str_split("TNFRSF8 TNFRSF4 CD27 TNFRSF25 TNFRSF14 CD28 SLAMF1 CD40LG TNFRSF9 TNFSF14 CD2 CD226 ICOS TNFRSF18", " ")[[1]]
    )
    my_cats_d <- rbindlist(lapply(names(my_cats), function(my_cat) {
      data.frame(cat = my_cat, Gene = my_cats[[my_cat]])
    }))
    my_cats_d <- my_cats_d[!duplicated(my_cats_d[['Gene']]),]
    my_cats_d$cat <- factor(my_cats_d$cat, names(my_cats))
    x <- left_join(x = my_cats_d, y = de_donor, by = "Gene")
    x <- x %>%
      filter(!is.na(ensembl_id)) %>%
      mutate(wt = abs(logFC) * -log10(P.Value)) %>%
      arrange(cat, -wt)
    de_x <- left_join(x = my_cats_d, y = de, by = "Gene")
    de_donor_x <- left_join(x =  my_cats_d, y = de_donor, by = "Gene")
    cluster_groups <- get_cluster_groups(analysis_name)
    if (length(cluster_groups) == length(unique(de_x$cluster))) {
      de_x$cluster_group <- cluster_groups[as.character(de_x$cluster)]
    }
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      x$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats1.pdf"),
      cluster_order = rev(cluster_order),
      title = "cats"
    )

    filename <- "data/colitis/CD4_DGE_lists.xlsx"
    x <- readxl::read_excel(filename, sheet = "supplemental-figure")[,1:2]
      colnames(x) <- c("Gene", "cat")
    x$cat <- str_replace(x$cat, "&", "and")
    x$cat <- str_split_fixed(x$cat, ",", 2)[,1]
    # x <- left_join(
    #   x = x,
    #   y = de_donor %>% select(Gene, ensembl_id),
    #   by = "Gene"
    # )
    x <- x[!duplicated(x[[1]]),]
    # # x <- x %>% group_by(cat) %>% mutate(wt = abs(logFC) * -log10(P.Value)) %>% arrange(-wt) %>% ungroup()
    # x <- x %>%
    #   filter(!is.na(ensembl_id)) %>%
    #   mutate(wt = abs(logFC) * -log10(P.Value)) %>%
    #   arrange(cat, -wt)
    de_x <- left_join(x = x, y = de, by = "Gene")
    de_donor_x <- left_join(x =  x, y = de_donor, by = "Gene")
    if (length(cluster_groups) == length(unique(de_x$cluster))) {
      de_x$cluster_group <- cluster_groups[as.character(de_x$cluster)]
    }
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      unique(de_x$ensembl_id),
      glue("{out_dir}/cat/heatmap_cats2.pdf"),
      cluster_order = rev(cluster_order),
      title = "cats"
    )

  }

  if (analysis_name == "a12_4_4_t4_cd8_1_2") {

    # selected_genes <- 
    filename <- "data/colitis/CD8_DGE_lists.xlsx"
    sheets <- readxl::excel_sheets(filename)
    for (sheet in sheets) {
      # x <- rbindlist(lapply(sheets, function(sheet) {
        x <- readxl::read_excel(filename, sheet = sheet)[,1:2]
        colnames(x) <- c("gene", "cat")
        # x
      # })) %>% unique()
      x <- left_join(
        x = x,
        y = de_donor %>% select(Gene, ensembl_id, P.Value, logFC, CI.L, CI.R),
        by = c("gene" = "Gene")
      )
      # x <- x %>% group_by(cat) %>% arrange(P.Value) %>% ungroup()
      # cat_ids <- split(x$ensembl_id, x$cat)
      x <- x %>% group_by(cat) %>%
        mutate(wt = abs(logFC) * -log10(P.Value)) %>%
        ungroup() %>%
        arrange(cat, -wt)
      x <- x %>%
        filter(
          ensembl_id %in% de$ensembl_id,
          ensembl_id %in% de_donor$ensembl_id
        )
      #
      # for (cat_name in names(cat_ids)) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          # cat_ids[[cat_name]],
          # glue("{out_dir}/heatmap_{safe(cat_name)}.pdf"),
          # title = cat_name
          x$ensembl_id,
          glue("{out_dir}/cat/heatmap_{safe(sheet)}.pdf"),
          title = sheet
        )
      # }
    }

    filename <- "data/colitis/CD8_DGE_lists.xlsx"
    x <- readxl::read_excel(filename, sheet = "supplement heatmaps")[,1:2]
      colnames(x) <- c("Gene", "cat")
    x$cat <- str_replace(x$cat, "&", "and")
    x$cat <- str_split_fixed(x$cat, ",", 2)[,1]
    x <- x[!duplicated(x[[1]]),]
    de_x <- left_join(x = x, y = de, by = "Gene")
    de_donor_x <- left_join(x =  x, y = de_donor, by = "Gene")
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      (de_donor_x %>% arrange(-logFC))$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats2.pdf"),
      cluster_order = rev(cluster_order),
      title = "cats"
    )

    my_cats <- list(
      "coinhibitory" = str_split("MYO7A CTLA4 MT1E LAIR2 HAVCR2 PDCD1 MT2A TOX2 CD160 CD38 BTLA LAG3 TOX LAIR1 ENTPD1 CD244 TIGIT", " ")[[1]],
      "costimulatory" = str_split("TNFRSF8 TNFRSF4 CD27 TNFRSF25 TNFRSF14 CD28 SLAMF1 CD40LG TNFRSF9 TNFSF14 CD2 CD226 ICOS TNFRSF18", " ")[[1]]
    )
    for (cat in names(my_cats)) {
      x <- left_join(x = data.frame(Gene = my_cats[[cat]]), y = de_donor, by = "Gene")
      plot_de_summary(
        pb_donor_meta,
        pb_donor,
        de_donor,
        de,
        x$ensembl_id,
        glue("{out_dir}/cat/heatmap_{safe(cat)}.pdf"),
        title = cat
      )
    }

    for (my_cluster in unique(de$cluster)) {
      n <- 25
      my_ens <- (
        de %>% filter(cluster == my_cluster) %>%
          filter(adj.P.Val < 0.05) %>%
          mutate(wt = abs(logFC) * -log10(P.Value)) %>%
          top_n(wt = wt, n = n)  %>%
          arrange(-logFC)
      )$ensembl_id
      if (length(my_ens) > 1) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          my_ens,
          glue("{out_dir}/cat/heatmap_top_genes_cluster{my_cluster}.pdf"),
          title = glue("Top {n} genes in cluster {my_cluster}")
        )
      }
    }

    # for (cat_name in names(cat_ids)) {
    #   pb_meta[['cat']] <- colMeans(pb[cat_ids[[cat_name]],])
    #   pb_meta$cluster_o <- factor(pb_meta$cluster, rev(naturalsort(unique(pb_meta$cluster))))
    #   p <- ggplot(pb_meta) +
    #     aes(x = cat, y = cluster_o, fill = case, group = case) +
    #     ggforestplot::geom_stripes() +
    #     geom_point(
    #       shape = 21, size = 3,
    #       position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE)
    #     ) +
    #     scale_fill_manual(values = pals::okabe()) +
    #     labs(x = bquote('Mean log'[2]~'CPM'), y = NULL, title = cat_name)
    #   my_ggsave(
    #     slug = glue("quasirandom_{safe(cat_name)}.pdf"),
    #     out_dir = out_dir,
    #     type = "pdf",
    #     plot = p,
    #     scale = 1, width = 5, height = 5, units = "in", dpi = 300
    #   )
    # }

    # # selected_genes <- 
    # filename <- "data/CD8_DGE_lists.xlsx"
    # sheets <- readxl::excel_sheets(filename)[1:3]
    # x <- rbindlist(lapply(sheets, function(sheet) {
    #   x <- readxl::read_excel(filename, sheet = sheet)[,1:2]
    #   colnames(x) <- c("gene", "cat")
    #   x
    # })) %>% unique()
    # x <- left_join(
    #   x = x,
    #   y = de %>% select(cluster, Gene, ensembl_id, P.Value, logFC, CI.L, CI.R),
    #   by = c("gene" = "Gene")
    # )
    # x <- x %>% group_by(cat) %>% arrange(P.Value) %>% ungroup()
    # cat_ids <- split(x$ensembl_id, x$cat)
    # x %<>% group_by(cat, cluster) %>%
    #   summarize(
    #     logFC = mean(logFC),
    #     CI.L = mean(CI.L),
    #     CI.R = mean(CI.R),
    #     .groups = "drop"
    #   )

    # x$cluster <- factor(x$cluster, rev(naturalsort(unique(x$cluster))))
    # p <- ggplot(x) +
    #   aes(x = logFC, y = cluster) +
    #   ggforestplot::geom_stripes() +
    #   geom_vline(xintercept = 0, size = 0.3) +
    #   geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), height = 0) +
    #   scale_x_continuous(
    #     labels = function(x) fractional::fractional(2 ^ x)
    #   ) +
    #   geom_point() +
    #   scale_fill_manual(values = pals::okabe()) +
    #   facet_wrap(~ cat, ncol = length(names(cat_ids))) +
    #   labs(x = 'FC', y = NULL) +
    #   theme(panel.spacing = unit(1, "lines"))
    # my_ggsave(
    #   slug = glue("errorbarh_cats.pdf"),
    #   out_dir = out_dir,
    #   type = "pdf",
    #   plot = p,
    #   scale = 1, width = 5, height = 5, units = "in", dpi = 300
    # )

    my_cats <- list(
      "cytokine" = str_split("IL10 IL26 IL17A IL22 EBI3 IL23R IL12RB1 IL1R2 STAT1 JAK3 IFNG GZMB GZMH CXCL13 CCL3 CCL4 CXCR6 ITGB2 ITGAX ITGA2 CADM1", " ")[[1]],
      "coinhibitory" = str_split("MYO7A CTLA4 MT1E LAIR2 HAVCR2 PDCD1 MT2A TOX2 CD160 CD38 BTLA LAG3 TOX LAIR1 ENTPD1 CD244 TIGIT", " ")[[1]],
      "costimulatory" = str_split("TNFRSF8 TNFRSF4 CD27 TNFRSF25 TNFRSF14 CD28 SLAMF1 CD40LG TNFRSF9 TNFSF14 CD2 CD226 ICOS TNFRSF18", " ")[[1]]
    )
    my_cats_d <- rbindlist(lapply(names(my_cats), function(my_cat) {
      data.frame(cat = my_cat, Gene = my_cats[[my_cat]])
    }))
    my_cats_d$cat <- factor(my_cats_d$cat, names(my_cats))
    x <- left_join(x = my_cats_d, y = de_donor, by = "Gene")
    x <- x %>%
      mutate(wt = abs(logFC) * -log10(P.Value)) %>%
      arrange(cat, -wt)
    de_x <- left_join(x = my_cats_d, y = de, by = "Gene")
    de_donor_x <- left_join(x =  my_cats_d, y = de_donor, by = "Gene")
    cluster_groups <- get_cluster_groups(analysis_name)
    if (length(cluster_groups) == length(unique(de_x$cluster))) {
      de_x$cluster_group <- cluster_groups[as.character(de_x$cluster)]
    }
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      x$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats1.pdf"),
      cluster_order = rev(cluster_order),
      title = "cats"
    )

  }

  if (analysis_name == "a12_4_4_b5_1_3") {

    filename <- "data/colitis/B_DGE_lists.xlsx"
    sheets <- readxl::excel_sheets(filename)
    for (sheet in sheets) {
      message(sheet)
      # x <- rbindlist(lapply(sheets, function(sheet) {
      x <- readxl::read_excel(filename, sheet = sheet)[,1:2]
      colnames(x) <- c("gene", "cat")
      # })) %>% unique()
      x <- left_join(
        x = x,
        y = de_donor %>% select(Gene, ensembl_id, P.Value, logFC, CI.L, CI.R),
        by = c("gene" = "Gene")
      )
      x$cat <- str_replace(x$cat, "&", "and")
      x$cat <- str_split_fixed(x$cat, ",", 2)[,1]
      x <- x %>% group_by(cat) %>% mutate(wt = abs(logFC) * -log10(P.Value)) %>% arrange(-wt) %>% ungroup()
      # cat_ids <- split(x$ensembl_id, x$cat)
      #
      # for (cat_name in names(cat_ids)) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          # cat_ids[[cat_name]],
          x$ensembl_id,
          glue("{out_dir}/cat/heatmap_{safe(sheet)}.pdf"),
          # title = cat_name
          title = sheet
        )
      # }
    }

    # my_cats <- list(
    #   "coinhibitory" = str_split("MYO7A CTLA4 MT1E LAIR2 HAVCR2 PDCD1 MT2A TOX2 CD160 CD38 BTLA LAG3 TOX LAIR1 ENTPD1 CD244 TIGIT", " ")[[1]],
    #   "costimulatory" = str_split("TNFRSF8 TNFRSF4 CD27 TNFRSF25 TNFRSF14 CD28 SLAMF1 CD40LG TNFRSF9 TNFSF14 CD2 CD226 ICOS TNFRSF18", " ")[[1]]
    # )
    # for (cat in names(my_cats)) {
    #   x <- left_join(x = data.frame(Gene = my_cats[[cat]]), y = de_donor, by = "Gene")
    #   plot_de_summary(
    #     pb_donor_meta,
    #     pb_donor,
    #     de_donor,
    #     de,
    #     x$ensembl_id,
    #     glue("{out_dir}/cat/heatmap_{safe(cat)}.pdf"),
    #     title = cat
    #   )
    # }

    for (my_cluster in unique(de$cluster)) {
      n <- 25
      my_ens <- (
        de %>% filter(cluster == my_cluster) %>%
          filter(adj.P.Val < 0.05) %>%
          mutate(wt = abs(logFC) * -log10(P.Value)) %>%
          top_n(wt = wt, n = n)  %>%
          arrange(-logFC)
      )$ensembl_id
      if (length(my_ens) > 1) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          my_ens,
          glue("{out_dir}/cat/heatmap_top_genes_cluster{my_cluster}.pdf"),
          title = glue("Top {n} genes in cluster {my_cluster}")
        )
      }
    }

  }

  if (analysis_name == "a12_4_4_m3_2") {

    filename <- "data/colitis/M_DGE_lists.xlsx"
    sheets <- readxl::excel_sheets(filename)
    for (sheet in sheets) {
      message(sheet)
      # x <- rbindlist(lapply(sheets, function(sheet) {
      x <- readxl::read_excel(filename, sheet = sheet)[,1:2]
      colnames(x) <- c("gene", "cat")
      # })) %>% unique()
      x <- left_join(
        x = x,
        y = de_donor %>% select(Gene, ensembl_id, P.Value, logFC, CI.L, CI.R),
        by = c("gene" = "Gene")
      )
      x$cat <- str_replace(x$cat, "&", "and")
      x$cat <- str_split_fixed(x$cat, ",", 2)[,1]
      x <- x %>% group_by(cat) %>%
        mutate(wt = abs(logFC) * -log10(P.Value)) %>%
        ungroup() %>%
        arrange(cat, -wt)
      # cat_ids <- split(x$ensembl_id, x$cat)
      #
      # for (cat_name in names(cat_ids)) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          # cat_ids[[cat_name]],
          x$ensembl_id,
          glue("{out_dir}/cat/heatmap_{safe(sheet)}.pdf"),
          # title = cat_name
          title = sheet
        )
      # }
    }

    # my_cats <- list(
    #   "coinhibitory" = str_split("MYO7A CTLA4 MT1E LAIR2 HAVCR2 PDCD1 MT2A TOX2 CD160 CD38 BTLA LAG3 TOX LAIR1 ENTPD1 CD244 TIGIT", " ")[[1]],
    #   "costimulatory" = str_split("TNFRSF8 TNFRSF4 CD27 TNFRSF25 TNFRSF14 CD28 SLAMF1 CD40LG TNFRSF9 TNFSF14 CD2 CD226 ICOS TNFRSF18", " ")[[1]]
    # )
    # for (cat in names(my_cats)) {
    #   x <- left_join(x = data.frame(Gene = my_cats[[cat]]), y = de_donor, by = "Gene")
    #   plot_de_summary(
    #     pb_donor_meta,
    #     pb_donor,
    #     de_donor,
    #     de,
    #     x$ensembl_id,
    #     glue("{out_dir}/cat/heatmap_{safe(cat)}.pdf"),
    #     title = cat
    #   )
    # }

    for (my_cluster in unique(de$cluster)) {
      n <- 25
      my_ens <- (
        de %>% filter(cluster == my_cluster) %>%
          filter(adj.P.Val < 0.05) %>%
          mutate(wt = abs(logFC) * -log10(P.Value)) %>%
          top_n(wt = wt, n = n)  %>%
          arrange(-logFC)
      )$ensembl_id
      n <- length(my_ens)
      if (n > 1) {
        plot_de_summary(
          pb_donor_meta,
          pb_donor,
          de_donor,
          de,
          my_ens,
          glue("{out_dir}/cat/heatmap_top_genes_cluster{my_cluster}.pdf"),
          title = glue("Top {n} genes in cluster {my_cluster}")
        )
      }
    }

  }

  if (analysis_name == "n3_2") {

    # my_cats <- list(
    #   "isg" = str_split("TMEM173 CGAS ISG15 OAS2 GBP4", " ")[[1]],
    #   "apoptosis" = str_split("TNFRSF10A CASP8 ZBP1 ERBB4", " ")[[1]],
    #   "metabolism" = str_split("IDO1 PCSK9", " ")[[1]],
    #   "antimicrobial" = str_split("DUOX2 NOS2", " ")[[1]],
    #   "solute" = str_split("CA4 AQP8 AQP7 SLC22A5 SLC22A23 SLC2A13", " ")[[1]],
    #   "neuron" = str_split("ROBO1 NRXN3", " ")[[1]],
    #   "tsignal" = str_split("CD274 NECTIN2 COL4A1 HLA-G HLA-DRA HHLA2", " ")[[1]],
    #   "cytokine" = str_split("STAT1 JAK3 IL15 IL15RA LIFR CR1", " ")[[1]],
    #   "adhesion" = str_split("CXCL11 CXCL1 ICAM1 CEACAM1 CD44 CXCL13 ZEB1 CEACAM7", " ")[[1]]
    # )

    my_cats <- list(
      "cytokine" = str_split("STAT1 STAT2 STAT3 STAT4 JAK2 JAK3 IL18 IL1A IL1B IL15 IL15RA IL20RB IL31RA OSMR CR1", " ")[[1]],
      "adhesion" = str_split("CXCL1 CXCL2 CXCL3 CXCL8 CXCL9 CXCL10 CXCL11 CXCL13 CEACAM1 CEACAM7 CEACAM6 ICAM1 ICAM2 ALCAM MCAM PECAM1 CD44 ZEB1", " ")[[1]],
      "signal" = str_split("CD274 NECTIN2 COL4A1 HLA-E HLA-G HLA-DRA HHLA2", " ")[[1]],
      "isg" = str_split("TMEM173 CGAS ISG15 OAS2 GBP4", " ")[[1]],
      "apoptosis" = str_split("TNFRSF10A CASP1 CASP8 ZBP1 ERBB4", " ")[[1]],
      "antimicrobial" = str_split("DUOXA2 DUOX2 NOS2", " ")[[1]],
      "metabolism" = str_split("IDO1 PCSK9 PCK1 PARP14 ASS1", " ")[[1]],
      "solute" = str_split("AQP7 AQP8 AQP11 CA4 SLC2A13 SLC22A5 SLC22A23", " ")[[1]],
      "neuron" = str_split("ROBO1 NRXN3 LSAMP", " ")[[1]],
      "hypoxia" = str_split("HIF1A SERPINE1 CXCL12 PTGS2 VEGFA", " ")[[1]]
    )
    my_cats_d <- rbindlist(lapply(names(my_cats), function(my_cat) {
      data.frame(cat = my_cat, Gene = my_cats[[my_cat]])
    }))
    my_cats_d$cat <- factor(my_cats_d$cat, names(my_cats))
    x <- left_join(x = my_cats_d, y = de_donor, by = "Gene")
    x <- x %>%
      mutate(wt = abs(logFC) * -log10(P.Value)) %>%
      arrange(cat, -wt)
    de_x <- left_join(x = my_cats_d, y = de, by = "Gene")
    de_donor_x <- left_join(x =  my_cats_d, y = de_donor, by = "Gene")
    cluster_groups <- get_cluster_groups(analysis_name)
    if (length(cluster_groups) == length(unique(de_x$cluster))) {
      de_x$cluster_group <- cluster_groups[as.character(de_x$cluster)]
    }
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      x$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats1.pdf"),
      cluster_order = rev(cluster_order),
      title = "cats"
    )

    #
    plot_de_summary_h(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      x$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats1_h.pdf"),
      cluster_order = cluster_order,
      title = "cats"
    )

    my_cats <- list(
      "g1 isg"   = "ICAM1 CXCL11 CXCL10 STAT1 GBP4 ISG15 TMEM173",
      "g2 mhc"   = "HLA-DRA CD274 HLA-G",
      "g3 antim" = "DUOX2 NOS2",
      "g4 water" = "AQP8 AQP7",
      "g5 apop" = "ZBP1 CASP1"
    )
    my_cats_d <- rbindlist(lapply(names(my_cats), function(my_cat) {
      data.frame(
        cat = my_cat,
        Gene = str_split(my_cats[[my_cat]], " ")[[1]]
      )
    }))
    # my_cats_d$Gene %intersect% ensembl_to_symbol
    my_cats_d$cat <- factor(my_cats_d$cat, names(my_cats))
    x <- left_join(x = my_cats_d, y = de_donor, by = "Gene")
    x <- x %>%
      mutate(wt = abs(logFC) * -log10(P.Value)) %>%
      arrange(cat, -wt)
    de_x <- left_join(x = my_cats_d, y = de, by = "Gene")
    de_donor_x <- left_join(x =  my_cats_d, y = de_donor, by = "Gene")
    cluster_groups <- get_cluster_groups(analysis_name)
    if (length(cluster_groups) == length(unique(de_x$cluster))) {
      de_x$cluster_group <- cluster_groups[as.character(de_x$cluster)]
    }
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      x$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats3.pdf"),
      cluster_order = rev(cluster_order),
      title = "cats"
    )

    for (cat in names(my_cats)) {
      x <- left_join(x = data.frame(Gene = strsplit(my_cats[[cat]], " ")[[1]]), y = de_donor, by = "Gene")
      plot_de_summary(
        pb_donor_meta,
        pb_donor,
        de_donor,
        de,
        x$ensembl_id,
        glue("{out_dir}/cat/heatmap_{safe(cat)}.pdf"),
        title = cat
      )
    }

    my_cats <- list(
      "isg" = str_split("TMEM173 CGAS ISG15 OAS2 GBP4", " ")[[1]],
      "apoptosis" = str_split("TNFRSF10A CASP8 ZBP1 ERBB4", " ")[[1]],
      "antimicrobial" = str_split("DUOXA2 DUOX2 NOS2", " ")[[1]],
      "metabolism" = str_split("IDO1 PCSK9 PCK1 PARP14 ASS1", " ")[[1]],
      "solute" = str_split("AQP7 AQP8 AQP11 CA4 SLC2A13 SLC22A5 SLC22A23", " ")[[1]],
      "neuron" = str_split("ROBO1 NRXN3 LSAMP", " ")[[1]],
      "signal" = str_split("CD274 NECTIN2 COL4A1 HLA-G HLA-DRA HHLA2", " ")[[1]]
    )
    x <- left_join(
      x =  rbindlist(lapply(names(my_cats), function(my_cat) {
        data.frame(cat = my_cat, Gene = my_cats[[my_cat]])
      })),
      y = de_donor,
      by = "Gene"
    )
    x <- x %>%
      mutate(wt = abs(logFC) * -log10(P.Value)) %>%
      arrange(cat, -wt)
    de_x <- left_join(
      x =  rbindlist(lapply(names(my_cats), function(my_cat) {
        data.frame(cat = my_cat, Gene = my_cats[[my_cat]])
      })),
      y = de,
      by = "Gene"
    )
    de_donor_x <- left_join(
      x =  rbindlist(lapply(names(my_cats), function(my_cat) {
        data.frame(cat = my_cat, Gene = my_cats[[my_cat]])
      })),
      y = de_donor,
      by = "Gene"
    )
    plot_de_summary(
      pb_donor_meta,
      pb_donor,
      de_donor_x,
      de_x,
      ensembl_ids = x$ensembl_id,
      out_file = glue("{out_dir}/cat/heatmap_cats2.pdf"),
      cluster_order = cluster_order,
      title = "cats"
    )
    #
    plot_de_summary_h(
      pb_donor_meta,
      pb_donor,
      de_donor,
      de_x,
      x$ensembl_id,
      glue("{out_dir}/cat/heatmap_cats2_h.pdf"),
      cluster_order = cluster_order,
      title = "cats"
    )

    for (cat in names(my_cats)) {
      x <- left_join(x = data.frame(Gene = my_cats[[cat]]), y = de_donor, by = "Gene")
      plot_de_summary(
        pb_donor_meta,
        pb_donor,
        de_donor,
        de,
        x$ensembl_id,
        glue("{out_dir}/cat/heatmap_{safe(cat)}.pdf"),
        title = cat
      )
    }

    a1$log2cpm <- do_log2cpm(a1$counts)

    my_genes <- c(
      "STAT1", "CXCL1", "CXCL11", "ICAM1", "GBP4", "CASP8", "DUOX2", "PCSK9",
      "AQP8", "ROBO1", "CD274", "HLA-DRA"
    )
    for (my_gene in my_genes) {
      my_id <- names(which(ensembl_to_symbol == my_gene))
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$log2cpm[my_id,]),
        group = append_n(a1$obs$case),
        bins = umap_n_bins,
        palette = "oslo"
      ) +
      facet_wrap(~ group) +
      labs(title = ensembl_to_symbol[my_id]) +
      theme(
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-{safe(ensembl_to_symbol[my_id])}"),
        out_dir = file.path(out_dir, "umap/de_cluster/slim"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 6.5,
        height = 3,
        units = "in",
        dpi = 300
      )
    }

  }

}

# }}}

a1_file <- "paper/pseudobulk_donor.qs"
if (file.exists(a1_file)) {
  a1 <- qread(a1_file)
} else {
  a1 <- lapply(analyses, function(analysis_name) {
    out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
    a1 <- qread(file.path(out_dir, "pseudobulk_donor.qs"))
    a1$meta$analysis <- analysis_name
    a1
  })
  common_genes <- lapply(a1, function(x) rownames(x$log2cpm))
  # common_genes <- Reduce(intersect, common_genes)
  common_genes <- Reduce(union, common_genes)
  a1 <- list(
    meta = rbindlist(lapply(a1, function(x) x$meta)),
    # log2cpm = do.call(cbind, lapply(a1, function(x) x$log2cpm[common_genes,]))
    log2cpm = do.call(cbind, lapply(a1, function(x) select_rows(x$log2cpm, common_genes)))
  )
  stopifnot(all(a1$meta$donor == colnames(a1$log2cpm)))
  analysis_name_pub <- c(
    "a12_4_4_t4_cd8_1_2"         = "Tissue CD8 T cells",
    "luoma_cd45_a5_tcell2_cd8_3" = "Luoma Tissue CD8 T cells",
    "a12_4_4_t4_cd4_2_2"         = "Tissue CD4 T cells",
    "luoma_cd45_a5_tcell2_cd4_3" = "Luoma Tissue CD4 T cells",
    "a12_4_4_m3_2"               = "Tissue Myeloid cells",
    "a12_4_4_b5_1_3"             = "Tissue B cells",
    "n3_2"                       = "Tissue epithelial nuclei",
    "blood2_tcell5_cd8_5"        = "Blood CD8 T cells",
    "blood2_tcell5_cd4_5"        = "Blood CD4 T cells",
    "blood2_myeloid5"            = "Blood Myeloid cells",
    "blood2_bcell5"              = "Blood B cells"
  )
  a1$meta$analysis_pub <- factor(
    analysis_name_pub[a1$meta$analysis], levels = rev(unname(analysis_name_pub))
  )
  a1$meta$tissue <- ifelse(str_detect(a1$meta$analysis, "blood"), "Blood", "Tissue")
  #
  # sample_info <- janitor::clean_names(read_excel(
  #     path = "data/luoma_villani_combined_colitis_samples2.xlsx",
  #     sheet = 3
  # ))
  # a1$meta <- left_join(
  #   x = a1$meta,
  #   y = sample_info %>% select(donor, pathway) %>% rename(drug = pathway),
  #   by = c("donor")
  # )
  a1$meta <- a1$meta %>% mutate(id = glue("{analysis_pub}|{donor}")) %>%
    select(-analysis) %>%
    rename(dataset = analysis_pub)
  colnames(a1$log2cpm) <- a1$meta$id
  head(a1$meta)
  a1$log2cpm[1:3,1:3]
  #
  qsave(a1, "paper/pseudobulk_donor.qs")
}


# Write files for sharing
fwrite(a1$meta, "paper/pseudobulk_donor-meta.tsv", sep = "\t")
writeMMgz(a1$log2cpm, "paper/pseudobulk_donor.mtx.gz")
barcodes_file <- "paper/pseudobulk_donor-barcodes.tsv"
writeLines(
  text = colnames(a1$log2cpm),
  con = barcodes_file
)
genes_file <- "paper/pseudobulk_donor-genes.tsv"
writeLines(
  # text = ensembl_to_symbol[rownames(a1$log2cpm)],
  text = rownames(a1$log2cpm),
  con = genes_file
)

de_donor <- fread("paper/de-case-vs-control.tsv.gz")
de_donor$analysis_pub <- factor(
  as.character(de_donor$analysis), levels = rev(unname(analysis_name_pub))
)
de_donor$tissue <- ifelse(str_detect(de_donor$analysis, "Blood"), "Blood", "Tissue")

# Figure 6
# my_genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "CD86", "CD80")
my_genes <- c("PDCD1", "CD274", "PDCD1LG2")
my_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in% my_genes])
for (my_id in my_ids) {
  a1$meta[[my_id]] <- as.numeric(a1$log2cpm[my_id,])
}
p <- ggplot(
  a1$meta %>% pivot_longer(starts_with("ENSG")) %>%
    filter(!str_detect(dataset, "Luoma")) %>%
    mutate(
      symbol = factor(ensembl_to_symbol[name], my_genes)
    )
) +
  aes(x = value, y = dataset, group = case, fill = case) +
  ggforestplot::geom_stripes() +
  # geom_point(
  #   shape = 21, size = 3, stroke = 0.2,
  #   position = position_quasirandom(groupOnX = FALSE, width = 0.1, dodge.width = 0.5)
  # ) +
  geom_boxplot(
    aes(x = value, y = dataset, group = paste(dataset, case), fill = case),
    size = 0.2
  ) +
  # stat_summary(shape = 21, stroke = 0.2, size = 0.8, fun.data = mean_se) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  # facet_row(vars(symbol), scales = "free_x") +
  facet_grid(tissue ~ symbol, scales = "free", space = "free_y") +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = pals::okabe(), name = NULL) +
  labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
  theme(
    strip.text      = element_text(face = "italic"),
    panel.spacing   = unit(1, "lines"),
    legend.position = "bottom"
  )
my_ggsave(
  "global-dots",
  out_dir = "results/a20",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 4 + 1.5 * length(my_ids),
  height = 4.5,
  units = "in",
  dpi = 300
)

for (my_id in my_ids) {
  # my_id <- my_ids[1]
  #
  d1 <- de_donor %>%
    filter(!str_detect(analysis, "Luoma"), cluster == "all cells", ensembl_id == my_id) %>%
    mutate(dataset = analysis)
  d1$dataset[d1$dataset == "Tissue epithelial cells"] <- "Tissue epithelial nuclei"
  #
  d2 <- a1$meta %>%
    filter(!str_detect(dataset, "Luoma")) %>%
    select(dataset, case, id, starts_with("ENSG")) %>%
    pivot_longer(starts_with("ENSG")) %>%
    filter(name == my_id) %>%
    mutate(tissue = ifelse(str_detect(dataset, "Blood"), "Blood", "Tissue"))
  #
  all_analyses <- naturalsort(
    union(as.character(unique(d1$analysis)), as.character(unique(a1$meta$dataset)))
  )
  d1$dataset <- factor(as.character(d1$dataset), all_analyses)
  d2$dataset<- factor(as.character(d2$dataset), all_analyses)
  d1 <- full_join(d1, d2 %>% select(dataset) %>% unique, by = c("dataset"))
  d1 <- d1 %>% filter(!is.na(dataset))
  d1$tissue = ifelse(str_detect(d1$dataset, "Blood"), "Blood", "Tissue")
  #
  p_error <- ggplot(d1) +
    aes(x = log_fc, xmin = ci_l, xmax = ci_r, y = dataset, color = adj_p_val < 0.05) +
    ggforestplot::geom_stripes(aes(xmin = -Inf, xmax = Inf)) +
    scale_y_discrete(expand = c(0, 0.5)) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_point() +
    geom_errorbarh(height = 0) +
    scale_color_manual(values = c("grey60", "black"), guide = "none") +
    facet_grid(rows = vars(tissue), scales = "free", space = "free_y") +
    scale_x_continuous(
      breaks = scales::pretty_breaks(3),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    labs(x = "FC", y = NULL, title = ensembl_to_symbol[my_id]) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "italic"),
      strip.text = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  #
  p_donor <- ggplot(d2) +
      ggforestplot::geom_stripes(aes(y = dataset)) +
      scale_y_discrete(expand = c(0, 0.5), position = "right") +
      # geom_point(
      #   mapping = aes(x = value, y = dataset, group = case, fill = case),
      #   position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE),
      #   shape = 21, size = 2.5, stroke = 0.2
      # ) +
      geom_boxploth(
        mapping = aes(x = value, y = dataset, group = paste(dataset, case), fill = case),
        size = 0.3
      ) +
      scale_fill_manual(values = pals::okabe(8)) +
      facet_grid(rows = vars(tissue), scales = "free", space = "free_y") +
      scale_x_continuous(
        breaks = scales::pretty_breaks(3)
      ) +
      # guides(fill = guide_legend(reverse = TRUE, title = NULL, override.aes = list(size = 5))) +
      labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
      theme(
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)
      )
  layout <- "AB"
  p <- p_error + p_donor + plot_layout(
    design = layout,
    widths = c(1, 1),
    heights = c(1, 1),
  )
  my_ggsave(
    glue("global-errorbars-{ensembl_to_symbol[my_id]}"),
    out_dir = "results/a20/global",
    type = "pdf",
    plot = p,
    scale = 1,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300
  )
}


analysis_name_pub <- c(
  "a12_4_4_t4_cd8_1_2"         = "Tissue CD8 T cells",
  "luoma_cd45_a5_tcell2_cd8_3" = "Luoma Tissue CD8 T cells",
  "a12_4_4_t4_cd4_2_2"         = "Tissue CD4 T cells",
  "luoma_cd45_a5_tcell2_cd4_3" = "Luoma Tissue CD4 T cells",
  "a12_4_4_m3_2"               = "Tissue Myeloid cells",
  "a12_4_4_b5_1_3"             = "Tissue B cells",
  "n3_2"                       = "Tissue epithelial nuclei",
  "blood2_tcell5_cd8_5"        = "Blood CD8 T cells",
  "blood2_tcell5_cd4_5"        = "Blood CD4 T cells",
  "blood2_myeloid5"            = "Blood Myeloid cells",
  "blood2_bcell5"              = "Blood B cells"
)


a1_file <- "paper/pseudobulk_donor_cluster.qs"
if (file.exists(a1_file)) {
  a1 <- qread(a1_file)
} else {
  a1 <- lapply(analyses, function(analysis_name) {
    out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
    a1 <- qread(file.path(out_dir, "pseudobulk_donor_cluster.qs"))
    a1$meta$analysis <- analysis_name
    a1
  })
  common_genes <- lapply(a1, function(x) rownames(x$log2cpm))
  # common_genes <- Reduce(intersect, common_genes)
  common_genes <- Reduce(union, common_genes)
  a1 <- list(
    meta = rbindlist(lapply(a1, function(x) x$meta)),
    # log2cpm = do.call(cbind, lapply(a1, function(x) x$log2cpm[common_genes,]))
    log2cpm = do.call(cbind, lapply(a1, function(x) select_rows(x$log2cpm, common_genes)))
  )
  colnames(a1$log2cpm) <- str_replace(colnames(a1$log2cpm), "factor\\(cluster\\)", "")
  colnames(a1$log2cpm) <- str_replace(colnames(a1$log2cpm), "factor\\(donor\\)", "")
  colnames(a1$log2cpm) <- str_replace(colnames(a1$log2cpm), ":", "|")
  stopifnot(all(sprintf("%s|%s", a1$meta$cluster, a1$meta$donor) == colnames(a1$log2cpm)))
  a1$meta$analysis_pub <- factor(
    analysis_name_pub[a1$meta$analysis], levels = rev(unname(analysis_name_pub))
  )
  a1$meta$tissue <- ifelse(str_detect(a1$meta$analysis, "blood"), "Blood", "Tissue")
  #
  # sample_info <- janitor::clean_names(read_excel(
  #     path = "data/luoma_villani_combined_colitis_samples2.xlsx",
  #     sheet = 3
  # ))
  # a1$meta <- left_join(
  #   x = a1$meta,
  #   y = sample_info %>% select(donor, pathway) %>% rename(drug = pathway),
  #   by = c("donor")
  # )
  a1$meta <- a1$meta %>% mutate(id = glue("{analysis_pub}|{cluster}|{donor}")) %>%
    select(-analysis) %>%
    rename(dataset = analysis_pub)
  colnames(a1$log2cpm) <- a1$meta$id
  head(a1$meta)
  a1$log2cpm[1:3,1:3]
  #
  qsave(a1, "paper/pseudobulk_donor_cluster.qs")
}

# Write files for sharing
fwrite(a1$meta, "paper/pseudobulk_donor_cluster-meta.tsv", sep = "\t")
writeMMgz(a1$log2cpm, "paper/pseudobulk_donor_cluster.mtx.gz")
barcodes_file <- "paper/pseudobulk_donor_cluster-barcodes.tsv"
writeLines(
  text = colnames(a1$log2cpm),
  con = barcodes_file
)
genes_file <- "paper/pseudobulk_donor_cluster-genes.tsv"
writeLines(
  # text = ensembl_to_symbol[rownames(a1$log2cpm)],
  text = rownames(a1$log2cpm),
  con = genes_file
)

# Figure 7
cluster_order <-  readRDS("results/a20/n3_2/figures/composition-case-vs-control/cluster_order.rds")
my_genes <- c("S1PR1", "CX3CL1", "HIF1A")
my_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in% my_genes])
a1$meta %<>% select(!starts_with("ENSG"))
for (my_id in my_ids) {
  a1$meta[[my_id]] <- as.numeric(a1$log2cpm[my_id,])
}
my_d <- a1$meta %>%
  pivot_longer(starts_with("ENSG")) %>%
  filter(dataset == "Tissue epithelial nuclei") %>%
  mutate(symbol = factor(ensembl_to_symbol[name], my_genes))
# my_clusters <- (
#   my_d %>%
#   filter(symbol == "CX3CL1") %>%
#   group_by(cluster) %>%
#   summarize(mean = mean(value)) %>%
#   arrange(mean)
# )$cluster
my_d <- my_d %>% mutate(cluster = factor(cluster, cluster_order))
cluster_groups <- get_cluster_groups("n3_2")
my_d$cluster_group <- cluster_groups[as.character(my_d$cluster)]
p <- ggplot(my_d) +
  aes(x = value, y = cluster, group = case, fill = case) +
  ggforestplot::geom_stripes() +
  # geom_point(
  #   shape = 21, size = 3, stroke = 0.2,
  #   position = position_quasirandom(groupOnX = FALSE, width = 0.1, dodge.width = 0.5)
  # ) +
  geom_boxplot(
    aes(x = value, y = cluster, group = paste(cluster, case), fill = case),
    size = 0.2, outlier.size = 0.1
  ) +
  # stat_summary(shape = 21, stroke = 0.2, size = 0.8, fun.data = mean_se) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  # facet_row(vars(symbol), scales = "free_x") +
  # facet_grid(tissue ~ symbol, scales = "free", space = "free_y") +
  facet_grid(rows = vars(cluster_group), cols = vars(symbol), scales = "free", space = "free_y") +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = pals::okabe(), name = NULL) +
  labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
  theme(
    strip.placement = "outside",
    strip.text.x    = element_text(face = "italic"),
    strip.text.y    = element_text(hjust = 0, angle = 0),
    panel.spacing   = unit(0.5, "lines"),
    legend.position = "none"
  )
my_ggsave(
  "boxplot-case-S1PR1-CX3CL1-HIF1A",
  out_dir = "results/a20/n3_2/figures/boxplot",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 4 + 1 * length(my_ids),
  height = 6.5,
  units = "in",
  dpi = 300
)
#
cluster_colors <- mpn65[seq(length(unique(my_d$cluster)))]
names(cluster_colors) <- seq(length(unique(my_d$cluster)))
p <- ggplot(my_d) +
  aes(x = value, y = cluster) +
  ggforestplot::geom_stripes() +
  geom_boxplot(
    aes(x = value, y = cluster, fill = cluster), size = 0.2, outlier.size = 0.1
  ) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  facet_grid(rows = vars(cluster_group), cols = vars(symbol), scales = "free", space = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  scale_y_discrete(position = "right") +
  labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
  theme(
    strip.placement = "outside",
    strip.text.x    = element_text(face = "italic"),
    strip.text.y    = element_text(hjust = 0, angle = 0),
    panel.spacing   = unit(0.5, "lines"),
    legend.position = "none"
  )
my_ggsave(
  "boxplot-S1PR1-CX3CL1-HIF1A",
  out_dir = "results/a20/n3_2/figures/boxplot",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 4 + 1 * length(my_ids),
  height = 6,
  units = "in",
  dpi = 300
)

# Figure 7
cluster_order <-  readRDS("results/a20/n3_2/figures/composition-case-vs-control/cluster_order.rds")
my_genes <- c("SPINK5")
my_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in% my_genes])
a1$meta %<>% select(!starts_with("ENSG"))
for (my_id in my_ids) {
  a1$meta[[my_id]] <- as.numeric(a1$log2cpm[my_id,])
}
my_d <- a1$meta %>%
  pivot_longer(starts_with("ENSG")) %>%
  filter(dataset == "Tissue epithelial nuclei") %>%
  mutate(symbol = factor(ensembl_to_symbol[name], my_genes))
# my_clusters <- (
#   my_d %>%
#   filter(symbol == "CX3CL1") %>%
#   group_by(cluster) %>%
#   summarize(mean = mean(value)) %>%
#   arrange(mean)
# )$cluster
my_d <- my_d %>% mutate(cluster = factor(cluster, cluster_order))
cluster_groups <- get_cluster_groups("n3_2")
my_d$cluster_group <- cluster_groups[as.character(my_d$cluster)]
p <- ggplot(my_d) +
  aes(x = value, y = cluster, group = case, fill = case) +
  ggforestplot::geom_stripes() +
  # geom_point(
  #   shape = 21, size = 3, stroke = 0.2,
  #   position = position_quasirandom(groupOnX = FALSE, width = 0.1, dodge.width = 0.5)
  # ) +
  geom_boxplot(
    aes(x = value, y = cluster, group = paste(cluster, case), fill = case),
    size = 0.2, outlier.size = 0.1
  ) +
  # stat_summary(shape = 21, stroke = 0.2, size = 0.8, fun.data = mean_se) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  # facet_row(vars(symbol), scales = "free_x") +
  # facet_grid(tissue ~ symbol, scales = "free", space = "free_y") +
  facet_grid(rows = vars(cluster_group), cols = vars(symbol), scales = "free", space = "free_y") +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = pals::okabe(), name = NULL) +
  labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
  theme(
    strip.placement = "outside",
    strip.text.x    = element_text(face = "italic"),
    strip.text.y    = element_text(hjust = 0, angle = 0),
    panel.spacing   = unit(0.5, "lines"),
    legend.position = "none"
  )
my_ggsave(
  "boxplot-case-S1PR1-CX3CL1-HIF1A",
  out_dir = "results/a20/n3_2/figures/boxplot",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 4 + 1 * length(my_ids),
  height = 6.5,
  units = "in",
  dpi = 300
)
#
cluster_colors <- mpn65[seq(length(unique(my_d$cluster)))]
names(cluster_colors) <- seq(length(unique(my_d$cluster)))
p <- ggplot(my_d) +
  aes(x = value, y = cluster) +
  ggforestplot::geom_stripes() +
  geom_boxplot(
    aes(x = value, y = cluster, fill = cluster), size = 0.2, outlier.size = 0.1
  ) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  facet_grid(rows = vars(cluster_group), cols = vars(symbol), scales = "free", space = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  scale_y_discrete(position = "right") +
  labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
  theme(
    strip.placement = "outside",
    strip.text.x    = element_text(face = "italic"),
    strip.text.y    = element_text(hjust = 0, angle = 0),
    panel.spacing   = unit(0.5, "lines"),
    legend.position = "none"
  )
my_ggsave(
  "boxplot-S1PR1-CX3CL1-HIF1A",
  out_dir = "results/a20/n3_2/figures/boxplot",
  type = "pdf",
  plot = p,
  scale = 1,
  width = 4 + 1 * length(my_ids),
  height = 6,
  units = "in",
  dpi = 300
)

# Differentially expressed genes in each cell cluster
########################################################################

analyses <- c(
  "a12_4_4_t4_cd8_1_2",
  "a12_4_4_t4_cd4_2_2",
  "a12_4_4_m3_2",
  "a12_4_4_b5_1_3",
  "n3_2"
  # "blood2_myeloid5",
  # "blood2_bcell5",
  # "blood2_tcell5_cd4_5"
  # "blood2_tcell5_cd8_5"
)
de <- rbindlist(lapply(analyses, function(analysis_name) {
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
  de_donor_file <- glue("{out_dir}/de_case-vs-control.tsv.gz")
  de <- fread(de_donor_file)
  de$celltype <- analysis_name
  return(de)
}))
de$celltype <- str_remove(de$celltype, "a12_4_4_")
de <- de %>% filter(celltype != "blood2")
to_celltype <- c("B", "M", "E", "CD4", "CD8", "Blood")
names(to_celltype) <- c("b5_1_3", "m3_2", "n3_2", "t4_cd4_2_2", "t4_cd8_1_2", "blood2")
de$celltype <- to_celltype[de$celltype]
#
out_dir <- "results/a20/de-case-vs-control"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
de$cluster <- sprintf("%s-%s", de$celltype, de$cluster)
#
de_top_genes <- unique(c(
  (
    de %>% group_by(cluster) %>%
    top_n(
      wt = (abs(logFC) * -log10(P.Value)),
      n = 10
    ) %>%
    filter(abs(logFC) > log2(1.5) & adj.P.Val < 0.05)
  )$ensembl_id
))
de_top <- de %>% filter(ensembl_id %in% de_top_genes)

{
  # "IL17A" %in% de_top$Gene
  #
  mat <- dcast(data = de_top, formula = ensembl_id ~ cluster, value.var = "logFC")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  mat_fdr <- dcast(data = de_top, formula = ensembl_id ~ cluster, value.var = "adj.P.Val")
  mat_fdr_rows <- mat_fdr[[1]]
  mat_fdr[[1]] <- NULL
  mat_fdr <- as.matrix(mat_fdr)
  mat_fdr[is.na(mat_fdr)] <- 0
  rownames(mat_fdr) <- ensembl_to_symbol[mat_fdr_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 1)
  )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    top_annotation = column_ha,
    right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    # cell_fun = function(j, i, x, y, width, height, fill) {
    #   if (mat_fdr[i, j] < 0.05) {
    #     grid.circle(
    #       x = x, y = y, r = unit("0.33", "mm"),
    #       gp = gpar(fill = "white", col = "black", lwd = 0.1)
    #     )
    #   }
    # },
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-top-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.05 + 2)
  draw(ht)
  dev.off()
}

{
  my_pairs <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  my_pairs <- pivot_longer(my_pairs, cols = c("g_target", "g_source"))
  my_genes <- unique(my_pairs$value)
  gene_to_cat <- with(
    my_pairs %>% select(cat, value) %>% unique(),
    split(cat, value)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  #
  de_select <- de %>% filter(Gene %in% my_genes)
  # de_select <- de_select %>%
  #   filter(
  #     Gene %in% (
  #       de_select %>% group_by(Gene) %>% summarize(min_fdr = min(adj.P.Val)) %>% filter(min_fdr < 0.05)
  #     )$Gene
  #   )
  # "IL17A" %in% de_select$Gene
  #
  mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "logFC")
  # mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(cluster) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    border_gp = gpar(col = "black", lty = 1),
    row_split = gene_to_cat[rownames(mat)],
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-select-log2fc-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.15 + 2)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(cluster) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(min(mat), max(mat), length.out = 9),
    RColorBrewer::brewer.pal(name = "BuPu", n = 9)
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Mean log2CPM",
    top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-select-mean-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.15 + 2)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  mat <- t(scale(t(mat)))
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(cluster) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "z",
    top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-select-z-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.15 + 2)
  draw(ht)
  dev.off()
}

{
  my_drug_targets <- readxl::read_excel("data/drug-targets.xlsx")
  my_genes <- unique(my_drug_targets$Target)
  gene_to_cat <- with(
    my_drug_targets %>% select(Drug, Target) %>% unique(),
    split(Drug, Target)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  #
  de_select <- de %>% filter(Gene %in% my_genes)
  # de_select <- de_select %>%
  #   filter(
  #     Gene %in% (
  #       de_select %>% group_by(Gene) %>% summarize(min_fdr = min(adj.P.Val)) %>% filter(min_fdr < 0.05)
  #     )$Gene
  #   )
  # "IL17A" %in% de_select$Gene
  #
  mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "logFC")
  # mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(cluster) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    border_gp = gpar(col = "black", lty = 1),
    # row_split = gene_to_cat[rownames(mat)],
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-drugtargets-log2fc-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.15 + 2)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(cluster) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(min(mat), max(mat), length.out = 9),
    RColorBrewer::brewer.pal(name = "BuPu", n = 9)
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Mean log2CPM",
    top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-drugtargets-mean-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.15 + 2)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ cluster, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  mat <- t(scale(t(mat)))
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(cluster) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(cluster) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  col_colors <- mpn65[1:length(unique(de$celltype))]
  names(col_colors) <- as.character(naturalsort(unique(de$celltype)))
  # names(col_colors) <- as.character(seq_along(col_colors))
  column_ha <- HeatmapAnnotation(
    cluster = sapply(strsplit(colnames(mat), "-", 2), "[", 1),
    col = list(cluster = col_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "z",
    top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    column_order = order(uwot::umap(X = t(mat), n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    # column_names_rot = 0,
    column_names_gp = gpar(angle = 60, fontsize = 8, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-drugtargets-z-cluster.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.12 + 2, height = nrow(mat) * 0.15 + 2)
  draw(ht)
  dev.off()
}

# Differentially expressed genes in each cell type (E, CD8, CD4, M, B)
########################################################################

de <- rbindlist(lapply(analyses, function(analysis_name) {
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
  de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
  de <- fread(de_donor_file)
  de$celltype <- analysis_name
  return(de)
}))
de$celltype <- str_remove(de$celltype, "a12_4_4_")
to_celltype <- c("B", "M", "E", "CD4", "CD8", "Blood")
names(to_celltype) <- c("b5_1_3", "m3_2", "n3_2", "t4_cd4_2_2", "t4_cd8_1_2", "blood2")
de$celltype <- to_celltype[de$celltype]
#
out_dir <- "results/a20/de-case-vs-control"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#
de_top_genes <- unique(c(
  (
    # de %>% group_by(celltype) %>%
    # filter(
    #   abs(logFC) > log2(2) & adj.P.Val < 0.01
    # ) %>% ungroup
    de %>% group_by(celltype) %>%
    top_n(
      wt = (abs(logFC) * -log10(P.Value)),
      n = 30
    ) %>%
    filter(abs(logFC) > log2(1.5) & adj.P.Val < 0.05)
  )$ensembl_id
))
de_top <- de %>% filter(ensembl_id %in% de_top_genes)
# "IL17A" %in% de_top$Gene

{
  mat <- dcast(data = de_top, formula = ensembl_id ~ celltype, value.var = "logFC")
  # mat <- dcast(data = de_top, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  mat_fdr <- dcast(data = de_top, formula = ensembl_id ~ celltype, value.var = "adj.P.Val")
  mat_fdr_rows <- mat_fdr[[1]]
  mat_fdr[[1]] <- NULL
  mat_fdr <- as.matrix(mat_fdr)
  mat_fdr[is.na(mat_fdr)] <- 0
  rownames(mat_fdr) <- ensembl_to_symbol[mat_fdr_rows]
  #
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  row_order <- order(uwot::umap(X = mat, n_components = 1)[,1])
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    col = col_fun,
    row_order = row_order,
    show_row_names = TRUE,
    row_names_gp = gpar(fontface = "italic"),
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-top-tall.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.5 + 2, height = nrow(mat) * 0.14 + 2)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_top, formula = ensembl_id ~ celltype, value.var = "logFC")
  # mat <- dcast(data = de_top, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  mat_fdr <- dcast(data = de_top, formula = ensembl_id ~ celltype, value.var = "adj.P.Val")
  mat_fdr_rows <- mat_fdr[[1]]
  mat_fdr[[1]] <- NULL
  mat_fdr <- as.matrix(mat_fdr)
  mat_fdr[is.na(mat_fdr)] <- 0
  rownames(mat_fdr) <- ensembl_to_symbol[mat_fdr_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_top %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  row_order <- order(uwot::umap(X = mat, n_components = 1)[,1])
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    # top_annotation = column_ha,
    right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    # cell_fun = function(j, i, x, y, width, height, fill) {
    #   if (mat_fdr[i, j] < 0.05) {
    #     grid.circle(
    #       x = x, y = y, r = unit("0.33", "mm"),
    #       gp = gpar(fill = "white", col = "black", lwd = 0.1)
    #     )
    #   }
    # },
    col = col_fun,
    # cluster_rows = FALSE,
    row_order = row_order,
    # row_split = 6,
    show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-top.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = ncol(mat) * 0.5 + 2, height = nrow(mat) * 0.08 + 2)
  draw(ht)
  dev.off()
}

{
  # my_pairs <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  # my_genes <- unique(c(my_pairs$g_target, my_pairs$g_source))
  my_pairs <- readxl::read_excel("data/curated-ligand-receptor.xlsx")
  my_pairs <- pivot_longer(my_pairs, cols = c("g_target", "g_source"))
  my_genes <- unique(my_pairs$value)
  gene_to_cat <- with(
    my_pairs %>% select(cat, value) %>% unique(),
    split(cat, value)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  #
  de_select <- de %>% filter(Gene %in% my_genes)
  # de_select <- de_select %>%
  #   filter(
  #     Gene %in% (
  #       de_select %>% group_by(Gene) %>% summarize(min_fdr = min(adj.P.Val)) %>% filter(min_fdr < 0.05)
  #     )$Gene
  #   )
  # "IL17A" %in% de_select$Gene
  #
  mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "logFC")
  # mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    # top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    border_gp = gpar(col = "black", lty = 1),
    row_split = gene_to_cat[rownames(mat)],
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-select-log2fc.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = 5, height = nrow(mat) * 0.15 + 1)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(min(mat), max(mat), length.out = 9),
    RColorBrewer::brewer.pal(name = "BuPu", n = 9)
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Mean log2CPM",
    # top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-select-mean.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = 5, height = nrow(mat) * 0.15 + 1)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  mat <- t(scale(t(mat)))
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "z",
    # top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-select-z.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = 5, height = nrow(mat) * 0.15 + 1)
  draw(ht)
  dev.off()
}


{
  my_drug_targets <- readxl::read_excel("data/drug-targets.xlsx")
  my_genes <- unique(my_drug_targets$Target)
  gene_to_cat <- with(
    my_drug_targets %>% select(Drug, Target) %>% unique(),
    split(Drug, Target)
  )
  gene_to_cat <- unlist(lapply(gene_to_cat, "[", 1))
  #
  de_select <- de %>% filter(Gene %in% my_genes)
  # de_select <- de_select %>%
  #   filter(
  #     Gene %in% (
  #       de_select %>% group_by(Gene) %>% summarize(min_fdr = min(adj.P.Val)) %>% filter(min_fdr < 0.05)
  #     )$Gene
  #   )
  # "IL17A" %in% de_select$Gene
  #
  mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "logFC")
  # mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Log2FC",
    # top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    border_gp = gpar(col = "black", lty = 1),
    # row_split = gene_to_cat[rownames(mat)],
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-drugtargets-log2fc.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = 5, height = nrow(mat) * 0.15 + 1)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(min(mat), max(mat), length.out = 9),
    RColorBrewer::brewer.pal(name = "BuPu", n = 9)
  )
  ht <- Heatmap(
    matrix = mat,
    name = "Mean log2CPM",
    # top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-drugtargets-mean.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = 5, height = nrow(mat) * 0.15 + 1)
  draw(ht)
  dev.off()
}

{
  mat <- dcast(data = de_select, formula = ensembl_id ~ celltype, value.var = "AveExpr")
  mat_rows <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  rownames(mat) <- ensembl_to_symbol[mat_rows]
  mat <- t(scale(t(mat)))
  #
  # my_genes <- c("IL17A", "CXCL13", "IL26", "IL10", "LINC02195", "HLA-DRA", "CTLA4")
  my_genes <- (
    de %>% group_by(celltype) %>%
      top_n(wt = (abs(logFC) * -log10(P.Value)), n = 5)
  )$Gene
  # my_genes <- (
  #   de_select %>% group_by(Gene) %>%
  #     mutate(wt = logFC / mean(logFC)) %>%
  #     ungroup %>% group_by(celltype) %>%
  #     top_n(wt = wt, n = 7)
  # )$Gene
  row_ha <- rowAnnotation(
    gene = anno_mark(
      at = match(my_genes, rownames(mat)),
      labels = my_genes,
      labels_gp = gpar(fontface = "italic")
    )
  )
  celltype_colors <- mpn65[1:length(unique(de$celltype))]
  names(celltype_colors) <- as.character(seq_along(celltype_colors))
  column_ha <- HeatmapAnnotation(
    celltype = as.character(colnames(mat)),
    col = list(celltype = celltype_colors)
  )
  col_fun <- circlize::colorRamp2(
    seq(-max(abs(range(mat))), max(abs(range(mat))), length.out = 11),
    rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
  )
  ht <- Heatmap(
    matrix = mat,
    name = "z",
    # top_annotation = column_ha,
    # right_annotation = row_ha,
    # col = rev(scico(pal = "", n = 20)),
    # col = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
    col = col_fun,
    # row_split = 6,
    row_order = order(uwot::umap(X = mat, n_components = 1)[,1]),
    row_names_gp = gpar(fontface = "italic"),
    # show_row_names = FALSE,
    # show_row_dend = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 0,
    # column_names_gp = gpar(angle = 0, fontsize = 12, just = "center"),
    show_column_dend = FALSE
  )
  pdf_file <- file.path(out_dir, "heatmap-de-drugtargets-z.pdf")
  message(pdf_file)
  unlink(pdf_file)
  pdf(pdf_file, width = 5, height = nrow(mat) * 0.15 + 1)
  draw(ht)
  dev.off()
}


# UMAP with major colors
########################################################################

analyses <- c(
  "a12_4_4_t4_cd8_1_2",
  "a12_4_4_t4_cd4_2_2",
  "a12_4_4_m3_2",
  "a12_4_4_b5_1_3"
)
cell_ids <- list()
for (analysis_name in analyses) {
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  cell_ids[[analysis_name]] <- a1$obs$cell
}

# main_analysis <- "a12_4_4_4"
# m1_file <- as.character(glue("results/a20/{main_analysis}/data/{main_analysis}.qs"))
# print_status(glue("Reading {m1_file}"))
# m1 <- qread(m1_file)
# print_status(glue("done"))
# out_dir <- glue("results/a20/{main_analysis}/figures")

main_analysis <- "a12_4_4_min_genes500_n_pcs20"
m1_file <- as.character(glue("results/{main_analysis}/data/{main_analysis}.rds"))
print_status(glue("Reading {m1_file}"))
m1 <- readRDS(m1_file)
print_status(glue("done"))
out_dir <- glue("results/{main_analysis}/figures")

#tcell_clusters <- c(1, 2, 4, 5, 6, 7, 9, 16, 17, 19, 20, 21)
#bcell_clusters <- c(8, 10, 15, 23, 24, 11, 3, 18, 12)
#plasma_clusters <- c(13)
##
#int3$obs$cluster_major <- "Other"
#int3$obs$cluster_major[int3$obs$leiden %in% tcell_clusters] <- "T cells"
#int3$obs$cluster_major[int3$obs$leiden %in% bcell_clusters] <- "B cells"
#int3$obs$cluster_major[int3$obs$leiden %in% plasma_clusters] <- "Plasma cells"
#int3$obs$cluster_major[int3$obs$leiden %in% c(22)] <- "Mast cells"
#int3$obs$cluster_major[int3$obs$leiden %in% c(14)] <- "Myeloid cells"

# m1$obs$cluster_major <- "Other"
m1$obs$cluster_major <- NA
m1$obs$cluster_major[m1$obs$cell %in% cell_ids[["a12_4_4_t4_cd8_1_2"]]] <- "CD8 T cells"
m1$obs$cluster_major[m1$obs$cell %in% cell_ids[["a12_4_4_t4_cd4_2_2"]]] <- "CD4 T cells"
m1$obs$cluster_major[m1$obs$cell %in% cell_ids[["a12_4_4_m3_2"]]] <- "Myeloid cells"
m1$obs$cluster_major[m1$obs$cell %in% cell_ids[["a12_4_4_b5_1_3"]]] <- "B cells"
m1$obs$cluster_major[m1$obs$leiden %in% c(1, 7, 12, 14)] <- "Plasma cells"
m1$obs$cluster_major[m1$obs$leiden %in% c(17)] <- "Mast cells"
table(m1$obs$cluster_major, useNA = "always")
#
set.seed(1)
ix_rand <- sample(nrow(m1$obs))
p <- plot_scattermore(
  x = m1$obs$UMAP1[ix_rand],
  y = m1$obs$UMAP2[ix_rand],
  group = m1$obs$cluster_major[ix_rand],
  group_colors = mpn65,
  group_labels = FALSE,
  group_legend = TRUE,
  alpha = 0.25,
  pixels = 500
)
my_ggsave(
  "umap-cluster_major",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.75,
  width = 10, height = 5,
  units = "in", dpi = 300
)

set.seed(1)
ix_rand <- sample(nrow(m1$obs))
p <- plot_scattermore(
  x = m1$obs$UMAP1[ix_rand],
  y = m1$obs$UMAP2[ix_rand],
  group = pub_ids[as.character(m1$obs$donor[ix_rand])],
  group_colors = mpn65,
  group_labels = FALSE,
  group_legend = TRUE,
  alpha = 0.25,
  pixels = 500
)
my_ggsave(
  "umap-donor-scattermore",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.75,
  width = 11, height = 5,
  units = "in", dpi = 300
)

x <- m1$obs %>%
  group_by(donor) %>%
  count(cluster_major)
x$cluster_major[is.na(x$cluster_major)] <- "Other"
x <- x %>% group_by(donor) %>% mutate(pct = 100 * n / sum(n))
# x$donor <- pub_ids[x$donor]
x$pub <- pub_ids[x$donor]
x$pub <- naturalfactor(x$pub)
x$pub <- factor(x$pub, rev(levels(x$pub)))
p1 <- ggplot(x) +
  aes(x = pct, y = pub, fill = cluster_major) +
  geom_colh() +
  scale_fill_manual(values = mpn65[c(1,2,3,4,5,7,6)]) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none"
  ) +
  labs(x = "Percent", y = NULL)
x_pub <- x %>% group_by(pub) %>% summarize(n = sum(n))
p2 <- ggplot(x) +
  geom_colh(
    data = x,
    mapping = aes(x = n, y = pub, fill = cluster_major)
  ) +
  geom_text(
    data = x_pub,
    mapping = aes(x = n, y = pub, label = comma(n, accuracy = 1)),
    hjust = 0, nudge_x = 100
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.25)), labels = label_number_si()
  ) +
  scale_fill_manual(name = NULL, values = mpn65[c(1,2,3,4,5,7,6)]) +
  labs(x = "Cells", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
my_ggsave(
  "bars-cluster_major",
  out_dir = out_dir,
  type = "pdf",
  plot = p1 + p2 + plot_annotation(title = "Composition of each donor's cells"),
  scale = 0.75,
  width = 11, height = 10,
  units = "in", dpi = 300
)

m1$log2cpm <- do_log2cpm(m1$counts, median(colSums(m1$counts)))
this_symbol <- "MS4A1"
my_symbols <- c("HPGDS", "XBP1", "MS4A1", "JCHAIN", "LYZ", "CD3D", "CD8A", "CD4")
for (this_symbol in my_symbols) {
  this_gene <- names(which(ensembl_to_symbol == this_symbol))
  p <- plot_hexgene(
    x            = m1$obs$UMAP1,
    y            = m1$obs$UMAP2,
    z            = as.numeric(m1$log2cpm[this_gene,]),
    bins         = 101,
    palette      = "davos",
    direction    = -1,
    use_quantile = FALSE,
    text = FALSE
  ) +
  labs(title = this_symbol) +
  theme(legend.position = "none")
  my_ggsave(
    slug = glue("umap-{this_symbol}"),
    out_dir = glue("{out_dir}/gene"),
    plot = p,
    type = "pdf",
    scale = 0.8, width = 3.5, height = 3, units = "in", dpi = 300
  )
}


# Blood main clusters
########################################################################

m1_file <- as.character(glue("results/blood/blood2/blood2.qs"))
print_status(glue("Reading {m1_file}"))
m1 <- qread(m1_file)
print_status(glue("done"))
out_dir <- glue("results/blood/blood2/figures")

m1$obs$cluster <- m1$obs$leiden0.933
blood_clusters <- list(
  "tcell" = c(2, 3, 4, 5, 6, 11, 12, 13, 14),
  "mait" = c(10),
  "nk" = c(9),
  "bcell" = c(8),
  "myeloid" = c(1, 7, 13, 15, 16, 17)
)
m1$obs$cluster_major <- "Other"
m1$obs$cluster_major[m1$obs$cluster %in% blood_clusters$tcell] <- "T cells"
m1$obs$cluster_major[m1$obs$cluster %in% blood_clusters$mait] <- "MAIT cells"
m1$obs$cluster_major[m1$obs$cluster %in% blood_clusters$nk] <- "NK cells"
m1$obs$cluster_major[m1$obs$cluster %in% blood_clusters$bcell] <- "B cells"
m1$obs$cluster_major[m1$obs$cluster %in% blood_clusters$myeloid] <- "Myeloid cells"
table(m1$obs$cluster_major, useNA = "always")


#
set.seed(1)
ix_rand <- sample(nrow(m1$obs))
p <- plot_scattermore(
  x = m1$obs$UMAP1[ix_rand],
  y = m1$obs$UMAP2[ix_rand],
  group = m1$obs$cluster_major[ix_rand],
  group_colors = mpn65,
  group_labels = FALSE,
  group_legend = TRUE,
  alpha = 0.25,
  pixels = 500
)
my_ggsave(
  "umap-cluster_major",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.75,
  width = 10, height = 5,
  units = "in", dpi = 300
)



p <- plot_umap_by_factor(m1$obs, "cluster")
my_ggsave(
  "umap-facet-by-cluster",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1,
  width = 10, height = 6,
  units = "in", dpi = 300
)

set.seed(1)
ix_rand <- sample(nrow(m1$obs))
p <- plot_scattermore(
  x = m1$obs$UMAP1[ix_rand],
  y = m1$obs$UMAP2[ix_rand],
  group = pub_ids[as.character(m1$obs$donor[ix_rand])],
  group_colors = mpn65,
  group_labels = FALSE,
  group_legend = TRUE,
  alpha = 0.25,
  pixels = 500
)
my_ggsave(
  "umap-donor-scattermore",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.75,
  width = 11, height = 5,
  units = "in", dpi = 300
)

x <- m1$obs %>%
  dplyr::group_by(donor) %>%
  dplyr::count(cluster_major)
x$cluster_major[is.na(x$cluster_major)] <- "Other"
x <- x %>% group_by(donor) %>% mutate(pct = 100 * n / sum(n))
# x$donor <- pub_ids[x$donor]
x$pub <- pub_ids[x$donor]
x$pub <- naturalfactor(x$pub)
x$pub <- factor(x$pub, rev(levels(x$pub)))
p1 <- ggplot(x) +
  aes(x = pct, y = pub, fill = cluster_major) +
  geom_colh() +
  scale_fill_manual(values = mpn65[c(1,2,3,4,5,7,6)]) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none"
  ) +
  labs(x = "Percent", y = NULL)
x_pub <- x %>% group_by(pub) %>% summarize(n = sum(n))
p2 <- ggplot(x) +
  geom_colh(
    data = x,
    mapping = aes(x = n, y = pub, fill = cluster_major)
  ) +
  geom_text(
    data = x_pub,
    mapping = aes(x = n, y = pub, label = comma(n, accuracy = 1)),
    hjust = 0, nudge_x = 100
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.3)), labels = label_number_si()
  ) +
  scale_fill_manual(name = NULL, values = mpn65[c(1,2,3,4,5,7,6)]) +
  labs(x = "Cells", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
my_ggsave(
  "bars-cluster_major",
  out_dir = out_dir,
  type = "pdf",
  plot = p1 + p2 + plot_annotation(title = "Composition of each donor's cells"),
  scale = 0.75,
  width = 9.5, height = length(unique(m1$obs$donor)) * 0.5,
  units = "in", dpi = 300
)


# cluster_major2 {{{

m1$obs$cluster_major2 <- "Other"
cell_ids <- list()
for (my_type in c("b", "cd4", "cd8", "myeloid")) {
  cell_ids[[my_type]] <- h5read(glue("paper/blood-{my_type}.h5ad"), "/obs/cell")
}
ix <- m1$obs$cell %in% cell_ids$b
m1$obs$cluster_major2[m1$obs$cell %in% cell_ids$b] <- "B cell"
m1$obs$cluster_major2[m1$obs$cell %in% cell_ids$cd4] <- "CD4 T"
m1$obs$cluster_major2[m1$obs$cell %in% cell_ids$cd8] <- "CD8 T/NK/GDT"
m1$obs$cluster_major2[m1$obs$cell %in% cell_ids$myeloid] <- "MNPs"
m1$obs %>% count(cluster_major2)
length(cell_ids$b)
length(cell_ids$cd8)
m1$obs %>% count(cluster_major, cluster_major2)
#
set.seed(1)
ix <- which(m1$obs$cluster_major2 != "Other")
ix <- sample(ix)
p <- plot_scattermore(
  x = m1$obs$UMAP1[ix],
  y = m1$obs$UMAP2[ix],
  group = m1$obs$cluster_major2[ix],
  group_colors = mpn65[c(1,2,3,5)],
  group_labels = FALSE,
  group_legend = TRUE,
  alpha = 0.25,
  pixels = 500
)
my_ggsave(
  "umap-cluster_major2",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.75,
  width = 10, height = 5,
  units = "in", dpi = 300
)
p <- plot_scattermore(
  x = m1$obs$UMAP1[ix],
  y = m1$obs$UMAP2[ix],
  group = pub_ids[as.character(m1$obs$donor[ix])],
  group_colors = mpn65,
  group_labels = FALSE,
  group_legend = TRUE,
  alpha = 0.25,
  pixels = 500
)
my_ggsave(
  "umap-donor-scattermore2",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 0.75,
  width = 11, height = 5,
  units = "in", dpi = 300
)

x <- m1$obs %>%
  filter(cluster_major2 != "Other") %>%
  dplyr::group_by(donor) %>%
  dplyr::count(cluster_major2)
x <- x %>% group_by(donor) %>% mutate(pct = 100 * n / sum(n))
# x$donor <- pub_ids[x$donor]
x$pub <- pub_ids[x$donor]
x$pub <- naturalfactor(x$pub)
x$pub <- factor(x$pub, rev(levels(x$pub)))
p1 <- ggplot(x) +
  aes(x = pct, y = pub, fill = cluster_major2) +
  geom_colh() +
  scale_fill_manual(values = mpn65[c(1,2,3,5,7,6)]) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none"
  ) +
  labs(x = "Percent", y = NULL)
x_pub <- x %>% group_by(pub) %>% summarize(n = sum(n))
p2 <- ggplot(x) +
  geom_colh(
    data = x,
    mapping = aes(x = n, y = pub, fill = cluster_major2)
  ) +
  geom_text(
    data = x_pub,
    mapping = aes(x = n, y = pub, label = comma(n, accuracy = 1)),
    hjust = 0, nudge_x = 100
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.3)), labels = label_number_si()
  ) +
  scale_fill_manual(name = NULL, values = mpn65[c(1,2,3,5,7,6)]) +
  labs(x = "Cells", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
my_ggsave(
  "bars-cluster_major2",
  out_dir = out_dir,
  type = "pdf",
  plot = p1 + p2 + plot_annotation(title = "Composition of each donor's cells"),
  scale = 0.75,
  width = 9.5, height = length(unique(m1$obs$donor)) * 0.5,
  units = "in", dpi = 300
)

# }}}


m1$log2cpm <- do_log2cpm(m1$counts, median(colSums(m1$counts)))

my_symbols <- c("HPGDS", "XBP1", "MS4A1", "JCHAIN", "LYZ", "CD3D", "CD8A", "CD4", "SH2D1B", "TRAV1-2")
for (this_symbol in my_symbols) {
  this_gene <- names(which(ensembl_to_symbol == this_symbol))
  p <- plot_hexgene(
    x            = m1$obs$UMAP1,
    y            = m1$obs$UMAP2,
    z            = as.numeric(m1$log2cpm[this_gene,]),
    bins         = 101,
    palette      = "davos",
    direction    = -1,
    use_quantile = FALSE,
    text = FALSE
  ) +
  labs(title = this_symbol) +
  theme(legend.position = "none")
  my_ggsave(
    slug = glue("umap-{this_symbol}"),
    out_dir = glue("{out_dir}/gene"),
    plot = p,
    type = "pdf",
    scale = 0.8, width = 3.5, height = 3, units = "in", dpi = 300
  )
}



##
#obs_file <- "results/a20/a12_4_4_4/data/umap-coords.tsv"
#obs <- fread(obs_file)
#out_dir <- file.path(dirname(dirname(obs_file)), "figures")
#cluster_colors <- mpn65
#names(cluster_colors) <- seq_along(cluster_colors)
#n_donors <- length(unique(obs$donor))
#set.seed(1)
#ix_rand <- sample(nrow(m1$obs))
#stopifnot(all(obs$cell[ix_rand] == m1$obs$cell[ix_rand]))
#p1 <- plot_scattermore(
#  x = obs$UMAP1_4[ix_rand],
#  y = obs$UMAP2_4[ix_rand],
#  # group = obs$leiden1.8,
#  # group_colors = cluster_colors,
#  group = m1$obs$cluster_major[ix_rand],
#  group_colors = mpn65[2:8],
#  # group_colors = pals::okabe(8)[2:8],
#  group_labels = FALSE,
#  group_legend = TRUE,
#  pixels = 1000,
#  alpha = 0.35
#) +
#labs(
#  title = glue(
#    "{length(unique(obs$leiden1.8))} clusters of {comma(nrow(obs))} cells from {n_donors} donors"
#  )
#)
#my_ggsave(
#  glue("umap-major-types"),
#  out_dir = out_dir,
#  plot = p1,
#  type = "pdf",
#  scale = 1, width = 9, height = 5, units = "in", dpi = 300
#)



# CD8
# Full figures with donor-level volcano, cluster-level boxplot, lanes
########################################################################

# analyses <- c(
#   "a12_4_4_t4_cd8_1_2",
#   "a12_4_4_t4_cd4_2_2",
#   "a12_4_4_m3_2",
#   "a12_4_4_b5_1_3"
# )

#{
#
#  analysis_name <- "a12_4_4_t4_cd8_1_2"
#  params <- list(
#    min_cells_in_cluster = 50,
#    min_percent_of_cells_with_gene = 5
#  )
#  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
#  print_status(glue("Reading {a1_file}"))
#  a1 <- qread(a1_file)
#  print_status(glue("done"))
#  # my_leiden <- "leiden0.933"
#  #
#  # my_leiden <- "leiden1.22"
#  # a1$obs$cluster <- a1$obs[[my_leiden]]
#  # a1$obs$cluster <- recluster_cd8_leiden122(a1$obs$cluster)
#  #
#  my_leiden <- "leiden1.51"
#  a1$obs$cluster <- a1$obs[[my_leiden]]
#  a1$obs$cluster <- recluster_cd8_leiden151(a1$obs$cluster)
#	#
#	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
#	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#	#
#	sample_info <- janitor::clean_names(read_excel(
#			path = "data/luoma_villani_combined_colitis_samples2.xlsx",
#			sheet = 3
#	))
#	donor_to_sex <- unlist(split(sample_info$sex, sample_info$donor))
#  #
#	a1$obs$sex <- donor_to_sex[a1$obs$donor]
#	a1$obs$drug <- factor(a1$obs$drug, c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4"))
#  ## Skip cell clusters that have too few cells.
#  #exclude_clusters <- (
#  #  a1$obs %>% count(cluster) %>% filter(n < params[["min_cells_in_cluster"]])
#  #)$leiden
#  ##
#  #a1$obs <- a1$obs[!a1$obs$cluster %in% exclude_clusters,]
#  #
#  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))
#
#
#  # Pseudobulk at the donor level
#  ########################################################################
#  y <- with(a1$obs, model.matrix(~ 0 + factor(donor)))
#  y <- as(y, "dgCMatrix")
#  # y <- sweep(y, 2, colSums(y), "/") # means
#  pb <- as(a1$counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#  #
#  pb_meta <- data.frame(donor = colnames(pb))
#  pb_meta <- as_tibble(pb_meta)
#  pb_meta %<>%
#    mutate(
#      donor = str_replace(donor, "factor\\(donor\\)", "")
#    )
#  colnames(pb) <- str_replace(colnames(pb), "factor\\(donor\\)", "")
#  pb_meta <- left_join(
#    pb_meta,
#    a1$obs %>%
#      select(
#        donor, case
#      ) %>%
#      group_by(donor, case) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#
#  # x <- table(a1$obs$cluster)
#  # min_percent <- 100 * min(x) * 0.25 / sum(x)
#  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
#  keep_ens <- rownames(pb)[rowMeans(pb) > 0.5]
#  #
#  des1 <- with(pb_meta, model.matrix(
#    ~ case
#  ))
#  ob <- as.matrix(pb[keep_ens,])
#  fit1 <- lmFit(object = ob, design = des1)
#  fit1 <- eBayes(fit1)
#  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#  de_donor <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE) %>%
#    dplyr::rename(Gene = ID)
#  de_donor$ensembl_id <- rownames(de_donor)
#  de_donor <- as_tibble(de_donor)
#  de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
#  data.table::fwrite(de_donor, de_donor_file, sep = "\t")
#  # de_donor$GeneName <- de_donor$Gene
#  write_de_xlsx(de_donor, file.path(out_dir, "de_donor_case-vs-control.xlsx"))
#
#
#  my_ids <- (
#    de_donor %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    arrange(-logFC)
#  )$ensembl_id[1:20]
#  # my_ids <- de_donor$ensembl_id[1:15]
#  d <- as_tibble(reshape2::melt(as.matrix(pb)[my_ids,]))
#  colnames(d) <- c("ensembl_id", "donor", "value")
#  d$GeneName <- factor(
#    as.character(unname(ensembl_to_symbol[as.character(d$ensembl_id)])),
#    rev(as.character(unname(ensembl_to_symbol[as.character(my_ids)])))
#  )
#  d$Gene <- as.integer(d$GeneName)
#  d <- left_join(d, pb_meta, by = "donor")
#  #
#  p_donor <- ggplot(d) +
#      annotate(
#        geom = "rect",
#        xmin = -Inf,
#        xmax = Inf,
#        ymin = seq(from = 1, to = length(my_ids), by = 2) - 0.5,
#        ymax = seq(from = 1, to = length(my_ids), by = 2) + 0.5,
#        alpha = 0.2
#      ) +
#      # geom_vline(xintercept = 0, size = 0.3) +
#      geom_point(
#        data = d,
#        mapping = aes(x = value, y = Gene, group = case, fill = case),
#        position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE),
#        shape = 21, size = 2, stroke = 0.2
#      ) +
#      scale_fill_manual(values = okabe(8)) +
#      # scale_x_continuous(
#      #   breaks = seq(-10, 10, by = 1),
#      #   labels = function(x) fractional::fractional(2 ^ x)
#      # ) +
#      scale_y_continuous(
#        expand = c(0, 0),
#        breaks = seq(1, max(d$Gene)),
#        labels = levels(d$GeneName)
#      ) +
#      guides(fill = guide_legend(reverse = TRUE, title = NULL, override.aes = list(size = 5))) +
#      # annotation_logticks(sides = "b", size = 0.3) +
#      expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#      labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
#      theme(
#        axis.text.y = element_text(face = "italic"),
#        legend.margin = margin(0, 0, 0, 0),
#        legend.box.margin = margin(0, 0, 0, 0)
#      )
#  my_ggsave(
#    "de_donor",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p_donor,
#    scale = 1,
#    width = 5,
#    height = length(my_ids) * 0.3 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
#
#  for (my_id in my_ids) {
#    p <- plot_hexgene(
#      x = a1$obs$UMAP1,
#      y = a1$obs$UMAP2,
#      z = as.numeric(a1$log2cpm[my_id,]),
#      group = append_n(a1$obs$case),
#      bins = 37,
#      legend = FALSE,
#      palette = "oslo"
#      # palette = "lajolla",
#      # direction = 1
#    ) +
#    facet_wrap(~ group) +
#    labs(title = ensembl_to_symbol[my_id]) +
#    theme(
#      panel.spacing = unit(0.5, "lines"),
#      plot.title = element_text(face = "italic")
#    )
#    my_ggsave(
#      glue("umap-{safe(ensembl_to_symbol[my_id])}"),
#      out_dir = out_dir,
#      type = "pdf",
#      plot = p,
#      # scale = 0.8,
#      scale = 0.6,
#      width = 7,
#      height = 4,
#      units = "in",
#      dpi = 300
#    )
#  }
#
#  p <- plot_relative_hex(
#    x = a1$obs$UMAP1,
#    y = a1$obs$UMAP2,
#    group = a1$obs$case == "Case",
#    bins = 91,
#    # palette = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
#    # palette = rev(RColorBrewer::brewer.pal(name = "PiYG", n = 11))
#    # palette = rev(RColorBrewer::brewer.pal(name = "Spectral", n = 11))
#    # palette = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 11))
#    palette = "roma",
#    # direction = -1
#  ) +
#  labs(title = "Enrichment of cells from Cases")
#  my_ggsave(
#    glue("umap-case-enrichment"),
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 4.2,
#    height = 3,
#    units = "in",
#    dpi = 300
#  )
#  
#  # Pseudobulk at the cluster level
#  ########################################################################
#  y <- with(a1$obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
#  y <- as(y, "dgCMatrix")
#  # y <- sweep(y, 2, colSums(y), "/") # means
#  pb <- as(a1$counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#  #
#  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
#  colnames(pb_meta) <- c("cluster", "donor")
#  pb_meta <- as_tibble(pb_meta)
#  pb_meta %<>%
#    mutate(
#      cluster = str_replace(cluster, "factor\\(cluster\\)", ""),
#      donor = str_replace(donor, "factor\\(donor\\)", "")
#    )
#  pb_meta <- left_join(
#    pb_meta,
#    a1$obs %>%
#      select(
#        donor, case
#      ) %>%
#      group_by(donor, case) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#
#  # x <- table(a1$obs$cluster)
#  # min_percent <- 100 * min(x) * 0.25 / sum(x)
#  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
#  keep_ens <- rownames(pb)[rowMeans(pb) > 0.5]
#  #
#  de <- rbindlist(
#    pblapply(sort(unique(pb_meta$cluster)), function(this_cluster) {
#      ix_cluster <- pb_meta$cluster == this_cluster
#      des1 <- with(pb_meta[ix_cluster,], model.matrix(
#        ~ case
#      ))
#      ob <- as.matrix(pb[keep_ens,ix_cluster])
#      fit1 <- lmFit(object = ob, design = des1)
#      fit1 <- eBayes(fit1)
#      fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#      res <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE)
#      res$cluster <- this_cluster
#      res <- res %>% rename(Gene = ID)
#      res$ensembl_id <- rownames(res)
#      return(res)
#    })
#  )
#  de$GeneName <- de$Gene
#  de$cluster <- naturalfactor(de$cluster)
#  de_file <- glue("{out_dir}/de_case-vs-control.tsv.gz")
#  data.table::fwrite(de, de_file, sep = "\t")
#  #
#  write_de_xlsx(
#    de %>% mutate_if(is.numeric, signif, 4) %>% select(-GeneName),
#    glue("{out_dir}/de_case-vs-control.xlsx"),
#    col = "cluster"
#  )
#
#  de_summary <- de %>% group_by(cluster) %>%
#    summarize(
#      n_up = sum(logFC > log2(1.5) & adj.P.Val < 0.05),
#      n_down = sum(-logFC > log2(1.5) & adj.P.Val < 0.05)
#    )
#  de_summary_file <- glue("{out_dir}/de_summary_case-vs-control.tsv")
#  data.table::fwrite(de_summary, de_summary_file, sep = "\t")
#
#  p <- ggplot(
#    de_summary %>% pivot_longer(cols = c("n_up", "n_down")) %>%
#      mutate(value = ifelse(name == "n_down", -value, value))
#  ) +
#    aes(x = value, y = cluster, fill = name) +
#    geom_colh() +
#    scale_y_discrete(limits = rev(levels(naturalfactor(de_summary$cluster)))) +
#    scale_fill_manual(
#      # values = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[c(9,3)],
#      values = okabe(2),
#      guide = "none"
#    ) +
#    scale_x_continuous(labels = abs) +
#    labs(x = NULL, y = NULL)
#  my_ggsave(
#    "bars-summary",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 2,
#    height = 3,
#    units = "in",
#    dpi = 300
#  )
#
#  # logFC within each cluster for the genes discovered at donor-level
#  # mat <- dcast(de %>% filter(ensembl_id %in% my_ids), ensembl_id ~ cluster, value.var = "logFC")
#  # mat_rows <- mat[[1]]
#  # mat <- as.matrix(mat[,2:ncol(mat)])
#  # rownames(mat) <- mat_rows
#
#  cluster_colors <- mpn65
#  names(cluster_colors) <- as.character(seq_along(cluster_colors))
#  x <- de %>% filter(ensembl_id %in% my_ids)
#  x$ensembl_id = factor(x$ensembl_id, rev(my_ids))
#  p_mat <- ggplot(x) +
#  geom_tile(
#    aes(x = cluster, y = ensembl_id, fill = logFC)
#  ) +
#  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
#  scale_y_discrete(
#    position = "l", name = NULL, expand = c(0, 0),
#    labels = function(x) ensembl_to_symbol[x]
#  ) +
#  scale_fill_gradientn(
#    colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))[6:11],
#    # colors = RColorBrewer::brewer.pal(name = "Reds", n = 9),
#    # limits = c(0, 1),
#    name = "log2FC", guide = guide_colorbar(barheight = 10), breaks = pretty_breaks(5)
#  ) +
#  theme(
#    axis.text.y = element_text(size = 14, face = "italic"),
#    legend.position = "right"
#  )
#  p_bar <- ggplot(de %>% count(cluster)) +
#  geom_tile(
#    aes(y = 1, x = cluster, fill = cluster)
#  ) +
#  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
#  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
#  scale_fill_manual(values = cluster_colors, guide = "none")
#  # p <- (
#  #       p_bar + theme(plot.margin = margin(r = 0))
#  #     ) / (
#  #       p_mat + theme(
#  #         axis.text.x = element_blank(),
#  #         axis.ticks.x = element_blank(),
#  #         axis.title.y = element_blank(),
#  #         plot.margin = margin(l = 0)
#  #       )
#  #     ) + plot_layout(heights = c(1, length(my_ids)))
#  p <- (
#      plot_spacer() + p_bar + theme(plot.margin = margin(l = 2))
#  ) / (
#    (
#      p_donor + theme(legend.position = "none")
#    ) + (
#      p_mat + theme(
#        axis.text.y = element_blank(),
#        axis.text.x = element_blank(),
#        axis.ticks.x = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.title.y = element_blank(),
#        plot.margin = margin(l = 0)
#      )
#    )
#  ) + plot_layout(heights = c(1, length(my_ids)))
#  fig_height <- length(my_ids) * 0.3 + 2
#  fig_width <- length(unique(de$cluster)) * 0.4 + 4
#  my_ggsave(
#    slug = "de_donor_heatmap",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
#  )
#
#
#  for (this_cluster in unique(pb_meta$cluster)) {
#    ix <- pb_meta$cluster == this_cluster
#    n_donors <- table(pb_meta$case[ix])
#    #
#    p <- plot_limma_volcano(de[de$cluster == this_cluster,]) +
#      labs(title = glue("Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
#    my_ggsave(
#      glue("volcano-cluster-{this_cluster}"),
#      out_dir = glue("{out_dir}"),
#      plot = p,
#      type = "pdf",
#      scale = 1, width = 5, height = 4, units = "in", dpi = 300
#    )
#    #
#    my_ids <- c(
#      (
#        de %>%
#        filter(cluster == this_cluster) %>%
#        filter(logFC > 0) %>%
#        top_n(n = 20, wt = -log10(P.Value) * logFC)
#      )$ensembl_id,
#      (
#        de %>%
#        filter(cluster == this_cluster) %>%
#        filter(logFC < 0) %>%
#        top_n(n = 20, wt = -log10(P.Value) * -logFC)
#      )$ensembl_id
#    )
#    # my_ids <- head((
#    #   de %>%
#    #   filter(cluster == this_cluster) %>%
#    #   filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
#    #   arrange(abs(logFC) + -log10(P.Value))
#    # )$ensembl_id, 50)
#    mat <- as.matrix(pb[my_ids,ix])
#    a_col <- as.data.frame(pb_meta)[ix,]
#    a_col <- left_join(a_col, unique(a1$obs[,c("donor","drug")]), by = "donor")
#    rownames(a_col) <- colnames(mat)
#    a_col <- a_col[,-1]
#    h1 <- dendsort(hclust(dist(mat), method = "complete"))
#    h2 <- dendsort(hclust(dist(t(mat)), method = "complete"))
#    a_col <- a_col[h2$order,]
#    mat <- mat[h1$order,h2$order]
#    #
#    a_colors <- list()
#    a_colors[["drug"]] <- okabe(8)[4:7]
#    names(a_colors[["drug"]]) <- c("None", "PD-1", "CTLA-4", "PD-1/CTLA-4")
#    a_colors[["case"]] <- okabe(2)
#    names(a_colors[["case"]]) <- c("Control", "Case")
#    #
#    heatmap_file <- glue("{out_dir}/heatmap-cluster-{this_cluster}.pdf")
#    message(heatmap_file)
#    pheatmap(
#      filename = heatmap_file,
#      width = ncol(mat) * 0.2 + 3,
#      height = nrow(mat) * 0.1 + 3,
#      mat = mat,
#      cluster_col = FALSE,
#      cluster_row = FALSE,
#      hclust_method = "average",
#      show_colnames = FALSE,
#      show_rownames = TRUE,
#      scale = "row",
#      # color = rev(scico::scico(palette = "roma", n = 20)),
#      color = rev(RColorBrewer::brewer.pal("RdBu", n = 11)),
#      border_color = NA,
#      labels_row = ensembl_to_symbol[rownames(mat)],
#      # labels_col = labels_col,
#      annotation_col = a_col[,c("drug","case")],
#      annotation_colors = a_colors,
#      main = glue("Cluster {this_cluster}")
#    )
#  }
#
#  plot_de_by_cluster <- function(d) {
#    ggplot() +
#      annotate(
#        geom = "rect",
#        xmin = -Inf,
#        xmax = Inf,
#        ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#        ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#        alpha = 0.2
#      ) +
#      geom_vline(xintercept = 0, size = 0.3) +
#      geom_point(
#        data = d,
#        mapping = aes(x = logFC, y = Gene)
#      ) +
#      geom_errorbarh(
#        data = d,
#        mapping = aes(xmin = CI.L, xmax = CI.R, y = Gene),
#        height = 0
#      ) +
#      scale_x_continuous(
#        breaks = seq(-10, 10, by = 1),
#        labels = function(x) fractional::fractional(2 ^ x)
#      ) +
#      scale_y_continuous(
#        expand = c(0, 0),
#        breaks = seq(1, max(d$Gene)),
#        labels = levels(d$GeneName)
#      ) +
#      annotation_logticks(sides = "b", size = 0.3) +
#      expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#      facet_grid(~ cluster)
#    # facet_wrap(~ cluster, scales = "free_x")
#  }
#
#  #
#  my_genes <- unique((
#    de %>%
#    group_by(cluster) %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    top_n(wt = -log10(P.Value), n = 5) %>%
#    select(GeneName, Gene, logFC, P.Value, adj.P.Val, cluster) %>%
#    ungroup() %>%
#    group_by(GeneName) %>%
#    mutate(min_p = min(P.Value)) %>%
#    ungroup() %>%
#    arrange(min_p)
#  )$GeneName)
#  #
#  d <- de %>% filter(GeneName %in% my_genes)
#  d$GeneName <- factor(d$GeneName, rev(my_genes))
#  d$Gene <- as.integer(d$GeneName)
#  d$cluster <- naturalfactor(d$cluster)
#  #
#  p <- plot_de_by_cluster(d) +
#    theme(
#      axis.title.y = element_blank(),
#      axis.text.y = element_text(face = "italic"),
#      panel.grid.major.x = element_line()
#    )
#  my_ggsave(
#    "lanes_case-vs-control",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = length(unique(de$cluster)) * 2.8 + 2,
#    height = length(my_genes) * 0.2 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  #
#  d <- de %>% filter(ensembl_id %in% my_ids, cluster %in% c(4, 6, 8, 9))
#  d$GeneName <- factor(d$GeneName)
#  d$Gene <- as.integer(d$GeneName)
#  d$cluster <- naturalfactor(d$cluster)
#  #
#  p <- plot_de_by_cluster(d) +
#    theme(
#      axis.title.y = element_blank(),
#      axis.text.y = element_text(face = "italic"),
#      panel.grid.major.x = element_line()
#    )
#  my_ggsave(
#    "lanes_case-vs-control-cluster-4-6-8-9",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = length(unique(d$cluster)) * 2.8 + 2,
#    height = length(unique(d$GeneName)) * 0.25 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#
#  # d <- rbindlist(lapply(my_ids, function(my_id) {
#  #   d <- a1$obs %>% select(cell, cluster, donor, case) %>% as_tibble
#  #   d$gene <- a1$counts[my_id,]
#  #   d <- d %>%
#  #     filter(cluster %in% c(4, 6, 8, 9)) %>%
#  #     group_by(cluster, donor, case) %>%
#  #     summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#  #   d$ensembl_id <- my_id
#  #   d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#  #   d
#  # }))
#  # d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  # d$Gene <- as.integer(d$GeneName)
#  #
#  # p <- ggplot() +
#  #   annotate(
#  #     geom = "rect",
#  #     xmin = 0,
#  #     xmax = Inf,
#  #     ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#  #     ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#  #     alpha = 0.2
#  #   ) +
#  #   geom_point(
#  #     data = d,
#  #     mapping = aes(x = percent, y = Gene, fill = case, group = case),
#  #     # position = position_quasirandom(dodge.width = 0.5, width = 0.3, groupOnX = FALSE),
#  #     shape = 21, size = 3, stroke = 0.2
#  #   ) +
#  #   scale_fill_manual(values = okabe(8), guide = "none") +
#  #   scale_x_log10() +
#  #   scale_y_continuous(
#  #     expand = c(0, 0),
#  #     breaks = seq(1, max(d$Gene)),
#  #     labels = levels(d$GeneName)
#  #   ) +
#  #   annotation_logticks(sides = "b", size = 0.3) +
#  #   expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#  #   labs(x = "Percent", y = NULL) +
#  #   facet_wrap(~ cluster, ncol = length(unique(d$cluster)))
#  # my_ggsave(
#  #   "percent",
#  #   out_dir = out_dir,
#  #   type = "pdf",
#  #   plot = p,
#  #   scale = 1,
#  #   width = length(unique(d$cluster)) * 1.5 + 1,
#  #   height = length(my_ids) * 0.3 + 1,
#  #   units = "in",
#  #   dpi = 300
#  # )
#
#  my_ids <- (
#    de_donor %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    arrange(-logFC)
#  )$ensembl_id[1:20]
#  d <- rbindlist(lapply(my_ids, function(my_id) {
#    d <- a1$obs %>% select(cell, cluster, case) %>% as_tibble
#    d$gene <- a1$counts[my_id,]
#    d <- d %>%
#      filter(cluster %in% c(4, 6, 8, 9)) %>%
#      group_by(cluster, case) %>%
#      summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#    d$ensembl_id <- my_id
#    d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#    d
#  }))
#  d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  d$Gene <- as.integer(d$GeneName)
#  #
#  p <- ggplot() +
#    annotate(
#      geom = "rect",
#      xmin = -Inf,
#      xmax = Inf,
#      ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#      ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#      alpha = 0.2
#    ) +
#    geom_vline(xintercept = seq(0, max(d$percent), by = 10), size = 0.3, color = "white") +
#    geom_colh(
#      data = d,
#      mapping = aes(x = percent, y = Gene, fill = case, group = case),
#      position = position_dodgev(),
#      size = 0.3
#    ) +
#    # geom_point(
#    #   data = d,
#    #   mapping = aes(x = percent, y = Gene, fill = case, group = case),
#    #   shape = 21, size = 3, stroke = 0.3
#    # ) +
#    # scale_x_log10() +
#    # annotation_logticks(sides = "b", size = 0.3) +
#    scale_fill_manual(values = okabe(8), guide = "none") +
#    scale_y_continuous(
#      expand = c(0, 0),
#      breaks = seq(1, max(d$Gene)),
#      labels = levels(d$GeneName)
#    ) +
#    expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#    labs(x = "Percent of cells with expression", y = NULL) +
#    facet_wrap(~ cluster, ncol = length(unique(d$cluster))) +
#    theme(axis.text.y = element_text(face = "italic"))
#  my_ggsave(
#    "percent",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    # width = 5,
#    width = length(unique(d$cluster)) * 1.5 + 1,
#    height = length(my_ids) * 0.3 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  #d <- rbindlist(lapply(my_ids, function(my_id) {
#  #  d <- a1$obs %>% select(cell, cluster, donor, case) %>% as_tibble
#  #  d$gene <- a1$counts[my_id,]
#  #  d <- d %>%
#  #    group_by(cluster, donor, case) %>%
#  #    summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#  #  d$ensembl_id <- my_id
#  #  d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#  #  d
#  #}))
#  #d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  ##
#  #mat <- dcast(d, ensembl_id ~ cluster + case + donor, value.var = "percent")
#  #mat_rownames <- mat[[1]]
#  #mat <- as.matrix(mat[,2:ncol(mat)])
#  #rownames(mat) <- mat_rownames
#  #mat <- mat[my_ids,]
#  #a_col <- as.data.frame(str_split_fixed(colnames(mat), "_", 3))
#  #colnames(a_col) <- c("cluster", "case", "donor")
#  #rownames(a_col) <- colnames(mat)
#  ##
#  #a_colors <- list()
#  #a_colors[["cluster"]] <- mpn65[1:length(unique(a_col$cluster))]
#  #names(a_colors[["cluster"]]) <- as.character(naturalsort(unique(a_col$cluster)))
#  #a_colors[["case"]] <- okabe(2)
#  #names(a_colors[["case"]]) <- c("Control", "Case")
#  ##
#  #heatmap_file <- glue("{out_dir}/heatmap-cluster-percent.pdf")
#  #pheatmap(
#  #  filename = heatmap_file,
#  #  width = ncol(mat) * 0.02 + 2,
#  #  height = nrow(mat) * 0.1 + 2,
#  #  mat = mat,
#  #  cluster_col = FALSE,
#  #  cluster_row = FALSE,
#  #  # hclust_method = "average",
#  #  show_colnames = FALSE,
#  #  # color = rev(scico::scico(palette = "davos", n = 20)),
#  #  # scale = "row",
#  #  # color = rev(scico::scico(palette = "roma", n = 20)),
#  #  color = rev(scico::scico(palette = "davos", n = 20)),
#  #  border_color = NA,
#  #  labels_row = ensembl_to_symbol[rownames(mat)],
#  #  # labels_col = labels_col,
#  #  annotation_col = a_col[,c("cluster","case")],
#  #  annotation_colors = a_colors
#  #)
#
#}


# CD4
# Full figures with donor-level volcano, cluster-level boxplot, lanes
########################################################################

#{
#
#  analysis_name <- "a12_4_4_t4_cd4_2_2"
#  params <- list(
#    min_cells_in_cluster = 50,
#    min_percent_of_cells_with_gene = 5
#  )
#  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
#  print_status(glue("Reading {a1_file}"))
#  a1 <- qread(a1_file)
#  print_status(glue("done"))
#  if (analysis_name == "a12_4_4_t4_cd8_1_2") {
#    my_leiden <- "leiden0.933"
#  }
#  if (analysis_name == "a12_4_4_t4_cd4_2_2") {
#    my_leiden <- "leiden0.933"
#  }
#  if (analysis_name == "a12_4_4_m3_2") {
#    my_leiden <- "leiden0.933"
#  }
#  if (analysis_name == "a12_4_4_b5_1_3") {
#    my_leiden <- "leiden0.933"
#  }
#  a1$obs$cluster <- a1$obs[[my_leiden]]
#	#
#	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
#	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#	#
#	sample_info <- janitor::clean_names(read_excel(
#			path = "data/luoma_villani_combined_colitis_samples2.xlsx",
#			sheet = 3
#	))
#	donor_to_sex <- unlist(split(sample_info$sex, sample_info$donor))
#  #
#	a1$obs$sex <- donor_to_sex[a1$obs$donor]
#	a1$obs$drug <- factor(a1$obs$drug, c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4"))
#  ## Skip cell clusters that have too few cells.
#  #exclude_clusters <- (
#  #  a1$obs %>% count(cluster) %>% filter(n < params[["min_cells_in_cluster"]])
#  #)$leiden
#  ##
#  #a1$obs <- a1$obs[!a1$obs$cluster %in% exclude_clusters,]
#  #
#  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))
#
#
#  # Pseudobulk at the donor level
#  ########################################################################
#  y <- with(a1$obs, model.matrix(~ 0 + factor(donor)))
#  y <- as(y, "dgCMatrix")
#  # y <- sweep(y, 2, colSums(y), "/") # means
#  pb <- as(a1$counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#  #
#  pb_meta <- data.frame(donor = colnames(pb))
#  pb_meta <- as_tibble(pb_meta)
#  pb_meta %<>%
#    mutate(
#      donor = str_replace(donor, "factor\\(donor\\)", "")
#    )
#  colnames(pb) <- str_replace(colnames(pb), "factor\\(donor\\)", "")
#  pb_meta <- left_join(
#    pb_meta,
#    a1$obs %>%
#      select(
#        donor, case
#      ) %>%
#      group_by(donor, case) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#
#  # x <- table(a1$obs$cluster)
#  # min_percent <- 100 * min(x) * 0.25 / sum(x)
#  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
#  keep_ens <- rownames(pb)[rowMeans(pb) > 0.5]
#  #
#  des1 <- with(pb_meta, model.matrix(
#    ~ case
#  ))
#  ob <- as.matrix(pb[keep_ens,])
#  fit1 <- lmFit(object = ob, design = des1)
#  fit1 <- eBayes(fit1)
#  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#  de_donor <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE) %>%
#    rename(Gene = ID)
#  de_donor$ensembl_id <- rownames(de_donor)
#  de_donor <- as_tibble(de_donor)
#  # de_donor$GeneName <- de_donor$Gene
#  de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
#  data.table::fwrite(de_donor, de_donor_file, sep = "\t")
#  write_de_xlsx(de_donor, file.path(out_dir, "de_donor_case-vs-control.xlsx"))
#
#  my_ids <- (
#    de_donor %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    arrange(-logFC)
#  )$ensembl_id[1:20]
#  # my_ids <- de_donor$ensembl_id[1:15]
#  d <- as_tibble(reshape2::melt(as.matrix(pb)[my_ids,]))
#  colnames(d) <- c("ensembl_id", "donor", "value")
#  d$GeneName <- factor(
#    as.character(unname(ensembl_to_symbol[as.character(d$ensembl_id)])),
#    rev(as.character(unname(ensembl_to_symbol[as.character(my_ids)])))
#  )
#  d$Gene <- as.integer(d$GeneName)
#  d <- left_join(d, pb_meta, by = "donor")
#  #
#  p <- ggplot(d) +
#      annotate(
#        geom = "rect",
#        xmin = -Inf,
#        xmax = Inf,
#        ymin = seq(from = 1, to = length(my_ids), by = 2) - 0.5,
#        ymax = seq(from = 1, to = length(my_ids), by = 2) + 0.5,
#        alpha = 0.2
#      ) +
#      # geom_vline(xintercept = 0, size = 0.3) +
#      geom_point(
#        data = d,
#        mapping = aes(x = value, y = Gene, group = case, fill = case),
#        position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE),
#        shape = 21, size = 2, stroke = 0.2
#      ) +
#      scale_fill_manual(values = okabe(8)) +
#      # scale_x_continuous(
#      #   breaks = seq(-10, 10, by = 1),
#      #   labels = function(x) fractional::fractional(2 ^ x)
#      # ) +
#      scale_y_continuous(
#        expand = c(0, 0),
#        breaks = seq(1, max(d$Gene)),
#        labels = levels(d$GeneName)
#      ) +
#      guides(fill = guide_legend(reverse = TRUE, title = NULL, override.aes = list(size = 5))) +
#      # annotation_logticks(sides = "b", size = 0.3) +
#      expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#      labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
#      theme(
#        axis.text.y = element_text(face = "italic"),
#        legend.margin = margin(0, 0, 0, 0),
#        legend.box.margin = margin(0, 0, 0, 0)
#      )
#  my_ggsave(
#    "de_donor",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 5,
#    height = length(my_ids) * 0.3 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
#
#  for (my_id in my_ids) {
#    p <- plot_hexgene(
#      x = a1$obs$UMAP1,
#      y = a1$obs$UMAP2,
#      z = as.numeric(a1$log2cpm[my_id,]),
#      group = append_n(a1$obs$case)
#    ) +
#    facet_wrap(~ group) +
#    labs(title = ensembl_to_symbol[my_id]) +
#    theme(
#      panel.spacing = unit(0.5, "lines"),
#      plot.title = element_text(face = "italic")
#    )
#    my_ggsave(
#      glue("umap-{safe(ensembl_to_symbol[my_id])}"),
#      out_dir = out_dir,
#      type = "pdf",
#      plot = p,
#      scale = 0.8,
#      width = 7,
#      height = 5,
#      units = "in",
#      dpi = 300
#    )
#  }
#  
#  # Pseudobulk at the cluster level
#  ########################################################################
#  y <- with(a1$obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
#  y <- as(y, "dgCMatrix")
#  # y <- sweep(y, 2, colSums(y), "/") # means
#  pb <- as(a1$counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#  #
#  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
#  colnames(pb_meta) <- c("cluster", "donor")
#  pb_meta <- as_tibble(pb_meta)
#  pb_meta %<>%
#    mutate(
#      cluster = str_replace(cluster, "factor\\(cluster\\)", ""),
#      donor = str_replace(donor, "factor\\(donor\\)", "")
#    )
#  pb_meta <- left_join(
#    pb_meta,
#    a1$obs %>%
#      select(
#        donor, case
#      ) %>%
#      group_by(donor, case) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#
#  # x <- table(a1$obs$cluster)
#  # min_percent <- 100 * min(x) * 0.25 / sum(x)
#  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
#  keep_ens <- rownames(pb)[rowMeans(pb) > 0.5]
#  #
#  de <- rbindlist(
#    pblapply(sort(unique(pb_meta$cluster)), function(this_cluster) {
#      ix_cluster <- pb_meta$cluster == this_cluster
#      des1 <- with(pb_meta[ix_cluster,], model.matrix(
#        ~ case
#      ))
#      ob <- as.matrix(pb[keep_ens,ix_cluster])
#      fit1 <- lmFit(object = ob, design = des1)
#      fit1 <- eBayes(fit1)
#      fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#      res <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE)
#      res$cluster <- this_cluster
#      res <- res %>% rename(Gene = ID)
#      res$ensembl_id <- rownames(res)
#      return(res)
#    })
#  )
#  de$GeneName <- de$Gene
#  de$cluster <- naturalfactor(de$cluster)
#  de_file <- glue("{out_dir}/de_case-vs-control.tsv.gz")
#  data.table::fwrite(de, de_file, sep = "\t")
#
#  de_summary <- de %>% group_by(cluster) %>%
#    summarize(
#      n_up = sum(logFC > log2(1.5) & adj.P.Val < 0.05),
#      n_down = sum(-logFC > log2(1.5) & adj.P.Val < 0.05)
#    )
#  de_summary_file <- glue("{out_dir}/de_summary_case-vs-control.tsv")
#  data.table::fwrite(de_summary, de_summary_file, sep = "\t")
#
#  p <- ggplot(
#    de_summary %>% pivot_longer(cols = c("n_up", "n_down")) %>%
#      mutate(value = ifelse(name == "n_down", -value, value))
#  ) +
#    aes(x = value, y = cluster, fill = name) +
#    geom_colh() +
#    scale_y_discrete(limits = rev(levels(naturalfactor(de_summary$cluster)))) +
#    scale_fill_manual(
#      # values = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[c(9,3)],
#      values = okabe(2),
#      guide = "none"
#    ) +
#    scale_x_continuous(labels = abs) +
#    labs(x = NULL, y = NULL)
#  my_ggsave(
#    "bars-summary",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 3,
#    height = 3,
#    units = "in",
#    dpi = 300
#  )
#
#  for (this_cluster in unique(pb_meta$cluster)) {
#    ix <- pb_meta$cluster == this_cluster
#    n_donors <- table(pb_meta$case[ix])
#    #
#    p <- plot_limma_volcano(de[de$cluster == this_cluster,]) +
#      labs(title = glue("Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
#    my_ggsave(
#      glue("volcano-cluster-{this_cluster}"),
#      out_dir = glue("{out_dir}"),
#      plot = p,
#      type = "pdf",
#      scale = 1, width = 5, height = 4, units = "in", dpi = 300
#    )
#    #
#    my_ids <- c(
#      (
#        de %>%
#        filter(cluster == this_cluster) %>%
#        filter(logFC > 0) %>%
#        top_n(n = 20, wt = -log10(P.Value) * logFC)
#      )$ensembl_id,
#      (
#        de %>%
#        filter(cluster == this_cluster) %>%
#        filter(logFC < 0) %>%
#        top_n(n = 20, wt = -log10(P.Value) * -logFC)
#      )$ensembl_id
#    )
#    # my_ids <- head((
#    #   de %>%
#    #   filter(cluster == this_cluster) %>%
#    #   filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
#    #   arrange(abs(logFC) + -log10(P.Value))
#    # )$ensembl_id, 50)
#    mat <- as.matrix(pb[my_ids,ix])
#    a_col <- as.data.frame(pb_meta)[ix,]
#    a_col <- left_join(a_col, unique(a1$obs[,c("donor","drug")]), by = "donor")
#    rownames(a_col) <- colnames(mat)
#    a_col <- a_col[,-1]
#    h1 <- dendsort(hclust(dist(mat), method = "complete"))
#    h2 <- dendsort(hclust(dist(t(mat)), method = "complete"))
#    a_col <- a_col[h2$order,]
#    mat <- mat[h1$order,h2$order]
#    #
#    a_colors <- list()
#    a_colors[["drug"]] <- okabe(8)[4:7]
#    names(a_colors[["drug"]]) <- c("None", "PD-1", "CTLA-4", "PD-1/CTLA-4")
#    a_colors[["case"]] <- okabe(2)
#    names(a_colors[["case"]]) <- c("Control", "Case")
#    #
#    heatmap_file <- glue("{out_dir}/heatmap-cluster-{this_cluster}.pdf")
#    message(heatmap_file)
#    pheatmap(
#      filename = heatmap_file,
#      width = ncol(mat) * 0.2 + 3,
#      height = nrow(mat) * 0.1 + 3,
#      mat = mat,
#      cluster_col = FALSE,
#      cluster_row = FALSE,
#      hclust_method = "average",
#      show_colnames = FALSE,
#      show_rownames = TRUE,
#      scale = "row",
#      # color = rev(scico::scico(palette = "roma", n = 20)),
#      color = colorRampPalette(rev(RColorBrewer::brewer.pal("RdBu", n = 11)))(51),
#      border_color = NA,
#      labels_row = ensembl_to_symbol[rownames(mat)],
#      # labels_col = labels_col,
#      annotation_col = a_col[,c("drug","case")],
#      annotation_colors = a_colors,
#      main = glue("Cluster {this_cluster}")
#    )
#  }
#
#  plot_de_by_cluster <- function(d) {
#    ggplot() +
#      annotate(
#        geom = "rect",
#        xmin = -Inf,
#        xmax = Inf,
#        ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#        ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#        alpha = 0.2
#      ) +
#      geom_vline(xintercept = 0, size = 0.3) +
#      geom_point(
#        data = d,
#        mapping = aes(x = logFC, y = Gene)
#      ) +
#      geom_errorbarh(
#        data = d,
#        mapping = aes(xmin = CI.L, xmax = CI.R, y = Gene),
#        height = 0
#      ) +
#      scale_x_continuous(
#        breaks = seq(-10, 10, by = 1),
#        labels = function(x) fractional::fractional(2 ^ x)
#      ) +
#      scale_y_continuous(
#        expand = c(0, 0),
#        breaks = seq(1, max(d$Gene)),
#        labels = levels(d$GeneName)
#      ) +
#      annotation_logticks(sides = "b", size = 0.3) +
#      expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#      facet_grid(~ cluster)
#    # facet_wrap(~ cluster, scales = "free_x")
#  }
#
#  #
#  my_genes <- unique((
#    de %>%
#    group_by(cluster) %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    top_n(wt = -log10(P.Value), n = 5) %>%
#    select(GeneName, Gene, logFC, P.Value, adj.P.Val, cluster) %>%
#    ungroup() %>%
#    group_by(GeneName) %>%
#    mutate(min_p = min(P.Value)) %>%
#    ungroup() %>%
#    arrange(min_p)
#  )$GeneName)
#  #
#  d <- de %>% filter(GeneName %in% my_genes)
#  d$GeneName <- factor(d$GeneName, rev(my_genes))
#  d$Gene <- as.integer(d$GeneName)
#  d$cluster <- naturalfactor(d$cluster)
#  #
#  p <- plot_de_by_cluster(d) +
#    theme(
#      axis.title.y = element_blank(),
#      axis.text.y = element_text(face = "italic"),
#      panel.grid.major.x = element_line()
#    )
#  my_ggsave(
#    "lanes_case-vs-control",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = length(unique(de$cluster)) * 2.8 + 2,
#    height = length(my_genes) * 0.2 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  my_clusters <- c(5, 6, 7, 8, 9)
#  #
#  d <- de %>% filter(ensembl_id %in% my_ids, cluster %in% my_clusters)
#  d$GeneName <- factor(d$GeneName)
#  d$Gene <- as.integer(d$GeneName)
#  d$cluster <- naturalfactor(d$cluster)
#  #
#  p <- plot_de_by_cluster(d) +
#    theme(
#      axis.title.y = element_blank(),
#      axis.text.y = element_text(face = "italic"),
#      panel.grid.major.x = element_line()
#    )
#  my_ggsave(
#    glue('lanes_case-vs-control-cluster-{paste(my_clusters, collapse = "-")}'),
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = length(unique(d$cluster)) * 2.8 + 2,
#    height = length(unique(d$GeneName)) * 0.25 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#
#  # d <- rbindlist(lapply(my_ids, function(my_id) {
#  #   d <- a1$obs %>% select(cell, cluster, donor, case) %>% as_tibble
#  #   d$gene <- a1$counts[my_id,]
#  #   d <- d %>%
#  #     filter(cluster %in% c(4, 6, 8, 9)) %>%
#  #     group_by(cluster, donor, case) %>%
#  #     summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#  #   d$ensembl_id <- my_id
#  #   d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#  #   d
#  # }))
#  # d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  # d$Gene <- as.integer(d$GeneName)
#  #
#  # p <- ggplot() +
#  #   annotate(
#  #     geom = "rect",
#  #     xmin = 0,
#  #     xmax = Inf,
#  #     ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#  #     ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#  #     alpha = 0.2
#  #   ) +
#  #   geom_point(
#  #     data = d,
#  #     mapping = aes(x = percent, y = Gene, fill = case, group = case),
#  #     # position = position_quasirandom(dodge.width = 0.5, width = 0.3, groupOnX = FALSE),
#  #     shape = 21, size = 3, stroke = 0.2
#  #   ) +
#  #   scale_fill_manual(values = okabe(8), guide = "none") +
#  #   scale_x_log10() +
#  #   scale_y_continuous(
#  #     expand = c(0, 0),
#  #     breaks = seq(1, max(d$Gene)),
#  #     labels = levels(d$GeneName)
#  #   ) +
#  #   annotation_logticks(sides = "b", size = 0.3) +
#  #   expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#  #   labs(x = "Percent", y = NULL) +
#  #   facet_wrap(~ cluster, ncol = length(unique(d$cluster)))
#  # my_ggsave(
#  #   "percent",
#  #   out_dir = out_dir,
#  #   type = "pdf",
#  #   plot = p,
#  #   scale = 1,
#  #   width = length(unique(d$cluster)) * 1.5 + 1,
#  #   height = length(my_ids) * 0.3 + 1,
#  #   units = "in",
#  #   dpi = 300
#  # )
#
#  d <- rbindlist(lapply(my_ids, function(my_id) {
#    d <- a1$obs %>% select(cell, cluster, case) %>% as_tibble
#    d$gene <- a1$counts[my_id,]
#    d <- d %>%
#      filter(cluster %in% my_clusters) %>%
#      group_by(cluster, case) %>%
#      summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#    d$ensembl_id <- my_id
#    d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#    d
#  }))
#  d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  d$Gene <- as.integer(d$GeneName)
#  #
#  p <- ggplot() +
#    annotate(
#      geom = "rect",
#      xmin = -Inf,
#      xmax = Inf,
#      ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#      ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#      alpha = 0.2
#    ) +
#    geom_vline(xintercept = seq(0, max(d$percent), by = 10), size = 0.3, color = "white") +
#    geom_colh(
#      data = d,
#      mapping = aes(x = percent, y = Gene, fill = case, group = case),
#      position = position_dodgev(),
#      size = 0.3
#    ) +
#    # geom_point(
#    #   data = d,
#    #   mapping = aes(x = percent, y = Gene, fill = case, group = case),
#    #   shape = 21, size = 3, stroke = 0.3
#    # ) +
#    # scale_x_log10() +
#    # annotation_logticks(sides = "b", size = 0.3) +
#    scale_fill_manual(values = okabe(8), guide = "none") +
#    scale_y_continuous(
#      expand = c(0, 0),
#      breaks = seq(1, max(d$Gene)),
#      labels = levels(d$GeneName)
#    ) +
#    expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#    labs(x = "Percent of cells with expression", y = NULL) +
#    facet_wrap(~ cluster, ncol = length(unique(d$cluster))) +
#    theme(axis.text.y = element_text(face = "italic"))
#  my_ggsave(
#    "percent",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    # width = 5,
#    width = length(unique(d$cluster)) * 1.7 + 1,
#    height = length(my_ids) * 0.3 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  #d <- rbindlist(lapply(my_ids, function(my_id) {
#  #  d <- a1$obs %>% select(cell, cluster, donor, case) %>% as_tibble
#  #  d$gene <- a1$counts[my_id,]
#  #  d <- d %>%
#  #    group_by(cluster, donor, case) %>%
#  #    summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#  #  d$ensembl_id <- my_id
#  #  d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#  #  d
#  #}))
#  #d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  ##
#  #mat <- dcast(d, ensembl_id ~ cluster + case + donor, value.var = "percent")
#  #mat_rownames <- mat[[1]]
#  #mat <- as.matrix(mat[,2:ncol(mat)])
#  #rownames(mat) <- mat_rownames
#  #mat <- mat[my_ids,]
#  #a_col <- as.data.frame(str_split_fixed(colnames(mat), "_", 3))
#  #colnames(a_col) <- c("cluster", "case", "donor")
#  #rownames(a_col) <- colnames(mat)
#  ##
#  #a_colors <- list()
#  #a_colors[["cluster"]] <- mpn65[1:length(unique(a_col$cluster))]
#  #names(a_colors[["cluster"]]) <- as.character(naturalsort(unique(a_col$cluster)))
#  #a_colors[["case"]] <- okabe(2)
#  #names(a_colors[["case"]]) <- c("Control", "Case")
#  ##
#  #heatmap_file <- glue("{out_dir}/heatmap-cluster-percent.pdf")
#  #pheatmap(
#  #  filename = heatmap_file,
#  #  width = ncol(mat) * 0.02 + 2,
#  #  height = nrow(mat) * 0.1 + 2,
#  #  mat = mat,
#  #  cluster_col = FALSE,
#  #  cluster_row = FALSE,
#  #  # hclust_method = "average",
#  #  show_colnames = FALSE,
#  #  # color = rev(scico::scico(palette = "davos", n = 20)),
#  #  # scale = "row",
#  #  # color = rev(scico::scico(palette = "roma", n = 20)),
#  #  color = rev(scico::scico(palette = "davos", n = 20)),
#  #  border_color = NA,
#  #  labels_row = ensembl_to_symbol[rownames(mat)],
#  #  # labels_col = labels_col,
#  #  annotation_col = a_col[,c("cluster","case")],
#  #  annotation_colors = a_colors
#  #)
#
#}


# Nuclei
# Full figures with donor-level volcano, cluster-level boxplot, lanes
########################################################################

#{
#
#  analysis_name <- "n3_2"
#  params <- list(
#    min_cells_in_cluster = 50,
#    min_percent_of_cells_with_gene = 5
#  )
#  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
#  print_status(glue("Reading {a1_file}"))
#  a1 <- qread(a1_file)
#  print_status(glue("done"))
#  my_leiden <- "leiden0.933"
#  # my_leiden <- "leiden1.08"
#  a1$obs$leiden <- a1$obs[[my_leiden]]
#  a1$obs$cluster <- a1$obs[[my_leiden]]
#  a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
#  a1$de <- presto::wilcoxauc(a1$log2cpm, a1$obs$leiden)
#	#
#	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/de-case-vs-control"))
#	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#	#
#	sample_info <- janitor::clean_names(read_excel(
#			path = "data/luoma_villani_combined_colitis_samples2.xlsx",
#			sheet = 3
#	))
#	donor_to_sex <- unlist(split(sample_info$sex, sample_info$donor))
#  #
#	a1$obs$sex <- donor_to_sex[a1$obs$donor]
#	a1$obs$drug <- factor(a1$obs$drug, c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4"))
#  ## Skip cell clusters that have too few cells.
#  #exclude_clusters <- (
#  #  a1$obs %>% count(cluster) %>% filter(n < params[["min_cells_in_cluster"]])
#  #)$leiden
#  ##
#  #a1$obs <- a1$obs[!a1$obs$cluster %in% exclude_clusters,]
#  #
#  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))
#
#
#  # Pseudobulk at the donor level
#  ########################################################################
#  y <- with(a1$obs, model.matrix(~ 0 + factor(donor)))
#  y <- as(y, "dgCMatrix")
#  # y <- sweep(y, 2, colSums(y), "/") # means
#  pb <- as(a1$counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#  #
#  pb_meta <- data.frame(donor = colnames(pb))
#  pb_meta <- as_tibble(pb_meta)
#  pb_meta %<>%
#    mutate(
#      donor = str_replace(donor, "factor\\(donor\\)", "")
#    )
#  colnames(pb) <- str_replace(colnames(pb), "factor\\(donor\\)", "")
#  pb_meta <- left_join(
#    pb_meta,
#    a1$obs %>%
#      select(
#        donor, case
#      ) %>%
#      group_by(donor, case) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#
#  # x <- table(a1$obs$cluster)
#  # min_percent <- 100 * min(x) * 0.25 / sum(x)
#  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
#  keep_ens <- rownames(pb)[rowMeans(pb) > 0.5]
#  #
#  des1 <- with(pb_meta, model.matrix(
#    ~ case
#  ))
#  ob <- as.matrix(pb[keep_ens,])
#  fit1 <- lmFit(object = ob, design = des1)
#  fit1 <- eBayes(fit1)
#  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#  de_donor <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE) %>%
#    rename(Gene = ID)
#  de_donor$ensembl_id <- rownames(de_donor)
#  de_donor <- as_tibble(de_donor)
#  de_donor_file <- glue("{out_dir}/de_donor_case-vs-control.tsv.gz")
#  message(glue("Writing {de_donor_file}"))
#  data.table::fwrite(de_donor, de_donor_file, sep = "\t")
#  write_de_xlsx(de_donor, file.path(out_dir, "de_donor_case-vs-control.xlsx"))
#
#  n_donors <- table(pb_meta$case)
#  p <- plot_limma_volcano(de_donor, n_text = 20) +
#    labs(title = glue("Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
#  my_ggsave(
#    glue("volcano-de_donor"),
#    out_dir = glue("{out_dir}"),
#    plot = p,
#    type = "pdf",
#    scale = 1, width = 5, height = 4, units = "in", dpi = 300
#  )
#
#  my_ids <- (
#    de_donor %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    arrange(-logFC)
#  )$ensembl_id[1:20]
#  # my_ids <- de_donor$ensembl_id[1:15]
#  d <- as_tibble(reshape2::melt(as.matrix(pb)[my_ids,]))
#  colnames(d) <- c("ensembl_id", "donor", "value")
#  d$GeneName <- factor(
#    as.character(unname(ensembl_to_symbol[as.character(d$ensembl_id)])),
#    rev(as.character(unname(ensembl_to_symbol[as.character(my_ids)])))
#  )
#  d$Gene <- as.integer(d$GeneName)
#  d <- left_join(d, pb_meta, by = "donor")
#  #
#  p <- ggplot(d) +
#      annotate(
#        geom = "rect",
#        xmin = -Inf,
#        xmax = Inf,
#        ymin = seq(from = 1, to = length(my_ids), by = 2) - 0.5,
#        ymax = seq(from = 1, to = length(my_ids), by = 2) + 0.5,
#        alpha = 0.2
#      ) +
#      # geom_vline(xintercept = 0, size = 0.3) +
#      geom_point(
#        data = d,
#        mapping = aes(x = value, y = Gene, group = case, fill = case),
#        position = position_quasirandom(dodge.width = 0.5, width = 0.2, groupOnX = FALSE),
#        shape = 21, size = 2, stroke = 0.2
#      ) +
#      scale_fill_manual(values = okabe(8)) +
#      # scale_x_continuous(
#      #   breaks = seq(-10, 10, by = 1),
#      #   labels = function(x) fractional::fractional(2 ^ x)
#      # ) +
#      scale_y_continuous(
#        expand = c(0, 0),
#        breaks = seq(1, max(d$Gene)),
#        labels = levels(d$GeneName)
#      ) +
#      guides(fill = guide_legend(reverse = TRUE, title = NULL, override.aes = list(size = 5))) +
#      # annotation_logticks(sides = "b", size = 0.3) +
#      expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#      labs(x = bquote("Log"[2]~"CPM"), y = NULL) +
#      theme(
#        axis.text.y = element_text(face = "italic"),
#        legend.margin = margin(0, 0, 0, 0),
#        legend.box.margin = margin(0, 0, 0, 0)
#      )
#  my_ggsave(
#    "lanes-de_donor",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 5,
#    height = length(my_ids) * 0.3 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  my_ids <- (
#    de_donor %>%
#    filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
#    arrange(-logFC)
#  )$ensembl_id#[1:80]
#  mat <- as.matrix(pb[my_ids,])
#  a_col <- as.data.frame(pb_meta)
#  a_col <- left_join(a_col, unique(a1$obs[,c("donor","drug")]), by = "donor")
#  rownames(a_col) <- colnames(mat)
#  a_col <- a_col[,-1]
#  h1 <- dendsort(hclust(dist(mat), method = "complete"))
#  h2 <- dendsort(hclust(dist(t(mat)), method = "complete"))
#  a_col <- a_col[h2$order,]
#  mat <- mat[h1$order,h2$order]
#  #
#  a_colors <- list()
#  a_colors[["drug"]] <- okabe(8)[4:7]
#  names(a_colors[["drug"]]) <- c("None", "PD-1", "CTLA-4", "PD-1/CTLA-4")
#  a_colors[["case"]] <- okabe(2)
#  names(a_colors[["case"]]) <- c("Control", "Case")
#  #
#  heatmap_file <- glue("{out_dir}/heatmap-de_donor-all.pdf")
#  pheatmap(
#    filename = heatmap_file,
#    width = ncol(mat) * 0.2 + 3,
#    height = nrow(mat) * 0.005 + 3,
#    mat = mat,
#    cluster_col = FALSE,
#    cluster_row = FALSE,
#    hclust_method = "average",
#    show_colnames = FALSE,
#    show_rownames = FALSE,
#    scale = "row",
#    # color = rev(scico::scico(palette = "roma", n = 20)),
#    color = rev(RColorBrewer::brewer.pal("RdBu", n = 11)),
#    border_color = NA,
#    labels_row = ensembl_to_symbol[rownames(mat)],
#    # labels_col = labels_col,
#    annotation_col = a_col[,c("drug","case")],
#    annotation_colors = a_colors,
#    title = glue("{nrow(mat)} genes")
#  )
#
#  my_ids <- (
#    de_donor %>%
#    filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
#    arrange(-logFC * -log10(P.Value))
#  )$ensembl_id[1:50]
#  mat <- as.matrix(pb[my_ids,])
#  a_col <- as.data.frame(pb_meta)
#  a_col <- left_join(a_col, unique(a1$obs[,c("donor","drug")]), by = "donor")
#  rownames(a_col) <- colnames(mat)
#  a_col <- a_col[,-1]
#  h1 <- dendsort(hclust(dist(mat), method = "complete"))
#  h2 <- dendsort(hclust(dist(t(mat)), method = "complete"))
#  a_col <- a_col[h2$order,]
#  mat <- mat[h1$order,h2$order]
#  #
#  a_colors <- list()
#  a_colors[["drug"]] <- okabe(8)[4:7]
#  names(a_colors[["drug"]]) <- c("None", "PD-1", "CTLA-4", "PD-1/CTLA-4")
#  a_colors[["case"]] <- okabe(2)
#  names(a_colors[["case"]]) <- c("Control", "Case")
#  #
#  heatmap_file <- glue("{out_dir}/heatmap-de_donor.pdf")
#  pheatmap(
#    filename = heatmap_file,
#    width = ncol(mat) * 0.2 + 3,
#    height = nrow(mat) * 0.1 + 3,
#    mat = mat,
#    cluster_col = FALSE,
#    cluster_row = FALSE,
#    hclust_method = "average",
#    show_colnames = FALSE,
#    show_rownames = TRUE,
#    scale = "row",
#    # color = rev(scico::scico(palette = "roma", n = 20)),
#    color = rev(RColorBrewer::brewer.pal("RdBu", n = 11)),
#    border_color = NA,
#    labels_row = ensembl_to_symbol[rownames(mat)],
#    # labels_col = labels_col,
#    annotation_col = a_col[,c("drug","case")],
#    annotation_colors = a_colors
#  )
#
#  for (my_id in my_ids) {
#    p <- plot_hexgene(
#      x = a1$obs$UMAP1,
#      y = a1$obs$UMAP2,
#      z = as.numeric(a1$log2cpm[my_id,]),
#      group = a1$obs$case
#    ) +
#    facet_wrap(~ group) +
#    labs(title = ensembl_to_symbol[my_id]) +
#    theme(
#      panel.spacing = unit(0.5, "lines"),
#      plot.title = element_text(face = "italic")
#    )
#    my_ggsave(
#      glue("umap-{safe(ensembl_to_symbol[my_id])}"),
#      out_dir = out_dir,
#      type = "pdf",
#      plot = p,
#      scale = 0.8,
#      width = 7,
#      height = 5,
#      units = "in",
#      dpi = 300
#    )
#  }
#
#  # UMAPs for best markers for each cluster
#  ########################################################################
#  a1$de <- a1$de %>%
#    mutate(symbol = ensembl_to_symbol[feature])
#  a1$de_top <- a1$de %>%
#    # dplyr::filter(!group %in% exclude_clusters) %>%
#    dplyr::group_by(group) %>%
#    dplyr::filter(pct_in > 10) %>%
#    # dplyr::filter(group %in% c(12, 14, 16, 17, 19, 25, 29)) %>%
#    dplyr::mutate(
#      rank1 = rank(100 * auc - pct_out),
#      rank2 = rank(logFC * (pct_in - pct_out)),
#    ) %>%
#    dplyr::top_n(n = 5, wt = (rank1 + rank2))
#  for (this_cluster in unique(a1$de_top$group)) {
#    ix_cluster <- a1$de_top$group == this_cluster
#    my_ids <- a1$de_top$feature[ix_cluster]
#    for (my_id in my_ids) {
#      p <- plot_hexgene(
#        x = a1$obs$UMAP1,
#        y = a1$obs$UMAP2,
#        z = as.numeric(a1$log2cpm[my_id,])
#      ) +
#      labs(title = ensembl_to_symbol[my_id]) +
#      theme(
#        panel.spacing = unit(0.5, "lines"),
#        plot.title = element_text(face = "italic")
#      )
#      my_ggsave(
#        glue("umap-{safe(ensembl_to_symbol[my_id])}"),
#        out_dir = glue("{out_dir}/cluster-{this_cluster}"),
#        type = "pdf",
#        plot = p,
#        scale = 0.8,
#        width = 4,
#        height = 5,
#        units = "in",
#        dpi = 300
#      )
#    }
#  }
#
#  
#  # Pseudobulk at the cluster level
#  ########################################################################
#  y <- with(a1$obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
#  y <- as(y, "dgCMatrix")
#  # y <- sweep(y, 2, colSums(y), "/") # means
#  pb <- as(a1$counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#  #
#  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
#  colnames(pb_meta) <- c("cluster", "donor")
#  pb_meta <- as_tibble(pb_meta)
#  pb_meta %<>%
#    mutate(
#      cluster = str_replace(cluster, "factor\\(cluster\\)", ""),
#      donor = str_replace(donor, "factor\\(donor\\)", "")
#    )
#  pb_meta <- left_join(
#    pb_meta,
#    a1$obs %>%
#      select(
#        donor, case
#      ) %>%
#      group_by(donor, case) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#
#  # x <- table(a1$obs$cluster)
#  # min_percent <- 100 * min(x) * 0.25 / sum(x)
#  # keep_ens <- a1$counts_stats$gene[a1$counts_stats$percent >= min_percent]
#  keep_ens <- rownames(pb)[rowMeans(pb) > 0.5]
#  #
#  de <- rbindlist(
#    pblapply(sort(unique(pb_meta$cluster)), function(this_cluster) {
#      ix_cluster <- pb_meta$cluster == this_cluster
#      des1 <- with(pb_meta[ix_cluster,], model.matrix(
#        ~ case
#      ))
#      ob <- as.matrix(pb[keep_ens,ix_cluster])
#      fit1 <- lmFit(object = ob, design = des1)
#      fit1 <- eBayes(fit1)
#      fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#      res <- topTable(fit1, coef = "caseCase", number = 1e6, confint = TRUE)
#      res$cluster <- this_cluster
#      res <- res %>% rename(Gene = ID)
#      res$ensembl_id <- rownames(res)
#      return(res)
#    })
#  )
#  de$GeneName <- de$Gene
#  de$cluster <- naturalfactor(de$cluster)
#  de_file <- glue("{out_dir}/de_case-vs-control.tsv.gz")
#  data.table::fwrite(de, de_file, sep = "\t")
#
#  de_summary <- de %>% group_by(cluster) %>%
#    summarize(
#      n_up = sum(logFC > log2(1.5) & adj.P.Val < 0.05),
#      n_down = sum(-logFC > log2(1.5) & adj.P.Val < 0.05)
#    )
#  de_summary_file <- glue("{out_dir}/de_summary_case-vs-control.tsv")
#  data.table::fwrite(de_summary, de_summary_file, sep = "\t")
#
#  for (this_cluster in unique(pb_meta$cluster)) {
#    ix <- pb_meta$cluster == this_cluster
#    n_donors <- table(pb_meta$case[ix])
#    #
#    p <- plot_limma_volcano(de[de$cluster == this_cluster,]) +
#      labs(title = glue("Case (n={n_donors['Case']}) vs Control (n={n_donors['Control']})"))
#    my_ggsave(
#      glue("volcano-cluster-{this_cluster}"),
#      out_dir = glue("{out_dir}"),
#      plot = p,
#      type = "pdf",
#      scale = 1, width = 5, height = 4, units = "in", dpi = 300
#    )
#    #
#    my_ids <- c(
#      (
#        de %>%
#        filter(cluster == this_cluster) %>%
#        filter(logFC > 0) %>%
#        top_n(n = 20, wt = -log10(P.Value) * logFC)
#      )$ensembl_id,
#      (
#        de %>%
#        filter(cluster == this_cluster) %>%
#        filter(logFC < 0) %>%
#        top_n(n = 20, wt = -log10(P.Value) * -logFC)
#      )$ensembl_id
#    )
#    # my_ids <- head((
#    #   de %>%
#    #   filter(cluster == this_cluster) %>%
#    #   filter(abs(logFC) > log2(2), adj.P.Val < 0.05) %>%
#    #   arrange(abs(logFC) + -log10(P.Value))
#    # )$ensembl_id, 50)
#    mat <- as.matrix(pb[my_ids,ix])
#    a_col <- as.data.frame(pb_meta)[ix,]
#    a_col <- left_join(a_col, unique(a1$obs[,c("donor","drug")]), by = "donor")
#    rownames(a_col) <- colnames(mat)
#    a_col <- a_col[,-1]
#    h1 <- dendsort(hclust(dist(mat), method = "complete"))
#    h2 <- dendsort(hclust(dist(t(mat)), method = "complete"))
#    a_col <- a_col[h2$order,]
#    mat <- mat[h1$order,h2$order]
#    #
#    a_colors <- list()
#    a_colors[["drug"]] <- okabe(8)[4:7]
#    names(a_colors[["drug"]]) <- c("None", "PD-1", "CTLA-4", "PD-1/CTLA-4")
#    a_colors[["case"]] <- okabe(2)
#    names(a_colors[["case"]]) <- c("Control", "Case")
#    #
#    heatmap_file <- glue("{out_dir}/heatmap-cluster-{this_cluster}.pdf")
#    message(heatmap_file)
#    pheatmap(
#      filename = heatmap_file,
#      width = ncol(mat) * 0.2 + 3,
#      height = nrow(mat) * 0.1 + 3,
#      mat = mat,
#      cluster_col = FALSE,
#      cluster_row = FALSE,
#      hclust_method = "average",
#      show_colnames = FALSE,
#      show_rownames = TRUE,
#      scale = "row",
#      # color = rev(scico::scico(palette = "roma", n = 20)),
#      color = colorRampPalette(rev(RColorBrewer::brewer.pal("RdBu", n = 11)))(51),
#      border_color = NA,
#      labels_row = ensembl_to_symbol[rownames(mat)],
#      # labels_col = labels_col,
#      annotation_col = a_col[,c("drug","case")],
#      annotation_colors = a_colors,
#      main = glue("Cluster {this_cluster}")
#    )
#  }
#
#  p <- ggplot(
#    de_summary %>% pivot_longer(cols = c("n_up", "n_down")) %>%
#      mutate(value = ifelse(name == "n_down", -value, value))
#  ) +
#    aes(x = value, y = cluster, fill = name) +
#    geom_colh() +
#    scale_y_discrete(limits = rev(levels(naturalfactor(de_summary$cluster)))) +
#    scale_fill_manual(
#      # values = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[c(9,3)],
#      values = okabe(2),
#      guide = "none"
#    ) +
#    scale_x_continuous(labels = abs) +
#    labs(x = NULL, y = NULL)
#  my_ggsave(
#    "bars-summary",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = 3,
#    height = length(unique(de_summary$cluster)) * 0.25,
#    units = "in",
#    dpi = 300
#  )
#
#  plot_de_by_cluster <- function(d) {
#    ggplot() +
#      annotate(
#        geom = "rect",
#        xmin = -Inf,
#        xmax = Inf,
#        ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#        ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#        alpha = 0.2
#      ) +
#      geom_vline(xintercept = 0, size = 0.3) +
#      geom_point(
#        data = d,
#        mapping = aes(x = logFC, y = Gene)
#      ) +
#      geom_errorbarh(
#        data = d,
#        mapping = aes(xmin = CI.L, xmax = CI.R, y = Gene),
#        height = 0
#      ) +
#      scale_x_continuous(
#        breaks = seq(-10, 10, by = 1),
#        labels = function(x) fractional::fractional(2 ^ x)
#      ) +
#      scale_y_continuous(
#        expand = c(0, 0),
#        breaks = seq(1, max(d$Gene)),
#        labels = levels(d$GeneName)
#      ) +
#      annotation_logticks(sides = "b", size = 0.3) +
#      expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#      facet_grid(~ cluster)
#    # facet_wrap(~ cluster, scales = "free_x")
#  }
#
#  #
#  my_genes <- unique((
#    de %>%
#    group_by(cluster) %>%
#    filter(logFC > log2(2), adj.P.Val < 0.05) %>%
#    top_n(wt = -log10(P.Value), n = 5) %>%
#    select(GeneName, Gene, logFC, P.Value, adj.P.Val, cluster) %>%
#    ungroup() %>%
#    group_by(GeneName) %>%
#    mutate(min_p = min(P.Value)) %>%
#    ungroup() %>%
#    arrange(min_p)
#  )$GeneName)
#  #
#  d <- de %>% filter(GeneName %in% my_genes)
#  d$GeneName <- factor(d$GeneName, rev(my_genes))
#  d$Gene <- as.integer(d$GeneName)
#  d$cluster <- naturalfactor(d$cluster)
#  #
#  p <- plot_de_by_cluster(d) +
#    theme(
#      axis.title.y = element_blank(),
#      axis.text.y = element_text(face = "italic"),
#      panel.grid.major.x = element_line()
#    )
#  my_ggsave(
#    "lanes_case-vs-control",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = length(unique(de$cluster)) * 2.8 + 2,
#    height = length(my_genes) * 0.2 + 1,
#    units = "in",
#    dpi = 300,
#    limitsize = FALSE
#  )
#
#  my_clusters <- c(5, 6, 7, 8, 9)
#  #
#  d <- de %>% filter(ensembl_id %in% my_ids, cluster %in% my_clusters)
#  d$GeneName <- factor(d$GeneName)
#  d$Gene <- as.integer(d$GeneName)
#  d$cluster <- naturalfactor(d$cluster)
#  #
#  p <- plot_de_by_cluster(d) +
#    theme(
#      axis.title.y = element_blank(),
#      axis.text.y = element_text(face = "italic"),
#      panel.grid.major.x = element_line()
#    )
#  my_ggsave(
#    glue('lanes_case-vs-control-cluster-{paste(my_clusters, collapse = "-")}'),
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    width = length(unique(d$cluster)) * 2.8 + 2,
#    height = length(unique(d$GeneName)) * 0.25 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#
#  # d <- rbindlist(lapply(my_ids, function(my_id) {
#  #   d <- a1$obs %>% select(cell, cluster, donor, case) %>% as_tibble
#  #   d$gene <- a1$counts[my_id,]
#  #   d <- d %>%
#  #     filter(cluster %in% c(4, 6, 8, 9)) %>%
#  #     group_by(cluster, donor, case) %>%
#  #     summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#  #   d$ensembl_id <- my_id
#  #   d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#  #   d
#  # }))
#  # d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  # d$Gene <- as.integer(d$GeneName)
#  #
#  # p <- ggplot() +
#  #   annotate(
#  #     geom = "rect",
#  #     xmin = 0,
#  #     xmax = Inf,
#  #     ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#  #     ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#  #     alpha = 0.2
#  #   ) +
#  #   geom_point(
#  #     data = d,
#  #     mapping = aes(x = percent, y = Gene, fill = case, group = case),
#  #     # position = position_quasirandom(dodge.width = 0.5, width = 0.3, groupOnX = FALSE),
#  #     shape = 21, size = 3, stroke = 0.2
#  #   ) +
#  #   scale_fill_manual(values = okabe(8), guide = "none") +
#  #   scale_x_log10() +
#  #   scale_y_continuous(
#  #     expand = c(0, 0),
#  #     breaks = seq(1, max(d$Gene)),
#  #     labels = levels(d$GeneName)
#  #   ) +
#  #   annotation_logticks(sides = "b", size = 0.3) +
#  #   expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#  #   labs(x = "Percent", y = NULL) +
#  #   facet_wrap(~ cluster, ncol = length(unique(d$cluster)))
#  # my_ggsave(
#  #   "percent",
#  #   out_dir = out_dir,
#  #   type = "pdf",
#  #   plot = p,
#  #   scale = 1,
#  #   width = length(unique(d$cluster)) * 1.5 + 1,
#  #   height = length(my_ids) * 0.3 + 1,
#  #   units = "in",
#  #   dpi = 300
#  # )
#
#  d <- rbindlist(lapply(my_ids, function(my_id) {
#    d <- a1$obs %>% select(cell, cluster, case) %>% as_tibble
#    d$gene <- a1$counts[my_id,]
#    d <- d %>%
#      filter(cluster %in% my_clusters) %>%
#      group_by(cluster, case) %>%
#      summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#    d$ensembl_id <- my_id
#    d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#    d
#  }))
#  d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  d$Gene <- as.integer(d$GeneName)
#  #
#  p <- ggplot() +
#    annotate(
#      geom = "rect",
#      xmin = -Inf,
#      xmax = Inf,
#      ymin = seq(from = 1, to = max(d$Gene), by = 2) - 0.5,
#      ymax = seq(from = 1, to = max(d$Gene), by = 2) + 0.5,
#      alpha = 0.2
#    ) +
#    geom_vline(xintercept = seq(0, max(d$percent), by = 10), size = 0.3, color = "white") +
#    geom_colh(
#      data = d,
#      mapping = aes(x = percent, y = Gene, fill = case, group = case),
#      position = position_dodgev(),
#      size = 0.3
#    ) +
#    # geom_point(
#    #   data = d,
#    #   mapping = aes(x = percent, y = Gene, fill = case, group = case),
#    #   shape = 21, size = 3, stroke = 0.3
#    # ) +
#    # scale_x_log10() +
#    # annotation_logticks(sides = "b", size = 0.3) +
#    scale_fill_manual(values = okabe(8), guide = "none") +
#    scale_y_continuous(
#      expand = c(0, 0),
#      breaks = seq(1, max(d$Gene)),
#      labels = levels(d$GeneName)
#    ) +
#    expand_limits(y = c(0.5, max(d$Gene) + 0.5)) +
#    labs(x = "Percent of cells with expression", y = NULL) +
#    facet_wrap(~ cluster, ncol = length(unique(d$cluster))) +
#    theme(axis.text.y = element_text(face = "italic"))
#  my_ggsave(
#    "percent",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1,
#    # width = 5,
#    width = length(unique(d$cluster)) * 1.7 + 1,
#    height = length(my_ids) * 0.3 + 1,
#    units = "in",
#    dpi = 300
#  )
#
#  #d <- rbindlist(lapply(my_ids, function(my_id) {
#  #  d <- a1$obs %>% select(cell, cluster, donor, case) %>% as_tibble
#  #  d$gene <- a1$counts[my_id,]
#  #  d <- d %>%
#  #    group_by(cluster, donor, case) %>%
#  #    summarize(percent = 100 * sum(gene > 0) / n(), .groups = "keep")
#  #  d$ensembl_id <- my_id
#  #  d$GeneName <- ensembl_to_symbol[d$ensembl_id]
#  #  d
#  #}))
#  #d$GeneName <- factor(d$GeneName, rev(ensembl_to_symbol[my_ids]))
#  ##
#  #mat <- dcast(d, ensembl_id ~ cluster + case + donor, value.var = "percent")
#  #mat_rownames <- mat[[1]]
#  #mat <- as.matrix(mat[,2:ncol(mat)])
#  #rownames(mat) <- mat_rownames
#  #mat <- mat[my_ids,]
#  #a_col <- as.data.frame(str_split_fixed(colnames(mat), "_", 3))
#  #colnames(a_col) <- c("cluster", "case", "donor")
#  #rownames(a_col) <- colnames(mat)
#  ##
#  #a_colors <- list()
#  #a_colors[["cluster"]] <- mpn65[1:length(unique(a_col$cluster))]
#  #names(a_colors[["cluster"]]) <- as.character(naturalsort(unique(a_col$cluster)))
#  #a_colors[["case"]] <- okabe(2)
#  #names(a_colors[["case"]]) <- c("Control", "Case")
#  ##
#  #heatmap_file <- glue("{out_dir}/heatmap-cluster-percent.pdf")
#  #pheatmap(
#  #  filename = heatmap_file,
#  #  width = ncol(mat) * 0.02 + 2,
#  #  height = nrow(mat) * 0.1 + 2,
#  #  mat = mat,
#  #  cluster_col = FALSE,
#  #  cluster_row = FALSE,
#  #  # hclust_method = "average",
#  #  show_colnames = FALSE,
#  #  # color = rev(scico::scico(palette = "davos", n = 20)),
#  #  # scale = "row",
#  #  # color = rev(scico::scico(palette = "roma", n = 20)),
#  #  color = rev(scico::scico(palette = "davos", n = 20)),
#  #  border_color = NA,
#  #  labels_row = ensembl_to_symbol[rownames(mat)],
#  #  # labels_col = labels_col,
#  #  annotation_col = a_col[,c("cluster","case")],
#  #  annotation_colors = a_colors
#  #)
#
#  # Markers
#  ########################################################################
#  out_dir <- glue("results/a20/{analysis_name}/figures/composition-case-vs-control")
#  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#	a1$obs$drug <- factor(a1$obs$drug, c("None", "CTLA-4", "PD-1", "PD-1/CTLA-4"))
#  a1$obs <- as_tibble(a1$obs)
#  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))
#  #
#  source("R/functions/composition.R")
#  print_status("do_masc()")
#  res3 <- do_masc(
#    a1$obs,
#    form1 = "is_cluster ~ 1 + (1|donor)",
#    form2 = "is_cluster ~ 1 + case + (1|donor)",
#    mc.cores = 8
#  )
#  print_status("done")
#  res3_coef <- summarize_masc(res3)
#  res3_coef$value <- str_replace(res3_coef$value, "case", "")
#
#  my_res <- 0.1
#  my_leiden <- run_leiden(
#    adj        = a1$knn$simil_cells,
#    resolution = my_res,
#    iterations = 20,
#    seed       = 42
#  )
#  a1$obs[[glue("leiden{signif(my_res, 3)}")]] <- my_leiden
#  n_donors <- length(unique(a1$obs$donor))
#  p1 <- plot_scattermore(
#    x = a1$obs$UMAP1,
#    y = a1$obs$UMAP2,
#    group = a1$obs$leiden0.1,
#    group_colors = mpn65,
#    pixels = 1000,
#    alpha = 0.35
#  ) +
#  labs(
#    title = glue(
#      "{length(unique(a1$obs$leiden))} clusters of {comma(nrow(a1$obs))} cells from {n_donors} donor{ifelse(n_donors > 1, 's', '')}"
#    )
#  )
#  my_ggsave(
#    "umap-clusters-leiden0.1",
#    out_dir = out_dir,
#    plot = p1,
#    type = "pdf",
#    scale = 1, width = 5.5, height = 5, units = "in", dpi = 300
#  )
#
#  a1$obs[["major_cluster"]] <- "Epithelial cells"
#  a1$obs[["major_cluster"]][a1$obs$leiden %in% c(7, 21, 22, 23)] <- "Mesenchymal cells"
#  a1$obs[["major_cluster"]][a1$obs$leiden %in% c(8)] <- "Stem cells"
#  n_donors <- length(unique(a1$obs$donor))
#  p1 <- plot_scattermore(
#    x            = a1$obs$UMAP1,
#    y            = a1$obs$UMAP2,
#    group        = a1$obs$major_cluster,
#    group_colors = mpn65,
#    pixels       = 1000,
#    alpha        = 0.35,
#    group_labels = FALSE,
#    group_legend = TRUE
#  )
#  my_ggsave(
#    "umap-major_clusters",
#    out_dir = out_dir,
#    plot = p1,
#    type = "pdf",
#    scale = 1, width = 7.5, height = 3, units = "in", dpi = 300
#  )
#
#  set.seed(1)
#  ix_rand <- sample(nrow(a1$obs))
#  p1 <- plot_scattermore(
#    x            = a1$obs$UMAP1[ix_rand],
#    y            = a1$obs$UMAP2[ix_rand],
#    group        = pub_ids[as.character(a1$obs$donor[ix_rand])],
#    # group        = as.character(a1$obs$donor[ix_rand]),
#    group_colors = mpn65,
#    pixels       = 1000,
#    alpha        = 0.35,
#    group_labels = FALSE,
#    group_legend = TRUE
#  )
#  my_ggsave(
#    "umap-donor",
#    out_dir = out_dir,
#    plot = p1,
#    type = "pdf",
#    scale = 1, width = 9.5, height = 4, units = "in", dpi = 300
#  )
#
#  x <- a1$obs %>%
#    group_by(donor) %>%
#    count(major_cluster)
#  x <- x %>% group_by(donor) %>% mutate(pct = 100 * n / sum(n))
#  x$pub <- pub_ids[as.character(x$donor)]
#  # x$pub <- x$donor
#  x$pub <- naturalfactor(x$pub)
#  x$pub <- factor(x$pub, rev(levels(x$pub)))
#  p1 <- ggplot(x) +
#    aes(x = pct, y = pub, fill = major_cluster) +
#    geom_colh() +
#    scale_fill_manual(values = mpn65[c(1,2,3,4,5,7,6)]) +
#    scale_x_continuous(expand = c(0, 0)) +
#    theme(
#      legend.position = "none"
#    ) +
#    labs(x = "Percent", y = NULL)
#  x_pub <- x %>% group_by(pub) %>% summarize(n = sum(n))
#  p2 <- ggplot(x) +
#    geom_colh(
#      data = x,
#      mapping = aes(x = n, y = pub, fill = major_cluster)
#    ) +
#    geom_text(
#      data = x_pub,
#      mapping = aes(x = n, y = pub, label = comma(n, accuracy = 1)),
#      hjust = 0, nudge_x = 100
#    ) +
#    scale_x_continuous(
#      expand = expansion(mult = c(0.01, 0.25)), labels = label_number_si()
#    ) +
#    scale_fill_manual(name = NULL, values = mpn65[c(1,2,3,4,5,7,6)]) +
#    labs(x = "Cells", y = NULL) +
#    theme(
#      axis.text.y = element_blank(),
#      axis.ticks.y = element_blank()
#    )
#  my_ggsave(
#    "bars-major_cluster",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p1 + p2 + plot_annotation(title = "Composition of each donor's cells"),
#    scale = 0.75,
#    width = 11, height = 10,
#    units = "in", dpi = 300
#  )
#
#  # - Epithelial: EPCAM, CDH1
#  # - STEM: LGR5, SMOC2
#  # - Mesenchymal: VIM, COL3A1
#
#  my_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in%
#    c("EPCAM", "CDH1", "LGR5", "SMOC2", "VIM", "COL3A1")])
#  for (my_id in my_ids) {
#    my_symbol <- ensembl_to_symbol[my_id]
#    p <- plot_hexgene(
#      x = a1$obs$UMAP1,
#      y = a1$obs$UMAP2,
#      z = as.numeric(a1$log2cpm[my_id,]),
#      bins         = 101,
#      palette      = "davos",
#      direction    = -1,
#      use_quantile = FALSE,
#      text = FALSE,
#      legend = FALSE
#    ) +
#    labs(title = my_symbol) +
#    theme(
#      plot.title = element_text(face = "italic")
#    )
#    my_ggsave(
#      # glue("umap-{safe(my_symbol)}-{my_contrast}"),
#      glue("umap-{safe(my_symbol)}"),
#      out_dir = glue("{out_dir}/umap"),
#      type = "pdf",
#      plot = p,
#      scale = 0.8,
#      width = 3.5,
#      height = 3,
#      units = "in",
#      dpi = 300
#    )
#  }
#
#  #
#  tsv_file1 <- glue("{out_dir}/masc_1-case-complete.tsv")
#  fwrite(res3_coef, tsv_file1, sep = "\t")
#  #
#  for (key in safe(unique(res3_coef$value))) {
#    tsv_file2 <- glue("{out_dir}/masc_1-case-{safe(key)}.tsv")
#    xx <- res3_coef %>% filter(value == key)
#    xx %>%
#      dplyr::mutate(
#        OR = exp(est),
#        OR_2.5 = exp(est_low),
#        OR_97.5 = exp(est_high),
#        lrt_fdr = p.adjust(lrt_p, method = "fdr")
#      ) %>%
#      dplyr::select(cluster, OR, OR_2.5, OR_97.5, p, lrt_p, lrt_fdr) %>%
#      dplyr::mutate_if(is.numeric, signif, 5) %>%
#      dplyr::arrange(lrt_p) %>%
#      write_tsv(tsv_file2)
#  }
#
#  n_clusters <- length(unique(a1$obs$leiden))
#  #
#  res3_coef$cluster <- factor(
#    as.character(res3_coef$cluster),
#    (
#      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
#    )$cluster
#  )
#  p1 <- ggplot(res3_coef %>% filter(value == "drugPD-1/CTLA-4")) +
#    aes(x = est, y = cluster) +
#    geom_vline(xintercept = 0, size = 0.3) +
#    ggforestplot::geom_stripes() +
#    geom_point(
#    ) +
#    geom_errorbarh(
#      mapping = aes(xmin = est_low, xmax = est_high),
#      height = 0
#    ) +
#    scale_x_continuous(
#      breaks = seq(-10, 10, by = 1),
#      labels = function(x) fractional::fractional(2 ^ x)
#    ) +
#    annotation_logticks(sides = "b", size = 0.3) +
#    expand_limits(y = c(0.5, max(as.integer(res3_coef$cluster)) + 0.5)) +
#    # facet_grid(~ value) +#, scales = "free_x") +
#    labs(title = "Case vs Control", y = NULL, x = "OR")
#  my_ggsave(
#    "lanes_1-case",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p1,
#    scale = 1,
#    width = 4,
#    height = length(unique(res3_coef$cluster)) * 0.25 + 1,
#    units = "in", dpi = 300
#  )
#  #
#  source("R/functions/okabe-ito.R")
#  a1$obs$cluster <- factor(
#    as.character(a1$obs$leiden),
#    (
#      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
#    )$cluster
#  )
#	p2 <- plot_composition_h(
#		d = a1$obs, fill = case, group = donor, x = cluster,
#    legend.position = "right"
#	) +
#	labs(title = glue("Composition of each donor (n = {length(unique(a1$obs$donor))})"))
#	my_ggsave(
#		"compositionh",
#		out_dir = out_dir,
#		type = "pdf",
#		plot = p2,
#		scale = 1, height = n_clusters * 0.5 + 1, width = 6, units = "in", dpi = 300
#	)
#  #
#  p_both <- (
#    p1 + 
#      labs(title = "Case vs Control", y = NULL, x = "OR")
#  ) + (
#    p2 +
#      theme(
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank()
#      )
#  )
#  p_both <- p_both + plot_layout(widths = c(2, 3))
#	my_ggsave(
#		"lanes_1-case-compositionh",
#		out_dir = out_dir,
#		type = "pdf",
#		plot = p_both,
#		scale = 1, height = n_clusters * 0.5 + 1, width = 10, units = "in", dpi = 300
#	)
#
#  a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))
#  a1$de <- presto::wilcoxauc(a1$log2cpm, a1$obs$leiden)
#
#  # Heatmap of best markers for each cluster
#  ########################################################################
#  cluster_colors <- mpn65
#  names(cluster_colors) <- as.character(seq_along(cluster_colors))
#  # Top markers for each cluster
#  a1$de <- a1$de %>%
#    mutate(symbol = ensembl_to_symbol[feature])
#  a1$de_top <- a1$de %>%
#    # dplyr::filter(!group %in% exclude_clusters) %>%
#    dplyr::group_by(group) %>%
#    dplyr::filter(pct_in > 10) %>%
#    # dplyr::filter(group %in% c(12, 14, 16, 17, 19, 25, 29)) %>%
#    dplyr::mutate(
#      rank1 = rank(100 * auc - pct_out),
#      rank2 = rank(logFC * (pct_in - pct_out)),
#    ) %>%
#    dplyr::top_n(n = 5, wt = (rank1 + rank2))
#  # these_genes <- unique(a1$de_top$feature)
#  #
#  # Subset the table
#  x <- a1$de %>% dplyr::filter(feature %in% these_genes)
#  x <- as.data.frame(dcast.data.table(
#    data = as.data.table(x), formula = symbol ~ group, value.var = "auc"
#  ))
#  rownames(x) <- x$symbol
#  x$symbol <- NULL
#  #
#  my_ens <- c()
#  for (my_cluster in rev(levels(res3_coef$cluster))) {
#    add_ens <- (
#      a1$de_top %>% filter(group == my_cluster) %>% arrange(rank1 + rank2)
#    )$feature %>% head(3)
#    my_ens <- c(my_ens, setdiff(add_ens, my_ens))
#  }
#  # # Order the heatmap
#  set.seed(1)
#  # xo <- seriation::seriate(as.matrix(x), method = "BEA_TSP")
#  xd <- a1$de %>% dplyr::filter(feature %in% my_ens)
#  # xd$symbol <- factor(xd$symbol, rownames(x)[xo[[1]]])
#  # xd$group <- factor(xd$group, colnames(x)[xo[[2]]])
#  xd$group <- factor(xd$group, levels(res3_coef$cluster))
#  xd$symbol <- factor(xd$symbol, ensembl_to_symbol[my_ens])
#  #
#  # Plot
#  p3 <- ggplot(xd) +
#  geom_tile(
#    aes(x = symbol, y = group, fill = auc)
#  ) +
#  scale_x_discrete(position = "t", name = NULL, expand = c(0, 0)) +
#  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
#  scale_fill_gradientn(
#    colors = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)),
#    limits = c(0, 1),
#    name = "AUC", guide = guide_colorbar(barwidth = 15), breaks = pretty_breaks(5)
#  ) +
#  theme(
#    axis.text.x = element_text(size = 14, face = "italic", angle = 60, hjust = 0),
#    legend.position = "bottom"
#  )
#  p4 <- ggplot(xd) +
#  geom_tile(
#    aes(x = 1, y = group, fill = group)
#  ) +
#  scale_x_discrete(position = "b", name = NULL, expand = c(0, 0)) +
#  scale_y_discrete(position = "l", name = NULL, expand = c(0, 0)) +
#  scale_fill_manual(values = cluster_colors, guide = "none")
#  p <- (
#    p4 + theme(plot.margin = margin(r = 0))
#  ) + (
#    p3 + theme(
#      axis.text.y = element_blank(),
#      axis.title.y = element_blank(),
#      axis.ticks.y = element_blank(),
#      plot.margin = margin(l = 0)
#    )
#  ) + plot_layout(widths = c(1, length(these_genes)))
#  fig_width  <- length(these_genes) * 0.3
#  fig_height <- length(unique(a1$de_top$group)) * 0.4 + 3
#  my_ggsave(
#    slug = "heatmap-top-markers_h",
#    out_dir = out_dir,
#    type = "pdf",
#    plot = p,
#    scale = 1, width = fig_width, height = fig_height, units = "in", dpi = 300
#  )
#  #
#  p <- (
#    p4 + theme(plot.margin = margin(r = 0))
#  ) + (
#    p3 + theme(
#      axis.text.y = element_blank(),
#      axis.title.y = element_blank(),
#      axis.ticks.y = element_blank(),
#      plot.margin = margin(l = 0)
#    )
#  ) + (
#    p1 + labs(title = "Case vs Control", y = NULL, x = "OR")
#  ) + (
#    p2 +
#      theme(
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank()
#      )
#  ) + plot_layout(widths = c(0.25, length(these_genes) * 0.2, 1.5, 3))
#	my_ggsave(
#		"full",
#		out_dir = out_dir,
#		# type = c("pdf", "png"),
#    type = "pdf",
#		plot = p,
#		scale = 1, height = n_clusters * 0.5 + 1, width = 8 + length(these_genes) * 0.3,
#    units = "in", dpi = 300
#	)
#
#}


# Drug analysis logistic OR, per-donor composition
# Figure 5
########################################################################

# analysis_name <- "a12_4_4_t4_cd8_1_2"
# analysis_name <- "a12_4_4_t4_cd4_2_2"
# analysis_name <- "a12_4_4_m3_2"
# analysis_name <- "a12_4_4_b5_1_3"
# analyses <- c(
#   "a12_4_4_t4_cd8_1_2",
#   "a12_4_4_t4_cd4_2_2",
#   "a12_4_4_m3_2",
#   "a12_4_4_b5_1_3",
#   "n3_2"
# )

for (analysis_name in analyses) {

  params <- list(
    min_cells_in_cluster = 50,
    min_percent_of_cells_with_gene = 5
  )
  a1_file <- as.character(glue("results/a20/{analysis_name}/data/{analysis_name}.qs"))
  print_status(glue("Reading {a1_file}"))
  a1 <- qread(a1_file)
  print_status(glue("done"))
  a1$obs$leiden <- as.character(a1$obs$leiden)
  a1$obs$cluster <- a1$obs$leiden
  out_dir <- glue("results/a20/{analysis_name}/figures")
  source("R/functions/composition.R")
  source("R/plot-composition.R")
	#
	out_dir <- as.character(glue("results/a20/{analysis_name}/figures/drug"))
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
  my_drug <- c("None", "PD-1", "PD-1/CTLA-4")
  a1$obs <- a1$obs %>% filter(drug %in% my_drug)
  a1$obs$drug <- factor(a1$obs$drug, my_drug)
  ix <- with(
    a1$obs,
    case == "Case" | (case == "Control" & drug == "PD-1") |
      (case == "Control" & drug == "None")
  )
  a1$obs <- a1$obs[ix,]
  # table(a1$obs$drug, a1$obs$case)
  # Skip cell clusters that have too few cells.
  exclude_clusters <- (
    a1$obs %>% count(leiden) %>% filter(n < params[["min_cells_in_cluster"]])
  )$leiden
  #
  a1$obs <- a1$obs[!a1$obs$leiden %in% exclude_clusters,]
  #
  a1$obs <- as_tibble(a1$obs)
  a1$obs$case <- factor(a1$obs$case, c("Control", "Case"))
  #
  a1$obs$case_drug <- str_replace_all(
    sprintf("%s_%s", a1$obs$case, a1$obs$drug),
    "[^[:alnum:] ]", ""
  )
  a1$obs$case_drug <- factor(a1$obs$case_drug, c("ControlNone", "ControlPD1", "CasePD1", "CasePD1CTLA4"))
  #
  a1$obs$case_drug3 <- as.character(a1$obs$case_drug)
  a1$obs$case_drug3[a1$obs$case_drug3 == "ControlNone"] <- "Control"
  a1$obs$case_drug3[a1$obs$case_drug3 == "ControlPD1"] <- "Control"
  a1$obs$case_drug3 <- factor(a1$obs$case_drug3, c("Control", "CasePD1", "CasePD1CTLA4"))
  table(a1$obs$case_drug3)
  a1$obs$sex[a1$obs$donor == "SIC_33"] <- "M"

  a1$counts <- a1$counts[,a1$obs$cell]
  a1$logcpm <- do_log1p_cpm(a1$counts)

  source("R/functions/composition.R")
  print_status("do_masc()")
  res3_file <- glue("{out_dir}/composition.qs")
  if (!file.exists(res3_file)) {
    res3 <- do_masc(
      a1$obs[a1$obs$case == "Case",],
      form1 = "is_cluster ~ 1 + sex + (1|donor)",
      form2 = "is_cluster ~ 1 + drug + sex + (1|donor)",
      # a1$obs,
      # form1 = "is_cluster ~ 1 + case + (1|donor)",
      # form2 = "is_cluster ~ 1 + case + drug + (1|donor)",
      mc.cores = 8
    )
    qsave(res3, res3_file)
  } else {
    res3 <- qread(res3_file)
  }
  print_status("done")
  res3_coef <- summarize_masc(res3)
  res3_coef$value <- str_replace(res3_coef$value, "case", "")
  res3_coef <- res3_coef %>% group_by(value) %>%
    mutate(lrt_fdr = p.adjust(lrt_p, method = "fdr")) %>%
    ungroup

  #
  tsv_file1 <- glue("{out_dir}/masc_1-case-drug-complete.tsv")
  fwrite(res3_coef, tsv_file1, sep = "\t")
  #
  for (key in unique(res3_coef$value)) {
    tsv_file2 <- glue("{out_dir}/masc_{safe(key)}.tsv")
    xx <- res3_coef %>% filter(value == key)
    xx %>%
      dplyr::mutate(
        OR = exp(est),
        OR_2.5 = exp(est_low),
        OR_97.5 = exp(est_high),
        lrt_fdr = p.adjust(lrt_p, method = "fdr")
      ) %>%
      dplyr::select(cluster, OR, OR_2.5, OR_97.5, p, lrt_p, lrt_fdr) %>%
      dplyr::mutate(label = glue("{safe(key)} cluster {cluster} {signif(OR,2)}-fold (95% CI {signif(OR_2.5,2)} to {signif(OR_97.5,2)}, P = {str_replace(signif(lrt_p, 2), '-0', '-')})")) %>%
      dplyr::mutate_if(is.numeric, signif, 5) %>%
      dplyr::arrange(lrt_p) %>%
      write_tsv(tsv_file2)
  }

  res3_coef$cluster <- factor(
    as.character(res3_coef$cluster),
    (
      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
    )$cluster
  )
  n_clusters <- length(unique(a1$obs$leiden))
  d_error <- res3_coef %>% filter(value == "drugPD-1/CTLA-4") %>%
    mutate(
      est      = log2(exp(est)),
      est_low  = log2(exp(est_low)),
      est_high = log2(exp(est_high))
    )
  vlines <- seq(-8, 8, by = 1)
  vlines <- vlines[vlines > min(d_error$est_low)]
  vlines <- vlines[vlines < max(d_error$est_high)]
  vlines <- vlines[vlines != 0]
  p_error <- ggplot(d_error) +
    aes(x = est, y = cluster) +
    ggforestplot::geom_stripes() +
		geom_vline(
			xintercept = vlines,
			size = 0.3, color = "white"
		) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_point(aes(color = lrt_fdr < 0.05)) +
    geom_errorbarh(
      mapping = aes(xmin = est_low, xmax = est_high, color = lrt_fdr < 0.05),
      height = 0
    ) +
    geom_text(
      data = d_error %>% filter(lrt_fdr < 0.05),
      mapping = aes(
        x = ifelse(est > 0, -Inf, Inf),
        hjust = ifelse(est > 0, 0, 1),
        y = cluster,
        label = str_replace(sprintf(" %s ", signif(lrt_p, 1)), "-0", "-")
      ),
      size = 5, color = "grey30"
    ) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = seq(-8, 8, by = 1),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    # annotation_logticks(sides = "b", size = 0.3) +
    # expand_limits(y = c(0.5, max(as.integer(res3_coef$cluster)) + 0.5)) +
    # facet_grid(~ value) +#, scales = "free_x") +
    labs(title = "Case PD-1/CTLA-4 vs Case PD-1", y = NULL, x = "OR")
    theme(plot.background = element_blank())
  # p_error <- ggplot(res3_coef %>% filter(value == "drugPD-1/CTLA-4")) +
  #   aes(x = est, y = cluster) +
  #   geom_vline(xintercept = 0, size = 0.3) +
  #   ggforestplot::geom_stripes() +
  #   geom_point(
  #   ) +
  #   geom_errorbarh(
  #     mapping = aes(xmin = est_low, xmax = est_high),
  #     height = 0
  #   ) +
  #   scale_x_continuous(
  #     breaks = seq(-10, 10, by = 1),
  #     labels = function(x) fractional::fractional(2 ^ x)
  #   ) +
  #   annotation_logticks(sides = "b", size = 0.3) +
  #   expand_limits(y = c(0.5, max(as.integer(res3_coef$cluster)) + 0.5)) +
  #   # facet_grid(~ value) +#, scales = "free_x") +
  #   labs(title = "Case PD-1/CTLA-4 vs Case PD-1", y = NULL, x = "OR")
  my_ggsave(
    "lanes_est",
    out_dir = out_dir,
    type = "pdf",
    plot = p_error,
    scale = 1,
    width = 4,
    height = length(unique(res3_coef$cluster)) * 0.25 + 1,
    units = "in", dpi = 300
  )

  #
  n_clusters <- length(unique(a1$obs$leiden))
  source("R/functions/okabe-ito.R")
  a1$obs$cluster <- factor(
    as.character(a1$obs$leiden),
    (
      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
    )$cluster
  )
	p_comp <- plot_composition_h(
		d = a1$obs, fill = case_drug, group = donor, x = cluster,
    legend.position = "right"
	) +
	labs(title = glue("Composition of each donor (n = {length(unique(a1$obs$donor))})"))
	my_ggsave(
		"compositionh",
		out_dir = out_dir,
		type = "pdf",
		plot = p_comp,
		scale = 1, height = n_clusters * 0.5 + 1, width = 6, units = "in", dpi = 300
	)
  #
  p_both <- (
    p_error + labs(title = "", y = NULL, x = "OR")
  ) + (
    p_comp +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
  my_widths <- c(1, 3)
  if (analysis_name == "blood2_myeloid5") {
    my_widths <- c(1, 1)
  }
  p_both <- p_both + plot_layout(widths = my_widths)
  my_width <- 10
	my_ggsave(
		"lanes_compositionh",
		out_dir = out_dir,
		type = "pdf",
		plot = p_both,
		scale = 1, height = n_clusters * 0.5 + 1, width = my_width, units = "in", dpi = 300
	)

  #
  n_clusters <- length(unique(a1$obs$leiden))
  source("R/functions/okabe-ito.R")
  a1$obs$cluster <- factor(
    as.character(a1$obs$leiden),
    (
      res3_coef %>% select(cluster, lrt_p) %>% unique %>% arrange(-lrt_p)
    )$cluster
  )
	p_comp3 <- plot_composition_h(
		d = a1$obs, fill = case_drug3, group = donor, x = cluster,
    legend.position = "right"
	) +
	labs(title = glue("Composition of each donor (n = {length(unique(a1$obs$donor))})"))
	my_ggsave(
		"compositionh-3",
		out_dir = out_dir,
		type = "pdf",
		plot = p_comp3,
		scale = 1, height = n_clusters * 0.5 + 1, width = 6, units = "in", dpi = 300
	)
  #
  p_both <- (
    p_error + labs(title = "", y = NULL, x = "OR")
  ) + (
    p_comp3 +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
  p_both <- p_both + plot_layout(widths = c(1, 3))
	my_ggsave(
		"lanes_compositionh-3",
		out_dir = out_dir,
		type = "pdf",
		plot = p_both,
		scale = 1, height = n_clusters * 0.5 + 1, width = 10, units = "in", dpi = 300
	)

  source("R/plot-composition.R")
  res3_drug <- res3_coef %>% filter(value == "drugPD-1/CTLA-4") %>% arrange(lrt_p)
  res3_drug$fdr <- p.adjust(res3_drug$lrt_p, method = "fdr")
  my_clusters <- as.character(( res3_drug %>% filter(fdr < 0.05))$cluster)
  if (analysis_name == "n3_2") {
    my_clusters <- as.character(( res3_drug %>% filter(lrt_p < 0.06))$cluster)
  }
  comp <- get_composition(a1$obs, fill = case_drug3, group = donor, x = cluster)
  res3_drug_top <- res3_drug %>% filter(cluster %in% my_clusters)
  res3_drug_top$cluster <- factor(res3_drug_top$cluster)
  # p_error_top <- ggplot(res3_drug_top) +
  #   aes(x = est, y = cluster) +
  #   geom_vline(xintercept = 0, size = 0.3) +
  #   ggforestplot::geom_stripes() +
  #   geom_point(
  #   ) +
  #   geom_errorbarh(
  #     mapping = aes(xmin = est_low, xmax = est_high),
  #     height = 0
  #   ) +
  #   scale_x_continuous(
  #     breaks = seq(-10, 10, by = 1),
  #     labels = function(x) fractional::fractional(2 ^ x)
  #   ) +
  #   annotation_logticks(sides = "b", size = 0.3) +
  #   expand_limits(y = c(0.5, max(as.integer(res3_drug_top$cluster)) + 0.5)) +
  #   # facet_grid(~ value) +#, scales = "free_x") +
  #   labs(title = "Case PD-1/CTLA-4 vs Case PD-1", y = NULL, x = "OR")
  res3_drug_top <- res3_drug_top %>%
    mutate(
      est      = log2(exp(est)),
      est_low  = log2(exp(est_low)),
      est_high = log2(exp(est_high))
    )
  vlines <- seq(-8, 8, by = 1)
  vlines <- vlines[vlines > min(res3_drug_top$est_low)]
  vlines <- vlines[vlines < max(res3_drug_top$est_high)]
  vlines <- vlines[vlines != 0]
  p_error_top <- ggplot(res3_drug_top) +
    aes(x = est, y = cluster) +
    ggforestplot::geom_stripes() +
		geom_vline(
			xintercept = vlines,
			size = 0.3, color = "white"
		) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_point(aes(color = lrt_fdr < 0.05)) +
    geom_errorbarh(
      mapping = aes(xmin = est_low, xmax = est_high, color = lrt_fdr < 0.05),
      height = 0
    ) +
    geom_text(
      # data = res3_drug_top %>% filter(lrt_fdr < 0.05),
      data = res3_drug_top,
      mapping = aes(
        x = ifelse(est > 0, -Inf, Inf),
        hjust = ifelse(est > 0, 0, 1),
        y = cluster,
        label = str_replace(sprintf(" %s ", signif(lrt_p, 1)), "-0", "-")
      ),
      size = 5, color = "grey30"
    ) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = seq(-8, 8, by = 1),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    labs(title = "Case PD-1/CTLA-4 vs Case PD-1", y = NULL, x = "OR")
    theme(plot.background = element_blank())
	p_comp3_top <- plot_composition_h(
		d = NULL,
    comp = comp %>% filter(cluster %in% my_clusters),
    fill = case_drug3,
    group = donor,
    x = cluster,
    legend.position = "right"
	) +
	labs(title = glue("Composition of each donor (n = {length(unique(a1$obs$donor))})"))
  p_both_top <- (
    p_error_top + labs(title = "", y = NULL, x = "OR")
  ) + (
    p_comp3_top +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
  my_widths <- c(1, 2)
  my_width <- 8
  if (analysis_name == "blood2_myeloid5") {
    my_widths <- c(1, 1)
    my_width <- 10
  }
  p_both_top <- p_both_top + plot_layout(widths = my_widths)
	my_ggsave(
		"compositionh-3_top",
		out_dir = out_dir,
		type = "pdf",
		plot = p_both_top,
		scale = 1, height = length(my_clusters) * 0.5 + 1, width = my_width, units = "in", dpi = 300
	)
  
  if (analysis_name == "n3_2") {
    # p <- plot_composition_h(
    #   d = NULL,
    #   fill = case_drug3,
    #   group = donor,
    #   x = cluster,
    #   legend.position = "right"
    # )
    p <- ggplot(comp %>% filter(cluster %in% c("13"))) +
      aes(x = freq, y = case_drug3, fill = case_drug3) +
      geom_boxplot(size = 0.3, outlier.shape = NA, alpha = 0.5) +
      geom_point(position = position_quasirandom(groupOnX = FALSE, width = 0.4),
        shape = 21, size = 3, stroke = 0.3) +
      scale_fill_manual(values = okabe(), guide = "none") +
      scale_y_discrete(position = "r") +
      scale_x_log10() +
      annotation_logticks(sides = "b", size = 0.3) +
      labs(x = NULL, y = NULL)
    my_ggsave(
      "compositionh-cluster13",
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1, height = 1.5, width = 4, units = "in", dpi = 300
    )
  }

  try({

    a1$obs$donor_leiden <- with(a1$obs, glue("{donor} {leiden}"))
    # ix <- a1$obs$donor_leiden %in% (a1$obs %>% count(donor_leiden) %>% filter(n > 5))$donor_leiden
    # ix <- ix & a1$obs$case == "Control"
    ix <- a1$obs$case == "Control"
    sum(ix)

    source("R/functions/composition.R")
    print_status("do_masc()")
    res4_file <- glue("{out_dir}/composition-control.qs")
    if (!file.exists(res4_file)) {
      # form1 <- "is_cluster ~ 1 + sex + (1|donor)"
      # form2 <- "is_cluster ~ 1 + drug + sex + (1|donor)"
      # if (str_detect(analysis_name, "^blood")) {
        form1 <- "is_cluster ~ 1 + (1|donor)"
        form2 <- "is_cluster ~ 1 + drug + (1|donor)"
      # }
      res4 <- do_masc(
        my_obs = a1$obs[ix,],
        form1 = form1,
        form2 = form2,
        mc.cores = 8
      )
      qsave(res4, res4_file)
    } else {
      res4 <- qread(res4_file)
    }
    print_status("done")
    res4_coef <- summarize_masc(res4)
    res4_coef$value <- str_replace(res4_coef$value, "case", "")
    res4_coef <- res4_coef %>% group_by(value) %>%
      mutate(lrt_fdr = p.adjust(lrt_p, method = "fdr")) %>%
      ungroup
    #
    tsv_file1 <- glue("{out_dir}/masc-control-complete.tsv")
    fwrite(res4_coef, tsv_file1, sep = "\t")
    #
    for (key in unique(res4_coef$value)) {
      tsv_file2 <- glue("{out_dir}/masc-control-{safe(key)}.tsv")
      xx <- res4_coef %>% filter(value == key)
      xx %>%
        dplyr::mutate(
          OR = exp(est),
          OR_2.5 = exp(est_low),
          OR_97.5 = exp(est_high),
          lrt_fdr = p.adjust(lrt_p, method = "fdr")
        ) %>%
        dplyr::select(cluster, OR, OR_2.5, OR_97.5, p, lrt_p, lrt_fdr) %>%
        dplyr::mutate(label = glue("{safe(key)} cluster {cluster} {signif(OR,2)}-fold (95% CI {signif(OR_2.5,2)} to {signif(OR_97.5,2)}, P = {str_replace(signif(lrt_p, 2), '-0', '-')})")) %>%
        dplyr::mutate_if(is.numeric, signif, 5) %>%
        dplyr::arrange(lrt_p) %>%
        write_tsv(tsv_file2)
    }

    exclude_clusters <- c(
      (
        a1$obs[ix,] %>% group_by(drug, leiden) %>%
        summarize(n_donors = length(unique(donor)), .groups = "drop") %>%
        arrange(n_donors) %>%
        filter(n_donors <= 1)
      )$leiden,
      (
      a1$obs[ix,] %>% group_by(drug, leiden) %>%
        group_by(leiden) %>% summarize(n_drug = length(unique(drug))) %>%
        filter(n_drug <= 1)
      )$leiden,
      (
        res4_coef %>% filter(value == "drugPD-1") %>% filter(sd > 3)
      )$cluster
    )

    res4_coef$cluster <- factor(
      as.character(res4_coef$cluster),
      (
        res4_coef %>% select(cluster, value, lrt_p) %>% filter(value == "drugPD-1") %>% arrange(-lrt_p)
      )$cluster
    )
    n_clusters <- length(unique(a1$obs$leiden))
    d_error <- res4_coef %>% filter(value == "drugPD-1") %>%
      mutate(
        est      = log2(exp(est)),
        est_low  = log2(exp(est_low)),
        est_high = log2(exp(est_high))
      )
    d_error$est[d_error$cluster %in% exclude_clusters] <- NA
    d_error$est_low[d_error$cluster %in% exclude_clusters] <- NA
    d_error$est_high[d_error$cluster %in% exclude_clusters] <- NA
    vlines <- seq(-8, 8, by = 1)
    vlines <- vlines[vlines > min(d_error$est_low)]
    vlines <- vlines[vlines < max(d_error$est_high)]
    vlines <- vlines[vlines != 0]
    p_error <- ggplot(d_error) +
      aes(x = est, y = cluster) +
      ggforestplot::geom_stripes() +
      geom_vline(
        xintercept = vlines,
        size = 0.3, color = "white"
      ) +
      geom_vline(xintercept = 0, size = 0.3) +
      geom_point(aes(color = lrt_fdr < 0.05)) +
      geom_errorbarh(
        mapping = aes(xmin = est_low, xmax = est_high, color = lrt_fdr < 0.05),
        height = 0
      ) +
      geom_text(
        data = d_error %>% filter(lrt_fdr < 0.05),
        mapping = aes(
          x = ifelse(est > 0, -Inf, Inf),
          hjust = ifelse(est > 0, 0, 1),
          y = cluster,
          label = str_replace(sprintf(" %s ", signif(lrt_p, 1)), "-0", "-")
        ),
        size = 5, color = "grey30"
      ) +
      scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "black"), guide = "none") +
      scale_x_continuous(
        breaks = seq(-8, 8, by = 1),
        labels = function(x) fractional::fractional(2 ^ x)
      ) +
      labs(title = "Control PD-1 vs Control None", y = NULL, x = "OR")
      theme(plot.background = element_blank())
    my_ggsave(
      "masc-control-errorbar",
      out_dir = out_dir,
      type = "pdf",
      plot = p_error,
      scale = 1,
      width = 4,
      height = length(unique(res4_coef$cluster)) * 0.25 + 1,
      units = "in", dpi = 300
    )

    #
    n_clusters <- length(unique(a1$obs$leiden))
    source("R/functions/okabe-ito.R")
    a1$obs$cluster <- factor(
      as.character(a1$obs$leiden),
      levels(res4_coef$cluster)
    )
    p_comp <- plot_composition_h(
      d = a1$obs[ix,], fill = case_drug, group = donor, x = cluster,
      legend.position = "right"
    ) +
    labs(title = glue("Composition of each donor (n = {length(unique(a1$obs$donor[ix]))})"))
    my_ggsave(
      "masc-control-compositionh",
      out_dir = out_dir,
      type = "pdf",
      plot = p_comp,
      scale = 1, height = n_clusters * 0.5 + 1, width = 6, units = "in", dpi = 300
    )
    #
    p_both <- (
      p_error + labs(title = "", y = NULL, x = "OR")
    ) + (
      p_comp +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
    )
    my_widths <- c(1, 2)
    if (analysis_name == "blood2_myeloid5") {
      my_widths <- c(1, 1)
    }
    p_both <- p_both + plot_layout(widths = my_widths)
    my_width <- 10
    my_ggsave(
      "masc-control-lanes",
      out_dir = out_dir,
      type = "pdf",
      plot = p_both,
      scale = 1, height = n_clusters * 0.3 + 1, width = my_width, units = "in", dpi = 300
    )
  })

  gene_pct <- 100 * rowSums(a1$counts > 0) / ncol(a1$counts)
  stopifnot(nrow(a1$obs) == ncol(a1$counts))
  y <- with(a1$obs, model.matrix(~ 0 + donor))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/")
  pb <- as(a1$counts %*% y, "dgCMatrix")
  pb <- pb[,colSums(pb) > 0]
  pb <- do_log1p_cpm(pb, median(colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "donor", "")
  # my_ens <- names(ensembl_to_symbol[ensembl_to_symbol == "LINC00278"])
  # sort(pb[my_ens,])
  #
  pb_meta <- tibble(donor = colnames(pb))
  #
  pb_meta <- left_join(
    pb_meta,
    a1$obs %>% group_by(case, case_drug, donor, drug, sex) %>% count() %>% ungroup,
    by = c("donor")
  )
  pb <- pb[,as.character(pb_meta$donor)]
  stopifnot(all(pb_meta$donor == colnames(pb)))
  pb_meta$case_drug3 <- as.character(pb_meta$case_drug)
  pb_meta$case_drug3[pb_meta$case_drug3 == "ControlNone"] <- "Control"
  pb_meta$case_drug3[pb_meta$case_drug3 == "ControlPD1"] <- "Control"
  pb_meta$case_drug3 <- factor(pb_meta$case_drug3, c("Control", "CasePD1", "CasePD1CTLA4"))
  #
  des1 <- with(pb_meta, model.matrix(
    ~ 0 + case_drug + sex
  ))
  des2 <- with(pb_meta, model.matrix(
    ~ 0 + case_drug3 + sex
  ))
  colnames(des1) <- str_replace_all(colnames(des1), "[^[:alnum:] ]", "")
  colnames(des1) <- str_replace(colnames(des1), "casedrug", "")
  cont <- makeContrasts(
    contrasts = list(
      "ControlPD1 - ControlNone",
      "CasePD1 - ControlPD1",
      "CasePD1CTLA4 - CasePD1"
    ),
    levels = des1
  )
  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
  stopifnot(nrow(des1) == nrow(pb_meta))
  #
  object <- as.matrix(pb)
  object <- object[rowMeans(object) > 0.5,]
  fit1 <- lmFit(object = object, design = des1)
  fit1 <- eBayes(fit1)
  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  #
  fit2 <- contrasts.fit(fit1, cont)
  fit2 <- eBayes(fit2)
  #
  res <- rbindlist(lapply(colnames(cont), function(this_contrast) {
    res <- topTable(fit2, coef = this_contrast, number = 1e6, confint = TRUE)
    res$contrast <- this_contrast
    colnames(res)[1] <- "Gene"
    res$ensembl_id <- rownames(res)
    return(res)
  }))
  res$percent <- gene_pct[res$ensembl_id]

  # res %>%
  #   filter(contrast == "CasePD1CTLA4vsCasePD1", percent > 2, adj.P.Val < 0.05) %>%
  #   arrange(P.Value) %>%
  #   head()

  # drug_auc <- rbindlist(lapply(colnames(cont), function(my_contrast) {
  #   ix_cells <- a1$obs$case_drug %in% as.character(str_split_fixed(my_contrast, "vs", 2))
  #   retval <- as.data.table(presto::wilcoxauc(a1$logcpm[,ix_cells], a1$obs$case_drug[ix_cells]))
  #   retval$contrast <- my_contrast
  #   retval
  # }))
  # x <- left_join(
  #   x = res,
  #   y = drug_auc %>% select(ensembl_id = feature, contrast, auc, pct_in, pct_out),
  #   by = c("ensembl_id", "contrast")
  # )
  #
  write_de_xlsx(
    res %>%
      mutate_if(is.numeric, signif, 4),
    glue("{out_dir}/de_donor.xlsx"),
    col = "contrast"
  )
  fwrite(res, glue("{out_dir}/de_donor.tsv.gz"), sep = "\t")
  #
  for (this_contrast in colnames(cont)) {
    p <- plot_limma_volcano(
      res %>% filter(percent > 2, contrast == this_contrast),
      n_text = 10, fdr = 0.1
    ) +
      labs(title = glue("{this_contrast}"))
    my_ggsave(
      glue("volcano-donor-{this_contrast}"),
      out_dir = glue("{out_dir}/volcano"),
      plot = p,
      type = "pdf",
      scale = 1, width = 5, height = 4, units = "in", dpi = 300
    )
  }

  contrast_pairs <- combn(colnames(cont), m = 2)
  for (i in seq(ncol(contrast_pairs))) {
    contrast1 <- contrast_pairs[1,i]
    contrast2 <- contrast_pairs[2,i]
    res2 <- as_tibble(inner_join(
      x = res %>% filter(contrast == contrast1),
      y = res %>% filter(contrast == contrast2),
      by = c("ensembl_id", "Gene")
    )) %>% filter(percent.x > 2) # %>% filter(adj.P.Val.x < 0.5 | adj.P.Val.y < 0.5)
    res2_text <- res2 %>% filter(
      (abs(logFC.x) > log2(1.5) & adj.P.Val.x < 0.1) |
      (abs(logFC.y) > log2(1.5) & adj.P.Val.y < 0.1)
    )
    p <- ggplot() +
      aes(logFC.x, logFC.y, label = Gene) +
      geom_hline(yintercept = 0, size = 0.3) +
      geom_vline(xintercept = 0, size = 0.3) +
      geom_abline(intercept = 0, slope = c(-1, 1), size = 0.3) +
      geom_point(data = res2, size = 0.5, alpha = 0.2) +
      geom_point(data = res2_text, size = 0.5, color = "red") +
      geom_text_repel(
        data = res2_text, size = 3, fontface = "italic"
      ) +
      labs(
        x = contrast1,
        y = contrast2,
        title = "Log2 Fold-changes, either contrast FC > 1.5 and FDR < 10%"
      )
    my_ggsave(
      slug = glue("point-donor-{contrast1}-{contrast2}"),
      out_dir = glue("{out_dir}"),
      type = "pdf",
      plot = p,
      scale = 1, units = "in", dpi = 300,
      width = 8,
      height = 6
    )
  }

  # my_ids <- c(
  #   "DUOX2", "DUOXA2", "CXCL1", "CXCL2", "ICAM1", "NOS2", "HIF1A", "ANXA1",
  #   "TNFAIP8", "SLC6A14", "ANO6"
  # )
  # my_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in% my_ids])

  for (my_contrast in colnames(cont)) {
    my_ids <- (
      res %>% filter(contrast == my_contrast) %>%
        top_n(n = 15, wt = -log10(P.Value) * abs(logFC))
    )$ensembl_id
    for (my_id in my_ids) {
      # my_ids <- names(ensembl_to_symbol)[which(ensembl_to_symbol %in% c("CD28", "ITGB2"))]
      my_symbol <- ensembl_to_symbol[my_id]
      ix_cells <- a1$obs$case_drug %in% as.character(str_split_fixed(my_contrast, "vs", 2))
      # ix_cells <- a1$obs$case_drug %in% unique(a1$obs$case_drug)
      n_case_drug <- length(unique(a1$obs$case_drug[ix_cells]))
      p <- plot_hexgene(
        x = a1$obs$UMAP1[ix_cells],
        y = a1$obs$UMAP2[ix_cells],
        z = as.numeric(a1$logcpm[my_id,ix_cells]),
        # palette = "batlow", # davos is better
        palette = "oslo",
        group = naturalfactor(append_n(a1$obs$case_drug[ix_cells]))
      ) +
      facet_wrap(~ group, ncol = n_case_drug) +
      labs(title = my_symbol) +
      theme(
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        # glue("umap-{safe(my_symbol)}-{my_contrast}"),
        glue("umap-{safe(my_symbol)}"),
        out_dir = glue("{out_dir}/umap-donor/{my_contrast}"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 3.5 * n_case_drug,
        height = 5,
        units = "in",
        dpi = 300
      )
      my_ggsave(
        # glue("umap-{safe(my_symbol)}-{my_contrast}"),
        glue("umap-{safe(my_symbol)}"),
        out_dir = glue("{out_dir}/umap-donor/{my_contrast}/slim"),
        type = "pdf",
        plot = p + theme(
          legend.position = "none",
          strip.text = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(face = "italic")
        ),
        scale = 0.8,
        width = 3.25 * n_case_drug,
        height = 3,
        units = "in",
        dpi = 300
      )
    }
  }

  # # CD8 T cells from tissue
  # my_symbols <- c("CD28", "ITGB2", "LINC00861", "GZMK", "KLRC2", "PXN", "SELL", "TNFRSF4")
  # my_ids <- ensembl_to_symbol[ensembl_to_symbol %in% my_symbols]
  # my_ids <- unname(setNames(names(my_ids), my_ids)[my_symbols])

  # Epithelial cells from tissue
  my_symbols <- c(
    "CXCL1", "DUOX2", "HIF1A", "ICAM1", "NOS2", "SLC6A14", "TNFAIP8", "ANO6", "GNLY"
  )
  my_ids <- ensembl_to_symbol[ensembl_to_symbol %in% my_symbols]
  my_ids <- unname(setNames(names(my_ids), my_ids)[my_symbols])

  pb_meta$case_drug3 <- as.character(pb_meta$case_drug)
  pb_meta$case_drug3[pb_meta$case_drug3 == "ControlNone"] <- "Control"
  pb_meta$case_drug3[pb_meta$case_drug3 == "ControlPD1"] <- "Control"
  pb_meta$case_drug3 <- factor(pb_meta$case_drug3, c("Control", "CasePD1", "CasePD1CTLA4"))
  #
  for (my_contrast in unique(res$contrast)) {
    # my_ids <- c(
    #   (
    #     res %>% filter(logFC > 0, contrast == my_contrast) %>%
    #       mutate(wt = -log10(P.Value) * abs(logFC)) %>%
    #       top_n(n = 15, wt = wt) %>%
    #       arrange(-wt)
    #   )$ensembl_id,
    #   (
    #     res %>% filter(logFC < 0, contrast == my_contrast) %>%
    #       mutate(wt = -log10(P.Value) * abs(logFC)) %>%
    #       top_n(n = 15, wt = wt) %>%
    #       arrange(-wt)
    #   )$ensembl_id
    # )
    x <- rbindlist(lapply(my_ids, function(this_ens) {
      cbind(pb_meta,
        ens = this_ens,
        gene = pb[this_ens,],
        symbol = ensembl_to_symbol[this_ens]
      )
    }))
    x$symbol <- factor(x$symbol, levels = ensembl_to_symbol[my_ids])
    my_res <- res %>% filter(contrast == my_contrast, ensembl_id %in% my_ids) %>%
      rename(symbol = Gene)
    my_res$symbol <- factor(my_res$symbol, levels(x$symbol))
    my_res$lab <- sprintf("%.2f-fold p=%s", 2^my_res$logFC, format.pval(my_res$P.Value, 1))
    p <- ggplot() +
      geom_quasirandom(
        data = x,
        mapping = aes(x = gene, y = case_drug, group = case_drug, color = case_drug),
        width = 0.4,
        size = 3,
        groupOnX = FALSE
      ) +
      geom_text(
        data = my_res,
        mapping = aes(label = lab),
        x = -Inf, y = Inf,
        hjust = -0.05, vjust = 1.2
      ) +
      # facet_wrap(~ symbol, scales = "free_x", ncol = length(my_ids) / 3) +
      facet_wrap(~ symbol, scales = "free_x", ncol = 1) +
      guides(color = "none") +
      scale_color_manual(values = c("grey60", "grey20", pals::okabe(8)[2:8])) +
      scale_y_discrete(position = "right") +
      scale_x_continuous(breaks = scales::pretty_breaks(3)) +
      labs(
        y = NULL, x = bquote("Log"[2]~"CPM"),
        title = glue("All cells: {my_contrast}")
      ) +
      theme(strip.text = element_text(face = "italic"), panel.spacing = unit(0.25, "lines"))
    my_ggsave(
      slug = glue("dots-donor-{my_contrast}"),
      out_dir = glue("{out_dir}/dots-donor"),
      type = "pdf",
      plot = p,
      scale = 1, units = "in", dpi = 300,
      # width = 1 + length(my_ids) * 0.8,
      # height = 7
      width = 4,
      height = length(my_ids) * 2,
      limitsize = FALSE
    )
    p <- ggplot() +
      geom_quasirandom(
        data = x,
        mapping = aes(x = gene, y = case_drug3, group = case_drug3, color = case_drug3),
        width = 0.4,
        size = 3,
        groupOnX = FALSE
      ) +
      geom_text(
        data = my_res,
        mapping = aes(label = lab),
        x = -Inf, y = Inf,
        hjust = -0.05, vjust = 1.2
      ) +
      # facet_wrap(~ symbol, scales = "free_x", ncol = length(my_ids) / 3) +
      facet_wrap(~ symbol, scales = "free_x", ncol = 1) +
      guides(color = "none") +
      scale_color_manual(values = pals::okabe(8)) +
      scale_y_discrete(position = "right") +
      scale_x_continuous(breaks = scales::pretty_breaks(3)) +
      labs(
        y = NULL, x = bquote("Log"[2]~"CPM"),
        title = glue("All cells: {my_contrast}")
      ) +
      theme(strip.text = element_text(face = "italic"), panel.spacing = unit(0.5, "lines"))
    my_ggsave(
      slug = glue("dots-donor-{my_contrast}-3"),
      out_dir = glue("{out_dir}/dots-donor"),
      type = "pdf",
      plot = p,
      scale = 1, units = "in", dpi = 300,
      # width = 1 + length(my_ids) * 0.8,
      # height = 7
      width = 4,
      height = length(my_ids) * 2,
      limitsize = FALSE
    )
  }

  # drug_auc <- presto::wilcoxauc(a1$logcpm, a1$obs$case_drug) %>% arrange(pval) %>%
  #   mutate(gene = ensembl_to_symbol[feature]) %>% as_tibble
  # my_ids <- (
  #   drug_auc %>% filter(group == "CasePD1CTLA4") %>% arrange(-auc) %>% select(-group, -statistic) %>% filter(pct_out < 90) %>% head(10)
  # )$feature

  # Summarize percent of cells expressing each gene in each donor, instead of logCPM
  stopifnot(nrow(a1$obs) == ncol(a1$counts))
  y <- with(a1$obs, model.matrix(~ 0 + donor))
  y <- as(y, "dgCMatrix")
  pb <- as((a1$counts > 0) %*% y, "dgCMatrix")
  colnames(pb) <- str_replace(colnames(pb), "donor", "")
  n_cells_per_donor <- a1$obs %>% count(donor)
  stopifnot(n_cells_per_donor$donor == colnames(pb))
  pb <- sweep(pb, 2, as.numeric(n_cells_per_donor$n), "/")
  pb <- as.matrix(log1p(pb))
  #
  pb_meta <- tibble(donor = colnames(pb))
  #
  pb_meta <- left_join(
    pb_meta,
    a1$obs %>% group_by(case, case_drug, donor, drug, sex) %>% count() %>% ungroup,
    by = c("donor")
  )
  pb <- pb[,as.character(pb_meta$donor)]
  stopifnot(all(pb_meta$donor == colnames(pb)))
  #
  des1 <- with(pb_meta, model.matrix(
    ~ 0 + case_drug + sex
  ))
  colnames(des1) <- str_replace_all(colnames(des1), "[^[:alnum:] ]", "")
  colnames(des1) <- str_replace(colnames(des1), "casedrug", "")
  cont <- makeContrasts(
    contrasts = list(
      "ControlPD1 - ControlNone",
      "CasePD1 - ControlPD1",
      "CasePD1CTLA4 - CasePD1"
    ),
    levels = des1
  )
  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
  stopifnot(nrow(des1) == nrow(pb_meta))
  #
  object <- as.matrix(pb)
  object <- object[rowMeans(object) > 0.5,]
  fit1 <- lmFit(object = object, design = des1)
  fit1 <- eBayes(fit1)
  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
  #
  fit2 <- contrasts.fit(fit1, cont)
  fit2 <- eBayes(fit2)
  #
  res <- rbindlist(lapply(colnames(cont), function(this_contrast) {
    res <- topTable(fit2, coef = this_contrast, number = 1e6, confint = TRUE)
    res$contrast <- this_contrast
    colnames(res)[1] <- "Gene"
    res$ensembl_id <- rownames(res)
    return(res)
  }))
  #
  fwrite(res, glue("{out_dir}/de_donor_percent.tsv.gz"), sep = "\t")
  write_de_xlsx(
    res %>%
      mutate_if(is.numeric, signif, 4),
    glue("{out_dir}/de_donor_percent.xlsx"),
    col = "contrast"
  )
  # res %>% filter(contrast == "CasePD1CTLA4vsCasePD1", logFC > 0) %>% head


  # Drug DE per cluster
  stopifnot(nrow(a1$obs) == ncol(a1$counts))
  y <- with(a1$obs, model.matrix(~ 0 + factor(leiden):donor))
  y <- as(y, "dgCMatrix")
  # y <- sweep(y, 2, colSums(y), "/")
  pb <- as(a1$counts %*% y, "dgCMatrix")
  pb <- pb[,colSums(pb) > 0]
  pb <- do_log1p_cpm(pb, median(colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "donor", "")
  colnames(pb) <- str_replace(colnames(pb), "factor\\(leiden\\)", "")
  colnames(pb) <- str_replace(colnames(pb), "source", "")
  #
  pb_meta <- tibble(leiden_donor = colnames(pb)) %>%
    mutate(leiden = str_split_fixed(leiden_donor, ":", 2)[,1]) %>%
    mutate(donor = str_split_fixed(leiden_donor, ":", 2)[,2])
  #
  pb_meta <- left_join(
    pb_meta,
    a1$obs %>% group_by(leiden, case, case_drug, donor, drug, sex) %>% count() %>% ungroup,
    by = c("leiden", "donor")
  )
  pb <- pb[,as.character(pb_meta$leiden_donor)]
  stopifnot(all(pb_meta$leiden_donor == colnames(pb)))

  cont <- makeContrasts(
    contrasts = list(
      "ControlPD1 - ControlNone",
      "CasePD1 - ControlPD1",
      "CasePD1CTLA4 - CasePD1"
    ),
    levels = des1
  )
  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
  res <- rbindlist(lapply(unique(pb_meta$leiden), function(this_cluster) {
    ix_cluster <- pb_meta$leiden == this_cluster
    des1 <- with(pb_meta[ix_cluster,], model.matrix(
      ~ 0 + case_drug + sex
    ))
    colnames(des1) <- str_replace_all(colnames(des1), "[^[:alnum:] ]", "")
    colnames(des1) <- str_replace(colnames(des1), "casedrug", "")
    stopifnot(nrow(des1) == nrow(pb_meta[ix_cluster,]))
    #
    object <- as.matrix(pb[,ix_cluster])
    object <- object[rowMeans(object) > 0.5,]
    fit1 <- lmFit(object = object, design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
    #
    fit2 <- contrasts.fit(fit1, cont)
    fit2 <- eBayes(fit2)
    #
    res <- rbindlist(lapply(colnames(cont), function(this_contrast) {
      res <- topTable(fit2, coef = this_contrast, number = 1e6, confint = TRUE)
      res$contrast <- this_contrast
      colnames(res)[1] <- "Gene"
      res$ensembl_id <- rownames(res)
      return(res)
    }))
    res$cluster <- this_cluster
    return(res)
  }))
  fwrite(res, glue("{out_dir}/de_cluster.tsv.gz"), sep = "\t")

  for (this_contrast in colnames(cont)) {
    write_de_xlsx(
      res %>% filter(contrast == this_contrast) %>%
        mutate_if(is.numeric, signif, 4),
      glue("{out_dir}/de_cluster_{this_contrast}.xlsx"),
      col = "cluster"
    )
  }

  for (this_contrast in colnames(cont)) {
    for (this_cluster in unique(res$cluster)) {
      p <- plot_limma_volcano(res %>% filter(cluster == this_cluster, contrast == this_contrast)) +
        labs(title = glue("Cluster {this_cluster}: {this_contrast}"))
      my_ggsave(
        glue("volcano-cluster-{this_cluster}-{this_contrast}"),
        out_dir = glue("{out_dir}/volcano-cluster"),
        plot = p,
        type = "pdf",
        scale = 1, width = 5, height = 4, units = "in", dpi = 300
      )
    }
  }

  for (this_contrast in colnames(cont)) {
    res_summary <- res %>% filter(contrast == this_contrast) %>%
      group_by(cluster) %>%
      summarize(
        n_up   = sum(logFC >  log2(1.5) & adj.P.Val < 0.05),
        n_down = sum(logFC < -log2(1.5) & adj.P.Val < 0.05)
      ) %>%
      arrange(as.numeric(cluster))
    fwrite(res_summary, glue("{out_dir}/de_cluster_summary_{this_contrast}.tsv"), sep = "\t")
    #
    p <- ggplot(
      res_summary %>% pivot_longer(cols = c("n_up", "n_down")) %>%
        mutate(value = ifelse(name == "n_down", -value, value))
    ) +
      aes(x = value, y = cluster, fill = name) +
      annotate(
        geom = "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = seq(from = 1, to = max(as.numeric(res_summary$cluster)), by = 2) - 0.5,
        ymax = seq(from = 1, to = max(as.numeric(res_summary$cluster)), by = 2) + 0.5,
        alpha = 0.2
      ) +
      geom_colh() +
      scale_y_discrete(limits = rev(levels(naturalfactor(res_summary$cluster)))) +
      scale_fill_manual(
        # values = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[c(9,3)],
        values = okabe(8)[c(6,7)],
        guide = "none"
      ) +
      # scale_x_continuous(labels = abs) +
      scale_x_continuous(labels = abs, breaks = scales::pretty_breaks(4)) +
      labs(x = NULL, y = NULL, title = this_contrast) +
      theme(plot.title = element_text(size = 8))
    my_ggsave(
      glue("bars-summary-{this_contrast}"),
      out_dir = out_dir,
      type = "pdf",
      plot = p,
      scale = 1,
      width = 2,
      height = nrow(res_summary) * 0.2 + 1,
      units = "in",
      dpi = 300
    )
  }

  # my_ids <- (
  #   res %>% filter(cluster == 1, contrast == "CasePD1CTLA4vsCasePD1") %>%
  #     top_n(n = 5, wt = -log10(P.Value) * abs(logFC))
  # )$ensembl_id
  # for (my_id in my_ids) {
  #   my_symbol <- ensembl_to_symbol[my_id]
  #   ix_cells <- a1$obs$case_drug %in% c("CasePD1", "CasePD1CTLA4")
  #   p <- plot_hexgene(
  #     x = a1$obs$UMAP1[ix_cells],
  #     y = a1$obs$UMAP2[ix_cells],
  #     z = as.numeric(a1$logcpm[my_id,ix_cells]),
  #     group = a1$obs$case_drug[ix_cells]
  #   ) +
  #   facet_wrap(~ group) +
  #   labs(title = my_symbol) +
  #   theme(
  #     panel.spacing = unit(0.5, "lines"),
  #     plot.title = element_text(face = "italic")
  #   )
  #   my_ggsave(
  #     glue("umap-{safe(my_symbol)}"),
  #     out_dir = glue("{out_dir}/umap/"),
  #     type = "pdf",
  #     plot = p,
  #     scale = 0.8,
  #     width = 7,
  #     height = 5,
  #     units = "in",
  #     dpi = 300
  #   )
  # }
  # my_ids <- (
  #   res %>% filter(cluster == 2, contrast == "ControlPD1vsControlNone") %>%
  #     top_n(n = 5, wt = -log10(P.Value) * abs(logFC))
  # )$ensembl_id
  # for (my_id in my_ids) {
  #   my_symbol <- ensembl_to_symbol[my_id]
  #   ix_cells <- a1$obs$case_drug %in% c("ControlPD1", "ControlNone")
  #   p <- plot_hexgene(
  #     x = a1$obs$UMAP1[ix_cells],
  #     y = a1$obs$UMAP2[ix_cells],
  #     z = as.numeric(a1$logcpm[my_id,ix_cells]),
  #     group = a1$obs$case_drug[ix_cells]
  #   ) +
  #   facet_wrap(~ group) +
  #   labs(title = my_symbol) +
  #   theme(
  #     panel.spacing = unit(0.5, "lines"),
  #     plot.title = element_text(face = "italic")
  #   )
  #   my_ggsave(
  #     glue("umap-{safe(my_symbol)}"),
  #     out_dir = glue("{out_dir}/umap/"),
  #     type = "pdf",
  #     plot = p,
  #     scale = 0.8,
  #     width = 7,
  #     height = 5,
  #     units = "in",
  #     dpi = 300
  #   )
  # }

  pb_meta$case_drug3 <- as.character(pb_meta$case_drug)
  pb_meta$case_drug3[pb_meta$case_drug3 == "ControlNone"] <- "Control"
  pb_meta$case_drug3[pb_meta$case_drug3 == "ControlPD1"] <- "Control"
  pb_meta$case_drug3 <- factor(pb_meta$case_drug3, c("Control", "CasePD1", "CasePD1CTLA4"))
  #
  for (my_contrast in unique(res$contrast)) {
    for (my_cluster in unique(pb_meta$leiden)) {
      # my_cluster <- 7
      # my_contrast <- "CasePD1CTLA4vsCasePD1"
      my_res <- res %>%
        filter(cluster == my_cluster, contrast == my_contrast) %>%
        top_n(n = 20, wt = -log10(P.Value) * abs(logFC))
      my_ids <- my_res$ensembl_id
      ix_cluster <- pb_meta$leiden == my_cluster
      x <- rbindlist(lapply(my_ids, function(this_ens) {
        cbind(pb_meta[ix_cluster,],
          ens = this_ens,
          gene = pb[this_ens,ix_cluster],
          symbol = ensembl_to_symbol[this_ens]
        )
      }))
      x$symbol <- factor(as.character(x$symbol), my_res$Gene)
      p <- ggplot() +
        geom_quasirandom(
          data = x,
          mapping = aes(x = gene, y = case_drug, group = case_drug, color = case_drug),
          width = 0.4,
          size = 3,
          groupOnX = FALSE
        ) +
        geom_text(
          data = my_res %>% rename(symbol = Gene) %>% mutate(symbol = factor(symbol)),
          mapping = aes(label = glue("{signif(2^logFC, 2)}-fold p={format.pval(P.Value, 1)}")),
          x = -Inf, y = Inf, hjust = 0, vjust = 1
        ) +
        # facet_wrap(~ symbol, scales = "free_x", ncol = length(my_ids)) +
        facet_wrap(facets = vars(symbol), scales = "free_x", ncol = 1) +
        guides(color = "none") +
        scale_color_manual(values = pals::okabe(8)) +
        scale_y_discrete(position = "right") +
        scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        labs(
          y = NULL, x = bquote("Log"[2]~"CPM"),
          title = glue("Cluster {my_cluster} - {my_contrast}")
        ) +
        theme(strip.text = element_text(face = "italic"), panel.spacing = unit(0.5, "lines"))
      my_ggsave(
        slug = glue("dots-cluster-{my_cluster}-{my_contrast}"),
        out_dir = glue("{out_dir}/dots-cluster"),
        type = "pdf",
        plot = p,
        scale = 1, units = "in", dpi = 300,
        # width = length(my_ids) * 2.1,
        # height = 3 
        width = 4,
        height = length(my_ids) * 2
      )
      p <- ggplot() +
        geom_quasirandom(
          data = x,
          mapping = aes(x = gene, y = case_drug3, group = case_drug3, color = case_drug3),
          width = 0.4,
          size = 3,
          groupOnX = FALSE
        ) +
        # facet_wrap(~ symbol, scales = "free_x", ncol = length(my_ids)) +
        geom_text(
          data = my_res %>% rename(symbol = Gene) %>% mutate(symbol = factor(symbol)),
          mapping = aes(label = glue("{signif(2^logFC, 2)}-fold p={format.pval(P.Value, 1)}")),
          x = -Inf, y = Inf, hjust = 0, vjust = 1
        ) +
        # facet_wrap(~ symbol, scales = "free_x", ncol = length(my_ids)) +
        facet_wrap(facets = vars(symbol), scales = "free_x", ncol = 1) +
        guides(color = "none") +
        scale_color_manual(values = pals::okabe(8)) +
        scale_y_discrete(position = "right") +
        scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        labs(
          y = NULL, x = bquote("Log"[2]~"CPM"),
          title = glue("Cluster {my_cluster} - {my_contrast}")
        ) +
        theme(strip.text = element_text(face = "italic"), panel.spacing = unit(0.5, "lines"))
      my_ggsave(
        slug = glue("dots-cluster-{my_cluster}-{my_contrast}-3"),
        out_dir = glue("{out_dir}/dots-cluster"),
        type = "pdf",
        plot = p,
        scale = 1, units = "in", dpi = 300,
        # width = length(my_ids) * 2.1,
        # height = 3 
        width = 4,
        height = length(my_ids) * 2
      )
    }
  }

  if (analysis_name == "a12_4_4_t4_cd4_2_2") {
    for (my_symbol in c("GNLY", "GZMA")) {
      my_id <- names(which(ensembl_to_symbol == my_symbol))
      p <- plot_hexgene(
        x = a1$obs$UMAP1,
        y = a1$obs$UMAP2,
        z = as.numeric(a1$logcpm[my_id,]),
        palette = "oslo",
        group = append_n(a1$obs$case_drug),
        bins = 51
      ) +
      facet_wrap(~ group) +
      labs(title = my_symbol) +
      theme(
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "italic")
      )
      my_ggsave(
        glue("umap-case_drug-{safe(my_symbol)}"),
        out_dir = glue("{out_dir}/umap-facet"),
        type = "pdf",
        plot = p,
        scale = 0.8,
        width = 7,
        height = 8,
        units = "in",
        dpi = 300
      )
    }
  }

if (FALSE) {
  # my_ids <- (
  #   res %>% filter(cluster == 2, contrast == "ControlPD1vsControlNone") %>%
  #     top_n(n = 5, wt = -log10(P.Value) * abs(logFC))
  # )$ensembl_id
  # my_ids <- (
  #   res %>% filter(cluster == 6, contrast == "CasePD1CTLA4vsCasePD1") %>%
  #     top_n(n = 5, wt = -log10(P.Value) * abs(logFC))
  # )$ensembl_id
  # ensembl_to_symbol[my_ids]
  for (my_contrast in unique(res$contrast)) {
    for (my_cluster in unique(pb_meta$leiden)) {
      # my_cluster <- 6
      # my_contrast <- "CasePD1CTLA4vsCasePD1"
      my_ids <- (
        res %>% filter(cluster == my_cluster, contrast == my_contrast) %>%
          top_n(n = 15, wt = -log10(P.Value) * abs(logFC))
      )$ensembl_id
      for (my_id in my_ids) {
        my_symbol <- ensembl_to_symbol[my_id]
        # ix_cells <- a1$obs$case_drug %in% c("CasePD1CTLA4", "CasePD1")
        ix_cells <- a1$obs$case_drug %in% as.character(str_split_fixed(my_contrast, "vs", 2))
        p <- plot_hexgene(
          x = a1$obs$UMAP1[ix_cells],
          y = a1$obs$UMAP2[ix_cells],
          z = as.numeric(a1$logcpm[my_id,ix_cells]),
          palette = "oslo",
          # palette = "batlow", # davos is better
          group = append_n(a1$obs$case_drug[ix_cells])
        ) +
        facet_wrap(~ group) +
        labs(title = my_symbol) +
        theme(
          panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(face = "italic")
        )
        my_ggsave(
          glue("umap-cluster-{my_cluster}-{safe(my_symbol)}"),
          out_dir = glue("{out_dir}/umap-cluster/{my_contrast}"),
          type = "pdf",
          plot = p,
          scale = 0.8,
          width = 7,
          height = 5,
          units = "in",
          dpi = 300
        )
        my_ggsave(
          glue("umap-cluster-{my_cluster}-{safe(my_symbol)}"),
          out_dir = glue("{out_dir}/umap-cluster/{my_contrast}/slim"),
          type = "pdf",
          plot = p +
          theme(
            legend.position = "none",
            strip.text = element_blank(),
            panel.spacing = unit(0.5, "lines"),
            plot.title = element_text(face = "italic")
          ),
          scale = 0.8,
          width = 6.5,
          height = 3,
          units = "in",
          dpi = 300
        )
      }
    }
  }
}

  # res %>% group_by(contrast) %>% summarize(n = sum(logFC > log2(2) & adj.P.Val < 0.05))

}


# B cells
# Ratio of IgG plasma cells (cluster 8) to IgA plasma cells (10, 1, 6, 7, 4, 9)
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
  if (analysis_name == "a12_4_4_b5_1_3") {
    a1$obs$leiden <- a1$obs$leiden0.933
  }
  out_dir <- glue("results/a20/{analysis_name}/figures")
  source("R/functions/composition.R")
  source("R/plot-composition.R")
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
  a1$obs$case <- factor(my_obs$case, c("Control", "Case"))

  x <- a1$obs %>%
    group_by(donor) %>%
    summarize(ratio = sum(leiden == 8) / sum(leiden %in% c(1, 4, 6, 7, 9, 10))) %>%
    left_join(a1$obs %>% select(donor, case) %>% unique, by = "donor")
  x_pval <- wilcox.test(ratio ~ case, x)$p.value
  p <- ggplot(x) +
    aes(y = case, x = ratio, fill = case) +
    geom_quasirandom(groupOnX = FALSE, shape = 21, size = 3, stroke = 0.2) +
    scale_fill_manual(values = pals::okabe()[1:2], guide = FALSE) +
    labs(
      x = "Ratio of IgG to IgA cells", y = NULL,
      subtitle = glue("Wilcoxon Rank Sum Test P = {signif(x_pval, 2)}")
    )
  my_ggsave(
    "igg-iga-ratio",
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 1,
    width = 5,
    height = 2,
    units = "in", dpi = 300
  )


# CD8 - Selected genes {{{
########################################################################

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
#
# my_leiden <- "leiden1.22"
# a1$obs$cluster <- a1$obs[[my_leiden]]
# a1$obs$cluster <- recluster_cd8_leiden122(a1$obs$cluster)
#
my_leiden <- "leiden1.51"
a1$obs$cluster <- a1$obs[[my_leiden]]
a1$obs$cluster <- recluster_cd8_leiden151(a1$obs$cluster)
#
out_dir <- as.character(glue("results/a20/{analysis_name}/figures/umap"))
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
a1$log2cpm <- do_log2cpm(a1$counts, median(colSums(a1$counts)))

selected_genes <- c("ITGB2", "ITGBA1", "ITGAE", "KLRG1", "CXCR3")
selected_ids <- names(ensembl_to_symbol[ensembl_to_symbol %in% selected_genes])

for (my_id in selected_ids) {
  p <- plot_hexgene(
    x = a1$obs$UMAP1,
    y = a1$obs$UMAP2,
    z = as.numeric(a1$log2cpm[my_id,]),
    bins = 47,
    palette = "oslo"
  ) +
  labs(title = ensembl_to_symbol[my_id]) +
  theme(
    # legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(face = "italic")
  )
  my_ggsave(
    glue("umap-{safe(ensembl_to_symbol[my_id])}"),
    out_dir = out_dir,
    type = "pdf",
    plot = p,
    scale = 0.8,
    width = 3.5,
    # height = 3,
    height = 5,
    units = "in",
    dpi = 300
  )
}

de_file <- glue("results/a20/{analysis_name}/figures/pseudobulk_de_ova.tsv.gz")
stopifnot(file.exists(de_file))
de <- fread(de_file)
de$cluster <- str_split_fixed(de$contrast, " ", 2)[,1]
de$cluster <- factor(de$cluster, rev(naturalsort(unique(de$cluster))))

plot_de_by_gene <- function(d) {
  ggplot(d) +
    aes(x = logFC, y = cluster) +
    ggforestplot::geom_stripes() +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_point() +
    geom_errorbarh(
      mapping = aes(xmin = CI.L, xmax = CI.R),
      height = 0
    ) +
    scale_x_continuous(
      breaks = seq(-10, 10, by = 1),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    # scale_y_continuous(
    #   expand = c(0, 0),
    #   breaks = seq(1, max(d$cluster))
    # ) +
    # annotation_logticks(sides = "b", size = 0.3) +
    facet_grid(~ Gene) +
    labs(y = NULL) +
    theme(
      strip.text = element_text(face = "italic"),
      panel.spacing = unit(0.5, "lines")
    )
}
p <- plot_de_by_gene(
  de %>% filter(ensembl_id %in% selected_ids)
)
my_ggsave(
  "selected-genes-1",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1,
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

# Pseudobulk at the cluster level
########################################################################
y <- with(a1$obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
y <- as(y, "dgCMatrix")
# y <- sweep(y, 2, colSums(y), "/") # means
pb <- as(a1$counts %*% y, "dgCMatrix")
pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
#
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
  a1$obs %>%
    select(
      donor, case
    ) %>%
    group_by(donor, case) %>%
    summarize_if(is.numeric, mean),
  by = "donor"
)
pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
stopifnot(nrow(pb_meta) == ncol(pb))

d <- rbindlist(lapply(selected_ids, function(selected_id) {
  retval <- pb_meta
  retval[['log2cpm']] <- as.numeric(pb[selected_id,])
  retval[['ensembl_id']] <- selected_id
  as.data.table(retval)
}))
d[['symbol']] <- ensembl_to_symbol[d$ensembl_id]
# d$cluster <- factor(d$cluster, rev(naturalsort(unique(d$cluster))))
d$cluster <- factor(d$cluster, rev(c(4,3,7,1,6,11,10,12,8,5,9,2)))

cluster_colors <- mpn65[seq_along(unique(d$cluster))]
names(cluster_colors) <- unique(d$cluster)
plot_by_gene <- function(d) {
  ggplot(d) +
    aes(x = log2cpm, y = cluster, fill = cluster) +
    ggforestplot::geom_stripes() +
    # geom_boxplot(size = 0.3) +
    # geom_quasirandom(groupOnX = FALSE, shape = 21, size = 2, stroke = 0.3) +
    stat_summary(
      fun.data = function(x) {
        retval <- quantile(x, probs = c(0.25, 0.5, 0.75))
        names(retval) <- c("ymin", "y", "ymax")
        retval
      },
      shape = 21, size = 0.5, stroke = 0.3
    ) +
    scale_fill_manual(values = cluster_colors) +
    scale_x_continuous(breaks = pretty_breaks(3)) +
    facet_grid(~ symbol) +
    labs(y = NULL) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "italic"),
      panel.spacing = unit(1, "lines")
    )
}
p <- plot_by_gene(d)
my_ggsave(
  "selected-genes-2",
  out_dir = out_dir,
  type = "pdf",
  plot = p,
  scale = 1,
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

# }}}

# control pd1 vs control none bars {{{

de <- fread("paper/de-contrasts.tsv.gz")

my_de <- de %>%
  filter(!str_detect(cluster, "-all$")) %>%
  filter(contrast == "ControlPD1 vs ControlNone") %>%
  group_by(analysis, cluster) %>%
  summarize(
    n_up = sum( log_fc > log2(1.5) & adj_p_val < 0.1 ),
    n_down = sum(-log_fc > log2(1.5) & adj_p_val < 0.1 )
  ) %>%
  ungroup()

out_dir <- "results/a20/de-controlpd1-vs-controlnone"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (this_analysis in unique(my_de$analysis)) {
  this_de <- my_de %>%
    filter(analysis == this_analysis) %>%
    arrange(naturalfactor(cluster))
  fwrite(this_de, glue("{out_dir}/{this_analysis}.tsv"), sep = "\t")
  my_breaks <- c(-max(this_de$n_down), max(this_de$n_up))
  if (all(my_breaks == c(0, 0))) {
    my_breaks <- c(-1, 1)
  }
  #
  p <- ggplot(
    this_de %>% pivot_longer(cols = c("n_up", "n_down")) %>%
      mutate(value = ifelse(name == "n_down", -value, value))
  ) +
    aes(x = value, y = cluster, fill = name) +
    # annotate(
    #   geom = "rect",
    #   xmin = -Inf,
    #   xmax = Inf,
    #   ymin = seq(from = 1, to = max(as.numeric(this_de$cluster)), by = 2) - 0.5,
    #   ymax = seq(from = 1, to = max(as.numeric(this_de$cluster)), by = 2) + 0.5,
    #   alpha = 0.2
    # ) +
    geom_stripes(even = "#ffffff", odd = "#eeeeee") +
    geom_vline(linewidth = 0.3, xintercept = 0) +
    geom_colh() +
    scale_y_discrete(
      labels = \(d) str_remove(d, "^.+-"),
      limits = rev(levels(naturalfactor(this_de$cluster)))
    ) +
    scale_fill_manual(
      # values = RColorBrewer::brewer.pal(name = "RdBu", n = 11)[c(9,3)],
      values = okabe(8)[c(3,4)],
      guide = "none"
    ) +
    # scale_x_continuous(labels = abs) +
    scale_x_continuous(labels = abs, breaks = pretty_breaks(4)(my_breaks)) +
    labs(x = NULL, y = NULL, title = this_analysis) +
    theme(plot.title = element_text(size = 8))
  my_ggsave(
    glue("{this_analysis}"),
    out_dir = file.path(out_dir, "bars"),
    type = "pdf",
    plot = p,
    scale = 1,
    width = 3,
    height = nrow(this_de) * 0.2 + 1,
    units = "in",
    dpi = 300
  )
}

# }}}

