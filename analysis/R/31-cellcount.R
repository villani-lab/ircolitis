
library(glue)
library(ggplot2)
library(ggbeeswarm)
library(rhdf5)
library(dplyr)
library(data.table)
library(pbapply)
library(rhdf5)
library(sitools)
library(ggforestplot)
library(naturalsort)
library(forcats)
library(scales)

source("R/functions/theme-kamil.R")
theme_set(theme_kamil)


out_dir <- "results/estimated-cells"

my_celltypes <- c("b", "cd4", "cd8", "myeloid")
my_celltype <- my_celltypes[1]

my_res <- rbindlist(pblapply(my_celltypes, function(my_celltype) {
  h5_file <- glue("paper/tissue-{my_celltype}.h5")
  my_obs <- h5read(h5_file, "obs")
  my_res <- my_obs %>%
    count(case, donor, leiden) %>%
    # group_by(donor) %>%
    # mutate(freq = n / sum(n)) %>%
    # ungroup() %>%
    mutate(celltype = my_celltype) %>%
    mutate(cluster = glue("{celltype}-{leiden}"))
  return(my_res)
}))
my_res %<>% group_by(donor) %>% mutate(freq = n / sum(n)) %>% ungroup

# Avg CD45+ Cells/Biopsy (x10^5)
d <- data.frame(
  avg_cd45_cells = c(
    0.1382, 0.5043, 0.1763, 0.3758, 0.2029, 0.3590, 0.4399, 0.2621, 0.8103,
    3.5280, 6.5957, 2.9240, 8.8366, 3.6426
  ) * 1e5,
  donor = c(
    "SIC_19", "SIC_33", "SIC_35", "SIC_53", "MC_2", "SIC_188", "SIC_196", "MC_9",
    "SIC_186", "SIC_140", "SIC_36", "SIC_40", "SIC_43", "SIC_71"
  )
)
d <- my_res %>% inner_join(d, by = "donor") %>%
  mutate(estimated_cells = freq * avg_cd45_cells)
d

d %>% group_by(case, donor) %>%
  summarize(
    n_cd8 = sum(n[celltype == "cd8"]),
    n_cd4 = sum(n[celltype == "cd4"])
  ) %>%
  group_by(case, donor) %>%
  summarize(
    n_tcell = n_cd8 + n_cd4
  )

for (my_celltype in unique(d$celltype)) {
  my_d <- d %>% filter(celltype == my_celltype) %>%
    mutate(leiden = fct_rev(naturalfactor(leiden)))
  my_breaks <- scales::log_breaks(n = 5)(my_d$estimated_cells)
  my_breaks <- my_breaks[2:(length(my_breaks)-1)]
  p <- ggplot(my_d) +
    aes(x = estimated_cells, y = leiden, fill = case) +
    ggforestplot::geom_stripes(even = "#ffffff", odd = "#eeeeee") +
    geom_vline(size = 0.5, color = "white", xintercept = my_breaks) +
    annotation_logticks(sides = "b", size = 0.3) +
    geom_boxplot(linewidth = 0.3, outlier.shape = NA, alpha = 0.5) +
    geom_point(
      position = position_quasirandom(dodge.width = 1), size = 1.5, stroke = 0.3, shape = 21
    ) +
    scale_x_log10(labels = f2si, breaks = my_breaks) +
    scale_fill_manual(
      name = NULL,
      values = c(Case = pals::okabe(3)[2], Control = "gray50"),
      guide = guide_legend(reverse = TRUE)
    ) +
    scale_y_discrete(position = "right") +
    labs(
      x = "Estimated number of cells per biopsy",
      y = NULL,
      title = glue("Estimated number of {my_celltype} cells")
    ) +
    theme(
      legend.position = "top"
    )
  p_file <- glue("{out_dir}/estimated-cells-{my_celltype}.pdf")
  dir.create(dirname(p_file), recursive = TRUE, showWarnings = FALSE)
  message(p_file)
  ggsave(p_file, p, width = 5, height = 5)
}

