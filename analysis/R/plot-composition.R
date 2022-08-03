
#' Given a dataframe with one row for each cell, summarize the composition of
#' each donor's cells across multiple clusters.
#
#' @param d A dataframe with one row for each cell, and multiple columns.
#' @param fill The name of the column that describes the type of each donor
#'   (e.g., case/control).
#' @param group The name of the donor column.
#' @param x The name of the cluster column.
#'
#' @examples
#'      get_composition(obs, fill = case, group = donor, x = leiden)
#'
get_composition <- function(d, fill, group, x) {
  d %>%
    dplyr::select({{ fill }}, {{ group }}, {{ x }}) %>%
    dplyr::group_by({{ fill }}, {{ group }}, {{ x }}) %>%
    dplyr::summarize(n = n(), .groups = "drop_last") %>%
    dplyr::mutate(freq = 100 * n / sum(n))
}

# get_composition <- function(d, groups=NULL) {
#   d %>%
#     group_by(across({{groups}})) %>%
#     summarize(n = n()) %>%  # Could also use tally() here
#     ungroup %>%
#     mutate(percent = 100 * n / sum(n))
# }


plot_composition <- function(d, fill, group, x, legend.position = "bottom") {
  n_x <- d %>% dplyr::select({{ x }}) %>% unique %>% nrow
  d <- d %>%
    dplyr::select({{ fill }}, {{ group }}, {{ x }}) %>%
    dplyr::group_by({{ fill }}, {{ group }}, {{ x }}) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(freq = 100 * n / sum(n))
  # d_t <- rbindlist(lapply(unique(d$leiden), function(this_leiden) {
  #   x <- d %>% dplyr::filter(leiden == this_leiden)
  #   if (length(table(x$case)) != 2) {
  #     return(
  #       tibble(
  #         statistic = NA, p.value = 1, method = NA,
  #         alternative = NA, leiden = this_leiden
  #       )
  #     )
  #   }
  #   # x <- broom::tidy(t.test(n ~ case, x))
  #   # x <- broom::tidy(t.test(freq ~ case, x))
  #   x <- broom::tidy(wilcox.test(log(freq) ~ case, x))
  #   x$leiden <- this_leiden
  #   return(x)
  # }))
  d_label <- d %>%
    dplyr::group_by({{ fill }}) %>%
    # summarise(n = sum(n)) %>%
    dplyr::summarise(n = length(unique({{ group }}))) %>%
    dplyr::mutate(label = sprintf("%s (n = %s)", {{ fill }}, n))
  # d_t <- inner_join(
  #   d_t,
  #   d %>% dplyr::group_by(leiden) %>% dplyr::mutate(y = mean(freq) + 4) %>% select(leiden, y) %>% unique,
  #   by = "leiden"
  # )
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
        y = c(10^(0:2), 0.5 * 10^(0:2))
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
      data = d %>% dplyr::group_by({{ x }}, {{ fill }}) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        x = factor({{ x }}), ymin = ymin, y = y, ymax = ymax,
        group = {{ fill }}, fill = {{ fill }}
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9),
      size = 0,
      alpha = 0
    ) +
    annotate(
      geom = "rect",
      ymin = 0, ymax = Inf,
      xmin = as.integer(seq(n_x)) - 0.5,
      xmax = as.integer(seq(n_x)) + 0.5,
      fill = rep(
        c(NA, "black"),
        length.out = n_x
      ),
      alpha = 0.1
    ) +
    geom_crossbar(
      data = d %>% dplyr::group_by({{ x }}, {{ fill }}) %>%
        dplyr::summarize(
          ymin = quantile(freq, 0.75),
          y = quantile(freq, 0.5),
          ymax = quantile(freq, 0.25)
        ),
      mapping = aes(
        x = factor({{ x }}), ymin = ymin, y = y, ymax = ymax,
        group = {{ fill }}, fill = {{ fill }} 
      ),
      color = "grey20",
      width = 0.5,
      position = position_dodge(width = 0.9),
      alpha = 0.3,
      size = 0.3
    ) +
    geom_point(
      data = d,
      mapping = aes(x = factor({{ x }}), y = freq, fill = {{ fill }}),
      position = ggbeeswarm::position_quasirandom(groupOnX = TRUE, dodge.width = 0.9),
      alpha = 0.9, shape = 21
    ) +
    # geom_text(
    #   data = d_t,
    #   mapping = aes(
    #     x = factor(leiden),
    #     y = 100,
    #     # y = y + 2,
    #     label = ifelse(
    #       # p.value < 0.05 / 12,
    #       p.value < 0.05 / length(unique(a2$obs$leiden)),
    #       scientific_10(p.value),
    #       ""
    #     )
    #   ),
    #   size = 5, parse = TRUE
    # ) +
    scale_fill_manual(
      name = NULL, values = Okabe_Ito, #pals::okabe(8)[1:8],
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
      # legend.position = "bottom",
      legend.position = legend.position,
      legend.background = element_blank()
      # panel.grid.major.y = element_line(size = 0.1, color = "#00000055"),
      # panel.grid.minor.y = element_line(size = 0.1, color = "#00000055")
    )
  p
}

plot_composition_h <- function(
  d = NULL, fill = NULL, group = NULL, x = NULL, legend.position = "bottom", comp = NULL
) {
  if (!is.null(comp)) {
    d <- comp
  } else {
    # d <- get_composition(d, fill, group, x)
    d <- d %>%
      dplyr::select({{ fill }}, {{ group }}, {{ x }}) %>%
      dplyr::group_by({{ fill }}, {{ group }}, {{ x }}) %>%
      dplyr::summarize(n = n(), .groups = "drop_last") %>%
      dplyr::mutate(freq = 100 * n / sum(n))
  }
  d_label <- d %>%
    dplyr::group_by({{ fill }}) %>%
    # summarise(n = sum(n)) %>%
    dplyr::summarise(n = length(unique({{ group }}))) %>%
    dplyr::mutate(label = sprintf("%s (n = %s)", {{ fill }}, n))
  vlines <- c(0.001, 0.01, 0.1, 1, 10, 100)
  vlines <- vlines[vlines > min(d$freq)]
	ggplot(d) +
		aes(y = factor({{x}}), x = freq, fill = {{fill}}) +
		ggforestplot::geom_stripes(
			odd = "#66666633", even = "#00000000"
		) +
		geom_vline(
			# xintercept = c(0.01, 0.1, 1, 10, 100),
			xintercept = vlines,
			# size = 0.3, alpha = 0.3
			size = 0.3, color = "white"
		) +
		geom_boxploth(
			outlier.shape = NA, coef = 0,
      alpha = 0.3,
      size = 0.3
		) +
		geom_point(
			position = ggbeeswarm::position_quasirandom(
				groupOnX = FALSE, dodge.width = -0.9
			),
			alpha = 0.9, shape = 21, size = 1.5, stroke = 0.3
		) +
		scale_fill_manual(
			name = NULL, values = pals::okabe(),
      labels = d_label$label
		) +
		scale_x_log10(
			name = "Percent",
			labels = function(x) signif(x, 3)
		) +
		annotation_logticks(side = "b", size = 0.3) +
		labs(y = NULL) +
		theme(
			# legend.position = "bottom",
			legend.position = legend.position,
			legend.background = element_blank()
		)
}

























  #fig_width <- 2 + length(unique(a2$obs$leiden)) * 0.7
  #my_ggsave(
  #  "composition-leiden",
  #  #out_dir = glue("figures/{analysis_name}"),
  #  out_dir = out_dir,
  #  type = "pdf",
  #  plot = p +
  #    labs(
  #      title = "Per-donor percent of cells in each cluster (Wilcoxon Rank Sum P-value)",
  #      y = NULL
  #    ),
  #  scale = 1, width = fig_width, height = 4, units = "in", dpi = 300
  #)
  #dd <- a2$obs %>%
  #  dplyr::group_by(donor) %>%
  #  dplyr::summarize(
  #    total_n = n(),
  #    n_1 = sum(leiden == 1)
  #  )
  #d <- left_join(d, dd, by = "donor")
  ## ddd <- fread("results/a12_4_4_min_genes500_n_pcs20/data/cells.tsv.gz")
  ## d <- left_join(
  ##   d,
  ##   ddd %>% group_by(donor) %>% summarize(n_plasma = sum(leiden == 1) / n()),
  ##   by = "donor"
  ## )
  #p <- ggplot(d) +
  #  annotation_logticks() +
  #  aes(x = total_n, y = freq, fill = case) +
  #  geom_point(shape = 21, size = 3) +
  #  facet_wrap(~ leiden, ncol = 5) +
  #  scale_fill_manual(
  #    name = NULL, values = pals::okabe(3)[c(2,3)]
  #  ) +
  #  scale_y_log10() +
  #  scale_x_log10(labels = scales::label_number_si())
  #fig_height <- 2 + ceiling(n_clusters / 5) * 0.7
  #my_ggsave(
  #  "composition-scatter",
  #  #out_dir = glue("figures/{analysis_name}"),
  #  out_dir = out_dir,
  #  type = "pdf",
  #  plot = p +
  #    labs(
  #      title = "Per-donor percent of cells in each cluster",
  #      y = "Percent",
  #      x = "Total cells from each donor"
  #    ),
  #  scale = 2, width = 6, height = fig_height, units = "in", dpi = 300
  #)
  #p <- ggplot(d) +
  #  annotation_logticks() +
  #  aes(x = n_1, y = freq, fill = case) +
  #  geom_point(shape = 21, size = 3) +
  #  facet_wrap(~ leiden, ncol = 5) +
  #  scale_fill_manual(
  #    name = NULL, values = pals::okabe(3)[c(2,3)]
  #  ) +
  #  scale_y_log10() +
  #  scale_x_log10(labels = scales::label_number_si())
  #fig_height <- 2 + ceiling(n_clusters / 5) * 0.7
  #my_ggsave(
  #  "composition-scatter-c1",
  #  #out_dir = glue("figures/{analysis_name}"),
  #  out_dir = out_dir,
  #  type = "pdf",
  #  plot = p +
  #    labs(
  #      title = "Per-donor percent of cells in each cluster",
  #      y = "Percent",
  #      x = "Total cells from cluster 1"
  #    ),
  #  scale = 2, width = 6, height = fig_height, units = "in", dpi = 300
  #)
