#' @importFrom ggplot2 theme_classic theme element_rect element_line
#'   element_blank element_text
#' @importFrom grid unit
theme_kamil <- theme_classic(base_family = "Helvetica") +
theme(
  axis.line             = element_blank(),
  axis.text             = element_text(size = 16),
  axis.ticks            = element_line(size = 0.4),
  axis.title            = element_text(size = 16),
  legend.background     = element_rect(colour = "transparent", fill = "transparent"),
  legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
  legend.text           = element_text(size = 16),
  legend.title          = element_text(size = 16),
  panel.background      = element_rect(fill = "transparent"),
  panel.border          = element_rect(size = 0.5, fill = NA),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  panel.spacing         = unit(2, "lines"),
  plot.background       = element_rect(fill = "transparent", color = NA),
  plot.subtitle         = element_text(size = 14),
  plot.title            = element_text(size = 16),
  strip.background      = element_blank(),
  strip.text            = element_text(size = 16)
)

scientific_10 <- function(x) {
  ifelse(x < 0.01,
    gsub("e", "%*%10^", scales::scientific_format(digits = 1)(x)),
    signif(x, 1)
  )
}
