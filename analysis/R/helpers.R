
safe <- function(x) {
  janitor::make_clean_names(x, case = "none")
}

print_status <- function(x) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  this_message <- glue("{timestamp}\t{x}")
  message(this_message)
}

cache <- function(filename, x, verbose = TRUE) {
  filename <- as.character(filename)
  if (file.exists(filename)) {
    if (verbose) {
      print_status(glue::glue("Reading cached result from '{filename}'"))
    }
    x <- readRDS(filename)
    if (verbose) {
      print_status(glue::glue("Done reading cached result from '{filename}'"))
    }
    return(x)
  }
  if (grepl("/", filename)) {
    if (verbose) {
      print_status(glue::glue("Creating directories '{dirname(filename)}'"))
    }
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  }
  if (verbose) {
    print_status(glue::glue("Writing result to '{filename}'"))
  }
  saveRDS(x, filename)
  if (verbose) {
    print_status(glue::glue("Done writing result to '{filename}'"))
  }
  return(x)
}

# Jupyter functions -----------------------------------------------------------

#' Print a figure to a PNG file and then display with IRdisplay::display_png()
#' This is necessary to use the showtext package.
# show_plot <- function(obj, width = NULL, height = NULL, units = "in", res = 120) {
#   out_dir <- "notebooks/figures"
#   dir.create(out_dir, showWarnings = FALSE)
# 
#   if (is.null(width)) {
#     width <- getOption("repr.plot.width", 3)
#   }
#   if (is.null(height)) {
#     height <- getOption("repr.plot.height", 3)
#   }
# 
#   filename <- tempfile(pattern = "jupyter-", tmpdir = out_dir)
#   file_png <- sprintf("%s.png", filename)
# 
#   png(file_png, width = width, height = height, units = units, res = res)
#   print(obj)
#   dev.off()
# 
#   IRdisplay::display_png(file = file_png)
# }

show_plot <- function(obj, name = NULL, width = NULL, height = NULL, units = "in", res = 300, optimize_png = TRUE, save_pdf = TRUE, out_dir = "figures") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(width)) {
    width <- getOption("repr.plot.width", 3)
  }
  if (is.null(height)) {
    height <- getOption("repr.plot.height", 3)
  }
  if (is.null(name)) {
    name <- substr(tempfile("", tmpdir = ""), 2, 6)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  # filename <- sprintf("%s/%s_%s", out_dir, name, timestamp)
  filename <- sprintf("%s/%s", out_dir, name)
  file_png <- sprintf("%s.png", filename)
  file_pdf <- sprintf("%s.pdf", filename)

  message(sprintf("\nWriting file %s", file_png))
  ggsave(
    filename = file_png,
    plot = obj, width = width, height = height, units = units, dpi = res
  )
  
  if (optimize_png && Sys.which("pngquant") != "") {
    opt_filename <- sprintf(
      "%s-fs8.png", substr(file_png, 1, nchar(file_png) - 4)
    )
    command <- sprintf(
      'pngquant --ext -fs8.png -- "%s" && mv -f "%s" "%s"',
      file_png, opt_filename, file_png
    )
    system(command)
  }

  if (save_pdf) {
    # pdf.options(encoding = 'ISOLatin2')
    ggsave(
      filename = file_pdf,
      plot = obj,
      device = cairo_pdf,
      width = width, height = height, units = units, dpi = res
    )
  }

  IRdisplay::display_png(file = file_png)
}

my_ggsave <- function(
  slug, out_dir = "figures", types = c("png", "pdf"), optimize = TRUE,
  use_cairo = FALSE, ...
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in types) {
    filename <- file.path(out_dir, glue("{slug}.{ext}"))
    message(glue("Writing {filename}"))
    if (ext == "pdf") {
      my_device <- "pdf"
      if (use_cairo) {
        my_device = cairo_pdf
      }
      ggsave(filename = filename, device = my_device, bg = "transparent", ...)
    } else if (ext == "png" && optimize) {
      ggsave(filename = filename, ...)
      optimize_png(filename)
    } else {
      ggsave(filename = filename, bg = "transparent", ...)
    }
  }
}

f2si <- function(number, unit = "") {
  ix_n <- which(!is.na(number) & number > 0)
  sifactor <- c(
    1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
    0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
    1e+24
  )
  pre <- c(
    " y", " z", " a", " f", " p", " n", " u", " m",
    "", " k", " M", " G", " T", " P", " E", " Z", " Y"
  )
  absolutenumber <- number[ix_n] * sign(number[ix_n])
  ix <- findInterval(absolutenumber, sifactor)
  sistring <- rep("0", length(number))
  if (length(ix) > 0) {
    sistring[ix_n] <- paste(number[ix_n] / sifactor[ix], pre[ix],
      sep = "",
      unit = unit
    )
  }
  else {
    sistring[ix_n] <- as.character(number[ix_n])
  }
  return(sistring)
}

# # FIXME This still does not work... :(
# # Thanks to tomelgin
# # https://stackoverflow.com/questions/11340444/is-there-an-r-function-to-format-number-using-unit-prefix
# f2si <- function(number, rounding = TRUE, digits = ifelse(rounding, NA, 6)) {
#   lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06, 
#     0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21, 
#     1e+24, 1e+27)
#   pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "k", 
#     "M", "G", "T", "P", "E", "Z", "Y", NA)
#   ix <- findInterval(number, lut)
#   if (ix>0 && ix<length(lut)) {# && lut[ix]!=1) {
#     if (rounding==T && !is.numeric(digits)) {
#       sistring <- paste(round(number/lut[ix]), pre[ix])
#     }
#     else if (rounding == T || is.numeric(digits)) {
#       sistring <- paste(signif(number/lut[ix], digits), pre[ix])
#     }
#     else {
#       sistring <- paste(number/lut[ix], pre[ix])
#     }
#   }
#   else {
#     sistring <- as.character(number)
#   }
#   return(sistring)
# }

plot_hexgene <- function(
  x, y, z, group = NULL,
  fun = "mean", bins = 101,
  legend = TRUE, hex = TRUE,
  palette = "davos", direction = -1,
  use_quantile = FALSE,
  italic = TRUE,
  text = TRUE,
  text_size = 7,
  color_title = NULL
) {
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(z))
  if (!is.null(group)) {
    stopifnot(length(x) == length(group))
  }
  z <- as.numeric(z)
  if (is.null(color_title)) {
    if (length(unique(z)) == 2) {
      color_title <- "Proportion"
    } else {
      color_title <- bquote("Mean log"[2]~"CPM")
    }
  }
  d <- data.frame(x, y, z)
  d_text <- data.frame(
    text = glue(
      "{scales::comma(sum(d$z > 0))} ({signif(100 * sum(d$z > 0) / length(d$z), 2)}%) cells"
    )
  )
  if (!is.null(group)) {
    d$group <- group
    d_text <- d %>%
      group_by(group) %>%
      summarize(
        text = glue(
          "{scales::comma(sum(z > 0))} ({signif(100 * sum(z > 0) / length(z), 2)}%) cells"
        ),
        .groups = "keep"
      )
  }
  stat_fun <- stat_summary_hex
  if (!hex) {
    stat_fun <- stat_summary_2d
  }
  plot_colors <- scico::scico(
    n = 20, palette = palette, direction = direction
  )[2:20]
  total_obs <- nrow(d)
  p <- ggplot(d) +
    stat_fun(
      mapping = aes(
        x = x,
        y = y,
        z = z
      ),
      # fun = fun, bins = bins
      bins = bins,
      fun = function(x) {
        mean(x)
      }
    ) +
    theme(
      legend.position = "bottom",
      # legend.position = c(0, 0),
      # legend.justification = c(0, 0),
      legend.background = element_blank(),
      plot.title = element_text(
        face = ifelse(italic, "italic", "plain"),
        size = 20,
        margin = unit(c(0, 0, 0, 0), "lines")
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  if (all(d$z == 0)) {
    p <- p + scale_fill_gradientn(
      colors = plot_colors[1],
      breaks = 0
    )
  } else if (use_quantile) {
    # qn <- quantile(z, c(0.01, 0.99), na.rm = TRUE)
    # qn01 <- rescale(c(qn, range(z)))
    # p <- p + scale_fill_gradientn(
    #   colors = plot_colors,
    #   values = c(0, seq(qn01[1], qn01[2], length.out = 17), 1)
    #   # breaks = scales::pretty_breaks(3)
    # )
    p <- p + scale_fill_gradientn(
      colors = plot_colors,
      breaks = function(x) {
        quantile(x, probs = seq(0, 1, length.out = length(plot_colors)))
      }
    )
  } else {
    p <- p + scale_fill_gradientn(
      colors = plot_colors,
      breaks = scales::pretty_breaks(3)
    )
  }
  if (legend) {
    p <- p + guides(
      fill = guide_colorbar(
        direction = "horizontal", title.position = "top",
        # title = bquote("Mean log"[2]~"CPM"),
        title = color_title,
        barwidth = 8.6
      )
    )
  } else {
    p <- p + guides(fill = "none")
  }
  if (text) {
    p <- p + 
      geom_text(
        data = d_text,
        mapping = aes(label = text),
        size = text_size,
        x = -Inf,
        y = -Inf,
        hjust = -0.05,
        vjust = -0.5
      ) +
    scale_y_continuous(expand = expansion(mult = c(0.22, 0.04)))
  }
  return(p)
}

plot_relative_hex <- function(
  x, y, group = NULL,
  fun = "mean", bins = 101,
  legend = TRUE, hex = TRUE,
  palette = "davos", direction = -1,
  italic = FALSE
) {
  stopifnot(length(x) == length(y))
  stopifnot(length(x) == length(group))
  color_title <- "Fold Enrichment"
  d <- data.frame(x, y, group)
  stat_fun <- stat_summary_hex
  if (!hex) {
    stat_fun <- stat_summary_2d
  }
  if (length(palette) == 1) {
    plot_colors <- scico::scico(
      n = 20, palette = palette, direction = direction
    )
  } else {
    plot_colors <- palette
  }
  total_obs <- nrow(d)
  d$group <- as.integer(d$group)# != d$group[1])
  bg_ratio <- sum(d$group) / length(d$group)
  p <- ggplot(d) +
    stat_fun(
      mapping = aes(
        x = x,
        y = y,
        z = group
      ),
      # fun = fun, bins = bins
      bins = bins,
      fun = function(x) {
        (sum(x) / length(x)) / bg_ratio
      }
    ) +
    scale_fill_gradientn(
      colors = plot_colors,
      breaks = scales::pretty_breaks(3),
      limits = function(a) {
        if (a[2] - 1 > 1 - a[1]) {
          return(c(1 + 1 - a[2], a[2]))
        }
        return(c(a[1], 1 + 1 - a[1]))
      }
    ) +
    guides(fill = guide_colorbar(title = "Fold", barheight = 10)) +
    theme(
      legend.position = "right",
      # legend.position = c(0, 0),
      # legend.justification = c(0, 0),
      legend.background = element_blank(),
      plot.title = element_text(
        face = ifelse(italic, "italic", "plain"),
        size = 20
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  return(p)
}

plot_scattermore <- function(
  x, y, group,
  group_colors,
  group_labels = TRUE,
  group_legend = FALSE,
  alpha = 0.25,
  pixels = 500,
  pointsize = 1
) {
  n_groups <- length(unique(group))
  group_colors <- rep_len(group_colors, length.out = n_groups)
  d <- data.frame(
    x = x, y = y, group = naturalsort::naturalfactor(group),
    stringsAsFactors = FALSE
  )
  p <- ggplot() +
    scattermore::geom_scattermost(
      xy        = d,
      color     = grDevices::adjustcolor(group_colors, alpha.f = alpha)[d$group],
      pointsize = pointsize,
      pixels    = c(pixels, pixels)
    ) +
    theme(
      # plot.title = element_text(size = 32),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  if (group_legend) {
    d_count <- naturalsort::naturalfactor(append_n(as.character(d[["group"]])))
    p <- p +
      geom_point(
        # data = d[,c("group"),drop=FALSE] %>% unique,
        data = data.frame(group = d_count),
        mapping = aes(x = 1, y = 1, fill = group),
        size = 0, alpha = 0
      ) +
      scale_fill_manual(name = NULL, values = group_colors) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21, size = 5, alpha = 1))
      )
  }
  if (group_labels) {
    group_text <- d %>%
      group_by(group) %>%
      summarize(x = mean(x), y = mean(y), .groups = 'drop')
    group_text <- rbind(group_text, group_text %>% mutate(group = ""))
    set.seed(42)
    p <- p +
      ggrepel::geom_text_repel(
        data = group_text,
        mapping = aes(x = x, y = y, label = group),
        size = 6, color = "black", bg.color = "#ffffff66"
      )
  }
  return(p)
}


plot_hexmix <- function(
  x, y, group, bins = 101,
  group_colors,
  group_labels = TRUE,
  group_legend = FALSE
  #hex = TRUE
) {
  n_groups <- length(unique(group))
  group_colors <- rep_len(group_colors, length.out = n_groups)
  colormat <- grDevices::col2rgb(group_colors)
  hex <- hexbin::hexbin(x, y, xbins = bins, IDs = TRUE)
  d <- data.frame(
    x = x, y = y, group = group, hex = hex@cID,
    stringsAsFactors = FALSE
  )
  x <- d %>%
    dplyr::group_by(hex) %>%
    dplyr::count(group) %>%
    tidyr::spread(group, n) %>%
    as.data.frame
  rownames(x) <- x$hex
  x$hex <- NULL
  x[is.na(x)] <- 0
  x <- as.matrix(x)
  x <- ((x / rowSums(x)) %*% t(colormat)) / 256
  hex_color <- rgb(red = x[,1], green = x[,2], blue = x[,3])
  x <- cbind(
    as.data.frame(x),
    data.frame(hexbin::hcell2xy(hex), hex = hex@cell, count = hex@count)
  )
  x$hex_color <- hex_color
  # TODO
  # geom_fun <- geom_hex
  # if (!hex) {
  #   geom_fun <- geom_tile
  # }
  p <- ggplot() +
    geom_hex(
      data = x,
      mapping = aes(x, y, fill = hex),
      stat = "identity", color = NA, fill = x$hex_color
    ) +
    theme(
      # plot.title = element_text(size = 32),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  if (group_legend) {
    p <- p +
      geom_point(
        data = d[,c("group"),drop=FALSE] %>% unique,
        mapping = aes(x = 1, y = 1, fill = group),
        size = 0
      ) +
      scale_fill_manual(name = NULL, values = group_colors) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21, size = 5))
      )
  }
  if (group_labels) {
    group_text <- d %>%
      group_by(group) %>%
      summarize(x = mean(x), y = mean(y))
    group_text <- rbind(group_text, group_text %>% mutate(group = ""))
    set.seed(42)
    p <- p +
      # shadowtext::geom_shadowtext(
      #   data = group_text,
      #   mapping = aes(x = x, y = y, label = group),
      #   size = 8, color = "black", bg.color = "#ffffff33"
      # )
      ggrepel::geom_text_repel(
        data = group_text,
        mapping = aes(x = x, y = y, label = group),
        size = 6, color = "black", bg.color = "#ffffff66"
      )
  }
  return(p)
}

plot_umap_by_factor <- function(obs, this_factor, palette = "oslo") {
  # The dataframe needs to have these columns.
  stopifnot(all(c("UMAP1", "UMAP2", this_factor) %in% colnames(obs)))
  my_factor <- obs[[this_factor]] 
  plots <- lapply(naturalsort::naturalsort(unique(my_factor)), function(my_value) {
    ggplot(obs) +
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
          z = my_factor == my_value,
          group = -1
        ),
        bins = 71, color = NA
      ) +
      # guides(fill = guide_colorbar(barheight = 10)) +
      scale_fill_gradientn(
        guide = "none",
        name = NULL,
        limits = c(0.01, 1),
        colors = scico::scico(
          n = 20, palette = palette,
          direction = ifelse(palette %in% c("grayC", "bilbao", "lajolla"), 1, -1)
        ),
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
          "{my_value} (n = {scales::comma(sum(my_factor == my_value))})"
        )
      )
  })
  plots[[length(plots)]] <- plots[[length(plots)]] + 
    guides(fill = guide_colorbar(barwidth = 10)) +
    theme(legend.position = "bottom")
  ncol <- which.min(abs(sapply(1:16, function(i) {
    i - (16 / 9) * (length(plots) / i)
  })))
  p <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = sprintf("Proportion of each hexagon occupied by each %s", this_factor)
    )
  return(p)
}

# Pure functions --------------------------------------------------------------

#' Select rows from a sparse matrix, but fill with zeros if the rows are
#' absent.
#'
#' @param m A sparse matrix.
#' @param rows A vector of rownames.
#' @return A sparse matrix with selected rows.
#' @examples
#' set.seed(42)
#' m <- Matrix::rsparsematrix(nrow = 3, ncol = 5, density = 0.3)
#' rownames(m) <- letters[1:3]
#' select_rows(m, c("f", "b", "c", "d", "e", "a"))
#' #> 6 x 5 sparse Matrix of class "dgCMatrix"
#' #> 
#' #> f .    .     .    . .
#' #> b .    0.63  .    . .
#' #> c .    .    -0.11 . 0.4
#' #> d .    .     .    . .
#' #> e .    .     .    . .
#' #> a 0.36 .     .    . .
select_rows <- function(m, rows, fill = 0) {
  common_rows <- intersect(rownames(m), rows)
  new_rows <- setdiff(rows, rownames(m))
  empty <- Matrix::Matrix(
    data = fill,
    nrow = length(new_rows),
    ncol = ncol(m),
    sparse = TRUE
  )
  rownames(empty) <- new_rows
  return(rbind(m, empty)[rows,])
}


#' Append "(n = X)" to the levels of a factor or character vector.
#'
#' @param x A factor or character vector.
#' @return The same vector with new levels.
append_n <- function(x) {
  stopifnot(is.character(x) || is.factor(x))
  tab <- table(x)
  new_x <- sprintf(
      "%s (n = %s)",
      names(tab),
      scales::comma(as.integer(tab), accuracy = 1)
  )
  names(new_x) <- names(tab)
  return(new_x[x])
}

# Return a dataframe where 'col1' and 'col2' are factors with levels in order.
seriate_dataframe <- function(
  d, x, y, z, fun.aggregate = NULL, method = "BEA_TSP"
) {
  mat <- as.data.frame(reshape2::dcast(
    data          = d,
    formula       = as.formula(sprintf("%s ~ %s", x, y)),
    value.var     = z,
    fun.aggregate = fun.aggregate
  ))
  rownames(mat) <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  if (any(mat < 0)) {
    mat <- mat + abs(min(mat))
  }
  mat_order <- seriation::seriate(mat, method = method)
  d[[x]] <- factor(as.character(d[[x]]), rownames(mat)[mat_order[[1]]])
  d[[y]] <- factor(as.character(d[[y]]), colnames(mat)[mat_order[[2]]])
  return(d)
}

order_dataframe <- function(data, formula, value.var, method = "BEA_TSP") {
  x <- dcast.data.table(
    data = as.data.table(data),
    formula = as.formula(formula),
    value.var = value.var
  )
  x[is.na(x)] <- 0
  x <- as.data.frame(x)
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  x <- as.matrix(x)
  xo <- seriation::seriate(x, method = method)
  return(xo)
}

#' Order the rows of a dataframe with UMAP
#'
umap_sorted_rows <- function(data, formula, value.var, ...) {
  x <- dcast.data.table(
    data = as.data.table(data),
    formula = as.formula(formula),
    value.var = value.var
  )
  x[is.na(x)] <- 0
  x <- as.data.frame(x)
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  x <- as.matrix(x)
  n_neighbors <- ifelse(nrow(x) >= 15, 15, nrow(x) - 1)
  xo <- rownames(x)[order(uwot::umap(x, n_components = 1, n_neighbors = n_neighbors, ...)[,1])]
  return(xo)
}

corner <- function(x) {
  nr <- min(5, nrow(x))
  nc <- min(5, ncol(x))
  x[1:nr, 1:nc]
}

# Read about Bessel's correction: https://en.wikipedia.org/wiki/Bessel%27s_correction
# The sd() function already uses (n-1) as the denominator.
# So, here we also use (n-1) as the denominator.
sem <- function(x) sd(x) / sqrt(length(x) - 1)

# This function can read bed narrowPeak files from ENCODE.
read_encode_bed <- function(file) {
  retval <- readr::read_tsv(
    file = file,
    col_names = FALSE,
    col_types = 'ciicicdddi'
  )
  colnames(retval) <- c(
    "chrom", "chromStart", "chromEnd",
    "name", "score", "strand",
    "signalValue", "pValue", "qValue",
    "peak"
  )
  GenomicRanges::GRanges(
    seqnames = retval$chrom,
    ranges = IRanges::IRanges(
      start = retval$chromStart,
      end = retval$chromEnd,
      names = retval$name
    ),
    strand = stringr::str_replace(retval$strand, "\\.", "*"),
    score = retval$score,
    signalValue = retval$signalValue,
    pValue = retval$pValue,
    qValue = retval$qValue
  )
}

sigfig <- function(vec, digits = 2){
  retval <- gsub("\\.$", "",
    formatC(
      x      = signif(x =vec, digits = digits),
      digits = digits,
      format = "fg",
      flag   = "#"
    )
  )
  ix_small <- which(vec < 1e-3)
  if (length(ix_small) != 0L) {
    retval[ix_small] <- gsub("\\.$", "",
      formatC(
        x      = signif(x = vec[ix_small], digits = digits),
        digits = digits,
        format = "g",
        flag   = "#"
      )
    )
  }
  # Convert 1e-06 to 1e-6
  retval <- gsub("e-0(\\d)$", "e-\\1", retval)
  retval
}

p2z <- function(pval) qnorm(1 - (pval / 2))

scale_rows <- function(x, ...) t(scale(t(x), ...))

sort_hclust <- function(...) { as.hclust(dendsort(as.dendrogram(...))) }

cosine_normalize <- function (X, MARGIN = 1) {
  sweep(X, MARGIN, apply(X, MARGIN, function(x) sqrt(sum(x^2))), "/")
}

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

trim_image <- function(filename) {
  new_filename <- file.path(
    dirname(filename),
    sprintf("convert_trim_%s", basename(filename))
  )
  command <- sprintf(
    'convert "%s" -trim "%s" && mv -f "%s" "%s"',
    filename, new_filename,
    new_filename, filename
  )
  system(command)
}

# pdf_to_png <- function(pdf_file, density = 120) {
#   x <- magick::image_read(pdf_file, density = glue::glue("{density}x{density}"))
#   png_file <- stringr::str_replace(pdf_file, "pdf$", "png")
#   magick::image_write(x, path = png_file, format = "png")
# }

add_pngs <- function(dir, recursive = TRUE) {
  # dir <- "/home/ks38/work/github.com/slowkow/colitis/results/int3_tcell2_cd8_2-1/figures"
  pdf_files <- file.path(dir, list.files(path = dir, pattern = ".pdf$", recursive = recursive))
  pb <- txtProgressBar(min = 0, max = length(pdf_files), initial = 0, style = 3)
  for (i in seq_along(pdf_files)) {
    setTxtProgressBar(pb, i)
    pdf_file <- pdf_files[i]
    png_dir <- file.path(dirname(pdf_file), "png")
    dir.create(png_dir, showWarnings = FALSE)
    png_file <- file.path(png_dir,
      paste(stringr::str_remove(basename(pdf_file), '.pdf$'), '.png', sep = '')
    )
    if (!file.exists(png_file)) {
      res <- 200
      cmd <- glue::glue('gs -sDEVICE=pngalpha -o {shQuote(png_file)} -r{res} {shQuote(pdf_file)}')
      out <- system(cmd, intern = TRUE)
    }
  }
}

pdf_to_png <- function(file, res=200) {
  stopifnot(file.exists(file))
  png <- paste(stringr::str_remove(file, '.pdf$'), '.png', sep = '')
  if (!file.exists(png)) {
    cmd <- glue::glue('gs -sDEVICE=pngalpha -o {shQuote(png)} -r{res} {shQuote(file)}')
    system(cmd)
  }
  png
}

optimize_png <- function(filename) {
  opt_filename <- sprintf(
    "%s-fs8.png", substr(filename, 1, nchar(filename) - 4)
  )
  command <- sprintf(
    'pngquant --ext -fs8.png -- "%s" && mv -f "%s" "%s"',
    filename, opt_filename, filename
  )
  system(command)
}

ggsave_optimize_png <- function(filename, overwrite = TRUE, ...) {
  if (!overwrite && file.exists(filename)) {
    message("Not overwriting ", filename)
  } else {
    message(filename)
  }
  ggplot2::ggsave(filename, ...)
  n <- nchar(filename)
  if (substr(filename, n - 3, n) == ".png") {
    optimize_png(filename)
  }
}

sort_hclust <- function(...) as.hclust(dendsort::dendsort(as.dendrogram(...)))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

low_density <- function(x, y, n = 100, q = 0.05) {
  kde <- MASS::kde2d(x, y, n = n)
  kx <- cut(x, kde$x, labels = FALSE, include.lowest = TRUE)
  ky <- cut(y, kde$y, labels = FALSE, include.lowest = TRUE)
  kz <- sapply(seq_along(kx), function(i) kde$z[kx[i], ky[i]])
  kz < quantile(kz, q)
}

extreme_idx <- function(xs, low = 0.005, high = 0.995) {
  qs <- quantile(xs, probs = c(low, high))
  xs < qs[1] | xs > qs[2]
}

# Given an Ensembl ID like "ENSG00.1", return "ENSG00".
drop_version <- function(xs) stringr::str_split_fixed(xs, "\\.", 2)[,1]

# Strictly standardized mean difference.
# https://en.wikipedia.org/wiki/Strictly_standardized_mean_difference
ssmd <- function(x, y, independent = TRUE, na.rm = FALSE) {
  # (mean1 - mean2) / sqrt(sd1 ^ 2 + sd2 ^ 2)
  num <- mean(x, na.rm = na.rm) - mean(y, na.rm = na.rm)
  if (!independent) {
    den <- sqrt(
      var(x, na.rm = na.rm) +
        var(y, na.rm = na.rm) -
        2 * cov(x, y)
    )
  } else {
    den <- sqrt(var(x, na.rm = na.rm) + var(y, na.rm = na.rm))
  }
  num / den
}

#' Compute a hypergeometric p-value for your gene set of interest relative to
#' a universe of genes that you have defined.
#' 
#' @param ids A vector with genes of interest.
#' @param universe A vector with all genes, including the genes of interest.
#' @param gene_sets A list of gene vectors.
#' @return A dataframe with one row for each gene set, ordered by p-value.
do_hyper <- function(ids, universe, gene_sets) {
  gene_set_universe <- unique(unlist(gene_sets))
  ids <- ids[ids %in% gene_set_universe]
  universe <- universe[universe %in% gene_set_universe]
  retval <- as.data.frame(t(sapply(gene_sets, function(gene_set) {
    x <- sum(ids %in% gene_set)
    n <- length(ids)
    X <- sum(universe %in% gene_set)
    N <- length(universe)
    if (x > 0 && X > 0) {
      pval <- phyper(x - 1, X, N - X, n, lower.tail = FALSE)
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = (x / n) / (X / N),
        "pval" = pval
      )
    } else {
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = 0,
        "pval" = 1
      )
    }
    retvec
  })))
  retval$Name <- names(gene_sets)
  retval <- retval[order(retval$pval),]
  retval$qval <- p.adjust(retval$pval, method = "fdr")
  retval
}

do_fisher <- function(ids, universe, gene_sets) {
  retval <- as.data.frame(t(sapply(gene_sets, function(gene_set) {
    x <- sum(ids %in% gene_set)
    n <- length(ids)
    X <- sum(universe %in% gene_set)
    N <- length(universe)
    if (x > 0 && X > 0) {
      fish <- fisher.test(
        x = matrix(c(x, n - x, X - x, N - X - (n - x)), 2),
        alternative = "greater"
      )
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = (x / n) / (X / N),
        "orlow" = fish$conf.int[1],
        "orhigh" = fish$conf.int[2],
        "or" = unname(fish$estimate),
        "pval" = fish$p.value
      )
    } else {
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = 1,
        "orlow" = 1,
        "orhigh" = 1,
        "oddsratio" = 1,
        "pval" = 1
      )
    }
    retvec
  })))
  retval$Name <- names(gene_sets)
  retval <- retval[order(retval$pval),]
  retval$qval <- p.adjust(retval$pval, method = "fdr")
  retval
}

bonf_alpha <- function(alpha, n) 1 - (1 - alpha) ^ (1 / n)

# Given a matrix, choose the top rows for each column and return an index.
make_idx <- function(x, n = 1) { 
  xs <- unique(as.vector(apply(x, 2, function(xs) {
    tail(names(sort(xs)), n = n)
  })))
  rownames(x) %in% xs
}

#

# Plotting functions ----------------------------------------------------------

pow_trans <- function(power = 2) {
  force(power)
  trans <- function(x) x ^ power
  inv <- function(x) x ^ (1 / power)
  scales::trans_new(paste0("power-", format(power)), trans, inv)
}
               
theme_clean <- function(base_size) {
  theme_classic() +
  theme(
    panel.spacing    = unit(2, "lines"),
    panel.border     = element_rect(size = 0.5, fill = NA),
    axis.ticks       = element_line(size = 0.4),
    axis.line        = element_blank(),
    plot.title       = element_text(size = 16),
    strip.background = element_blank(),
    strip.text       = element_text(size = 16),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 16),
    axis.text        = element_text(size = 16),
    axis.title       = element_text(size = 16)
  )
}
  
# Overwrite default draw_colnames in the pheatmap package.©
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

# https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
get_breaks <- function(test, paletteLength) {
  myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
  # use floor and ceiling to deal with even/odd length palette lengths
  c(
    seq(min(test), 0, length.out = ceiling(paletteLength / 2) + 1), 
    seq(max(test) / paletteLength, max(test), length.out = floor(paletteLength / 2))
  )
}

# limma -----------------------------------------------------------------------

lmCompare <- function(object, fit1, fit2) {
  # Residuals
  res1 <- limma::residuals.MArrayLM(fit1, object)
  res2 <- limma::residuals.MArrayLM(fit2, object)
  # Residual sums of squares
  rss1 <- apply(res1, 1, function(xs) sum(xs ^ 2))
  rss2 <- apply(res2, 1, function(xs) sum(xs ^ 2))
  # Degrees of freedom
  df1 <- fit2$rank - fit1$rank
  df2 <- ncol(object) - fit2$rank
  # Number of parameters in each model
  if (df1 <= 0) {
    stop("fit2 must have more parameters than fit1")
  }
  # F statistics and p values
  #fs <- (rss1 - rss2) / (rss2 / df1)
  fs <- (rss1 - rss2) / rss2
  fs <- fs * df2 / df1
  ps <- pf(q = fs, df1 = df1, df2 = df2, lower.tail = FALSE)
  data.frame(
    "fstat" = fs,
    "pval"  = ps
  )
}

plot_limma_volcano <- function(
  dat, lfc = log2(1.5), fdr = 0.05, n_text = 10, classes = NULL,
  title = "", subtitle = "", text_ids = NULL,
  fdr_line = TRUE, seed = NA,
  expand_x = expansion(mult = 0.1),
  expand_y = expansion(mult = c(0.01, 0.05)),
  text_size = 4
) {
  # fdrs <- p.adjust(p = dat$P.Value, method = "fdr")
  fdrs <- dat$adj.P.Val
  p_threshold <- -log10(dat$P.Value[tail(which(fdrs < fdr), 1)])
  #  
  # n_genes <- sum(abs(dat$logFC) >= lfc & -log10(dat$P.Value) >= p_threshold)
  n_genes <- sum(abs(dat$logFC) >= lfc & fdrs < fdr)
  # 
  p <- ggplot(mapping = aes(x = logFC, y = -log10(P.Value))) +
    scattermore::geom_scattermore(
        data = dat,
        pointsize = 1,
        pixels = c(1224, 1024),
        color = "grey40"
    ) +
    geom_point(
      data = subset(dat, abs(logFC) >= lfc & adj.P.Val <= fdr),
      size = 0.4, color = "red"
    ) +
    scale_x_continuous(
      # breaks = scales::pretty_breaks(8),
      breaks = seq(-10, 10, by = 1),
      labels = function(x) fractional::fractional(2 ^ x),
      expand = expand_x
    ) +
    scale_y_continuous(expand = expand_y) +
    geom_vline(xintercept = 0, size = 0.3) +
    geom_vline(
      xintercept = c(-lfc, lfc),
      color = "grey70",
      size = 0.3, linetype = 2
    ) +
    labs(
      # x = bquote("Log"[2]~"Fold-Change"),
      x = "Fold-Change",
      y = bquote("-Log"[10]~"P"),
      title = title,
      subtitle = sprintf(
        "%s gene%s (%.2f FC, %s%% FDR)",
        comma(n_genes), ifelse(n_genes == 1, "", "s"), 2^lfc, signif(fdr * 100, 1)
      )
    )
  # if (fdr_line) {
  #   p <- p +
  #     geom_hline(
  #       yintercept = p_threshold,
  #       color = "grey70"
  #     ) +
  #     annotate(
  #       geom = "text",
  #       x = min(dat$logFC),
  #       y = p_threshold,
  #       hjust = 0,
  #       vjust = -0.3,
  #       color = "grey50",
  #       label = sprintf("%s%% FDR", signif(100 * fdr, 2))
  #     )
  # }
  # 
  if (!is.null(classes)) {
    # %1 percent in the tails
    q50 <- quantile(dat$logFC, c(0.005, 0.995))
    p <- p + annotate(
      geom = "text",
      fontface = "bold",
      # x = log2(c(1 / (lfc * 1.2), lfc * 1.2)),
      # x = range(dat$logFC),
      x = range(dat$logFC[dat$logFC > q50[1] & dat$logFC < q50[2]]),
      y = max(-log10(dat$P.Value)) * 1.15,
      label = classes,
      size = 5
    )
  }
  # 
  if (n_text > 0) {
    text_dat <- subset(dat, abs(logFC) >= lfc & adj.P.Val <= fdr)
    text_dat1 <- text_dat %>%
        filter(logFC > 0) %>%
        top_n(n = n_text / 2, wt = -log10(P.Value) * logFC)
    text_dat2 <- text_dat %>%
        filter(logFC < 0) %>%
        top_n(n = n_text / 2, wt = -log10(P.Value) * -logFC)
    text_dat$mylabel <- ""
    ix <- text_dat$Gene %in% c(text_dat1$Gene, text_dat2$Gene)
    text_dat$mylabel[ix] <- text_dat$Gene[ix]
    p <- p + 
      geom_point(
        data = text_dat %>% filter(mylabel != ""),
        mapping = aes(x = logFC, y = -log10(P.Value)),
        shape = 21, size = 1.25, color = "black", fill = "red"
      ) +
      geom_text_repel(
        data = text_dat,
        mapping = aes(x = logFC, y = -log10(P.Value), label = mylabel),
        # nudge_x = ifelse(text_dat$logFC > 0, 0.5, -0.5),
        size = text_size, fontface = "italic", seed = seed
      )
  } else if (!is.null(text_ids)) {
    text_dat <- dat[text_ids,]
    p <- p + 
      geom_point(
        data = text_dat,
        mapping = aes(x = logFC, y = -log10(P.Value)),
        shape = 21, size = 1.25, color = "black", fill = "red"
      ) +
      geom_text_repel(
        data = text_dat,
        mapping = aes(x = logFC, y = -log10(P.Value), label = Gene),
        size = text_size, fontface = "italic", seed = seed
      )
  }
  p <- p + annotate(
    geom = "text",
    x = -Inf,
    y = 0,
    label = glue("{comma(nrow(dat))} genes"),
    size = 5,
    hjust = -0.25, vjust = -0.5
  )
  return(p)
}


plot_limma_ma <- function(dat) {
  ggplot() +
    geom_point(
      data = dat,
      mapping = aes(x = AveExpr, y = logFC),
      size = 0.1
    ) +
    geom_point(
      data = subset(dat, abs(logFC) >= log2(1.5) & adj.P.Val <= 0.05),
      mapping = aes(x = AveExpr, y = logFC),
      size = 1, color = "red"
    ) +
    geom_text_repel(
      data = subset(dat, abs(logFC) >= log2(4) & adj.P.Val <= 0.05),
      mapping = aes(x = AveExpr, y = logFC, label = Gene),
      size = 3, color = "red"
    ) +
    geom_hline(
      yintercept = log2(c(1.5, 1, 0.6666)),
      color = "grey70"
    ) +
    scale_y_continuous(
      labels = function(x) signif(2^x, 2),
      breaks = pretty_breaks(5)
    ) +
    labs(
      x = bquote("Mean Log"[2]~"TPM"),
      y = bquote("Fold-Change")
    )
}

#' Create a quantile-quantile plot with ggplot2.
#'
#' Assumptions:
#'   - Expected P values are uniformly distributed.
#'   - Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Vector of p-values.
#' @param ci Size of the confidence interval, 95% by default.
#' @return A ggplot2 plot.
#' @examples
#' library(ggplot2)
#' gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

########################################################################
# github.com/immunogenomics/presto

compute_hash <- function(data_df, vars_use) {
    base <- 1
    hash <- rep(0, nrow(data_df))
    for (varname in vars_use) {
        vals <- factor(data.frame(data_df)[, varname, drop = TRUE])
        nlevel <- nlevels(vals)
        hash <- hash + (as.integer(vals) - 1) * base
        base <- base * nlevel
    }
    return(hash)
}

collapse_counts <- function(meta_data, varnames) {
    ## give each unique row a hash value for indexing
    hash <- compute_hash(meta_data, varnames)
    idx_keep <- which(!is.na(hash))
    hash <- hash[idx_keep]
    hash <- factor(sprintf('sample_%d', as.integer(hash)))
    meta_data <- meta_data[idx_keep, ]
    ## one hot encoded design matrix, sample level
    design_collapsed <- data.frame(meta_data)[, varnames, drop = FALSE] %>% 
        cbind(sample_id = hash) %>% 
        unique()
    design_collapsed <- data.table(meta_data)[
        , varnames, drop = FALSE, with = FALSE
    ][
        , sample_id := hash
    ][
        , N := .N, by = sample_id
    ] %>% 
    unique()
}
