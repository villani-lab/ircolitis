# composition.R
# Kamil Slowikowski

pacman::p_load(
  broom,
  broom.mixed,
  RhpcBLASctl,
  data.table,
  lme4
)

do_lm <- function(
  my_obs, form1, form2,
  mc.cores = 1, threads = 2, verbose = TRUE
) {
  RhpcBLASctl::blas_set_num_threads(threads)
  parallel::mclapply(
    mc.cores = mc.cores,
    X = unique(my_obs$cluster),
    FUN = function(this_cluster) {
      #message(this_cluster)
      mod0 <- NULL
      mod1 <- NULL
      mod_a <- NULL
      try({
        my_d <- my_obs %>%
        filter(cluster == this_cluster)
        mod0 <- lm(
          data    = my_d,
          formula = as.formula(form1)
        )
        mod1 <- lm(
          data    = my_d,
          formula = as.formula(form2)
        )
        mod_a <- anova(mod0, mod1)
      })
      return(list(
        "mod0" = mod0, "mod1" = mod1, "mod_a" = mod_a,
        "cluster" = this_cluster
      ))
    }
  )
}

#' Convert each result in the list returned by `do_lm()` into a table.
tidy_lm <- function(x) {
  lapply(x, function(y) {
    lapply(c("mod0", "mod1", "mod_a"), function(z) {
      retval <- suppressWarnings(tidy(y[[z]]))
      retval[["cluster"]] <- y[["cluster"]]
      retval
    })
  })
}

#' Summarize the results of `do_lm()` as a table.
summarize_lm <- function(res) {
  rbindlist(lapply(seq_along(res), function(i) {
    if (!is.null(res[[i]]$mod1)) {
			#message(i)
      x <- summary(res[[i]]$mod1)$coef
      x <- data.frame(
        coef = rownames(x),
        est = x[,1],
        sd = x[,2],
        t = x[,3],
        p = x[,4]
      )
      ci <- as.data.frame(confint(res[[i]]$mod1, method = "Wald"))
      colnames(ci) <- c("est_low", "est_high")
			ci$coef <- rownames(ci)
			x <- left_join(x, ci, by = "coef")
      x$lrt_p <- res[[i]]$mod_a[2,"Pr(>F)"]
      x <- pivot_longer(x, "coef")
      x$cluster <- res[[i]]$cluster
      x$name <- NULL
      return(x)
    }
  }))
}

#' Repeatedly call `do_lm()` with a shuffled version of one factor.
#' Useful to determine if p-values are calibrated.
repeat_lm <- function(
  x, form, shuffle_col = NA, n_iter = 10, mc.cores = 4, seed = 42
) {
  if (is.na(shuffle_col) || !shuffle_col %in% colnames(x)) {
    stop("shuffle_col should be the name of one column in x")
  }
  set.seed(seed)
  for (i in seq(n_iter)) {
    mycol <- sprintf("%s%s", shuffle_col, i)
    x[[mycol]] <- sample(x[[shuffle_col]])
  }
  my_res <- mclapply(
    X = seq(n_iter),
    mc.cores = mc.cores,
    FUN = function(i) {
      retval <- NULL
      form2 <- sprintf("%s + %s%s", form, shuffle_col, i)
      try({
        retval <- do_lm(x, form1 = form, form2 = form2, mc.cores = 1, threads = 2)
        retval <- tidy_lm(retval)
      })
      return(retval)
    }
  )
  # Reduce the glmer objects to summary tables.
  retval <- lapply(1:3, function(j) {
    rbindlist(lapply(seq_along(my_res), function(i) {
      items <- lapply(my_res[[i]], "[[", j)
      items <- items[!is.null(items)]
      items <- items[sapply(items, ncol) > 1]
      retval <- rbindlist(items)
      retval$iter <- i
      retval
    }))
  })
  names(retval) <- c("mod0", "mod1", "anova")
  retval
}

do_masc <- function(
  my_obs, clusters = NULL, form1, form2,
  mc.cores = 1, threads = 2, verbose = TRUE
) {
  RhpcBLASctl::blas_set_num_threads(threads)
  my_obs_clusters <- clusters
  if (is.null(my_obs_clusters)) {
    my_obs_clusters <- sort(unique(my_obs$leiden))
  }
  stopifnot(all(my_obs_clusters %in% my_obs$leiden))
  if (verbose) {
    print_status(glue("Running MASC with {nrow(my_obs)} cells and {length(unique(my_obs$leiden))} clusters"))
  }
  res <- parallel::mclapply(
    mc.cores = mc.cores,
    X = my_obs_clusters,
    FUN = function(this_cluster) {
      if (verbose) {
        print_status(glue("Testing cluster {this_cluster}"))
      }
      my_d <- my_obs %>%
        mutate(
          is_cluster = leiden == this_cluster
        )
      mod0 <- lme4::glmer(
        data    = my_d,
        # formula = is_cluster ~ 1 + chemistry + facs_sorting + (1|donor),
        formula = as.formula(form1),
        family  = binomial,
        nAGQ    = 2,
        control = lme4::glmerControl(optimizer = "bobyqa")
      )
      mod1 <- lme4::glmer(
        data    = my_d,
        # formula = is_cluster ~ is_case + chemistry + facs_sorting + (1|donor),
        formula = as.formula(form2),
        family  = binomial,
        nAGQ    = 2,
        control = lme4::glmerControl(optimizer = "bobyqa")
      )
      mod_a <- anova(mod0, mod1)
      list(
        "mod0" = mod0, "mod1" = mod1, "mod_a" = mod_a,
        "cluster" = this_cluster
      )
    }
  )
  if (verbose) {
    print_status(glue("Finished MASC"))
  }
  return(res)
}

#' Convert each result in the list returned by `do_masc()` into a table.
tidy_masc <- function(x) {
  lapply(x, function(y) {
    lapply(c("mod0", "mod1", "mod_a"), function(z) {
      retval <- suppressWarnings(tidy(y[[z]]))
      retval[["cluster"]] <- y[["cluster"]]
      retval
    })
  })
}

#' Summarize the results of `do_masc()` as a table.
summarize_masc <- function(res) {
  rbindlist(lapply(seq_along(res), function(i) {
    x <- summary(res[[i]]$mod1)$coef
    x <- data.frame(
      coef = rownames(x),
      est = x[,1],
      sd = x[,2],
      z = x[,3],
      p = x[,4]
    )
    ci <- confint.merMod(res[[i]]$mod1, method = "Wald")
    ci <- ci[2:nrow(ci),]
    colnames(ci) <- c("est_low", "est_high")
    x <- cbind(x, ci)
    x$lrt_p <- res[[i]]$mod_a[2,"Pr(>Chisq)"]
    x <- pivot_longer(x, "coef")
    x$cluster <- res[[i]]$cluster
    x$name <- NULL
    return(x)
  }))
}

#' Repeatedly call `do_masc()` with a shuffled version of one factor. 
#' Useful to determine if p-values are calibrated.
repeat_masc <- function(
  x, form, shuffle_col = NA, n_iter = 10, mc.cores = 4, seed = 42
) {
  if (is.na(shuffle_col) || !shuffle_col %in% colnames(x)) {
    stop("shuffle_col should be the name of one column in x")
  }
  set.seed(seed)
  for (i in seq(n_iter)) {
    mycol <- sprintf("%s%s", shuffle_col, i)
    x[[mycol]] <- sample(x[[shuffle_col]])
  }
  my_res <- mclapply(X = seq(n_iter), mc.cores = mc.cores, FUN = function(i) {
    form2 <- sprintf("%s + %s%s", form, shuffle_col, i)
    retval <- do_masc(x, form1 = form, form2 = form2, mc.cores = 1, threads = 2)
    tidy_masc(retval)
  })
  # Reduce the glmer objects to summary tables.
  retval <- lapply(1:3, function(j) {
    rbindlist(lapply(seq_along(my_res), function(i) {
      retval <- rbindlist(lapply(my_res[[i]], "[[", j))
      retval$iter <- i
      retval
    }))
  })
  names(retval) <- c("mod0", "mod1", "anova")
  retval
}


plot_masc <- function(xx, max_fdr = 0.1) {
  xx$fdr <- p.adjust(xx$p)
  xx$lrt_fdr <- p.adjust(xx$lrt_p)
  xx$est <- exp(xx$est)
  xx$est_low <- exp(xx$est_low)
  xx$est_high <- exp(xx$est_high)
  ggplot() +
    geom_vline(xintercept = 1, size = 0.2) +
    geom_errorbarh(
      data = xx,#[xx$lrt_fdr < max_fdr,],
      mapping = aes(y = -log10(p), xmin = est_low, xmax = est_high),
      height = 0
    ) +
    geom_point(
      data = xx,
      mapping = aes(x = est, y = -log10(p), color = lrt_fdr < max_fdr),
      size = 3
    ) +
    geom_text_repel(
      data = xx[xx$lrt_fdr < max_fdr,],
      mapping = aes(x = est, y = -log10(p), label = cluster),
      size = 5
    ) +
    scale_color_manual(
      values = c("grey30", "red"),
      labels = c(
        glue("FDR > {max_fdr}"),
        glue("FDR < {max_fdr}")
      ),
      name = NULL
    ) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(x = "Odds Ratio", y = bquote("-Log"[10]~"P"))
}


