run_analysis <- function(obs, counts, exclude_genes, mito_genes, params) {
  n_cells <- ncol(counts)
  default_params <- list(
    analysis_name = "analysis",
    min_percent   = 100 * (50 / n_cells),
    loess_span    = 0.01,
    n_pcs         = "mcv",
    max_pcs       = 50,
    equal_pcs     = TRUE,
    n_harmony     = 0,
    harmony_vars  = c("batch"),
    n_knn         = round(sqrt(n_cells)),
    leiden_res    = 1.3,
    leiden_iter   = 20,
    umap_spread   = 1,
    umap_min_dist = 0.05,
    log_file      = ""
  )
  for (key in names(params)) {
    if (key %in% names(default_params)) {
      # stopifnot(class(default_params[[key]]) == class(params[[key]]))
      # if (class(default_params[[key]]) != class(params[[key]])) {
      #   stop(glue("class for param {key} does not match"))
      # }
      default_params[[key]] <- params[[key]]
    }
  }
  do_analysis(
    obs           = obs,
    counts        = counts,
    exclude_genes = exclude_genes,
    mito_genes    = mito_genes,
    analysis_name = default_params[["analysis_name"]],
    min_percent   = default_params[["min_percent"]],
    loess_span    = default_params[["loess_span"]],
    n_pcs         = default_params[["n_pcs"]],
    max_pcs       = default_params[["max_pcs"]],
    equal_pcs     = default_params[["equal_pcs"]],
    n_harmony     = default_params[["n_harmony"]],
    harmony_vars  = default_params[["harmony_vars"]],
    n_knn         = default_params[["n_knn"]],
    leiden_res    = default_params[["leiden_res"]],
    leiden_iter   = default_params[["leiden_iter"]],
    umap_spread   = default_params[["umap_spread"]],
    umap_min_dist = default_params[["umap_min_dist"]],
    log_file      = default_params[["log_file"]]
  )
}

#' Main function for running a complete analysis
#' @param obs
#' @param log2cpm
#' @param min_percent
#' @param n_genes Number of robust genes to use for PCA.
#' @param n_pcs
#' @param n_harmony Max number of Harmony iterations.
#' @param n_knn
#' @param leiden_res
#' @param leiden_iter
do_analysis <- function(
  obs, counts,
  exclude_genes = NULL,
  analysis_name = "",
  mito_genes    = NULL,
  loess_span    = 0.05,
  min_percent   = 0.5,
  # n_genes       = 2000,
  n_pcs         = 25,
  max_pcs       = 50,
  equal_pcs     = TRUE,
  harmony_vars  = c("channel"),
  n_harmony     = 20,
  n_knn         = 30,
  leiden_res    = 1.3,
  leiden_iter   = 10,
  dist_method   = "euclidean",
  umap_spread   = 1,
  umap_min_dist = 0.01,
  umap_threads  = 8,
  umap_epochs   = 800,
  random_seed   = 42,
  log_file      = ""
) {
  # browser()
  if (ncol(counts) != nrow(obs)) {
    stop("Columns in counts must match rows in obs")
  }
  cat("", file = log_file, append = FALSE)
  cat(analysis_name, file = log_file, append = FALSE)
  log_message(glue::glue(
    "Computing cell statistics for {scales::comma(ncol(counts))} cells"
  ), log_file)
  obs$n_counts    <- Matrix::colSums(counts)
  obs$n_features  <- Matrix::colSums(counts > 0)
  ix_mito         <- which(rownames(counts) %in% mito_genes)
  obs$mito_counts <- colSums(counts[ix_mito,])
  obs$mito_pct    <- 100 * obs$mito_counts / obs$n_counts
  log_message(glue::glue(
    "Computing gene statistics for {scales::comma(nrow(counts))} genes"
  ), log_file)
  my_loess <- do_loess(counts, exclude_genes, loess_span, min_percent)
  log_message(glue::glue(
    "{scales::comma(sum(my_loess$counts_stats$include))} genes have expression in > {signif(min_percent, 3)}% of cells"
  ), log_file)
  ix_genes <- which(my_loess$counts_stats$residuals > 0)
  log_message(glue::glue(
    "{scales::comma(length(ix_genes))} genes have residual variance > 0 for the model log10(sd) ~ log10(mean)"
  ), log_file)
  # Take the top 80 percent
  x <- my_loess$counts_stats$residuals[ix_genes]
  ix_genes <- ix_genes[x > quantile(x, 0.2)]
  # ix_genes <- ix_genes[x > quantile(x, 0.1)]
  # Log2CPM
  log2cpm <- do_log2cpm(counts, total = median(colSums(counts)))
  # Molecular cross validation
  mcv <- NA
  if (n_pcs == "mcv") {
    log_message(glue::glue(
      "Running MCV with {scales::comma(length(ix_genes))} genes and {scales::comma(ncol(log2cpm))} cells, max PCs {max_pcs}"
    ), log_file)
    set.seed(random_seed)
    mcv <- do_mcv_pca(
      counts   = counts,
      ix_genes = ix_genes,
      max_pcs  = max_pcs,
      mcv_reps = 5
    )
    n_pcs <- as.integer(mcv$optimal_k[1])
    log_message(glue::glue("MCV found {n_pcs} PCs is optimal"), log_file)
  }
  if (n_pcs > 0) {
    log_message(glue::glue(
      "Running PCA with {scales::comma(length(ix_genes))} genes and {scales::comma(ncol(log2cpm))} cells, keep {n_pcs} PCs"
    ), log_file)
    set.seed(random_seed)
    pca <- RSpectra::svds(
      A    = t(log2cpm[ix_genes,]),
      k    = n_pcs,
      opts = list(
        center = TRUE,
        scale  = TRUE,
        # scale = counts_stats$sd[ix_genes] - counts_stats$residuals[ix_genes],
        maxitr = 2000,
        tol    = 1e-10
      )
    )
    if (equal_pcs) {
      pca$x <- pca$u
    } else {
      # By default, weight PCs by their eigenvalues.
      pca$x <- pca$u %*% diag(pca$d)
    }
    pca$genes <- rownames(log2cpm)[ix_genes]
    if (n_harmony > 0) {
      log_message(glue::glue("Running Harmony with {n_pcs} PCs"), log_file)
      hm <- harmony::HarmonyMatrix(
        data_mat         = t(pca$x),
        meta_data        = obs,
        vars_use         = harmony_vars,
        max.iter.harmony = n_harmony,
        do_pca           = FALSE,
        return_object    = TRUE
      )
      pca_h <- as.matrix(t(hm$Z_corr))
      log_message(glue::glue(
        "Finding {n_knn} nearest neighbors with {dist_method} on {n_pcs} harmonized PCs"
      ), log_file)
    } else {
      hm <- NULL
      pca_h <- pca$x
      log_message(glue::glue(
        "Finding {n_knn} nearest neighbors with {dist_method} on {n_pcs} PCs"
      ), log_file)
    }
    colnames(pca_h) <- sprintf("PC%s", 1:ncol(pca_h))
    rownames(pca_h) <- colnames(log2cpm)
    #
    # Add PCA coordinates
    obs <- cbind(obs, pca_h)
    # KNN
    knn <- do_knn(pca_h, n_knn, dist_method)
    # snn <- compute_snn(1.0 * (knn$simil_cells > 0))
    if (n_harmony > 0) {
      log_message(glue::glue("Running UMAP with KNN from harmonized PCs"), log_file)
    } else {
      log_message(glue::glue("Running UMAP with KNN from PCs"), log_file)
    }
    # Add UMAP coordinates
    pca_h_umap <- uwot::umap(
      X         = NULL,
      init      = pca_h[,1:2],
      nn_method = knn$nn_method,
      # nn_method = list(dist = 1 / snn, idx = 1.0 * (snn > 0)),
      spread    = umap_spread,
      min_dist  = umap_min_dist,
      n_threads = umap_threads,
      n_epochs  = umap_epochs
    )
    obs$UMAP1 <- pca_h_umap[,1]
    obs$UMAP2 <- pca_h_umap[,2]
  } else {
    # No PCA
    X <- as.matrix(scale(t(log2cpm[ix_genes,])))
    knn <- do_knn(X, n_knn, dist_method)
    # snn <- compute_snn(1.0 * (knn$simil_cells > 0))
    log_message(glue::glue("Running UMAP with {length(ix_genes)} features"), log_file)
    # Add UMAP coordinates
    log2cpm_umap <- uwot::umap(
      X         = NULL,
      nn_method = knn$nn_method,
      # nn_method = list(dist = 1 / snn, idx = 1.0 * (snn > 0)),
      spread    = umap_spread,
      min_dist  = umap_min_dist,
      n_threads = umap_threads,
      n_epochs  = umap_epochs
    )
    obs$UMAP1 <- log2cpm_umap[,1]
    obs$UMAP2 <- log2cpm_umap[,2]
  }
  # leiden_res <- c(0.5, 1.0, 1.5)
  # leiden_res <- 1.3
  for (my_res in leiden_res) {
    log_message(glue::glue(
      "Running Leiden community detection on KNN with resolution {signif(my_res, 3)} and {leiden_iter} iterations"
    ), log_file)
    # leiden <- run_leiden(
    #   adj        = knn$simil_cells,
    #   resolution = my_res,
    #   iterations = leiden_iter
    # )
    # my_leiden <- as.integer(1 + unname(t(leiden)[,1])[3:ncol(leiden)])
    # obs[[glue("leiden{my_res}")]] <- my_leiden
    my_leiden <- run_leiden(
      adj        = knn$simil_cells,
      resolution = my_res,
      iterations = leiden_iter,
      seed       = random_seed
    )
    obs[[glue("leiden{signif(my_res, 3)}")]] <- my_leiden
  }
  # Pick the last one as the "main" clustering
  obs[["leiden"]] <- my_leiden
  #
  # # Try to guess the best resolution automatically # 1 hour for 25k cells
  #   log_message(glue::glue(
  #     "Using MCV to find optimal Leiden resolution"
  #   ), log_file)
  #   mcv_leiden <- do_mcv_leiden(
  #     counts, ix_genes,
  #     n_pcs       = n_pcs,
  #     n_knn       = n_knn,
  #     res_range   = leiden_res,
  #     leiden_reps = leiden_iter,
  #     mcv_reps    = 5
  #   )
  #   log_message(glue::glue(
  #     "Optimal resolution is {mcv_leiden$optimal_res[1]}"
  #   ), log_file)
  #
  # log_message(glue::glue(
  #   "Running Leiden community detection on KNN with resolution {leiden_res} and {leiden_iter} iterations"
  # ), log_file)
  # leiden <- run_leiden(
  #   adj        = knn$simil_cells,
  #   resolution = leiden_res,
  #   iterations = leiden_iter
  # )
  # obs$leiden <- as.integer(1 + unname(t(leiden)[,1])[3:ncol(leiden)])
  #
  # Build SNN
  # log_message(glue::glue("Building rank-based SNN from KNN"), log_file)
  # snn                 <- scran:::build_snn_rank(knn$index)
  # g                   <- igraph::make_graph(edges = snn[[1]])
  # igraph::E(g)$weight <- snn[[2]]
  # g                   <- igraph::simplify(g, edge.attr.comb = "first")
  # snn_adj             <- igraph::as_adjacency_matrix(
  #   g, attr = "weight", sparse = TRUE
  # )
  # log_message(glue::glue(
  #   "Running Leiden community detection on SNN with resolution {leiden_res} and {leiden_iter} iterations"
  # ), log_file)
  # obs$leiden_snn <- run_leiden(
  #   adj        = snn_adj,
  #   resolution = leiden_res,
  #   iterations = leiden_iter
  # )
  log_message("done", log_file)
  if (n_pcs > 0) {
    list(
      obs           = obs,
      counts        = counts,
      ix_genes      = ix_genes,
      exclude_genes = exclude_genes,
      counts_stats  = my_loess$counts_stats,
      fit           = my_loess$fit,
      mcv           = mcv,
      pca           = pca,
      pca_h         = pca_h,
      hm            = hm,
      knn           = knn,
      # leiden        = leiden,
      # snn           = snn,
      de            = presto::wilcoxauc(log2cpm, obs$leiden)
    )
  } else {
    list(
      obs           = obs,
      counts        = counts,
      ix_genes      = ix_genes,
      exclude_genes = exclude_genes,
      counts_stats  = my_loess$counts_stats,
      fit           = my_loess$fit,
      # pca           = pca,
      # pca_h         = pca_h,
      # hm            = hm,
      knn           = knn,
      # leiden        = leiden,
      # snn           = snn,
      de            = presto::wilcoxauc(log2cpm, obs$leiden)
    )
  }
}

log_message <- function(x, log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  this_message <- glue::glue("{timestamp}\t{x}")
  if (log_file != "") {
    cat(this_message, file = log_file, append = TRUE)
    cat("\n", file = log_file, append = TRUE)
  }
  message(this_message)
}

do_cpm <- function(A, total = 1e4) {
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  A@x <- total * A@x
  return(A)
}

do_log2cpm <- function(A, total = 1e4) {
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  A@x <- total * A@x
  A@x <- log2(A@x + 1)
  return(A)
}

do_loess <- function(counts, exclude_genes, loess_span, min_percent) {
  d <- data.table::data.table(
    mean    = Matrix::rowMeans(counts),
    sd      = proxyC::rowSds(counts),
    percent = 100 * Matrix::rowSums(counts > 0) / ncol(counts)
  )
  d$gene <- rownames(counts)
  d$exclude <- d$gene %in% exclude_genes
  d$include <- (
    d$percent >= min_percent & !d$exclude
  )
  ix_include <- d$include
  # Need to be careful that this looks good
  fit <- loess(
    formula = log10(sd) ~ log10(mean),
    data    = d[ix_include,],
    span    = loess_span,
    degree  = 2
  )
  # plot(fit$x, fit$residuals)
  d$fitted <- NA
  d$fitted[ix_include] <- fit$fitted
  d$residuals <- NA
  d$residuals[ix_include] <- fit$residuals
  d$rank <- NA
  d$rank[ix_include] <- (
    rank(rank(-fit$residuals) + rank(-fit$y / fit$fitted))
  )
  # p <- ggplot() +
  #   stat_binhex(
  #     data = d,
  #     mapping = aes(log10(mean), residuals),
  #     # mapping = aes(log10(mean), log10(sd)),
  #     # mapping = aes(log10(mean), fitted),
  #     size = 1, bins = 131
  #   ) +
  #   # geom_point(
  #   #   data = d,
  #   #   mapping = aes(log10(mean), fitted),
  #   #   size = 0.1
  #   # ) +
  #   geom_hline(yintercept = 0, size = 0.3, linetype = 2) +
  #   scale_fill_gradientn(
  #     colors = scico::scico(20)[4:20],
  #     trans = "log10",
  #     breaks = scales::log_breaks(7)
  #   ) +
  #   guides(
  #     fill = guide_colorbar(barheight = 10)
  #   ) +
  #   labs(
  #     x = bquote("Log"[10]~"Mean"),
  #     y = bquote("Log"[10]~"SD (resid.)"),
  #     title = glue::glue("Select {comma(sum(d$include))} of {comma(nrow(d))} total genes"),
  #     fill = "Genes"
  #   )
  # p
  return(list(counts_stats = d, fit = fit))
}

# snn <- compute_snn(1.0 * (knn$simil_cells > 0))
compute_snn <- function(adj, prune = 1/20) {
  snn <- Matrix::tcrossprod(adj)
  k <- nrow(adj)
  snn@x[ (snn@x / (k + (k - snn@x)) < prune)] <- 0
  snn
}

do_knn <- function(pca_h, n_knn, dist_method = "euclidean") {
  if (dist_method == "euclidean") {
    knn <- BiocNeighbors::findKNN(
      X       = pca_h,
      k       = n_knn,
      BNPARAM = BiocNeighbors::HnswParam(),
      BPPARAM = BiocParallel::MulticoreParam(workers = 4)
    )
    simil_cells <- Matrix::sparseMatrix(
      i = rep(1:nrow(knn$index), n_knn),
      j = as.integer(knn$index),
      x = as.numeric(knn$distance)
    )
    return(list(
      simil_cells = simil_cells,
      nn_method = list(
        idx = knn$index,
        dist = knn$distance
      )
    ))
  } else {
    simil_cells <- proxyC::simil(
      x      = Matrix::Matrix(pca_h, sparse = TRUE),
      margin = 1,
      method = dist_method, # correlation, cosine
      rank   = n_knn + 1
    )
    d <- summary(simil_cells)
    d <- d[order(d$j, d$i),]
    nn_method <- list(
      idx = matrix(d$i, ncol = n_knn + 1, byrow = TRUE)[,2:n_knn],
      dist = matrix(1 - d$x, ncol = n_knn + 1, byrow = TRUE)[,2:n_knn]
    )
    return(list(
      simil_cells = simil_cells,
      nn_method = nn_method
    ))
  }
}

  # a2$obs$cluster <- factor(a2$obs$leiden)
  # pc_lm <- rbindlist(lapply(seq(ncol(a2$pca$x)), function(my_pc) {
  #   my_form <- glue::glue("PC{my_pc} ~ 0 + cluster")
  #   res <- broom::tidy(summary(lm(as.formula(my_form), a2$obs)))
  #   # res <- broom::tidy(anova(aov(as.formula(my_form), a2$obs)))[1,]
  #   res$pc <- my_pc
  #   res
  # }))

#' @param x Sparse matrix (e.g. gene expression)
#' @param clusters Factor of cluster assignments
cluster_means <- function(x, clusters) {
  stopifnot(ncol(x) == length(clusters))
  clusters <- factor(clusters)
  y <- model.matrix(~ 0 + clusters)
  colnames(y) <- substr(colnames(y), 9, nchar(colnames(y)))
  y <- sweep(y, 2, colSums(y), "/")
  x %*% y
}

##' @importFrom limma lmFit eBayes makeContrasts contrasts.fit topTable
#pseudobulk_markers <- function(counts, cluster, donor) {
#  # Make a pseudobulk matrix
#  y <- model.matrix(~ 0 + factor(cluster):factor(donor))
#  pb <- as(counts %*% y, "dgCMatrix")
#  pb <- do_log2cpm(pb, median(colSums(pb)))
#  #
#  library(limma)
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
#    sample_info %>%
#      select(donor, case, chemistry, qubit_library_quantification_ng_ul) %>%
#      group_by(donor, case, chemistry) %>%
#      summarize_if(is.numeric, mean),
#    by = "donor"
#  )
#  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
#  stopifnot(nrow(pb_meta) == ncol(pb))
#  # Test all pairs of clusters
#  print_status("Finding marker genes with pseudobulk")
#  pb_meta$x <- factor(pb_meta$cluster)
#  des1 <- with(pb_meta, model.matrix(~ 0 + x))
#  fit1 <- lmFit(object = as.matrix(pb[rowMeans(pb) > 0.5,]), design = des1)
#  fit1 <- eBayes(fit1)
#  fit1$genes <- ensembl_to_symbol[rownames(fit1$coefficients)]
#  cluster_pairs <- t(combn(levels(pb_meta$x), 2))
#  cont <- makeContrasts(contrasts = lapply(seq(nrow(cluster_pairs)), function(i) {
#    glue("x{cluster_pairs[i,1]} - x{cluster_pairs[i,2]}")
#  }), levels = des1)
#  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
#  fit2 <- contrasts.fit(fit1, cont)
#  fit2 <- eBayes(fit2)
#  res <- rbindlist(lapply(colnames(cont), function(this_coef) {
#    x <- topTable(fit2, coef = this_coef, number = nrow(fit1$coefficients))
#    this_coef <- str_replace_all(this_coef, "x", "")
#    this_coef <- str_replace(this_coef, "vs", " vs ")
#    x$coef <- this_coef
#    x
#  }))
#}

##' @param x matrix of objects (rows) and attributes (columns)
#get_knn <- function(x, n_neighbors = 30) {
#  n_objects <- nrow(x)
#  n_neighbors <- n_neighbors + 1 # self is a neighbor
#  idx <- RANN::nn2(x, k = n_neighbors)$nn.idx[,2:(n_neighbors)]
#  Matrix::sparseMatrix(rep(1:n_objects, n_neighbors - 1), as.integer(idx))
#}

# get_snn <- function(knn_adj, prune_snn = NULL) {
#   n_knn <- sum(knn_adj[1,])
#   knn_adj <- as(knn_adj, "dgCMatrix")
#   snn <- Matrix::tcrossprod(knn_adj)
#   # snn@x <- snn@x / (2 * knn_n - snn@x)
#   # if (is.null(prune_snn)) {
#   #   snn@x[snn@x < prune_snn] <- 0
#   # }
#   snn
# }


#' @param adj adjacency matrix
run_leiden <- function(adj, resolution = 1.3, iterations = 3, seed = 1, invert_dist = FALSE) {
  stopifnot(
    is(adj, "dgCMatrix") || is(adj, "dgTMatrix")
  )
  d <- summary(adj)
  if (invert_dist) {
    d$x <- 1 / d$x
  }
  stopifnot(all(d$x > 0))
  # Sort the i,j columns.
  d$i <- as.integer(d$i - 1)
  d$j <- as.integer(d$j - 1)
  ix <- d$i > d$j
  d$i2 <- d$i
  d$i2[ix] <- d$j[ix]
  d$j2 <- d$j
  d$j2[ix] <- d$i[ix]
  d <- d[order(d$i2, d$j2),]
  d <- unique(d[,c("i2", "j2","x")])
  # Save to a temporary file.
  infile <- tempfile()
  data.table::fwrite(d, infile, sep = "\t", col.names = FALSE)
  py <- reticulate::py_run_string(convert = FALSE, code = "
def leiden(infile, resolution, iterations, seed = 1):
    import numpy as np
    from igraph import Graph
    import leidenalg
    g = Graph.Read_Ncol(infile, directed=False)
    vertex_names = np.array([int(x) for x in g.vs.get_attribute_values('name')])
    o = np.argsort([int(x) for x in vertex_names])
    partition = leidenalg.find_partition(
        graph = g,
        weights = 'weight',
        partition_type = leidenalg.RBConfigurationVertexPartition,
        resolution_parameter = float(resolution),
        n_iterations = int(iterations),
        seed = int(seed)
    )
    groups = np.array(partition.membership)[o]
    return groups
")
  clusters <- 1 + reticulate::py_to_r(py$leiden(infile, resolution, iterations, seed))
  return(clusters)
}

##' @param adj adjacency matrix
#run_leiden_old <- function(adj, resolution = 1.3, iterations = 3, invert_dist = FALSE) {
#  # An adjacency matrix should only have zeros and ones.
#  #stopifnot(is(adj, "ngCMatrix"))
#  stopifnot(
#    is(adj, "dgCMatrix") || is(adj, "dgTMatrix")
#  )
#  d <- summary(adj)
#  if (invert_dist) {
#    d$x <- 1 / d$x
#  }
#  stopifnot(all(d$x > 0))
#  d$i <- as.integer(d$i - 1)
#  d$j <- as.integer(d$j - 1)
#  ix <- d$i > d$j
#  d$i2 <- d$i
#  d$i2[ix] <- d$j[ix]
#  d$j2 <- d$j
#  d$j2[ix] <- d$i[ix]
#  d <- d[order(d$i2, d$j2),]
#  d <- unique(d[,c("i2", "j2","x")])
#  # Make temporary files
#  infile <- tempfile()
#  outfile <- tempfile()
#  data.table::fwrite(d, infile, sep = "\t", col.names = FALSE)
#  print(glue::glue("infile = {infile}"))
#  print(glue::glue("outfile = {outfile}"))
#  # browser()
#  # Java version does not seem to work...
#  # resolution <- 1.3
#  # algorithm <- "Leiden"
#  # iterations <- 3
#  # my_seed <- 42
#  # randomness <- 0.1
#  # cmd <- glue::glue(
#  #   "java -jar /Users/kamil/src/RunNetworkClustering.jar --seed {my_seed} --randomness {randomness} -r {resolution} -i {iterations} -a {algorithm} -o {outfile} {infile}"
#  # )
#  # Python version seems to be ok
#  # resolution <- 1.3
#  # iterations <- 3
#  cmd <- glue::glue(
#    "python scripts/run-leiden.py {infile} {outfile} --resolution {resolution} --iterations {iterations}"
#  )
#  print(glue::glue("Running system command: {cmd}"))
#  system(cmd, wait = TRUE)
#  # out <- data.table::fread(outfile, header = FALSE)$V1 + 1
#  # outfile <- "~/testing.txt"
#  out <- data.table::fread(outfile, header = FALSE)
#  # unlink(infile)
#  # unlink(outfile)
#  return(out)
#}

#' Normalize each row of a sparse matrix.
#'
#' Subtract the mean from each row. Then divide each row by the standard
#' deviation.
#'
#' @param X A sparse matrix.
#' @param axis Scale rows if axis=1, or columns if axis=0.
scale_data <- function(X, axis = 1) {
  if (axis == 0) {
    X <- t(X)
  }
  X_mean <- Matrix::rowMeans(X)
  X_std <- proxyC::rowSds(X)
  X <- as.matrix(X - X_mean)
  X <- X / X_std
  X[is.na(X)] <- 0
  if (axis == 0) {
    X <- t(X)
  }
  return(X)
}

#' Molecular Cross Validation (MCV) in R
#'
#' @param counts Matrix of genes (rows) and cells (columns).
#' @param selected Index of which genes to include.
#' @param max_pcs Maximum number of PCs to test.
#'
do_mcv <- function(counts, ix_genes, max_pcs = 50) {
	#
  X1 <- counts
  X2 <- counts
  X1@x <- as.numeric(rbinom(n = length(counts@x), size = counts@x, prob = 0.5))
	stopifnot(all(counts@x >= X1@x))
  X2@x <- counts@x - X1@x
  stopifnot(all(X1@x + X2@x == counts@x))
  X1 <- do_log2cpm(X1, total = median(colSums(X1)))
  X2 <- do_log2cpm(X2, total = median(colSums(X2)))
	#
	ix_genes <- ix_genes[rowSums(X1[ix_genes,]) > 0]
	# X1 <- t(X1[ix_genes,])
	# X2 <- t(X2[ix_genes,])
	X1 <- t(scale_data(X1[ix_genes,]))
	X2 <- t(scale_data(X2[ix_genes,]))
	#
  pca <- RSpectra::svds(
    A    = X1,
    k    = max_pcs,
    opts = list(
      center = FALSE, # can be true
      scale  = FALSE, # should be false
      maxitr = 2000,
      tol    = 1e-10
    )
  )
  pca$genes <- rownames(X1)[ix_genes]
  pca$x <- pca$u %*% diag(pca$d)
	# X1 <- scale(X1) # this doesn't seem to work with scale=TRUE above
	# X2 <- scale(X2)
  # X2[is.na(X2)] <- 0
	k_range <- c(seq(2, 10, 1), seq(11, 29, 2), seq(30, max_pcs, 2))
  k_range <- k_range[k_range <= max_pcs]
  mcv_loss <- numeric(length(k_range))
  rec_loss <- numeric(length(k_range))
	mean_squared_error <- function(x, y) mean( (x - y) ^ 2 , na.rm = TRUE)
	for (i in seq_along(k_range)) {
		k <- k_range[i]
    U <- pca$x[,seq(k)]
    V <- pca$v[,seq(k)]
		reconstruction <- U %*% t(V)
    mcv_loss[i] <- mean_squared_error(reconstruction, X2)
    rec_loss[i] <- mean_squared_error(reconstruction, X1)
	}
	plot(k_range, mcv_loss)
  #
  optimal_k <- k_range[which.min(mcv_loss)]
  #
	list(
		optimal_k = optimal_k,
    k_range   = k_range,
    mcv_loss  = mcv_loss,
    rec_loss  = rec_loss
	)
}

#' PCA with Molecular Cross Validation (MCV)
#'
#' @param counts Matrix of genes (rows) and cells (columns).
#' @param ix_genes Integer vector with an index of which genes to include.
#' @param max_pcs Maximum number of PCs to test.
#' @param mcv_reps Number of times to repeat the MCV algorithm. 
#'
do_mcv_pca <- function(
  counts, ix_genes,
  max_pcs = 100,
  mcv_reps = 3
) {
  # pb <- progress::progress_bar$new(
  #   format = "  Running MCV with pca [:bar] :percent :elapsedfull eta: :eta",
  #   total = mcv_reps * (max_pcs - 1)
  # )
  # res <- lapply(seq(mcv_reps), function(mcv_rep) {
  res <- pbmcapply::pbmclapply(
    X = seq(mcv_reps),
    mc.cores = 4,
    FUN = function(mcv_rep) {
    #
    X1 <- counts
    X2 <- counts
    X1@x <- as.numeric(rbinom(n = length(counts@x), size = counts@x, prob = 0.5))
    stopifnot(all(counts@x >= X1@x))
    X2@x <- counts@x - X1@x
    stopifnot(all(X1@x + X2@x == counts@x))
    X1 <- do_log2cpm(X1, total = median(colSums(X1)))
    X2 <- do_log2cpm(X2, total = median(colSums(X2)))
    #
    ix_genes2 <- ix_genes[rowSums(X1[ix_genes,]) > 0]
    X1 <- t(scale_data(X1[ix_genes2,]))
    X2 <- t(scale_data(X2[ix_genes2,]))
    #
    pca <- RSpectra::svds(
      A    = X1,
      k    = max_pcs,
      opts = list(
        center = FALSE,
        scale  = FALSE,
        maxitr = 2000,
        tol    = 1e-10
      )
    )
    pca$x <- pca$u %*% diag(pca$d)
    k_range <- seq(2, max_pcs) #c(seq(2, 10, 1), seq(11, 29, 2), seq(30, max_pcs, 2))
    k_range <- k_range[k_range <= max_pcs]
    mcv_loss <- numeric(length(k_range))
    rec_loss <- numeric(length(k_range))
    mean_squared_error <- function(x, y) mean( (x - y) ^ 2 , na.rm = TRUE)
    for (i in seq_along(k_range)) {
      # pb$tick()
      k <- k_range[i]
      U <- pca$x[,seq(k)]
      V <- pca$v[,seq(k)]
      reconstruction <- U %*% t(V)
      mcv_loss[i] <- mean_squared_error(reconstruction, X2)
      rec_loss[i] <- mean_squared_error(reconstruction, X1)
    }
    optimal_k <- k_range[which.min(mcv_loss)]
    list(
      optimal_k = optimal_k,
      k_range   = k_range,
      mcv_loss  = mcv_loss,
      rec_loss  = rec_loss
    )
  })
  retval <- do.call(rbind, lapply(seq_along(res), function(i) {
    x <- res[[i]]
    data.frame(
      rep = i,
      k = x$k_range,
      mcv_loss = x$mcv_loss,
      rec_loss = x$rec_loss
    )
  }))
  retval$optimal_k <- names(which.min(
    lapply(split(retval$mcv_loss, retval$k), mean)
  ))
  return(retval)
}

#' K-Means with Molecular Cross Validation (MCV)
#'
#' @param counts Matrix of genes (rows) and cells (columns).
#' @param selected Index of which genes to include.
#' @param max_pcs Maximum number of PCs to test.
#'
do_mcv_kmeans <- function(
  counts, selected,
  n_pcs = 7, max_k = 10,
  kmeans_reps = 3, mcv_reps = 3
) {
  pb <- progress::progress_bar$new(
    format = "  Running MCV with kmeans [:bar] :percent :elapsedfull eta: :eta",
    total = mcv_reps * kmeans_reps * (max_k - 1)
  )
  res <- lapply(seq(mcv_reps), function(mcv_rep) {
    #
    X1 <- counts
    X2 <- counts
    X1@x <- as.numeric(rbinom(n = length(counts@x), size = counts@x, prob = 0.5))
    stopifnot(all(counts@x >= X1@x))
    X2@x <- counts@x - X1@x
    stopifnot(all(X1@x + X2@x == counts@x))
    #
    selected2 <- selected[rowSums(X1[selected,]) > 0]
    X1 <- do_log2cpm(X1, total = median(colSums(X1)))
    X2 <- do_log2cpm(X2, total = median(colSums(X2)))
    #
    X1 <- t(scale_data(X1[selected2,]))
    X2 <- t(scale_data(X2[selected2,]))
    #
    pca1 <- RSpectra::svds(
      A    = X1,
      k    = n_pcs,
      opts = list(
        center = FALSE,
        scale  = FALSE,
        maxitr = 2000,
        tol    = 1e-10
      )
    )
    pca1$x <- pca1$u %*% diag(pca1$d)
    #
    k_range <- rep(seq(2, max_k), kmeans_reps)
    mcv_loss <- numeric(length(k_range))
    rec_loss <- numeric(length(k_range))
    #
    for (i in seq_along(k_range)) {
      pb$tick()
      k <- k_range[i]
      km <- ClusterR::KMeans_rcpp(
        data = pca1$x,
        clusters = k,
        seed = i,
        num_init = 3
      )
      centroids <- matrix(ncol = ncol(X1), nrow = k)
      for (j in seq(k)) {
        centroids[j,] <- colMeans(X1[km$clusters == j,,drop=FALSE])
      }
      rec_loss[i] <- mean(unlist(lapply(seq(k), function(j) {
        as.numeric(sweep(
          x      = X1[km$clusters == j,,drop=FALSE],
          MARGIN = 2,
          STATS  = centroids[j,],
          FUN    = "-"
        )) ^ 2
      })))
      mcv_loss[i] <- mean(unlist(lapply(seq(k), function(j) {
        as.numeric(sweep(
          x      = X2[km$clusters == j,,drop=FALSE],
          MARGIN = 2,
          STATS  = centroids[j,],
          FUN    = "-"
        )) ^ 2
      })))
    }
    list(
      k_range   = k_range,
      mcv_loss  = mcv_loss,
      rec_loss  = rec_loss
    )
  })
  retval <- do.call(rbind, lapply(seq_along(res), function(i) {
    x <- res[[i]]
    data.frame(
      rep = i,
      k = x$k_range,
      mcv_loss = x$mcv_loss,
      rec_loss = x$rec_loss
    )
  }))
  retval$optimal_k <- as.numeric(names(which.min(
    lapply(split(retval$mcv_loss, retval$k), mean)
  )))
  return(retval)
}

#' Molecular Cross Validation (MCV) for Leiden
#'
#' @param counts Matrix of genes (rows) and cells (columns).
#' @param selected Index of which genes to include.
#' @param max_pcs Maximum number of PCs to test.
#'
do_mcv_leiden <- function(
  counts, ix_genes,
  n_pcs = 7, n_knn = 30,
  res_range = seq(1, 2, length.out = 5),
  leiden_reps = 3, mcv_reps = 3
) {
  pb <- progress::progress_bar$new(
    format = "  Running MCV with leiden [:bar] :percent :elapsedfull eta: :eta",
    total = mcv_reps * leiden_reps * length(res_range)
  )
  res <- lapply(seq(mcv_reps), function(mcv_rep) {
    #
    X1 <- counts
    X2 <- counts
    X1@x <- as.numeric(rbinom(n = length(counts@x), size = counts@x, prob = 0.5))
    stopifnot(all(counts@x >= X1@x))
    X2@x <- counts@x - X1@x
    stopifnot(all(X1@x + X2@x == counts@x))
    X1 <- do_log2cpm(X1, total = median(colSums(X1)))
    X2 <- do_log2cpm(X2, total = median(colSums(X2)))
    #
    ix_genes <- ix_genes[rowSums(X1[ix_genes,]) > 0]
    X1 <- t(scale_data(X1[ix_genes,]))
    X2 <- t(scale_data(X2[ix_genes,]))
    #
    pca1 <- RSpectra::svds(
      A    = X1,
      k    = n_pcs,
      opts = list(
        center = FALSE,
        scale  = FALSE,
        maxitr = 2000,
        tol    = 1e-10
      )
    )
    pca1$x <- pca1$u %*% diag(pca1$d)
    knn1 <- do_knn(pca1$x, n_knn = n_knn, dist_method = "euclidean")
    #
    k_range <- rep(res_range, leiden_reps)
    mcv_loss <- numeric(length(k_range))
    rec_loss <- numeric(length(k_range))
    #
    for (i in seq_along(k_range)) {
      pb$tick()
      k <- k_range[i]
      clusters <- run_leiden(
        adj = knn1$simil_cells,
        resolution = k,
        iterations = 3,
        seed = i
      )
      centroids <- matrix(ncol = ncol(X1), nrow = max(clusters))
      for (j in seq(max(clusters))) {
        centroids[j,] <- colMeans(X1[clusters == j,,drop=FALSE])
      }
      rec_loss[i] <- mean(unlist(lapply(seq(max(clusters)), function(j) {
        as.numeric(sweep(
          x      = X1[clusters == j,,drop=FALSE],
          MARGIN = 2,
          STATS  = centroids[j,],
          FUN    = "-"
        )) ^ 2
      })))
      mcv_loss[i] <- mean(unlist(lapply(seq(max(clusters)), function(j) {
        as.numeric(sweep(
          x      = X2[clusters == j,,drop=FALSE],
          MARGIN = 2,
          STATS  = centroids[j,],
          FUN    = "-"
        )) ^ 2
      })))
    }
    list(
      k_range   = k_range,
      mcv_loss  = mcv_loss,
      rec_loss  = rec_loss
    )
  })
  retval <- do.call(rbind, lapply(seq_along(res), function(i) {
    x <- res[[i]]
    data.frame(
      rep = i,
      res = x$k_range,
      mcv_loss = x$mcv_loss,
      rec_loss = x$rec_loss
    )
  }))
  retval$optimal_res <- as.numeric(names(which.min(
    lapply(split(retval$mcv_loss, retval$res), mean)
  )))
  return(retval)
}

#' Molecular Cross Validation (MCV) for Leiden (parallel)
#'
#' @param counts Matrix of genes (rows) and cells (columns).
#' @param selected Index of which genes to include.
#' @param max_pcs Maximum number of PCs to test.
#'
do_mcv_leiden2 <- function(
  counts, ix_genes,
  n_pcs = 7, res_range = seq(1, 2, length.out = 5),
  leiden_reps = 3, mcv_reps = 3
) {
  pb <- progress::progress_bar$new(
    format = "  Running MCV with leiden [:bar] :percent :elapsedfull eta: :eta",
    total = mcv_reps
  )
  res <- lapply(
    X = seq(mcv_reps),
    FUN = function(mcv_rep) {
    pb$tick()
    #
    X1 <- counts
    X2 <- counts
    X1@x <- as.numeric(rbinom(n = length(counts@x), size = counts@x, prob = 0.5))
    stopifnot(all(counts@x >= X1@x))
    X2@x <- counts@x - X1@x
    stopifnot(all(X1@x + X2@x == counts@x))
    X1 <- do_log2cpm(X1, total = median(colSums(X1)))
    X2 <- do_log2cpm(X2, total = median(colSums(X2)))
    #
    ix_genes <- ix_genes[rowSums(X1[ix_genes,]) > 0]
    X1 <- t(scale_data(X1[ix_genes,]))
    X2 <- t(scale_data(X2[ix_genes,]))
    #
    pca1 <- RSpectra::svds(
      A    = X1,
      k    = n_pcs,
      opts = list(
        center = FALSE,
        scale  = FALSE,
        maxitr = 2000,
        tol    = 1e-10
      )
    )
    pca1$x <- pca1$u %*% diag(pca1$d)
    knn1 <- do_knn(pca1$x, n_knn = 30, dist_method = "euclidean")
    #
    k_range <- rep(res_range, leiden_reps)
    #
    losses <- mclapply(
      mc.cores = 4,
      X = seq_along(k_range),
      FUN = function(i) {
        k <- k_range[i]
        clusters <- run_leiden(
          adj = knn1$simil_cells,
          resolution = k,
          iterations = 3,
          seed = i
        )
        centroids <- matrix(ncol = ncol(X1), nrow = max(clusters))
        for (j in seq(max(clusters))) {
          centroids[j,] <- colMeans(X1[clusters == j,,drop=FALSE])
        }
        list(
          rec_loss = mean(unlist(lapply(seq(max(clusters)), function(j) {
            as.numeric(sweep(
              x      = X1[clusters == j,,drop=FALSE],
              MARGIN = 2,
              STATS  = centroids[j,],
              FUN    = "-"
            )) ^ 2
          }))),
          mcv_loss = mean(unlist(lapply(seq(max(clusters)), function(j) {
            as.numeric(sweep(
              x      = X2[clusters == j,,drop=FALSE],
              MARGIN = 2,
              STATS  = centroids[j,],
              FUN    = "-"
            )) ^ 2
          })))
        )
      }
    )
    list(
      k_range   = k_range,
      mcv_loss  = sapply(losses, "[[", "mcv_loss"),
      rec_loss  = sapply(losses, "[[", "rec_loss")
    )
  })
  retval <- do.call(rbind, lapply(seq_along(res), function(i) {
    x <- res[[i]]
    data.frame(
      rep = i,
      res = x$k_range,
      mcv_loss = x$mcv_loss,
      rec_loss = x$rec_loss
    )
  }))
  retval$optimal_res <- as.numeric(names(which.min(
    lapply(split(retval$mcv_loss, retval$res), mean)
  )))
  return(retval)
}




##' Main function for running a complete analysis
##' @param obs
##' @param log2cpm
##' @param min_percent
##' @param n_genes Number of robust genes to use for PCA.
##' @param n_pcs
##' @param n_harmony Max number of Harmony iterations.
##' @param n_knn
##' @param leiden_res
##' @param leiden_iter
#do_analysis_old <- function(
#  obs, counts,
#  exclude_genes = NULL,
#  mito_genes    = NULL,
#  loess_span    = 0.05,
#  min_percent   = 0.5,
#  # n_genes       = 2000,
#  n_pcs         = 30,
#  harmony_vars  = c("channel"),
#  n_harmony     = 20,
#  n_knn         = 30,
#  leiden_res    = 1.3,
#  leiden_iter   = 10,
#  umap_spread   = 1,
#  umap_min_dist = 0.01,
#  random_seed   = 42
#) {
#  # browser()
#  if (ncol(counts) != nrow(obs)) {
#    stop("Columns in counts must match rows in obs")
#  }
#  print_status(glue::glue(
#    "Computing cell statistics for {scales::comma(ncol(counts))} cells"
#  ))
#  obs$n_counts    <- Matrix::colSums(counts)
#  obs$n_features  <- Matrix::colSums(counts > 0)
#  ix_mito         <- which(rownames(counts) %in% mito_genes)
#  obs$mito_counts <- colSums(counts[ix_mito,])
#  obs$mito_pct    <- 100 * obs$mito_counts / obs$n_counts
#  print_status(glue::glue(
#    "Computing gene statistics for {scales::comma(nrow(counts))} genes"
#  ))
#  counts_stats <- data.table::data.table(
#    mean    = Matrix::rowMeans(counts),
#    sd      = proxyC::rowSds(counts),
#    percent = 100 * Matrix::rowSums(counts > 0) / counts@Dim[1]
#  )
#  counts_stats$gene <- rownames(counts)
#  counts_stats$exclude <- counts_stats$gene %in% exclude_genes
#  counts_stats$include <- (
#    counts_stats$percent >= min_percent & !counts_stats$exclude
#  )
#  print_status(glue::glue(
#    "{scales::comma(sum(counts_stats$include))} genes have expression in > {signif(min_percent, 3)}% of cells"
#  ))
#  ix_include <- counts_stats$include
#  # Need to be careful that this looks good
#  fit <- loess(
#    formula = log10(sd) ~ log10(mean),
#    data    = counts_stats[ix_include,],
#    span    = loess_span,
#    degree  = 2
#  )
#  # plot(fit$x, fit$residuals)
#  counts_stats$fitted <- NA
#  counts_stats$fitted[ix_include] <- fit$fitted
#  counts_stats$residuals <- NA
#  counts_stats$residuals[ix_include] <- fit$residuals
#  counts_stats$rank <- NA
#  counts_stats$rank[ix_include] <- (
#    rank(rank(-fit$residuals) + rank(-fit$y / fit$fitted))
#  )
#  # This is what I've been using...
#  # ix_genes <- which(with(counts_stats, include & rank <= n_genes))
#  # This probably works better
#  ix_genes <- which(counts_stats$residuals > 0)
#  # plot(fit$x, fit$residuals)
#  # ix_genes <- which(
#  #   counts_stats$residuals > median(counts_stats$residuals, na.rm = TRUE)
#  # )
#  print_status(glue::glue(
#    "{scales::comma(length(ix_genes))} genes have residual variance > 0 for the model log10(sd) ~ log10(mean)"
#  ))
#  # Take the top 80 percent
#  x <- counts_stats$residuals[ix_genes]
#  ix_genes <- ix_genes[x > quantile(x, 0.2)]
#  # Log2CPM
#  log2cpm <- do_log2cpm(counts, total = median(colSums(counts)))
#  if (n_pcs > 0) {
#    print_status(glue::glue(
#      "Running PCA with {scales::comma(length(ix_genes))} genes and {scales::comma(ncol(log2cpm))} cells"
#    ))
#    set.seed(random_seed)
#    pca <- RSpectra::svds(
#      A    = t(log2cpm[ix_genes,]),
#      k    = n_pcs,
#      opts = list(
#        center = TRUE,
#        scale  = TRUE,
#        maxitr = 2000,
#        tol    = 1e-10
#      )
#    )
#    pca$genes <- rownames(log2cpm)[ix_genes]
#    #
#    if (n_harmony > 0) {
#      print_status(glue::glue("Running Harmony with {n_pcs} PCs"))
#      hm <- harmony::HarmonyMatrix(
#        data_mat         = t(pca$u),
#        meta_data        = obs,
#        vars_use         = harmony_vars,
#        max.iter.harmony = n_harmony,
#        do_pca           = FALSE,
#        return_object    = TRUE
#      )
#      pca_h <- as.matrix(t(hm$Z_corr))
#      print_status(glue::glue(
#        "Finding {n_knn} nearest neighbors with harmonized PCs"
#      ))
#    } else {
#      hm <- NULL
#      pca_h <- pca$u
#      print_status(glue::glue(
#        "Finding {n_knn} nearest neighbors with PCs"
#      ))
#    }
#    colnames(pca_h) <- sprintf("PC%s", 1:ncol(pca_h))
#    rownames(pca_h) <- colnames(log2cpm)
#    #
#    # Add PCA coordinates
#    obs <- cbind(obs, pca_h)
#    # KNN
#    knn <- BiocNeighbors::findKNN(
#      X       = pca_h,
#      k       = n_knn,
#      BNPARAM = BiocNeighbors::HnswParam(),
#      BPPARAM = BiocParallel::MulticoreParam(workers = 4)
#    )
#    knn$adj <- Matrix::sparseMatrix(
#      i = rep(1:nrow(knn$index), n_knn),
#      j = as.integer(knn$index)
#    )
#    if (n_harmony > 0) {
#      print_status(glue::glue("Running UMAP with harmonized PCs"))
#    } else {
#      print_status(glue::glue("Running UMAP with PCs"))
#    }
#    # Add UMAP coordinates
#    pca_h_umap <- uwot::umap(
#      X        = pca_h,
#      spread   = umap_spread,
#      min_dist = umap_min_dist,
#      n_threads = 8
#    )
#    # pca_h_umap <- uwot::umap(
#    #   X         = NULL,
#    #   spread    = umap_spread,
#    #   min_dist  = umap_min_dist,
#    #   nn_method = list(
#    #     "idx"  = knn$index,
#    #     "dist" = knn$distance
#    #   )
#    # )
#    obs$UMAP1 <- pca_h_umap[,1]
#    obs$UMAP2 <- pca_h_umap[,2]
#    simil_cells <- proxyC::simil(
#      Matrix::Matrix(pca_h, sparse = TRUE),
#      margin = 1, method = "correlation", rank = n_knn + 1
#    )
#  } else {
#    # No PCA
#    X <- as.matrix(scale(t(log2cpm[ix_genes,])))
#    knn <- BiocNeighbors::findKNN(
#      X       = X,
#      k       = n_knn,
#      BNPARAM = BiocNeighbors::HnswParam(),
#      BPPARAM = BiocParallel::MulticoreParam(workers = 4)
#    )
#    knn$adj <- Matrix::sparseMatrix(
#      i = rep(1:nrow(knn$index), n_knn),
#      j = as.integer(knn$index)
#    )
#    print_status(glue::glue("Running UMAP with {length(ix_genes)} features"))
#    # Add UMAP coordinates
#    log2cpm_umap <- uwot::umap(
#      X        = X,
#      spread   = umap_spread,
#      min_dist = umap_min_dist,
#      n_threads = 8
#    )
#    obs$UMAP1 <- log2cpm_umap[,1]
#    obs$UMAP2 <- log2cpm_umap[,2]
#    simil_cells <- proxyC::simil(
#      X, margin = 1, method = "correlation", rank = n_knn + 1
#    )
#  }
#  # obs$leiden_knn <- run_leiden(
#  #   # adj        = knn$adj,
#  #   adj        = simil_cells,
#  #   resolution = leiden_res,
#  #   iterations = leiden_iter
#  # )
#  # obs$leiden <- obs$leiden_knn
#  print_status(glue::glue(
#    "Running Leiden community detection on KNN with resolution {leiden_res} and {leiden_iter} iterations"
#  ))
#  leiden <- run_leiden(
#    adj        = simil_cells,
#    resolution = leiden_res,
#    iterations = leiden_iter
#  )
#  obs$leiden <- as.integer(1 + unname(t(leiden)[,1])[3:ncol(leiden)])
#  # Build SNN
#  # print_status(glue::glue("Building rank-based SNN from KNN"))
#  # snn                 <- scran:::build_snn_rank(knn$index)
#  # g                   <- igraph::make_graph(edges = snn[[1]])
#  # igraph::E(g)$weight <- snn[[2]]
#  # g                   <- igraph::simplify(g, edge.attr.comb = "first")
#  # snn_adj             <- igraph::as_adjacency_matrix(
#  #   g, attr = "weight", sparse = TRUE
#  # )
#  # print_status(glue::glue(
#  #   "Running Leiden community detection on SNN with resolution {leiden_res} and {leiden_iter} iterations"
#  # ))
#  # obs$leiden_snn <- run_leiden(
#  #   adj        = snn_adj,
#  #   resolution = leiden_res,
#  #   iterations = leiden_iter
#  # )
#  print_status("done")
#  if (n_pcs > 0) {
#    list(
#      obs           = obs,
#      counts        = counts,
#      ix_genes      = ix_genes,
#      exclude_genes = exclude_genes,
#      counts_stats  = counts_stats,
#      fit           = fit,
#      pca           = pca,
#      pca_h         = pca_h,
#      hm            = hm,
#      knn           = knn,
#      simil_cells   = simil_cells,
#      leiden        = leiden,
#      # snn           = snn,
#      de            = presto::wilcoxauc(log2cpm, obs$leiden)
#    )
#  } else {
#    list(
#      obs           = obs,
#      counts        = counts,
#      ix_genes      = ix_genes,
#      exclude_genes = exclude_genes,
#      counts_stats  = counts_stats,
#      fit           = fit,
#      # pca           = pca,
#      # pca_h         = pca_h,
#      # hm            = hm,
#      knn           = knn,
#      simil_cells   = simil_cells,
#      leiden        = leiden,
#      # snn           = snn,
#      de            = presto::wilcoxauc(log2cpm, obs$leiden)
#    )
#  }
#}
