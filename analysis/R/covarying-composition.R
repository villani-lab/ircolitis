# Copyright 2021 Kamil Slowikowski
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE


#' Summarize genes by cluster and donor.
#'
#' @param log2cpm Sparse matrix of single-cell RNA-seq gene expression (rows are genes).
#' @param cluster Character vector of cluster names.
#' @param donor Character vector of donor names.
#' @param min_donors Genes are kept if they are expressed in at least this many
#' donors, for each cluster.
#' @param min_cells Genes are kept if they are expressed in at least this many
#' cells, for each donor for each cluster.
#' @return A data.table with one row for each gene-pair, for each cluster-pair.
matrix_to_table <- function(
  log2cpm, cluster, donor,
  min_donors = 10,
  min_cells = 25,
  ensembl_to_symbol = NULL
) {
  stopifnot(is(log2cpm,"dgCMatrix"))
  stopifnot(!is.null(ensembl_to_symbol))
  stopifnot(length(cluster) == length(donor))
  stopifnot(length(cluster) == ncol(log2cpm))
  message(glue("log2cpm matrix has {nrow(log2cpm)} genes, {ncol(log2cpm)} samples"))
  obs <- data.table(cluster = cluster, donor = donor)
  # d <- as.data.table(summary(log2cpm[omni_ens,]))
  d <- as.data.table(summary(log2cpm))
  d$symbol  <- ensembl_to_symbol[rownames(log2cpm)[d$i]]
  d$cluster <- obs$cluster[d$j]
  d$donor   <- obs$donor[d$j]
  y <- obs %>% group_by(cluster, donor) %>% count() %>% as.data.table
  d <- y[d, on = c("cluster", "donor")]
  d <- d[
    ,
    .(
      mean = mean(c(x, rep(0, n[1] - length(x)))),
      percent = sum(x > 0) / n[1]
    ),
    by = .(cluster, donor, symbol)
  ]
  setkey(d, symbol, cluster, donor)
  # The gene is robustly detected this many donors for each cluster.
  d_nonzero <- d[,.(nonzero = sum(mean > 0)),by = .(cluster, symbol)]
  robust_genes <- d_nonzero[nonzero >= min_donors]$symbol %>% unique
  message(glue("{length(robust_genes)} robust genes detected in >= {min_donors} donors"))
  d <- d[symbol %in% robust_genes]
  # Add zeros when we don't have anything.
  d <- tidyr::complete(
    d, cluster, donor, symbol, fill = list(mean = 0, percent = 0)
  ) %>% as.data.table
  # Add the number of cells each donor has in each cluster
  d <- left_join(
    x = d,
    y = obs %>% count(donor, cluster),
    by = c("cluster", "donor")
  ) %>% as.data.table
  d <- d[!is.na(d$n)]
  d <- d[d$n >= min_cells]
  return(d)
}

#' Product analysis.
#'
#' @param counts Sparse matrix of single-cell RNA-seq gene expression (rows are genes).
#' @param cluster Character vector of cluster names.
#' @param donor Character vector of donor names.
#' @param min_donors Genes are kept if they are expressed in at least this many
#' donors, for each cluster.
#' @param min_cells Genes are kept if they are expressed in at least this many
#' cells, for each donor for each cluster.
#' @return A data.table with one row for each gene-pair, for each cluster-pair.
gene_gene_products <- function(
  counts, cluster, donor, case,
  source_ens, target_ens,
  min_donors = 10,
  min_cells = 25,
  ensembl_to_symbol = NULL
) {
  stopifnot(is(counts, "dgCMatrix"))
  stopifnot(!is.null(ensembl_to_symbol))
  stopifnot(length(cluster) == length(donor))
  stopifnot(length(case) == length(donor))
  stopifnot(length(cluster) == ncol(counts))
  stopifnot(length(source_ens) == length(target_ens))
  stopifnot(all(source_ens %in% rownames(counts)))
  stopifnot(all(target_ens %in% rownames(counts)))
  all_ens <- unique(c(source_ens, target_ens))
  message(glue("counts matrix has {nrow(counts)} genes, {ncol(counts)} samples"))
  my_obs <- data.table(cluster = cluster, donor = donor, case = case)
  # Pseudobulk at the cluster level
  ########################################################################
  message("creating pseudobulk log2cpm matrix")
  y <- with(my_obs, model.matrix(~ 0 + factor(cluster):factor(donor)))
  y <- as(y, "dgCMatrix")
  pb <- as(counts %*% y, "dgCMatrix")
  pb <- do_log2cpm(pb, median(Matrix::colSums(pb)))
  colnames(pb) <- str_replace(colnames(pb), "factor\\(cluster\\)", "")
  colnames(pb) <- str_replace(colnames(pb), "factor\\(donor\\)", "")
  pb <- pb[all_ens,]
  my_counts <- counts[all_ens,]
  #
  pb_meta <- str_split_fixed(colnames(pb), ":", 2)
  colnames(pb_meta) <- c("cluster", "donor")
  pb_meta <- as_tibble(pb_meta)
  pb_meta <- left_join(pb_meta, unique(my_obs[,c("donor","case")]), by = "donor")
  pb_meta$case <- factor(pb_meta$case, c("Control", "Case"))
  stopifnot(nrow(pb_meta) == ncol(pb))
  d <- as.data.table(summary(pb))
  d$cluster <- pb_meta$cluster[d$j]
  d$donor <- pb_meta$donor[d$j]
  d$ens <- rownames(pb)[d$i]
  pb <- as.matrix(pb)
  my_pairs <- data.frame(source_ens = source_ens, target_ens = target_ens)
  message(glue("testing {nrow(my_pairs)} gene pairs"))
  #
  # Test one gene pair at a time
  res <- pblapply(seq(nrow(my_pairs)), function(omni_i) {
    e1 <- my_pairs$source_ens[omni_i]
    e2 <- my_pairs$target_ens[omni_i]
    g1 <- ensembl_to_symbol[e1]
    g2 <- ensembl_to_symbol[e2]
    #
    x <- cbind.data.frame(pb_meta, value = pb[e1,])
    y <- cbind.data.frame(pb_meta, value = pb[e2,])
    stopifnot(nrow(x) == nrow(y))
    xx <- data.table::dcast(as.data.table(x), donor ~ cluster, value.var = "value")
    xx_donor <- xx$donor
    xx <- as.matrix(xx[,2:ncol(xx)])
    rownames(xx) <- xx_donor
    yy <- data.table::dcast(as.data.table(y), donor ~ cluster, value.var = "value")
    yy_donor <- yy$donor
    yy <- as.matrix(yy[,2:ncol(yy)])
    rownames(yy) <- yy_donor
    # Remove NAs
    xx[is.na(xx)] <- 0
    yy[is.na(yy)] <- 0
    #
    # Add logCPM of gene1 to logCPM of gene2
    # res <- t(do.call(cbind, lapply(1:ncol(xx), function(i) xx[,i] + yy)))
    res <- t(do.call(cbind, lapply(1:ncol(xx), function(i) (xx[,i] + yy) / 2)))
    #
    # Dataframe with one row for each pair of cell clusters
    d_res <- expand.grid(c2 = colnames(xx), c1 = colnames(xx))
    d_res$g1 <- g1
    d_res$g2 <- g2
    #
    res_meta <- left_join(
      x = data.frame(donor = colnames(res)),
      y = unique(pb_meta[,c("donor","case")]),
      by = "donor"
    )
    des1 <- with(res_meta, model.matrix(~ case))
    # 
    fit1 <- limma::lmFit(object = res, design = des1)
    fit1 <- limma::eBayes(fit1)
    limma_res <- limma::topTable(
      fit     = fit1,
      coef    = "caseCase",
      sort.by = "none",
      number  = nrow(res),
      confint = TRUE
    )
    limma_res <- as.data.table(cbind.data.frame(d_res, limma_res))
    #
    return(list("limma" = limma_res, "product" = res))
  })
  #
  products <- do.call(rbind, lapply(res, "[[", "product"))
  limma_res <- rbindlist(lapply(res, "[[", "limma")) %>% select(-ID)
  # Drop duplicates
  limma_res <- limma_res %>%
    rowwise() %>%
    mutate(c_rank = diff(rank(c(c1, c2)))) %>%
    mutate(
      id = ifelse(c_rank == 1, glue("{c1} {g1} {c2} {g2}"), glue("{c2} {g2} {c1} {g1}"))
    ) %>%
    select(-c_rank) %>%
    ungroup()
  ix_dup <- duplicated(limma_res$id)
  products <- products[!ix_dup,]
  limma_res <- limma_res[!ix_dup,]
  rownames(products) <- limma_res$id
  limma_res <- as.data.table(limma_res)
  #
  d_n <- d[, .(donors = sum(x > 0)), by = .(cluster, ens)]
  d_n$cluster <- as.character(d_n$cluster)
  d_n$symbol <- ensembl_to_symbol[d_n$ens]
  #
  res2 <- d_n[limma_res, on = c("symbol" = "g1", "cluster" = "c1")]
  res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "g2")] 
  colnames(res3) <- c(
    "c1", "e1", "donors1", "g1", "c2", "e2", "donors2", "g2", "logFC",
    "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "id"
  )
  limma_res <- res3
  rm(res2)
  rm(res3)
  #
  limma_res %<>% arrange(P.Value)
  limma_res <- limma_res[donors1 >= min_donors & donors2 >= min_donors,]
  # Get statistics for each cluster, regardless of donors.
  d2 <- as.data.table(summary(my_counts))
  d2$ens <- rownames(my_counts)[d2$i]
  d2$cluster <- my_obs$cluster[d2$j]
  y <- my_obs %>% group_by(cluster) %>% count() %>% as.data.table
  d2 <- y[d2, on = c("cluster")]
  d2 <- d2[
    ,
    .(
      percent = sum(x > 0) / n[1]
    ),
    by = .(cluster, ens)
  ]
  limma_res <- left_join(
    x = limma_res,
    y = d2 %>% select(c1 = cluster, e1 = ens, p1 = percent),
    by = c("c1", "e1")
  )
  limma_res <- left_join(
    x = limma_res,
    y = d2 %>% select(c2 = cluster, e2 = ens, p2 = percent),
    by = c("c2", "e2")
  )
  limma_res %<>%
    group_by(c1, c2) %>%
    mutate(fdr = p.adjust(P.Value, method = "BH")) %>%
    ungroup()
  return(list(
    pb = pb,
    pb_meta = pb_meta,
    # products = products,
    limma = limma_res
  ))
}

#' Product analysis.
#'
#' @param log2cpm Sparse matrix of single-cell RNA-seq gene expression (rows are genes).
#' @param cluster Character vector of cluster names.
#' @param donor Character vector of donor names.
#' @param min_donors Genes are kept if they are expressed in at least this many
#' donors, for each cluster.
#' @param min_cells Genes are kept if they are expressed in at least this many
#' cells, for each donor for each cluster.
#' @return A data.table with one row for each gene-pair, for each cluster-pair.
gene_gene_products2 <- function(
  log2cpm, cluster, donor,
  min_donors = 10,
  min_cells = 25,
  source_genes, target_genes,
  ensembl_to_symbol = NULL
) {
  stopifnot(is(log2cpm,"dgCMatrix"))
  stopifnot(!is.null(ensembl_to_symbol))
  stopifnot(length(cluster) == length(donor))
  stopifnot(length(cluster) == ncol(log2cpm))
  stopifnot(length(source_genes) == length(target_genes))
  message(glue("log2cpm matrix has {nrow(log2cpm)} genes, {ncol(log2cpm)} samples"))
  obs <- data.table(cluster = cluster, donor = donor)
  # d <- as.data.table(summary(log2cpm[omni_ens,]))
  d <- as.data.table(summary(log2cpm))
  d$symbol  <- ensembl_to_symbol[rownames(log2cpm)[d$i]]
  d$cluster <- obs$cluster[d$j]
  d$donor   <- obs$donor[d$j]
  y <- obs %>% group_by(cluster, donor) %>% count() %>% as.data.table
  d <- y[d, on = c("cluster", "donor")]
  d <- d[
    ,
    .(
      mean = mean(c(x, rep(0, n[1] - length(x)))),
      percent = sum(x > 0) / n[1]
    ),
    by = .(cluster, donor, symbol)
  ]
  setkey(d, symbol, cluster, donor)
  # The gene is robustly detected this many donors for each cluster.
  d_nonzero <- d[,.(nonzero = sum(mean > 0)),by = .(cluster, symbol)]
  robust_genes <- d_nonzero[nonzero >= min_donors]$symbol %>% unique
  message(glue("{length(robust_genes)} robust genes detected in >= {min_donors} donors"))
  d <- d[symbol %in% robust_genes]
  # Add zeros when we don't have anything.
  d <- tidyr::complete(
    d, cluster, donor, symbol, fill = list(mean = 0, percent = 0)
  ) %>% as.data.table
  ix <- source_genes %in% robust_genes & target_genes %in% robust_genes
  omni <- data.table(
    source_genesymbol = source_genes[ix],
    target_genesymbol = target_genes[ix]
  )
  # Add the number of cells each donor has in each cluster
  d <- left_join(
    x = d,
    y = obs %>% count(donor, cluster),
    by = c("cluster", "donor")
  ) %>% as.data.table
  d <- d[!is.na(d$n)]
  d <- d[d$n >= min_cells]
  # Compute gene-gene products
  res <- pblapply(seq(nrow(my_omni)), function(omni_i) {
    # g1 <- "ADAM15"
    # g2 <- "ITGAV"
    g1 <- my_omni$source_genesymbol[omni_i]
    g2 <- my_omni$target_genesymbol[omni_i]
    #
    x <- d[symbol == g1]
    y <- d[symbol == g2]
    stopifnot(nrow(x) == nrow(y))
    xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
    # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
    xx_donor <- xx$donor
    xx <- as.matrix(xx[,2:ncol(xx)])
    rownames(xx) <- xx_donor
    yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
    # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
    yy_donor <- yy$donor
    yy <- as.matrix(yy[,2:ncol(yy)])
    rownames(yy) <- yy_donor
    # Remove NAs
    xx[is.na(xx)] <- 0
    yy[is.na(yy)] <- 0
    xx <- log(xx)
    yy <- log(yy)
    xx[!is.finite(xx)] <- 0
    yy[!is.finite(yy)] <- 0
    #
    # res <- t(do.call(cbind, lapply(1:ncol(xx), function(i) xx[,i] * yy)))
    res <- t(do.call(cbind, lapply(1:ncol(xx), function(i) xx[,i] + yy)))
    # res <- log(res)
    # res[!is.finite(res)] <- NA
    #
    d_res <- expand.grid(c2 = colnames(xx), c1 = colnames(xx))
    d_res$g1 <- g1
    d_res$g2 <- g2
    #
    res_meta <- left_join(
      x = data.frame(donor = colnames(res)),
      y = unique(pb_meta[,c("donor","case")]),
      by = "donor"
    )
    des1 <- with(res_meta, model.matrix(~ case))
    # 
    fit1 <- limma::lmFit(object = res, design = des1)
    fit1 <- limma::eBayes(fit1)
    limma_res <- limma::topTable(fit1, coef = "caseCase", sort.by = "none", number = 1e3, confint = TRUE)
    limma_res <- cbind.data.frame(d_res, limma_res)
    #
    return(list("limma" = limma_res, "product" = res))
  })
  #
  products <- do.call(rbind, lapply(res, "[[", "product"))
  limma_res <- rbindlist(lapply(res, "[[", "limma")) %>% select(-ID)
  limma_res$id <- with(limma_res, glue("{c1}_{c2}_{g1}_{g2}"))
  rownames(products) <- limma_res$id
  #
  d_n <- d[, .(donors = sum(mean > 0)), by = .(cluster, symbol)]
  d_n$cluster <- as.character(d_n$cluster)
  #
  res2 <- d_n[limma_res, on = c("symbol" = "g1", "cluster" = "c1")]
  res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "g2")] 
  colnames(res3) <- c("c1", "g1", "donors1", "c2", "g2", "donors2", "logFC", "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "id")
  limma_res <- res3
  rm(res2)
  rm(res3)
  #
  limma_res %<>% arrange(P.Value)
  limma_res <- limma_res[donors1 >= min_donors & donors2 >= min_donors,]
  # Get statistics for each cluster, regardless of donors.
  d2 <- as.data.table(summary(log2cpm))
  d2$symbol  <- ensembl_to_symbol[rownames(log2cpm)[d2$i]]
  d2$cluster <- obs$cluster[d2$j]
  y <- obs %>% group_by(cluster) %>% count() %>% as.data.table
  d2 <- y[d2, on = c("cluster")]
  d2 <- d2[
    ,
    .(
      mean = mean(c(x, rep(0, n[1] - length(x)))),
      percent = sum(x > 0) / n[1]
    ),
    by = .(cluster, symbol)
  ]
  limma_res <- left_join(
    x = limma_res,
    y = d2 %>% select(c1 = cluster, g1 = symbol, p1 = percent, m1 = mean),
    by = c("c1", "g1")
  )
  limma_res <- left_join(
    x = limma_res,
    y = d2 %>% select(c2 = cluster, g2 = symbol, p2 = percent, m2 = mean),
    by = c("c2", "g2")
  )
  limma_res %<>% group_by(c1, c2) %>% mutate(fdr = p.adjust(P.Value, method = "BH"))
  return(list(products = products, limma = limma_res))
}

#' Covarying composition analysis.
#'
#' @param log2cpm Sparse matrix of single-cell RNA-seq gene expression (rows are genes).
#' @param cluster Character vector of cluster names.
#' @param donor Character vector of donor names.
#' @param min_donors Genes are kept if they are expressed in at least this many
#' donors, for each cluster.
#' @param min_cells Genes are kept if they are expressed in at least this many
#' cells, for each donor for each cluster.
#' @return A data.table with one row for each gene-pair, for each cluster-pair.
covarying_composition <- function(
  log2cpm, cluster, donor,
  min_donors = 10,
  min_cells = 25,
  source_genes, target_genes,
  ensembl_to_symbol = NULL
) {
  stopifnot(is(log2cpm,"dgCMatrix"))
  stopifnot(!is.null(ensembl_to_symbol))
  stopifnot(length(cluster) == length(donor))
  stopifnot(length(cluster) == ncol(log2cpm))
  stopifnot(length(source_genes) == length(target_genes))
  message(glue("log2cpm matrix has {nrow(log2cpm)} genes, {ncol(log2cpm)} samples"))
  obs <- data.table(cluster = cluster, donor = donor)
  # d <- as.data.table(summary(log2cpm[omni_ens,]))
  d <- as.data.table(summary(log2cpm))
  d$symbol  <- ensembl_to_symbol[rownames(log2cpm)[d$i]]
  d$cluster <- obs$cluster[d$j]
  d$donor   <- obs$donor[d$j]
  y <- obs %>% group_by(cluster, donor) %>% count() %>% as.data.table
  d <- y[d, on = c("cluster", "donor")]
  d <- d[
    ,
    .(
      mean = mean(c(x, rep(0, n[1] - length(x)))),
      percent = sum(x > 0) / n[1]
    ),
    by = .(cluster, donor, symbol)
  ]
  setkey(d, symbol, cluster, donor)
  # The gene is robustly detected this many donors for each cluster.
  d_nonzero <- d[,.(nonzero = sum(mean > 0)),by = .(cluster, symbol)]
  robust_genes <- d_nonzero[nonzero >= min_donors]$symbol %>% unique
  message(glue("{length(robust_genes)} robust genes detected in >= {min_donors} donors"))
  d <- d[symbol %in% robust_genes]
  # Add zeros when we don't have anything.
  d <- tidyr::complete(
    d, cluster, donor, symbol, fill = list(mean = 0, percent = 0)
  ) %>% as.data.table
  ix <- source_genes %in% robust_genes & target_genes %in% robust_genes
  omni <- data.table(
    source_genesymbol = source_genes[ix],
    target_genesymbol = target_genes[ix]
  )
  # Add the number of cells each donor has in each cluster
  d <- left_join(
    x = d,
    y = obs %>% count(donor, cluster),
    by = c("cluster", "donor")
  ) %>% as.data.table
  d <- d[!is.na(d$n)]
  d <- d[d$n >= min_cells]
  # Compute spearman correlations
  res <- pblapply(seq(nrow(omni)), function(omni_i) {
    # g1 <- "ADAM15"
    # g2 <- "ITGAV"
    g1 <- omni$source_genesymbol[omni_i]
    g2 <- omni$target_genesymbol[omni_i]
    x <- d[symbol == g1]
    y <- d[symbol == g2]
    stopifnot(nrow(x) == nrow(y))
    xx <- data.table::dcast(x, donor ~ cluster, value.var = "percent")
    # xx <- data.table::dcast(x, donor ~ cluster, value.var = "mean")
    xx_donor <- xx$donor
    xx <- as.matrix(xx[,2:ncol(xx)])
    rownames(xx) <- xx_donor
    yy <- data.table::dcast(y, donor ~ cluster, value.var = "percent")
    # yy <- data.table::dcast(y, donor ~ cluster, value.var = "mean")
    yy_donor <- yy$donor
    yy <- as.matrix(yy[,2:ncol(yy)])
    rownames(yy) <- yy_donor
    # Remove zeros
    xx[xx == 0] <- NA
    yy[yy == 0] <- NA
    sp <- Hmisc::rcorr(xx, yy, "spearman")
    sp_r <- sp$r[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
    sp_p <- sp$P[seq(ncol(xx)),ncol(xx) + seq(ncol(xx))]
    ut <- upper.tri(sp_r, diag = FALSE)
    lt <- lower.tri(sp_r, diag = FALSE)
    res <- rbindlist(list(
      data.table(
        g1 = g1,
        g2 = g2,
        c1 = rownames(sp_r)[row(sp_r)[ut]],
        c2 = rownames(sp_r)[col(sp_r)[ut]],
        estimate = sp_r[ut],
        p.value = sp_p[ut]
      ),
      data.table(
        g1 = g1,
        g2 = g2,
        c1 = rownames(sp_r)[row(sp_r)[lt]],
        c2 = rownames(sp_r)[col(sp_r)[lt]],
        estimate = sp_r[lt],
        p.value = sp_p[lt]
      )
    ))
  })
  res <- rbindlist(res)
  res <- res[!is.na(res$estimate)]
  res$p.value[res$p.value == 0] <- 1
  d_n <- d[, .(donors = sum(mean > 0)), by = .(cluster, symbol)]
  d_n$cluster <- as.character(d_n$cluster)
  res2 <- d_n[res, on = c("symbol" = "g1", "cluster" = "c1")]
  res3 <- d_n[res2, on = c("cluster" = "c2", "symbol" = "g2")] 
  colnames(res3) <- c("c1", "g1", "donors1", "c2", "g2", "donors2", "estimate", "p.value")
  res <- res3
  rm(res2)
  rm(res3)
  res %<>% arrange(p.value)
  res <- unique(res)
  res <- res[donors1 >= min_donors & donors2 >= min_donors,]
  # Get statistics for each cluster, regardless of donors.
  d2 <- as.data.table(summary(log2cpm))
  d2$symbol  <- ensembl_to_symbol[rownames(log2cpm)[d2$i]]
  d2$cluster <- obs$cluster[d2$j]
  y <- obs %>% group_by(cluster) %>% count() %>% as.data.table
  d2 <- y[d2, on = c("cluster")]
  d2 <- d2[
    ,
    .(
      mean = mean(c(x, rep(0, n[1] - length(x)))),
      percent = sum(x > 0) / n[1]
    ),
    by = .(cluster, symbol)
  ]
  res <- left_join(
    x = res,
    y = d2 %>% select(c1 = cluster, g1 = symbol, p1 = percent, m1 = mean),
    by = c("c1", "g1")
  )
  res <- left_join(
    x = res,
    y = d2 %>% select(c2 = cluster, g2 = symbol, p2 = percent, m2 = mean),
    by = c("c2", "g2")
  )
  res %<>% group_by(c1, c2) %>% mutate(fdr = p.adjust(p.value, method = "BH"))
  return(list(percents = d, correlations = res))
}

#' Covarying composition analysis.
#'
#' @param log2cpm Sparse matrix of single-cell RNA-seq gene expression (rows are genes).
#' @param cluster Character vector of cluster names.
plot_covarying_genes <- function(
  d, g1, g2, c1, c2, value.var = "percent", out_dir = "results/a20/communication/cc",
  title = TRUE,
  scale = 0.8
) {
  x <- d[cluster == c1][symbol == g1]
  y <- d[cluster == c2][symbol == g2]
  both_donors <- intersect(x$donor, y$donor)
  y <- y[donor %in% both_donors]
  x <- x[donor %in% both_donors]
  stopifnot(nrow(x) == nrow(y))
  stopifnot(all(x$donor == y$donor))
  x[[value.var]][x[[value.var]] == 0] <- NA
  y[[value.var]][y[[value.var]] == 0] <- NA
  x_t <- list(estimate = 0, pvalue_label = 1)
  try({
    x_t <- broom::tidy(
      suppressWarnings(cor.test(x[[value.var]], y[[value.var]], method = "spearman", exact = FALSE))
    )
    x_t$pvalue_label <- signif(x_t$p.value, 2) %>% str_replace("e-0", "e-")
  })
  x$g1_pct <- x[[value.var]]
  x$g2_pct <- y[[value.var]]
  x$g1_pct[is.na(x$g1_pct)] <- 0
  x$g2_pct[is.na(x$g2_pct)] <- 0
  #
  x_case <- unique(x[!is.na(x$g1_pct) & !is.na(x$g2_pct),c("donor","case")])
  xc <- table(x_case$case)
  xn <- sprintf("%s (n = %s)", names(xc), xc)
  case_to_n <- xn
  names(case_to_n) <- names(xc)
  x$case_n <- case_to_n[x$case]
  #
  c_name <- function(x) {
    retval <- x
    if (str_detect(x, "^T")) {
      retval <- glue("CD4T-{str_remove(x, '^T')}")
    }
    if (str_detect(x, "^CT")) {
      retval <- glue("CD8T-{str_remove(x, '^CT')}")
    }
    if (str_detect(x, "^E")) {
      retval <- glue("E-{str_remove(x, '^E')}")
    }
    retval
  }
  p <- ggplot(x %>% as.data.frame) + 
    geom_point(
      mapping = aes(x = g1_pct, y = g2_pct, fill = case_n),
      size = 4, shape = 21, stroke = 0.3
    ) +
    labs(
      x = glue("<i>{g1}</i> {c_name(c1)}"),
      y = glue("<i>{g2}</i> {c_name(c2)}")
    ) +
    scale_fill_manual(name = NULL, values = pals::okabe()[2:1]) +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    scale_y_continuous(breaks = pretty_breaks(4)) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    theme_kamil +
    theme(
      # axis.title.x = element_text(face = "italic"),
      # axis.title.y = element_text(face = "italic")
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown()
    )
  if (title) {
    p <- p + labs(
      subtitle = glue("rho = {signif(x_t$estimate, 2)}, P = {x_t$pvalue_label}")
    )
  }
  my_filename <- glue("{out_dir}/cc_{value.var}_{g1}_{g2}_{c1}_{c2}.pdf")
  message(my_filename)
  ggsave(
    filename = my_filename,
    plot = p,
    scale = scale, width = 6, height = 3, units = "in", dpi = 300
  )
}

plot_communication <- function(
  percents, correlations,
  log_percent = FALSE,
  out_dir = "results/a20/communication/cc"
) {
  for (i in seq(nrow(correlations))) {
    g1 <- correlations$g1[i]
    g2 <- correlations$g2[i]
    c1 <- correlations$c1[i]
    c2 <- correlations$c2[i]
    x <- percents[cluster == c1][symbol == g1]# %>% filter(cluster == c1, symbol == g1)
    y <- percents[cluster == c2][symbol == g2]#%>% filter(cluster == c2, symbol == g2)
    both_donors <- intersect(x$donor, y$donor)
    y <- y[donor %in% both_donors]
    x <- x[donor %in% both_donors]
    stopifnot(nrow(x) == nrow(y))
    stopifnot(all(x$donor == y$donor))
    x$percent[x$percent == 0] <- NA
    y$percent[y$percent == 0] <- NA
    x_t <- broom::tidy(
      suppressWarnings(cor.test(x$percent, y$percent, method = "spearman", exact = FALSE))
    )
    x_t$pvalue_label <- signif(x_t$p.value, 2) %>% str_replace("e-0", "e-")
    x$g1_pct <- x$percent
    x$g2_pct <- y$percent
    x$g1_pct[is.na(x$g1_pct)] <- 0
    x$g2_pct[is.na(x$g2_pct)] <- 0
    #
    x_case <- unique(x[!is.na(x$g1_pct) & !is.na(x$g2_pct),c("donor","case")])
    xc <- table(x_case$case)
    xn <- sprintf("%s (n = %s)", names(xc), xc)
    case_to_n <- xn
    names(case_to_n) <- names(xc)
    x$case_n <- case_to_n[x$case]
    #
    p <- ggplot(x %>% as.data.frame) + 
      geom_point(
        mapping = aes(x = g1_pct, y = g2_pct, color = case_n),
        size = 3
      ) +
      labs(
        # x = g1, y = g2, title = glue("{c1} and {c2}"),
        x = glue("<i>{g1}</i>+ {c1}"),
        y = glue("<i>{g2}</i>+ {c2}"),
        subtitle = glue("rho = {signif(x_t$estimate, 2)}, P = {x_t$pvalue_label}")
      ) +
      # scale_x_continuous(labels = function(x) signif(100 * x, 2)) +
      scale_color_manual(name = NULL, values = pals::okabe()[2:1]) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      theme_kamil +
      theme(
        # axis.title.x = element_text(face = "italic"),
        # axis.title.y = element_text(face = "italic")
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        axis.ticks = element_line(size = 0.3)
      )
    if (log_percent) {
      p <- p +
      scale_x_continuous(
        labels = function(x) signif(100 * x, 2), trans = "log10",
        breaks = log_breaks()
      ) +
      scale_y_continuous(
        labels = function(x) signif(100 * x, 2), trans = "log10",
        breaks = log_breaks()
      ) +
      annotation_logticks(sides = "bl", size = 0.3)
    } else {
      p <- p +
      scale_x_continuous(labels = function(x) signif(100 * x, 2)) +
      scale_y_continuous(labels = function(x) signif(100 * x, 2))
    }
    my_ggsave(
      glue("cc_{c1}_{c2}_{g1}_{g2}"),
      out_dir = out_dir,
      types = "pdf",
      plot = p,
      scale = 1, width = 5.5, height = 3, units = "in", dpi = 300
    )
  }
}

plot_communication2 <- function(
  percents, correlations,
  log_percent = TRUE,
  out_dir = "results/a20/communication/cc"
) {
  for (i in seq(nrow(correlations))) {
    g1 <- correlations$g1[i]
    g2 <- correlations$g2[i]
    c1 <- correlations$c1[i]
    c2 <- correlations$c2[i]
    x <- percents[cluster == c1][symbol == g1]# %>% filter(cluster == c1, symbol == g1)
    y <- percents[cluster == c2][symbol == g2]#%>% filter(cluster == c2, symbol == g2)
    both_donors <- intersect(x$donor, y$donor)
    y <- y[donor %in% both_donors]
    x <- x[donor %in% both_donors]
    stopifnot(nrow(x) == nrow(y))
    stopifnot(all(x$donor == y$donor))
    x$percent[x$percent == 0] <- NA
    y$percent[y$percent == 0] <- NA
    x_lm <- data.frame(x = x$percent, y = y$percent, case = x$case)
    x_t <- rbindlist(lapply(c("Case", "Control"), function(this_case) {
      retval <- broom::tidy(
        suppressWarnings(cor.test(
          x$percent[x$case == this_case],
          y$percent[y$case == this_case],
          method = "spearman",
          exact = FALSE
        ))
      )
      retval$case <- this_case
      retval
    }))
    x_t$pvalue_label <- signif(x_t$p.value, 1) %>% str_replace("e-0", "e-")
    x_t$label <-  with(x_t, glue(
      "rho = {signif(estimate, 2)}, P = {pvalue_label}"
    ))
    x$g1_pct <- x$percent
    x$g2_pct <- y$percent
    x$g1_pct[is.na(x$g1_pct)] <- 0
    x$g2_pct[is.na(x$g2_pct)] <- 0
    #
    x_case <- unique(x[!is.na(x$g1_pct) & !is.na(x$g2_pct),c("donor","case")])
    xc <- table(x_case$case)
    xn <- sprintf("%s (n = %s)", names(xc), xc)
    case_to_n <- xn
    names(case_to_n) <- names(xc)
    x$case_n <- case_to_n[x$case]
    x_t$case_n <- case_to_n[x_t$case]
    #
    p <- ggplot() +
      geom_point(
        data = x %>% as.data.frame,
        # mapping = aes(x = g1_pct, y = g2_pct, color = case_n, size = g1_pct == 0 | g2_pct == 0)
        mapping = aes(x = g1_pct, y = g2_pct, color = case_n),
        size = 3
      ) +
      # scale_size_manual(values = c(3, 0.66), guide = "none") +
      geom_label(
        data    = x_t,
        mapping = aes(x = 0, y = Inf, label = label),
        hjust   = -0.07,
        vjust   = 1.3,
        size    = 5,
        label.padding = unit(0.1, "lines"),
        alpha = 0.5,
        label.size = NA
      ) +
      scale_color_manual(name = NULL, values = pals::okabe()[2:1]) +
      guides(
        # color = guide_legend(override.aes = list(size = 5))
        color = "none"
      ) +
      labs(
        # x = g1, y = g2, title = glue("{c1} and {c2}")
        x = glue("<i>{g1}</i>+ {c1}"),
        y = glue("<i>{g2}</i>+ {c2}")
      ) +
      theme_kamil +
      theme(
        panel.spacing = unit(1, "lines"),
        # axis.title.x = element_text(face = "italic"),
        # axis.title.y = element_text(face = "italic")
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        axis.ticks = element_line(size = 0.3)
      ) +
      facet_wrap(~ case_n)
    if (log_percent) {
      p <- p +
      scale_x_continuous(
        labels = function(x) signif(100 * x, 2), trans = "log10",
        breaks = log_breaks()
      ) +
      scale_y_continuous(
        labels = function(x) signif(100 * x, 2), trans = "log10",
        breaks = log_breaks()
      ) +
      annotation_logticks(sides = "bl", size = 0.3)
    } else {
      p <- p +
      scale_x_continuous(labels = function(x) signif(100 * x, 2)) +
      scale_y_continuous(labels = function(x) signif(100 * x, 2))
    }
    my_ggsave(
      glue("cc_{g1}_{g2}_{c1}_{c2}"),
      out_dir = out_dir,
      types = "pdf",
      plot = p,
      scale = 1, width = 6, height = 3, units = "in", dpi = 300
    )
  }
}

