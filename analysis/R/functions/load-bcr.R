#!/usr/bin/env Rscript
# load-bcr.R
# Kamil Slowikowski

pacman::p_load(
  data.table,
  Biostrings,
  jsonlite
)

# @param filename Path to "filtered_contig_annotations.csv"
munge_bcr_contigs <- function(filename) {
  d <- data.table::fread(filename)
  d <- d[d$productive,]
  x <- d[
    chain == "IGH",
    list(
      IGHV        = toString(v_gene),
      IGHD        = toString(d_gene),
      IGHJ        = toString(j_gene),
      IGHC        = toString(c_gene),
      IGH_cdr3    = toString(cdr3),
      IGH_cdr3_nt = toString(cdr3_nt),
      IGH_length  = toString(length),
      IGH_reads   = toString(reads),
      IGH_umis    = toString(umis)
    ),
    by = barcode
  ]
  y <- d[
    chain == "IGK",
    list(
      IGKV        = toString(v_gene),
      IGKD        = toString(d_gene),
      IGKJ        = toString(j_gene),
      IGKC        = toString(c_gene),
      IGK_cdr3    = toString(cdr3),
      IGK_cdr3_nt = toString(cdr3_nt),
      IGK_length  = toString(length),
      IGK_reads   = toString(reads),
      IGK_umis    = toString(umis)
    ),
    by = barcode
  ]
  res <- merge.data.table(x, y, by = 'barcode', all = TRUE)
  rm(x)
  rm(y)
  y <- d[
    chain == "IGL",
    list(
      IGLV        = toString(v_gene),
      IGLD        = toString(d_gene),
      IGLJ        = toString(j_gene),
      IGLC        = toString(c_gene),
      IGL_cdr3    = toString(cdr3),
      IGL_cdr3_nt = toString(cdr3_nt),
      IGL_length  = toString(length),
      IGL_reads   = toString(reads),
      IGL_umis    = toString(umis)
    ),
    by = barcode
  ]
  res <- merge.data.table(res, y, by = 'barcode', all = TRUE)
  return(res)
}

load_bcr <- function(vdj_files) {
  bcr <- data.table::rbindlist(pbapply::pblapply(vdj_files, function(vdj_file) {
    x <- munge_bcr_contigs(vdj_file)
    channel <- vdj_file %>% dirname %>% basename
    x$channel <- channel
    # x$cell <- paste(channel, x$barcode %>% str_remove("-\\d$"), sep = "_")
    x$cell <- paste(channel, x$barcode, sep = "|")
    return(x)
  }))
  # bcr$donor <- str_extract(bcr$channel, "^[^_]+_\\d+")
  # Multiple chains... probably not good?
  sum(str_detect(bcr$IGH_cdr3, ","), na.rm = TRUE) / nrow(bcr)
  sum(str_detect(bcr$IGK_cdr3, ","), na.rm = TRUE) / nrow(bcr)
  sum(str_detect(bcr$IGL_cdr3, ","), na.rm = TRUE) / nrow(bcr)
  # Make a new TCR data table without the invalid entries.
  # Eventually get the sequences from here:
  # http://www.imgt.org/vquest/refseqh.html#:~:text=The%20IMGT%2FV%2DQUEST%20reference,one%20sequence%20for%20each%20allele
  # bcr <- bcr[
  #   !str_detect(IGH_cdr3, ",") &
  #   !str_detect(IGK_cdr3, ",") &
  #   !str_detect(IGL_cdr3, ",")
  # ]
  return(bcr)
}

