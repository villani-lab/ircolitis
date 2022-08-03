#!/usr/bin/env Rscript
# load-tcr.R
# Kamil Slowikowski

library(data.table)
library(Biostrings)
library(jsonlite)

# @param filename Path to "filtered_contig_annotations.csv"
munge_tcr_contigs <- function(filename) {
  d <- data.table::fread(filename)
  d <- d[d$productive,]
  x <- d[
    chain == "TRA",
    list(
      TRAV        = toString(v_gene),
      TRAD        = toString(d_gene),
      TRAJ        = toString(j_gene),
      TRAC        = toString(c_gene),
      TRA_cdr3    = toString(cdr3),
      TRA_cdr3_nt = toString(cdr3_nt),
      TRA_length  = toString(length),
      TRA_reads   = toString(reads),
      TRA_umis    = toString(umis)
    ),
    by = barcode
  ]
  y <- d[
    chain == "TRB",
    list(
      TRBV        = toString(v_gene),
      TRBD        = toString(d_gene),
      TRBJ        = toString(j_gene),
      TRBC        = toString(c_gene),
      TRB_cdr3    = toString(cdr3),
      TRB_cdr3_nt = toString(cdr3_nt),
      TRB_length  = toString(length),
      TRB_reads   = toString(reads),
      TRB_umis    = toString(umis)
    ),
    by = barcode
  ]
  merge.data.table(x, y, by = 'barcode', all = TRUE)
}

load_tcr <- function(vdj_files) {
  tcr <- data.table::rbindlist(pbapply::pblapply(vdj_files, function(vdj_file) {
    x <- munge_tcr_contigs(vdj_file)
    channel <- vdj_file %>% dirname %>% basename
    x$channel <- channel
    # x$cell <- paste(channel, x$barcode %>% str_remove("-\\d$"), sep = "_")
    x$cell <- paste(channel, x$barcode, sep = "|")
    return(x)
  }))
  # tcr$donor <- str_extract(tcr$channel, "^[^_]+_\\d+")

  # range(sapply(tcr$TRB_cdr3, nchar))
  # #> 5 40

  # Trim the first 4 bases and last 1 base.
  tcr$TRB_cdr3_trim <- substr(tcr$TRB_cdr3, 4, nchar(tcr$TRB_cdr3))
  tcr$TRB_cdr3_trim <- substr(tcr$TRB_cdr3_trim, 1, nchar(tcr$TRB_cdr3_trim) - 1)
  tcr$TRA_cdr3_trim <- substr(tcr$TRA_cdr3, 4, nchar(tcr$TRA_cdr3))
  tcr$TRA_cdr3_trim <- substr(tcr$TRA_cdr3_trim, 1, nchar(tcr$TRA_cdr3_trim) - 1)

  tcr$TRB_cdr3[1:5]
  tcr$TRB_cdr3_trim[1:5]

  # About 2% of TCRs have two beta chains
  # sum(str_detect(tcr$TRB_cdr3, ",")) / nrow(tcr)

  # Thanks to Neal Smith.
  cdr <- jsonlite::read_json(path = "data/human_cdrs_for_10X.json")
  # unique(unlist(lapply(cdr, names)))

  # cdr <- rbindlist(lapply(names(cdr), function(i) data.table(
  #   cdr = i,
  #   gene = names(cdr[[i]]),
  #   sequence = cdr[[i]]
  # )))

  tcr$TRAV2 <- str_replace_all(tcr$TRAV, "/", "")
  tcr$TRBV2 <- str_replace_all(tcr$TRBV, "/", "")

  # Make a new TCR data table without the invalid entries.
  tcr <- tcr[
    !is.na(TRAV) & !is.na(TRBV) &
    !is.na(TRB_cdr3_trim) & !is.na(TRA_cdr3_trim) &
    !str_detect(TRB_cdr3_trim, ",") & !str_detect(TRA_cdr3_trim, ",") &
    !str_detect(TRAV, ",") & !str_detect(TRBV, ",")
  ]
  stopifnot(all(tcr$TRAV2 %in% names(cdr[["cdr2.5"]])))
  stopifnot(all(tcr$TRBV2 %in% names(cdr[["cdr1"]])))
  tcr[["TRAV_cdr1"]]  <- unlist(cdr[["cdr1"]][tcr$TRAV2])
  tcr[["TRAV_cdr2"]]  <- unlist(cdr[["cdr2"]][tcr$TRAV2])
  tcr[["TRAV_cdr25"]] <- unlist(cdr[["cdr2.5"]][tcr$TRAV2])
  tcr[["TRBV_cdr1"]]  <- unlist(cdr[["cdr1"]][tcr$TRBV2])
  tcr[["TRBV_cdr2"]]  <- unlist(cdr[["cdr2"]][tcr$TRBV2])
  tcr[["TRBV_cdr25"]] <- unlist(cdr[["cdr2.5"]][tcr$TRBV2])

  return(tcr)
}

