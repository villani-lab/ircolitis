#' Write differential expression results to an Excel file.
#' @param d Dataframe with differential expression results from presto::wilcoxauc()
#' @param fname Filename of the output Excel file.
write_de <- function(d, fname, n = 500) {
  wb <- openxlsx::createWorkbook()
  #fname <- "analysis/nuclei/louvain_de.xlsx"
  unlink(fname)
  for (this_group in sort(unique(d$group))) {
    openxlsx::addWorksheet(wb, as.character(this_group))
    x <- d %>%
      dplyr::filter(group == this_group) %>%
      dplyr::top_n(n = n, wt = abs(0.5 - auc)) %>%
      dplyr::arrange(-auc) %>%
      # top_n(n = 300, wt = -logFC * (pct_in - pct_out)) %>%
      as.data.frame
    openxlsx::writeDataTable(
      wb,
      as.character(this_group),
      x = x, rowNames = FALSE, tableStyle = "TableStyleLight1"
    )
  }
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}

write_ava_xlsx <- function(d, fname, n = 500) {
  wb <- openxlsx::createWorkbook()
  #fname <- "analysis/nuclei/louvain_de.xlsx"
  unlink(fname)
  for (this_coef in sort(unique(d$coef))) {
    openxlsx::addWorksheet(wb, as.character(this_coef))
    x <- d %>%
      dplyr::filter(coef == this_coef) %>%
      dplyr::top_n(n = n, wt = -log10(P.Value)) %>%
      arrange(P.Value) %>%
      as.data.frame
    openxlsx::writeDataTable(
      wb,
      as.character(this_coef),
      x = x, rowNames = FALSE, tableStyle = "TableStyleLight1"
    )
  }
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)
}
