#' @param x A sparse matrix from the Matrix package.
#' @param file A filename that ends in ".gz".
writeMMgz <- function(x, file) {
  stopifnot(
    is(x, "dgCMatrix") | is(x, "ngCMatrix")
  )
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  con <- gzfile(file)
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    con
  )
  close(con)
  data.table::fwrite(
    x         = Matrix::summary(x),
    file      = file,
    append    = TRUE,
    sep       = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}

