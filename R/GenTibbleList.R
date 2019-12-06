#' Generate a List of Tibbles from a List of Gene Expression Sets
#'
#' @param ExpressionSetList a list of [Biobase::ExpressionSet()]
#'
#' @author Shixiang Wang
#' @return a `list`
#' @export
#'
#' @examples
#' data("example_data")
#' GenTibbleList(list(example_data))
GenTibbleList <- function(ExpressionSetList) {
  stopifnot(is.list(ExpressionSetList), class(ExpressionSetList[[1]]) == "ExpressionSet")

  res <- list()
  i <- 1
  for (gset in ExpressionSetList) {
    eset <- exprs(gset)
    # find gene symbol column
    fdata <- fData(gset)
    symbol_col <- grep("^gene.?symbol", colnames(fdata), value = TRUE, ignore.case = TRUE)

    if (length(symbol_col) == 0) {
      message("Find nothing about gene symbol in fData, try search it...")
      symbol_col2 <- grep("^gene_assignment", colnames(fdata), value = TRUE, ignore.case = TRUE)
      message("Find ", symbol_col2)

      message("Processing...")
      strlist <- strsplit(fdata[, symbol_col2], split = " // ")
      rowname <- sapply(strlist, function(x) trimws(x[2]))
      rownames(eset) <- rowname

      # stop("Something wrong with your fData of input List, please check it")
    }
    if (length(symbol_col) > 1) {
      warning("Multiple columns of fData match gene symbol, only use the first one")
      symbol_col <- symbol_col[1]
      rownames(eset) <- fdata[, symbol_col]
    } else if (length(symbol_col) == 1) {
      rownames(eset) <- fdata[, symbol_col]
    }



    # remove duplicate rows, keep the one with biggest mean value
    eset %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      mutate(
        mean_expr = rowMeans(.[, -1], na.rm = TRUE),
        rowname = sub("^(\\w+)\\..*", "\\1", rowname)
      ) %>%
      arrange(rowname, desc(mean_expr)) %>%
      distinct(rowname, .keep_all = TRUE) %>%
      select(-mean_expr) %>%
      as.tibble() -> res[[i]]


    i <- i + 1
  }

  message("Done.")
  return(res)
}
