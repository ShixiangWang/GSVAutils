#' Apply GSVA Method to a Expression Matrix List
#'
#' @param GeneSets a `data.frame` containing groups and genes.
#' @param group_col the column name for groups.
#' @param gene_col the column name for genes.
#' @param ExprMatList a `list` of expression matrix represented with `tibble`.
#' See also [GenTibbleList()].
#' @inheritParams GSVA::gsva
#' @param ... other arguments passing to [GSVA::gsva()]
#'
#' @return a result `list`
#' @export
#'
#' @examples
#' data("example_gsets")
#' data("example_data")
#'
#' ExprList <- GenTibbleList(list(example_data))
#' res <- ApplyGSVA(example_gsets,
#'   group_col = "Cell_type",
#'   gene_col = "Symbol", ExprMatList = ExprList, method = "gsva"
#' )
ApplyGSVA <- function(GeneSets, group_col, gene_col, ExprMatList,
                      method = c("ssgsea", "gsva", "zscore", "plage"),
                      kcdf = c("Gaussian", "Poisson"), ...) {
  stopifnot(inherits(GeneSets, "tbl_df") &
    inherits(group_col, "character") &
    inherits(gene_col, "character") &
    inherits(ExprMatList, "list"))

  method <- match.arg(method)
  kcdf <- match.arg(kcdf)
  i <- 1
  resList <- list()
  groups <- names(table(GeneSets[, group_col]))
  gset_list <- lapply(groups, function(x) {
    GeneSets[GeneSets[, group_col] == x, gene_col] %>%
      unlist() %>%
      as.character()
  })

  names(gset_list) <- groups

  if (length(ExprMatList) > 1) {
    for (expr_mat in ExprMatList) {
      if (!inherits(expr_mat, "tbl_df")) {
        stop("All elements of ExprMatList should be class tibble!")
      }
      expr_mat <- as.data.frame(expr_mat)
      rownames(expr_mat) <- expr_mat[, 1]
      expr_mat <- expr_mat[, -1] %>% as.matrix()

      res <- gsva(expr = expr_mat, gset.idx.list = gset_list, method = method, kcdf = kcdf, ...)
      res <- as.data.frame(t(res))

      resList[[i]] <- res
      names(resList)[i] <- names(ExprMatList)[i]
      i <- i + 1
    }
  } else {
    for (expr_mat in ExprMatList) {
      if (!inherits(expr_mat, "tbl_df")) {
        stop("All elements of ExprMatList should be class tibble!")
      }
      expr_mat <- as.data.frame(expr_mat)
      rownames(expr_mat) <- expr_mat[, 1]
      expr_mat <- expr_mat[, -1] %>% as.matrix()

      res <- gsva(expr = expr_mat, gset.idx.list = gset_list, method = method, kcdf = kcdf, ...)
      res <- as.data.frame(t(res))
      resList[[i]] <- res
      names(resList)[i] <- names(ExprMatList)[i]
      i <- i + 1
    }
  }
  return(resList)
}
