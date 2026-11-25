#' Standard normal variate (SNV)
#'
#' Centres and scales each row (spectrum) of a matrix: each row ends with mean
#' 0 and standard deviation 1 (computed with \eqn{n - 1}).
#'
#' @param X numeric matrix, one spectrum per row.
#' @return The matrix after SNV.
#' @examples
#' X <- matrix(rnorm(20), 4)
#' rowMeans(snv(X))
#' @export
snv <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  m <- rowMeans(X)
  s <- apply(X, 1, stats::sd)
  (X - m) / s
}

#' Preprocess a list of blocks
#'
#' Applies, in this order and block by block: optional SNV (by row), column
#' centring, optional column scaling, and scaling of the block to unit total
#' inertia (\eqn{\|X_k\|_F = 1}), the default normalisation of the package.
#'
#' @param blocks a list of numeric matrices (or a single matrix) with the same
#'   number of rows.
#' @param center logical, centre the columns.
#' @param scale logical or logical vector (one value per block): scale the
#'   columns to unit standard deviation.
#' @param snv logical or logical vector (one value per block): apply SNV to the
#'   rows first.
#' @param unit_inertia logical, scale each block to unit total inertia.
#' @return A named list of preprocessed matrices.
#' @examples
#' blocks <- list(A = matrix(rnorm(40), 10), B = matrix(rnorm(30), 10))
#' sapply(cf_preprocess(blocks), function(X) sum(X^2))
#' @export
cf_preprocess <- function(blocks, center = TRUE, scale = FALSE, snv = FALSE,
                          unit_inertia = TRUE) {
  if (is.matrix(blocks) || is.data.frame(blocks)) blocks <- list(blocks)
  K <- length(blocks)
  if (is.null(names(blocks))) names(blocks) <- paste0("Block", seq_len(K))
  scale <- rep_len(scale, K)
  use_snv <- rep_len(snv, K)
  n <- unique(vapply(blocks, NROW, integer(1)))
  if (length(n) != 1) stop("all blocks must have the same number of rows")
  out <- vector("list", K)
  for (k in seq_len(K)) {
    X <- as.matrix(blocks[[k]])
    storage.mode(X) <- "double"
    if (anyNA(X)) stop("block ", names(blocks)[k], " contains missing values")
    if (use_snv[k]) X <- snv(X)
    if (center) X <- sweep(X, 2, colMeans(X))
    if (scale[k]) {
      s <- apply(X, 2, stats::sd)
      s[s == 0] <- 1
      X <- sweep(X, 2, s, "/")
    }
    if (unit_inertia) X <- X / sqrt(sum(X^2))
    out[[k]] <- X
  }
  names(out) <- names(blocks)
  out
}
