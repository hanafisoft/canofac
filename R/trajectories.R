#' Trajectory of each block as the power q varies
#'
#' Fits the q-family on a grid of powers and returns, for each block, the
#' normalised pseudo-singular values \eqn{\delta_{kh} = d_{kh}/\sigma_1(X_k)},
#' which measure how close the reference component \eqn{u_h} is to the dominant
#' direction of block \eqn{k}. The dimensions are followed from the largest
#' value of q: at each q, they are matched to those obtained at the largest q.
#'
#' @inheritParams canofac
#' @param qgrid numeric vector of powers (all \code{>= 1}).
#' @return An object of class \code{"cf_trajectories"}: a list with \code{q} (sorted grid), \code{delta} and \code{d}
#'   (arrays \eqn{K \times R \times} length of the grid), \code{sigma1}, \code{criterion} and \code{settings}.
#' @examples
#' set.seed(1)
#' blocks <- cf_preprocess(lapply(c(4, 6, 8), function(p) matrix(rnorm(15 * p), 15)))
#' tr <- cf_trajectories(blocks, R = 2, qgrid = c(1, 2, 4, 8))
#' tr$delta[, , "q=4"]
#' plot(tr)
#' @export
cf_trajectories <- function(blocks, R = 2, qgrid = c(1, 1.5, 1.9, 2, 2.1, 3, 4, 6, 8, 10),
                            axes = c("free", "orthonormal"),
                            extraction = c("simultaneous", "sequential"),
                            tol = 1e-10, maxit = 5000, nstart = NULL) {
  axes <- match.arg(axes)
  extraction <- match.arg(extraction)
  qgrid <- sort(unique(qgrid))
  if (any(qgrid < 1)) stop("all values of qgrid must be >= 1")
  K <- length(blocks)
  if (is.null(names(blocks))) names(blocks) <- paste0("Block", seq_len(K))
  sigma1 <- vapply(blocks, function(X) svd(X, nu = 0, nv = 0)$d[1], numeric(1))
  labs <- paste0("q=", qgrid)
  dn <- list(names(blocks), paste0("Dim", seq_len(R)), labs)
  d <- array(NA_real_, c(K, R, length(qgrid)), dn)
  crit <- stats::setNames(numeric(length(qgrid)), labs)
  Uref <- NULL
  for (i in rev(seq_along(qgrid))) {
    fit <- .cf_fit(blocks, R, qgrid[i], axes, extraction, tol, maxit,
                   .nstart(nstart, qgrid[i], axes))
    if (is.null(Uref)) {
      Uref <- fit$U
      pi <- seq_len(R)
    } else {
      pi <- .match_dims(fit$U, Uref)
    }
    d[, , i] <- fit$d[, pi, drop = FALSE]
    crit[i] <- sum(abs(fit$d)^qgrid[i])
  }
  out <- list(q = qgrid, delta = sweep(d, 1, sigma1, "/"), d = d,
              sigma1 = sigma1, criterion = crit,
              settings = list(R = R, axes = axes, extraction = extraction))
  class(out) <- "cf_trajectories"
  out
}
