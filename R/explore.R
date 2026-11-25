#' Multimodality of the criterion: criterion reached from random starts
#'
#' Runs the algorithm with unconstrained axes from \code{nstart} random
#' centred orthonormal starts, drawn by the common generator of canofac with
#' the seed \code{seed} (same starts as the Python version), and from the
#' deterministic CPCA start, and counts the starts that reach the best
#' criterion to within \code{rtol} in relative terms, and the number of
#' distinct optima (article on pseudo-singular values, multimodality study).
#' In sequential extraction, each step starts from the same start.
#'
#' @inheritParams canofac
#' @param nstart number of random starts.
#' @param seed seed of the common generator.
#' @param rtol relative tolerance to identify the best criterion and distinct
#'   optima.
#' @return An object of class \code{"cf_multistart"}: \code{criteria}
#'   (random starts), \code{deterministic} (CPCA start), \code{best},
#'   \code{share_best} (\%), \code{n_optima}, \code{optima},
#'   \code{deterministic_best} and the settings. French help:
#'   \code{cf_help("cf_multistart", "fr")}.
#' @examples
#' bl <- cf_preprocess(lapply(c(4, 6, 5), function(p) matrix(rnorm(15 * p), 15)))
#' ms <- cf_multistart(bl, R = 2, q = 4, nstart = 20)
#' ms
#' plot(ms)
#' @export
cf_multistart <- function(blocks, R = 2, q = 4, nstart = 100,
                          extraction = c("simultaneous", "sequential"),
                          seed = 2026, rtol = 1e-5, tol = 1e-10, maxit = 5000) {
  extraction <- match.arg(extraction)
  blocks <- lapply(blocks, function(X) { X <- as.matrix(X); storage.mode(X) <- "double"; X })
  n <- nrow(blocks[[1]])
  r0 <- if (extraction == "simultaneous") R else 1
  starts <- .starts(n, r0, nstart + 1, seed = seed)
  run_sim <- function(U0) sum(.free_simultaneous(blocks, R, q, tol, maxit, U0 = U0)$d^q)
  run_seq <- function(U0) {
    Xc <- blocks
    tot <- 0
    for (h in seq_len(R)) {
      f <- .free_simultaneous(Xc, 1, q, tol, maxit, U0 = U0)
      tot <- tot + sum(f$d^q)
      u <- f$U[, 1]
      Xc <- lapply(Xc, function(X) X - outer(u, drop(crossprod(u, X))))
    }
    tot
  }
  run <- if (extraction == "simultaneous") run_sim else run_seq
  crits <- vapply(starts, run, numeric(1))
  det <- run(NULL)
  best <- max(c(crits, det))
  reach <- crits >= best * (1 - rtol)
  vals <- sort(crits, decreasing = TRUE)
  optima <- vals[1]
  for (v in vals[-1]) if (abs(v - optima[length(optima)]) > rtol * abs(best)) optima <- c(optima, v)
  out <- list(criteria = crits, deterministic = det, best = best,
              share_best = 100 * mean(reach), n_optima = length(optima), optima = optima,
              deterministic_best = det >= best * (1 - rtol),
              settings = list(R = R, q = q, nstart = nstart, extraction = extraction, seed = seed))
  class(out) <- "cf_multistart"
  out
}

#' @export
print.cf_multistart <- function(x, ...) {
  cat(.tr("multistart_title"), "\n", sep = "")
  cat(sprintf(.tr("multistart_head"), x$share_best, x$n_optima,
              if (x$deterministic_best) .tr("yes") else .tr("no")))
  invisible(x)
}

#' @export
plot.cf_multistart <- function(x, ...) {
  cr <- sort(x$criteria, decreasing = TRUE)
  graphics::plot(seq_along(cr), cr, pch = 19, cex = 0.6, col = .palette(1),
                 xlab = .tr("start"), ylab = .tr("criterion_lab"),
                 main = .tr("multistart_title"), cex.main = 0.95, ...)
  graphics::abline(h = x$deterministic, lty = 2, col = .palette(2)[2])
  invisible(x)
}

#' Angles between the components of two fits
#'
#' Matches the components of \code{fit2} to those of \code{fit1} (optimal
#' assignment) and gives the angle (degrees) between each pair of matched
#' components, and the largest principal angle between the two subspaces
#' (article on pseudo-singular values: sequential versus simultaneous
#' extraction).
#'
#' @param fit1,fit2 \code{canofac} objects with the same number of rows.
#' @return A list (printed): \code{components} (angles) and \code{subspace}.
#' @examples
#' bl <- cf_preprocess(lapply(c(4, 6, 5), function(p) matrix(rnorm(15 * p), 15)))
#' cf_angles(canofac(bl, R = 2, q = 4), canofac(bl, R = 2, q = 4, extraction = "sequential"))
#' @export
cf_angles <- function(fit1, fit2) {
  U1 <- fit1$U; U2 <- fit2$U
  pm <- .match_dims(U2, U1)
  cs <- pmin(1, abs(colSums(U1 * U2[, pm, drop = FALSE])))
  comp <- stats::setNames(acos(cs) * 180 / base::pi, paste0("Dim", seq_along(cs)))
  sv <- pmin(1, svd(crossprod(U1, U2), nu = 0, nv = 0)$d)
  sub <- acos(min(sv)) * 180 / base::pi
  cat(.tr("angles_head"), "\n", sep = "")
  print(round(comp, 1))
  cat(sprintf(.tr("subspace_angle"), sub))
  invisible(list(components = comp, subspace = sub))
}

#' Pseudo-singular values of several fits
#'
#' Table of the pseudo-singular values of several fits, the dimensions being
#' matched to those of the first fit (article on pseudo-singular values,
#' tables by method). The names of the classical methods are given
#' automatically.
#'
#' @param fits list of \code{canofac} objects.
#' @param names optional names of the fits.
#' @param digits number of digits printed.
#' @param match logical: match the dimensions to those of the first fit
#'   (default); with \code{FALSE}, each fit keeps its own order (decreasing
#'   contribution to its criterion), as in the tables of the article.
#' @return A matrix (printed): one row per block, one column per fit and
#'   dimension.
#' @examples
#' bl <- cf_preprocess(lapply(c(4, 6, 5), function(p) matrix(rnorm(15 * p), 15)))
#' cf_compare(list(canofac(bl, R = 2, q = 2), canofac(bl, R = 2, q = 4)))
#' @export
cf_compare <- function(fits, names = NULL, digits = 2, match = TRUE) {
  if (is.null(names))
    names <- vapply(fits, function(f) {
      m <- .method_name(f$settings$q, f$settings$axes, f$settings$extraction)
      if (nzchar(m)) m else sprintf("q = %g", f$settings$q)
    }, character(1))
  ref <- fits[[1]]$U
  R <- ncol(ref)
  out <- do.call(cbind, lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    pm <- if (match) .match_dims(f$U, ref) else seq_len(R)
    D <- f$d[, pm, drop = FALSE]
    colnames(D) <- paste0(names[i], " \u00b7 Dim", seq_len(R))
    D
  }))
  cat(.tr("compare_table_head"), "\n", sep = "")
  print(round(out, digits))
  invisible(out)
}
