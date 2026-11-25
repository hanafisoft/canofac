#' Test along the number of dimensions (EXPERIMENTAL choice of R)
#'
#' \strong{Experimental.} The choice of the number of dimensions \eqn{R} is not
#' part of the published articles (pseudo-singular values; LOBO test): a
#' warning is issued at each call, and the printed result and the plot say so.
#' The models are re-estimated for each \eqn{r} and are not nested for
#' simultaneous methods (with INDORT, \eqn{u_1} may rotate by more than 30
#' degrees between \eqn{R = 1} and \eqn{R = 2}); the successive tests are not
#' corrected for multiplicity. Model \eqn{r} is tested with the seed
#' \code{seed + r - 1}. French help: \code{cf_help("cf_test_path", "fr")}.
#'
#' Fits and tests the models with \eqn{r = 1, \dots, R_{max}} dimensions
#' (each model estimated anew, so that in simultaneous extraction the
#' dimensions of model \eqn{r} are estimated together), with
#' \code{\link{cf_test}} for each model. By default
#' \eqn{R_{max} = \min(K + c, n - 1)}: at most \eqn{c} directions linking
#' several blocks and one dimension of its own for each block.
#'
#' The procedure does not stop at the first model containing a dimension
#' without signal, because a signal may appear in a larger model (in
#' sequential extraction, the order of the dimensions follows their
#' contribution to the criterion, not their significance; in simultaneous
#' extraction, a dimension of model \eqn{r} may mix two structures that model
#' \eqn{r + 1} separates). The chosen number of dimensions is the largest
#' \eqn{r} whose dimensions are all retained: with the LOBO
#' statistic (default), every component of the other blocks is linked to at
#' least one block (\eqn{\min_k p_{kh} \le \alpha / K}); with the
#' re-estimation statistic, every dimension has at least one linked block (see \code{\link{cf_test}}). Since this choice uses the same data as the tests, the
#' p-values of the chosen model are slightly optimistic: report the whole
#' profile.
#'
#' @inheritParams cf_test
#' @param Rmax largest number of dimensions; by default
#'   \code{min(K + c, n - 1)}.
#' @param c assumed maximal number of directions linking several blocks, used for the default
#'   \code{Rmax}.
#' @return An object of class \code{"cf_test_path"}: \code{tests} (list of
#'   \code{cf_test} objects, one per model), \code{profile} (data frame with,
#'   for each \eqn{r}, the number of retained dimensions, of dimensions linked to all the blocks, to several and to one block, and
#'   whether all dimensions are retained), \code{R} (chosen number of
#'   dimensions, 0 if none), \code{Rmax} and \code{experimental = TRUE}.
#' @examples
#' set.seed(3)
#' u <- scale(rnorm(20))[, 1]
#' blocks <- cf_preprocess(lapply(c(5, 8, 6), function(p) outer(u, rnorm(p)) +
#'                                  matrix(rnorm(20 * p), 20)))
#' pth <- suppressWarnings(cf_test_path(blocks, q = 2, B = 49, seed = 1))
#' pth
#' plot(pth)
#' @export
cf_test_path <- function(blocks, Rmax = NULL, c = 2, q = 2,
                         axes = c("free", "orthonormal"),
                         extraction = c("simultaneous", "sequential"),
                         scheme = c("rows", "groups", "within", "variables", "components"),
                         groups = NULL, B = 999, alpha = 0.05, statistic = c("lobo", "refit", "fixed"),
                         matching = TRUE, seed = NULL,
                         tol = 1e-10, maxit = 5000, nstart = NULL,
                         subspace = c("inertia", "equal"),
                         weighting = c("none", "mfa", "statis")) {
  warning(.tr("experimental_warning"), call. = FALSE)
  statistic <- match.arg(statistic)
  if (statistic == "fixed") statistic <- "lobo"   # alias (canofac < 0.6.0)
  subspace <- match.arg(subspace)
  weighting <- match.arg(weighting)
  axes <- match.arg(axes)
  extraction <- match.arg(extraction)
  scheme <- match.arg(scheme)
  K <- length(blocks)
  n <- NROW(blocks[[1]])
  if (is.null(Rmax)) Rmax <- min(K + c, n - 1)
  Rmax <- min(Rmax, n - 1)
  seed <- .draw_seed(seed)
  tests <- vector("list", Rmax)
  for (r in seq_len(Rmax)) {
    tests[[r]] <- cf_test(blocks, R = r, q = q, axes = axes,
                          extraction = extraction, scheme = scheme,
                          groups = groups, B = B,
                          alpha = alpha, statistic = statistic,
                          matching = matching, seed = seed + r - 1,
                          tol = tol, maxit = maxit, nstart = nstart,
                          subspace = subspace, weighting = weighting)
  }
  names(tests) <- paste0("R=", seq_len(Rmax))
  profile <- data.frame(
    r = seq_len(Rmax),
    retained = vapply(tests, function(t) sum(t$retained), numeric(1)),
    all = vapply(tests, function(t) sum(t$link == "all"), numeric(1)),
    several = vapply(tests, function(t) sum(t$link == "several"), numeric(1)),
    one = vapply(tests, function(t) sum(t$link == "one"), numeric(1)),
    all_retained = vapply(tests, function(t) all(t$retained), logical(1)),
    row.names = NULL)
  ok <- which(profile$all_retained)
  out <- list(tests = tests, profile = profile,
              R = if (length(ok)) max(ok) else 0L, Rmax = Rmax, experimental = TRUE,
              settings = list(c = c, q = q, axes = axes, extraction = extraction,
                              scheme = scheme, B = B, alpha = alpha,
                              statistic = statistic, seed = seed))
  class(out) <- "cf_test_path"
  out
}

#' @export
print.cf_test_path <- function(x, ...) {
  s <- x$settings
  cat(.tr("experimental_title"), "\n", .tr("experimental_note"), "\n", sep = "")
  cat(.tr("title_path"), "\n", sep = "")
  cat(sprintf(.tr("settings_test"), s$q, s$axes, s$extraction, s$scheme, s$B, s$alpha))
  cat(sprintf("  Rmax = %d\n\n", x$Rmax))
  print(x$profile, row.names = FALSE)
  cat(sprintf(.tr("chosen_R"), x$R))
  cat("\n", .tr("experimental_title"), "\n", sep = "")
  invisible(x)
}

#' @rdname plot.canofac
#' @export
plot.cf_test_path <- function(x, ...) {
  Rmax <- x$Rmax
  K <- nrow(x$tests[[1]]$p)
  M <- matrix(NA_real_, Rmax, Rmax)          # rows: model r, columns: dimension h
  Ret <- matrix(FALSE, Rmax, Rmax)           # retained dimensions (bold)
  for (r in seq_len(Rmax)) {
    M[r, seq_len(r)] <- apply(x$tests[[r]]$p, 2, min)
    Ret[r, seq_len(r)] <- x$tests[[r]]$retained
  }
  Z <- -log10(M)
  cols <- grDevices::hcl.colors(20, "YlOrRd", rev = TRUE)
  graphics::image(seq_len(Rmax), seq_len(Rmax), t(Z), col = cols,
                  xlab = .tr("dim_h"), ylab = .tr("n_dims"),
                  main = paste0(.tr("experimental_title"), "\n", .tr("path_title")),
                  cex.main = 0.9, ...)
  for (r in seq_len(Rmax)) for (h in seq_len(r)) {
    graphics::text(h, r, sprintf("%.1f", 100 * M[r, h]), cex = 0.7,
                   font = if (Ret[r, h]) 2 else 1)
  }
  graphics::abline(h = x$R + 0.5, lty = 2)
  invisible(x)
}
