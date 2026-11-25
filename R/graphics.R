#' Graphics module of canofac
#'
#' Base-graphics functions for the results of the package. All labels follow
#' the language set by \code{\link{cf_lang}} (English or French).
#' \itemize{
#'   \item \code{cf_plot_scores}: individuals on two reference components,
#'     coloured by group, with a normal ellipse per group;
#'   \item \code{cf_plot_saliences}: table of the saliences \eqn{d_{kh}^2}
#'     (or of \eqn{d_{kh}}, or of \eqn{\delta_{kh}}), blocks in rows and
#'     dimensions in columns, drawn as a heat map with the values;
#'   \item \code{cf_plot_trajectories}: trajectory of each block in the
#'     \eqn{(\delta_{k1}, \delta_{k2})} plane as \eqn{q} varies, with markers
#'     at the values of \eqn{q} of the classical methods (\eqn{q = 1},
#'     \eqn{q = 2}, \eqn{q = 4}) and at the largest \eqn{q};
#'   \item \code{cf_plot_test}: p-values of \code{\link{cf_test}} (blocks in
#'     rows, components or dimensions in columns, plus the subspace test),
#'     linked blocks in bold and the outcome (all, several, one or no block
#'     linked) above each column, that of the subspace in capitals;
#'   \item \code{cf_plot_compare}: pseudo-singular values of several fits on
#'     the same blocks (for example \eqn{q = 1, 2, 4}), one group of bars per
#'     block, the dimensions of each fit being matched to those of the first.
#' }
#' The \code{plot} methods of \code{canofac}, \code{cf_trajectories} and
#' \code{cf_test} objects call these functions.
#'
#' @param x a \code{canofac} object (\code{cf_plot_scores},
#'   \code{cf_plot_saliences}), a \code{cf_trajectories} object
#'   (\code{cf_plot_trajectories}) or a \code{cf_test} object
#'   (\code{cf_plot_test}).
#' @param dims the two dimensions to plot.
#' @param group optional factor used to colour the individuals.
#' @param ellipse logical, draw a normal ellipse for each group (at least
#'   three individuals).
#' @param level coverage of the ellipses.
#' @param labels logical or character vector: label the individuals
#'   (\code{cf_plot_scores}) or the blocks (\code{cf_plot_trajectories}).
#' @param main optional title; a default title in the current language is
#'   used when \code{NULL}.
#' @param ... further graphical parameters.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(1)
#' u <- scale(rnorm(30))[, 1]
#' g <- factor(ifelse(u > 0, "high", "low"))
#' blocks <- cf_preprocess(lapply(c(5, 8, 12), function(p) outer(u, rnorm(p)) +
#'                                  matrix(rnorm(30 * p), 30)))
#' fit <- canofac(blocks, R = 2, q = 2)
#' cf_plot_scores(fit, group = g)
#' cf_plot_saliences(fit)
#' tr <- cf_trajectories(blocks, R = 2, qgrid = c(1, 1.5, 2, 3, 4, 6))
#' cf_plot_trajectories(tr)
#' fits <- list("q = 1" = canofac(blocks, 2, 1), CPCA = fit,
#'              INDORT = canofac(blocks, 2, 4))
#' cf_plot_compare(fits, dim = 1)
#' old <- cf_lang("fr")
#' cf_plot_saliences(fit)
#' cf_lang(old)
#' @name cf_graphics
NULL

# Normal ellipse of a set of 2-d points.
.ellipse_points <- function(P, level = 0.95, n = 100) {
  m <- colMeans(P)
  e <- eigen(stats::cov(P), symmetric = TRUE)
  r <- sqrt(stats::qchisq(level, 2))
  th <- seq(0, 2 * pi, length.out = n)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), 2)
  t(m + r * A %*% rbind(cos(th), sin(th)))
}

#' @rdname cf_graphics
#' @export
cf_plot_scores <- function(x, dims = c(1, 2), group = NULL, ellipse = !is.null(group),
                           level = 0.95, labels = FALSE, main = NULL, ...) {
  if (!inherits(x, "canofac")) stop("x must be a canofac object")
  U <- x$U[, dims, drop = FALSE] * sqrt(nrow(x$U))
  col <- "black"
  if (!is.null(group)) {
    group <- as.factor(group)
    pal <- .palette(nlevels(group))
    col <- pal[as.integer(group)]
  }
  lim <- range(U)
  if (ellipse && !is.null(group)) {
    E <- lapply(levels(group), function(l) {
      P <- U[group == l, , drop = FALSE]
      if (nrow(P) >= 3) .ellipse_points(P, level) else NULL
    })
    lim <- range(lim, unlist(E))
  }
  graphics::plot(U[, 1], U[, 2], pch = 19, col = col, asp = 1,
                 xlim = lim, ylim = lim,
                 xlab = paste(.tr("ref_comp"), dims[1]),
                 ylab = paste(.tr("ref_comp"), dims[2]),
                 main = if (is.null(main)) .tr("scores_title") else main, ...)
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey50")
  if (ellipse && !is.null(group))
    for (i in seq_along(E)) if (!is.null(E[[i]]))
      graphics::lines(E[[i]], col = pal[i], lwd = 1.5)
  if (!isFALSE(labels)) {
    lab <- if (is.character(labels)) labels else
      if (!is.null(rownames(x$U))) rownames(x$U) else seq_len(nrow(U))
    graphics::text(U[, 1], U[, 2], lab, pos = 3, cex = 0.7, col = col)
  }
  if (!is.null(group))
    graphics::legend("topright", levels(group), pch = 19, col = pal,
                     bty = "n", cex = 0.8)
  invisible(x)
}

# Heat map of a matrix with its values; rows from top to bottom.
.heat <- function(M, txt, main, bold = NULL, zlim = NULL, top = NULL,
                  palette = "YlOrRd", cex = 0.8) {
  nr <- nrow(M); nc <- ncol(M)
  if (is.null(zlim)) zlim <- range(M, finite = TRUE)
  cols <- grDevices::hcl.colors(20, palette, rev = TRUE)
  mar_left <- 1 + 0.55 * max(nchar(rownames(M)), 4)
  op <- graphics::par(mar = c(3, mar_left, if (is.null(top)) 4 else 6, 1))
  on.exit(graphics::par(op))
  Z <- t(M[rev(seq_len(nr)), , drop = FALSE])
  graphics::image(seq_len(nc), seq_len(nr), Z, col = cols, zlim = zlim,
                  axes = FALSE, xlab = "", ylab = "")
  graphics::title(main = main, line = if (is.null(top)) 1.5 else 4, cex.main = 0.95)
  graphics::axis(1, at = seq_len(nc), labels = colnames(M), tick = FALSE, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(nr), labels = rev(rownames(M)), las = 1,
                 tick = FALSE, cex.axis = 0.8)
  if (!is.null(top))
    graphics::axis(3, at = seq_len(nc), labels = top, tick = FALSE, cex.axis = 0.75,
                   las = 1, padj = 0)
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    f <- if (!is.null(bold) && isTRUE(bold[i, j])) 2 else 1
    graphics::text(j, nr - i + 1, txt[i, j], font = f, cex = cex)
  }
  graphics::box()
}

#' @rdname cf_graphics
#' @param stat \code{"saliences"} (\eqn{d_{kh}^2}), \code{"d"} or
#'   \code{"delta"}.
#' @param digits number of decimals of the printed values.
#' @export
cf_plot_saliences <- function(x, stat = c("saliences", "d", "delta"), digits = 2,
                              main = NULL) {
  if (!inherits(x, "canofac")) stop("x must be a canofac object")
  stat <- match.arg(stat)
  M <- switch(stat, saliences = x$saliences, d = x$d, delta = x$delta)
  colnames(M) <- paste(.tr("dim"), seq_len(ncol(M)))
  txt <- matrix(formatC(M, format = "f", digits = digits), nrow(M))
  if (is.null(main))
    main <- switch(stat, saliences = .tr("sal_title"), d = .tr("d_title"),
                   delta = .tr("delta_title"))
  .heat(M, txt, main, zlim = c(0, max(M, 1e-12)), palette = "Blues")
  invisible(x)
}

# Name of the classical method for a value of q, axes and extraction.
.method_name <- function(q, axes, extraction) {
  if (axes == "free") {
    if (q == 1 && extraction == "sequential") return("MAXBET")
    if (q == 2) return("CPCA")
    if (q == 4) return(if (extraction == "sequential") "ACCPS" else "INDORT")
  } else {
    if (q == 1 && extraction == "simultaneous") return("simultaneous MAXBET")
    if (q == 2) return(if (extraction == "sequential") "MCOA" else "simultaneous MCOA")
  }
  ""
}

#' @rdname cf_graphics
#' @param mark values of \eqn{q} marked on the trajectories, in addition to the
#'   smallest and the largest value of the grid.
#' @param legend logical, add a legend of the markers.
#' @export
cf_plot_trajectories <- function(x, dims = c(1, 2), labels = TRUE, mark = c(2, 4),
                                 main = NULL, legend = TRUE, ...) {
  if (!inherits(x, "cf_trajectories")) stop("x must be a cf_trajectories object")
  D <- x$delta
  K <- dim(D)[1]
  q <- x$q
  nq <- length(q)
  col <- .palette(K)
  axes <- if (is.null(x$settings)) "free" else x$settings$axes
  extr <- if (is.null(x$settings)) "simultaneous" else x$settings$extraction
  # markers: smallest q (open circle), marked q (square, triangle, ...), largest q (dot)
  mq <- c(q[1], intersect(mark, q[-c(1, nq)]), q[nq])
  pch <- c(1, rep_len(c(15, 17, 18, 8, 3, 4), length(mq) - 2), 19)
  lab <- vapply(mq, function(v) {
    m <- .method_name(v, axes, extr)
    if (nzchar(m)) sprintf("q = %g (%s)", v, m) else sprintf("q = %g", v)
  }, character(1))
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
                 xlab = bquote(delta[.(dims[1])]),
                 ylab = bquote(delta[.(dims[2])]),
                 main = if (is.null(main)) .tr("traj_title") else main, ...)
  graphics::abline(0, 1, lty = 3, col = "grey60")
  for (k in seq_len(K)) {
    graphics::lines(D[k, dims[1], ], D[k, dims[2], ], col = col[k], lwd = 1.5)
    for (i in seq_along(mq)) {
      j <- match(mq[i], q)
      graphics::points(D[k, dims[1], j], D[k, dims[2], j], pch = pch[i], col = col[k])
    }
    if (!isFALSE(labels)) {
      nm <- if (is.character(labels)) labels[k] else dimnames(D)[[1]][k]
      graphics::text(D[k, dims[1], nq], D[k, dims[2], nq], nm,
                     pos = 4, cex = 0.7, col = col[k])
    }
  }
  if (legend)
    graphics::legend("bottomright", lab, pch = pch, bty = "n", cex = 0.8)
  invisible(x)
}

#' @rdname cf_graphics
#' @param subspace logical, add a column with the p-values of the subspace test.
#' @export
cf_plot_test <- function(x, digits = 1, subspace = TRUE, main = NULL) {
  if (!inherits(x, "cf_test")) stop("x must be a cf_test object")
  P <- x$p
  fixed <- x$settings$statistic %in% c("lobo", "fixed")
  colnames(P) <- paste(.tr(if (fixed) "comp" else "dim"), seq_len(ncol(P)))
  bold <- x$linked
  top <- .link_label(x$link)
  if (subspace) {
    P <- cbind(P, x$p_subspace)
    colnames(P)[ncol(P)] <- if (.lang() == "fr") "Sous-espace" else "Subspace"
    bold <- cbind(bold, if (is.null(x$linked_subspace)) FALSE else x$linked_subspace)
    top <- c(top, if (is.null(x$link_subspace)) "" else
      toupper(.link_label(x$link_subspace)))
  }
  txt <- matrix(formatC(100 * P, format = "f", digits = digits), nrow(P))
  zmax <- -log10(1 / (x$settings$B + 1))
  if (is.null(main)) main <- .tr("pvalues_title")
  .heat(-log10(P), txt, main, bold = bold, zlim = c(0, zmax), top = top)
  invisible(x)
}

#' @rdname cf_graphics
#' @param fits a named list of \code{canofac} objects fitted on the same
#'   blocks (\code{cf_plot_compare}).
#' @param dim dimension shown (\code{cf_plot_compare}), in the numbering of
#'   the first fit.
#' @export
cf_plot_compare <- function(fits, dim = 1, main = NULL, ...) {
  if (!all(vapply(fits, inherits, logical(1), "canofac")))
    stop("fits must be a list of canofac objects")
  if (is.null(names(fits)))
    names(fits) <- vapply(fits, function(f) {
      m <- .method_name(f$settings$q, f$settings$axes, f$settings$extraction)
      if (nzchar(m)) m else sprintf("q = %g", f$settings$q)
    }, character(1))
  ref <- fits[[1]]$U
  M <- vapply(fits, function(f) f$d[, .match_dims(f$U, ref)[dim]],
              numeric(nrow(fits[[1]]$d)))
  M <- matrix(M, ncol = length(fits), dimnames = list(rownames(fits[[1]]$d), names(fits)))
  graphics::barplot(t(M), beside = TRUE, las = 2, col = .palette(ncol(M)),
                    ylab = .tr("pvs"), ylim = c(0, max(M, 1e-12) * 1.15),
                    main = if (is.null(main)) sprintf(.tr("compare_title"), dim) else main,
                    legend.text = colnames(M),
                    args.legend = list(bty = "n", cex = 0.8, horiz = TRUE, x = "top"), ...)
  invisible(fits)
}


#' Permutation distributions of the subspace statistic
#'
#' One histogram per block of the permuted values of the subspace statistic
#' of \code{\link{cf_test}} (restricted RV by default), with the observed
#' value as a vertical line. Also \code{plot(x, type = "null")}.
#'
#' @param x a \code{cf_test} object (not computed with \code{pvalue = "gamma"}).
#' @param breaks number of bins.
#' @param main optional title.
#' @return \code{x}, invisibly.
#' @export
cf_plot_null <- function(x, breaks = 30, main = NULL) {
  if (!inherits(x, "cf_test")) stop("x must be a cf_test object")
  if (length(x$null) == 0) stop("no permutation distribution stored (pvalue = 'gamma')")
  nm <- names(x$null)
  # grid of at most 4 columns, so that 10 blocks (Friday) still fit in a standard figure
  nc <- min(4, length(nm))
  op <- graphics::par(mfrow = c(ceiling(length(nm) / nc), nc), oma = c(0, 0, 2, 0),
                      mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.6, 0), cex = 0.8)
  on.exit(graphics::par(op))
  for (b in nm) {
    d <- x$null[[b]]
    graphics::hist(d$permuted, breaks = breaks, col = "#9DB9E6", border = "white",
                   xlim = range(c(d$permuted, d$observed)),
                   main = sprintf("%s  p = %.3g", b, x$p_subspace[[b]]), cex.main = 0.9,
                   xlab = .tr("statistic"), ylab = .tr("frequency"))
    graphics::abline(v = d$observed, col = "#B4491F", lwd = 2)
  }
  graphics::mtext(if (is.null(main)) .tr("null_title") else main, outer = TRUE, cex = 0.9)
  invisible(x)
}
