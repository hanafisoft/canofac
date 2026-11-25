#' Plots for canofac objects
#'
#' Base graphics for the results of \code{\link{canofac}},
#' \code{\link{cf_trajectories}} and \code{\link{cf_test}}.
#'
#' @name plot.canofac
NULL

.palette <- function(n) {
  base <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
            "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
  rep_len(base, n)
}

#' @param x an object of class \code{"canofac"}, \code{"cf_trajectories"} or
#'   \code{"cf_test"}.
#' @param type for \code{canofac} objects: \code{"scores"} (individuals on the reference
#'   components, see \code{\link{cf_plot_scores}}), \code{"d"} (pseudo-singular values),
#'   \code{"saliences"} (\eqn{d^2}) or \code{"delta"}
#'   (\eqn{d_{kh}/\sigma_1(X_k)}), or \code{"loadings"} (axes \eqn{V_k}
#'   plotted against the variables; for spectra, against the wavelengths or
#'   wavenumbers read from the column names, with a reversed axis when they
#'   decrease). For \code{cf_test} objects: \code{"pvalues"} (map of the
#'   p-values, \code{\link{cf_plot_test}}) or \code{"null"} (permutation
#'   distributions of the subspace statistic, \code{\link{cf_plot_null}}).
#' @param dims the two dimensions to plot (scores and trajectories).
#' @param group optional factor used to colour the individuals.
#' @param block for \code{type = "loadings"}: names or indices of the blocks
#'   to plot (one panel per block); all blocks by default.
#' @param labels logical or character vector: label the individuals.
#' @param ... further graphical parameters.
#' @rdname plot.canofac
#' @export
plot.canofac <- function(x, type = c("scores", "d", "saliences", "delta", "loadings"),
                         dims = c(1, 2), group = NULL, labels = FALSE,
                         block = NULL, ...) {
  type <- match.arg(type)
  if (type == "loadings") {
    if (is.null(block)) block <- seq_along(x$V)
    if (is.character(block)) block <- match(block, names(x$V))
    dims <- dims[dims <= ncol(x$U)]
    op <- graphics::par(mfrow = grDevices::n2mfrow(length(block)))
    on.exit(graphics::par(op))
    for (k in block) {
      V <- x$V[[k]][, dims, drop = FALSE]
      wl <- suppressWarnings(as.numeric(rownames(V)))
      spectral <- !is.null(rownames(V)) && !anyNA(wl)
      xx <- if (spectral) wl else seq_len(nrow(V))
      xlim <- range(xx)
      if (spectral && length(xx) > 1 && xx[1] > xx[length(xx)]) xlim <- rev(xlim)
      cols <- .palette(length(dims))
      lab_axis <- !spectral && nrow(V) <= 30 && !is.null(rownames(V))
      graphics::matplot(xx, V, xaxt = if (lab_axis) "n" else "s", type = if (spectral || nrow(V) > 30) "l" else "h",
                        lty = 1, lwd = if (spectral) 1.2 else 3, col = cols,
                        xlim = xlim,
                        xlab = if (spectral) .tr("wavelength") else .tr("variable"),
                        ylab = .tr("loading"), main = names(x$V)[k], ...)
      graphics::abline(h = 0, lty = 3, col = "grey50")
      if (lab_axis)
        graphics::axis(1, at = xx, labels = rownames(V), las = 2, cex.axis = 0.7)
      graphics::legend("topright", paste(.tr("dim"), dims), col = cols, lty = 1,
                       bty = "n", cex = 0.8)
    }
    return(invisible(x))
  }
  if (type == "scores") {
    cf_plot_scores(x, dims = dims, group = group, labels = labels, ...)
  } else {
    M <- switch(type, d = x$d, saliences = x$saliences, delta = x$delta)
    ylab <- switch(type, d = .tr("pvs"),
                   saliences = .tr("salience"), delta = expression(delta[kh]))
    graphics::barplot(t(M), beside = TRUE, las = 2,
                      col = .palette(ncol(M)), ylab = ylab,
                      legend.text = colnames(M),
                      args.legend = list(bty = "n", cex = 0.8), ...)
  }
  invisible(x)
}

#' @rdname plot.canofac
#' @export
plot.cf_trajectories <- function(x, dims = c(1, 2), labels = TRUE, ...) {
  cf_plot_trajectories(x, dims = dims, labels = labels, ...)
}

#' @rdname plot.canofac
#' @export
plot.cf_test <- function(x, type = c("pvalues", "null"), ...) {
  type <- match.arg(type)
  if (type == "null") return(cf_plot_null(x, ...))
  cf_plot_test(x, ...)
}
