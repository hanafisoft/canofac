#' Rose plot of the pseudo-singular values
#'
#' Draws, in each panel, three concentric circles describing the
#' pseudo-singular values \eqn{d_{kh}} of one dimension across the \eqn{K}
#' blocks (one sector per block):
#' \itemize{
#'   \item inner rose: squared deviations \eqn{(d_{kh} - \bar d_h)^2}, whose
#'     mean over the blocks is the variance \eqn{V_h};
#'   \item middle rose: the values \eqn{d_{kh}} themselves, whose mean is
#'     \eqn{\bar d_h}; it carries the block names;
#'   \item outer ring: composition \eqn{\bar d_h^2 / (\bar d_h^2 + V_h)}, the
#'     share of the mean (orange) and of the variance (teal) in
#'     \eqn{\overline{d^2_h}}. A component close to a common structure has a
#'     ring dominated by the mean.
#' }
#' With one \code{canofac} fit, there is one panel per dimension. With a named
#' list of fits (for example \code{q = 1, 2, 4}), there is one panel per fit,
#' for the dimension \code{dim}, the dimensions of each fit being matched to
#' those of the first one. The scales of the roses are common to all panels.
#'
#' @param x a \code{canofac} object or a named list of \code{canofac} objects
#'   fitted on the same blocks.
#' @param dim dimension shown when \code{x} is a list.
#' @param stat \code{"d"} (pseudo-singular values) or \code{"saliences"}
#'   (\eqn{d^2}).
#' @param block_names optional names of the blocks.
#' @param labels logical, print the values on the roses.
#' @param legend logical, add a legend of the blocks under the panels.
#' @param main optional overall title.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(1)
#' u <- scale(rnorm(20))[, 1]
#' blocks <- cf_preprocess(lapply(c(5, 8, 12, 6), function(p) outer(u, rnorm(p)) +
#'                                  matrix(rnorm(20 * p), 20)))
#' cf_rose(canofac(blocks, R = 2, q = 4))
#' fits <- list("q = 1" = canofac(blocks, 2, 1), "CPCA" = canofac(blocks, 2, 2),
#'              "INDORT" = canofac(blocks, 2, 4))
#' cf_rose(fits, dim = 1)
#' @export
cf_rose <- function(x, dim = 1, stat = c("d", "saliences"), block_names = NULL,
                    labels = TRUE, legend = TRUE, main = NULL) {
  stat <- match.arg(stat)
  val <- function(f) if (stat == "d") f$d else f$saliences
  if (inherits(x, "canofac")) {
    M <- val(x)
    panels <- lapply(seq_len(ncol(M)), function(h) M[, h])
    titles <- paste(.tr("dimension"), seq_len(ncol(M)))
  } else {
    if (!all(vapply(x, inherits, logical(1), "canofac")))
      stop("x must be a canofac object or a list of canofac objects")
    ref <- x[[1]]$U
    panels <- lapply(x, function(f) {
      pi <- .match_dims(f$U, ref)
      val(f)[, pi[dim]]
    })
    titles <- if (!is.null(names(x))) names(x) else paste(.tr("fit"), seq_along(x))
    titles <- paste0(titles, " - ", tolower(.tr("dimension")), " ", dim)
  }
  K <- length(panels[[1]])
  if (is.null(block_names)) {
    block_names <- names(panels[[1]])
    if (is.null(block_names)) block_names <- paste0("B", seq_len(K))
  }
  colors <- .block_colors(K)
  dev2 <- lapply(panels, function(v) (v - mean(v))^2)
  vmax <- max(unlist(panels), 1e-12)
  dmax <- max(unlist(dev2), 1e-12)
  sym <- if (stat == "d") "d" else "d^2"

  op <- graphics::par(mfrow = c(1, length(panels)), mar = c(1, 1, 4, 1),
                      oma = c(if (legend) 3 else 0, 0, if (is.null(main)) 0 else 2, 0))
  on.exit(graphics::par(op))
  for (i in seq_along(panels)) {
    v <- panels[[i]]
    graphics::plot(NA, xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45), asp = 1,
                   axes = FALSE, xlab = "", ylab = "")
    .rose(dev2[[i]], dmax, r0 = 0, r1 = 0.32, colors = colors,
          label = labels, cex = 0.55)
    .rose(v, vmax, r0 = 0.40, r1 = 0.72, colors = colors, label = labels,
          cex = 0.55, block_names = block_names)
    .ring(mean(v), mean(dev2[[i]]), r_in = 0.84, r_out = 0.92, label_r = 1.08)
    graphics::title(main = sprintf("%s\n%s = %.3f  V = %.3f", titles[i], .tr("rose_mean"),
                                   mean(v), mean(dev2[[i]])), cex.main = 0.9)
  }
  if (!is.null(main)) graphics::mtext(main, side = 3, outer = TRUE, line = 0.5)
  if (legend) {
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    graphics::legend("bottom", legend = block_names, fill = colors,
                     horiz = K <= 6, ncol = if (K <= 6) K else ceiling(K / 2),
                     bty = "n", cex = 0.8, xpd = TRUE)
  }
  graphics::mtext(sprintf(.tr("rose_bottom"), sym, sym), side = 1, outer = TRUE,
                  line = if (legend) 2 else 0, cex = 0.7)
  invisible(x)
}

.block_colors <- function(K) grDevices::rainbow(K, s = 0.65, v = 0.85)

.sector <- function(theta1, theta2, r_in, r_out, col, n = 80) {
  th <- seq(theta1, theta2, length.out = n)
  graphics::polygon(c(r_out * cos(th), r_in * cos(rev(th))),
                    c(r_out * sin(th), r_in * sin(rev(th))), col = col, border = NA)
}

.ring <- function(m, V, r_in, r_out, label_r, cex = 0.75,
                  col_mean = "#d97706", col_var = "#0f766e") {
  frac <- if (m^2 + V > 0) m^2 / (m^2 + V) else 1
  split <- frac * 2 * pi
  .sector(0, split, r_in, r_out, col_mean)
  if (frac < 1) .sector(split, 2 * pi, r_in, r_out, col_var)
  graphics::text(label_r * cos(split / 2), label_r * sin(split / 2),
                 sprintf("%.0f%%", 100 * frac), col = col_mean, font = 2, cex = cex)
  if (frac < 1) {
    mid <- (split + 2 * pi) / 2
    graphics::text(label_r * cos(mid), label_r * sin(mid),
                   sprintf("%.0f%%", 100 * (1 - frac)), col = col_var, font = 2, cex = cex)
  }
}

.rose <- function(values, vmax, r0, r1, colors, label = TRUE, cex = 0.6,
                  fmt = "%.2f", block_names = NULL, name_cex = 0.6, name_r = 1.14) {
  K <- length(values)
  width <- (2 * pi / K) * 0.88
  theta <- (seq_len(K) - 1) * 2 * pi / K
  for (k in seq_len(K)) {
    r_k <- r0 + (r1 - r0) * max(values[k], 0) / vmax
    .sector(theta[k] - width / 2, theta[k] + width / 2, r0, r_k, colors[k])
    if (label)
      graphics::text(r_k * cos(theta[k]) * 1.08, r_k * sin(theta[k]) * 1.08,
                     sprintf(fmt, values[k]), cex = cex)
  }
  if (!is.null(block_names))
    for (k in seq_len(K))
      graphics::text(r1 * name_r * cos(theta[k]), r1 * name_r * sin(theta[k]),
                     block_names[k], cex = name_cex, font = 2, col = colors[k])
}
