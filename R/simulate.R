#' Simulated multiblock data of the article on the LOBO test
#'
#' Simulates \eqn{K} blocks with the scenarios of the simulation study of the
#' LOBO test. Two orthonormal centred directions \eqn{u_0}, \eqn{u_1} of the
#' individuals; in each block, noise with eigenvalues proportional to
#' \eqn{1/j} (\eqn{j = 1, \dots, p}) in a random orthonormal basis, plus a
#' signal: \code{"S0"} none (H0); \code{"SA"} \eqn{u_0} in all the blocks;
#' \code{"SB"} \eqn{u_0} in all the blocks and \eqn{u_1} in block 1 only;
#' \code{"SC"} \eqn{u_0} in blocks 1 and 2 only; \code{"SD"} \eqn{u_0} in all
#' the blocks and \eqn{u_1} (intensity \code{theta2}) in blocks 1 and 2;
#' \code{"SE"} \eqn{u_0} in blocks 1 and 2 and \eqn{u_1} (\code{theta2}) in all
#' the blocks. Intensities are multiplied by \eqn{\sqrt{n}}.
#'
#' All the draws come from the common generator of canofac and the orthonormal
#' bases are polar factors (no QR sign convention): the same seed gives the
#' same blocks in R and in Python.
#'
#' @param scenario one of \code{"S0"}, \code{"SA"}, \code{"SB"}, \code{"SC"},
#'   \code{"SD"}, \code{"SE"}.
#' @param n number of individuals.
#' @param p number of variables per block.
#' @param K number of blocks.
#' @param theta intensity of the first direction.
#' @param theta2 intensity of the second direction (SD, SE).
#' @param seed seed of the common generator.
#' @param preprocess logical: centre and bring each block to unit total inertia.
#' @return A named list of \eqn{K} matrices.
#' @examples
#' bl <- cf_simulate("SC", n = 20, p = 10, K = 4, theta = 1.8, seed = 3)
#' cf_test(bl, R = 2, B = 99, seed = 1)
#' @export
cf_simulate <- function(scenario = c("SA", "S0", "SB", "SC", "SD", "SE"), n = 20, p = 10, K = 4,
                        theta = 1.6, theta2 = 1.5, seed = 1, preprocess = TRUE) {
  scenario <- match.arg(scenario)
  g <- .gen_new(seed)
  G <- matrix(.gen_normal(g, n * 2), n, 2)
  U <- .polar(sweep(G, 2, colMeans(G)))
  u0 <- U[, 1]; u1 <- U[, 2]
  lam <- 1 / seq_len(p)
  lam <- lam * p / sum(lam)
  t1 <- theta * sqrt(n); t2 <- theta2 * sqrt(n)
  out <- vector("list", K)
  for (k in seq_len(K)) {
    Q <- .polar(matrix(.gen_normal(g, p * p), p, p))
    E <- sweep(matrix(.gen_normal(g, n * p), n, p), 2, sqrt(lam), "*") %*% t(Q)
    S <- matrix(0, n, p)
    a <- .gen_normal(g, p); a <- a / sqrt(sum(a^2))
    if (scenario %in% c("SA", "SB", "SD") || (scenario %in% c("SC", "SE") && k <= 2))
      S <- S + t1 * outer(u0, a)
    if (scenario %in% c("SB", "SD", "SE")) {
      b <- .gen_normal(g, p); b <- b / sqrt(sum(b^2))
      if (scenario == "SB" && k == 1) S <- S + t1 * outer(u1, b)
      if (scenario == "SD" && k <= 2) S <- S + t2 * outer(u1, b)
      if (scenario == "SE") S <- S + t2 * outer(u1, b)
    }
    out[[k]] <- S + E
  }
  names(out) <- paste0("Block", seq_len(K))
  if (preprocess) out <- cf_preprocess(out)
  out
}
