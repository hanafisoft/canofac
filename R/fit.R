#' Canonical factorisation of multiblock data in the q-family
#'
#' Estimates the reference components (scores) \eqn{U}, shared by all blocks, the axes
#' (loadings) \eqn{V_k} and the pseudo-singular values \eqn{d_{kh}} of the
#' canonical factorisation \eqn{X_k = U D_k V_k' + E_k} by maximising
#' \deqn{\sum_k \sum_h |d_{kh}|^q.}
#' The family includes MAXBET (\code{q = 1}), CPCA (\code{q = 2}, unconstrained
#' axes), MCOA (\code{q = 2}, orthonormal axes, sequential extraction), ACCPS /
#' CCSWA / ComDim (\code{q = 4}, unconstrained axes, sequential) and INDORT
#' (\code{q = 4}, unconstrained axes, simultaneous). STATIS and multiple factor
#' analysis (MFA) are CPCA (\code{q = 2}, unconstrained axes, simultaneous) on
#' weighted blocks: use \code{weighting = "statis"} or \code{weighting = "mfa"}.
#'
#' The blocks are used as given: preprocess them first with
#' \code{\link{cf_preprocess}} (by default, column centring and unit total
#' inertia).
#'
#' @param blocks list of numeric matrices with the same number of rows,
#'   already preprocessed.
#' @param R number of dimensions.
#' @param q power applied to the pseudo-singular values, \code{q >= 1}.
#' @param axes \code{"free"} (unconstrained axes of unit norm) or
#'   \code{"orthonormal"} (\eqn{V_k'V_k = I}).
#' @param extraction \code{"simultaneous"} or \code{"sequential"} (deflation).
#' @param tol relative tolerance on the criterion.
#' @param maxit maximum number of iterations.
#' @param nstart number of starts for unconstrained axes: the deterministic
#'   CPCA start plus \code{nstart - 1} random starts, the best solution being
#'   kept. For \code{q > 2} the criterion can be multimodal (several local
#'   maxima, one per block able to dominate a dimension), hence the default
#'   \code{NULL}: 10 starts if \code{q > 2} and \code{axes = "free"}, 1
#'   otherwise. The random starts are drawn by the common generator of
#'   canofac with the seed 12345 (\code{\link{cf_permutations}}), so the
#'   procedure is deterministic, gives the same starts as the Python version,
#'   and leaves R's global random number generator untouched. Ignored for
#'   orthonormal axes.
#' @param weighting weighting of the blocks before the fit: \code{"none"}
#'   (default, the blocks as given), \code{"mfa"} (each block divided by its
#'   first singular value, as in MFA) or \code{"statis"} (each block brought to
#'   \eqn{\|X_kX_k'\|_F = 1}, then its scalar product matrix weighted by the
#'   first eigenvector \eqn{\alpha} of the matrix \eqn{C_{kl} = tr(W_kW_l)}, as
#'   in STATIS). STATIS and MFA are obtained with \code{q = 2},
#'   \code{axes = "free"} and \code{extraction = "simultaneous"}.
#' @return An object of class \code{"canofac"}: a list with \code{U}
#'   (\eqn{n \times R}), \code{V} (list of \eqn{p_k \times R} matrices),
#'   \code{d} (\eqn{K \times R} pseudo-singular values), \code{saliences}
#'   (\eqn{d^2}), \code{delta} (\eqn{d_{kh}/\sigma_1(X_k)}), \code{criterion},
#'   \code{block_weights} (relative weights of the scalar product matrices of
#'   the blocks as given, summing to 1) and the settings. With a weighting,
#'   \code{d} and \code{saliences} are computed on the weighted blocks;
#'   \code{delta} does not depend on the weighting.
#' @references Bougeard S., Peltier C., Jaillais B., Boulet J.-C., Hanafi M.
#'   (2024). Benchmarking multiblock methods with canonical factorization.
#'   Chemometrics and Intelligent Laboratory Systems, 254, 105240.
#' @examples
#' set.seed(1)
#' u <- scale(rnorm(20))[, 1]
#' blocks <- lapply(c(5, 8, 12), function(p) outer(u, rnorm(p)) +
#'                    matrix(rnorm(20 * p), 20))
#' blocks <- cf_preprocess(blocks)
#' fit <- canofac(blocks, R = 2, q = 4)          # INDORT
#' fit
#' plot(fit)
#' plot(fit, type = "saliences")
#' statis <- canofac(blocks, R = 2, q = 2, weighting = "statis")
#' statis$block_weights
#' @export
canofac <- function(blocks, R = 2, q = 2,
                    axes = c("free", "orthonormal"),
                    extraction = c("simultaneous", "sequential"),
                    tol = 1e-10, maxit = 5000, nstart = NULL,
                    weighting = c("none", "mfa", "statis")) {
  axes <- match.arg(axes)
  extraction <- match.arg(extraction)
  weighting <- match.arg(weighting)
  nstart <- .nstart(nstart, q, axes)
  if (is.matrix(blocks)) blocks <- list(blocks)
  if (q < 1) stop("q must be >= 1")
  K <- length(blocks)
  if (is.null(names(blocks))) names(blocks) <- paste0("Block", seq_len(K))
  blocks <- lapply(blocks, function(X) { X <- as.matrix(X); storage.mode(X) <- "double"; X })
  n <- nrow(blocks[[1]])
  if (R > n - 1) stop("R must be smaller than n")
  wb <- .weight_blocks(blocks, weighting)
  blocks <- wb$blocks

  res <- .cf_fit(blocks, R, q, axes, extraction, tol, maxit, nstart)
  dimnames(res$d) <- list(names(blocks), paste0("Dim", seq_len(R)))
  names(res$V) <- names(blocks)
  rownames(res$U) <- rownames(blocks[[1]])
  for (k in seq_len(K)) {
    rownames(res$V[[k]]) <- colnames(blocks[[k]])
    colnames(res$V[[k]]) <- paste0("Dim", seq_len(R))
  }
  colnames(res$U) <- paste0("Dim", seq_len(R))
  sigma1 <- vapply(blocks, function(X) svd(X, nu = 0, nv = 0)$d[1], numeric(1))
  res$saliences <- res$d^2
  res$delta <- res$d / sigma1
  res$sigma1 <- sigma1
  res$criterion <- sum(abs(res$d)^q)
  res$block_weights <- wb$weights
  res$call <- match.call()
  res$settings <- list(R = R, q = q, axes = axes, extraction = extraction,
                       nstart = nstart, weighting = weighting)
  class(res) <- "canofac"
  res
}

# Core estimation (no checks) -------------------------------------------------
.nstart <- function(nstart, q, axes) {
  if (is.null(nstart)) nstart <- if (axes == "free" && q > 2) 10L else 1L
  max(1L, as.integer(nstart))
}

# Weighting of the blocks (STATIS, MFA). Returns the weighted blocks and the
# relative weights of the scalar product matrices of the blocks as given.
.wnorm <- function(X) {
  if (ncol(X) > nrow(X)) sqrt(sum(tcrossprod(X)^2)) else sqrt(sum(crossprod(X)^2))
}

.weight_blocks <- function(blocks, weighting = "none") {
  K <- length(blocks)
  if (weighting == "none" || K == 0)
    return(list(blocks = blocks, weights = stats::setNames(rep(1 / max(K, 1), K), names(blocks))))
  if (weighting == "mfa") {
    s1 <- vapply(blocks, function(X) svd(X, nu = 0, nv = 0)$d[1], numeric(1))
    out <- lapply(seq_len(K), function(k) blocks[[k]] / s1[k])
    w <- 1 / s1^2
  } else {
    nw <- vapply(blocks, .wnorm, numeric(1))
    Y <- lapply(seq_len(K), function(k) blocks[[k]] / sqrt(nw[k]))
    W <- lapply(Y, tcrossprod)
    C <- matrix(0, K, K)
    for (k in seq_len(K)) for (l in k:K) C[k, l] <- C[l, k] <- sum(W[[k]] * W[[l]])
    a <- abs(eigen(C, symmetric = TRUE)$vectors[, 1])
    out <- lapply(seq_len(K), function(k) sqrt(a[k]) * Y[[k]])
    w <- a / nw
  }
  names(out) <- names(blocks)
  list(blocks = out, weights = stats::setNames(w / sum(w), names(blocks)))
}

# Random orthonormal centred starts, drawn by the common generator of canofac
# (same starts in R and Python; R's global generator is not touched)
.starts <- function(n, R, nstart, seed = 12345) {
  if (nstart <= 1) return(list())
  g <- .gen_new(seed)
  lapply(seq_len(nstart - 1), function(i) {
    U <- matrix(.gen_normal(g, n * R), n, R)
    .polar(sweep(U, 2, colMeans(U)))
  })
}

# Unconstrained axes: CPCA start plus random starts, best criterion kept
.free_best <- function(blocks, R, q, tol, maxit, nstart = 1) {
  best <- .free_simultaneous(blocks, R, q, tol, maxit)
  cbest <- sum(best$d^q)
  for (U0 in .starts(nrow(blocks[[1]]), R, nstart)) {
    f <- .free_simultaneous(blocks, R, q, tol, maxit, U0 = U0)
    cf <- sum(f$d^q)
    if (cf > cbest * (1 + 1e-9)) { best <- f; cbest <- cf }
  }
  best
}

.cf_fit <- function(blocks, R, q, axes, extraction, tol = 1e-10, maxit = 5000,
                    nstart = 1) {
  K <- length(blocks)
  if (extraction == "sequential") {
    Xc <- blocks
    n <- nrow(blocks[[1]])
    U <- matrix(0, n, R)
    V <- lapply(blocks, function(X) matrix(0, ncol(X), R))
    D <- matrix(0, K, R)
    for (h in seq_len(R)) {
      u <- .free_best(Xc, 1, q, tol, maxit, nstart)$U[, 1]
      for (k in seq_len(K)) {
        a <- drop(crossprod(Xc[[k]], u))
        D[k, h] <- sqrt(sum(a^2))
        V[[k]][, h] <- if (D[k, h] > 0) a / D[k, h] else 0
      }
      U[, h] <- u
      for (k in seq_len(K)) {
        if (axes == "free") {
          Xc[[k]] <- Xc[[k]] - outer(u, drop(crossprod(u, Xc[[k]])))
        } else {
          v <- V[[k]][, h]
          Xc[[k]] <- Xc[[k]] - outer(drop(Xc[[k]] %*% v), v)
        }
      }
    }
    return(list(U = U, V = V, d = D))
  }

  if (axes == "free") {
    fit <- .free_best(blocks, R, q, tol, maxit, nstart)
    U <- fit$U
    D <- fit$d
    V <- lapply(seq_len(K), function(k) {
      Vk <- crossprod(blocks[[k]], U)
      sweep(Vk, 2, ifelse(D[k, ] > 0, D[k, ], 1), "/")
    })
  } else {
    start <- .cf_fit(blocks, R, 2, "orthonormal", "sequential", tol, maxit)  # MCOA
    fit <- .orth_simultaneous(blocks, R, q, tol, maxit, start$U, start$V)
    U <- fit$U
    V <- fit$V
    D <- fit$d
    if (q == 1) {                      # identification constraint
      M <- Reduce(`+`, lapply(seq_len(K), function(k)
        crossprod(U, blocks[[k]] %*% V[[k]])))
      M <- (M + t(M)) / 2
      Q <- eigen(M, symmetric = TRUE)$vectors
      U <- U %*% Q
      V <- lapply(V, function(Vk) Vk %*% Q)
      D <- .d_axes(blocks, U, V)
    }
    for (k in seq_len(K)) {            # positive pseudo-singular values
      sg <- ifelse(D[k, ] < 0, -1, 1)
      V[[k]] <- sweep(V[[k]], 2, sg, "*")
      D[k, ] <- D[k, ] * sg
    }
  }
  ord <- order(-colSums(abs(D)^q))
  list(U = U[, ord, drop = FALSE],
       V = lapply(V, function(Vk) Vk[, ord, drop = FALSE]),
       d = D[, ord, drop = FALSE])
}

.free_simultaneous <- function(blocks, R, q, tol, maxit, U0 = NULL, eps = 1e-12) {
  U <- if (is.null(U0)) svd(do.call(cbind, blocks), nu = R, nv = 0)$u else U0
  W <- lapply(blocks, tcrossprod)
  crit_old <- -Inf
  for (it in seq_len(maxit)) {
    D <- .d_free(blocks, U)
    crit <- sum(D^q)
    if (abs(crit - crit_old) <= tol * max(1, abs(crit))) break
    crit_old <- crit
    G <- matrix(0, nrow(U), R)
    for (k in seq_along(W))
      G <- G + sweep(W[[k]] %*% U, 2, pmax(D[k, ], eps)^(q - 2), "*")
    U <- .polar(G)
  }
  list(U = U, d = .d_free(blocks, U), iterations = it)
}

.orth_simultaneous <- function(blocks, R, q, tol, maxit, U0, V0) {
  U <- U0
  V <- V0
  K <- length(blocks)
  crit_old <- -Inf
  for (it in seq_len(maxit)) {
    D <- .d_axes(blocks, U, V)
    crit <- sum(abs(D)^q)
    if (abs(crit - crit_old) <= tol * max(1, abs(crit))) break
    crit_old <- crit
    Wt <- abs(D)^(q - 1) * sign(D)
    U <- .polar(Reduce(`+`, lapply(seq_len(K), function(k)
      sweep(blocks[[k]] %*% V[[k]], 2, Wt[k, ], "*"))))
    D <- .d_axes(blocks, U, V)
    Wt <- abs(D)^(q - 1) * sign(D)
    V <- lapply(seq_len(K), function(k)
      .polar(sweep(crossprod(blocks[[k]], U), 2, Wt[k, ], "*")))
  }
  list(U = U, V = V, d = .d_axes(blocks, U, V), iterations = it)
}
