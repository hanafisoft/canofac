#' Common random number generator of canofac (identical in R and Python)
#'
#' All the random draws of canofac (permutations of \code{\link{cf_test}},
#' random starts of \code{\link{canofac}}) come from a generator of the
#' package, the same in R and in Python: the Park-Miller "minimal standard"
#' generator (multiplier 48271, modulus \eqn{2^{31} - 1}). All products stay
#' below \eqn{2^{53}} and are exact in double precision, so that R and Python
#' produce the same sequence for the same seed. Permutations are drawn by the
#' Fisher-Yates algorithm; permutation \eqn{b} uses the draws
#' \eqn{b(n - 1) + 1, \dots, (b + 1)(n - 1)} of the stream (skip-ahead), so
#' that the result does not depend on how the computation is vectorised. The
#' global random number generator of R (\code{set.seed}) is neither used nor
#' modified, except to draw a seed when none is given.
#'
#' @param B number of permutations.
#' @param n size of each permutation.
#' @param seed integer seed.
#' @return \code{cf_permutations}: a \eqn{B \times n} integer matrix, one
#'   permutation of \code{1:n} per row (Python returns \code{0:(n-1)}).
#' @seealso \code{\link{cf_test}}; French help: \code{cf_help("generateur", "fr")}.
#' @examples
#' cf_permutations(2, 5, seed = 3)
#' @export
cf_permutations <- function(B, n, seed = 1) {
  .gen_perms(.gen_new(seed), B, n)
}

.MINSTD_M <- 2147483647
.MINSTD_A <- 48271

# (a * b) mod M for doubles a, b in [0, M), exact: b is split into 16-bit halves.
.mulmod <- function(a, b) {
  bh <- floor(b / 65536)
  bl <- b - bh * 65536
  ((((a * bh) %% .MINSTD_M) * 65536) %% .MINSTD_M + a * bl) %% .MINSTD_M
}

.powmod <- function(a, e) {
  r <- 1
  a <- a %% .MINSTD_M
  while (e > 0) {
    if (e %% 2 == 1) r <- .mulmod(r, a)
    a <- .mulmod(a, a)
    e <- e %/% 2
  }
  r
}

# A generator is an environment holding the state (a double in [1, M - 1]).
.gen_new <- function(seed) {
  g <- new.env(parent = emptyenv())
  g$state <- (as.numeric(seed) %% (.MINSTD_M - 1)) + 1
  g
}

.gen_uniform <- function(g, size) {
  out <- numeric(size)
  x <- g$state
  for (i in seq_len(size)) {
    x <- (x * .MINSTD_A) %% .MINSTD_M
    out[i] <- x / .MINSTD_M
  }
  g$state <- x
  out
}

.gen_normal <- function(g, size) {
  m <- (size + 1) %/% 2
  u <- matrix(.gen_uniform(g, 2 * m), ncol = 2, byrow = TRUE)
  r <- sqrt(-2 * log(u[, 1]))
  z <- as.vector(t(cbind(r * cos(2 * pi * u[, 2]), r * sin(2 * pi * u[, 2]))))
  z[seq_len(size)]
}

# B permutations of 1:n (B x n matrix), from the current state of g.
.gen_perms <- function(g, B, n) {
  B <- as.integer(B); n <- as.integer(n)
  perm <- matrix(rep(seq_len(n), each = B), B, n)
  if (n < 2 || B < 1) return(perm)
  step <- .powmod(.MINSTD_A, n - 1)
  x <- numeric(B)
  s <- g$state
  for (b in seq_len(B)) {
    x[b] <- s
    s <- .mulmod(s, step)
  }
  rows <- seq_len(B)
  for (i in (n - 1):1) {                    # 0-based position i, column i + 1
    x <- (x * .MINSTD_A) %% .MINSTD_M
    j <- pmin(floor(x / .MINSTD_M * (i + 1)), i)
    tmp <- perm[cbind(rows, i + 1)]
    perm[cbind(rows, i + 1)] <- perm[cbind(rows, j + 1)]
    perm[cbind(rows, j + 1)] <- tmp
  }
  g$state <- s
  perm
}

# Seed to use: the one given, or one drawn by R's generator.
.draw_seed <- function(seed) {
  if (!is.null(seed)) return(as.numeric(seed))
  as.numeric(sample.int(.Machine$integer.max - 1L, 1L))
}
