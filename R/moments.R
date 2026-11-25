#' Exact permutation moments and gamma approximation
#'
#' For doubly centred symmetric \eqn{n \times n} matrices \eqn{A} and \eqn{B}
#' and a uniformly random permutation matrix \eqn{\Pi}, the mean and the
#' variance of \eqn{T = tr(A \Pi B \Pi')} are exact (Kazi-Aoual et al., 1995;
#' Josse et al., 2008):
#' \eqn{E(T) = tr(A) tr(B) / (n - 1)}, and \eqn{Var(T)} is an explicit
#' function of \eqn{n} and, for \eqn{A} and \eqn{B}, of the trace, the trace
#' of the square and the sum of the squared diagonal terms. All the statistics
#' of the LOBO test have this form with \eqn{A = W_k}: \eqn{M = u_h u_h'} for
#' \eqn{S_{kh}^2}, \eqn{P} for \eqn{T_k^2}, \eqn{P W_{(-k)} P} for the
#' restricted RV. In particular \eqn{E(S_{kh}^2) = 1/(n - 1)} for a block of
#' unit inertia, the chance level of a salience.
#'
#' \code{gamma_pvalue} returns the upper-tail p-value of \eqn{T} under the
#' gamma law with these two moments, used by \code{cf_test(pvalue = "gamma")}
#' and \code{"hybrid"}.
#'
#' @param A,B doubly centred symmetric matrices (for example \eqn{W_k} and
#'   \eqn{U U'}).
#' @param obs observed value of \eqn{tr(A B)}.
#' @return \code{perm_moments}: a vector \code{c(mean, var)};
#'   \code{gamma_pvalue}: a p-value.
#' @references Kazi-Aoual F., Hitier S., Sabatier R., Lebreton J.-D. (1995).
#'   Refined approximations to permutation tests for multivariate inference.
#'   Computational Statistics and Data Analysis, 20, 643-656.
#'
#'   Josse J., Pagès J., Husson F. (2008). Testing the significance of the RV
#'   coefficient. Computational Statistics and Data Analysis, 53, 82-91.
#' @seealso \code{\link{cf_test}}; French help: \code{cf_help("approximation", "fr")}.
#' @examples
#' X <- scale(matrix(rnorm(60), 20), scale = FALSE)
#' u <- svd(scale(matrix(rnorm(80), 20), scale = FALSE))$u[, 1:2]
#' perm_moments(tcrossprod(X), tcrossprod(u))
#' @export
perm_moments <- function(A, B) {
  n <- nrow(A)
  s <- function(M) {
    t <- sum(diag(M)); t2 <- sum(M * M); s2 <- sum(diag(M)^2)
    c(t, t2, t * t / t2, s2 / t2)
  }
  a <- s(A); b <- s(B)
  m <- a[1] * b[1] / (n - 1)
  ca <- n * (n + 1) * a[4] - (n - 1) * (a[3] + 2)
  cb <- n * (n + 1) * b[4] - (n - 1) * (b[3] + 2)
  v <- 2 * (n - 1 - a[3]) * (n - 1 - b[3]) / ((n + 1) * (n - 1)^2 * (n - 2)) *
    (1 + (n - 3) / (2 * n * (n - 1)) * ca * cb / ((n - 1 - a[3]) * (n - 1 - b[3])))
  c(mean = m, var = v * a[2] * b[2])
}

#' @rdname perm_moments
#' @export
gamma_pvalue <- function(obs, A, B) {
  if (sum(diag(B)) <= 0) return(1)
  mv <- perm_moments(A, B)
  m <- mv[[1]]; v <- mv[[2]]
  if (!is.finite(v) || v <= 0 || m <= 0) return(1)
  stats::pgamma(obs, shape = m * m / v, scale = v / m, lower.tail = FALSE)
}
