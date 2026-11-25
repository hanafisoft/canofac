# Internal helpers -----------------------------------------------------------

# Polar factor of G (closest matrix with orthonormal columns).
.polar <- function(G) {
  s <- svd(G)
  s$u %*% t(s$v)
}

# Principal coordinates of a block: same W_k = X X', cost independent of p_k.
.compress <- function(X) {
  s <- svd(X)
  r <- sum(s$d > s$d[1] * 1e-12)
  s$u[, seq_len(r), drop = FALSE] %*% diag(s$d[seq_len(r)], r)
}

# Pseudo-singular values for unconstrained axes: d_kh = ||X_k' u_h||.
.d_free <- function(blocks, U) {
  matrix(unlist(lapply(blocks, function(X) sqrt(colSums(crossprod(X, U)^2)))),
         nrow = length(blocks), byrow = TRUE)
}

# Pseudo-singular values for given axes: d_kh = u_h' X_k v_kh.
.d_axes <- function(blocks, U, V) {
  matrix(unlist(lapply(seq_along(blocks), function(k)
    colSums(U * (blocks[[k]] %*% V[[k]])))),
    nrow = length(blocks), byrow = TRUE)
}

# Optimal assignment maximizing the total similarity (Hungarian algorithm,
# O(R^3)). C is a square similarity matrix; returns a, a permutation such that
# row i is assigned to column a[i] and sum(C[cbind(i, a)]) is maximal.
# Arrays u, v, p, way, minv, used are indexed from 0 in the algorithm and
# stored with an offset of 1.
.lsap_max <- function(C) {
  n <- nrow(C)
  if (n == 1L) return(1L)
  a <- max(C) - C
  u <- numeric(n + 1L); v <- numeric(n + 1L)
  p <- integer(n + 1L); way <- integer(n + 1L)
  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 0L
    minv <- rep(Inf, n + 1L)
    used <- rep(FALSE, n + 1L)
    repeat {
      used[j0 + 1L] <- TRUE
      i0 <- p[j0 + 1L]
      delta <- Inf
      j1 <- 0L
      for (j in seq_len(n)) {
        if (!used[j + 1L]) {
          cur <- a[i0, j] - u[i0 + 1L] - v[j + 1L]
          if (cur < minv[j + 1L]) {
            minv[j + 1L] <- cur
            way[j + 1L] <- j0
          }
          if (minv[j + 1L] < delta) {
            delta <- minv[j + 1L]
            j1 <- j
          }
        }
      }
      for (j in 0:n) {
        if (used[j + 1L]) {
          u[p[j + 1L] + 1L] <- u[p[j + 1L] + 1L] + delta
          v[j + 1L] <- v[j + 1L] - delta
        } else {
          minv[j + 1L] <- minv[j + 1L] - delta
        }
      }
      j0 <- j1
      if (p[j0 + 1L] == 0L) break
    }
    repeat {
      j1 <- way[j0 + 1L]
      p[j0 + 1L] <- p[j1 + 1L]
      j0 <- j1
      if (j0 == 0L) break
    }
  }
  ans <- integer(n)
  for (j in seq_len(n)) ans[p[j + 1L]] <- j
  ans
}

# Match the dimensions of U to those of Uref.
# Returns pi such that dimension pi[j] of U is associated with dimension j of
# Uref, maximizing sum_j |Uref[, j]' U[, pi[j]]| (optimal assignment, exact
# for any number of dimensions).
.match_dims <- function(U, Uref) {
  C <- abs(crossprod(Uref, U))
  .lsap_max(C)
}

# Effective rank D2 = (tr W)^2 / tr(W^2).
.effective_rank <- function(X) {
  l <- svd(X, nu = 0, nv = 0)$d^2
  sum(l)^2 / sum(l^2)
}
