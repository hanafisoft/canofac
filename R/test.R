#' Permutation test of the links between each block and the reference components
#'
#' Tests, block by block, whether each block is linked to the reference
#' components of the other blocks beyond what a block independent of them
#' (\code{scheme = "rows"}) or without structure (\code{scheme = "variables"})
#' would produce, and designates the blocks linked to the subspace of the
#' reference components and to each of its components: all the blocks,
#' several, one or none.
#'
#' \strong{LOBO statistic} (\code{statistic = "lobo"}, default; leave-one-block-out, \code{"fixed"} in canofac < 0.6.0).
#' For each block \eqn{k}, the reference components \eqn{U^{(-k)}} are
#' estimated once, by the same procedure, on the \eqn{K - 1} other blocks. The
#' statistics are \eqn{S_{kh} = \|X_k' u_h^{(-k)}\|} for each component and
#' \eqn{T_k = \|X_k' U^{(-k)}\|_F} for the subspace. Only block \eqn{k} is
#' permuted, and nothing is re-estimated: the test is exact under the null
#' hypothesis of the scheme, needs no matching of the dimensions and is fast.
#' All the blocks are linked to the subspace when \eqn{\max_k p_k \le \alpha}
#' (intersection-union test, \code{link_subspace = "all"}); otherwise the
#' linked blocks are those whose subspace p-value, adjusted by Holm's method
#' over the \eqn{K} blocks, is at most \eqn{\alpha}: \code{"several"},
#' \code{"one"} or \code{"none"}. The p-values \eqn{p_{kh}} bear on the
#' components of the other blocks, which differ from one tested block to
#' another and whose ordering is not stable; they and the linked blocks derived
#' from them dimension by dimension are descriptive. The
#' \eqn{h}th component of the other blocks is \emph{linked} (\code{retained})
#' when \eqn{\min_k p_{kh} \le \alpha / K}. This statistic was proposed by
#' J.-M. Galharret.
#'
#' \strong{Subspace statistic} (\code{subspace}). \code{"inertia"} (default,
#' restricted RV) uses
#' \eqn{T_k^w = tr(W_k P^{(-k)} W_{(-k)} P^{(-k)})}, where
#' \eqn{P^{(-k)} = U^{(-k)}U^{(-k)'}} and \eqn{W_{(-k)}} is the sum of the
#' scalar product matrices of the other blocks, weighted as in the method: it
#' weights the directions of the subspace by the inertia of the other blocks
#' (for CPCA, MFA and STATIS, \eqn{T_k^w = \sum_h \lambda_h^{(-k)} S_{kh}^2}) and
#' is the numerator of the RV coefficient between \eqn{X_k} and the other
#' blocks projected on \eqn{U^{(-k)}} (restricted RV coefficient).
#' \code{"equal"} uses \eqn{T_k^2 = \sum_h S_{kh}^2}, which gives the \eqn{R}
#' components the same weight. Both statistics use the same subspace, null
#' hypothesis and permutations, and are exact; \code{"inertia"} is at least as
#' powerful, and clearly more so for small samples. Both p-values are returned
#' (\code{p_subspace_equal}, \code{p_subspace_inertia}); \code{p_subspace} and
#' the linked blocks use the chosen one. The dimension-wise statistics
#' \eqn{S_{kh}} are the same for both.
#'
#' \strong{Weighting of the blocks} (\code{weighting}). With \code{"statis"} or
#' \code{"mfa"}, the components of the other blocks are those of STATIS or MFA
#' applied to the \eqn{K - 1} other blocks (STATIS weights recomputed on these
#' blocks); use \code{q = 2}, \code{axes = "free"} and
#' \code{extraction = "simultaneous"}. The p-values do not depend on the scale of
#' the tested block.
#'
#' \strong{Links are not labels.} Being linked to a component is not carrying
#' a structure common to the blocks: a component may mix directions that the
#' blocks carry differently, and a block may be linked to it through a small
#' share only; two directions carried by two blocks, distinct but correlated
#' through the individuals, also link these blocks. Labeling parts of the
#' retained space as common, shared, private or residual requires per-block
#' thresholds and a change of basis (CoPRes,
#' \url{https://github.com/hanafisoft/CoPRes}).
#'
#' \strong{Re-estimation statistic} (\code{statistic = "refit"}). The
#' pseudo-singular value \eqn{d_{kh}} of the full analysis, recomputed by
#' re-estimating the whole model on each permuted data set (behaviour of
#' canofac < 0.4.0). The test is exact under the global null hypothesis of the
#' scheme, but not dimension by dimension: a permuted block keeps the
#' structures that linked it to the other blocks, now carried by this block
#' alone, which attract a dimension of the re-estimated model and hide a
#' second link. Use it as a global test (subspace) or, with
#' \code{scheme = "variables"}, to recognise directions carried by one block
#' without link to the others. With \code{matching = TRUE}, the
#' dimensions are matched to the solution on the other blocks; the p-value of
#' observed dimension 1 is then slightly liberal and that of the last
#' dimension slightly conservative.
#'
#' \strong{Permutation of the principal components} (\code{scheme =
#' "components"}). The scores of each principal component of the tested block
#' are permuted independently, and the block is rebuilt with its loadings.
#' Unlike variable permutation, the scheme is invariant under rotation of the
#' variables and preserves the spectrum and the loadings of the block (for
#' spectra, the continuity of each loading); it randomizes the orientation of
#' each component. Its null hypothesis, principal components independent of
#' each other and of the other blocks, is harder to interpret; it is an
#' alternative to variable permutation for blocks whose variables are ordered
#' and correlated, such as spectra, when the question is whether the block
#' carries a signal of its own.
#'
#' The blocks must be preprocessed with \code{\link{cf_preprocess}} (column
#' centring and unit total inertia, which both schemes leave unchanged). With
#' \code{scheme = "rows"}, blocks with more variables than rows are replaced by
#' their principal coordinates, which changes neither the estimation nor the
#' test.
#'
#' \strong{Exchangeability groups} (\code{scheme = "groups"} or
#' \code{"within"}, argument \code{groups}). With replicates, batches or other
#' groups of rows, free permutation of the rows is liberal; \code{"groups"}
#' permutes whole groups (of equal sizes), \code{"within"} permutes the rows
#' inside each group.
#'
#' \strong{p-values without permutation} (\code{pvalue}). The statistics are
#' \eqn{T = tr(W_k \Pi M \Pi')} with \eqn{M} fixed, whose permutation mean and
#' variance are exact (Kazi-Aoual et al., 1995; Josse et al., 2008; see
#' \code{\link{perm_moments}}). \code{"gamma"} fits a gamma law with these
#' two moments, without permutation: in the simulations of the companion
#' article, it gives the same decisions as the permutation test in 99\% of the
#' cases but is slightly liberal (5.2\% rejections under H0 instead of
#' 4.8\%). \code{"hybrid"} replaces the gamma p-value by the permutation
#' p-value when it lies between \eqn{\alpha / (c K)} and \eqn{c \alpha}
#' (\code{hybrid_c} = 3): same decisions as the permutation test while
#' permuting about one test in ten, and very small p-values beyond
#' \eqn{1/(B + 1)}. Both need \code{statistic = "lobo"} and
#' \code{scheme = "rows"}; \code{pmethod} says which method gave each p-value.
#'
#' \strong{Comparison tests} (\code{compare = TRUE}). The full RV test (block
#' against the concatenation of the other blocks) and PROTEST (Procrustes
#' statistic between the first \eqn{R} principal components of the block and
#' of the other blocks) are computed with the same permutations. They are
#' comparisons, not part of the LOBO test: the full RV test is more powerful on
#' average, but it says whether the block is linked to the other blocks, not to
#' the retained components.
#'
#' \strong{Random draws.} The permutations come from the common generator of
#' canofac (\code{\link{cf_permutations}}): with the same \code{seed}, R and
#' Python give the same permutations, hence the same p-values. Order of the
#' draws: block 1, then block 2, ..., \code{B} permutations per block (for
#' \code{"variables"}, \code{B} x \eqn{p_k} permutations, variable by
#' variable; for \code{"within"}, group by group).
#'
#' @inheritParams canofac
#' @param scheme \code{"rows"} (permute the rows of the tested block),
#'   \code{"groups"} (permute whole groups of rows, of equal sizes),
#'   \code{"within"} (permute the rows within each group),
#'   \code{"variables"} (permute each variable of the tested block
#'   independently) or \code{"components"} (permute the scores of each
#'   principal component of the tested block independently).
#' @param groups group of each row, for \code{scheme = "groups"} or
#'   \code{"within"}.
#' @param B number of permutations.
#' @param alpha level.
#' @param statistic \code{"lobo"} (LOBO statistic, leave-one-block-out, default; \code{"fixed"} is an alias) or
#'   \code{"refit"} (re-estimation statistic).
#' @param matching logical, for \code{statistic = "refit"}: match the
#'   dimensions to the solution on the other blocks.
#' @param seed seed of the common generator (same permutations in R and
#'   Python); if \code{NULL}, a seed is drawn by R's generator and kept in
#'   \code{settings$seed}.
#' @param subspace subspace statistic of the LOBO test:
#'   \code{"inertia"} (\eqn{T_k^w}, restricted RV, weights by the inertia of
#'   the other blocks, default) or \code{"equal"} (\eqn{T_k}, equal weights).
#' @param pvalue \code{"permutation"} (default), \code{"gamma"} or
#'   \code{"hybrid"}.
#' @param hybrid_c width of the zone of the hybrid mode (3).
#' @param compare logical: also compute the full RV test and PROTEST.
#' @param weighting weighting of the blocks, \code{"none"}, \code{"mfa"} or
#'   \code{"statis"} (see \code{\link{canofac}}).
#' @return An object of class \code{"cf_test"}: \code{p} (\eqn{K \times R}
#'   p-values; components of the other blocks for \code{"lobo"}, observed
#'   dimensions for \code{"refit"}), \code{p_holm} (adjusted over the blocks),
#'   \code{linked} (\eqn{K \times R} logical, Holm-adjusted p-value at most
#'   \eqn{\alpha}) and \code{link} (\code{"all"}, \code{"several"},
#'   \code{"one"} or \code{"none"}, dimension by dimension),
#'   \code{p_subspace}, \code{linked_subspace} and \code{link_subspace} (the
#'   same for the subspace, with the chosen subspace statistic),
#'   \code{p_subspace_equal} and \code{p_subspace_inertia} (both subspace
#'   p-values, LOBO statistic only), \code{all_linked} (logical, one per dimension:
#'   \eqn{\max_k p_{kh} \le \alpha}), \code{retained} (logical, one per
#'   dimension: linked component for \code{"lobo"}; at least one linked block
#'   for \code{"refit"}), \code{refs} (the \eqn{U^{(-k)}}), \code{fit} (the observed
#'   \code{canofac} fit), \code{p_compare} (full RV and PROTEST p-values, if
#'   \code{compare = TRUE}), \code{pmethod} (\code{"perm"} or \code{"gamma"}
#'   for each p-value, components then subspace), \code{null} (permutation
#'   distributions of the subspace statistic, for \code{plot(x, type = "null")})
#'   and the settings, including the seed. French help:
#'   \code{cf_help("cf_test", "fr")}.
#' @examples
#' set.seed(2)
#' u <- scale(rnorm(20))[, 1]
#' blocks <- cf_preprocess(lapply(c(5, 8, 12), function(p)
#'   outer(u, rnorm(p)) + matrix(rnorm(20 * p), 20)))
#' tt <- cf_test(blocks, R = 2, q = 2, B = 99, seed = 1)
#' tt
#' plot(tt)
#' cbind(inertia = tt$p_subspace_inertia, equal = tt$p_subspace_equal)
#' @export
cf_test <- function(blocks, R = 2, q = 2,
                    axes = c("free", "orthonormal"),
                    extraction = c("simultaneous", "sequential"),
                    scheme = c("rows", "groups", "within", "variables", "components"),
                    groups = NULL,
                    B = 999, alpha = 0.05, statistic = c("lobo", "refit", "fixed"),
                    matching = TRUE, seed = NULL,
                    tol = 1e-10, maxit = 5000, nstart = NULL,
                    subspace = c("inertia", "equal"),
                    weighting = c("none", "mfa", "statis"),
                    pvalue = c("permutation", "gamma", "hybrid"), hybrid_c = 3,
                    compare = FALSE) {
  axes <- match.arg(axes)
  extraction <- match.arg(extraction)
  scheme <- match.arg(scheme)
  statistic <- match.arg(statistic)
  if (statistic == "fixed") statistic <- "lobo"   # alias (canofac < 0.6.0)
  subspace <- match.arg(subspace)
  weighting <- match.arg(weighting)
  pvalue <- match.arg(pvalue)
  if (pvalue != "permutation" && (statistic != "lobo" || scheme != "rows"))
    stop("pvalue = 'gamma' or 'hybrid' needs statistic = 'lobo' and scheme = 'rows'")
  if (scheme %in% c("groups", "within") && is.null(groups))
    stop("scheme = 'groups' or 'within' needs groups")
  if (compare && scheme %in% c("variables", "components"))
    stop("compare = TRUE needs a row-type scheme ('rows', 'groups', 'within')")
  # the restricted RV statistic is defined for the LOBO statistic only
  if (statistic == "refit") subspace <- "equal"
  nstart <- .nstart(nstart, q, axes)
  seed <- .draw_seed(seed)
  gen <- .gen_new(seed)
  K <- length(blocks)
  if (is.null(names(blocks))) names(blocks) <- paste0("Block", seq_len(K))
  blocks <- lapply(blocks, function(X) { X <- as.matrix(X); storage.mode(X) <- "double"; X })
  n <- nrow(blocks[[1]])
  if (pvalue != "permutation" && !all(vapply(blocks, .is_centred, logical(1))))
    stop("the gamma approximation needs column-centred blocks (see cf_preprocess)")
  work <- blocks
  if (scheme %in% c("rows", "groups", "within", "components"))
    work <- lapply(blocks, function(X) if (ncol(X) > n - 1) .compress(X) else X)
  pcs <- NULL
  if (scheme == "components")
    pcs <- lapply(work, function(X) {
      s <- svd(X)
      list(scores = sweep(s$u, 2, s$d, "*"), loadings = s$v)
    })
  gidx <- NULL
  if (scheme %in% c("groups", "within")) {
    if (length(groups) != n) stop("groups must have one value per row")
    lev <- unique(as.character(groups))                 # order of first appearance
    gidx <- lapply(lev, function(l) which(as.character(groups) == l))
    if (scheme == "groups" && length(unique(lengths(gidx))) != 1)
      stop("scheme = 'groups' needs groups of equal sizes")
  }
  # Row permutations (B x n) for the row-type schemes, drawn in the same order as in Python
  row_perms <- function(Bk) {
    if (scheme == "rows") return(.gen_perms(gen, Bk, n))
    perm <- matrix(rep(seq_len(n), each = Bk), Bk, n)
    if (scheme == "groups") {
      gp <- .gen_perms(gen, Bk, length(gidx))
      G <- do.call(rbind, gidx)                          # groups x size
      for (g in seq_along(gidx)) perm[, gidx[[g]]] <- G[gp[, g], , drop = FALSE]
      return(perm)
    }
    for (ig in gidx) {                                   # within
      pg <- .gen_perms(gen, Bk, length(ig))
      perm[, ig] <- matrix(ig[pg], Bk)
    }
    perm
  }
  # Permuted versions of block k for the other schemes (list of Bk matrices)
  other_perms <- function(k, Bk) {
    if (scheme == "components") {
      sc <- pcs[[k]]$scores
      r <- ncol(sc)
      P <- .gen_perms(gen, Bk * r, n)
      return(lapply(seq_len(Bk), function(b) {
        S <- vapply(seq_len(r), function(j) sc[P[(b - 1) * r + j, ], j], numeric(n))
        tcrossprod(matrix(S, n), pcs[[k]]$loadings)
      }))
    }
    X <- blocks[[k]]
    p <- ncol(X)
    P <- .gen_perms(gen, Bk * p, n)
    lapply(seq_len(Bk), function(b) {
      Xp <- vapply(seq_len(p), function(j) X[P[(b - 1) * p + j, ], j], numeric(n))
      cf_preprocess(list(matrix(Xp, n)))[[1]]
    })
  }
  match_on <- matching && K > 1

  substat <- function(D, k) {
    if (axes == "orthonormal" && q == 1) sum(D[k, ]) else sum(D[k, ]^2)
  }
  M <- function(bl) .cf_fit(.weight_blocks(bl, weighting)$blocks, R, q, axes,
                            extraction, tol, maxit, nstart)

  P <- matrix(1, K, R, dimnames = list(names(blocks), paste0("Dim", seq_len(R))))
  Psub <- stats::setNames(rep(1, K), names(blocks))
  Psub_eq <- Psub_in <- NULL
  pmethod <- matrix("perm", K, R + 1, dimnames = list(names(blocks), c(paste0("Comp", seq_len(R)), "Subspace")))
  null <- if (pvalue != "gamma") list() else NULL
  refs <- vector("list", K)
  names(refs) <- names(blocks)
  comp_p <- if (compare) matrix(1, K, 2, dimnames = list(names(blocks), c("RV", "PROTEST"))) else NULL
  if (statistic == "lobo") {
    Psub_eq <- Psub_in <- Psub
    for (k in seq_len(K)) {
      Xk <- work[[k]]
      if (K > 1) {
        others <- .weight_blocks(work[-k], weighting)$blocks
        Uref <- .cf_fit(others, R, q, axes, extraction, tol, maxit, nstart)$U
        # U' W_(-k) U: inertia of the other blocks in the subspace
        Cw <- Reduce(`+`, lapply(others, function(X) tcrossprod(crossprod(Uref, X))))
      } else {
        Uref <- matrix(0, n, R)
        Cw <- matrix(0, R, R)
      }
      refs[[k]] <- Uref
      A <- crossprod(Xk, Uref)
      AtA <- crossprod(A)
      s_obs <- sqrt(diag(AtA))
      t_obs <- sqrt(sum(s_obs^2))
      tw_obs <- sum(AtA * Cw)
      Xr <- if (compare && K > 1) do.call(cbind, work[-k]) else NULL
      need <- rep(TRUE, R + 1)
      pg <- NULL
      if (pvalue != "permutation") {
        Wk <- tcrossprod(Xk)
        pg <- c(vapply(seq_len(R), function(h) gamma_pvalue(s_obs[h]^2, Wk, tcrossprod(Uref[, h])), numeric(1)),
                gamma_pvalue(t_obs^2, Wk, tcrossprod(Uref)),
                gamma_pvalue(tw_obs, Wk, Uref %*% Cw %*% t(Uref)))
        lo <- alpha / (hybrid_c * K); hi <- hybrid_c * alpha
        if (pvalue == "gamma") {
          need <- rep(FALSE, R + 1)
        } else {
          zone <- pg >= lo & pg <= hi
          need <- c(zone[seq_len(R)], zone[if (subspace == "equal") R + 1 else R + 2])
          need_eq_in <- zone[R + 1:2]
        }
      }
      run_perm <- pvalue == "permutation" || any(need) || (compare && pvalue != "gamma")
      if (run_perm) {
        if (scheme %in% c("rows", "groups", "within")) {
          Pr <- row_perms(B)
          Ap <- lapply(seq_len(B), function(b) crossprod(Xk[Pr[b, ], , drop = FALSE], Uref))
        } else {
          Pr <- NULL
          Ap <- lapply(other_perms(k, B), function(Xp) crossprod(Xp, Uref))
        }
        Mp <- lapply(Ap, crossprod)
        s <- t(vapply(Mp, function(m) sqrt(diag(m)), numeric(R)))
        if (R == 1) s <- matrix(s, ncol = 1)
        tt <- sqrt(rowSums(s^2))
        tw <- vapply(Mp, function(m) sum(m * Cw), numeric(1))
        p_comp <- (1 + colSums(sweep(s, 2, s_obs - 1e-12, ">="))) / (B + 1)
        p_eq <- (1 + sum(tt >= t_obs - 1e-12)) / (B + 1)
        p_in <- (1 + sum(tw >= tw_obs - 1e-12 * abs(tw_obs))) / (B + 1)
        null[[names(blocks)[k]]] <- list(observed = if (subspace == "inertia") tw_obs else t_obs,
                                         permuted = if (subspace == "inertia") tw else tt)
        if (compare && !is.null(Pr) && K > 1) {
          Wk <- tcrossprod(Xk)
          Wr <- tcrossprod(Xr)
          rv_obs <- sum(Wk * Wr)
          rv_null <- vapply(seq_len(B), function(b) sum(Wk[Pr[b, ], Pr[b, ]] * Wr), numeric(1))
          Zk <- .pcs(Xk, R); Zr <- .pcs(Xr, R)
          pr_obs <- sum(svd(crossprod(Zk, Zr), nu = 0, nv = 0)$d)
          pr_null <- vapply(seq_len(B), function(b)
            sum(svd(crossprod(Zk[Pr[b, ], , drop = FALSE], Zr), nu = 0, nv = 0)$d), numeric(1))
          comp_p[k, ] <- c((1 + sum(rv_null >= rv_obs * (1 - 1e-12))) / (B + 1),
                           (1 + sum(pr_null >= pr_obs - 1e-12)) / (B + 1))
        }
      }
      if (compare && pvalue == "gamma" && K > 1) {
        Wk <- tcrossprod(Xk); Wr <- tcrossprod(Xr)
        comp_p[k, ] <- c(gamma_pvalue(sum(Wk * Wr), Wk, Wr), NA)
      }
      if (pvalue == "permutation") {
        P[k, ] <- p_comp
        Psub_eq[k] <- p_eq
        Psub_in[k] <- p_in
      } else {
        P[k, ] <- pg[seq_len(R)]
        Psub_eq[k] <- pg[R + 1]
        Psub_in[k] <- pg[R + 2]
        pmethod[k, ] <- "gamma"
        if (pvalue == "hybrid" && run_perm) {
          for (h in seq_len(R)) if (need[h]) { P[k, h] <- p_comp[h]; pmethod[k, h] <- "perm" }
          if (need_eq_in[1]) Psub_eq[k] <- p_eq
          if (need_eq_in[2]) Psub_in[k] <- p_in
          if (need[R + 1]) pmethod[k, R + 1] <- "perm"
        }
      }
    }
    Psub <- if (subspace == "inertia") Psub_in else Psub_eq
    colnames(P) <- paste0("Comp", seq_len(R))
  } else {
  obs <- M(work)
  for (k in seq_len(K)) {
    if (match_on) {
      Uref <- M(work[-k])$U
      pi_obs <- .match_dims(obs$U, Uref)
      s_obs <- obs$d[k, pi_obs]
    } else {
      s_obs <- obs$d[k, ]
    }
    sub_obs <- substat(obs$d, k)
    cnt <- numeric(R)
    cnt_sub <- 0
    if (match_on) refs[[k]] <- Uref
    if (scheme %in% c("rows", "groups", "within")) {
      Pr <- row_perms(B); Xs <- NULL
    } else {
      Pr <- NULL; Xs <- other_perms(k, B)
    }
    permd <- numeric(B)
    for (b in seq_len(B)) {
      wp <- work
      wp[[k]] <- if (!is.null(Pr)) work[[k]][Pr[b, ], , drop = FALSE] else Xs[[b]]
      sol <- M(wp)
      sv <- if (match_on) sol$d[k, .match_dims(sol$U, Uref)] else sol$d[k, ]
      cnt <- cnt + (sv >= s_obs - 1e-12)
      permd[b] <- substat(sol$d, k)
      cnt_sub <- cnt_sub + (permd[b] >= sub_obs - 1e-12)
    }
    p_ref <- (1 + cnt) / (B + 1)
    if (match_on) {
      P[k, pi_obs] <- p_ref
    } else {
      P[k, ] <- p_ref
    }
    Psub[k] <- (1 + cnt_sub) / (B + 1)
    null[[names(blocks)[k]]] <- list(observed = sub_obs, permuted = permd)
  }
  }
  fit <- canofac(blocks, R, q, axes, extraction, tol, maxit, nstart,
                 weighting = weighting)
  # blocks linked to each dimension: Holm adjustment over the K blocks
  P_holm <- apply(P, 2, stats::p.adjust, method = "holm")
  if (K == 1) P_holm <- matrix(P_holm, 1, dimnames = dimnames(P))
  dimnames(P_holm) <- dimnames(P)
  linked <- P_holm <= alpha
  all_linked <- apply(P, 2, max) <= alpha
  link <- .link_code(all_linked, colSums(linked))
  names(link) <- colnames(P)
  # subspace: same rules on the subspace p-values
  linked_sub <- stats::p.adjust(Psub, method = "holm") <= alpha
  names(linked_sub) <- names(Psub)
  link_sub <- .link_code(max(Psub) <= alpha, sum(linked_sub))
  # retained: for "lobo", the component of the other blocks is linked
  # (min_k p <= alpha/K, i.e. at least one block linked after Holm); for
  # "refit", at least one linked block
  retained <- if (statistic == "lobo") apply(P, 2, min) <= alpha / K else
    link != "none"
  out <- list(p = P, p_holm = P_holm, p_subspace = Psub,
              p_subspace_equal = Psub_eq, p_subspace_inertia = Psub_in,
              linked = linked, link = link,
              linked_subspace = linked_sub, link_subspace = link_sub,
              all_linked = all_linked, retained = retained, refs = refs,
              fit = fit, p_compare = comp_p, pmethod = pmethod, null = null,
              settings = list(R = R, q = q, axes = axes, extraction = extraction,
                              scheme = scheme, B = B, alpha = alpha,
                              statistic = statistic, matching = matching,
                              nstart = nstart, subspace = subspace,
                              weighting = weighting, pvalue = pvalue,
                              hybrid_c = hybrid_c, seed = seed,
                              n_groups = if (is.null(gidx)) 0L else length(gidx)))
  class(out) <- "cf_test"
  out
}

# Outcome of the test for a dimension or the subspace: all the blocks linked
# (max_k p <= alpha), otherwise several, one or no block linked after Holm.
.link_code <- function(all_linked, n_linked) {
  ifelse(all_linked, "all",
         ifelse(n_linked >= 2, "several", ifelse(n_linked == 1, "one", "none")))
}

# First r principal components of a block, scaled to unit Frobenius norm (PROTEST).
.pcs <- function(X, r) {
  s <- svd(X, nu = r, nv = 0)
  Z <- sweep(s$u[, seq_len(r), drop = FALSE], 2, s$d[seq_len(r)], "*")
  Z / sqrt(sum(Z^2))
}

.is_centred <- function(X, tol = 1e-8) {
  max(abs(colMeans(X))) <= tol * max(1, max(abs(X)))
}
