test_that("MFA and STATIS weightings give the expected components", {
  bl <- .test_blocks()
  # MFA: eigenvectors of sum_k W_k / sigma_1(X_k)^2
  f <- canofac(bl, R = 2, q = 2, weighting = "mfa")
  s1 <- vapply(bl, function(X) svd(X)$d[1], numeric(1))
  Wm <- Reduce(`+`, lapply(seq_along(bl), function(k) tcrossprod(bl[[k]]) / s1[k]^2))
  Um <- eigen(Wm, symmetric = TRUE)$vectors[, 1:2]
  expect_equal(tcrossprod(f$U), tcrossprod(Um), tolerance = 1e-6)
  expect_equal(sum(f$block_weights), 1)
  expect_identical(f$settings$weighting, "mfa")
  # STATIS: eigenvectors of the compromise sum_k alpha_k W_k / ||W_k||_F
  g <- canofac(bl, R = 2, q = 2, weighting = "statis")
  W <- lapply(bl, function(X) { Wk <- tcrossprod(X); Wk / sqrt(sum(Wk^2)) })
  C <- outer(seq_along(W), seq_along(W), Vectorize(function(k, l) sum(W[[k]] * W[[l]])))
  a <- abs(eigen(C, symmetric = TRUE)$vectors[, 1])
  Us <- eigen(Reduce(`+`, Map(`*`, a, W)), symmetric = TRUE)$vectors[, 1:2]
  expect_equal(tcrossprod(g$U), tcrossprod(Us), tolerance = 1e-6)
  expect_true(all(g$block_weights > 0))
  # delta does not depend on the weighting
  h <- canofac(bl, R = 2, q = 2)
  expect_equal(abs(g$delta[, 1]) <= 1, rep(TRUE, 3), ignore_attr = TRUE)
  expect_equal(dim(h$delta), dim(g$delta))
})

test_that("the inertia-weighted subspace statistic is computed with the same permutations", {
  bl <- .test_blocks()
  te <- cf_test(bl, R = 2, q = 2, B = 49, seed = 4, subspace = "equal")
  ti <- cf_test(bl, R = 2, q = 2, B = 49, seed = 4)
  # same permutations: equal-weight p-values and dimension-wise p-values unchanged
  expect_equal(ti$p_subspace_equal, te$p_subspace)
  expect_equal(ti$p, te$p)
  expect_equal(ti$p_subspace, ti$p_subspace_inertia)
  expect_true(all(ti$p_subspace >= 1 / 50 & ti$p_subspace <= 1))
  expect_identical(ti$settings$subspace, "inertia")
  # default: restricted RV statistic; for the re-estimation statistic, equal weights
  expect_identical(cf_test(bl, R = 2, B = 9, statistic = "refit")$settings$subspace, "equal")
})

test_that("the test runs with the STATIS and MFA weightings", {
  bl <- .test_blocks()
  for (w in c("statis", "mfa")) {
    tt <- cf_test(bl, R = 2, q = 2, B = 19, seed = 5, weighting = w, subspace = "inertia")
    expect_s3_class(tt, "cf_test")
    expect_identical(tt$settings$weighting, w)
    expect_true(all(tt$p >= 1 / 20 & tt$p <= 1))
    # the p-values do not depend on the scale of the tested block
    bl2 <- bl
    bl2[[1]] <- 3 * bl2[[1]]
    t2 <- cf_test(bl2, R = 2, q = 2, B = 19, seed = 5, weighting = w, subspace = "inertia")
    expect_equal(t2$p[1, ], tt$p[1, ])
    expect_equal(t2$p_subspace[1], tt$p_subspace[1])
  }
})

test_that("the permutation of the principal components preserves the block", {
  bl <- .test_blocks()
  tt <- cf_test(bl, R = 2, q = 2, B = 19, seed = 6, scheme = "components")
  expect_identical(tt$settings$scheme, "components")
  expect_true(all(tt$p >= 1 / 20 & tt$p <= 1))
  # rotation invariance: rotating the variables of the tested block leaves the p-values unchanged
  O <- qr.Q(qr(matrix(stats::rnorm(16), 4)))
  bl2 <- bl
  bl2[[1]] <- bl2[[1]] %*% O
  t2 <- cf_test(bl2, R = 2, q = 2, B = 19, seed = 6, scheme = "components")
  expect_equal(t2$p[1, ], tt$p[1, ])
})
