test_that("pseudo-singular values match the reference implementation", {
  bl <- .test_blocks()
  for (e in .expected) {
    fit <- canofac(bl, R = 2, q = e$q, axes = e$axes, extraction = e$extraction,
                   tol = 1e-14, maxit = 1e5)
    expect_equal(unname(fit$d), e$d, tolerance = 1e-5,
                 label = paste(e$q, e$axes, e$extraction))
    expect_equal(fit$criterion, e$crit, tolerance = 1e-5)
  }
})

test_that("preprocessing gives centred blocks of unit total inertia", {
  bl <- .test_blocks()
  for (X in bl) {
    expect_equal(sum(X^2), 1)
    expect_equal(unname(colMeans(X)), rep(0, ncol(X)), tolerance = 1e-12)
  }
  S <- snv(matrix(1:12 + 0.5 * (1:12)^2, 3))
  expect_equal(rowMeans(S), rep(0, 3), tolerance = 1e-12)
  expect_equal(apply(S, 1, sd), rep(1, 3))
})

test_that("q = 2 with free axes is CPCA, sequential = simultaneous", {
  bl <- .test_blocks()
  s <- canofac(bl, 2, 2, "free", "sequential")
  m <- canofac(bl, 2, 2, "free", "simultaneous")
  expect_equal(s$criterion, m$criterion, tolerance = 1e-8)
  ev <- eigen(Reduce(`+`, lapply(bl, tcrossprod)), symmetric = TRUE)$values
  expect_equal(unname(colSums(m$d^2)), ev[1:2], tolerance = 1e-8)
})

test_that("MCOA saliences are bounded by the saliences on the initial tables", {
  bl <- .test_blocks()
  f <- canofac(bl, 2, 2, "orthonormal", "simultaneous")
  a <- t(vapply(bl, function(X) colSums(crossprod(X, f$U)^2), numeric(2)))
  expect_true(all(f$d^2 <= a + 1e-10))
})

test_that("common components are orthonormal", {
  bl <- .test_blocks()
  for (ax in c("free", "orthonormal")) for (ex in c("simultaneous", "sequential")) {
    f <- canofac(bl, 2, 4, ax, ex)
    expect_equal(unname(crossprod(f$U)), diag(2), tolerance = 1e-8)
  }
})

test_that("multiple starts keep the best solution and leave the RNG untouched", {
  bl <- .test_blocks()
  set.seed(42); a <- runif(1)
  set.seed(42); f10 <- canofac(bl, 2, 4, nstart = 10); b <- runif(1)
  expect_identical(a, b)
  f1 <- canofac(bl, 2, 4, nstart = 1)
  expect_gte(f10$criterion, f1$criterion - 1e-10)
  expect_equal(f10$settings$nstart, 10L)
  expect_equal(canofac(bl, 2, 2)$settings$nstart, 1L)
})
