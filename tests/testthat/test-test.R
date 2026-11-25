test_that("the permutation test returns valid p-values", {
  bl <- .test_blocks()
  t1 <- cf_test(bl, R = 2, q = 2, B = 19, seed = 1)
  expect_s3_class(t1, "cf_test")
  expect_true(all(t1$p >= 1 / 20 & t1$p <= 1))
  expect_length(t1$all_linked, 2)
  t2 <- cf_test(bl, R = 2, q = 4, scheme = "variables", B = 19, seed = 1)
  expect_true(all(t2$p_subspace >= 1 / 20 & t2$p_subspace <= 1))
})

test_that("the test is close to its level under H0", {
  skip_on_cran()
  set.seed(10)
  rej <- replicate(40, {
    bl <- cf_preprocess(lapply(c(4, 6, 5), function(p) matrix(rnorm(20 * p), 20)))
    cf_test(bl, R = 1, q = 2, B = 39)$p[1, 1] <= 0.05
  })
  expect_lt(mean(rej), 0.25)
})

test_that("trajectories have the expected shape", {
  bl <- .test_blocks()
  tr <- cf_trajectories(bl, R = 2, qgrid = c(1, 2, 4))
  expect_equal(dim(tr$delta), c(3, 2, 3))
  expect_true(all(tr$delta <= 1 + 1e-10))
})

test_that("the test path returns a profile and a chosen R", {
  bl <- .test_blocks()
  expect_warning(pth <- cf_test_path(bl, Rmax = 3, q = 2, B = 9, seed = 1), "experimental")
  expect_s3_class(pth, "cf_test_path")
  expect_true(pth$experimental)
  expect_output(print(pth), "EXPERIMENTAL")
  expect_equal(nrow(pth$profile), 3)
  expect_true(pth$R >= 0 && pth$R <= 3)
})

test_that("the linked blocks are given for each dimension and the subspace", {
  bl <- .test_blocks()
  outcomes <- c("all", "several", "one", "none")
  t1 <- cf_test(bl, R = 2, q = 2, B = 19, seed = 2, statistic = "refit")
  expect_true(all(t1$link %in% outcomes))
  expect_equal(dim(t1$linked), c(3, 2))
  # re-estimation statistic: retained iff at least one block is linked
  expect_identical(unname(t1$retained), unname(t1$link != "none"))
  expect_true(all(t1$retained[t1$all_linked]))
  expect_identical(unname(t1$link[t1$all_linked]), rep("all", sum(t1$all_linked)))
  t2 <- cf_test(bl, R = 2, q = 2, B = 19, seed = 2)
  expect_identical(t2$settings$statistic, "lobo")
  expect_true(t2$link_subspace %in% outcomes)
  expect_length(t2$linked_subspace, 3)
  # LOBO statistic: a component is linked iff min_k p <= alpha / K
  expect_identical(unname(t2$retained), unname(apply(t2$p, 2, min) <= 0.05 / 3))
  expect_length(t2$refs, 3)
})

test_that("the LOBO statistic matches its definition", {
  bl <- .test_blocks()
  t2 <- cf_test(bl, R = 2, q = 2, B = 9, seed = 3)
  Uref <- canofac(bl[-1], R = 2, q = 2)$U
  expect_equal(abs(crossprod(t2$refs[[1]], Uref)), diag(2), tolerance = 1e-6)
})

test_that("the LOBO test detects a second link", {
  skip_on_cran()
  set.seed(11)
  n <- 50
  U <- qr.Q(qr(scale(matrix(rnorm(2 * n), n), scale = FALSE)))
  bl <- cf_preprocess(lapply(1:3, function(k) {
    X <- 25 * outer(U[, 1], rnorm(20) / sqrt(20)) + matrix(rnorm(n * 20), n)
    if (k < 3) X <- X + 20 * outer(U[, 2], rnorm(20) / sqrt(20))
    X
  }))
  tf <- cf_test(bl, R = 2, q = 2, B = 99, seed = 1)
  # blocks 1 and 2 carry the second structure
  expect_true(all(tf$p[1:2, 2] <= 0.05))
})
