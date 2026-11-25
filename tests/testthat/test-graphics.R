test_that("the optimal matching is exact and agrees with exhaustive search", {
  perms <- function(R) {
    if (R == 1) return(matrix(1L, 1, 1))
    P <- perms(R - 1)
    do.call(rbind, lapply(seq_len(R), function(i)
      cbind(i, matrix(setdiff(seq_len(R), i)[P], nrow = nrow(P)))))
  }
  set.seed(4)
  for (it in 1:60) {
    R <- sample(1:6, 1)
    C <- abs(matrix(rnorm(R * R), R))
    a <- canofac:::.lsap_max(C)
    expect_setequal(a, seq_len(R))
    P <- perms(R)
    best <- max(apply(P, 1, function(p) sum(C[cbind(seq_len(R), p)])))
    expect_equal(sum(C[cbind(seq_len(R), a)]), best, tolerance = 1e-12)
  }
  # identity when U = Uref, permutation recovered otherwise
  U <- qr.Q(qr(matrix(rnorm(40), 10, 4)))
  expect_equal(canofac:::.match_dims(U, U), 1:4)
  expect_equal(canofac:::.match_dims(U[, c(3, 1, 4, 2)], U), c(2, 4, 1, 3))
})

test_that("the language option switches the labels", {
  old <- cf_lang("fr")
  on.exit(cf_lang(old))
  expect_identical(cf_lang(), "fr")
  expect_identical(canofac:::.link_label(c("all", "none")), c("tous les blocs", "aucun bloc"))
  cf_lang("en")
  expect_identical(canofac:::.link_label("several"), "several blocks")
  expect_error(cf_lang("de"))
})

test_that("the printed results are available in both languages", {
  bl <- .test_blocks()
  fit <- canofac(bl, R = 2, q = 2)
  tt <- cf_test(bl, R = 2, q = 2, B = 9, seed = 1)
  old <- cf_lang("fr")
  on.exit(cf_lang(old))
  expect_output(print(fit), "Pseudo-valeurs")
  expect_output(print(tt), "Test par permutation")
  cf_lang("en")
  expect_output(print(fit), "Pseudo-singular values")
  expect_output(print(tt), "Permutation test")
})

test_that("the graphics functions run", {
  bl <- .test_blocks()
  fit <- canofac(bl, R = 2, q = 2)
  tr <- cf_trajectories(bl, R = 2, qgrid = c(1, 2, 4, 6))
  tt <- cf_test(bl, R = 2, q = 2, B = 9, seed = 1)
  g <- rep(c("a", "b", "c"), 4)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (l in c("en", "fr")) {
    old <- cf_lang(l)
    expect_invisible(cf_plot_scores(fit, group = g))
    expect_invisible(cf_plot_saliences(fit))
    expect_invisible(cf_plot_saliences(fit, stat = "delta"))
    expect_invisible(cf_plot_trajectories(tr))
    expect_invisible(cf_plot_test(tt))
    expect_invisible(cf_plot_compare(list(canofac(bl, 2, 1), fit, canofac(bl, 2, 4))))
    expect_invisible(plot(fit))
    expect_invisible(plot(tr))
    expect_invisible(plot(tt))
    cf_lang(old)
  }
})
