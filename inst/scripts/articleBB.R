# ==============================================================================
# canofac : reproduction de l'application et d'une étude de simulation réduite de
#   l'article BB « Which blocks are linked to the reference components? A
#   leave-one-block-out permutation test for multiblock dimension reduction »
# Données : potato (package multiblock, ou fichier : cf_potato("potato.csv")).
# Les données simulées viennent de cf_simulate() (générateur commun) ; l'étude
#   publiée a utilisé un autre générateur : taux de rejet équivalents aux
#   fluctuations d'échantillonnage près. Même script en Python :
#   python/scripts/articleBB.py.
# ==============================================================================
library(canofac)
out <- "articleBB_results"
dir.create(out, showWarnings = FALSE)
B <- 999
nrep <- 100                  # 0 : pas de simulation
cf_lang("en")                # "fr" pour des sorties en français

methods <- list(MAXBET = list(q = 1, extraction = "sequential"), CPCA = list(q = 2),
                INDORT = list(q = 4, nstart = 10), MFA = list(q = 2, weighting = "mfa"),
                STATIS = list(q = 2, weighting = "statis"))

## ---- Application : potato -----------------------------------------------------
bl <- cf_potato()            # ou cf_potato("potato.csv")
rows <- list()
for (m in names(methods)) {
  t <- do.call(cf_test, c(list(blocks = bl, R = 2, B = B, seed = 1, compare = TRUE), methods[[m]]))
  rows[[m]] <- data.frame(method = m, t(t$p_subspace), linked = t$link_subspace,
                          min_p_k1 = min(t$p[, 1]), min_p_k2 = min(t$p[, 2]),
                          comp1 = t$retained[1], comp2 = t$retained[2],
                          RV = t(t$p_compare[, "RV"]), PROTEST = t(t$p_compare[, "PROTEST"]),
                          check.names = FALSE)
  if (m == "CPCA") {
    print(t)
    pdf(file.path(out, "figBB_potato_null.pdf"), 10, 3.5); plot(t, type = "null"); dev.off()
    pdf(file.path(out, "figBB_potato_pvalues.pdf"), 6, 4); plot(t); dev.off()
  }
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE, digits = 3)
write.csv(tab, file.path(out, "tableBB_potato.csv"), row.names = FALSE)

## Statistique de réestimation (comparaison)
tr <- cf_test(bl, R = 2, q = 2, B = min(B, 199), seed = 1, statistic = "refit")
cat("\nRe-estimation statistic (CPCA), subspace p-values:", round(tr$p_subspace, 3), "->", tr$link_subspace, "\n")

## Robustesse au choix de R
prof <- do.call(rbind, lapply(c("MAXBET", "CPCA", "INDORT"), function(m) do.call(rbind, lapply(1:5, function(R) {
  t <- do.call(cf_test, c(list(blocks = bl, R = R, B = B, seed = 1), methods[[m]]))
  data.frame(method = m, R = R, t(t$p_subspace),
             linked_components = paste(which(t$retained), collapse = ""), check.names = FALSE)
}))))
print(prof, row.names = FALSE)
write.csv(prof, file.path(out, "tableBB_potato_R.csv"), row.names = FALSE)

## p-values sans permutation (CPCA, composante 2)
g <- cf_test(bl, R = 2, q = 2, seed = 1, pvalue = "gamma")
h <- cf_test(bl, R = 2, q = 2, B = 19999, seed = 1, pvalue = "hybrid")
cat("\nPotato, CPCA: min_k p_k2 by permutation / gamma / hybrid:",
    round(tab["CPCA", "min_p_k2"], 4), round(min(g$p[, 2]), 4), round(min(h$p[, 2]), 4), "\n")

## ---- Simulations (réduites) ---------------------------------------------------
if (nrep > 0) {
  configs <- list(list("S0", 0), list("SA", 1.1), list("SA", 1.6), list("SB", 1.8), list("SC", 1.8),
                  list("SD", 1.8), list("SE", 1.8))
  carriers <- list(S0 = integer(0), SA = 1:4, SB = 1:4, SC = 1:2, SD = 1:4, SE = 1:4)
  sim <- list()
  for (i in seq_along(configs)) for (n in c(20, 100)) {
    scen <- configs[[i]][[1]]; th <- configs[[i]][[2]]
    res <- lapply(seq_len(nrep), function(rep) {
      b <- cf_simulate(scen, n = n, p = 10, K = 4, theta = th,
                       seed = 100000 * (n == 100) + 1000 * (i - 1) + rep)
      t <- cf_test(b, R = 2, q = 2, B = 99, seed = rep, compare = TRUE)
      tgt <- seq_len(4) %in% carriers[[scen]]
      rule <- function(p) if (max(p) <= 0.05) rep(TRUE, 4) else p.adjust(p, "holm") <= 0.05
      c(rv = mean(t$p_subspace_inertia <= 0.05), eq = mean(t$p_subspace_equal <= 0.05),
        RV = mean(t$p_compare[, "RV"] <= 0.05), PROTEST = mean(t$p_compare[, "PROTEST"] <= 0.05),
        exact_rv = all(rule(t$p_subspace_inertia) == tgt), exact_eq = all(rule(t$p_subspace_equal) == tgt))
    })
    r <- round(100 * colMeans(do.call(rbind, res)), 1)
    sim[[length(sim) + 1]] <- data.frame(scenario = scen, theta = th, n = n, t(r))
    print(tail(sim, 1)[[1]], row.names = FALSE)
  }
  sim <- do.call(rbind, sim)
  write.csv(sim, file.path(out, "tableBB_simulations.csv"), row.names = FALSE)
}
