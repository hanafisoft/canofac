# ==============================================================================
# canofac : reproduction des tableaux et figures de l'article A
#   « Pseudo-singular values: a family of power criteria for multiblock dimension
#     reduction »
# Données : potato, blocs bruts NIRraw, FIDraw, CPMGraw (multiblock), breast.TCGA (mixOmics),
#   Friday (ade4).
# Sorties : tableaux (CSV) et figures (PDF) dans le dossier `out`.
# Remarque : les départs aléatoires viennent du générateur commun de canofac (graine
#   2026), celui de l'article : R et Python donnent les mêmes tableaux. Même script en
#   Python : python/scripts/articleA.py.
# ==============================================================================
library(canofac)
out <- "articleA_results"
dir.create(out, showWarnings = FALSE)
cf_lang("en")                 # "fr" pour des sorties en français

data <- list()
if (requireNamespace("multiblock", quietly = TRUE)) data[["Potato"]] <- list(blocks = cf_potato(blocks = "A"))
if (requireNamespace("mixOmics", quietly = TRUE)) data[["breast.TCGA"]] <- cf_tcga()
if (requireNamespace("ade4", quietly = TRUE)) data[["Friday"]] <- list(blocks = cf_friday())
cat("datasets:", paste(names(data), collapse = ", "), "\n")

## ---- Tableau : multimodalité (100 départs aléatoires, graine 2026) -------------
mm <- do.call(rbind, lapply(names(data), function(nm) {
  do.call(rbind, lapply(c("simultaneous", "sequential"), function(ex) {
    cells <- vapply(c(1, 2, 4, 6, 10), function(q) {
      ms <- cf_multistart(data[[nm]]$blocks, R = 2, q = q, nstart = 100, extraction = ex, seed = 2026)
      paste0(round(ms$share_best), if (ms$n_optima > 1) paste0(" (", ms$n_optima, ")") else "",
             if (ms$deterministic_best) "" else "*")
    }, character(1))
    data.frame(dataset = nm, extraction = ex, t(setNames(cells, paste0("q=", c(1, 2, 4, 6, 10)))),
               check.names = FALSE)
  }))
}))
print(mm, row.names = FALSE)
write.csv(mm, file.path(out, "tableA_multimodality.csv"), row.names = FALSE)

## ---- Figure : trajectoires selon q (meilleur de 20 départs pour q > 2) ----------
qgrid <- c(1, 1.25, 1.5, 1.75, 1.9, 2, 2.1, 2.5, 3, 4, 6, 8, 10)
for (nm in names(data)) {
  pdf(file.path(out, paste0("figA_traj_", strsplit(nm, " ")[[1]][1], ".pdf")), 6, 6)
  plot(cf_trajectories(data[[nm]]$blocks, R = 2, qgrid = qgrid, nstart = 20), main = nm)
  dev.off()
}

## ---- Tableaux : pseudo-valeurs par méthode (ordre propre à chaque méthode) -----
for (nm in names(data)) {
  bl <- data[[nm]]$blocks
  fits <- list(canofac(bl, 2, 1), canofac(bl, 2, 1, axes = "orthonormal"), canofac(bl, 2, 2),
               canofac(bl, 2, 4, nstart = 10))
  cat("\n", nm, "\n", sep = "")
  tab <- cf_compare(fits, names = c("q=1, unconstrained", "simultaneous MAXBET", "CPCA", "INDORT"),
                    match = FALSE)
  write.csv(tab, file.path(out, paste0("tableA_methods_", strsplit(nm, " ")[[1]][1], ".csv")))
}

## ---- Tableau : extraction séquentielle et simultanée (Friday) ------------------
if (!is.null(data$Friday)) {
  bl <- data$Friday$blocks
  ss <- do.call(rbind, lapply(list(list("q=1, unconstrained", 1, "free"), list("q=2, unconstrained", 2, "free"),
                                   list("q=2, orthonormal", 2, "orthonormal"), list("q=4, unconstrained", 4, "free")),
                              function(s) {
    fs <- canofac(bl, 2, s[[2]], s[[3]], "sequential")
    fm <- canofac(bl, 2, s[[2]], s[[3]], "simultaneous")
    a <- cf_angles(fs, fm)
    data.frame(criterion = s[[1]], angle = round(max(a$components), 1),
               crit_seq = round(fs$criterion, 3), crit_sim = round(fm$criterion, 3))
  }))
  print(ss, row.names = FALSE)
  write.csv(ss, file.path(out, "tableA_seq_sim_Friday.csv"), row.names = FALSE)
}

## ---- Figure : individus de breast.TCGA (CPCA et INDORT) -----------------------
if (!is.null(data$breast.TCGA)) {
  pdf(file.path(out, "figA_indiv_tcga.pdf"), 11, 5.5)
  op <- par(mfrow = c(1, 2))
  cf_plot_scores(canofac(data$breast.TCGA$blocks, 2, 2), group = data$breast.TCGA$subtype, main = "CPCA")
  cf_plot_scores(canofac(data$breast.TCGA$blocks, 2, 4, nstart = 20), group = data$breast.TCGA$subtype,
                 main = "INDORT")
  par(op)
  dev.off()
}
