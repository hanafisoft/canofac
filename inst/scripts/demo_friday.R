# ==============================================================================
# canofac : exemple complet sur les données Friday (1987)
#   16 étangs, 91 espèces de macro-invertébrés réparties en 10 blocs (groupes taxonomiques)
#   Prétraitement : centrage par colonne, inertie totale unité par bloc
# Pré-requis : install.packages("ade4") ; package canofac installé
# ==============================================================================

library(canofac)

## Option : enregistrer toutes les figures dans un PDF (mettre TRUE)
en_pdf <- FALSE
if (en_pdf) pdf("friday_canofac.pdf", width = 9, height = 6)

## ---- 1. Données -------------------------------------------------------------
blocks <- cf_friday()                         # 10 blocs prétraités
sapply(blocks, ncol)                          # nombre d'espèces par bloc
sapply(blocks, function(X) sum(X^2))          # inertie totale = 1

## ---- 2. Ajustements de la famille en q ---------------------------------------
R <- 2
fit_cpca   <- canofac(blocks, R = R, q = 2)                             # CPCA
fit_indort <- canofac(blocks, R = R, q = 4)                             # INDORT
fit_accps  <- canofac(blocks, R = R, q = 4, extraction = "sequential")  # ACCPS
fit_mcoa   <- canofac(blocks, R = R, q = 2, axes = "orthonormal",
                      extraction = "sequential")                        # MCOA
fit_maxbet <- canofac(blocks, R = R, q = 1)                             # q = 1

summary(fit_indort)          # d, saliences (parts d'inertie), delta

## Comparaison des critères séquentiel / simultané (cf. tableau de l'article)
c(ACCPS = fit_accps$criterion, INDORT = fit_indort$criterion)

## ---- 3. Carte des étangs sur les composantes de référence ----------------------
plot(fit_indort, labels = TRUE, main = "INDORT (q = 4) : étangs")
plot(fit_cpca,   labels = TRUE, main = "CPCA (q = 2) : étangs")

## ---- 4. Pseudo-valeurs singulières, saliences, delta ------------------------
op <- par(mfrow = c(1, 3), mar = c(8, 4, 3, 1))
plot(fit_indort, type = "d",         main = "d_kh")
plot(fit_indort, type = "saliences", main = "Saliences d_kh^2")
plot(fit_indort, type = "delta",     main = "delta_kh = d_kh / sigma_1")
par(op)

## ---- 5. Axes (loadings) de quelques groupes ---------------------------------
plot(fit_indort, type = "loadings",
     block = c("Hemiptera", "Malacostraca", "Mollusca", "Oligochaeta"))

## ---- 6. Roses des pseudo-valeurs singulières --------------------------------
cf_rose(fit_indort, main = "INDORT (q = 4)")          # un panneau par dimension

fits <- list("MAXBET (q = 1)" = fit_maxbet,
             "CPCA (q = 2)"   = fit_cpca,
             "INDORT (q = 4)" = fit_indort)
cf_rose(fits, dim = 1, main = "Dimension 1 selon la méthode")
cf_rose(fits, dim = 2, main = "Dimension 2 selon la méthode")

## ---- 7. Trajectoires des groupes quand q varie ------------------------------
tr <- cf_trajectories(blocks, R = R,
                      qgrid = c(1, 1.5, 1.9, 2, 2.1, 3, 4, 6, 8, 10))
plot(tr, main = "Friday : trajectoires selon q")
round(tr$delta[, , "q=4"], 2)

## ---- 8. Test : blocs liés aux composantes des autres blocs ---------------
## B = 999 : quelques minutes ; réduire à 199 pour un essai rapide
test2 <- cf_test(blocks, R = R, q = 2, B = 999, seed = 1)
test4 <- cf_test(blocks, R = R, q = 4, B = 999, seed = 1)
test2
test4
plot(test2, main = "Test, CPCA (q = 2)")
plot(test4, main = "Test, INDORT (q = 4)")

## Permutation des variables (bloc = bruit ?) au lieu des lignes
test2v <- cf_test(blocks, R = R, q = 2, scheme = "variables", B = 999, seed = 1)
test2v

## ---- 9. Choix du nombre de dimensions ---------------------------------------
## Modèles r = 1, ..., Rmax, Rmax = min(K + 2, n - 1) = 12 pour Friday.
## Long : chaque modèle est testé bloc par bloc ; B = 199 pour un essai.
pth <- cf_test_path(blocks, q = 2, B = 199, seed = 1)
pth                      # profil : dimensions retenues, liées à tous les blocs, à plusieurs, à un seul
plot(pth)                # carte des min_k p (%) par modèle et par dimension
pth$R                    # nombre de dimensions choisi

if (en_pdf) dev.off()
