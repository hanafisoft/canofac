# ==============================================================================
# canofac, pas à pas / step by step
#
# Ce script parcourt toute la démarche de l'article sur un exemple :
#   1. données et prétraitement        6. choix du nombre de dimensions
#   2. un premier ajustement (CPCA)    7. test bloc par bloc
#   3. la famille en q                 8. croisement des deux schémas
#   4. trajectoires selon q            9. carte des individus (breast.TCGA)
#   5. départs multiples (q > 2)      10. export des figures
#
# Les résultats (impressions et graphiques) sont en français ou en anglais
# selon cf_lang("fr") ou cf_lang("en").
#
# Pré-requis : install.packages("ade4") pour Friday ; pour breast.TCGA (étape 9),
# BiocManager::install("mixOmics").
# Lancement : source(system.file("scripts", "canofac_step_by_step.R",
#                                 package = "canofac"), echo = TRUE)
# ==============================================================================

library(canofac)

## ---- 0. Langue des résultats / language of the results ---------------------
cf_lang("fr")          # "en" pour des libellés en anglais / "en" for English

## Mettre TRUE pour enregistrer toutes les figures dans un PDF (étape 10)
en_pdf <- FALSE
if (en_pdf) pdf("canofac_pas_a_pas.pdf", width = 9, height = 6.5)

## Nombre de permutations : 999 pour des résultats définitifs (quelques minutes),
## 199 pour un essai rapide.
B <- 199


## ---- 1. Données et prétraitement ---------------------------------------------
## Friday (1987) : 16 étangs, 91 espèces réparties en 10 blocs (groupes taxonomiques).
## cf_friday() renvoie les blocs déjà prétraités : centrage des colonnes, puis
## inertie totale unité pour chaque bloc.
blocks <- cf_friday()
sapply(blocks, ncol)                    # nombre de variables p_k de chaque bloc
sapply(blocks, function(X) sum(X^2))    # inertie totale de chaque bloc : 1

## Pour vos propres données : une liste de matrices ayant les mêmes lignes
## (mêmes individus, même ordre), puis cf_preprocess().
##   mes_blocs <- cf_preprocess(list(A = XA, B = XB, C = XC))            # centrage + inertie 1
##   mes_blocs <- cf_preprocess(list(...), scale = TRUE)                  # variables d'unités différentes
##   mes_blocs <- cf_preprocess(list(...), snv = TRUE)                    # spectres (SNV)
## (voir ?cf_preprocess pour les options exactes)


## ---- 2. Un premier ajustement : la CPCA (q = 2) -------------------------------
## Chaque bloc est résumé par les composantes de référence U (communes à tous
## les blocs par construction) ; d_kh mesure ce que la dimension h reproduit du
## bloc k ; la salience d_kh^2 est la part d'inertie du bloc portée par h.
R <- 2
fit_cpca <- canofac(blocks, R = R, q = 2)
fit_cpca                                   # pseudo-valeurs singulières d_kh
summary(fit_cpca)                          # d_kh, saliences, delta_kh

cf_plot_saliences(fit_cpca)                # tableau des saliences
cf_plot_saliences(fit_cpca, stat = "delta")# delta_kh = d_kh / sigma_1(X_k)
cf_plot_scores(fit_cpca, labels = TRUE)    # étangs sur les composantes de référence


## ---- 3. La famille en q ---------------------------------------------------------
## Le critère sum_k sum_h |d_kh|^q retrouve les méthodes classiques :
##   q = 1 : MAXBET ; q = 2 : CPCA (axes libres), MCOA (axes orthonormés) ;
##   q = 4 : ACCPS (séquentiel), INDORT (simultané).
## En extraction simultanée, q < 2 égalise les saliences d'un bloc entre les
## dimensions, q > 2 les concentre sur une dimension.
fit_q1     <- canofac(blocks, R = R, q = 1)
fit_indort <- canofac(blocks, R = R, q = 4)                              # INDORT
fit_accps  <- canofac(blocks, R = R, q = 4, extraction = "sequential")  # ACCPS
fit_mcoa   <- canofac(blocks, R = R, q = 2, axes = "orthonormal",
                      extraction = "sequential")                         # MCOA
fit_smcoa  <- canofac(blocks, R = R, q = 2, axes = "orthonormal")       # MCOA simultanée

## Pseudo-valeurs singulières par méthode, dimension 1 puis 2
fits <- list("q = 1" = fit_q1, "CPCA" = fit_cpca, "INDORT" = fit_indort)
cf_plot_compare(fits, dim = 1)
cf_plot_compare(fits, dim = 2)

## Roses : moyenne et dispersion des d_kh sur les blocs, dimension par dimension
cf_rose(fits, dim = 1)

## Extraction séquentielle et simultanée : critères atteints
c(ACCPS = fit_accps$criterion, INDORT = fit_indort$criterion)
c(MCOA = fit_mcoa$criterion, "MCOA simultanée" = fit_smcoa$criterion)


## ---- 4. Trajectoires selon q -----------------------------------------------------
## Pour chaque bloc, le point (delta_k1, delta_k2) quand q va de 1 à 10.
## Marqueurs : q = 1, q = 2 (CPCA), q = 4 (INDORT) et q = 10.
tr <- cf_trajectories(blocks, R = R, qgrid = c(1, 1.5, 1.9, 2, 2.1, 3, 4, 6, 8, 10))
cf_plot_trajectories(tr)
round(tr$delta[, , "q=4"], 2)


## ---- 5. Départs multiples pour q > 2 --------------------------------------------
## À partir de q = 4, le critère peut avoir plusieurs maxima locaux. Par
## défaut, canofac() utilise 10 départs (CPCA + 9 départs aléatoires fixés) si
## q > 2 ; comparer avec un seul départ :
fit_1 <- canofac(blocks, R = R, q = 4, nstart = 1)
fit_20 <- canofac(blocks, R = R, q = 4, nstart = 20)
c(un_depart = fit_1$criterion, vingt_departs = fit_20$criterion)


## ---- 6. Choix du nombre de dimensions ------------------------------------------
## Modèles à r = 1, ..., Rmax dimensions, chacun réestimé et testé ; on retient
## le plus grand r dont toutes les composantes des autres blocs sont liées
## (min_k p_kh <= alpha / K). Le test à référence fixe est rapide.
pth <- cf_test_path(blocks, Rmax = 4, q = 2, B = B, seed = 1)
pth                     # profil des modèles
plot(pth)               # min_k p_kh (%) par modèle et composante ; liées en gras
pth$R                   # nombre de dimensions choisi


## ---- 7. Test bloc par bloc ---------------------------------------------------------
## Statistique à référence fixe (par défaut) : pour chaque bloc k, les composantes
## U^(-k) sont estimées une fois sur les autres blocs ; on compare le lien entre
## le bloc k et ces composantes à celui obtenu en permutant les lignes du bloc k.
## Blocs liés au sous-espace :
##   tous les blocs   si max_k p_k <= alpha ;
##   sinon, les blocs dont la p-value ajustée par Holm est <= alpha
##   (plusieurs, un seul ou aucun).
## Être lié n'est pas une étiquette : une composante peut mélanger des
## directions que les blocs portent différemment. Pour étiqueter des parties
## de l'espace retenu (commune, partagée, propre, résiduelle), voir CoPRes
## (https://github.com/hanafisoft/CoPRes).
## Les p-values par composante des autres blocs sont descriptives ; une
## composante est liée si min_k p_kh <= alpha / K.
test2 <- cf_test(blocks, R = R, q = 2, scheme = "rows", B = B, seed = 1)
test2
cf_plot_test(test2)     # p-values (%), blocs liés en gras, issue du sous-espace en capitales

test4 <- cf_test(blocks, R = R, q = 4, scheme = "rows", B = B, seed = 1)
cf_plot_test(test4)

## Les valeurs stockées ne dépendent pas de la langue :
test2$link_subspace      # "all", "several", "one" ou "none"
test2$linked_subspace    # blocs liés au sous-espace des autres blocs
test2$retained           # composantes des autres blocs liées


## ---- 8. Croiser les deux schémas ------------------------------------------------
## Une direction portée par un seul bloc n'est pas dans U^(-k) : la
## statistique à référence fixe ne la voit pas. On la reconnaît
## en croisant la statistique de réestimation avec permutation des variables
## (H0(k) = « le bloc k est du bruit ») et le test à référence fixe : le bloc
## n'est pas du bruit dans cette direction, mais elle n'est pas liée aux autres.
test2v <- cf_test(blocks, R = R, q = 2, scheme = "variables", B = B, seed = 1,
                  statistic = "refit")
cf_plot_test(test2v)
data.frame(reestimation_variables = test2v$link)

## Pour comparaison : la statistique de réestimation avec permutation des lignes
## perd le second lien quand une première structure liant les blocs est forte.
test2r <- cf_test(blocks, R = R, q = 2, scheme = "rows", B = B, seed = 1,
                  statistic = "refit")
round(100 * cbind(reference_fixe = test2$p_subspace, reestimation = test2r$p_subspace), 1)


## ---- 9. Carte des individus avec groupes : breast.TCGA --------------------------
## Nécessite mixOmics (Bioconductor). Ellipses à 95 % par sous-type.
if (requireNamespace("mixOmics", quietly = TRUE)) {
  tc <- cf_tcga()
  fit_tc2 <- canofac(tc$blocks, R = 2, q = 2)
  fit_tc4 <- canofac(tc$blocks, R = 2, q = 4)
  op <- par(mfrow = c(1, 2))
  cf_plot_scores(fit_tc2, group = tc$subtype, main = "CPCA (q = 2)")
  cf_plot_scores(fit_tc4, group = tc$subtype, main = "INDORT (q = 4)")
  par(op)
  cf_plot_saliences(fit_tc4)
}


## ---- 10. Export des figures -------------------------------------------------------
if (en_pdf) dev.off()

## Revenir aux libellés anglais / back to English labels
## cf_lang("en")
