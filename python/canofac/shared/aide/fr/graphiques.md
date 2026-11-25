# Graphiques

| Fonction | Contenu |
|---|---|
| `plot(fit)` / `cf_plot_scores(fit, group = ...)` | individus sur les composantes de référence (× √n), ellipses à 95 % par groupe |
| `plot(fit, type = "d" / "saliences" / "delta")` | barres par bloc et dimension |
| `plot(fit, type = "loadings")` | axes de chaque bloc ; en spectres si les noms de colonnes sont des longueurs d'onde |
| `cf_plot_saliences(fit, stat = ...)` | carte K × R des saliences, pseudo-valeurs ou delta |
| `plot(traj)` / `cf_plot_trajectories()` | trajectoire de chaque bloc dans le plan (delta_k1, delta_k2) selon q |
| `cf_plot_compare(fits, dim)` | pseudo-valeurs de plusieurs méthodes sur une dimension |
| `cf_rose(fit)` | rose des pseudo-valeurs : écarts à la moyenne, valeurs, part moyenne / variance |
| `plot(test)` / `cf_plot_test()` | carte des p-values (%) ; blocs liés en gras ; issue en haut |
| `plot(test, type = "null")` | distributions de permutation du sous-espace ; trait : valeur observée |
| `plot(path)` | choix de R (**expérimental**) : min_k p_kh par modèle et dimension |
| `plot(ms)` | multimodalité : critère atteint par chaque départ |

Les titres, axes et légendes suivent `cf_lang()`. En R, graphiques de base (`graphics`) ; en Python,
matplotlib (chaque fonction renvoie les axes).
