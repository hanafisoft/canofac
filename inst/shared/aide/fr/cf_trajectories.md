# cf_trajectories() — trajectoires des blocs selon q

```
R       cf_trajectories(blocks, R = 2, qgrid = c(1, 1.5, 1.9, 2, 2.1, 3, 4, 6, 8, 10),
                        axes = "free", extraction = "simultaneous", tol, maxit, nstart = NULL)
Python  cf_trajectories(blocks, R=2, qgrid=(1, 1.5, 1.9, 2, 2.1, 3, 4, 6, 8, 10), ...)
```

Pour chaque valeur de q, ajuste la méthode et range les pseudo-valeurs. Les dimensions sont appariées d'une
valeur de q à la suivante en partant de la plus grande (affectation optimale), pour suivre chaque bloc dans
le plan `(delta_k1, delta_k2)` (article A, figure des trajectoires). Le saut à q = 2 (bascule entre régime
d'égalisation et régime de concentration) apparaît comme une discontinuité.

Valeur : `q`, `d` et `delta` (tableaux K × R × nombre de q), `sigma1`, `criterion`.
Graphique : `plot()` / `cf_plot_trajectories()` ; marqueurs à q = 1, 2 (CPCA), 4 (INDORT) et la plus grande valeur.
