# cf_compare() — pseudo-valeurs de plusieurs méthodes

```
R       cf_compare(fits, names = NULL, digits = 2, match = TRUE)
Python  cf_compare(fits, names=None, digits=2, match=True)
```

Tableau des pseudo-valeurs singulières de plusieurs ajustements (article A, tableaux des pseudo-valeurs par
méthode). Avec `match = TRUE`, les dimensions sont appariées à celles du premier ajustement ; avec
`match = FALSE`, chaque ajustement garde son propre ordre, comme dans les tableaux de l'article. Les noms des méthodes classiques (CPCA, INDORT, MAXBET,
ACOM…) sont donnés automatiquement. Graphique : `cf_plot_compare()`.
