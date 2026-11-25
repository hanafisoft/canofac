# cf_multistart() — multimodalité du critère

```
R       cf_multistart(blocks, R = 2, q = 4, nstart = 100, extraction = "simultaneous", seed = 2026, rtol = 1e-5)
Python  cf_multistart(blocks, R=2, q=4, nstart=100, extraction="simultaneous", seed=2026, rtol=1e-5)
```

Lance l'algorithme (axes sans contrainte) depuis `nstart` départs aléatoires orthonormés centrés, tirés par
le générateur commun avec la graine `seed`, et depuis le départ déterministe CPCA. Compte la part des départs
qui atteignent le meilleur critère à `rtol` près, et le nombre d'optima distincts (article A, étude de
multimodalité). En extraction séquentielle, chaque étape part du même départ.

Valeur : `criteria` (départs aléatoires), `deterministic`, `best`, `share_best` (%), `n_optima`, `optima`,
`deterministic_best`. Graphique : `plot()`.

Les pourcentages de l'article A ont été obtenus avec un autre générateur ; avec le générateur commun, ils
peuvent différer légèrement, sans changer les conclusions (multimodalité dès q = 4 sur breast.TCGA).
