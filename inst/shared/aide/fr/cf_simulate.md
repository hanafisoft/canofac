# cf_simulate() — données simulées de l'article sur le test LOBO

```
R       cf_simulate(scenario = "SA", n = 20, p = 10, K = 4, theta = 1.6, theta2 = 1.5, seed = 1, preprocess = TRUE)
Python  cf_simulate(scenario="SA", n=20, p=10, K=4, theta=1.6, theta2=1.5, seed=1, preprocess=True)
```

Deux directions orthonormées centrées u0, u1 des individus ; dans chaque bloc, un bruit de valeurs propres
proportionnelles à 1/j dans une base orthonormée aléatoire, plus un signal :

| Scénario | Signal |
|---|---|
| S0 | aucun (H0) |
| SA | u0 dans tous les blocs (intensité `theta`) |
| SB | u0 dans tous les blocs, u1 dans le bloc 1 seulement |
| SC | u0 dans les blocs 1 et 2 seulement |
| SD | u0 dans tous les blocs, u1 (`theta2`) dans les blocs 1 et 2 |
| SE | u0 dans les blocs 1 et 2, u1 (`theta2`) dans tous les blocs |

Les intensités sont multipliées par √n. Tous les tirages viennent du générateur commun et les bases
orthonormées sont des facteurs polaires : la même graine donne les mêmes blocs en R et en Python.
Les jeux de l'étude publiée ont été tirés avec un autre générateur ; les taux de rejet obtenus avec
`cf_simulate()` sont équivalents aux fluctuations d'échantillonnage près.
