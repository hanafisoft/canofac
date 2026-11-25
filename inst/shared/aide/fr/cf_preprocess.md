# cf_preprocess() — prétraitement des blocs

```
R       cf_preprocess(blocks, center = TRUE, scale = FALSE, snv = FALSE, unit_inertia = TRUE)
Python  cf_preprocess(blocks, center=True, scale=False, snv=False, unit_inertia=True)
```

Dans l'ordre : SNV des lignes (spectres : chaque spectre centré et réduit par sa moyenne et son écart-type),
centrage des colonnes, réduction des colonnes (variables d'unités différentes), puis inertie totale unité
`||X_k||^2 = 1`. `scale` et `snv` acceptent une valeur par bloc. Cette normalisation est celle des deux
articles ; le test LOBO avec l'approximation gamma suppose des colonnes centrées.

`snv(X)` applique la SNV seule.
