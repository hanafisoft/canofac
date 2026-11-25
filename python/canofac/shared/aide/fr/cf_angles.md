# cf_angles() — angles entre deux ajustements

```
R       cf_angles(fit1, fit2)
Python  cf_angles(fit1, fit2)
```

Apparie les composantes de `fit2` à celles de `fit1` (affectation optimale), puis donne l'angle (degrés)
entre chaque paire de composantes et le plus grand angle principal entre les deux sous-espaces. Sert à
comparer les extractions séquentielle et simultanée (article A, tableau séquentiel / simultané) ou deux
méthodes.

Valeur : `components` (R angles), `subspace`.
