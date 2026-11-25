# cf_test_path() — choix de R par tests successifs (EXPÉRIMENTAL)

> **Fonction expérimentale.** Le choix du nombre de dimensions R ne fait pas partie des articles publiés
> (A et BB). Un avertissement est émis à chaque appel, et l'affichage du résultat le rappelle.

```
R       cf_test_path(blocks, Rmax = NULL, c = 2, q = 2, ..., B = 999, alpha = 0.05, seed = NULL)
Python  cf_test_path(blocks, Rmax=None, c=2, q=2, ..., B=999, alpha=0.05, seed=None)
```

Pour r = 1, …, Rmax (par défaut min(K + c, n − 1)), le modèle à r dimensions est réestimé puis testé par
`cf_test()` (graine `seed + r − 1` pour le modèle r). R est le plus grand r dont toutes les composantes des
autres blocs sont liées (min_k p_kh ≤ alpha / K).

## Limites

- Les modèles sont réestimés pour chaque r : ils ne sont pas emboîtés pour les méthodes simultanées
  (avec INDORT, u_1 peut tourner de plus de 30° entre R = 1 et R = 2).
- Le choix utilise les mêmes données que les tests : les p-values du modèle retenu sont optimistes.
- Les tests successifs ne sont pas corrigés pour la multiplicité.

Présenter **tout le profil** (`profile`, `plot()`), pas seulement R. Pour une conclusion qui ne dépende pas de
R, utiliser plutôt le test du sous-espace (RV restreint), peu sensible à un R un peu trop grand.
