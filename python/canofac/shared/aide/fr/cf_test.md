# cf_test() — test LOBO : quels blocs sont liés aux composantes de référence ?

```
R       cf_test(blocks, R = 2, q = 2, axes = "free", extraction = "simultaneous",
                scheme = "rows", groups = NULL, B = 999, alpha = 0.05, statistic = "lobo",
                matching = TRUE, seed = NULL, tol, maxit, nstart = NULL, subspace = "inertia",
                weighting = "none", pvalue = "permutation", hybrid_c = 3, compare = FALSE)
Python  cf_test(blocks, R=2, q=2, ..., pvalue="permutation", hybrid_c=3, compare=False)
```

## Principe (article BB)

Pour chaque bloc k :

1. **Retirer le bloc k** et appliquer la même procédure M (normalisation, méthode, départs, ordre des
   dimensions) aux K − 1 autres blocs : on obtient les composantes de référence `U^(-k)`, **fixées**.
2. **Mesurer le lien** entre le bloc k et `U^(-k)` : par composante, `S_kh = ||X_k' u_h^(-k)||` (descriptif) ;
   pour le sous-espace, le **RV restreint** `T_k^w = tr(W_k P W_(-k) P)` (défaut, `subspace = "inertia"`),
   qui pondère chaque composante par l'inertie des autres blocs, ou `T_k^2 = tr(W_k P)` (`"equal"`).
3. **Permuter** les lignes du bloc k, B fois ; `U^(-k)` ne change pas.
4. **Comparer** : p-value = (1 + nombre de valeurs permutées ≥ valeur observée) / (B + 1).

Décisions : **tous les blocs liés** si max_k p_k ≤ alpha (test intersection–union) ; sinon les blocs liés
après ajustement de **Holm**. Une composante des autres blocs est liée si min_k p_kh ≤ alpha / K.

## Arguments propres au test

- `scheme` : `"rows"` (lignes du bloc) ; `"groups"` (groupes entiers de lignes, de même taille : réplicats) ;
  `"within"` (lignes à l'intérieur de chaque groupe) ; `"variables"` (chaque variable indépendamment :
  H0 « le bloc est du bruit ») ; `"components"` (chaque composante principale du bloc indépendamment).
- `groups` : groupe de chaque ligne, pour `"groups"` et `"within"`. Avec des individus structurés
  (réplicats, lots), la permutation libre des lignes est trop libérale : permuter des groupes entiers.
- `B` : nombre de permutations ; au minimum B ≥ K/alpha − 1 ; **10 000 à 20 000 près d'un seuil corrigé**
  (B ≥ 1,96² α′(1 − α′)/ε² pour connaître p à ± ε près autour du seuil α′).
- `statistic` : `"lobo"` (défaut ; `"fixed"` est un alias) ou `"refit"` (réestimation du modèle sur chaque
  jeu permuté, avec appariement `matching` : comparaison seulement, moins puissante pour la 2e direction).
- `seed` : graine du **générateur commun** : même graine, mêmes permutations et mêmes p-values en R et en
  Python (`cf_help("generateur", "fr")`).
- `pvalue` : `"permutation"` (défaut), `"gamma"` ou `"hybrid"` (`cf_help("approximation", "fr")`).
- `compare` : ajoute le **RV complet** (bloc contre la concaténation des autres blocs) et **PROTEST**
  (R premières composantes principales), avec les mêmes permutations. Ce sont des comparaisons, hors du test
  LOBO : le RV complet est plus puissant en moyenne, mais il dit si le bloc est lié aux autres blocs, pas aux
  composantes retenues.

## Valeur

`p` (K × R, par composante), `p_holm`, `p_subspace`, `p_subspace_equal`, `p_subspace_inertia`, `linked`,
`link` (`"all"`, `"several"`, `"one"`, `"none"`), `linked_subspace`, `link_subspace`, `all_linked`,
`retained`, `refs` (les `U^(-k)`), `p_compare`, `pmethod` (`"perm"` ou `"gamma"` pour chaque p-value),
`null` (distributions de permutation du sous-espace), `fit`, `settings`. Les valeurs stockées ne dépendent
pas de la langue.

## Remarques

- Permuter les lignes de X_k ou celles de la référence entière donne le même test.
- Avec U^(-k) orthonormée, T_k^2 est la co-inertie entre X_k et U^(-k), et l'inertie expliquée par la RDA de
  X_k sur U^(-k) (mêmes p-values) ; T_k^w devient le numérateur du RV complet quand R = rang(X_(-k)).
- Le test du sous-espace ne dépend pas de la rotation des composantes ; les tests par composante si.
- Graphiques : `plot(test)` (carte des p-values, blocs liés en gras) et `plot(test, type = "null")`
  (distributions de permutation, trait : valeur observée).
