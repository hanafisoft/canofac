# p-values sans permutation : approximation gamma et mode hybride

Argument `pvalue` de `cf_test()` (statistique LOBO, permutation des lignes).

Les statistiques s'écrivent `T = tr(W_k Pi M Pi')` avec M fixe (`u_h u_h'` par composante, `P` pour T_k^2,
`P W_(-k) P` pour le RV restreint). Leur moyenne et leur variance sous permutation sont exactes
(Kazi-Aoual et al., 1995 ; Josse et al., 2008) :

- E(T) = tr(W_k) tr(M) / (n − 1) ; en particulier E(S_kh^2) = 1/(n − 1) pour un bloc d'inertie 1, la
  « salience de hasard » ;
- Var(T) : fonction explicite de n et, pour W_k et M, de la trace, de la trace du carré et de la somme des
  carrés des termes diagonaux.

## Les trois modes

- `"permutation"` (défaut) : test exact, référence.
- `"gamma"` : loi gamma ajustée sur les deux moments exacts, **sans permutation**. Sur les simulations de
  l'article BB : même décision que la permutation dans 99 % des cas, mais **légèrement libérale** (5,2 % de
  rejets sous H0 au lieu de 4,8 %) ; sur potato, la composante 2 basculerait (0,014 au lieu de 0,018 pour un
  seuil de 0,0167).
- `"hybrid"` : p-value gamma, remplacée par la p-value de permutation quand elle tombe entre
  alpha / (hybrid_c K) et hybrid_c × alpha (hybrid_c = 3 par défaut). Sur les simulations de BB : mêmes
  décisions que la permutation (désaccords 0 % à 5 %, 0,04 % au seuil alpha/K) en ne permutant qu'environ
  un test sur dix ; très petites p-values possibles, au-delà de 1/(B + 1).

L'approximation normale est écartée (trop libérale). `pmethod` indique, pour chaque p-value, la méthode
utilisée. Les fonctions `perm_moments()` et `gamma_pvalue()` donnent les moments et la p-value.
