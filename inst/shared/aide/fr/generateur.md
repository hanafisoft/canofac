# Générateur aléatoire commun (argument seed)

Toutes les tirages aléatoires de canofac (permutations du test, départs aléatoires de l'estimation) viennent
d'un générateur propre au package, **identique en R et en Python** : le générateur « minimal standard » de
Park et Miller (multiplicateur 48271, module 2^31 − 1). Tous les calculs restent exacts en double précision
(R) et en entiers 64 bits (Python). Les permutations sont tirées par l'algorithme de Fisher–Yates ; la
permutation b utilise les tirages b(n − 1) + 1 à (b + 1)(n − 1) du flux, si bien que le résultat ne dépend pas
de la vectorisation.

Conséquences :

- avec la même graine `seed`, `cf_test()` donne **les mêmes p-values en R et en Python** ;
- le générateur global du langage (`set.seed()` en R, `numpy.random` en Python) n'est ni utilisé ni modifié ;
- sans `seed`, une graine est tirée par le générateur du langage (le résultat change d'un appel à l'autre)
  et elle est gardée dans `settings$seed` pour reproduire le calcul ;
- les départs aléatoires de `canofac()` (q > 2) utilisent la graine 12345 : mêmes ajustements en R et en Python.

Ordre des tirages dans `cf_test()` : bloc 1 puis bloc 2…, B permutations par bloc ; pour `"variables"`,
B × p permutations (variable par variable) ; pour `"within"`, groupe par groupe.

`cf_generator(seed)` (Python) et `cf_permutations(B, n, seed)` (R et Python) donnent accès aux tirages.
