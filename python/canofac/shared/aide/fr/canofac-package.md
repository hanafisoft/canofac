# canofac : factorisation canonique des données multiblocs

**canofac** réunit le code des deux articles :

- **article A** — *Pseudo-singular values: a family of power criteria for multiblock dimension reduction* :
  famille de critères en puissance q sur les pseudo-valeurs singulières (MAXBET, ACP consensus CPCA, ACOM,
  ACCPS/ComDim, INDORT ; STATIS et AFM par pondération des blocs), trajectoires selon q, multimodalité,
  extraction séquentielle et simultanée ;
- **article BB** — *Which blocks are linked to the reference components? A leave-one-block-out permutation
  test for multiblock dimension reduction* : test de permutation LOBO.

Le package existe en **R** et en **Python**, avec les mêmes fonctions, les mêmes arguments, les mêmes
résultats, les mêmes libellés (anglais ou français, `cf_lang()`) et le **même générateur aléatoire** :
avec la même graine (`seed`), les permutations, donc les p-values, sont identiques dans les deux langages.

## Fonctions

| Tâche | Fonction | Aide |
|---|---|---|
| Prétraitement (centrage, réduction, SNV, inertie unité) | `cf_preprocess()` | `cf_help("cf_preprocess", "fr")` |
| Estimation d'une méthode de la famille | `canofac()` | `cf_help("canofac", "fr")` |
| Trajectoires des blocs selon q | `cf_trajectories()` | `cf_help("cf_trajectories", "fr")` |
| Multimodalité (départs multiples) | `cf_multistart()` | `cf_help("cf_multistart", "fr")` |
| Angles entre deux ajustements | `cf_angles()` | `cf_help("cf_angles", "fr")` |
| Pseudo-valeurs de plusieurs méthodes | `cf_compare()` | `cf_help("cf_compare", "fr")` |
| Test LOBO | `cf_test()` | `cf_help("cf_test", "fr")` |
| p-values sans permutation (gamma, hybride) | `cf_test(pvalue = ...)` | `cf_help("approximation", "fr")` |
| Données simulées (scénarios de l'article BB) | `cf_simulate()` | `cf_help("cf_simulate", "fr")` |
| Choix de R (**expérimental**) | `cf_test_path()` | `cf_help("cf_test_path", "fr")` |
| Graphiques | `cf_plot_*()`, `cf_rose()`, `plot()` | `cf_help("graphiques", "fr")` |
| Données | `cf_read_blocks()`, `cf_write_blocks()`, `cf_potato()`, `cf_friday()`, `cf_tcga()` | `cf_help("donnees", "fr")` |
| Générateur aléatoire commun | argument `seed` | `cf_help("generateur", "fr")` |
| Langue des sorties | `cf_lang()` | `cf_help("cf_lang", "fr")` |

## Reproduire les articles

Scripts `articleA` et `articleBB`, en R (`system.file("scripts", package = "canofac")`) et en Python
(`python/scripts`).

## R et Python

En R, les blocs sont une liste nommée de matrices ; en Python, un dictionnaire `nom -> tableau`
(numpy ou pandas). `TRUE`/`FALSE` en R s'écrivent `True`/`False` en Python ; les résultats sont des listes
en R et des objets à attributs en Python (`test$p_subspace` en R, `test.p_subspace` en Python).

## Auteur

Mohamed Hanafi (mohamed.hanafi@oniris-nantes.fr), auteur et mainteneur. Licence GPL-3.
