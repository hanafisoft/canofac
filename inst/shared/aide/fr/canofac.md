# canofac() — factorisation canonique, famille de critères en puissance q

Ajuste `X_k = U D_k V_k' + E_k` en maximisant `somme_k somme_h |d_kh|^q`, où `d_kh` est la pseudo-valeur
singulière du bloc k sur la dimension h (article A).

## Usage

```
R       canofac(blocks, R = 2, q = 2, axes = "free", extraction = "simultaneous",
                tol = 1e-10, maxit = 5000, nstart = NULL, weighting = "none")
Python  canofac(blocks, R=2, q=2, axes="free", extraction="simultaneous",
                tol=1e-10, maxit=5000, nstart=None, weighting="none")
```

## Arguments

- `blocks` : blocs prétraités (mêmes individus en lignes), voir `cf_preprocess()`.
- `R` : nombre de dimensions.
- `q` : puissance appliquée aux pseudo-valeurs, `q >= 1`.
- `axes` : `"free"` (axes sans contrainte, de norme 1) ou `"orthonormal"` (axes orthonormés dans chaque bloc).
- `extraction` : `"simultaneous"` (les R dimensions ensemble) ou `"sequential"` (une à une, avec déflation).
- `tol`, `maxit` : tolérance relative sur le critère, nombre maximal d'itérations.
- `nstart` : nombre de départs pour les axes sans contrainte : le départ déterministe (CPCA) et
  `nstart - 1` départs aléatoires tirés par le générateur commun (graine 12345, mêmes départs en R et en
  Python) ; le meilleur critère est retenu. Par défaut 10 si `q > 2`, 1 sinon.
- `weighting` : `"none"`, `"mfa"` (AFM : bloc divisé par sa première valeur singulière) ou `"statis"`
  (STATIS : `||W_k|| = 1` puis poids du compromis).

## Méthodes retrouvées

| q | Axes sans contrainte | Axes orthonormés |
|---|---|---|
| 1 | MAXBET (séquentiel) | MAXBET simultané (contrainte d'identification) |
| 2 | CPCA (séquentiel = simultané) ; AFM, STATIS par pondération | ACOM séquentielle ; ACOM simultanée |
| 4 | ACCPS (séquentiel) ; INDORT (simultané) | — |

## Valeur

- `U` : composantes de référence (n × R), orthonormées ;
- `V` : axes de chaque bloc ;
- `d` : pseudo-valeurs singulières (K × R) ; `saliences` = `d^2` (part d'inertie du bloc, inertie totale 1) ;
- `delta` = `d / sigma_1(X_k)` (entre 0 et 1) ; `sigma1` ; `criterion` ; `block_weights` ; `settings`.

`print` (R) ou l'affichage de l'objet (Python) donne les pseudo-valeurs ; `summary` y ajoute les saliences et
les `delta`. Les dimensions sont rangées par contribution décroissante au critère.

## Exemple

```
R       bl <- cf_preprocess(list(A = XA, B = XB, C = XC)); fit <- canofac(bl, R = 2, q = 4); summary(fit)
Python  bl = cf_preprocess({"A": XA, "B": XB, "C": XC}); fit = canofac(bl, R=2, q=4); fit.summary()
```
