# canofac (Python)

Canonical factorisation for multiblock data — Python version of the R package
[canofac](https://github.com/hanafisoft/canofac), with the same functions, arguments, results, labels
(English or French) and random draws: with the same `seed`, the permutation test gives the same p-values in
R and in Python.

*Factorisation canonique des données multiblocs — version Python du package R canofac : mêmes fonctions,
mêmes arguments, mêmes résultats, mêmes libellés (anglais ou français) et mêmes tirages aléatoires.
Aide en français : `cf_help("canofac-package", "fr")`.*

## Contents

- **Article A** — pseudo-singular values and the power-q family: `canofac()` (MAXBET, CPCA, MCOA, ACCPS,
  INDORT; STATIS and MFA by weighting), `cf_trajectories()`, `cf_multistart()`, `cf_angles()`, `cf_compare()`.
- **Article BB** — leave-one-block-out (LOBO) permutation test: `cf_test()`, with the restricted RV by
  default, schemes `rows`, `groups`, `within`, `variables`, `components`, comparison with the full RV and
  PROTEST (`compare=True`), and p-values by permutation, gamma approximation or hybrid (`pvalue=`).
- `cf_test_path()`: **experimental** choice of R (not part of the articles; flagged at each call).
- Plots (matplotlib): `fit.plot()`, `test.plot()`, `test.plot("null")`, `cf_plot_*()`, `cf_rose()`.
- Language: `cf_lang("fr")`. Help: `help(cf_test)` (English), `cf_help("cf_test", "fr")` (French).

## Installation

```
pip install ./python          # from the repository
```

Requires numpy, scipy, pandas, matplotlib.

## Example

```python
import canofac as cf
blocks = cf.cf_potato("potato.csv")                 # Chemical, Compression, Sensory (article BB)
raw = cf.cf_potato("potato.csv", blocks="A")        # NIRraw, FIDraw, CPMGraw (article A)
fit = cf.canofac(blocks, R=2, q=2)                  # CPCA
fit.summary()
test = cf.cf_test(blocks, R=2, q=2, B=999, seed=1)  # LOBO test, restricted RV
print(test)
test.plot(); test.plot("null")
cf.cf_lang("fr"); print(test)
```

The datasets are not distributed with canofac; export them from R with
`cf_write_blocks(cf_potato(preprocess = FALSE), "potato.csv")` (see `cf_help("donnees", "fr")`).

## Tests

```
python3 sync_shared.py      # copy the files shared with R (labels, French help, reference values)
python3 -m pytest tests
```

The reference values (`inst/shared/reference/reference_values.json`) are checked by the Python and the R
tests.

Author and maintainer: Mohamed Hanafi. Licence GPL-3.
