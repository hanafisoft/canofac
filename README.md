# canofac — Canonical Factorisation for Multiblock Data

<!-- badges: start -->
[![R-CMD-check](https://github.com/hanafisoft/canofac/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hanafisoft/canofac/actions/workflows/R-CMD-check.yaml)
[![python](https://github.com/hanafisoft/canofac/actions/workflows/python.yaml/badge.svg)](https://github.com/hanafisoft/canofac/actions/workflows/python.yaml)
[![DOI](https://zenodo.org/badge/1472117148.svg)](https://doi.org/10.5281/zenodo.24540840)
<!-- badges: end -->

`canofac` estimates exploratory multiblock methods written in the canonical
factorisation \(X_k = U D_k V_k' + E_k\) (Bougeard et al., 2024), where \(U\)
holds the reference components (scores), shared by all blocks, \(V_k\) the axes
(loadings) of block \(k\) and \(D_k\) the pseudo-singular values \(d_{kh}\).

* **q-family**: maximisation of \(\sum_k \sum_h |d_{kh}|^q\), which includes
  MAXBET (q = 1), CPCA and MCOA (q = 2), ACCPS/CCSWA/ComDim and INDORT (q = 4),
  with unconstrained or orthonormal axes and sequential or simultaneous
  extraction (`canofac()`); STATIS and multiple factor analysis through the
  weighting of the blocks (`weighting = "statis"` or `"mfa"`).
* **Trajectories** of each block as q varies (`cf_trajectories()`).
* **Exact LOBO (leave-one-block-out) permutation test**, block by block
  (`cf_test()`): each block is compared with the reference components
  estimated on the other blocks, with a subspace statistic that weights the
  components by the inertia of the other blocks (restricted RV, default) or
  equally (`subspace = "equal"`); the test
  designates the blocks linked to the
  subspace and to each component (all, several, one or none); choice of the
  number of dimensions (`cf_test_path()`). Being linked is not a label: to
  label parts of the retained space as common, shared, private or residual,
  see [CoPRes](https://github.com/hanafisoft/CoPRes).
* **Multiple starts** for q > 2 (`nstart`), since the criterion can be
  multimodal.
* **Graphics module**: individuals on the reference components with an
  ellipse per group (`cf_plot_scores()`), table of the saliences
  (`cf_plot_saliences()`), annotated trajectories
  (`cf_plot_trajectories()`), p-values of the test with the linked blocks
  (`cf_plot_test()`), comparison of methods
  (`cf_plot_compare()`), rose plots (`cf_rose()`), loadings drawn as spectra.
* **Labels in English or French**: `cf_lang("fr")` switches the printed and
  plotted results to French, `cf_lang("en")` back to English.
* Default preprocessing: column centring and unit total inertia of each block
  (`cf_preprocess()`), with optional SNV and column scaling.

## Installation

```r
# install.packages("remotes")
remotes::install_github("hanafisoft/canofac", build_vignettes = TRUE)
```

The examples use the Friday data of the package **ade4**
(`install.packages("ade4")`) and the breast.TCGA data of the Bioconductor
package **mixOmics** (`BiocManager::install("mixOmics")`).

## Example

```r
library(canofac)
cf_lang("en")                               # or cf_lang("fr")
blocks <- cf_friday()                       # requires the ade4 package
fit <- canofac(blocks, R = 2, q = 4)        # INDORT
summary(fit)
cf_plot_saliences(fit)
test <- cf_test(blocks, R = 2, q = 4, B = 999, seed = 1)
test
cf_plot_test(test)
statis <- cf_test(blocks, R = 2, q = 2, B = 999, seed = 1,
                  weighting = "statis")
statis
```

## Step-by-step script

`canofac_step_by_step.R` runs the whole analysis, commented in French:
preprocessing, q-family, trajectories, multiple starts, choice of the number
of dimensions, test, crossing of the two permutation schemes, maps of the
individuals and export of the figures.

```r
file.edit(system.file("scripts", "canofac_step_by_step.R", package = "canofac"))
```

`demo_friday.R`, in the same folder, gives a shorter analysis of the Friday
data.

## Data used in the articles

| Data | Access | Loader |
|---|---|---|
| potato (Thygesen et al., 2001) | package **multiblock**, or a CSV file | `cf_potato()` (article BB), `cf_potato(blocks = "A")` (article A) |
| breast.TCGA | package **mixOmics** | `cf_tcga()` |
| Friday (1987) | package **ade4** | `cf_friday()` |

## Citation

Hanafi M. (2026). canofac: Canonical Factorisation for Multiblock Data (R and
Python). Zenodo. <https://doi.org/10.5281/zenodo.24540840> — in R: `citation("canofac")`.
The article presenting the methods is in preparation;
please also cite Bougeard et al. (2024), *Chemometrics and Intelligent
Laboratory Systems*, 254, 105240.

## Author and contact

Mohamed Hanafi (author and maintainer).

Contact: mohamed.hanafi@oniris-nantes.fr, mohammed44hanafi@gmail.com.
Questions and bug reports: <https://github.com/hanafisoft/canofac/issues>.

## Licence

GPL-3
