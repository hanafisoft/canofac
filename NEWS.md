# canofac 0.7.0

Version covering the two articles (pseudo-singular values and the power-q family;
LOBO permutation test), released together for R and Python.

* **Python version** (directory `python/`): same functions, arguments, results,
  labels and random draws as the R package; matplotlib plots; pytest tests.
* **Common random number generator** (`cf_permutations()`): Park-Miller MINSTD
  generator and Fisher-Yates permutations, identical in R and Python. With the
  same `seed`, `cf_test()` gives the same permutations, hence the same p-values,
  in both languages; the random starts of `canofac()` (q > 2) are the same too.
  R's global generator is no longer used (only to draw a seed when none is
  given; the seed is kept in `settings$seed`). **The p-values of a given seed
  differ from those of canofac 0.6.0**, and so may the random starts.
* `cf_test()`:
  - `pvalue = "gamma"` (gamma approximation with the exact permutation mean and
    variance, no permutation; slightly liberal) and `pvalue = "hybrid"`
    (permutation only when the gamma p-value lies between alpha/(c K) and
    c alpha, `hybrid_c = 3`); `pmethod` records the method of each p-value.
  - schemes `"groups"` (whole groups of rows, e.g. replicates) and `"within"`
    (rows within groups), argument `groups`.
  - `compare = TRUE`: full RV test and PROTEST with the same permutations.
  - `null`: permutation distributions of the subspace statistic;
    `plot(x, type = "null")` and `cf_plot_null()`.
* New functions for the article on pseudo-singular values: `cf_multistart()`
  (multimodality), `cf_angles()` (angles between two fits), `cf_compare()`
  (pseudo-singular values of several methods).
* `perm_moments()` and `gamma_pvalue()`: exact permutation moments of
  tr(A P B P').
* Data: one of the former data loaders is removed, with every reference to its
  dataset. `cf_potato()` (potato data of the package multiblock;
  `blocks = "BB"`, the default, for the article on the LOBO test, `blocks = "A"`
  for the raw NIR, FID and CPMG blocks of the article on pseudo-singular values); `cf_read_blocks()` and `cf_write_blocks()`: CSV exchange format
  shared with Python.
* `cf_test_path()` (choice of R) is **experimental**: a warning at each call,
  and the printed result and the plot say so; it now accepts `groups`, and
  model r uses the seed `seed + r - 1`.
* Help in two languages: `cf_help(topic, "fr")` (French pages of
  `inst/shared/aide/fr`, shared with Python) and `cf_help(topic, "en")`;
  French vignette "Prise en main de canofac".
* Labels read from `inst/shared/labels.csv`, shared with Python.
* Reference values shared with Python (`inst/shared/reference`), checked by the
  tests of both versions (`tests/testthat/test-crosscheck.R`).
* Scripts reproducing the articles: `inst/scripts/articleA.R`,
  `inst/scripts/articleBB.R` (and their Python versions in `python/scripts`).

# canofac 0.6.0

* The fixed-reference test is renamed the LOBO test (leave-one-block-out):
  the reference components used to test a block are estimated without it.
  `statistic = "lobo"` is the new default name; `statistic = "fixed"` remains
  accepted as an alias. Printed labels and documentation are updated.

# canofac 0.5.2

* New permutation scheme `scheme = "components"` in `cf_test()` and
  `cf_test_path()`: the scores of each principal component of the tested block
  are permuted independently, and the block is rebuilt with its loadings.
  Invariant under rotation of the variables, it preserves the spectrum and the
  loadings of the block; it is an alternative to variable permutation for
  ordered and correlated variables such as spectra.

# canofac 0.5.1

* The restricted RV statistic T_k^w (`subspace = "inertia"`) becomes the
  default subspace statistic of `cf_test()` and `cf_test_path()`: same
  subspace, null hypothesis and permutations as T_k, exact, at least as
  powerful and clearly more so for small samples. Use `subspace = "equal"`
  for the equal-weight statistic of canofac 0.4.0 and 0.5.0. With
  `statistic = "refit"`, the subspace statistic is the sum of the squared
  pseudo-singular values (equal weights), as before.

# canofac 0.5.0

* STATIS and multiple factor analysis (MFA): new argument `weighting`
  (`"none"`, `"mfa"`, `"statis"`) in `canofac()`, `cf_test()` and
  `cf_test_path()`. MFA divides each block by its first singular value;
  STATIS brings each block to ||X_k X_k'||_F = 1 and weights its scalar
  product matrix by the first eigenvector of C_kl = tr(W_k W_l). Use them with
  `q = 2`, `axes = "free"`, `extraction = "simultaneous"`. In `cf_test()`, the
  components of the other blocks are those of the weighted method applied to
  the K - 1 other blocks (STATIS weights recomputed). `canofac()` returns
  `block_weights`.
* New subspace statistic for the fixed-reference test: `subspace = "inertia"`,
  T_k^w = tr(W_k P W_(-k) P), which weights the components of the other
  blocks by their inertia (numerator of the RV coefficient between the block
  and the other blocks projected on the same components). Same subspace,
  null hypothesis and permutations as T_k (`subspace = "equal"`, default),
  exact as well, more powerful for small samples. `cf_test()` returns both
  subspace p-values, `p_subspace_equal` and `p_subspace_inertia`;
  `p_subspace` and the linked blocks use the chosen statistic.
* `print()` shows the subspace statistic and the weighting.

# canofac 0.4.0

* New default statistic of `cf_test()` and `cf_test_path()`:
  `statistic = "fixed"`, the fixed-reference statistic proposed by
  J.-M. Galharret. For each block k, the reference components U^(-k) are
  estimated once on the other blocks; the statistics are
  S_kh = ||X_k' u_h^(-k)|| and T_k = ||X_k' U^(-k)||_F, and only block k is
  permuted, without re-estimation. The test is exact, needs no matching and
  is much faster; unlike the re-estimation statistic, it keeps its power for
  a second link between the blocks.
* Vocabulary of links (breaking change): the test no longer classifies the
  dimensions as common, partially common, private or residual, since being
  linked to a component is not a label (a component may mix directions that
  the blocks carry differently; labeling parts of the retained space is done
  by CoPRes, <https://github.com/hanafisoft/CoPRes>). Elements of `cf_test`
  objects: `carriers` becomes `linked`, `structure` becomes `link` (values
  `"all"`, `"several"`, `"one"`, `"none"`), `common` becomes `all_linked`,
  and the subspace is described by `linked_subspace` and `link_subspace`.
  Columns of the profile of `cf_test_path()`: `common`, `partially_common`,
  `private` become `all`, `several`, `one`.
* The linked blocks are given for the subspace; the dimension-wise
  p-values bear on the components of the other blocks and are descriptive.
  `retained` means that the component of the other blocks is linked to at
  least one block (min_k p <= alpha / K); `cf_test_path()` chooses R with this
  rule.
* The previous statistic remains available with `statistic = "refit"`: exact
  test of the global null hypothesis, used for the subspace and, with
  `scheme = "variables"`, to recognise directions carried by one block
  without link to the others.
* `print()` and `cf_plot_test()` show the blocks linked to the subspace.

# canofac 0.3.0

* Categories of the test renamed: `"specific"` becomes `"private"` (French:
  propre) and `"noise"` becomes `"residual"` (French: résidus); the column
  `specific` of the profile of `cf_test_path()` becomes `private`.
* Repository: <https://github.com/hanafisoft/canofac>; author and
  maintainer Mohamed Hanafi;
  `inst/CITATION`; R CMD check on Linux, Windows and macOS with GitHub Actions.
* Language of the results: `cf_lang("fr")` or `cf_lang("en")` (option
  `canofac.lang`) sets the language of the printed results, plot titles, axis
  labels, legends and names of the categories of the test. The values stored
  in the objects are not translated.
* Graphics module: `cf_plot_scores()` (individuals on the reference components,
  normal ellipse per group), `cf_plot_saliences()` (table of the saliences,
  pseudo-singular values or delta), `cf_plot_trajectories()` (markers at
  q = 1, 2, 4 and at the largest q, with the names of the methods),
  `cf_plot_test()` (p-values, carrier blocks in bold, category of each
  dimension, subspace test) and `cf_plot_compare()` (pseudo-singular values of
  several fits). `plot()` methods now call these functions; `plot()` of a
  `cf_test` object shows the table of p-values instead of a bar plot.
* Matching of the dimensions: optimal assignment (Hungarian algorithm), exact
  for any number of dimensions; previously exhaustive search up to 7
  dimensions and greedy assignment beyond.
* `cf_trajectories()` stores its settings (axes, extraction).
* New script `inst/scripts/canofac_step_by_step.R`.
* Documentation of the level of the test with matching updated from the new
  simulation study.

# canofac 0.2.0

* Retention rule of `cf_test()`: a dimension is now retained when it is a
  common structure (max_k p <= alpha) or has at least one carrier block
  (min_k p <= alpha / K), i.e. when it is not noise. Previously only the second
  condition was used, so that a dimension common to all blocks with all
  p-values between alpha / K and alpha was not retained. `cf_test_path()`
  chooses R with this rule, and its plot shows the retained dimensions in bold.
* Vocabulary (breaking change): structure labels are now `"common"`,
  `"partially_common"`, `"specific"`, `"noise"`; element `global` of
  `cf_test` objects is renamed `common`; columns `global` and `partial` of
  `cf_test_path` profiles are renamed `common` and `partially_common`. The
  components U are called reference components in the documentation.

# canofac 0.1.0

* `nstart`: multiple starts for unconstrained axes (default 10 when q > 2),
  since the criterion can be multimodal from q = 4 (e.g. breast.TCGA).
* `cf_test()`: documentation of the effect of dimension matching on the level.
* Structure labels renamed `"global"`, `"partial"`, `"specific"`, `"noise"`;
  element `common` of `cf_test` objects renamed `global`, column `common` of
  `cf_test_path` profiles renamed `global`.

* First version: q-family estimation, trajectories, permutation test,
  preprocessing and data loaders.

* plot() methods for canofac, cf_trajectories and cf_test objects.
* plot(fit, type = "loadings"): axes of each block, as spectra when the
  column names are wavelengths or wavenumbers.
* cf_rose(): rose plot of the pseudo-singular values (inner: squared
  deviations, middle: values, outer ring: mean / variance composition).
* cf_test_path(): tests of the models with r = 1, ..., Rmax dimensions,
  Rmax = min(K + c, n - 1) by default, and choice of R as the largest model
  whose dimensions are all retained.
* cf_test() describes each dimension as common, partial, specific or none,
  with the blocks carrying it (Holm adjustment over the blocks).
