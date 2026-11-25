"""canofac: canonical factorisation for multiblock data (Python version).

Exploratory multiblock methods written in the canonical factorisation
X_k = U D_k V_k' (Bougeard et al., 2024): the power-q family of criteria on the
pseudo-singular values (MAXBET, CPCA, MCOA, ACCPS/ComDim, INDORT; STATIS and
MFA through the weighting of the blocks) - article A - and the
leave-one-block-out (LOBO) permutation test - article BB. Same functions,
arguments, results, labels (English or French) and random draws as the R
package canofac.

Main functions: canofac, cf_preprocess, cf_trajectories, cf_multistart,
cf_angles, cf_compare, cf_test, cf_test_path (experimental), cf_lang, cf_help,
plots cf_plot_* and cf_rose, data cf_read_blocks / cf_write_blocks and loaders.

Author and maintainer: Mohamed Hanafi. Licence GPL-3.
Help in French: cf_help("canofac", "fr").
"""
from .lang import cf_lang
from .preprocess import cf_preprocess, snv
from .fit import canofac, Canofac
from .explore import cf_trajectories, cf_multistart, cf_angles, cf_compare, method_name
from .test import cf_test, CfTest, holm
from .path import cf_test_path
from .moments import perm_moments, gamma_pvalue
from .data import cf_read_blocks, cf_write_blocks, cf_potato, cf_friday, cf_tcga
from .help import cf_help
from .simulate import cf_simulate
from ._rng import Generator as cf_generator, cf_permutations
from .graphics import (cf_plot_scores, cf_plot_saliences, cf_plot_trajectories, cf_plot_test,
                       cf_plot_compare, cf_rose)

__version__ = "0.7.0"
__all__ = ["cf_lang", "cf_preprocess", "snv", "canofac", "Canofac", "cf_trajectories", "cf_multistart",
           "cf_angles", "cf_compare", "method_name", "cf_test", "CfTest", "holm", "cf_test_path",
           "perm_moments", "gamma_pvalue", "cf_read_blocks", "cf_write_blocks", "cf_potato", "cf_friday",
           "cf_tcga", "cf_help", "cf_simulate", "cf_generator", "cf_permutations", "cf_plot_scores", "cf_plot_saliences",
           "cf_plot_trajectories", "cf_plot_test", "cf_plot_compare", "cf_rose"]
