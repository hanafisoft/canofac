"""Exact permutation moments of T = tr(A P B P') and gamma approximation of its p-value.

For doubly centred symmetric n x n matrices A and B and a uniformly random
permutation matrix P (Kazi-Aoual et al., 1995; Josse et al., 2008):

    E(T)   = tr(A) tr(B) / (n - 1),
    Var(T) = explicit function of n and, for A and B, of the trace, the trace of
             the square and the sum of the squared diagonal terms.

The gamma approximation fits a gamma law with these two moments. In the
simulations of article BB it gives the same decisions as the permutation test
in 99 % of the cases but is slightly liberal (rejection rate 5.2 % instead of
4.8 % under H0); the hybrid mode permutes only near the thresholds.
"""
import numpy as np
from scipy import stats


def _s(M):
    t, t2, s2 = np.trace(M), np.sum(M * M), np.sum(np.diag(M) ** 2)
    return t, t2, t * t / t2, s2 / t2


def perm_moments(A, B):
    """Mean and variance of tr(A P B P') over the permutations P (A, B doubly centred)."""
    n = A.shape[0]
    ta, ta2, ba, da = _s(A)
    tb, tb2, bb, db = _s(B)
    m = ta * tb / (n - 1)
    ca = n * (n + 1) * da - (n - 1) * (ba + 2)
    cb = n * (n + 1) * db - (n - 1) * (bb + 2)
    v = 2 * (n - 1 - ba) * (n - 1 - bb) / ((n + 1) * (n - 1) ** 2 * (n - 2)) \
        * (1 + (n - 3) / (2 * n * (n - 1)) * ca * cb / ((n - 1 - ba) * (n - 1 - bb)))
    return m, v * ta2 * tb2


def gamma_pvalue(obs, A, B):
    """Upper-tail p-value of tr(A B) under the gamma law with the exact permutation moments."""
    if np.trace(B) <= 0:
        return 1.0
    m, v = perm_moments(A, B)
    if not np.isfinite(v) or v <= 0 or m <= 0:
        return 1.0
    return float(stats.gamma.sf(obs, m * m / v, scale=v / m))


def is_centred(X, tol=1e-8):
    return np.max(np.abs(X.mean(axis=0))) <= tol * max(1.0, np.max(np.abs(X)))
