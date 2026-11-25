"""Simulated multiblock data of the article on the LOBO test (scenarios S0, SA to SE)."""
import numpy as np
from ._rng import Generator
from ._core import polar
from .preprocess import cf_preprocess

SCENARIOS = ("S0", "SA", "SB", "SC", "SD", "SE")


def cf_simulate(scenario="SA", n=20, p=10, K=4, theta=1.6, theta2=1.5, seed=1, preprocess=True):
    """Simulate K blocks with the scenarios of the article on the LOBO test.

    Two orthonormal centred directions u0, u1 of the individuals; in each block, noise with
    eigenvalues proportional to 1/j (j = 1, ..., p) in a random orthonormal basis, plus signal:

    - S0: no signal (H0);
    - SA: u0 in all the blocks (intensity theta);
    - SB: u0 in all the blocks, u1 in block 1 only (theta);
    - SC: u0 in blocks 1 and 2 only;
    - SD: u0 in all the blocks, u1 in blocks 1 and 2 (intensity theta2);
    - SE: u0 in blocks 1 and 2, u1 in all the blocks (theta2).

    The intensities are multiplied by sqrt(n). All the draws come from the common
    generator, and the orthonormal bases are polar factors (no QR sign convention):
    the same seed gives the same blocks in R and in Python.

    Returns
    -------
    dict
        Blocks ``Block1``, ..., ``BlockK`` (centred, unit total inertia if ``preprocess``).
    """
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    g = Generator(seed)
    G = g.normal(n * 2).reshape((n, 2), order="F")
    U = polar(G - G.mean(axis=0))
    u0, u1 = U[:, 0], U[:, 1]
    lam = 1.0 / np.arange(1, p + 1)
    lam = lam * p / lam.sum()
    t1, t2 = theta * np.sqrt(n), theta2 * np.sqrt(n)
    out = {}
    for k in range(K):
        Q = polar(g.normal(p * p).reshape((p, p), order="F"))
        E = (g.normal(n * p).reshape((n, p), order="F") * np.sqrt(lam)) @ Q.T
        S = np.zeros((n, p))
        a = g.normal(p)
        a /= np.linalg.norm(a)
        if scenario in ("SA", "SB", "SD") or (scenario in ("SC", "SE") and k < 2):
            S += t1 * np.outer(u0, a)
        if scenario in ("SB", "SD", "SE"):
            b = g.normal(p)
            b /= np.linalg.norm(b)
            if scenario == "SB" and k == 0:
                S += t1 * np.outer(u1, b)
            if scenario == "SD" and k < 2:
                S += t2 * np.outer(u1, b)
            if scenario == "SE":
                S += t2 * np.outer(u1, b)
        out[f"Block{k + 1}"] = S + E
    return cf_preprocess(out) if preprocess else out
