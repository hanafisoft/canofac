"""Exploration of the q-family (article A): trajectories, multiple starts, angles, comparison of methods."""
import numpy as np
import pandas as pd
from ._core import cf_fit, free_best, match_dims, nstart_default, weight_blocks
from .fit import _as_blocks
from .lang import _tr


class CfTrajectories:
    """Result of :func:`cf_trajectories`: ``q``, ``d`` and ``delta`` (K x R x len(q)), ``sigma1``, ``criterion``."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __repr__(self):
        lines = [_tr("traj_title")]
        for i, q in enumerate(self.q):
            df = pd.DataFrame(self.delta[:, :, i], index=self.block_names,
                              columns=[f"Dim{h + 1}" for h in range(self.delta.shape[1])])
            lines.append(f"q = {q:g}\n" + df.round(2).to_string())
        return "\n".join(lines)

    def plot(self, **kw):
        from . import graphics
        return graphics.plot_trajectories(self, **kw)


def cf_trajectories(blocks, R=2, qgrid=(1, 1.5, 1.9, 2, 2.1, 3, 4, 6, 8, 10), axes="free",
                    extraction="simultaneous", tol=1e-10, maxit=5000, nstart=None):
    """Pseudo-singular values of each block as q varies.

    The dimensions are matched from one value of q to the next, starting from
    the largest q (optimal assignment), so that each block can be followed in
    the plane (delta_k1, delta_k2).

    Parameters
    ----------
    blocks : dict, list or array
        Preprocessed blocks.
    R : int
        Number of dimensions.
    qgrid : sequence of float
        Values of q (>= 1).
    axes, extraction, tol, maxit, nstart
        As in :func:`canofac`.

    Returns
    -------
    CfTrajectories
        ``q``, ``d``, ``delta`` (arrays K x R x len(qgrid)), ``sigma1``, ``criterion``.
    """
    qgrid = sorted(set(float(x) for x in qgrid))
    if min(qgrid) < 1:
        raise ValueError("all values of qgrid must be >= 1")
    names, arrs, _, _ = _as_blocks(blocks)
    K = len(arrs)
    sigma1 = np.array([np.linalg.svd(X, compute_uv=False)[0] for X in arrs])
    d = np.full((K, R, len(qgrid)), np.nan)
    crit = np.zeros(len(qgrid))
    Uref = None
    for i in reversed(range(len(qgrid))):
        fit = cf_fit(arrs, R, qgrid[i], axes, extraction, tol, maxit, nstart_default(nstart, qgrid[i], axes))
        if Uref is None:
            Uref, pi = fit["U"], np.arange(R)
        else:
            pi = match_dims(fit["U"], Uref)
        d[:, :, i] = fit["d"][:, pi]
        crit[i] = np.sum(np.abs(fit["d"]) ** qgrid[i])
    return CfTrajectories(q=np.array(qgrid), d=d, delta=d / sigma1[:, None, None], sigma1=sigma1,
                          criterion=crit, block_names=names,
                          settings=dict(R=R, axes=axes, extraction=extraction))


class CfMultistart:
    """Result of :func:`cf_multistart`."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __repr__(self):
        return (_tr("multistart_title") + "\n" +
                _tr("multistart_head") % (self.share_best, self.n_optima,
                                          _tr("yes") if self.deterministic_best else _tr("no")))

    def plot(self, **kw):
        from . import graphics
        return graphics.plot_multistart(self, **kw)


def cf_multistart(blocks, R=2, q=4, nstart=100, extraction="simultaneous", seed=2026, rtol=1e-5,
                  tol=1e-10, maxit=5000):
    """Multimodality of the criterion: criterion reached from random starts (unconstrained axes).

    Runs the algorithm from ``nstart`` random centred orthonormal starts drawn by
    the common generator, and from the deterministic CPCA start, and counts the
    starts that reach the best criterion to within ``rtol`` in relative terms
    (article A, multimodality study).

    Returns
    -------
    CfMultistart
        ``criteria`` (random starts), ``deterministic`` (CPCA start), ``best``,
        ``share_best`` (%), ``n_optima`` (distinct values at ``rtol``),
        ``deterministic_best`` (bool).
    """
    from ._core import random_starts, free_simultaneous
    names, arrs, _, _ = _as_blocks(blocks)
    n = arrs[0].shape[0]

    def run(U0, Xs, r):
        if U0 is None:
            return np.sum(free_simultaneous(Xs, r, q, tol, maxit)[1] ** q)
        return np.sum(free_simultaneous(Xs, r, q, tol, maxit, U0=U0)[1] ** q)

    starts = random_starts(n, R if extraction == "simultaneous" else 1, nstart + 1, seed=seed)
    if extraction == "simultaneous":
        crits = np.array([run(U0, arrs, R) for U0 in starts])
        det = run(None, arrs, R)
    else:                                              # sequential: same start for each step
        def seq(U0):
            Xc = [X.copy() for X in arrs]
            tot = 0.0
            for h in range(R):
                if U0 is None:
                    U, D, _ = free_simultaneous(Xc, 1, q, tol, maxit)
                else:
                    U, D, _ = free_simultaneous(Xc, 1, q, tol, maxit, U0=U0)
                tot += np.sum(D ** q)
                u = U[:, 0]
                Xc = [X - np.outer(u, u @ X) for X in Xc]
            return tot
        crits = np.array([seq(U0) for U0 in starts])
        det = seq(None)
    best = max(crits.max(), det)
    reach = crits >= best * (1 - rtol)
    vals = np.sort(crits)[::-1]
    optima = [vals[0]]
    for v in vals[1:]:
        if abs(v - optima[-1]) > rtol * abs(best):
            optima.append(v)
    return CfMultistart(criteria=crits, deterministic=det, best=best, share_best=100 * reach.mean(),
                        n_optima=len(optima), optima=np.array(optima),
                        deterministic_best=bool(det >= best * (1 - rtol)),
                        settings=dict(R=R, q=q, nstart=nstart, extraction=extraction, seed=seed))


def cf_angles(fit1, fit2):
    """Angles (degrees) between the matched components of two fits, and largest principal angle between their subspaces.

    Returns
    -------
    dict
        ``components`` (array of R angles), ``subspace`` (float).
    """
    U1, U2 = fit1.U, fit2.U
    pi = match_dims(U2, U1)
    c = np.abs(np.sum(U1 * U2[:, pi], axis=0)).clip(0, 1)
    comp = np.degrees(np.arccos(c))
    s = np.linalg.svd(U1.T @ U2, compute_uv=False).clip(0, 1)
    sub = float(np.degrees(np.arccos(s.min())))
    print(_tr("angles_head"))
    print(pd.Series(comp.round(1), index=[f"Dim{h + 1}" for h in range(len(comp))]).to_string())
    print(_tr("subspace_angle") % sub, end="")
    return {"components": comp, "subspace": sub}


def method_name(q, axes, extraction):
    """Name of the classical method for (q, axes, extraction), or ''."""
    if axes == "free":
        if q == 1 and extraction == "sequential":
            return "MAXBET"
        if q == 2:
            return "CPCA"
        if q == 4:
            return "ACCPS" if extraction == "sequential" else "INDORT"
    else:
        if q == 1 and extraction == "simultaneous":
            return "simultaneous MAXBET"
        if q == 2:
            return "MCOA" if extraction == "sequential" else "simultaneous MCOA"
    return ""


def cf_compare(fits, names=None, digits=2, match=True):
    """Pseudo-singular values of several fits (article A, tables by method).

    With ``match=True`` the dimensions are matched to those of the first fit; with ``match=False`` each fit
    keeps its own order (decreasing contribution to its criterion), as in the tables of the article.

    Returns
    -------
    pandas.DataFrame
        One row per block, one column per method and dimension (``"CPCA · Dim1"``...).
    """
    if names is None:
        names = []
        for f in fits:
            m = method_name(f.settings["q"], f.settings["axes"], f.settings["extraction"])
            names.append(m if m else f"q = {f.settings['q']:g}")
    ref = fits[0].U
    cols = {}
    for nm, f in zip(names, fits):
        pi = match_dims(f.U, ref) if match else np.arange(ref.shape[1])
        for h in range(ref.shape[1]):
            cols[f"{nm} · Dim{h + 1}"] = f.d[:, pi[h]]
    df = pd.DataFrame(cols, index=fits[0].block_names)
    print(_tr("compare_table_head"))
    print(df.round(digits).to_string())
    return df
