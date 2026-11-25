"""Estimation of the q-family of criteria in canonical factorisation."""
import numpy as np
import pandas as pd
from ._core import cf_fit, weight_blocks, nstart_default
from .lang import _tr
from .preprocess import _names_of


def _as_blocks(blocks):
    """dict name -> float array, from a dict, list, tuple or single array (data frames accepted)."""
    if not isinstance(blocks, (dict, list, tuple)):
        blocks = [blocks]
    names, arrs = _names_of(blocks)
    cols = [list(X.columns) if hasattr(X, "columns") else None for X in arrs]
    rows = arrs[0].index.tolist() if hasattr(arrs[0], "index") else None
    arrs = [np.asarray(X, dtype=float) for X in arrs]
    n = {X.shape[0] for X in arrs}
    if len(n) != 1:
        raise ValueError("all blocks must have the same number of rows")
    return names, arrs, cols, rows


class Canofac:
    """Result of :func:`canofac`.

    Attributes
    ----------
    U : ndarray (n, R)
        Reference components (scores), orthonormal.
    V : dict name -> ndarray (p_k, R)
        Axes (loadings) of each block, of unit norm.
    d : ndarray (K, R)
        Pseudo-singular values d_kh.
    saliences : ndarray (K, R)
        Saliences d_kh^2.
    delta : ndarray (K, R)
        d_kh / sigma_1(X_k).
    sigma1 : ndarray (K,)
        First singular value of each (weighted) block.
    criterion : float
        Value of sum_k sum_h |d_kh|^q.
    block_weights : ndarray (K,)
        Relative weights of the blocks (STATIS, MFA; equal otherwise).
    block_names : list of str
    settings : dict
        R, q, axes, extraction, nstart, weighting.
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def table(self, what="d"):
        """K x R table (pandas DataFrame) of ``"d"``, ``"saliences"`` or ``"delta"``."""
        M = getattr(self, what)
        return pd.DataFrame(M, index=self.block_names,
                            columns=[f"Dim{h + 1}" for h in range(M.shape[1])])

    def __repr__(self):
        s = self.settings
        out = _tr("title_fit") + "\n"
        out += _tr("settings_fit") % (s["q"], s["axes"], s["extraction"], s["R"], len(self.block_names))
        if s.get("weighting", "none") != "none":
            out += (_tr("weighting_line") % s["weighting"].upper()) + "\n"
        out += _tr("criterion") % f"{self.criterion:.6g}"
        out += _tr("pvs_head") + "\n" + self.table("d").round(3).to_string()
        return out

    def summary(self):
        """Print the pseudo-singular values, saliences and delta_kh; return them as data frames."""
        print(repr(self))
        print("\n" + _tr("sal_head"))
        print(self.table("saliences").round(3).to_string())
        print("\n" + _tr("delta_head"))
        print(self.table("delta").round(3).to_string())
        return {w: self.table(w) for w in ("d", "saliences", "delta")}

    def plot(self, type="scores", **kw):
        """Plot the fit: ``type`` = "scores", "d", "saliences", "delta" or "loadings"."""
        from . import graphics
        return graphics.plot_fit(self, type=type, **kw)


def canofac(blocks, R=2, q=2, axes="free", extraction="simultaneous", tol=1e-10, maxit=5000,
            nstart=None, weighting="none"):
    """Canonical factorisation by the power-q family of criteria.

    Fits ``X_k = U D_k V_k' + E_k`` by maximising ``sum_k sum_h |d_kh|^q``,
    where ``d_kh`` is the pseudo-singular value of block k on dimension h.
    q = 1 gives MAXBET (sequential, unconstrained axes) or simultaneous MAXBET
    (orthonormal axes); q = 2 gives CPCA (unconstrained axes) or MCOA
    (orthonormal axes, sequential or simultaneous); q = 4 gives ACCPS
    (sequential) or INDORT (simultaneous). STATIS and multiple factor
    analysis are CPCA on weighted blocks (``weighting``).

    Parameters
    ----------
    blocks : dict, list or array
        Blocks with the same rows, preprocessed (see :func:`cf_preprocess`).
    R : int
        Number of dimensions.
    q : float
        Power applied to the pseudo-singular values, q >= 1.
    axes : {"free", "orthonormal"}
        Unconstrained axes of unit norm, or orthonormal axes in each block.
    extraction : {"simultaneous", "sequential"}
        Simultaneous fit or one dimension at a time with deflation.
    tol, maxit : float, int
        Relative tolerance on the criterion and maximum number of iterations.
    nstart : int, optional
        Number of starts for unconstrained axes: the deterministic start
        (CPCA) and ``nstart - 1`` random starts drawn by the common generator
        (seed 12345, same starts in R and Python); best criterion kept.
        Default 10 if q > 2 with unconstrained axes, else 1.
    weighting : {"none", "mfa", "statis"}
        Weighting of the blocks before the fit.

    Returns
    -------
    Canofac
        See :class:`Canofac`. French help: ``cf_help("canofac", "fr")``.

    References
    ----------
    Hanafi M., Bougeard S., Galharret J.-M., Van Deun K. Pseudo-singular values:
    a family of power criteria for multiblock dimension reduction (article A).
    """
    if axes not in ("free", "orthonormal"):
        raise ValueError("axes must be 'free' or 'orthonormal'")
    if extraction not in ("simultaneous", "sequential"):
        raise ValueError("extraction must be 'simultaneous' or 'sequential'")
    if q < 1:
        raise ValueError("q must be >= 1")
    names, arrs, cols, rows = _as_blocks(blocks)
    n = arrs[0].shape[0]
    if R > n - 1:
        raise ValueError("R must be smaller than n")
    nstart = nstart_default(nstart, q, axes)
    wb, w = weight_blocks(arrs, weighting)
    res = cf_fit(wb, R, q, axes, extraction, tol, maxit, nstart)
    sigma1 = np.array([np.linalg.svd(X, compute_uv=False)[0] for X in wb])
    d = res["d"]
    V = {nm: Vk for nm, Vk in zip(names, res["V"])}
    return Canofac(U=res["U"], V=V, d=d, saliences=d ** 2, delta=d / sigma1[:, None],
                   sigma1=sigma1, criterion=float(np.sum(np.abs(d) ** q)), block_weights=w,
                   block_names=names, variable_names=dict(zip(names, cols)), row_names=rows,
                   settings=dict(R=R, q=q, axes=axes, extraction=extraction, nstart=nstart,
                                 weighting=weighting))
