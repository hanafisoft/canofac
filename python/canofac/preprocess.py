"""Preprocessing of the blocks."""
import numpy as np


def snv(X):
    """Standard normal variate: centre and scale each row (spectrum) by its own mean and standard deviation."""
    X = np.asarray(X, dtype=float)
    m = X.mean(axis=1, keepdims=True)
    s = X.std(axis=1, ddof=1, keepdims=True)
    return (X - m) / s


def _names_of(blocks):
    if isinstance(blocks, dict):
        return list(blocks.keys()), list(blocks.values())
    blocks = list(blocks)
    return [f"Block{k + 1}" for k in range(len(blocks))], blocks


def cf_preprocess(blocks, center=True, scale=False, snv=False, unit_inertia=True):
    """Preprocess the blocks: SNV, column centring, scaling, unit total inertia.

    Parameters
    ----------
    blocks : dict of name -> array, list of arrays, or a single array
        Blocks with the same number of rows (individuals). Pandas data frames
        are accepted.
    center : bool
        Centre each column.
    scale : bool or list of bool
        Scale each column by its standard deviation (one value per block or
        a single value), for variables measured in different units.
    snv : bool or list of bool
        Apply the standard normal variate to the rows of the block (spectra).
    unit_inertia : bool
        Bring each block to unit total inertia, ``||X_k||_F^2 = 1``.

    Returns
    -------
    dict
        The preprocessed blocks, as a dict name -> array, or data frame when
        the block was a data frame (block names kept from a dict, else
        ``Block1``, ``Block2``, ...).
    """
    if not isinstance(blocks, (dict, list, tuple)):
        blocks = [blocks]
    names, arrs = _names_of(blocks)
    K = len(arrs)
    scale = list(scale) if isinstance(scale, (list, tuple)) else [scale] * K
    use_snv = list(snv) if isinstance(snv, (list, tuple)) else [snv] * K
    ns = {np.asarray(X).shape[0] for X in arrs}
    if len(ns) != 1:
        raise ValueError("all blocks must have the same number of rows")
    out = {}
    for k, (nm, X) in enumerate(zip(names, arrs)):
        X = np.asarray(X, dtype=float).copy()
        if X.ndim == 1:
            X = X[:, None]
        if np.isnan(X).any():
            raise ValueError(f"block {nm} contains missing values")
        if use_snv[k]:
            X = globals()["snv"](X)
        if center:
            X = X - X.mean(axis=0)
        if scale[k]:
            s = X.std(axis=0, ddof=1)
            s[s == 0] = 1.0
            X = X / s
        if unit_inertia:
            X = X / np.sqrt(np.sum(X ** 2))
        if hasattr(arrs[k], "columns"):                  # keep the names of a data frame
            import pandas as pd
            X = pd.DataFrame(X, index=arrs[k].index, columns=arrs[k].columns)
        out[nm] = X
    return out
