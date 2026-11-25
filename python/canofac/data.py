"""Reading and writing multiblock data, and loaders for the datasets of the articles.

Exchange format (same in R): one CSV file, first column the row names, then
the variables named ``<block>.x.<variable>`` (the layout of the potato data of
the R package multiblock); columns without ``.x.`` are extra variables (for
example a group). In R, ``cf_write_blocks(cf_friday(preprocess = FALSE), "friday.csv")``
exports the datasets distributed with ade4, mixOmics or multiblock for use in Python.
"""
import os
import numpy as np
import pandas as pd
from .preprocess import cf_preprocess

FRIDAY_SIZES = {"Hemiptera": 11, "Odonata": 7, "Trichoptera": 13, "Ephemeroptera": 4, "Coleoptera": 13,
                "Diptera": 22, "Hydracarina": 4, "Malacostraca": 3, "Mollusca": 8, "Oligochaeta": 6}


def cf_read_blocks(file, blocks=None):
    """Read blocks from a CSV file in the exchange format.

    Parameters
    ----------
    file : str
        CSV file: first column row names, variables ``<block>.x.<variable>``.
    blocks : list of str, optional
        Blocks to keep, in this order (default: all, in file order).

    Returns
    -------
    (dict, pandas.DataFrame)
        The blocks (name -> data frame) and the extra columns (may be empty).
    """
    d = pd.read_csv(file, index_col=0)
    names, cols, extra = [], {}, []
    for c in d.columns:
        if ".x." in c:
            b, v = c.split(".x.", 1)
            if b not in cols:
                names.append(b)
                cols[b] = []
            cols[b].append(c)
        else:
            extra.append(c)
    keep = blocks if blocks is not None else names
    missing = [b for b in keep if b not in cols]
    if missing:
        raise ValueError(f"blocks not found: {missing}")
    out = {}
    for b in keep:
        X = d[cols[b]].copy()
        X.columns = [c.split(".x.", 1)[1] for c in cols[b]]
        out[b] = X
    return out, d[extra]


def cf_write_blocks(blocks, file, extra=None):
    """Write blocks (dict name -> array or data frame) to a CSV file in the exchange format."""
    parts = []
    for b, X in blocks.items():
        X = X if isinstance(X, pd.DataFrame) else pd.DataFrame(np.asarray(X))
        Y = X.copy()
        Y.columns = [f"{b}.x.{c}" for c in X.columns]
        parts.append(Y.reset_index(drop=True))
    out = pd.concat(parts, axis=1)
    first = next(iter(blocks.values()))
    if isinstance(first, pd.DataFrame):
        out.index = first.index
    if extra is not None:
        for c in extra:
            out[c] = np.asarray(extra[c])
    out.to_csv(file)


POTATO_SETS = {
    "BB": ("Chemical", "Compression", "Sensory"),
    "A": ("NIRraw", "FIDraw", "CPMGraw"),
}


def cf_potato(file, blocks="BB", preprocess=True):
    """Potato data (Thygesen et al., 2001; R package multiblock): 26 samples, 9 blocks.

    ``blocks``: a list of block names, or ``"BB"`` (Chemical, Compression,
    Sensory: article on the LOBO test) or ``"A"`` (NIRraw, FIDraw, CPMGraw,
    measured on the raw tubers: article on pseudo-singular values). In the
    multiblock version, each NIR spectrum is SNV-transformed and the variables of
    Chemical and Compression are already scaled; preprocessing: centring and
    unit total inertia.
    """
    if isinstance(blocks, str):
        blocks = POTATO_SETS[blocks]
    bl, _ = cf_read_blocks(file, list(blocks))
    return cf_preprocess(bl) if preprocess else bl


def cf_friday(file, preprocess=True):
    """Friday (1987) data: 16 ponds, 91 species in 10 taxonomic blocks (article A).

    ``file``: exchange-format CSV (``cf_write_blocks`` in R from ``cf_friday(preprocess = FALSE)``),
    or the faunal table of ade4 ``friday87$fau`` as a CSV (16 rows, 91 columns in block order).
    Preprocessing: centring and unit total inertia.
    """
    d = pd.read_csv(file)
    if d.columns[0].startswith("Unnamed") or d.columns[0] == "" or d.iloc[:, 0].dtype == object:
        d = d.set_index(d.columns[0])
    if any(".x." in c for c in d.columns):
        bl, _ = cf_read_blocks(file)
    else:
        if d.shape[0] != 16 and d.shape[1] == 16:
            d = d.T
        bl, j = {}, 0
        for b, p in FRIDAY_SIZES.items():
            bl[b] = d.iloc[:, j:j + p]
            j += p
    return cf_preprocess(bl) if preprocess else bl


def cf_tcga(file, preprocess=True):
    """breast.TCGA training data (mixOmics): 150 tumours, mRNA, miRNA, proteins (article A).

    ``file``: exchange-format CSV with an extra column ``subtype``. Preprocessing:
    centring, scaling, unit total inertia. Returns (blocks, subtype).
    """
    bl, extra = cf_read_blocks(file)
    sub = extra["subtype"].to_numpy() if "subtype" in extra else None
    return (cf_preprocess(bl, scale=True) if preprocess else bl), sub
