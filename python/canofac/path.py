"""EXPERIMENTAL: choice of the number of dimensions R by successive tests.

The choice of R is not part of the published articles (A and BB). For each
r = 1, ..., Rmax the model is re-estimated and tested; the models are not
nested for simultaneous methods, and the p-values of the chosen model are
optimistic, since the choice uses the same data as the tests.
"""
import warnings
import numpy as np
import pandas as pd
from .test import cf_test
from .lang import _tr
from .fit import _as_blocks


class CfTestPath:
    """Result of :func:`cf_test_path` (experimental): ``tests``, ``profile``, ``R``, ``Rmax``."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __repr__(self):
        s = self.settings
        out = [_tr("experimental_title"), _tr("experimental_note").rstrip("\n"), "", _tr("title_path"),
               (_tr("settings_test") % (s["q"], s["axes"], s["extraction"], s["scheme"], s["B"],
                                         s["alpha"])).rstrip("\n"),
               f"  Rmax = {self.Rmax}", "", self.profile.to_string(index=False),
               (_tr("chosen_R") % self.R).rstrip("\n"), "", _tr("experimental_title")]
        return "\n".join(out)

    def plot(self, **kw):
        from . import graphics
        return graphics.plot_path(self, **kw)


def cf_test_path(blocks, Rmax=None, c=2, q=2, axes="free", extraction="simultaneous", scheme="rows",
                 groups=None, B=999, alpha=0.05, statistic="lobo", matching=True, seed=None, tol=1e-10,
                 maxit=5000, nstart=None, subspace="inertia", weighting="none"):
    """EXPERIMENTAL. Tests of the models with r = 1, ..., Rmax dimensions and choice of R.

    R is chosen as the largest r for which every component of the other blocks
    is linked (min_k p_kh <= alpha / K). This function is experimental: the
    choice of R is not part of the published articles, the models are
    re-estimated for each r (not nested for simultaneous methods), and the
    p-values of the chosen model are optimistic. Report the whole profile.
    A warning is issued at each call and the printed result says so.

    Parameters
    ----------
    Rmax : int, optional
        Default min(K + c, n - 1).
    Other parameters
        As in :func:`cf_test`.

    Returns
    -------
    CfTestPath
    """
    warnings.warn(_tr("experimental_warning"), UserWarning, stacklevel=2)
    names, arrs, _, _ = _as_blocks(blocks)
    K, n = len(arrs), arrs[0].shape[0]
    if Rmax is None:
        Rmax = min(K + c, n - 1)
    Rmax = min(Rmax, n - 1)
    if seed is None:
        seed = int(np.random.default_rng().integers(1, 2 ** 31 - 2))
    tests = {}
    for r in range(1, Rmax + 1):
        tests[r] = cf_test(dict(zip(names, arrs)), R=r, q=q, axes=axes, extraction=extraction, scheme=scheme,
                           groups=groups, B=B, alpha=alpha, statistic=statistic, matching=matching,
                           seed=seed + r - 1, tol=tol, maxit=maxit, nstart=nstart, subspace=subspace,
                           weighting=weighting)
    prof = pd.DataFrame({
        "r": list(tests),
        "retained": [int(t.retained.sum()) for t in tests.values()],
        "all": [int(np.sum(t.link == "all")) for t in tests.values()],
        "several": [int(np.sum(t.link == "several")) for t in tests.values()],
        "one": [int(np.sum(t.link == "one")) for t in tests.values()],
        "all_retained": [bool(t.retained.all()) for t in tests.values()]})
    ok = prof.loc[prof.all_retained, "r"]
    return CfTestPath(tests=tests, profile=prof, R=int(ok.max()) if len(ok) else 0, Rmax=Rmax,
                      experimental=True,
                      settings=dict(c=c, q=q, axes=axes, extraction=extraction, scheme=scheme, B=B, alpha=alpha,
                                    statistic=statistic, seed=seed))
