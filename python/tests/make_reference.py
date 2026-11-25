"""Generate the reference values shared by the R and Python tests (inst/shared/reference/reference_values.json).

The R tests (tests/testthat/test-crosscheck.R) and the Python tests (tests/test_canofac.py) must both
retrieve these values: same generator, same algorithms, same starts.
Usage: python3 tests/make_reference.py <output.json>
"""
import sys, json
import numpy as np
sys.path.insert(0, ".")
import canofac as cf
from canofac._rng import Generator


def ref_blocks():
    n = 12

    def mk(k, p):
        i = np.arange(1, n + 1)[:, None] * np.ones((1, p))
        j = np.ones((n, 1)) * np.arange(1, p + 1)[None, :]
        return np.sin(0.37 * i * j + k) + 0.5 * np.cos(1.3 * i + 0.7 * j * k) + 0.3 * k * np.sin(0.9 * i) * np.cos(0.4 * j)
    return cf.cf_preprocess({"A": mk(1, 4), "B": mk(2, 6), "C": mk(3, 9)})


def noise_blocks():
    """Blocks drawn by the common generator (seed 99, column-major), weakly linked: informative p-values."""
    g = Generator(99)
    n = 15
    f = g.normal(n)
    out = {}
    for name, p, a in (("A", 4, 0.8), ("B", 5, 0.5), ("C", 6, 0.0)):
        X = g.normal(n * p).reshape((n, p), order="F")
        X[:, 0] = X[:, 0] + a * f
        out[name] = X
    return cf.cf_preprocess(out)


SETTINGS = [(1, "free", "simultaneous"), (2, "free", "simultaneous"), (4, "free", "simultaneous"),
            (4, "free", "sequential"), (2, "orthonormal", "sequential"), (2, "orthonormal", "simultaneous"),
            (1, "orthonormal", "simultaneous"), (3, "free", "simultaneous")]
TESTS = [dict(name="rows_inertia", kw=dict(q=2)), dict(name="rows_equal", kw=dict(q=2, subspace="equal")),
         dict(name="rows_q4", kw=dict(q=4)), dict(name="variables", kw=dict(q=2, scheme="variables")),
         dict(name="components", kw=dict(q=2, scheme="components")),
         dict(name="groups", kw=dict(q=2, scheme="groups", groups=[1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5])),
         dict(name="within", kw=dict(q=2, scheme="within", groups=[1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5])),
         dict(name="statis", kw=dict(q=2, weighting="statis")), dict(name="mfa", kw=dict(q=2, weighting="mfa")),
         dict(name="gamma", kw=dict(q=2, pvalue="gamma")), dict(name="hybrid", kw=dict(q=2, pvalue="hybrid")),
         dict(name="compare", kw=dict(q=2, compare=True))]


def main(out):
    g = Generator(1)
    ref = {"rng": {"uniform_seed1": g.uniform(5).tolist(),
                   "permutations_seed11_B3_n7": (cf.cf_permutations(3, 7, seed=11) + 1).tolist(),
                   "normal_seed12345": Generator(12345).normal(4).tolist()}}
    bl = ref_blocks()
    ref["fits"] = []
    for q, ax, ex in SETTINGS:
        f = cf.canofac(bl, 2, q, ax, ex)
        ref["fits"].append(dict(q=q, axes=ax, extraction=ex, d=f.d.tolist(), criterion=f.criterion))
    ref["tests"] = []
    nb = noise_blocks()
    ref["noise_blocks_first_row"] = {k: v[0].tolist() for k, v in nb.items()}
    for t in TESTS:
        r = cf.cf_test(nb, R=2, B=99, seed=1, **t["kw"])
        kw = {k: (list(v) if isinstance(v, (list, tuple)) else v) for k, v in t["kw"].items()}
        ent = dict(name=t["name"], kw=kw, B=99, seed=1, p=r.p.tolist(), p_subspace=r.p_subspace.tolist(),
                   p_subspace_equal=r.p_subspace_equal.tolist(), p_subspace_inertia=r.p_subspace_inertia.tolist(),
                   link=list(r.link), link_subspace=r.link_subspace)
        if r.p_compare is not None:
            ent["p_compare"] = r.p_compare.values.tolist()
        ref["tests"].append(ent)
    sim = cf.cf_simulate("SD", n=12, p=5, K=3, theta=1.8, seed=7)
    ref["simulate_SD"] = {k: v[:2].tolist() for k, v in sim.items()}
    W = bl["A"] @ bl["A"].T
    u = np.linalg.svd(np.hstack([bl["B"], bl["C"]]), full_matrices=False)[0][:, :2]
    m, v = cf.perm_moments(W, u @ u.T)
    ref["moments"] = dict(mean=m, var=v)
    json.dump(ref, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
