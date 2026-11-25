"""Tests of the Python version of canofac; the reference values are shared with the R tests."""
import os, sys, json, warnings
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import canofac as cf
from canofac._rng import Generator
from make_reference import ref_blocks, noise_blocks, SETTINGS, TESTS

REF = json.load(open(os.path.join(os.path.dirname(cf.__file__), "shared", "reference", "reference_values.json")))
SLOW = {(1, "free", "simultaneous"), (2, "orthonormal", "simultaneous")}   # flat criteria: slow convergence of d


def test_generator_reference():
    assert np.allclose(Generator(1).uniform(5), REF["rng"]["uniform_seed1"], rtol=0, atol=1e-15)
    assert (cf.cf_permutations(3, 7, seed=11) + 1).tolist() == REF["rng"]["permutations_seed11_B3_n7"]
    assert np.allclose(Generator(12345).normal(4), REF["rng"]["normal_seed12345"], atol=1e-12)


def test_generator_batch_equals_sequential():
    g1, g2 = Generator(5), Generator(5)
    P = g1.permutations(7, 11)
    Q = np.array([g2.permutation(11) for _ in range(7)])
    assert (P == Q).all() and g1.state == g2.state
    assert all(sorted(p) == list(range(11)) for p in P)


@pytest.mark.parametrize("i", range(len(SETTINGS)))
def test_fits_reference(i):
    q, ax, ex = SETTINGS[i]
    e = REF["fits"][i]
    f = cf.canofac(ref_blocks(), 2, q, ax, ex)
    assert abs(f.criterion - e["criterion"]) <= 1e-7
    tol = 2e-3 if (q, ax, ex) in SLOW else 1e-6
    assert np.allclose(f.d, e["d"], atol=tol)


@pytest.mark.parametrize("i", range(len(TESTS)))
def test_tests_reference(i):
    t, e = TESTS[i], REF["tests"][i]
    r = cf.cf_test(noise_blocks(), R=2, B=99, seed=1, **t["kw"])
    assert np.allclose(r.p, e["p"], atol=1e-10)
    assert np.allclose(r.p_subspace, e["p_subspace"], atol=1e-10)
    assert list(r.link) == e["link"] and r.link_subspace == e["link_subspace"]


def test_moments_against_permutations():
    bl = ref_blocks()
    W = bl["A"] @ bl["A"].T
    u = np.linalg.svd(np.hstack([bl["B"], bl["C"]]), full_matrices=False)[0][:, :2]
    M = u @ u.T
    m, v = cf.perm_moments(W, M)
    P = cf.cf_permutations(20000, 12, seed=3)
    T = np.array([np.sum(W[np.ix_(p, p)] * M) for p in P])
    assert abs(T.mean() - m) < 4 * T.std() / np.sqrt(len(T))
    assert abs(T.var() / v - 1) < 0.05


def test_level_under_h0():
    rej = []
    for s in range(40):
        g = Generator(1000 + s)
        bl = cf.cf_preprocess([g.normal(20 * p).reshape((20, p), order="F") for p in (4, 6, 5)])
        rej.append(cf.cf_test(bl, R=1, B=39, seed=s).p[0, 0] <= 0.05)
    assert np.mean(rej) < 0.25


def test_language_and_values():
    bl = ref_blocks()
    t = cf.cf_test(bl, R=2, B=19, seed=2)
    old = cf.cf_lang("fr")
    txt = repr(t)
    cf.cf_lang(old)
    assert "Blocs" in txt and set(t.link) <= {"all", "several", "one", "none"}


def test_path_is_flagged_experimental():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p = cf.cf_test_path(ref_blocks(), Rmax=2, B=9, seed=1)
    assert any("experimental" in str(x.message) for x in w)
    assert repr(p).startswith("EXPERIMENTAL") and p.experimental


def test_gamma_hybrid_and_errors():
    bl = noise_blocks()
    g = cf.cf_test(bl, R=2, seed=1, pvalue="gamma")
    assert (g.pmethod == "gamma").all() and np.all((g.p > 0) & (g.p <= 1))
    with pytest.raises(ValueError):
        cf.cf_test(bl, R=2, pvalue="gamma", scheme="variables")
    with pytest.raises(ValueError):
        cf.cf_test(bl, R=2, scheme="groups")


def test_trajectories_and_tools(tmp_path):
    bl = ref_blocks()
    tr = cf.cf_trajectories(bl, R=2, qgrid=[1, 2, 4])
    assert tr.delta.shape == (3, 2, 3) and np.all(tr.delta <= 1 + 1e-10)
    ms = cf.cf_multistart(bl, R=2, q=4, nstart=5)
    assert 0 < ms.share_best <= 100
    a = cf.cf_angles(cf.canofac(bl, 2, 2), cf.canofac(bl, 2, 2, extraction="sequential"))
    assert a["subspace"] < 1e-4
    f = str(tmp_path / "b.csv")
    cf.cf_write_blocks(bl, f)
    back, _ = cf.cf_read_blocks(f)
    assert np.allclose(np.asarray(back["B"]), bl["B"])


def test_plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bl = ref_blocks()
    f = cf.canofac(bl, 2, 2)
    t = cf.cf_test(bl, R=2, B=19, seed=1)
    for fn in (lambda: f.plot(), lambda: f.plot("loadings"), lambda: cf.cf_plot_saliences(f),
               lambda: cf.cf_trajectories(bl, qgrid=[1, 2, 4]).plot(), lambda: t.plot(), lambda: t.plot("null"),
               lambda: cf.cf_rose(f), lambda: cf.cf_plot_compare([f, cf.canofac(bl, 2, 4)])):
        fn()
        plt.close("all")


def test_simulate_reference():
    sim = cf.cf_simulate("SD", n=12, p=5, K=3, theta=1.8, seed=7)
    for k, v in REF["simulate_SD"].items():
        assert np.allclose(sim[k][:2], v, atol=1e-10)
