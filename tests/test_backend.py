"""The backend layer: registry, auto selection, layout and manifest."""
import json

import numpy as np
import pytest

import fmmgen
from fmmgen import Backend, Kernel, Problem, backends, generate, register_backend
from fmmgen.backend import CartesianBackend
from fmmgen.kernel import nterms, power_law, source_indices
from fmmgen.vector import stokeslet
from kernel_util import HAVE_GCC, Lib


def test_registry_and_explicit_backend():
    assert "cartesian" in backends()
    pr = Problem(power_law(1))
    assert fmmgen.backend._resolve(pr, "cartesian", [2, 3])[3].name == "cartesian"
    with pytest.raises(ValueError):
        fmmgen.backend.get_backend("nope")
    with pytest.raises(ValueError, match="not being generated"):
        fmmgen.backend._resolve(pr, {"default": "cartesian", 9: "cartesian"}, [2, 3])


class _Cheap(Backend):
    name = "cheap"

    def supports(self, problem):
        return problem.source_order == 0, "charges only"

    def cost(self, problem, order):
        return 1.0 if order < 6 else 9.0  # the best backend depends on the order


class _Dear(CartesianBackend):
    name = "dear"

    def cost(self, problem, order):
        return 5.0


@pytest.fixture
def two_more_backends():
    register_backend(_Cheap())
    register_backend(_Dear())
    yield
    fmmgen.backend._REGISTRY.pop("cheap")
    fmmgen.backend._REGISTRY.pop("dear")


def test_recommend_is_advice_per_order(two_more_backends):
    rec = fmmgen.recommend(Problem(power_law(1)), [4, 8])
    assert rec[4][0][0] == "cheap" and rec[8][0][0] == "dear"  # ranking changes with order
    assert [n for n, _, _ in rec[4]][-1] == "cartesian"  # unknown cost goes last
    dip = fmmgen.recommend(Problem(power_law(1), source_order=1), [4])
    assert dip[4][-1][0] == "cheap" and dip[4][-1][2] == "charges only"


def test_per_order_mapping_and_unsupported_mix(two_more_backends, tmp_path):
    pr = Problem(power_law(1), source_order=1)
    with pytest.raises(ValueError, match="charges only"):
        fmmgen.backend._resolve(pr, {"default": "dear", 3: "cheap"}, [2, 3, 4])
    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path))
    generate(power_law(1), 3, "same", backend={"default": "dear"}, orders=[1, 2, 3], **kw)  # one backend: fine

    class NoPrefix(Backend):
        name = "noprefix"

        def supports(self, problem):
            return True, ""

    register_backend(NoPrefix())
    try:
        with pytest.raises(NotImplementedError, match="prefixed"):
            generate(power_law(1), 4, "mix", backend={"default": "cartesian", 4: "noprefix"}, orders=[2, 4], **kw)
    finally:
        fmmgen.backend._REGISTRY.pop("noprefix")


def test_mixed_build_manifest(two_more_backends, tmp_path):
    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path))
    man = generate(power_law(1), 4, "mix", backend={"default": "cartesian", 4: "dear"}, orders=[2, 3, 4], **kw)
    assert man.backend == "mixed" and man.by_order == {2: "cartesian", 3: "cartesian", 4: "dear"}
    data = json.loads((tmp_path / "mix.json").read_text())
    parts = {p["prefix"]: p for p in data["parts"]}
    assert set(parts) == {"cartesian_", "dear_"} and parts["dear_"]["orders"] == [4]
    assert "cartesian_M2L_3" in parts["cartesian_"]["functions"] and "dear_M2L_4" in parts["dear_"]["functions"]
    assert (tmp_path / "mix_cartesian.c").exists() and (tmp_path / "mix_dear.c").exists()


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
def test_backends_coexist_in_one_library(two_more_backends, tmp_path):
    """Two backends' operators link into one shared library and both run."""
    import ctypes
    import subprocess

    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path))
    generate(power_law(1), 4, "mix", backend={"default": "cartesian", 4: "dear"}, orders=[3, 4], **kw)
    so = tmp_path / "mix.so"
    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", "-I", str(tmp_path), str(tmp_path / "mix_cartesian.c"),
                    str(tmp_path / "mix_dear.c"), "-o", str(so), "-lm"], check=True, capture_output=True)
    lib = ctypes.CDLL(str(so))
    dp = ctypes.POINTER(ctypes.c_double)
    out = {}
    for fn, p in (("cartesian_S2M_3", 3), ("dear_S2M_4", 4)):
        f = getattr(lib, fn)
        f.argtypes = [ctypes.c_double] * 3 + [dp, dp]
        f.restype = None
        M = np.zeros(nterms(3, p))
        q = np.array([1.0])
        f(0.1, 0.2, 0.3, q.ctypes.data_as(dp), M.ctypes.data_as(dp))
        out[fn] = M
    assert out["cartesian_S2M_3"][0] == 1.0 and out["dear_S2M_4"][0] == 1.0
    assert len(out["cartesian_S2M_3"]) != len(out["dear_S2M_4"])


def test_prefix_applies_to_macros_and_functions(tmp_path):
    generate(power_law(1), 2, "pre", backend="cartesian", orders=[1, 2], prefix="lap_", src_dir=str(tmp_path),
             include_dir=str(tmp_path))
    h = (tmp_path / "pre.h").read_text()
    assert "lap_M2L_2" in h and "lap_P2P_batch" in h and "FMMGEN_LAP_MAXORDER 2" in h
    assert "void M2L_2" not in h
    with pytest.raises(ValueError):
        generate(power_law(1), 2, "bad", prefix="x-y", orders=[1, 2], src_dir=str(tmp_path), include_dir=str(tmp_path))


@pytest.mark.parametrize("s,nch", [(0, 1), (1, 2), (2, 3)])
def test_cartesian_layout_matches_formulas(s, nch):
    k = power_law(1)
    lay = CartesianBackend().layout(Problem(k, s, nch), 5)
    assert lay.multipole == len(source_indices(3, 5, s)) * nch
    assert lay.local == nterms(3, 5 - s) * nch
    assert lay.output == 4 * nch


def test_generate_writes_manifest(tmp_path):
    man = generate(power_law(1), 4, "ops", source_order=1, nchannels=2, src_dir=str(tmp_path),
                   include_dir=str(tmp_path))
    data = json.loads((tmp_path / "ops.json").read_text())
    assert data["backend"] == "cartesian" and data["orders"] == [2, 3, 4] and data["nchannels"] == 2
    assert (tmp_path / "ops.c").exists() and "M2L_3" in data["functions"]
    assert data["layouts"]["3"]["local"] == nterms(3, 2) * 2 and man.tiers == ["operator"]


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
def test_layout_matches_compiled_operators(tmp_path):
    """The manifest sizes are the ones the compiled operators actually read and write."""
    k = power_law(1)
    p, s, nch = 4, 1, 2
    lay = CartesianBackend().layout(Problem(k, s, nch), p)
    lib = Lib(k, tmp_path, [p], source_order=s, nchannels=nch)
    M = np.zeros(lay.multipole)
    S = np.ones(lay.source)
    lib.call(f"S2M_{p}", np.array([0.1, 0.2, 0.3]), S, M)
    L = np.zeros(lay.local)
    lib.call(f"M2L_{p}", np.array([3.0, 0.5, -0.5]), M, L)
    F = np.zeros(lay.output)
    lib.call(f"L2P_{p}", np.array([0.1, 0.0, 0.0]), L, F)
    assert np.all(np.isfinite(F)) and np.any(F != 0)


def test_generate_vector_kernel(tmp_path):
    man = generate(stokeslet(), 3, "st", src_dir=str(tmp_path), include_dir=str(tmp_path))
    assert len(man.groups) == 1 and man.groups[0]["nchannels"] == 3 and man.groups[0]["source_order"] == 2
    assert (tmp_path / "st_g0.c").exists() and (tmp_path / "st.json").exists()


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
def test_vector_groups_link_into_one_program(tmp_path):
    """All groups of a split vector kernel (6 groups) must not define the same symbols."""
    import subprocess

    import sympy as sp

    x, y, z = X = sp.symbols("x y z")
    r = sp.sqrt(x**2 + y**2 + z**2)
    Xa = np.array(X, dtype=object)
    G = np.eye(3, dtype=int) / r + np.einsum("i,j->ij", Xa, Xa) / r**3
    vk = Kernel(G, method="symbolic")
    assert len(vk.groups) == 6
    man = generate(vk, 2, "st", orders=[2], src_dir=str(tmp_path), include_dir=str(tmp_path))
    names = [fn for g in man.groups for fn in g["functions"]]
    assert len(names) == len(set(names)) and "st_g3_M2L_2" in names
    cfiles = sorted(str(p) for p in tmp_path.glob("st_g*.c"))
    subprocess.run(["gcc", "-O0", "-shared", "-fPIC", "-I", str(tmp_path), *cfiles, "-o", str(tmp_path / "st.so"),
                    "-lm"], check=True, capture_output=True)


def test_generate_vector_kernel_fortran(tmp_path):
    # a Fortran group is its own module, so no function prefix is needed (or allowed)
    man = generate(stokeslet(), 2, "st", orders=[2], language="fortran", src_dir=str(tmp_path))
    assert man.groups[0]["files"] == ["st_g0.f90"] and (tmp_path / "st_g0.f90").exists()
