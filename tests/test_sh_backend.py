"""The spherical-harmonic backend through `generate`: what it accepts, what it writes, and that it links and
runs next to the Cartesian backend in one library."""
import ctypes
import re
import subprocess

import numpy as np
import pytest

import fmmgen
from fmmgen import Kernel, Problem, generate
from fmmgen.backend import _is_laplace
from fmmgen.kernel import power_law, yukawa
from fmmgen.sh import generate_sh_code
from kernel_util import HAVE_GCC, cluster, direct


def test_only_the_laplace_charge_problem_is_supported(tmp_path):
    import sympy as sp

    x, y, z = sp.symbols("x y z")
    assert _is_laplace(power_law(1)) and not _is_laplace(yukawa(1.0))
    assert not _is_laplace(Kernel(2 / sp.sqrt(x**2 + y**2 + z**2), (x, y, z)))  # 1/r only, not a multiple of it
    kw = dict(backend="spherical", orders=[2], src_dir=str(tmp_path), include_dir=str(tmp_path))
    for what, args, msg in (
        ("yukawa", dict(kernel=yukawa(1.0)), "Laplace"),
        ("dipoles", dict(kernel=power_law(1), source_order=1), "charges only"),
        ("two channels", dict(kernel=power_law(1), nchannels=2), "charges only"),
        ("potential only", dict(kernel=power_law(1), potential=False), "potential and the field"),
    ):
        with pytest.raises(ValueError, match=msg):
            generate(args.pop("kernel"), 3, "s", **args, **kw)
    with pytest.raises(ValueError, match="planar or compressed"):
        generate(power_law(1), 3, "s", planar=True, **kw)
    with pytest.raises(ValueError, match="language"):
        generate(power_law(1), 3, "s", language="fortran", **{k: v for k, v in kw.items()})
    assert not list(tmp_path.glob("s.*"))  # nothing is written for a refused request
    # recommend lists it for Laplace charges and says why not otherwise
    rows = {n: why for n, _, why in fmmgen.recommend(Problem(power_law(1)), [4])[4]}
    assert rows["spherical"] == "" and rows["cartesian"] == ""
    assert "Laplace" in {n: why for n, _, why in fmmgen.recommend(Problem(yukawa(1.0)), [4])[4]}["spherical"]


def test_manifest_lists_what_the_header_declares(tmp_path):
    man = generate(power_law(1), 4, "s", backend="spherical", orders=[2, 4], m2l="rotation",
                   src_dir=str(tmp_path), include_dir=str(tmp_path))
    header = (tmp_path / "s.h").read_text()
    declared = set(re.findall(r"void (\w+)\(", header))
    assert set(man.functions) <= declared and {"S2M_4", "M2Lrot_setup_2", "P2P_batch", "M2L"} <= set(man.functions)
    assert "S2M_3" not in man.functions  # only the requested degrees
    assert man.layouts[4] == {"multipole": 25, "local": 25, "source": 1, "output": 4}
    assert man.options["m2l"] == "rotation" and man.options["rotation_table_size"][4] > 0
    assert man.files == ["s.h", "s.c"] and man.dependencies == ["libm"]
    assert "M2Lrot_4(x, y, z, M, L);" in (tmp_path / "s.c").read_text()  # M2L_4 calls the rotation M2L


def test_generator_orders_prefix_and_validation(tmp_path):
    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path))
    generate_sh_code(7, "p", orders=[2, 5], prefix="sh_", **kw)
    header = (tmp_path / "p.h").read_text()
    assert "void sh_M2L_2(" in header and "void sh_M2L_5(" in header and "void sh_M2L_3(" not in header
    assert "#define FMMGEN_SH_MAXORDER 7" in header and "FMMGEN_SH_MULTIPOLESIZE[]" in header
    assert "FMMGEN_MAXORDER" not in header and "void M2L(" not in header
    wrapper = re.search(r"void sh_M2L\(.*?\n\}\n", (tmp_path / "p.c").read_text(), re.S).group(0)
    assert "case 2:" in wrapper and "case 5:" in wrapper and "case 3:" not in wrapper
    for bad, msg in ((dict(prefix="a-b"), "identifier"), (dict(orders=[]), "orders"), (dict(orders=[0, 2]), "orders"),
                     (dict(orders=[2, 7]), "exceed")):
        with pytest.raises(ValueError, match=msg):
            generate_sh_code(7, "q", **{**kw, **bad})


class _Shared:
    def __init__(self, path):
        self.lib = ctypes.CDLL(str(path))

    def call(self, fn, xyz, *arrays):
        f = getattr(self.lib, fn)
        dp = ctypes.POINTER(ctypes.c_double)
        f.argtypes = [ctypes.c_double] * 3 + [dp] * len(arrays)
        f.restype = None
        f(*xyz, *[a.ctypes.data_as(dp) for a in arrays])


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
def test_spherical_and_cartesian_run_in_one_library(tmp_path):
    """Cartesian at low order and spherical harmonics at high order, one shared library, both accurate."""
    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path))
    man = generate(power_law(1), 6, "mix", backend={"default": "cartesian", 6: "spherical"}, orders=[3, 6], **kw)
    assert man.by_order == {3: "cartesian", 6: "spherical"}
    cfiles = [str(tmp_path / f"mix_{b}.c") for b in ("cartesian", "spherical")]
    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", "-I", str(tmp_path), *cfiles, "-o", str(tmp_path / "mix.so"),
                    "-lm"], check=True, capture_output=True)
    lib = _Shared(tmp_path / "mix.so")
    rng = np.random.default_rng(1)
    cs, ct = np.zeros(3), np.array([3.0, 1.0, -2.0])
    src, tgt, q = cluster(rng, cs, 20, 0.5, 3), cluster(rng, ct, 6, 0.5, 3), rng.standard_normal(20)
    phi_ref, _ = direct(power_law(1), tgt, src, q)
    errs = {}
    for prefix, p, size in (("cartesian_", 3, 20), ("spherical_", 6, 49)):
        M, L = np.zeros(size), np.zeros(size)
        for y_, qq in zip(src, q):
            lib.call(f"{prefix}S2M_{p}", cs - y_, np.array([qq]), M)
        lib.call(f"{prefix}M2L_{p}", ct - cs, M, L)
        phi = np.zeros(6)
        for i, t in enumerate(tgt):
            F = np.zeros(4)
            lib.call(f"{prefix}L2P_{p}", t - ct, L, F)
            phi[i] = F[0]
        errs[prefix] = np.linalg.norm(phi - phi_ref) / np.linalg.norm(phi_ref)
    assert errs["cartesian_"] < 5e-2 and errs["spherical_"] < 1e-3 and errs["spherical_"] < errs["cartesian_"]
