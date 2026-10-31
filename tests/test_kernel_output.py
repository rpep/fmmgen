"""What the generated text looks like, for every library kernel: valid Fortran 90,
no integer literal too big for C, and the same numbers from Fortran and C++."""
import os
import re
import shutil
import subprocess

import pytest

from fmmgen.kernel import (gaussian_screened_coulomb, generate_kernel_code, nterms, power_law, shell_indices,
                           source_indices, yukawa)
from test_kernel_library import LIB

gfortran = shutil.which("gfortran")
if gfortran is None and os.environ.get("FMMGEN_REQUIRE_FORTRAN"):
    pytest.fail("FMMGEN_REQUIRE_FORTRAN is set but gfortran was not found", pytrace=False)


def test_no_integer_literal_overflows_c(tmp_path):
    # Coefficients that sympy folds into exact integers can exceed `long`; C then
    # rejects or silently changes the literal.
    bad = []
    for i, (name, make) in enumerate(LIB.items()):
        generate_kernel_code(make(), 5, f"k{i}", language="c", src_dir=str(tmp_path), include_dir=str(tmp_path))
        text = (tmp_path / f"k{i}.c").read_text()
        if re.search(r"(?<![\w.])\d{19,}(?![\w.])", text):
            bad.append(name)
    assert not bad, bad


@pytest.mark.skipif(gfortran is None, reason="gfortran not available")
def test_fortran_is_valid_for_every_library_kernel(tmp_path):
    bad = {}
    for i, (name, make) in enumerate(LIB.items()):
        mod = f"k{i}"  # library names such as r^-1 are not valid module names
        generate_kernel_code(make(), 3, mod, language="fortran", src_dir=str(tmp_path))
        r = subprocess.run([gfortran, "-std=f2008", "-fsyntax-only", f"{mod}.f90"], cwd=tmp_path,
                           capture_output=True, text=True)
        if r.returncode:
            bad[name] = r.stderr.strip().splitlines()[:3]
    assert not bad, bad


CHAIN_CASES = [
    ("erf", lambda: gaussian_screened_coulomb(0.7), 0, 1),  # Piecewise radial recurrence -> merge()
    ("dipole2", lambda: yukawa(1.0), 1, 2),  # order-1 sources, two channels: the array layouts differ
    ("quad", lambda: power_law(3), 2, 1),
]


def _chain_driver(lang, p, s, nch):
    """Same calls in Fortran and C++: S2M, M2M, M2L, L2L, L2P and P2P, printing every array."""
    d = 3
    shell = len(shell_indices(d, s)) * nch
    nm, nl, nf = len(source_indices(d, p, s)) * nch, nterms(d, p - s) * nch, 4 * nch  # F: phi and 3 field terms/channel
    pts = [(0.3, -0.2, 0.1), (-0.5, 0.4, 0.2), (2.0, -1.5, 3.0), (0.4, 0.1, -0.3), (0.2, -0.1, 0.3), (1.1, 0.7, -0.9)]
    if lang == "f90":
        c = [", ".join(f"{v}d0" for v in pt) for pt in pts]
        return f"""program d
use k_f90
real(wp) :: S({shell}), M({nm}), Ms({nm}), L({nl}), Ls({nl}), F({nf}), F2({nf})
integer :: i
do i = 1, {shell}
  S(i) = sin(0.9d0*i) + 0.3d0
end do
M = 0; Ms = 0; L = 0; Ls = 0; F = 0; F2 = 0
call S2M({c[0]}, S, M, {p})
call M2M({c[1]}, M, Ms, {p})
call M2L({c[2]}, Ms, L, {p})
call L2L({c[3]}, L, Ls, {p})
call L2P({c[4]}, Ls, F, {p})
call P2P({c[5]}, S, F2)
write(*,'(es24.16)') Ms, L, Ls, F, F2
end program
"""
    c = [", ".join(str(v) for v in pt) for pt in pts]
    return f"""#include <cstdio>
#include <cmath>
#include "k_cpp.h"
int main() {{
  double S[{shell}], M[{nm}] = {{0}}, Ms[{nm}] = {{0}}, L[{nl}] = {{0}}, Ls[{nl}] = {{0}};
  double F[{nf}] = {{0}}, F2[{nf}] = {{0}};
  for (int i = 1; i <= {shell}; i++) S[i-1] = std::sin(0.9*i) + 0.3;
  S2M({c[0]}, S, M, {p});  M2M({c[1]}, M, Ms, {p});  M2L({c[2]}, Ms, L, {p});
  L2L({c[3]}, L, Ls, {p});  L2P({c[4]}, Ls, F, {p});  P2P({c[5]}, S, F2);
  for (double v : Ms) std::printf("%24.16e\\n", v);
  for (double v : L) std::printf("%24.16e\\n", v);
  for (double v : Ls) std::printf("%24.16e\\n", v);
  for (double v : F) std::printf("%24.16e\\n", v);
  for (double v : F2) std::printf("%24.16e\\n", v);
}}
"""


@pytest.mark.skipif(gfortran is None or shutil.which("g++") is None, reason="needs gfortran and g++")
@pytest.mark.parametrize("name,make,s,nch", CHAIN_CASES, ids=[c[0] for c in CHAIN_CASES])
def test_fortran_and_cpp_agree(tmp_path, name, make, s, nch):
    p = 4
    for lang, ext in (("fortran", "f90"), ("c++", "cpp")):
        generate_kernel_code(make(), p, f"k_{ext}", language=lang, orders=[p], source_order=s, nchannels=nch,
                             src_dir=str(tmp_path), include_dir=str(tmp_path))
        (tmp_path / f"d.{ext}").write_text(_chain_driver(ext, p, s, nch))
    subprocess.run([gfortran, "-O1", "k_f90.f90", "d.f90", "-o", "df"], cwd=tmp_path, check=True)
    subprocess.run(["g++", "-O1", "k_cpp.cpp", "d.cpp", "-o", "dc"], cwd=tmp_path, check=True)
    f90, cpp = (list(map(float, subprocess.run([f"./{e}"], cwd=tmp_path, capture_output=True, text=True,
                                               check=True).stdout.split())) for e in ("df", "dc"))
    assert len(f90) == len(cpp) and len(cpp) > 8
    scale = max(abs(v) for v in cpp)
    assert max(abs(a - b) for a, b in zip(f90, cpp)) < 1e-12 * scale


def test_fortran_printer_maps_conditionals_and_constants():
    import sympy as sp

    from fmmgen.printers import FortranCodePrinter

    pr = FortranCodePrinter()
    t = sp.Symbol("t")
    # three branches nest into merge(); a comparison on its own is a real 0/1 (CSE temporaries are real)
    pw = sp.Piecewise((1, t > 2), (2, t > 1), (3, True))
    assert pr.doprint(pw) == "merge(1.0d0, merge(2.0d0, 3.0d0, t > 1.0d0), t > 2.0d0)"
    assert pr.doprint(t > 2) == "merge(1.0d0, 0.0d0, t > 2.0d0)"
    assert pr.doprint(sp.Piecewise((5, True))) == "5.0d0"
    # a constant is a literal, not a `parameter` declaration, and extra-precision floats are shortened
    assert pr.doprint(sp.pi) == "3.1415926535897932d0"
    assert pr.doprint(sp.Float("0.1234567890123456789012345678901234567890", 60)) == "0.12345678901234568d0"


def test_header_macros_agree_with_generate_code(tmp_path):
    # a driver loops `order < FMMGEN_MAXORDER` over either generator's output: the conventions must match
    import fmmgen

    for sub in ("old", "new"):
        (tmp_path / sub).mkdir()
    fmmgen.generate_code(4, "ops", src_dir=str(tmp_path / "old"), include_dir=str(tmp_path / "old"), CSE=True,
                         minpow=11, cache_dir=None)
    generate_kernel_code(power_law(1), 3, "ops", src_dir=str(tmp_path / "new"), include_dir=str(tmp_path / "new"))

    def macros(d):
        text = (tmp_path / d / "ops.h").read_text()
        return {k: re.search(rf"#define FMMGEN_{k} (\d+)", text).group(1)
                for k in ("MINORDER", "MAXORDER", "SOURCEORDER", "SOURCESIZE", "OUTPUTSIZE")}

    assert macros("old") == macros("new") == {"MINORDER": "1", "MAXORDER": "4", "SOURCEORDER": "0",
                                              "SOURCESIZE": "1", "OUTPUTSIZE": "4"}
