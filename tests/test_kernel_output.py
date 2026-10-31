"""What the generated text looks like, for every library kernel: valid Fortran 90,
no integer literal too big for C, and the same numbers from Fortran and C++."""
import os
import re
import shutil
import subprocess

import pytest

from fmmgen.kernel import gaussian_screened_coulomb, generate_kernel_code
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


@pytest.mark.skipif(gfortran is None or shutil.which("g++") is None, reason="needs gfortran and g++")
def test_fortran_and_cpp_agree_for_erf(tmp_path):
    # erf has a Piecewise radial recurrence: the Fortran printer maps it to merge().
    k, p, n = gaussian_screened_coulomb(0.7), 6, 84
    for lang, ext in (("fortran", "f90"), ("c++", "cpp")):
        generate_kernel_code(k, p, f"gsc_{ext}", language=lang, orders=[p], src_dir=str(tmp_path),
                             include_dir=str(tmp_path))
    (tmp_path / "d.f90").write_text(f"""program d
use gsc_f90
real(wp) :: M({n}), L({n})
integer :: i
do i = 1, {n}
  M(i) = sin(0.37d0*i)
end do
L = 0
call M2L(2.0d0, -1.5d0, 3.0d0, M, L, {p})
do i = 1, {n}
  write(*,'(es24.16)') L(i)
end do
end program
""")
    (tmp_path / "d.cpp").write_text(f"""#include <cstdio>
#include <cmath>
#include "gsc_cpp.h"
int main() {{ double M[{n}], L[{n}] = {{0}}; for (int i = 1; i <= {n}; i++) M[i-1] = std::sin(0.37*i);
  M2L(2.0, -1.5, 3.0, M, L, {p}); for (int i = 0; i < {n}; i++) std::printf("%24.16e\\n", L[i]); }}
""")
    subprocess.run([gfortran, "-O1", "gsc_f90.f90", "d.f90", "-o", "df"], cwd=tmp_path, check=True)
    subprocess.run(["g++", "-O1", "gsc_cpp.cpp", "d.cpp", "-o", "dc"], cwd=tmp_path, check=True)
    a, b = (list(map(float, subprocess.run([f"./{e}"], cwd=tmp_path, capture_output=True, text=True,
                                           check=True).stdout.split())) for e in ("df", "dc"))
    assert len(a) == len(b) == n
    assert max(abs(x - y) / (abs(y) + 1e-300) for x, y in zip(a, b)) < 1e-12
