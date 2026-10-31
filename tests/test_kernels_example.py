"""examples/kernels: generate a Yukawa kernel, build the C++ and Fortran drivers and check that both
converge with the order and print the same errors."""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "kernels"
ROOT = EXAMPLE.parent.parent

pytestmark = pytest.mark.skipif(shutil.which("g++") is None or shutil.which("gfortran") is None,
                                reason="needs g++ and gfortran")


def test_cpp_and_fortran_drivers_agree(tmp_path):
    work = tmp_path / "kernels"
    shutil.copytree(EXAMPLE, work)
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    out = {}
    for target in ("run_cpp", "run_f90"):
        r = subprocess.run(["make", target, f"PYTHON={sys.executable}"], cwd=work, env=env, check=True,
                           capture_output=True, text=True)
        out[target] = [float(x) for x in re.findall(r"^\s+\d+\s+([0-9.]+[eE][-+]?\d+)\s*$", r.stdout, re.M)]
    cpp, f90 = out["run_cpp"], out["run_f90"]
    assert len(cpp) == len(f90) == 4
    assert cpp[-1] < 1e-3 * cpp[0]
    assert f90 == pytest.approx(cpp, rel=1e-6)
