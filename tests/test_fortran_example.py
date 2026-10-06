"""End-to-end test of the Fortran example in examples/fortran.

Generates a low-order operator module, builds the example driver with
gfortran (OpenMP on), runs it and checks the result is physically sensible:
errors fall with expansion order, variants that must agree do agree, and the
answer does not depend on the thread count.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import fmmgen

gfortran = shutil.which("gfortran")
DRIVER = Path(__file__).resolve().parent.parent / "examples" / "fortran" / "fmm_example.f90"

# CI sets this so a missing compiler fails the job instead of skipping it.
if gfortran is None and os.environ.get("FMMGEN_REQUIRE_FORTRAN"):
    pytest.fail("FMMGEN_REQUIRE_FORTRAN is set but gfortran was not found", pytrace=False)

pytestmark = pytest.mark.skipif(gfortran is None, reason="gfortran not available")


@pytest.fixture(scope="module")
def exe(tmp_path_factory):
    d = tmp_path_factory.mktemp("fortran_example")
    # Orders 1..3: enough to see convergence, quick to generate.
    fmmgen.generate_code(4, "operators", CSE=True, minpow=11, harmonic_derivs=True,
                         compress=True, planar=True, language="fortran", src_dir=str(d))
    subprocess.run([gfortran, "-O1", "-fopenmp", "-Wno-unused-dummy-argument",
                    "operators.f90", str(DRIVER), "-o", "main"],
                   cwd=d, check=True)
    return d


def run(exe_dir, *args, threads=1):
    out = subprocess.run(["./main", "--nparticles", "1000", *args], cwd=exe_dir, check=True,
                         capture_output=True, text=True,
                         env=dict(os.environ, OMP_NUM_THREADS=str(threads))).stdout
    errs, order = {}, None
    for line in out.splitlines():
        m = re.match(r"Order (\d+)", line)
        if m:
            order = int(m.group(1))
        if line.startswith("L2 errs"):
            errs[order] = [float(v) for v in line.split("=")[1].split(",") if v.strip()]
    return errs


@pytest.mark.parametrize("args", [
    ["--type", "0"],
    ["--type", "1"],
    ["--type", "0", "--compress", "1"],
    ["--type", "1", "--dim", "2"],
    ["--type", "0", "--dim", "2"],
], ids=["fmm", "barnes-hut", "fmm-compressed", "barnes-hut-planar", "fmm-planar"])
def test_error_decreases_with_order(exe, args):
    errs = run(exe, *args)
    assert sorted(errs) == [1, 2, 3]
    pot = [errs[p][0] for p in (1, 2, 3)]
    assert pot[0] > pot[1] > pot[2]
    assert pot[2] < 5e-2


def test_fmm_reaches_expected_accuracy(exe):
    # Exact-field comparison, so a wrong operator or tree bug shows up here.
    assert run(exe, "--type", "0")[3][0] < 1e-3


def test_compressed_matches_full(exe):
    full, comp = run(exe, "--type", "0"), run(exe, "--type", "0", "--compress", "1")
    for p in (1, 2, 3):
        assert comp[p] == pytest.approx(full[p], rel=1e-6)


@pytest.mark.parametrize("args", [["--type", "0"], ["--type", "1"]], ids=["fmm", "barnes-hut"])
def test_thread_count_does_not_change_result(exe, args):
    assert run(exe, *args, threads=1) == run(exe, *args, threads=3)


def test_input_file_round_trip(exe):
    """Reading back the generated particle file reproduces the generated run."""
    generated = run(exe, "--type", "0")
    particles = exe / "particles_n_1000.txt"
    assert particles.exists()
    assert run(exe, "--type", "0", "--input", str(particles)) == generated
