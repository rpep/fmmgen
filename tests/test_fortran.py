import shutil
import subprocess

import pytest

import fmmgen

gfortran = shutil.which("gfortran")


def test_fortran_rejects_unsupported_options(tmp_path):
    with pytest.raises(NotImplementedError):
        fmmgen.generate_code(3, "ops", language="fortran", atomic=True, src_dir=str(tmp_path))


@pytest.mark.skipif(gfortran is None, reason="gfortran not available")
def test_fortran_compiles_as_strict_f2003(tmp_path):
    fmmgen.generate_code(4, "ops", CSE=True, minpow=11, harmonic_derivs=True,
                         language="fortran", src_dir=str(tmp_path))
    src = (tmp_path / "ops.f90").read_text()
    assert "implicit none" in src
    subprocess.run([gfortran, "-std=f2003", "-Wall", "-Werror", "-c", "ops.f90"],
                   cwd=tmp_path, check=True)


@pytest.mark.skipif(gfortran is None, reason="gfortran not available")
def test_fortran_compressed_and_planar_compile(tmp_path):
    fmmgen.generate_code(4, "ops", CSE=True, minpow=11, harmonic_derivs=True,
                         compress=True, planar=True,
                         language="fortran", src_dir=str(tmp_path))
    src = (tmp_path / "ops.f90").read_text()
    for name in ("FMMGEN_MULTIPOLESIZE", "FMMGEN_PLANAR_LOCALSIZE", "M2Lc", "M2Lxy", "P2P_batchxy"):
        assert name in src
    # Planar kernels ignore z by design, so unused dummies are expected.
    subprocess.run([gfortran, "-std=f2003", "-Wall", "-Wno-unused-dummy-argument", "-Werror", "-c", "ops.f90"],
                   cwd=tmp_path, check=True)
