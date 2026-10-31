"""
pytest configuration for fmmgen tests.

This file picks a compiler for pyximport/Cython compilation on macOS,
where Apple's system clang does not support -fopenmp. Linux runners
(including CI) already ship a gcc/g++ with OpenMP support, so CC/CXX
are left alone there.
"""
import os
import platform
import shutil

if platform.system() == "Darwin" and "CC" not in os.environ:
    for gxx in ("g++-15", "g++-14", "g++-13", "g++-12", "g++-11"):
        if shutil.which(gxx):
            os.environ["CC"] = gxx
            os.environ["CXX"] = gxx
            break


import pytest  # noqa: E402


@pytest.fixture(autouse=True)
def _private_temp_dir(tmp_path, monkeypatch):
    """The generation cache defaults to a directory under the system temp dir; give every
    test its own so results never leak between tests or into the user's real cache."""
    import tempfile

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
