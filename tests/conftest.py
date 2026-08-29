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
