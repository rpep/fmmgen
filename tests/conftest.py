"""
pytest configuration for fmmgen tests.

This file configures the test environment to use g++-15 for compiling
Cython extensions via pyximport, ensuring OpenMP support works correctly.
"""
import os

# Configure compiler for pyximport/Cython compilation
# This is needed because Apple's clang doesn't support -fopenmp flag
os.environ['CC'] = 'g++-15'
os.environ['CXX'] = 'g++-15'
