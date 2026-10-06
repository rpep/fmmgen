# Fortran 90 example

Fortran counterpart of the C++ example in the parent directory, with the same
algorithms, options and output files:

- FMM (`--type 0`) and Barnes-Hut (`--type 1`), plus the direct reference solution
- 3D octree with the full operators, or the harmonic-compressed ones (`--compress 1`)
- 2D planar quadtree (`--dim 2`), which always uses the planar operators
- OpenMP over exclusive targets and level-synchronous M2M/L2L; results are identical for any thread count
- `--input FILE` to read particles instead of generating them (see below)
- `refresh_sources` to reuse a tree with new source strengths

```
make                                   # runs example.py, then builds ./main
./main --nparticles 5000 --type 0
./main --help
```

`operators.f90` is generated free-form Fortran 90 (one module, `implicit none`) and
builds with `-std=f2003 -Wall -Wextra`. OpenMP is optional: drop `-fopenmp` from `FFLAGS`.
`example.py` must keep `compress=True` and `planar=True`, because the variant is chosen
at run time (the C++ example selects it with `#ifdef`; Fortran has no equivalent).
`compress`, `planar` and `source_order` work; `cython`, `atomic`, `gpu` and single
precision do not.

Input particles: `--input FILE` reads one particle per line, comma separated -- the
positions (2 or 3) followed by the source strengths, the format of the
`particles_n_<N>.txt` file both programs write. The particle count comes from the file
and the file is never rewritten. Running both programs on the same file removes the
random-number difference: fields from the C++ and Fortran examples then agree to
about 1e-15 relative to the largest field value.

Differences from the C++ example:

- Without `--input`, particles come from `random_number`, so errors are similar but not identical.
- The program also prints the L2 relative error next to the mean relative error.
- The C++ program pins threads to physical cores on Linux. Fortran has no portable way
  to do that: set `OMP_PROC_BIND=close OMP_PLACES=cores` instead.
