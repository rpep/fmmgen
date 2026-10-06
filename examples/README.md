# Examples

- This directory: the C++ example (OpenMP FMM and Barnes-Hut). Build it with `make`; it can
  be tweaked by editing `example.py`.
- `fortran/`: the same program in Fortran 90, with its own README.

Source Order:- 0 for Charges, 1 for Dipoles, 2 for Quadrupoles, etc.

Pass `--input FILE` to read particles (one per line: positions then source strengths, comma
separated) instead of generating them; the file format is the one written to
`particles_n_<N>.txt`. Both examples read it, so the same file can be given to each to compare
results.
