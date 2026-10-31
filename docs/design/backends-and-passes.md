# Backends and passes: design draft

Status: draft for review. Nothing here is implemented beyond the backend registry,
`Problem`, manifests and function prefixes (`fmmgen/backend.py`).

## Goals

- A kernel is a specification; a backend decides how it is evaluated (Cartesian Taylor,
  spherical harmonics, kernel-independent, FFT-based M2L, ...).
- Backends coexist in one program (function prefixes) and can be chosen per order.
- FFT-based M2L is a first-class case, not an exception to the interface.
- Reference tree builders ship as examples; the library defines the contract they fill.

## Operator tiers

1. **Operators**: per-displacement functions (`s2m`, `m2m`, `m2l`, `l2l`, `l2p`, `m2p`, `p2p`).
   Cartesian and spherical-harmonic backends emit these.
2. **Passes**: whole-tree functions. This is the primary contract for tree code.

| Pass | Does | Notes |
|---|---|---|
| `upward` | S2M at leaves, then M2M children to parent | level by level, leaves first |
| `m2l` | far interactions (V-list) | batched per level, accumulated per target box |
| `downward` | L2L parent to children, then L2P at leaves | level by level, root first |
| `near` | P2P over U-list, plus M2P/P2L for W/X lists on adaptive trees | |

A backend with per-offset operators gets generated default passes (loops over the operators,
OpenMP over boxes). A backend with its own transforms implements the passes directly and may emit no
per-offset operators at all.

## Tree description

A plain structure of flat arrays. Anyone's tree builder can fill it; `check_tree` validates it.

- Particles in Morton order, structure of arrays; leaves own contiguous ranges
  (`body_offset`, `nbodies`).
- Per box: `level`, `parent`, `child[2^D]` (or -1), integer lattice coordinates `ijk[D]` at its level,
  centre and half-width.
- Interaction lists in compressed-row form grouped by target box, so one thread owns one target's
  accumulator: `U` (near, leaf-leaf), `V` (same-level far, M2L), and for adaptive trees `W` and `X`.
- Each V entry carries an `offset_class` (an index into a table of integer offset vectors for its level),
  so FFT backends reuse one transformed kernel per class. A tree whose V-list has arbitrary offsets
  (dual-tree traversal with an opening angle) is valid for per-offset backends and flagged as
  `uniform_v = false`; FFT backends reject it with a clear error.

## Context

Backends that need setup declare it in the manifest and implement:

```
ctx = <prefix>_init(config)   // order, root width, max level, threads, kernel parameters
<prefix>_free(ctx)
```

`init` builds FFT plans, per-level kernel transforms (one set for homogeneous kernels, scaled per
level; one per level otherwise), rotation or SVD tables and workspaces. All passes take `ctx`.

## Expansion layout is the backend's business

The tree code sees arrays of doubles with sizes from the manifest: `multipole` and `local` length per
box at order p (`Layout`). It must not index into them. A backend may use multipoles, equivalent densities
on surface grids, or spherical-harmonic coefficients.

## What an FFT backend looks like

`m2l(ctx, level, tree, in, out)`:

1. forward transform each source box's expansion onto the FFT grid (once per source box);
2. for each target box, multiply-accumulate over its V-list in Fourier space with the transformed
   kernel of each offset class;
3. one inverse transform per target box; scatter into the local expansions.

The pass owns its workspace (from `ctx`), is thread-safe across target boxes, and declares the FFT
library in the manifest `dependencies`. Planned candidates: kernel-independent (equivalent surfaces,
pseudo-inverse tables); Cartesian Taylor with the lattice-times-index FFT (prototyped outside this repository).

## Manifest additions

`tiers` (`operator`, `pass`), `runtime_files` (source shipped with the backend, not generated),
`dependencies`, `needs_context`, `requires_uniform_v`, `homogeneity` handling, and for mixed builds
the existing `by_order` and `parts`.

## Open questions

- Batch granularity for FFT passes: per level, or per tile of target boxes for cache and memory.
- GPU: the same pass contract with device pointers, or a separate tier.
- Precision: single and double, and whether FFT backends keep their own precision option.
- Complex-valued kernels (Helmholtz): channel pairs, or native complex arrays.
- How `recommend` works from a tolerance, since order means different accuracy per backend.
