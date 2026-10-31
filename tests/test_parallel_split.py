"""
generate_code(jobs=..., split=..., cache_dir=...): faster generation and
splittable output that must not change what is generated.
"""
import glob
from pathlib import Path
import re

import pytest

import fmmgen

ORDER = 4
COMMON = dict(
    precision="double", CSE=True, cython=False, potential=True, field=True,
    source_order=0, atomic=False, minpow=11, harmonic_derivs=True, compress=True,
    planar=True, language="c++",
)


def _emit(tmp_path, sub, **kw):
    d = tmp_path / sub
    d.mkdir()
    fmmgen.generate_code(ORDER, "ops", include_dir=str(d), src_dir=str(d), **{**COMMON, **kw})
    return d


def _functions(text):
    """Map function name -> full definition text, for the generated .cpp text."""
    parts = re.split(r"\n(?=void )", text)
    return {re.match(r"void (\w+)\(", p).group(1): p.strip() for p in parts if p.startswith("void ")}


def test_jobs_output_identical(tmp_path, monkeypatch):
    # cache_dir=None: with the default cache the second run would be a cache hit and never reach the pool
    import multiprocessing as mp

    pools = []
    real = mp.get_context

    def spy(*args, **kw):
        pools.append(1)
        return real(*args, **kw)

    monkeypatch.setattr(mp, "get_context", spy)
    a = _emit(tmp_path, "serial", cache_dir=None)
    assert not pools
    b = _emit(tmp_path, "parallel", jobs=2, cache_dir=None)
    assert pools, "jobs=2 did not use a worker pool"
    assert (a / "ops.h").read_text() == (b / "ops.h").read_text()
    assert (a / "ops.cpp").read_text() == (b / "ops.cpp").read_text()


def test_jobs_is_ignored_for_fortran(tmp_path, caplog):
    kw = dict(CSE=True, minpow=11, language="fortran", cache_dir=None)
    for sub in ("ser", "par"):
        (tmp_path / sub).mkdir()
    fmmgen.generate_code(3, "ops", src_dir=str(tmp_path / "ser"), **kw)
    with caplog.at_level("WARNING"):
        fmmgen.generate_code(3, "ops", src_dir=str(tmp_path / "par"), jobs=2, **kw)
    assert "jobs is ignored" in caplog.text
    assert (tmp_path / "ser" / "ops.f90").read_text() == (tmp_path / "par" / "ops.f90").read_text()


def test_split_has_same_functions(tmp_path):
    a = _emit(tmp_path, "whole")
    b = _emit(tmp_path, "split", split=True)
    assert (a / "ops.h").read_text() == (b / "ops.h").read_text()
    whole = _functions((a / "ops.cpp").read_text())
    parts = sorted(glob.glob(str(b / "ops_*.cpp")))
    # One file per operator, named after it.
    wrappers = set(_functions((b / "ops.cpp").read_text()))
    assert sorted(Path(fn).stem[len("ops_"):] for fn in parts) == sorted(set(whole) - wrappers)
    assert "ops_M2L_3.cpp" in {Path(fn).name for fn in parts}
    pieces = {}
    for fn in [str(b / "ops.cpp")] + parts:
        got = _functions(open(fn).read())
        assert not (set(got) & set(pieces)), "a function is defined in two files"
        pieces.update(got)
    assert pieces == whole
    # Dispatch wrappers stay in the main file; every file includes the header.
    assert "void M2L(" in (b / "ops.cpp").read_text()
    for fn in parts:
        assert '#include "ops.h"' in open(fn).read()


def test_split_removes_stale_files(tmp_path):
    d = _emit(tmp_path, "d", split=True)
    (d / "ops_part0.cpp").write_text("// old layout\n")
    (d / "ops_notes.cpp").write_text("// not an operator file\n")
    fmmgen.generate_code(ORDER, "ops", include_dir=str(d), src_dir=str(d), split=False, **COMMON)
    assert sorted(p.name for p in d.glob("*.cpp")) == ["ops.cpp", "ops_notes.cpp"]


def test_cache_reuses_orders(tmp_path):
    cache = tmp_path / "cache"
    a = _emit(tmp_path, "a", cache_dir=str(cache))
    n = len(list(cache.glob("*.json")))
    assert n > 0
    b = _emit(tmp_path, "b", cache_dir=str(cache), jobs=2)
    assert len(list(cache.glob("*.json"))) == n
    assert (a / "ops.cpp").read_text() == (b / "ops.cpp").read_text()


def test_split_rejects_cython(tmp_path):
    with pytest.raises(NotImplementedError):
        fmmgen.generate_code(2, "x", include_dir=str(tmp_path), src_dir=str(tmp_path),
                             split=True, **{**COMMON, "cython": True})


def test_default_cache_is_private_and_used(tmp_path):
    import stat

    a = _emit(tmp_path, "a", **{})  # cache_dir defaults to "auto"
    cache = [d for d in tmp_path.iterdir() if d.name.startswith("fmmgen-")]
    assert len(cache) == 1 and not stat.S_IMODE(cache[0].stat().st_mode) & 0o077
    assert list(cache[0].glob("*.json"))
    b = _emit(tmp_path, "b")
    assert (a / "ops.cpp").read_text() == (b / "ops.cpp").read_text()
    assert not list(cache[0].glob("*.tmp"))


def test_cache_off_and_corrupt_entries(tmp_path):
    off = _emit(tmp_path, "off", cache_dir=None)
    assert not [d for d in tmp_path.iterdir() if d.name.startswith("fmmgen-")]
    cache = tmp_path / "cache"
    a = _emit(tmp_path, "a", cache_dir=str(cache))
    for f in cache.glob("*.json"):
        f.write_text("{not json")  # damaged entries are regenerated, not fatal
    b = _emit(tmp_path, "b", cache_dir=str(cache))
    assert (a / "ops.cpp").read_text() == (b / "ops.cpp").read_text() == (off / "ops.cpp").read_text()


def test_foreign_cache_dir_is_ignored(tmp_path):
    import os

    if not hasattr(os, "getuid"):
        pytest.skip("POSIX only")
    d = tmp_path / f"fmmgen-{os.getuid()}"
    d.mkdir(mode=0o777)
    os.chmod(d, 0o777)  # writable by others: must not be used
    _emit(tmp_path, "a")
    assert not list(d.glob("*.json"))


def test_fortran_ignores_warm_cache(tmp_path):
    kw = dict(CSE=True, minpow=11, language="fortran", cache_dir=str(tmp_path / "cache"))
    for d in ("first", "second"):
        (tmp_path / d).mkdir()
        fmmgen.generate_code(3, "ops", src_dir=str(tmp_path / d), **kw)
    first, second = (tmp_path / d / "ops.f90" for d in ("first", "second"))
    assert "subroutine M2L(" in second.read_text()  # the dispatch wrappers need the printer's signatures
    assert first.read_text() == second.read_text()
    assert not (tmp_path / "cache").exists()
