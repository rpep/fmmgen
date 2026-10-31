"""generate_code(orders=..., prefix=...): emit exactly some degrees, under a prefix, without changing them."""
import re

import pytest

import fmmgen

KW = dict(CSE=True, minpow=11, compress=True, planar=True, harmonic_derivs=True, language="c", cache_dir=None)


def _functions(text):
    parts = re.split(r"\n(?=void )", text)
    return {re.match(r"void (\w+)\(", p).group(1): p.strip() for p in parts if p.startswith("void ")}


def test_orders_and_prefix_select_and_rename_without_changing_the_operators(tmp_path):
    for sub in ("all", "some"):
        (tmp_path / sub).mkdir()
    fmmgen.generate_code(5, "ops", src_dir=str(tmp_path / "all"), include_dir=str(tmp_path / "all"), **KW)
    fmmgen.generate_code(5, "ops", orders=[2, 4], prefix="lap_", src_dir=str(tmp_path / "some"),
                         include_dir=str(tmp_path / "some"), **KW)
    whole = _functions((tmp_path / "all" / "ops.c").read_text())
    some = _functions((tmp_path / "some" / "ops.c").read_text())
    assert all(n.startswith("lap_") for n in some)
    # no degree 3, and a degree-4 operator is the same text as in the full build, up to its name
    assert not any(re.search(r"_3$", n) for n in some) and "lap_M2L_4" in some and "lap_M2Lxy_2" in some
    assert some["lap_M2Lc_4"].replace("lap_", "") == whole["M2Lc_4"]
    # the dispatch wrapper switches over the emitted degrees only
    wrapper = some["lap_M2L"]
    assert "case 2:" in wrapper and "case 4:" in wrapper and "case 3:" not in wrapper
    header = (tmp_path / "some" / "ops.h").read_text()
    assert "#define FMMGEN_LAP_MAXORDER 5" in header and "#define FMMGEN_LAP_PLANAR 1" in header
    assert "FMMGEN_MAXORDER" not in header


def test_orders_and_prefix_are_validated(tmp_path):
    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path), **KW)
    with pytest.raises(ValueError, match="orders must be"):
        fmmgen.generate_code(5, "a", orders=[0], **kw)  # below the first degree
    with pytest.raises(ValueError, match="exceed"):
        fmmgen.generate_code(4, "a", orders=[2, 4], **kw)
    with pytest.raises(ValueError, match="identifier"):
        fmmgen.generate_code(3, "a", prefix="bad-name", **kw)
    with pytest.raises(ValueError, match="plain C"):
        fmmgen.generate_code(3, "a", prefix="p_", **{**kw, "language": "fortran"})
