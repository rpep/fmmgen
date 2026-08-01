"""
generate_code(compress=...) wiring.

Numerical equivalence of the compressed operators is covered symbolically in
test_harmonic_operators.py; what is checked here is that the option reaches
the emitted code intact -- sizes, symbol names, header contract -- and that
turning it off changes nothing.
"""
import re

import pytest

import fmmgen
from fmmgen.harmonic import Nkeep
from fmmgen.utils import Nterms

ORDER = 4
COMMON = dict(
    precision="double", CSE=True, cython=False, potential=True, field=True,
    source_order=0, atomic=False, minpow=11, language="c",
)


def _emit(tmp_path, name, **kw):
    fmmgen.generate_code(
        ORDER, name, include_dir=str(tmp_path), src_dir=str(tmp_path), **{**COMMON, **kw}
    )
    return (tmp_path / f"{name}.h").read_text(), (tmp_path / f"{name}.c").read_text()


def test_option_off_changes_nothing(tmp_path):
    """
    The whole point of the option: with it off, no compressed code path is
    entered and the output matches a build that never knew about it.

    Byte-identity against the pre-change generator is verified out of band
    (regenerate at an earlier commit and diff); what is pinned here is that
    compress=False is stable and carries none of the compressed markers.
    """
    h1, c1 = _emit(tmp_path, "off1", compress=False, harmonic_derivs=True)
    h2, c2 = _emit(tmp_path, "off2", compress=False, harmonic_derivs=True)
    assert h1.replace("off1", "X") == h2.replace("off2", "X")
    assert c1.replace("off1", "X") == c2.replace("off2", "X")

    assert "FMMGEN_COMPRESSED" not in h1
    assert "FMMGEN_MULTIPOLESIZE" not in h1
    assert "M2Lc_" not in h1 and "M2Lc_" not in c1


@pytest.mark.parametrize("harmonic_derivs", [False, True])
def test_compressed_declares_retained_sizes(tmp_path, harmonic_derivs):
    """Sizes must be (p+1)^2, and the header must say so -- callers can no
    longer derive them from Nterms."""
    h, _ = _emit(tmp_path, "on", compress=True, harmonic_derivs=harmonic_derivs)

    assert "#define FMMGEN_COMPRESSED 1" in h
    sizes = [int(v) for v in
             re.search(r"FMMGEN_MULTIPOLESIZE\[\] = \{([^}]*)\}", h).group(1).split(",")]
    assert sizes == [Nkeep(o) for o in range(ORDER + 1)]
    assert sizes[ORDER] == (ORDER + 1) ** 2 < Nterms(ORDER)


def test_compressed_output_is_a_superset(tmp_path):
    """
    Both variants live in one header, so the compressed output must contain
    everything the uncompressed output does, plus a c-suffixed counterpart for
    each expansion operator -- and nothing else. A caller can then hold both
    and switch at runtime with no second header and no redefinition.

    P2P must NOT be duplicated: it never touches an expansion, so the two
    variants would emit identical arithmetic.
    """
    h_u, _ = _emit(tmp_path, "u", compress=False)
    h_c, _ = _emit(tmp_path, "c", compress=True)

    def symbols(h):
        return set(re.findall(r"^void (\w+)\(", h, re.M))

    su, sc = symbols(h_u), symbols(h_c)
    assert su, "no symbols emitted"
    assert su < sc, "compressed output must be a strict superset"

    def suffixed(s):
        base, _, rest = s.partition("_")
        return f"{base}c_{rest}" if rest else f"{base}c"

    expansion_ops = {s for s in su if not s.startswith("P2P")}
    assert sc - su == {suffixed(s) for s in expansion_ops}
    assert not any(s.startswith("P2Pc") for s in sc), "P2P should not be duplicated"


@pytest.mark.parametrize("s", [0, 1])
def test_source_orders_zero_and_one(tmp_path, s):
    """
    Dipole sources (s=1) are the micromagnetics case and must work. The
    multipole then spans s <= |n| <= p and the local |m| <= p - s, so the two
    size tables genuinely differ -- which is the thing to check.
    """
    h, _ = _emit(tmp_path, f"s{s}", compress=True, source_order=s)

    def table(nm):
        return [int(v) for v in
                re.search(rf"{nm}\[\] = \{{([^}}]*)\}}", h).group(1).split(",")]

    msz, lsz = table("FMMGEN_MULTIPOLESIZE"), table("FMMGEN_LOCALSIZE")
    assert msz == [Nkeep(o, s) for o in range(ORDER + 1)]
    assert lsz == [Nkeep(o - s) for o in range(ORDER + 1)]
    if s:
        assert msz != lsz, "multipole and local sizes must differ at s > 0"
    assert f"#define FMMGEN_SOURCEORDER {s}" in h


def test_quadrupole_sources_rejected(tmp_path):
    """At degree 2 the source shell compresses 6 -> 5, so S2M would have to
    project the caller's moment array. Refused rather than guessed."""
    with pytest.raises(NotImplementedError, match="source_order"):
        fmmgen.generate_code(
            6, "quad", include_dir=str(tmp_path), src_dir=str(tmp_path),
            **{**COMMON, "source_order": 2}, compress=True,
        )
