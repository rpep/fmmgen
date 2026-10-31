"""
Spherical-harmonic (solid-harmonic) Laplace FMM operators, numpy reference.

Kernel 1/|r|, charge sources, potential and field targets. Operator names,
argument lists and sign conventions follow the generated Cartesian C
operators: all displacement arguments are "target minus source" style,

    S2M(x, S, M)    x = centre - source              S = [q]
    M2M(x, M, Ms)   x = parent centre - child centre
    M2L(x, M, L)    x = local centre - multipole centre
    L2L(x, L, Ls)   x = child centre - parent centre
    L2P(x, L, F)    x = target - local centre         F = [phi, -dphi/dx_i]
    M2P(x, M, F)    x = target - multipole centre
    P2P(x, S, F)    x = target - source

and every operator accumulates (+=) into its output.

Mathematics (complex form, Y-normalisation chosen so translations are plain
convolutions; n = degree, m = order, -n <= m <= n):

    R_n^m(r) = r^n P_n^m(cos t) e^{i m f} / (n+m)!         regular solid harmonic
    I_n^m(r) = (n-m)! P_n^m(cos t) e^{i m f} / r^(n+1)     irregular solid harmonic
    1/|r - r'| = sum_{n,m} conj(R_n^m(r')) I_n^m(r)        (|r'| < |r|)

    M_n^m = sum_s q_s conj(R_n^m(s - c))      phi(t) = sum M_n^m I_n^m(t - c)
    L_n^m                                     phi(t) = sum L_n^m R_n^m(t - c)

    M2M: Ms_n^m = sum_{j,k} conj(R_j^k(c - c')) M_{n-j}^{m-k}
    L2L: Ls_n^m = sum_{j,k} L_j^k R_{j-n}^{k-m}(c_child - c_parent)
    M2L: L_n^m  = (-1)^(n+m) sum_{j,k} M_j^k I_{j+n}^{k-m}(c_L - c_M)

Real-valued sources give M_n^{-m} = (-1)^m conj(M_n^m) (same for L), so only
m >= 0 is stored, as two reals (Re, Im) per (n, m>0). Layout inside degree n
is real-SH order mu = -n..n:

    slot n^2 + n + m   holds Re(c_n^m)   (m = 0..n)
    slot n^2 + n - m   holds Im(c_n^m)   (m = 1..n)

so a degree-p expansion has (p+1)^2 reals (same as the trace-free Cartesian
basis) and every degree is a contiguous real-SH block, which is what the
rotation-based M2L needs.

Two M2L variants:
  * plain (`m2l_plain`): the O(p^4) double sum over (j,k).
  * rotation (`m2l_rotation`): rotate the multipole to the z axis, shift along
    z (a single sum per order m), rotate back: O(p^3). Rotation matrices of
    real spherical harmonics are built by the Ivanic-Ruedenberg recurrence
    (J. Phys. Chem. 100, 6342 (1996), with the 1998 erratum), which needs only
    the 3x3 rotation matrix and has no pole singularity. Rotation-based
    translation: Cheng, Greengard, Rokhlin, J. Comput. Phys. 155, 468 (1999);
    White and Head-Gordon, J. Chem. Phys. 101, 6593 (1994).

Truncation: trunc="box" keeps all (j, n) <= p in M2L (classic); trunc="total"
keeps j + n <= p (same asymptotic error order, about a third of the work).
"""
from functools import lru_cache
from math import factorial, pi, sqrt

import numpy as np


def sizes(p):
    """(len(M), len(L)) for expansion degree p."""
    return (p + 1) ** 2, (p + 1) ** 2


def planar_keep(p):
    """Full-layout indices of the entries a planar expansion keeps, in planar order.

    Sources, targets and centres in the z = 0 plane: the solid harmonics with n + |m| odd vanish there, and so do
    the multipole and local coefficients built from them. A planar array holds the n + |m| even entries, n + 1 of
    degree n, at n(n+1)/2 + (n+k)/2 for the signed index k (k > 0 the real part of order k, k < 0 the imaginary
    part of order |k|)."""
    return [n * n + n + k for n in range(p + 1) for k in range(-n, n + 1, 2)]


def planar_size(p):
    """Length of a planar expansion array of degree p: (p+1)(p+2)/2."""
    return (p + 1) * (p + 2) // 2


def ire(n, m):
    return n * n + n + m


def iim(n, m):
    return n * n + n - m


# ---------------------------------------------------------------------------
# Solid harmonics, flat real layout
# ---------------------------------------------------------------------------
def _solid(p, x, y, z, irregular):
    rr = x * x + y * y + z * z
    w = x + 1j * y
    A = {}
    if irregular:
        A[(0, 0)] = 1.0 / sqrt(rr) + 0j
    else:
        A[(0, 0)] = 1.0 + 0j
    for m in range(p + 1):
        if m > 0:
            if irregular:
                A[(m, m)] = -(2 * m - 1) * w / rr * A[(m - 1, m - 1)]
            else:
                A[(m, m)] = -w / (2 * m) * A[(m - 1, m - 1)]
        if m + 1 <= p:
            A[(m + 1, m)] = ((2 * m + 1) * z / rr if irregular else z) * A[(m, m)]
        for n in range(m + 2, p + 1):
            if irregular:
                A[(n, m)] = ((2 * n - 1) * z * A[(n - 1, m)] - ((n - 1) ** 2 - m * m) * A[(n - 2, m)]) / rr
            else:
                A[(n, m)] = ((2 * n - 1) * z * A[(n - 1, m)] - rr * A[(n - 2, m)]) / (n * n - m * m)
    out = np.zeros((p + 1) ** 2)
    for n in range(p + 1):
        out[ire(n, 0)] = A[(n, 0)].real
        for m in range(1, n + 1):
            out[ire(n, m)] = A[(n, m)].real
            out[iim(n, m)] = A[(n, m)].imag
    return out


def solid_R(p, x, y, z):
    """Regular solid harmonics R_n^m(x,y,z), n <= p, flat layout."""
    return _solid(p, x, y, z, False)


def solid_I(p, x, y, z):
    """Irregular solid harmonics I_n^m(x,y,z), n <= p, flat layout."""
    return _solid(p, x, y, z, True)


# ---------------------------------------------------------------------------
# Bilinear term lists: out[o] += c * A[a] * B[b]
# ---------------------------------------------------------------------------
def cget(n, m):
    """Complex coefficient c_n^m (any sign of m) as ((s, idx) | None, (s, idx) | None)
    giving the real and imaginary part as signed references into the flat array."""
    if n < 0 or abs(m) > n:
        return None, None
    if m >= 0:
        return (1, ire(n, m)), ((1, iim(n, m)) if m > 0 else None)
    mm = -m
    sg = -1 if mm % 2 else 1
    return (sg, ire(n, mm)), (-sg, iim(n, mm))


def _cmul(terms, out, coef, A, B, conjA=False):
    """terms += complex product coef * A * B written into out = (o_re, o_im)."""
    (Ar, Ai), (Br, Bi) = A, B
    if conjA and Ai is not None:
        Ai = (-Ai[0], Ai[1])
    o_re, o_im = out

    def add(o, c, a, b):
        if o is not None and a is not None and b is not None:
            terms.append((o, c * a[0] * b[0], a[1], b[1]))

    add(o_re, coef, Ar, Br)
    add(o_re, -coef, Ai, Bi)
    add(o_im, coef, Ar, Bi)
    add(o_im, coef, Ai, Br)


def _out(n, m):
    return ire(n, m), (iim(n, m) if m > 0 else None)


def _pack(terms):
    if not terms:
        z = np.zeros(0, dtype=int)
        return z, np.zeros(0), z, z
    acc = {}
    for o, c, a, b in terms:
        acc[(o, a, b)] = acc.get((o, a, b), 0.0) + c
    keys = [k for k, v in acc.items() if v != 0.0]
    o = np.array([k[0] for k in keys], dtype=int)
    a = np.array([k[1] for k in keys], dtype=int)
    b = np.array([k[2] for k in keys], dtype=int)
    c = np.array([acc[k] for k in keys])
    return o, c, a, b


@lru_cache(maxsize=None)
def terms_S2M(p):
    """M[o] += c * R(x)[a] * S[0]."""
    t = []
    for n in range(p + 1):
        sg = -1.0 if n % 2 else 1.0
        t.append((ire(n, 0), sg, ire(n, 0), 0))
        for m in range(1, n + 1):
            t.append((ire(n, m), sg, ire(n, m), 0))
            t.append((iim(n, m), -sg, iim(n, m), 0))
    return _pack(t)


@lru_cache(maxsize=None)
def terms_M2M(p):
    """Ms[o] += c * R(x)[a] * M[b]."""
    t = []
    for n in range(p + 1):
        for m in range(n + 1):
            for j in range(n + 1):
                sg = -1.0 if j % 2 else 1.0
                for k in range(max(-j, m - (n - j)), min(j, m + n - j) + 1):
                    _cmul(t, _out(n, m), sg, cget(j, k), cget(n - j, m - k), conjA=True)
    return _pack(t)


@lru_cache(maxsize=None)
def terms_L2L(p):
    """Ls[o] += c * R(x)[a] * L[b]."""
    t = []
    for n in range(p + 1):
        for m in range(n + 1):
            for j in range(n, p + 1):
                for k in range(max(-j, m - (j - n)), min(j, m + j - n) + 1):
                    _cmul(t, _out(n, m), 1.0, cget(j - n, k - m), cget(j, k))
    return _pack(t)


@lru_cache(maxsize=None)
def terms_M2L(p, trunc="box"):
    """L[o] += c * I(x)[a] * M[b]; I needed to degree 2p (box) or p (total)."""
    t = []
    for n in range(p + 1):
        for m in range(n + 1):
            sg = -1.0 if (n + m) % 2 else 1.0
            for j in range(p + 1 - (n if trunc == "total" else 0)):
                for k in range(-j, j + 1):
                    _cmul(t, _out(n, m), sg, cget(j + n, k - m), cget(j, k))
    return _pack(t)


def degree_I_M2L(p, trunc="box"):
    return p if trunc == "total" else 2 * p


def _eval_terms(P, A, B, nout):
    o, c, a, b = P
    return np.bincount(o, weights=c * A[a] * B[b], minlength=nout)


@lru_cache(maxsize=None)
def terms_eval(p, kind):
    """kind 'L2P' (A = R, B = L) or 'M2P' (A = I, B = M). F[o] += c*A[a]*B[b].
    A = (R or I) array of degree p (L2P) or p+1 (M2P)."""
    t = []
    L = kind == "L2P"
    for n in range(p + 1):
        # potential: c_n^0 A_n^0 + 2 Re sum_{m>0} c_n^m A_n^m
        t.append((0, 1.0, ire(n, 0), ire(n, 0)))
        for m in range(1, n + 1):
            t.append((0, 2.0, ire(n, m), ire(n, m)))
            t.append((0, -2.0, iim(n, m), iim(n, m)))
        # dz: L2P F3 = -sum_m L_n^m R_{n-1}^m ; M2P F3 = +sum_m M_n^m I_{n+1}^m
        nn = n - 1 if L else n + 1
        sgz = -1.0 if L else 1.0
        if nn >= 0:
            t.append((3, sgz, ire(nn, 0), ire(n, 0)))
            for m in range(1, min(n, nn) + 1):
                t.append((3, 2.0 * sgz, ire(nn, m), ire(n, m)))
                t.append((3, -2.0 * sgz, iim(nn, m), iim(n, m)))
        # T = sum_m c_n^m A_{nn}^{m-1};  F1 = Re T, F2 = -Im T
        if nn >= 0:
            for m in range(-n, n + 1):
                tt = []
                _cmul(tt, (1, 2), 1.0, cget(nn, m - 1), cget(n, m))
                for o, c, a, b in tt:
                    t.append((1, c, a, b) if o == 1 else (2, -c, a, b))
    return _pack(t)


# ---------------------------------------------------------------------------
# Operators (numpy reference)
# ---------------------------------------------------------------------------
def S2M(p, x, y, z, S, M):
    R = solid_R(p, x, y, z)
    o, c, a, _ = terms_S2M(p)
    M += np.bincount(o, weights=c * R[a] * S[0], minlength=len(M))


def M2M(p, x, y, z, M, Ms):
    R = solid_R(p, x, y, z)
    Ms += _eval_terms(terms_M2M(p), R, M, len(Ms))


def L2L(p, x, y, z, L, Ls):
    R = solid_R(p, x, y, z)
    Ls += _eval_terms(terms_L2L(p), R, L, len(Ls))


def m2l_plain(p, x, y, z, M, L, trunc="box"):
    I = solid_I(degree_I_M2L(p, trunc), x, y, z)
    L += _eval_terms(terms_M2L(p, trunc), I, M, len(L))


def L2P(p, x, y, z, L, F):
    R = solid_R(p, x, y, z)
    F += _eval_terms(terms_eval(p, "L2P"), R, L, 4)


def M2P(p, x, y, z, M, F):
    I = solid_I(p + 1, x, y, z)
    F += _eval_terms(terms_eval(p, "M2P"), I, M, 4)


def P2P(x, y, z, S, F):
    r2 = x * x + y * y + z * z
    ri = 1.0 / sqrt(r2)
    F[0] += ri * S[0]
    t = ri ** 3 * S[0]
    F[1] += t * x
    F[2] += t * y
    F[3] += t * z


# ---------------------------------------------------------------------------
# Rotation-based M2L
# ---------------------------------------------------------------------------
def rotation_to_z(x, y, z):
    """3x3 rotation Q with Q (x,y,z)/|.| = (0,0,1). Returns (Q, rho)."""
    rho = sqrt(x * x + y * y + z * z)
    rxy = sqrt(x * x + y * y)
    if rxy < 1e-14 * rho:
        ca, sa = 1.0, 0.0
    else:
        ca, sa = x / rxy, y / rxy
    cb, sb = z / rho, rxy / rho
    Q = np.array([[cb * ca, cb * sa, -sb], [-sa, ca, 0.0], [sb * ca, sb * sa, cb]])
    return Q, rho


def ir_coefficients(p):
    """Ivanic-Ruedenberg u, v, w tables: dict l -> (u, v, w) arrays [m+l, n+l]."""
    out = {}
    for l in range(2, p + 1):
        u = np.zeros((2 * l + 1, 2 * l + 1))
        v = np.zeros_like(u)
        w = np.zeros_like(u)
        for m in range(-l, l + 1):
            d = 1 if m == 0 else 0
            for n in range(-l, l + 1):
                den = (2 * l) * (2 * l - 1) if abs(n) == l else (l + n) * (l - n)
                u[m + l, n + l] = sqrt((l + m) * (l - m) / den)
                v[m + l, n + l] = 0.5 * sqrt((1 + d) * (l + abs(m) - 1) * (l + abs(m)) / den) * (1 - 2 * d)
                w[m + l, n + l] = -0.5 * sqrt((l - abs(m) - 1) * (l - abs(m)) / den) * (1 - d)
        out[l] = (u, v, w)
    return out


def ir_rotations(Q, p, coeffs=None):
    """Real-SH rotation matrices D_l (list, l = 0..p) with S_l(Q r) = D_l S_l(r),
    mu = -l..l order, orthonormal real harmonics, Ivanic-Ruedenberg recurrence."""
    if coeffs is None:
        coeffs = ir_coefficients(p)
    Pm = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    R1 = Pm @ Q @ Pm.T
    Ds = [np.array([[1.0]]), R1]

    def P(i, l, a, b, Rp):
        ri1, rim1, ri0 = R1[i + 1, 2], R1[i + 1, 0], R1[i + 1, 1]
        if b == -l:
            return ri1 * Rp[a + l - 1, 0] + rim1 * Rp[a + l - 1, 2 * l - 2]
        if b == l:
            return ri1 * Rp[a + l - 1, 2 * l - 2] - rim1 * Rp[a + l - 1, 0]
        return ri0 * Rp[a + l - 1, b + l - 1]

    for l in range(2, p + 1):
        Rp = Ds[l - 1]
        u, v, w = coeffs[l]
        Rl = np.zeros((2 * l + 1, 2 * l + 1))
        for m in range(-l, l + 1):
            for n in range(-l, l + 1):
                val = 0.0
                if u[m + l, n + l] != 0:
                    val += u[m + l, n + l] * P(0, l, m, n, Rp)
                if v[m + l, n + l] != 0:
                    if m == 0:
                        V = P(1, l, 1, n, Rp) + P(-1, l, -1, n, Rp)
                    elif m > 0:
                        dd = 1 if m == 1 else 0
                        V = P(1, l, m - 1, n, Rp) * sqrt(1 + dd) - P(-1, l, -m + 1, n, Rp) * (1 - dd)
                    else:
                        dd = 1 if m == -1 else 0
                        V = P(1, l, m + 1, n, Rp) * (1 - dd) + P(-1, l, -m - 1, n, Rp) * sqrt(1 + dd)
                    val += v[m + l, n + l] * V
                if w[m + l, n + l] != 0:
                    if m > 0:
                        W = P(1, l, m + 1, n, Rp) + P(-1, l, -m - 1, n, Rp)
                    else:
                        W = P(1, l, m - 1, n, Rp) - P(-1, l, -m + 1, n, Rp)
                    val += w[m + l, n + l] * W
                Rl[m + l, n + l] = val
        Ds.append(Rl)
    return Ds[: p + 1]


def slot_scale(p):
    """Per-slot factor mapping stored (Re, Im) to orthonormal real-SH coefficients:
    mu = 0 -> 1, mu = +m -> sqrt2 (-1)^m (Re), mu = -m -> -sqrt2 (-1)^m (Im)."""
    s = np.ones((p + 1) ** 2)
    for n in range(p + 1):
        for m in range(1, n + 1):
            sg = -1.0 if m % 2 else 1.0
            s[ire(n, m)] = sqrt(2) * sg
            s[iim(n, m)] = -sqrt(2) * sg
    return s


def basis_scales(p):
    """(sb, sa): stored M -> orthonormal real-SH coefficients is  c = sb * M;
    orthonormal coefficients of a local expansion -> stored L is  L = c / sa.
    sb = slot_scale * b_nm,  sa = slot_scale * a_nm  with
    a_nm = kappa_n / sqrt((n+m)!(n-m)!),  b_nm = kappa_n sqrt((n+m)!(n-m)!),
    kappa_n = sqrt(4 pi / (2n+1)) (R_n^m = a_nm r^n Y_n^m, I_n^m = b_nm Y_n^m / r^(n+1))."""
    sc = slot_scale(p)
    a = np.zeros((p + 1) ** 2)
    b = np.zeros((p + 1) ** 2)
    for n in range(p + 1):
        kap = sqrt(4 * pi / (2 * n + 1))
        for m in range(-n, n + 1):
            f = sqrt(factorial(n + m) * factorial(n - m))
            a[ire(n, m)] = kap / f
            b[ire(n, m)] = kap * f
    return sc * b, sc * a


def zshift_table(p):
    """Z[n, j, m] = kappa_n/kappa_j (-1)^(n+m) a-b scale, without (j+n)!/rho^(j+n+1):
    a_nm = kappa_n/sqrt((n+m)!(n-m)!), b_jm = kappa_j sqrt((j+m)!(j-m)!);
    returns (a_nm/b_jm) (-1)^(n+m) (j+n)! (the rho powers are applied at run time)."""
    Z = np.zeros((p + 1, p + 1, p + 1))
    for n in range(p + 1):
        for j in range(p + 1):
            for m in range(min(n, j) + 1):
                a = sqrt((2 * j + 1) / (2 * n + 1)) / sqrt(factorial(n + m) * factorial(n - m))
                b = sqrt(factorial(j + m) * factorial(j - m))
                Z[n, j, m] = (-1.0 if (n + m) % 2 else 1.0) * a / b * factorial(j + n)
    return Z


def rotate_blocks(Ds, v, transpose=False):
    out = np.zeros_like(v)
    for n, D in enumerate(Ds):
        sl = slice(n * n, (n + 1) ** 2)
        out[sl] = (D.T if transpose else D) @ v[sl]
    return out


def m2l_rotation(p, x, y, z, M, L, trunc="box", Ds=None):
    """O(p^3) M2L through rotation to the z axis. Ds (optional): precomputed
    rotation matrices for this direction (see ir_rotations)."""
    Q, rho = rotation_to_z(x, y, z)
    if Ds is None:
        Ds = ir_rotations(Q, p)
    sb, sa = basis_scales(p)
    cM = rotate_blocks(Ds, M * sb)
    Z = zshift_table(p)
    N = np.arange(2 * p + 2)
    rp = 1.0 / rho ** (N + 1.0)
    cL = np.zeros_like(cM)
    for n in range(p + 1):
        for m in range(n + 1):
            jmax = p - n if trunc == "total" else p
            for j in range(m, jmax + 1):
                zz = Z[n, j, m] * rp[j + n]
                cL[ire(n, m)] += zz * cM[ire(j, m)]
                if m > 0:
                    cL[iim(n, m)] += zz * cM[iim(j, m)]
    L += rotate_blocks(Ds, cL, transpose=True) / sa


def M2L(p, x, y, z, M, L, variant="rotation", trunc="box"):
    if variant == "plain":
        m2l_plain(p, x, y, z, M, L, trunc)
    else:
        m2l_rotation(p, x, y, z, M, L, trunc)
