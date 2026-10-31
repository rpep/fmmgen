"""
C code generation for the spherical-harmonic Laplace FMM operators.

The mathematics, layouts and sign conventions are documented in
`fmmgen.sh_math` (the numpy reference). This module emits

    S2M_p  M2M_p  M2L_p  L2L_p  L2P_p  M2P_p  P2P  P2P_batch

for p = pmin..order-1, with the same argument lists as the Cartesian operators
from `fmmgen.generate_code` (expansion arrays have (p+1)^2 reals, see
`fmmgen.sh_math.sizes`; the header carries FMMGEN_MULTIPOLESIZE/LOCALSIZE).
`order` follows the generate_code convention: degrees p = pmin..order-1 are
emitted and FMMGEN_MAXORDER = order is an exclusive bound (so the example tree
driver, which loops `order < FMMGEN_MAXORDER`, runs unchanged).

In addition to M2L_p the file always provides both M2L variants:

    M2Lplain_p(x,y,z,M,L)       O(p^4) double sum
    M2Lrot_p(x,y,z,M,L)         O(p^3) rotate-shift-rotate, rotation built on the fly
    M2Lrot_setup_p(x,y,z,Dt)    build the rotation tables for one direction
    M2Lrot_apply_p(Dt,M,L)      apply them (Dt holds 1/rho then the rotation blocks;
                                its length is FMMGEN_ROTSIZE[p]); in a uniform tree
                                the same table serves every box pair with the same
                                offset, so this is the realistic per-pair cost

`m2l="plain"|"rotation"` selects what M2L_p calls (default plain: with the rotation built on the
fly it only wins above p of about 20, use the setup/apply pair). `trunc="box"` keeps all
(j, n) <= p in M2L, `trunc="total"` keeps j + n <= p.

unroll=False: compact loop code, runtime-p generic routines plus static const
coefficient tables (Z-shift, Ivanic-Ruedenberg u/v/w, basis scales); scales to
any order (code size grows like p^3 for the tables only).
unroll=True: straight-line code with every index and coefficient resolved at
generation time, like the Cartesian output. Code size grows like p^5 for the
plain M2L; `max_unroll` (default 10) refuses larger orders. The rotation
builder (Ivanic-Ruedenberg recurrence) stays a loop in both modes.
"""
import os

import numpy as np

from . import sh_math as sm

DEFAULT_MAX_UNROLL = 10

_RESTRICT = """\
#if defined(__cplusplus)
#  define FMMGEN_RESTRICT __restrict
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#  define FMMGEN_RESTRICT restrict
#else
#  define FMMGEN_RESTRICT
#endif
"""


def rot_size(p):
    """Doubles in one rotation table: 1 (1/rho) + sum_{n<=p} (2n+1)^2."""
    return 1 + sum((2 * n + 1) ** 2 for n in range(p + 1))


def _num(v):
    r = repr(float(v))
    return r


def _carray(ctype, name, vals, per=6):
    items = [_num(v) if ctype == "double" else str(int(v)) for v in vals]
    lines = []
    for i in range(0, len(items), per):
        lines.append("    " + ", ".join(items[i:i + per]))
    return f"static const {ctype} {name}[{len(items)}] = {{\n" + ",\n".join(lines) + "\n};\n"


# ---------------------------------------------------------------------------
# Tables (loop variant)
# ---------------------------------------------------------------------------
def _tables(P, full=True):
    out = []
    if full:
        sb, sa = sm.basis_scales(P)
        out.append(_carray("double", "SH_SB", sb))
        out.append(_carray("double", "SH_SAI", 1.0 / sa))
        Z = sm.zshift_table(P)
        # index (n*(P+1) + j)*(P+1) + m
        out.append(_carray("double", "SH_ZT", Z.reshape(-1)))
    # IR coefficients by (l, |m|, |n|)
    off = [0] * (P + 2)
    U, V, W = [], [], []
    for l in range(P + 1):
        off[l] = len(U)
        if l < 2:
            continue
        for am in range(l + 1):
            for an in range(l + 1):
                m, n = am, an
                d = 1 if m == 0 else 0
                den = (2 * l) * (2 * l - 1) if n == l else (l + n) * (l - n)
                U.append(np.sqrt((l + m) * (l - m) / den))
                V.append(0.5 * np.sqrt((1 + d) * (l + m - 1) * (l + m) / den) * (1 - 2 * d))
                W.append(-0.5 * np.sqrt((l - m - 1) * (l - m) / den) * (1 - d))
    off[P + 1] = len(U)
    out.append(_carray("int", "SH_IROFF", off[: P + 2]))
    out.append(_carray("double", "SH_IRU", U))
    out.append(_carray("double", "SH_IRV", V))
    out.append(_carray("double", "SH_IRW", W))
    roff = [0]
    for l in range(P + 1):
        roff.append(roff[-1] + (2 * l + 1) ** 2)
    out.append(_carray("int", "SH_ROTOFF", roff))
    return "\n".join(out)


_LOOP_CODE = r"""
#define SH_NS ((SH_PMAX + 1) * (SH_PMAX + 1))
#define SH_IRE(n, m) ((n) * (n) + (n) + (m))
#define SH_IIM(n, m) ((n) * (n) + (n) - (m))
#define SH_MAX(a, b) ((a) > (b) ? (a) : (b))
#define SH_MIN(a, b) ((a) < (b) ? (a) : (b))

/* Solid harmonics, flat layout. irr = 0: regular R_n^m, irr = 1: irregular I_n^m. */
static void sh_solid(int p, int irr, double x, double y, double z, double *A) {
    double r2 = x * x + y * y + z * z;
    double rinv2 = 1.0 / r2;
    A[0] = irr ? sqrt(rinv2) : 1.0;
    for (int m = 0; m <= p; m++) {
        if (m > 0) {
            double pre = A[SH_IRE(m - 1, m - 1)];
            double pim = (m > 1) ? A[SH_IIM(m - 1, m - 1)] : 0.0;
            double f = irr ? -(2 * m - 1) * rinv2 : -1.0 / (2 * m);
            A[SH_IRE(m, m)] = f * (x * pre - y * pim);
            A[SH_IIM(m, m)] = f * (x * pim + y * pre);
        }
        if (m + 1 <= p) {
            double f = irr ? (2 * m + 1) * z * rinv2 : z;
            A[SH_IRE(m + 1, m)] = f * A[SH_IRE(m, m)];
            if (m > 0) A[SH_IIM(m + 1, m)] = f * A[SH_IIM(m, m)];
        }
        for (int n = m + 2; n <= p; n++) {
            double a = (2 * n - 1) * z;
            double b = irr ? (double)((n - 1) * (n - 1) - m * m) : r2;
            double c = irr ? rinv2 : 1.0 / (double)(n * n - m * m);
            A[SH_IRE(n, m)] = (a * A[SH_IRE(n - 1, m)] - b * A[SH_IRE(n - 2, m)]) * c;
            if (m > 0) A[SH_IIM(n, m)] = (a * A[SH_IIM(n - 1, m)] - b * A[SH_IIM(n - 2, m)]) * c;
        }
    }
}

/* Expand to full complex rows k = -n..n: F[n * stride + off + k], off >= p. */
static void sh_expand(int p, int stride, int off, const double *A, double *Fr, double *Fi) {
    for (int n = 0; n <= p; n++) {
        double *fr = Fr + n * stride + off;
        double *fi = Fi + n * stride + off;
        fr[0] = A[SH_IRE(n, 0)];
        fi[0] = 0.0;
        for (int m = 1; m <= n; m++) {
            double re = A[SH_IRE(n, m)], im = A[SH_IIM(n, m)];
            double sg = (m & 1) ? -1.0 : 1.0;
            fr[m] = re;
            fi[m] = im;
            fr[-m] = sg * re;
            fi[-m] = -sg * im;
        }
    }
}

static void sh_s2m(int p, double x, double y, double z, double q, double *M) {
    double R[SH_NS];
    sh_solid(p, 0, x, y, z, R);
    for (int n = 0; n <= p; n++) {
        double s = (n & 1) ? -q : q;
        M[SH_IRE(n, 0)] += s * R[SH_IRE(n, 0)];
        for (int m = 1; m <= n; m++) {
            M[SH_IRE(n, m)] += s * R[SH_IRE(n, m)];
            M[SH_IIM(n, m)] -= s * R[SH_IIM(n, m)];
        }
    }
}

static void sh_m2m(int p, double x, double y, double z, const double *M, double *Ms) {
    const int S = 2 * p + 1;
    double R[SH_NS];
    double Rr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Ri[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    double Mr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Mi[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    sh_solid(p, 0, x, y, z, R);
    sh_expand(p, S, p, R, Rr, Ri);
    sh_expand(p, S, p, M, Mr, Mi);
    for (int n = 0; n <= p; n++) {
        for (int m = 0; m <= n; m++) {
            double sr = 0.0, si = 0.0;
            for (int j = 0; j <= n; j++) {
                double sg = (j & 1) ? -1.0 : 1.0;
                int lo = SH_MAX(-j, m - (n - j)), hi = SH_MIN(j, m + n - j);
                const double *rr = Rr + j * S + p, *ri = Ri + j * S + p;
                const double *mr = Mr + (n - j) * S + p + m, *mi = Mi + (n - j) * S + p + m;
                double ar = 0.0, ai = 0.0;
                for (int k = lo; k <= hi; k++) {
                    /* conj(R_j^k) * M_{n-j}^{m-k} */
                    ar += rr[k] * mr[-k] + ri[k] * mi[-k];
                    ai += rr[k] * mi[-k] - ri[k] * mr[-k];
                }
                sr += sg * ar;
                si += sg * ai;
            }
            Ms[SH_IRE(n, m)] += sr;
            if (m > 0) Ms[SH_IIM(n, m)] += si;
        }
    }
}

static void sh_l2l(int p, double x, double y, double z, const double *L, double *Ls) {
    const int S = 2 * p + 1;
    double R[SH_NS];
    double Rr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Ri[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    double Lr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Li[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    sh_solid(p, 0, x, y, z, R);
    sh_expand(p, S, p, R, Rr, Ri);
    sh_expand(p, S, p, L, Lr, Li);
    for (int n = 0; n <= p; n++) {
        for (int m = 0; m <= n; m++) {
            double sr = 0.0, si = 0.0;
            for (int j = n; j <= p; j++) {
                int lo = SH_MAX(-j, m - (j - n)), hi = SH_MIN(j, m + j - n);
                const double *lr = Lr + j * S + p, *li = Li + j * S + p;
                const double *rr = Rr + (j - n) * S + p - m, *ri = Ri + (j - n) * S + p - m;
                for (int k = lo; k <= hi; k++) {
                    /* L_j^k * R_{j-n}^{k-m} */
                    sr += lr[k] * rr[k] - li[k] * ri[k];
                    si += lr[k] * ri[k] + li[k] * rr[k];
                }
            }
            Ls[SH_IRE(n, m)] += sr;
            if (m > 0) Ls[SH_IIM(n, m)] += si;
        }
    }
}

/* O(p^4) M2L. total = 0: all j, n <= p (I needed to degree 2p); total = 1: j + n <= p. */
static void sh_m2l_plain(int p, int total, double x, double y, double z, const double *M, double *L) {
    const int D = total ? p : 2 * p;
    const int SI = 2 * D + 1, SM = 2 * p + 1;
    double I[(2 * SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    double Ir[(2 * SH_PMAX + 1) * (4 * SH_PMAX + 1)], Ii[(2 * SH_PMAX + 1) * (4 * SH_PMAX + 1)];
    double Mr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Mi[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    sh_solid(D, 1, x, y, z, I);
    sh_expand(D, SI, D, I, Ir, Ii);
    sh_expand(p, SM, p, M, Mr, Mi);
    for (int n = 0; n <= p; n++) {
        for (int m = 0; m <= n; m++) {
            double sr = 0.0, si = 0.0;
            int jmax = total ? p - n : p;
            for (int j = 0; j <= jmax; j++) {
                const double *mr = Mr + j * SM + p, *mi = Mi + j * SM + p;
                const double *ir = Ir + (j + n) * SI + D - m, *ii = Ii + (j + n) * SI + D - m;
                double ar = 0.0, ai = 0.0;
                for (int k = -j; k <= j; k++) {
                    ar += mr[k] * ir[k] - mi[k] * ii[k];
                    ai += mr[k] * ii[k] + mi[k] * ir[k];
                }
                sr += ar;
                si += ai;
            }
            double sg = ((n + m) & 1) ? -1.0 : 1.0;
            L[SH_IRE(n, m)] += sg * sr;
            if (m > 0) L[SH_IIM(n, m)] += sg * si;
        }
    }
}

/* Potential and field from a local (c = L, A = R, nn = n-1) or multipole
   (c = M, A = I, nn = n+1) expansion. isL = 1 for the local case. */
static void sh_eval(int p, int isL, const double *c, const double *A, double *F) {
    double phi = 0.0, fz = 0.0, t1 = 0.0, t2 = 0.0;
    for (int n = 0; n <= p; n++) {
        int nn = isL ? n - 1 : n + 1;
        phi += c[SH_IRE(n, 0)] * A[SH_IRE(n, 0)];
        for (int m = 1; m <= n; m++)
            phi += 2.0 * (c[SH_IRE(n, m)] * A[SH_IRE(n, m)] - c[SH_IIM(n, m)] * A[SH_IIM(n, m)]);
        if (nn < 0) continue;
        fz += c[SH_IRE(n, 0)] * A[SH_IRE(nn, 0)];
        for (int m = 1; m <= SH_MIN(n, nn); m++)
            fz += 2.0 * (c[SH_IRE(n, m)] * A[SH_IRE(nn, m)] - c[SH_IIM(n, m)] * A[SH_IIM(nn, m)]);
        /* T = sum_m c_n^m A_nn^{m-1} */
        for (int m = 1; m <= SH_MIN(n, nn + 1); m++) {
            double cr = c[SH_IRE(n, m)], ci = c[SH_IIM(n, m)];
            double ar = A[SH_IRE(nn, m - 1)], ai = (m > 1) ? A[SH_IIM(nn, m - 1)] : 0.0;
            t1 += cr * ar - ci * ai;
            t2 += cr * ai + ci * ar;
        }
        for (int a = 0; a <= SH_MIN(n, nn - 1); a++) {
            double cr = c[SH_IRE(n, a)], ci = (a > 0) ? c[SH_IIM(n, a)] : 0.0;
            double ar = A[SH_IRE(nn, a + 1)], ai = A[SH_IIM(nn, a + 1)];
            t1 -= cr * ar - ci * ai;
            t2 += cr * ai + ci * ar;
        }
    }
    F[0] += phi;
    F[1] += t1;
    F[2] -= t2;
    F[3] += isL ? -fz : fz;
}

static void sh_l2p(int p, double x, double y, double z, const double *L, double *F) {
    double R[SH_NS];
    sh_solid(p, 0, x, y, z, R);
    sh_eval(p, 1, L, R, F);
}

static void sh_m2p(int p, double x, double y, double z, const double *M, double *F) {
    double I[(SH_PMAX + 2) * (SH_PMAX + 2)];
    sh_solid(p + 1, 1, x, y, z, I);
    sh_eval(p, 0, M, I, F);
}

/* ---- rotation-based M2L ------------------------------------------------- */
/* Dt layout: Dt[0] = 1/rho, then blocks l = 0..p row-major (2l+1)x(2l+1) at 1 + SH_ROTOFF[l],
   real orthonormal spherical harmonics, S_l(Q r) = D_l S_l(r), Q d/|d| = z. */
static inline double sh_ir_P(const double *R1, int i, int l, int a, int b, const double *Rp) {
    double ri1 = R1[(i + 1) * 3 + 2], rim1 = R1[(i + 1) * 3 + 0], ri0 = R1[(i + 1) * 3 + 1];
    int w = 2 * l - 1;
    const double *row = Rp + (a + l - 1) * w;
    if (b == -l) return ri1 * row[0] + rim1 * row[w - 1];
    if (b == l) return ri1 * row[w - 1] - rim1 * row[0];
    return ri0 * row[b + l - 1];
}

static void sh_rot_setup(int p, double x, double y, double z, double *Dt) {
    double rho = sqrt(x * x + y * y + z * z);
    double rxy = sqrt(x * x + y * y);
    double ca = 1.0, sa = 0.0;
    if (rxy >= 1e-14 * rho) { ca = x / rxy; sa = y / rxy; }
    double cb = z / rho, sb = rxy / rho;
    double Q[9] = {cb * ca, cb * sa, -sb, -sa, ca, 0.0, sb * ca, sb * sa, cb};
    static const int perm[3] = {1, 2, 0}; /* (y, z, x) */
    double R1[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) R1[i * 3 + j] = Q[perm[i] * 3 + perm[j]];
    Dt[0] = 1.0 / rho;
    double *D = Dt + 1;
    D[0] = 1.0;
    if (p >= 1)
        for (int i = 0; i < 9; i++) D[SH_ROTOFF[1] + i] = R1[i];
    for (int l = 2; l <= p; l++) {
        const double *Rp = D + SH_ROTOFF[l - 1];
        double *Rl = D + SH_ROTOFF[l];
        const int W = 2 * l + 1;
        const double *tu = SH_IRU + SH_IROFF[l], *tv = SH_IRV + SH_IROFF[l], *tw = SH_IRW + SH_IROFF[l];
        for (int m = -l; m <= l; m++) {
            int am = m < 0 ? -m : m;
            for (int n = -l; n <= l; n++) {
                int an = n < 0 ? -n : n;
                int t = am * (l + 1) + an;
                double u = tu[t], v = tv[t], w = tw[t], val = 0.0;
                if (u != 0.0) val += u * sh_ir_P(R1, 0, l, m, n, Rp);
                if (v != 0.0) {
                    double V;
                    if (m == 0) V = sh_ir_P(R1, 1, l, 1, n, Rp) + sh_ir_P(R1, -1, l, -1, n, Rp);
                    else if (m > 0) {
                        if (m == 1) V = sh_ir_P(R1, 1, l, 0, n, Rp) * 1.4142135623730951;
                        else V = sh_ir_P(R1, 1, l, m - 1, n, Rp) - sh_ir_P(R1, -1, l, -m + 1, n, Rp);
                    } else {
                        if (m == -1) V = sh_ir_P(R1, -1, l, 0, n, Rp) * 1.4142135623730951;
                        else V = sh_ir_P(R1, 1, l, m + 1, n, Rp) + sh_ir_P(R1, -1, l, -m - 1, n, Rp);
                    }
                    val += v * V;
                }
                if (w != 0.0) {
                    double Wv;
                    if (m > 0) Wv = sh_ir_P(R1, 1, l, m + 1, n, Rp) + sh_ir_P(R1, -1, l, -m - 1, n, Rp);
                    else Wv = sh_ir_P(R1, 1, l, m - 1, n, Rp) - sh_ir_P(R1, -1, l, -m + 1, n, Rp);
                    val += w * Wv;
                }
                Rl[(m + l) * W + (n + l)] = val;
            }
        }
    }
}

static void sh_rot_apply_blocks(int p, const double *D, const double *in, double *out, int transpose) {
    out[0] = in[0];
    for (int n = 1; n <= p; n++) {
        const int W = 2 * n + 1;
        const double *B = D + SH_ROTOFF[n];
        const double *x = in + n * n;
        double *y = out + n * n;
        if (!transpose) {
            for (int a = 0; a < W; a++) {
                const double *row = B + a * W;
                double s = 0.0;
                for (int b = 0; b < W; b++) s += row[b] * x[b];
                y[a] = s;
            }
        } else {
            for (int a = 0; a < W; a++) y[a] = 0.0;
            for (int b = 0; b < W; b++) {
                const double *row = B + b * W;
                double xb = x[b];
                for (int a = 0; a < W; a++) y[a] += row[a] * xb;
            }
        }
    }
}

static void sh_m2l_rot_apply(int p, int total, const double *Dt, const double *M, double *L) {
    const double rinv = Dt[0];
    const double *D = Dt + 1;
    double c[SH_NS], cm[SH_NS], cl[SH_NS], rp[SH_PMAX + 2];
    const int P1 = SH_PMAX + 1;
    for (int i = 0; i < (p + 1) * (p + 1); i++) c[i] = SH_SB[i] * M[i];
    sh_rot_apply_blocks(p, D, c, cm, 0);
    rp[0] = 1.0;
    for (int k = 1; k <= p + 1; k++) rp[k] = rp[k - 1] * rinv;
    for (int j = 0; j <= p; j++)
        for (int i = j * j; i < (j + 1) * (j + 1); i++) cm[i] *= rp[j + 1];
    for (int n = 0; n <= p; n++) {
        int jmax = total ? p - n : p;
        for (int m = 0; m <= n; m++) {
            double sr = 0.0, si = 0.0;
            for (int j = m; j <= jmax; j++) {
                double z = SH_ZT[(n * P1 + j) * P1 + m];
                sr += z * cm[SH_IRE(j, m)];
                if (m > 0) si += z * cm[SH_IIM(j, m)];
            }
            cl[SH_IRE(n, m)] = sr * rp[n];
            if (m > 0) cl[SH_IIM(n, m)] = si * rp[n];
        }
    }
    sh_rot_apply_blocks(p, D, cl, c, 1);
    for (int i = 0; i < (p + 1) * (p + 1); i++) L[i] += SH_SAI[i] * c[i];
}
"""


_PLANAR_CODE = r"""
/* ---- planar variants: sources, targets and centres in the z = 0 plane --------------------
   The solid harmonics with n + |m| odd vanish at z = 0, and so do the multipole and local
   coefficients built from them, so an expansion keeps only the n + |m| even entries,
   n(n+1)/2 + (n+k)/2 for the signed index k (k > 0: real part of order k, k < 0: imaginary
   part of order |k|): (p+1)(p+2)/2 of them. The z field is identically zero. */
#define SHP_I(n, k) ((n) * ((n) + 1) / 2 + ((n) + (k)) / 2)
#define SHP_RE(n, m) SHP_I(n, m)
#define SHP_IM(n, m) SHP_I(n, -(m))
/* first k >= lo with k = j (mod 2) */
#define SHP_START(lo, j) ((lo) + ((((lo) - (j)) & 1) ? 1 : 0))

static void shp_scatter(int p, const double *Ac, double *A) {
    for (int i = 0; i < (p + 1) * (p + 1); i++) A[i] = 0.0;
    for (int n = 0; n <= p; n++)
        for (int k = -n; k <= n; k += 2) A[n * n + n + k] = Ac[SHP_I(n, k)];
}

static void shp_gather_add(int p, const double *A, double *Ac) {
    for (int n = 0; n <= p; n++)
        for (int k = -n; k <= n; k += 2) Ac[SHP_I(n, k)] += A[n * n + n + k];
}

static void sh_s2m_planar(int p, double x, double y, double q, double *Mc) {
    double M[SH_NS] = {0};
    sh_s2m(p, x, y, 0.0, q, M);
    shp_gather_add(p, M, Mc);
}

static void sh_l2p_planar(int p, double x, double y, const double *Lc, double *F) {
    double R[SH_NS], L[SH_NS];
    sh_solid(p, 0, x, y, 0.0, R);
    shp_scatter(p, Lc, L);
    sh_eval(p, 1, L, R, F);
}

static void sh_m2p_planar(int p, double x, double y, const double *Mc, double *F) {
    double I[(SH_PMAX + 2) * (SH_PMAX + 2)], M[SH_NS];
    sh_solid(p + 1, 1, x, y, 0.0, I);
    shp_scatter(p, Mc, M);
    sh_eval(p, 0, M, I, F);
}

/* The O(p^4) operators: only the n + m even outputs, and the inner sum over k skips the entries that
   are zero (k = j mod 2). */
static void sh_m2m_planar(int p, double x, double y, const double *Mc, double *Msc) {
    const int S = 2 * p + 1;
    double R[SH_NS], M[SH_NS];
    double Rr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Ri[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    double Mr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Mi[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    sh_solid(p, 0, x, y, 0.0, R);
    shp_scatter(p, Mc, M);
    sh_expand(p, S, p, R, Rr, Ri);
    sh_expand(p, S, p, M, Mr, Mi);
    for (int n = 0; n <= p; n++) {
        for (int m = n & 1; m <= n; m += 2) {
            double sr = 0.0, si = 0.0;
            for (int j = 0; j <= n; j++) {
                double sg = (j & 1) ? -1.0 : 1.0;
                int lo = SH_MAX(-j, m - (n - j)), hi = SH_MIN(j, m + n - j);
                const double *rr = Rr + j * S + p, *ri = Ri + j * S + p;
                const double *mr = Mr + (n - j) * S + p + m, *mi = Mi + (n - j) * S + p + m;
                double ar = 0.0, ai = 0.0;
                for (int k = SHP_START(lo, j); k <= hi; k += 2) {
                    ar += rr[k] * mr[-k] + ri[k] * mi[-k];
                    ai += rr[k] * mi[-k] - ri[k] * mr[-k];
                }
                sr += sg * ar;
                si += sg * ai;
            }
            Msc[SHP_RE(n, m)] += sr;
            if (m > 0) Msc[SHP_IM(n, m)] += si;
        }
    }
}

static void sh_l2l_planar(int p, double x, double y, const double *Lc, double *Lsc) {
    const int S = 2 * p + 1;
    double R[SH_NS], L[SH_NS];
    double Rr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Ri[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    double Lr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Li[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    sh_solid(p, 0, x, y, 0.0, R);
    shp_scatter(p, Lc, L);
    sh_expand(p, S, p, R, Rr, Ri);
    sh_expand(p, S, p, L, Lr, Li);
    for (int n = 0; n <= p; n++) {
        for (int m = n & 1; m <= n; m += 2) {
            double sr = 0.0, si = 0.0;
            for (int j = n; j <= p; j++) {
                int lo = SH_MAX(-j, m - (j - n)), hi = SH_MIN(j, m + j - n);
                const double *lr = Lr + j * S + p, *li = Li + j * S + p;
                const double *rr = Rr + (j - n) * S + p - m, *ri = Ri + (j - n) * S + p - m;
                for (int k = SHP_START(lo, j); k <= hi; k += 2) {
                    sr += lr[k] * rr[k] - li[k] * ri[k];
                    si += lr[k] * ri[k] + li[k] * rr[k];
                }
            }
            Lsc[SHP_RE(n, m)] += sr;
            if (m > 0) Lsc[SHP_IM(n, m)] += si;
        }
    }
}

static void sh_m2l_plain_planar(int p, int total, double x, double y, const double *Mc, double *Lc) {
    const int D = total ? p : 2 * p;
    const int SI = 2 * D + 1, SM = 2 * p + 1;
    double I[(2 * SH_PMAX + 1) * (2 * SH_PMAX + 1)], M[SH_NS];
    double Ir[(2 * SH_PMAX + 1) * (4 * SH_PMAX + 1)], Ii[(2 * SH_PMAX + 1) * (4 * SH_PMAX + 1)];
    double Mr[(SH_PMAX + 1) * (2 * SH_PMAX + 1)], Mi[(SH_PMAX + 1) * (2 * SH_PMAX + 1)];
    sh_solid(D, 1, x, y, 0.0, I);
    shp_scatter(p, Mc, M);
    sh_expand(D, SI, D, I, Ir, Ii);
    sh_expand(p, SM, p, M, Mr, Mi);
    for (int n = 0; n <= p; n++) {
        for (int m = n & 1; m <= n; m += 2) {
            double sr = 0.0, si = 0.0;
            int jmax = total ? p - n : p;
            for (int j = 0; j <= jmax; j++) {
                const double *mr = Mr + j * SM + p, *mi = Mi + j * SM + p;
                const double *ir = Ir + (j + n) * SI + D - m, *ii = Ii + (j + n) * SI + D - m;
                double ar = 0.0, ai = 0.0;
                for (int k = -j; k <= j; k += 2) {
                    ar += mr[k] * ir[k] - mi[k] * ii[k];
                    ai += mr[k] * ii[k] + mi[k] * ir[k];
                }
                sr += ar;
                si += ai;
            }
            Lc[SHP_RE(n, m)] += sr;
            if (m > 0) Lc[SHP_IM(n, m)] += si;
        }
    }
}
"""

_P2P_PLANAR = """double Rinv = 1.0 / sqrt(x*x + y*y);
    double Ftmp0 = Rinv*Rinv*Rinv*S[0];
    F[0] += Rinv*S[0];
    F[1] += Ftmp0*x;
    F[2] += Ftmp0*y;"""

_P2P_BATCH_PLANAR = """void P2P_batchxy(double tx, double ty, double tz,
    const double * FMMGEN_RESTRICT sx,
    const double * FMMGEN_RESTRICT sy,
    const double * FMMGEN_RESTRICT sz,
    const double * FMMGEN_RESTRICT S,
    size_t begin, size_t end, double * FMMGEN_RESTRICT F) {
    double facc0 = 0.0, facc1 = 0.0, facc2 = 0.0;
    (void)tz; (void)sz;
#pragma omp simd reduction(+:facc0,facc1,facc2)
    for (size_t u = begin; u < end; u++) {
        double x = tx - sx[u];
        double y = ty - sy[u];
        double Rinv = 1.0 / sqrt(x*x + y*y);
        double Ftmp0 = Rinv*Rinv*Rinv*S[u];
        facc0 += Rinv*S[u];
        facc1 += Ftmp0*x;
        facc2 += Ftmp0*y;
    }
    F[0] += facc0;
    F[1] += facc1;
    F[2] += facc2;
}
"""


def _split_code():
    code = _LOOP_CODE
    i1 = code.index("/* Solid harmonics")
    i2 = code.index("/* ---- rotation-based M2L")
    i3 = code.index("static void sh_rot_apply_blocks")
    macros = code[:i1]
    body = code[i1:i2]
    rot_setup = code[i2:i3]
    rot_apply = code[i3:]
    return macros, body, rot_setup, rot_apply


# ---------------------------------------------------------------------------
# Unrolled emitters
# ---------------------------------------------------------------------------
def _emit_solid(name, p, irr):
    """Straight-line statements filling `name` (flat, degree p) with R or I."""
    ire, iim = sm.ire, sm.iim
    L = [f"double {name}[{(p + 1) ** 2}];"]
    L.append("const double r2 = x*x + y*y + z*z;")
    if irr:
        L.append("const double ri2 = 1.0/r2;")
        L.append(f"{name}[0] = sqrt(ri2);")
    else:
        L.append(f"{name}[0] = 1.0;")
    A = name
    for m in range(p + 1):
        if m > 0:
            f = f"{-(2 * m - 1):.1f}*ri2" if irr else _num(-1.0 / (2 * m))
            pre = f"{A}[{ire(m - 1, m - 1)}]"
            if m > 1:
                pim = f"{A}[{iim(m - 1, m - 1)}]"
                L.append(f"{A}[{ire(m, m)}] = {f}*(x*{pre} - y*{pim});")
                L.append(f"{A}[{iim(m, m)}] = {f}*(x*{pim} + y*{pre});")
            else:
                L.append(f"{A}[{ire(m, m)}] = {f}*(x*{pre});")
                L.append(f"{A}[{iim(m, m)}] = {f}*(y*{pre});")
        if m + 1 <= p:
            f = f"{2 * m + 1:.1f}*z*ri2" if irr else "z"
            L.append(f"{A}[{ire(m + 1, m)}] = {f}*{A}[{ire(m, m)}];")
            if m > 0:
                L.append(f"{A}[{iim(m + 1, m)}] = {f}*{A}[{iim(m, m)}];")
        for n in range(m + 2, p + 1):
            for slot in ([ire] + ([iim] if m > 0 else [])):
                a = f"{A}[{slot(n - 1, m)}]"
                b = f"{A}[{slot(n - 2, m)}]"
                if irr:
                    L.append(f"{A}[{slot(n, m)}] = ({2 * n - 1:.1f}*z*{a} - {(n - 1) ** 2 - m * m:.1f}*{b})*ri2;")
                else:
                    L.append(f"{A}[{slot(n, m)}] = ({2 * n - 1:.1f}*z*{a} - r2*{b})*{_num(1.0 / (n * n - m * m))};")
    return L


def _term_str(c, a, b, An, Bn):
    body = f"{An}[{a}]*{Bn}[{b}]"
    if c == 1.0:
        return "+", body
    if c == -1.0:
        return "-", body
    s = "+" if c > 0 else "-"
    return s, f"{_num(abs(c))}*{body}"


def _emit_terms(P, outn, An, Bn, chunk=8):
    o, c, a, b = P
    groups = {}
    for oo, cc, aa, bb in zip(o.tolist(), c.tolist(), a.tolist(), b.tolist()):
        groups.setdefault(oo, []).append((cc, aa, bb))
    L = []
    for oo in sorted(groups):
        ts = groups[oo]
        for i in range(0, len(ts), chunk):
            parts = [_term_str(cc, aa, bb, An, Bn) for cc, aa, bb in ts[i:i + chunk]]
            s = ""
            for k, (sg, body) in enumerate(parts):
                s += (("-" if sg == "-" else "") if k == 0 else f" {sg} ") + body
            L.append(f"{outn}[{oo}] += {s};")
    return L


def _emit_rot_apply(p, trunc):
    ire, iim = sm.ire, sm.iim
    N = (p + 1) ** 2
    sb, sa = sm.basis_scales(p)
    Z = sm.zshift_table(p)
    L = ["const double rinv = Dt[0];", "const double *D = Dt + 1;",
         f"double c[{N}], cm[{N}], cl[{N}], rp[{p + 2}];"]
    L.append("rp[0] = 1.0;")
    for k in range(1, p + 2):
        L.append(f"rp[{k}] = rp[{k - 1}]*rinv;")
    for i in range(N):
        L.append(f"c[{i}] = {_num(sb[i])}*M[{i}];")
    L.append("cm[0] = c[0];")
    roff = [0]
    for n in range(p + 1):
        roff.append(roff[-1] + (2 * n + 1) ** 2)
    for n in range(1, p + 1):
        W = 2 * n + 1
        for a in range(W):
            s = " + ".join(f"D[{roff[n] + a * W + b}]*c[{n * n + b}]" for b in range(W))
            L.append(f"cm[{n * n + a}] = {s};")
    for j in range(p + 1):
        for i in range(j * j, (j + 1) ** 2):
            L.append(f"cm[{i}] *= rp[{j + 1}];")
    for n in range(p + 1):
        jmax = p - n if trunc == "total" else p
        for m in range(n + 1):
            for slot in ([ire] + ([iim] if m > 0 else [])):
                s = " + ".join(f"{_num(Z[n, j, m])}*cm[{slot(j, m)}]" for j in range(m, jmax + 1)) or "0.0"
                if n > 0:
                    L.append(f"cl[{slot(n, m)}] = rp[{n}]*({s});")
                else:
                    L.append(f"cl[{slot(n, m)}] = {s};")
    L.append("c[0] = cl[0];")
    for n in range(1, p + 1):
        W = 2 * n + 1
        for a in range(W):
            s = " + ".join(f"D[{roff[n] + b * W + a}]*cl[{n * n + b}]" for b in range(W))
            L.append(f"c[{n * n + a}] = {s};")
    for i in range(N):
        L.append(f"L[{i}] += {_num(1.0 / sa[i])}*c[{i}];")
    return L


_SIG = "double x, double y, double z, double * FMMGEN_RESTRICT {a}, double * FMMGEN_RESTRICT {b}"


def _fn(name, a, b, body):
    return f"void {name}(" + _SIG.format(a=a, b=b) + ") {\n    " + "\n    ".join(body) + "\n}\n"


def _decl(name, a, b):
    return f"void {name}(" + _SIG.format(a=a, b=b) + ");"


_ARGS = {"S2M": ("S", "M"), "M2M": ("M", "Ms"), "L2L": ("L", "Ls"), "L2P": ("L", "F"),
         "M2P": ("M", "F"), "M2Lplain": ("M", "L")}


def unrolled_function_bodies(p, trunc):
    """Return {funcname: body lines} for the unrolled variants at order p."""
    out = {}
    out[f"S2M_{p}"] = _emit_solid("R", p, False) + _emit_terms(sm.terms_S2M(p), "M", "R", "S")
    out[f"M2M_{p}"] = _emit_solid("R", p, False) + _emit_terms(sm.terms_M2M(p), "Ms", "R", "M")
    out[f"L2L_{p}"] = _emit_solid("R", p, False) + _emit_terms(sm.terms_L2L(p), "Ls", "R", "L")
    D = sm.degree_I_M2L(p, trunc)
    out[f"M2Lplain_{p}"] = _emit_solid("I", D, True) + _emit_terms(sm.terms_M2L(p, trunc), "L", "I", "M")
    out[f"L2P_{p}"] = _emit_solid("R", p, False) + _emit_terms(sm.terms_eval(p, "L2P"), "F", "R", "L")
    out[f"M2P_{p}"] = _emit_solid("I", p + 1, True) + _emit_terms(sm.terms_eval(p, "M2P"), "F", "I", "M")
    return out


_P2P_BATCH = """void P2P_batch(double tx, double ty, double tz,
    const double * FMMGEN_RESTRICT sx,
    const double * FMMGEN_RESTRICT sy,
    const double * FMMGEN_RESTRICT sz,
    const double * FMMGEN_RESTRICT S,
    size_t begin, size_t end, double * FMMGEN_RESTRICT F) {
    double facc0 = 0.0, facc1 = 0.0, facc2 = 0.0, facc3 = 0.0;
#pragma omp simd reduction(+:facc0,facc1,facc2,facc3)
    for (size_t u = begin; u < end; u++) {
        double x = tx - sx[u];
        double y = ty - sy[u];
        double z = tz - sz[u];
        double Rinv = 1.0 / sqrt(x*x + y*y + z*z);
        double Ftmp0 = Rinv*Rinv*Rinv*S[u];
        facc0 += Rinv*S[u];
        facc1 += Ftmp0*x;
        facc2 += Ftmp0*y;
        facc3 += Ftmp0*z;
    }
    F[0] += facc0;
    F[1] += facc1;
    F[2] += facc2;
    F[3] += facc3;
}
"""


def generate_sh_code(order, name, unroll=False, m2l="plain", trunc="box", language="c",
                     include_dir=None, src_dir=None, pmin=1, max_unroll=DEFAULT_MAX_UNROLL, orders=None,
                     prefix="", planar=False):
    """Write <name>.h and <name>.c (or .cpp) with the SH operators for p = pmin..order-1
    (the same convention as fmmgen.generate_code: FMMGEN_MAXORDER = order is exclusive).

    orders: an explicit list of degrees to emit instead of pmin..order-1 (they need not be contiguous;
    `order` must then exceed the largest). prefix: prepended to every function and to the header macros
    (FMMGEN_<PREFIX>MAXORDER ...), so the operators can be linked into one program with another
    backend's.

    planar: also emit S2Mxy_p, M2Mxy_p, M2Lxy_p, L2Lxy_p, L2Pxy_p, M2Pxy_p, P2Pxy and P2P_batchxy for charges
    and targets confined to the z = 0 plane (every displacement then has z = 0, and z is ignored). An expansion
    keeps only its n + |m| even entries, (p+1)(p+2)/2 of them (see `sh_math.planar_keep`), the O(p^4) operators
    do about a quarter of the work, and the z field is exactly zero. They are loops whatever `unroll` is, and
    M2Lxy_p is the plain double sum: the rotation-based M2L has no planar form, because rotating an in-plane
    displacement onto the z axis destroys the parity structure.

    Returns a report dict (paths, orders, functions, generated size)."""
    if language not in ("c", "c++"):
        raise ValueError("language must be 'c' or 'c++'")
    if m2l not in ("rotation", "plain"):
        raise ValueError("m2l must be 'rotation' or 'plain'")
    if trunc not in ("box", "total"):
        raise ValueError("trunc must be 'box' or 'total'")
    if prefix and not prefix.replace("_", "a").isalnum():
        raise ValueError("prefix must be a valid C identifier fragment")
    if orders is not None:
        orders = sorted(set(orders))
        if not orders or orders[0] < 1:
            raise ValueError("orders must be a non-empty list of degrees >= 1")
        pmin = orders[0]
        if order <= orders[-1]:
            raise ValueError("order must exceed the largest of orders (FMMGEN_MAXORDER is exclusive)")
    pmax = order - 1 if orders is None else orders[-1]
    if pmax < pmin:
        raise ValueError("order must exceed pmin (orders pmin..order-1 are generated)")
    if unroll and pmax > max_unroll:
        raise ValueError(
            f"unroll=True with p up to {pmax} > max_unroll={max_unroll}: straight-line code grows like p^5 "
            "(see fmmgen/sh.py); use unroll=False or raise max_unroll explicitly")
    include_dir = include_dir or "."
    src_dir = src_dir or include_dir
    ext = "c" if language == "c" else "cpp"
    orders = list(range(pmin, pmax + 1)) if orders is None else orders
    pf = prefix
    total = 1 if trunc == "total" else 0

    macros, body, rot_setup, rot_apply = _split_code()
    src = [f'#include "{name}.h"', "#include <math.h>", "#include <stddef.h>",
           f"#define SH_PMAX {pmax}", _tables(pmax, full=not unroll), macros]
    if not unroll or planar:
        src.append(body)
    if planar:
        src.append(_PLANAR_CODE)
    src.append(rot_setup)
    if not unroll:
        src.append(rot_apply)

    sizes_list = ", ".join(str((p + 1) ** 2) for p in range(order + 1))
    m = "FMMGEN_" + (prefix.upper().rstrip("_") + "_" if prefix else "")
    hdr = ["#pragma once", "#include <stddef.h>", _RESTRICT.rstrip(),
           f"#define {m}MINORDER {pmin}", f"#define {m}MAXORDER {order}",
           f"#define {m}SOURCEORDER 0", f"#define {m}SOURCESIZE 1", f"#define {m}OUTPUTSIZE 4",
           f"#define {m}SPHERICAL 1",
           f"/* (p+1)^2 reals per expansion; index by order p. {m}ROTSIZE: doubles per rotation table. */",
           f"static const size_t {m}MULTIPOLESIZE[] = {{{sizes_list}}};",
           f"static const size_t {m}LOCALSIZE[] = {{{sizes_list}}};",
           f"static const size_t {m}ROTSIZE[] = {{" + ", ".join(str(rot_size(p)) for p in range(order + 1)) + "};"]
    if planar:
        psizes = ", ".join(str(sm.planar_size(p)) for p in range(order + 1))
        hdr += [f"#define {m}PLANAR 1",
                "/* planar arrays: (p+1)(p+2)/2 reals, the n + |m| even entries; index by order p. */",
                f"static const size_t {m}PLANAR_MULTIPOLESIZE[] = {{{psizes}}};",
                f"static const size_t {m}PLANAR_LOCALSIZE[] = {{{psizes}}};"]

    names = [("S2M", "S", "M"), ("M2M", "M", "Ms"), ("M2L", "M", "L"), ("L2L", "L", "Ls"),
             ("L2P", "L", "F"), ("M2P", "M", "F")]
    extra = [("M2Lplain", "M", "L"), ("M2Lrot", "M", "L")]
    pnames = [("S2Mxy", "S", "M"), ("M2Mxy", "M", "Ms"), ("M2Lxy", "M", "L"), ("L2Lxy", "L", "Ls"),
              ("L2Pxy", "L", "F"), ("M2Pxy", "M", "F")] if planar else []
    for p in orders:
        for nm, a, b in names + extra + pnames:
            hdr.append(_decl(f"{pf}{nm}_{p}", a, b))
        hdr.append(f"void {pf}M2Lrot_setup_{p}(double x, double y, double z, double * FMMGEN_RESTRICT Dt);")
        hdr.append(f"void {pf}M2Lrot_apply_{p}(const double * FMMGEN_RESTRICT Dt, "
                   "const double * FMMGEN_RESTRICT M, double * FMMGEN_RESTRICT L);")
    hdr.append(_decl(f"{pf}P2P", "S", "F"))
    hdr.append(f"void {pf}P2P_batch(double tx, double ty, double tz, const double * FMMGEN_RESTRICT sx, "
               "const double * FMMGEN_RESTRICT sy, const double * FMMGEN_RESTRICT sz, "
               "const double * FMMGEN_RESTRICT S, size_t begin, size_t end, double * FMMGEN_RESTRICT F);")
    if planar:
        hdr.append(_decl(f"{pf}P2Pxy", "S", "F"))
        hdr.append(f"void {pf}P2P_batchxy(double tx, double ty, double tz, const double * FMMGEN_RESTRICT sx, "
                   "const double * FMMGEN_RESTRICT sy, const double * FMMGEN_RESTRICT sz, "
                   "const double * FMMGEN_RESTRICT S, size_t begin, size_t end, double * FMMGEN_RESTRICT F);")
    for nm, a, b in names + extra + pnames:
        hdr.append(f"void {pf}{nm}(double x, double y, double z, double * FMMGEN_RESTRICT {a}, "
                   f"double * FMMGEN_RESTRICT {b}, int order);")

    setup_fn = ("void {P}M2Lrot_setup_{p}(double x, double y, double z, double * FMMGEN_RESTRICT Dt) "
                "{{\n    sh_rot_setup({p}, x, y, z, Dt);\n}}\n")
    apply_sig = ("void {P}M2Lrot_apply_{p}(const double * FMMGEN_RESTRICT Dt, const double * FMMGEN_RESTRICT M, "
                 "double * FMMGEN_RESTRICT L) {{\n    ")
    for p in orders:
        if unroll:
            for fname, lines in unrolled_function_bodies(p, trunc).items():
                a, b = _ARGS[fname.rsplit("_", 1)[0]]
                src.append(_fn(pf + fname, a, b, lines))
            src.append(setup_fn.format(p=p, P=pf))
            src.append(apply_sig.format(p=p, P=pf) + "\n    ".join(_emit_rot_apply(p, trunc)) + "\n}\n")
        else:
            src.append(_fn(f"{pf}S2M_{p}", "S", "M", [f"sh_s2m({p}, x, y, z, S[0], M);"]))
            src.append(_fn(f"{pf}M2M_{p}", "M", "Ms", [f"sh_m2m({p}, x, y, z, M, Ms);"]))
            src.append(_fn(f"{pf}L2L_{p}", "L", "Ls", [f"sh_l2l({p}, x, y, z, L, Ls);"]))
            src.append(_fn(f"{pf}L2P_{p}", "L", "F", [f"sh_l2p({p}, x, y, z, L, F);"]))
            src.append(_fn(f"{pf}M2P_{p}", "M", "F", [f"sh_m2p({p}, x, y, z, M, F);"]))
            src.append(_fn(f"{pf}M2Lplain_{p}", "M", "L", [f"sh_m2l_plain({p}, {total}, x, y, z, M, L);"]))
            src.append(setup_fn.format(p=p, P=pf))
            src.append(apply_sig.format(p=p, P=pf) + f"sh_m2l_rot_apply({p}, {total}, Dt, M, L);\n}}\n")
        src.append(_fn(f"{pf}M2Lrot_{p}", "M", "L",
                       [f"double Dt[{rot_size(p)}];", f"{pf}M2Lrot_setup_{p}(x, y, z, Dt);",
                        f"{pf}M2Lrot_apply_{p}(Dt, M, L);"]))
        alias = "M2Lrot" if m2l == "rotation" else "M2Lplain"
        src.append(_fn(f"{pf}M2L_{p}", "M", "L", [f"{pf}{alias}_{p}(x, y, z, M, L);"]))
        if planar:
            for nm, a, b, call in (("S2Mxy", "S", "M", f"sh_s2m_planar({p}, x, y, S[0], M);"),
                                   ("M2Mxy", "M", "Ms", f"sh_m2m_planar({p}, x, y, M, Ms);"),
                                   ("M2Lxy", "M", "L", f"sh_m2l_plain_planar({p}, {total}, x, y, M, L);"),
                                   ("L2Lxy", "L", "Ls", f"sh_l2l_planar({p}, x, y, L, Ls);"),
                                   ("L2Pxy", "L", "F", f"sh_l2p_planar({p}, x, y, L, F);"),
                                   ("M2Pxy", "M", "F", f"sh_m2p_planar({p}, x, y, M, F);")):
                src.append(_fn(f"{pf}{nm}_{p}", a, b, ["(void)z;", call]))

    src.append(_fn(f"{pf}P2P", "S", "F", [
        "double Rinv = 1.0 / sqrt(x*x + y*y + z*z);",
        "double Ftmp0 = Rinv*Rinv*Rinv*S[0];",
        "F[0] += Rinv*S[0];", "F[1] += Ftmp0*x;", "F[2] += Ftmp0*y;", "F[3] += Ftmp0*z;"]))
    src.append(_P2P_BATCH.replace("void P2P_batch(", f"void {pf}P2P_batch(", 1))
    if planar:
        src.append(_fn(f"{pf}P2Pxy", "S", "F", ["(void)z;"] + _P2P_PLANAR.split("\n    ")))
        src.append(_P2P_BATCH_PLANAR.replace("void P2P_batchxy(", f"void {pf}P2P_batchxy(", 1))
    for nm, a, b in names + extra + pnames:
        cases = "".join(f"        case {p}:\n            {pf}{nm}_{p}(x, y, z, {a}, {b});\n            break;\n"
                        for p in orders)
        src.append(f"void {pf}{nm}(double x, double y, double z, double * FMMGEN_RESTRICT {a}, "
                   f"double * FMMGEN_RESTRICT {b}, int order) {{\n    switch (order) {{\n{cases}    }}\n}}\n")

    os.makedirs(include_dir, exist_ok=True)
    os.makedirs(src_dir, exist_ok=True)
    hpath = os.path.join(include_dir, f"{name}.h")
    cpath = os.path.join(src_dir, f"{name}.{ext}")
    with open(hpath, "w") as f:
        f.write("\n".join(hdr) + "\n")
    text = "\n".join(src)
    with open(cpath, "w") as f:
        f.write(text)
    functions = [f"{pf}{nm}_{p}" for p in orders for nm in ("S2M", "M2M", "M2L", "L2L", "L2P", "M2P", "M2Lplain",
                                                              "M2Lrot", "M2Lrot_setup", "M2Lrot_apply")]
    if planar:
        functions += [f"{pf}{nm}_{p}" for p in orders for nm in ("S2Mxy", "M2Mxy", "M2Lxy", "L2Lxy", "L2Pxy", "M2Pxy")]
        functions += [f"{pf}P2Pxy", f"{pf}P2P_batchxy"]
    functions += [f"{pf}{nm}" for nm in ("P2P", "P2P_batch", "S2M", "M2M", "M2L", "L2L", "L2P", "M2P", "M2Lplain",
                                         "M2Lrot")]
    if planar:
        functions += [f"{pf}{nm}xy" for nm in ("S2M", "M2M", "M2L", "L2L", "L2P", "M2P")]
    return {"header": hpath, "source": cpath, "orders": orders, "unroll": unroll, "functions": functions,
            "lines": text.count("\n"), "bytes": len(text)}
