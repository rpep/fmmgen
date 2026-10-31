"""Pluggable backends: how a kernel is turned into FMM operators.

A `Kernel` (or `VectorKernel`) says *what* is summed; a `Backend` decides *how*
(Cartesian Taylor, spherical harmonics, kernel-independent, ...). The tree code
needs only the manifest that `generate` writes next to the sources: function
names, array sizes per order, channels, dependencies.

    from fmmgen import Kernel, generate
    man = generate(Kernel(...), 8, "ops", backend="cartesian", language="c")
    fmmgen.recommend(Problem(kernel), [4, 8, 12])      # advice per order; you decide

Operator tiers
--------------
* "operator": per-displacement calls (`S2M_p(x, y, z, S, M)`, `M2L_p(x, y, z, M, L)`, ...),
  as the Cartesian output always had.
* "batch": `M2L_batch`-style calls over many pairs, for backends that amortise a
  transform across offsets (FFT or rotation based). Declared in the manifest;
  no backend implements it yet, but the contract is reserved so adding one does
  not change this module.
"""
from dataclasses import dataclass, field
import json
import os

__all__ = ["Problem", "Layout", "Manifest", "Backend", "CartesianBackend", "SphericalBackend", "register_backend",
           "backends", "get_backend", "recommend", "generate"]


@dataclass(frozen=True)
class Problem:
    """What to generate: a scalar kernel with sources of one order and `nchannels` channels."""

    kernel: object
    source_order: int = 0
    nchannels: int = 1
    potential: bool = True
    field: bool = True


@dataclass(frozen=True)
class Layout:
    """Array lengths (in doubles) at expansion order p, so tree code never assumes a basis."""

    multipole: int
    local: int
    source: int  # coefficients of one source (all channels)
    output: int  # entries of F for one target (all channels)


@dataclass
class Manifest:
    backend: str
    kernel: str
    language: str
    orders: list
    source_order: int
    nchannels: int
    tiers: list
    files: list
    functions: list
    layouts: dict
    dependencies: list = field(default_factory=list)
    groups: list = field(default_factory=list)  # one sub-manifest per scalar group of a vector kernel
    prefix: str = ""  # prepended to every function name
    by_order: dict = field(default_factory=dict)  # order -> backend name (parts of a mixed build)
    parts: list = field(default_factory=list)  # one sub-manifest per backend of a mixed build
    options: dict = field(default_factory=dict)  # backend-specific facts a driver may need (see the backend)

    def to_dict(self):
        return {k: v for k, v in self.__dict__.items()}

    def write(self, path):
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2, sort_keys=True)


class Backend:
    """Protocol. Subclasses set `name`, `tiers` and implement the four methods."""

    name = "backend"
    tiers = ("operator",)
    dependencies = ()
    supports_prefix = False  # can emit functions under a prefix, so it can coexist with other backends

    def supports(self, problem):
        """(bool, reason)."""
        raise NotImplementedError

    def cost(self, problem, order):
        """Relative cost estimate used by `recommend` (lower is better), or None if unknown."""
        return None

    def layout(self, problem, p):
        raise NotImplementedError

    def generate(self, problem, order, name, orders=None, language="c", **options):
        """Write the sources; return a Manifest."""
        raise NotImplementedError


class CartesianBackend(Backend):
    """Taylor/Cartesian operators with the derivative-tensor recurrence (`generate_kernel_code`)."""

    name = "cartesian"
    supports_prefix = True

    def supports(self, problem):
        k = problem.kernel
        if not 1 <= k.dim <= 3:
            return False, "dimension must be 1, 2 or 3"
        return True, ""

    def layout(self, problem, p):
        from .kernel import nterms, shell_indices, source_indices

        k, s, nch = problem.kernel, problem.source_order, problem.nchannels
        if p < s:
            raise ValueError("order must be >= source_order")
        out = (int(problem.potential) + (k.dim if problem.field else 0)) * nch
        return Layout(len(source_indices(k.dim, p, s)) * nch, nterms(k.dim, p - s) * nch,
                      len(shell_indices(k.dim, s)) * nch, out)

    def generate(self, problem, order, name, orders=None, language="c", prefix="", **options):
        from .kernel import generate_kernel_code

        orders = list(orders) if orders is not None else list(range(problem.source_order + 1, order + 1))
        generate_kernel_code(problem.kernel, order, name, orders=orders, language=language,
                             source_order=problem.source_order, nchannels=problem.nchannels,
                             potential=problem.potential, field=problem.field, prefix=prefix, **options)
        ext = {"c": ("h", "c"), "c++": ("h", "cpp"), "fortran": (None, "f90")}[language]
        files = [f"{name}.{e}" for e in ext if e]
        ops = ["S2M", "M2M", "M2L", "L2L", "L2P", "M2P"]
        functions = [f"{prefix}{o}_{p}" for p in orders for o in ops] + [f"{prefix}P2P", f"{prefix}P2P_batch"]
        return Manifest(self.name, problem.kernel.name, language, orders, problem.source_order, problem.nchannels,
                        list(self.tiers), files, functions,
                        {p: self.layout(problem, p).__dict__ for p in orders}, list(self.dependencies), prefix=prefix)


def _is_laplace(kernel):
    """True if the scalar kernel is exactly 1/|r| in three dimensions (checked numerically)."""
    if kernel.dim != 3:
        return False
    for pt in ((0.7, -0.4, 0.9), (-1.3, 0.2, 2.1)):
        value = float(kernel.expr.subs(dict(zip(kernel.variables, pt))))
        if abs(value * (sum(c * c for c in pt) ** 0.5) - 1.0) > 1e-12:
            return False
    return True


class SphericalBackend(Backend):
    """Spherical-harmonic operators for the Laplace kernel 1/r (see `fmmgen.sh`).

    (p+1)^2 real coefficients per expansion instead of the Cartesian (p+1)(p+2)(p+3)/6, an O(p^4) M2L
    or an O(p^3) rotation-based one (`m2l="rotation"`, or the cached `M2Lrot_setup_p`/`M2Lrot_apply_p`
    pair), and accuracy that keeps improving to 1e-8 and beyond, which the Cartesian operators reach only
    at orders that take hours to generate. At moderate accuracy (about 3e-4 and looser) the Laplace-tuned
    Cartesian operators are faster, so which backend is better depends on the order: see `recommend`.

    Options: `m2l` ("plain" or "rotation"), `trunc` ("box" or "total"), `unroll`, `max_unroll`; see
    `fmmgen.sh.generate_sh_code`. Cartesian-only options (CSE, minpow, horner) are ignored.
    """

    name = "spherical"
    supports_prefix = True
    dependencies = ("libm",)
    _IGNORED = ("CSE", "minpow", "horner", "report")

    def supports(self, problem):
        if not _is_laplace(problem.kernel):
            return False, "only the three-dimensional Laplace kernel 1/r"
        if problem.source_order != 0 or problem.nchannels != 1:
            return False, "charges only (source_order 0, one channel)"
        if not (problem.potential and problem.field):
            return False, "always computes the potential and the field"
        return True, ""

    def layout(self, problem, p):
        if p < 1:
            raise ValueError("order must be >= 1")
        return Layout((p + 1) ** 2, (p + 1) ** 2, 1, 4)

    def generate(self, problem, order, name, orders=None, language="c", prefix="", **options):
        from .sh import generate_sh_code, rot_size

        ok, why = self.supports(problem)
        if not ok:
            raise ValueError(f"spherical backend: {why}")
        orders = sorted(orders) if orders is not None else list(range(1, order + 1))
        opts = {k: v for k, v in options.items() if k not in self._IGNORED}
        report = generate_sh_code(max(orders) + 1, name, orders=orders, language=language, prefix=prefix, **opts)
        files = [os.path.basename(report["header"]), os.path.basename(report["source"])]
        info = {"m2l": opts.get("m2l", "plain"), "trunc": opts.get("trunc", "box"), "unroll": report["unroll"],
                "rotation_table_size": {p: rot_size(p) for p in orders}}
        return Manifest(self.name, problem.kernel.name, language, orders, 0, 1, list(self.tiers), files,
                        report["functions"], {p: self.layout(problem, p).__dict__ for p in orders},
                        list(self.dependencies), prefix=prefix, options=info)


_REGISTRY = {}


def register_backend(backend):
    _REGISTRY[backend.name] = backend
    return backend


def backends():
    return dict(_REGISTRY)


def get_backend(name):
    try:
        return _REGISTRY[name]
    except KeyError:
        raise ValueError(f"unknown backend {name!r}; registered: {sorted(_REGISTRY)}") from None


register_backend(CartesianBackend())
register_backend(SphericalBackend())


def _resolve(problem, backend, orders):
    """{order: Backend} from backend = name | {order or "default": name}; unsupported pairs are errors."""
    spec = backend if isinstance(backend, dict) else {"default": backend}
    unknown = [k for k in spec if k != "default" and k not in orders]
    if unknown:
        raise ValueError(f"backend given for orders {unknown} that are not being generated {list(orders)}")
    out = {}
    for p in orders:
        name = spec.get(p, spec.get("default"))
        if name is None:
            raise ValueError(f"no backend for order {p}: give one or a 'default'")
        b = get_backend(name)
        ok, why = b.supports(problem)
        if not ok:
            raise ValueError(f"backend {name!r} does not support this problem at order {p}: {why}")
        out[p] = b
    return out


def recommend(problem, orders):
    """Advice, never a silent choice: for each order, the supporting backends ranked by their cost
    estimate (cheapest first, unknown costs last), as [(name, cost_or_None, reason)]. The caller passes the
    chosen mapping to `generate(backend={order: name, ...})`; the best backend can change with order."""
    table = {}
    for p in orders:
        rows = []
        for b in _REGISTRY.values():
            ok, why = b.supports(problem)
            rows.append((b.name, b.cost(problem, p) if ok else None, "" if ok else why, ok))
        rows.sort(key=lambda r: (not r[3], r[1] is None, r[1] if r[1] is not None else 0))
        table[p] = [(n, c, w) for n, c, w, _ in rows]
    return table


def generate(kernel, order, name, backend="cartesian", source_order=0, nchannels=1, orders=None, language="c",
             manifest=True, **options):
    """Generate operators for `kernel`; returns the Manifest.

    backend is a registered name, or a mapping {order: name, "default": name} to use different backends
    at different orders (see `recommend`). Backends coexist: in a mixed build each backend writes its
    own orders to {name}_{backend}.* with every function prefixed "{backend}_" (e.g. cartesian_M2L_4), and
    the manifest lists the parts and which backend serves which order. Backends that cannot emit
    prefixed functions cannot be mixed.
    A VectorKernel is lowered to its scalar groups, each generated as {name}_g{g}.
    The manifest is also written to {name}.json (in src_dir if given) unless manifest=False.
    """
    from .vector import VectorKernel

    def one(pr, nm, base=""):
        if language == "fortran":
            base = ""  # a Fortran group is namespaced by its module, which is named after `nm`
        ords = list(orders) if orders is not None else list(range(pr.source_order + 1, order + 1))
        chosen = _resolve(pr, backend, ords)
        opts = {k: v for k, v in options.items() if k not in ("potential", "field")}
        names = list(dict.fromkeys(b.name for b in chosen.values()))
        if len(names) == 1:
            opts["prefix"] = base + opts.get("prefix", "")
            return chosen[ords[0]].generate(pr, order, nm, orders=ords, language=language, **opts)
        # mixed build: every backend writes its own orders under its own prefix, so all of it links together
        if "prefix" in opts:
            raise ValueError("prefix cannot be combined with a mixed-backend build (prefixes are the backend names)")
        parts = []
        for bn in names:
            b = get_backend(bn)
            if not b.supports_prefix:
                raise NotImplementedError(f"backend {bn!r} cannot emit prefixed functions, so it cannot be mixed")
            sub = [q for q in ords if chosen[q].name == bn]
            parts.append(b.generate(pr, order, f"{nm}_{bn}", orders=sub, language=language,
                                    prefix=f"{base}{bn}_", **opts))
        return Manifest("mixed", parts[0].kernel, language, ords, pr.source_order, pr.nchannels, ["operator"],
                        [f for m in parts for f in m.files], [], {}, sorted({d for m in parts for d in m.dependencies}),
                        by_order={q: chosen[q].name for q in ords}, parts=[m.to_dict() for m in parts])

    pot, fld = options.get("potential", True), options.get("field", True)
    if isinstance(kernel, VectorKernel):
        # each group is its own operator set: distinct files and distinct function prefixes, so the
        # groups of one vector kernel link into a single program
        subs = [one(Problem(g.kernel, g.source_order, g.nchannels, pot, fld), f"{name}_g{i}", f"{name}_g{i}_")
                for i, g in enumerate(kernel.groups)]
        man = Manifest(subs[0].backend, kernel.contract, language, subs[0].orders, 0, 0, subs[0].tiers,
                       [f for m in subs for f in m.files], [], {}, [], groups=[m.to_dict() for m in subs])
    else:
        man = one(Problem(kernel, source_order, nchannels, pot, fld), name)
    if manifest:
        d = options.get("src_dir") or options.get("include_dir") or "."
        man.write(os.path.join(d, f"{name}.json"))
    return man
