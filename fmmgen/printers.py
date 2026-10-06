from sympy.printing.c import C99CodePrinter as C99Base
from sympy.printing.cxx import CXX11CodePrinter as CXX11Base
from sympy.printing.fortran import FCodePrinter as FBase
import logging
import re
import sympy as sp
from sympy import cse
from sympy.polys.polyfuncs import horner as sp_horner
from fmmgen.opts import basic as opts
from sympy import count_ops


logger = logging.getLogger(name="fmmgen")


def integral_exponent(e):
    """Return int(e) if e is an integer-valued real number, else None.

    Phi_derivatives builds R as (dx**2 + dy**2 + dz**2)**(0.5) using a Python
    float, so every exponent in the derived expressions is a sympy Float rather
    than an Integer, and Float(3.0).is_integer is False. Testing is_integer
    directly therefore disabled the pow -> multiplication replacement below
    entirely: minpow silently did nothing and 167 pow() calls survived in the
    generated operators, two of them per P2P call.

    Note the caller needs a genuine Python int, since it indexes range().
    """
    if e.is_Integer:
        return int(e)
    if e.is_Number and e.is_real and not e.is_infinite:
        f = float(e)
        if f == int(f):
            return int(f)
    return None


class CCodePrinter(C99Base):
    def __init__(self, settings={}, minpow=False):
        super(C99Base, self).__init__(settings)
        self.minpow = minpow

    def _print_Pow(self, expr):
        if self.minpow:
            n = integral_exponent(expr.exp)
            if n is not None and 0 < n <= self.minpow:
                base = self._print(expr.base)
                return "(" + "*".join([base] * n) + ")"

            elif n is not None and -self.minpow <= n < 0:
                base = self._print(expr.base)
                return "(1 / (" + "*".join([base] * abs(n)) + "))"
            else:
                return super()._print_Pow(expr)
        else:
            return super()._print_Pow(expr)


class CXXCodePrinter(CXX11Base):
    def __init__(self, settings={}, minpow=False):
        super(CXX11Base, self).__init__(settings)
        self.minpow = minpow

    def _print_Pow(self, expr):
        if self.minpow:
            n = integral_exponent(expr.exp)
            if n is not None and 0 < n <= self.minpow:
                base = self._print(expr.base)
                return "(" + "*".join([base] * n) + ")"

            elif n is not None and -self.minpow <= n < 0:
                base = self._print(expr.base)
                return "(1 / (" + "*".join([base] * abs(n)) + "))"
            else:
                return super()._print_Pow(expr)
        else:
            return super()._print_Pow(expr)


class FortranCodePrinter(FBase):
    """Free-form Fortran 90 expression printer.

    Differences from stock sympy that matter here:

    - Array elements print as a single 1-based index, `M(4)`, not `M(4, 1)`:
      every generated array is a flat vector that the driver addresses
      through assumed-size dummies.
    - Line wrapping is left to FortranStatementWriter, which has to count
      continuation lines (Fortran 90 allows only 39 of them).
    - Floats are `d0` literals, so this printer is double precision only.
    """

    def __init__(self, settings=None, minpow=False):
        settings = dict(settings or {})
        settings.setdefault("source_format", "free")
        settings.setdefault("standard", 90)
        super().__init__(settings)
        self.minpow = minpow
        # Hook for P2P_batch: S(k) there is the k-th moment of the u-th source,
        # which lives at S(source_size*(u-1) + k).
        self.moment_stride = None

    def _print_Integer(self, expr):
        # Coefficients reach 1e11 at order 12, past default-integer range
        # (Fortran rejects the literal; C silently promotes it). Real
        # literals everywhere also rule out integer division.
        return f"{int(expr)}.0d0"

    def _print_MatrixElement(self, expr):
        name = self._print(expr.parent)
        rows, cols = expr.parent.shape
        assert cols == 1, "generated arrays are column vectors"
        i = int(expr.i)
        if self.moment_stride is not None and name == "S":
            return f"S({self.moment_stride}*(u-1)+{i + 1})"
        return f"{name}({i + 1})"

    def _print_Pow(self, expr):
        if self.minpow:
            n = integral_exponent(expr.exp)
            if n is not None and 0 < n <= self.minpow:
                base = self._print(expr.base)
                return "(" + "*".join([base] * n) + ")"
            elif n is not None and -self.minpow <= n < 0:
                base = self._print(expr.base)
                return "(1.0d0/(" + "*".join([base] * abs(n)) + "))"
        # Phi_derivatives builds exponents as Python floats (see
        # integral_exponent), which would print as `x**3.0d0`, a pow() call.
        # An integer exponent lets the compiler expand it to multiplications.
        n = integral_exponent(expr.exp)
        if n is not None and n != 0:
            return f"{self.parenthesize(expr.base, 1000)}**({n})"
        return super()._print_Pow(expr)


class FortranStatementWriter:
    """Turns `target = target + expr` into wrapped free-form Fortran lines.

    A long Add is split into several statements, each at most `max_chars`
    printed characters, rather than one statement with hundreds of
    continuation lines: Fortran 90 allows 39 and Fortran 2003 255, and a
    compiler that enforces either limit would otherwise reject the P2P and
    M2L bodies at high order.
    """

    LINE = 90

    def __init__(self, printer, max_chars=2400):
        self.printer = printer
        self.max_chars = max_chars

    def _wrap(self, text, indent="  "):
        # Spaces only ever separate tokens in sympy's output, so breaking at
        # one can never split a literal such as 1.0d-9.
        out = []
        cur = indent
        for piece in text.split(" "):
            if len(cur) + len(piece) + 1 > self.LINE and cur.strip():
                out.append(cur.rstrip() + " &")
                cur = indent + "    "
            cur += piece + " "
        out.append(cur.rstrip())
        return "\n".join(out)

    def assign(self, target, expr, operator="="):
        """Statements for `target op expr`; operator is "=" or "+=" ."""
        terms = list(expr.args) if expr.is_Add else [expr]
        chunks, cur, size = [], [], 0
        for t in terms:
            n = len(self.printer.doprint(t))
            if cur and size + n > self.max_chars:
                chunks.append(cur)
                cur, size = [], 0
            cur.append(t)
            size += n
        if cur:
            chunks.append(cur)
        if not chunks:
            chunks = [[expr]]

        lines = []
        for k, chunk in enumerate(chunks):
            rhs = self.printer.doprint(sp.Add(*chunk, evaluate=False) if len(chunk) > 1 else chunk[0])
            # sympy's own wrapper may already have inserted "&" breaks.
            rhs = " ".join(part.strip().rstrip("&").strip() for part in rhs.split("\n"))
            if operator == "+=" or k > 0:
                lines.append(self._wrap(f"{target} = {target} + ({rhs})"))
            else:
                lines.append(self._wrap(f"{target} = {rhs}"))
        return "\n".join(lines) + "\n"


language_mapping = {
    "c": CCodePrinter,
    "c++": CXXCodePrinter,
    "fortran": FortranCodePrinter,
}


class SymbolIterator:
    def __init__(self, name):
        self.name = name
        self.num = 0

    def __iter__(self):
        return self

    def __next__(self):
        num = self.num
        self.num += 1
        return sp.Symbol(self.name + "tmp" + str(num))


class FunctionPrinter:
    def __init__(self, language="c", precision="double", debug=True, gpu=False, minpow=False, horner=False):
        logger.info(f'Function Printer created with precision "{precision}"')

        self.gpu = gpu
        if self.gpu:
            logger.info("Writing CUDA __device__ functions is enabled")

        if not debug:
            logger.info("CSE is enabled")
        else:
            logger.info("CSE is disabled")
        self.debug = debug
        self.horner = horner
        if self.horner:
            logger.info("Horner-form preprocessing is enabled")

        try:
            if minpow:
                self.printer = language_mapping[language](minpow=minpow)
            else:
                self.printer = language_mapping[language]()
        except KeyError:
            raise ValueError("Language not supported")

        self.language = language
        self.precision = precision
        assert self.precision in ["float", "double"]
        if language == "fortran":
            if precision != "double":
                raise NotImplementedError("Fortran output supports precision='double' only")
            if gpu:
                raise NotImplementedError("Fortran output does not support gpu=True")
            self.statements = FortranStatementWriter(self.printer)
        # Argument lists of every Fortran routine generated so far, in order,
        # as {name: [(argname, kind)]} with kind one of "scalar", "in", "inout".
        # The writer builds the order-dispatch wrappers from this.
        self.signatures = {}

    def _array(
        self,
        name,
        matrix,
        allocate=False,
        operator="=",
        atomic=False,
        ignore_symbols=[],
        coords=None,
        horner=None,
    ):
        opscount = 0
        code = ""

        # Horner-form preprocessing, entry by entry, before CSE ever sees the
        # matrix. This only pays off for entries that are genuinely
        # polynomials in the coordinates/Rinv with array-element coefficients
        # (S2M, M2M, L2L, L2P, M2P, the M2L derivative array D) -- measured
        # 20-45% fewer post-CSE ops at order 7-8 for those. It is a waste, not
        # just a no-op, on M2L's own M[i]*D[j] contraction: there is no
        # coordinate polynomial left to factor at that point (D is an opaque
        # placeholder array), and sympy's horner() still pays for a full
        # multivariate poly conversion over every M/D symbol in play --
        # measured 26s for zero benefit at order 7. That is why generate()'s
        # M2L call passes horner=False explicitly rather than relying on
        # this being a no-op.
        use_horner = self.horner if horner is None else horner
        if use_horner:
            matrix = sp.Matrix([sp_horner(e) if e.free_symbols else e for e in matrix])

        # R/Rinv's definition below must only reference coordinates the
        # function actually has as parameters. Every operator used to take
        # x, y AND z, so hardcoding all three was never wrong -- the planar
        # (2D-plane) 2-argument P2P/P2P_batch variants are the first
        # functions generated without z in scope at all, and referencing it
        # anyway is a compile error, not a warning.
        if coords is None:
            coords = sp.symbols("x y z")

        # Rinv is emitted as 1.0/sqrt(...) rather than pow(..., -0.5).
        #
        # An earlier note here read: "Testing on Godbolt with GCC 9.1 and ICPC
        # shows that pow(x, 0.5) generates fewer instructions than sqrt(x), so
        # will swap." That is misleading: instruction count at the call site is
        # not cost, because a libm call is one instruction that dispatches to
        # hundreds. Checked again with -O3:
        #
        #   g++-15  : pow(x,-0.5) emits a real libm call;
        #             1.0/sqrt(x) emits hardware fsqrt, no call.
        #   clang++ : both forms fold to hardware fsqrt.
        #
        # So the explicit form is a large win on GCC and a no-op on clang, i.e.
        # it cannot regress. Measured 2.5-2.9x on the P2P operator, which is
        # dominated by this single expression. ICPC not retested.
        #
        # R is left as sqrt() in case it gets used in expansions in future.

        r_squared = " + ".join(f"{s}*{s}" for s in coords)
        if sp.symbols("R") in matrix.free_symbols:
            code += f"{self.precision} R = sqrt({r_squared});\n"
        if sp.symbols("Rinv") in matrix.free_symbols:
            code += f"{self.precision} Rinv = 1.0 / sqrt({r_squared});\n"

        if allocate:
            code += f"{self.precision} {name}[{len(matrix)}];\n"

        if not self.debug:
            # print('Printing with CSE')
            iterator = SymbolIterator(name)
            # print(f'ignoring {name} in cse')
            # Stock sympy.cse. fmmgen previously vendored a patched copy of
            # sympy's tree_cse/cse (fmmgen/cse.py) to add `ignore` and
            # `light_ignore` filtering, but neither was doing anything:
            #
            #  - inserting the light_ignore loop between the `ignore` loop and
            #    its `for...else` detached them, so breaking out of the
            #    `ignore` loop no longer skipped elimination;
            #  - `light_ignore` compared against whole expressions, but bare
            #    Symbols return early as atoms, so it never matched. Its only
            #    effect was a debug print on every subexpression.
            #
            # Verified byte-identical generated output across p = 1..11 with
            # stock cse and no ignore arguments, so the vendored copy and both
            # parameters were removed. This also drops a dependency on sympy
            # internals (opt_cse, preprocess_for_cse, Unevaluated, ...).
            sub_expressions, rmatrix = cse(
                matrix,
                optimizations=opts,
                symbols=iterator,
            )

            rmatrix = sp.Matrix(rmatrix)
            for i, (var, sub_expr) in enumerate(sub_expressions):
                opscount += count_ops(sub_expr)
                code += f"{self.precision} " + self.printer.doprint(sub_expr, assign_to=var) + "\n"

            opscount += count_ops(rmatrix)
            tmp = self.printer.doprint(rmatrix, assign_to=name).replace("=", operator)

        else:
            # print('Printing without CSE')
            opscount += count_ops(matrix)
            tmp = self.printer.doprint(matrix, assign_to=name).replace("=", operator)

        if atomic:
            lines = tmp.split("\n")
            for line in lines:
                code += "#pragma omp atomic\n"
                code += line + "\n"
        else:
            code += tmp + "\n"
        return code, opscount

    def _generate_body(self, LHS, RHS, internal=[], operator="=", atomic=False, ignore=[], coords=None, horner=None):
        # Find the reduced RHS equation.
        opscount = 0
        logger.debug(f"Generating body for LHS = {str(LHS)}")
        code = ""

        # `horner` here overrides only the top-level LHS/RHS array below, not
        # `internal` arrays (e.g. M2L's D): those are always genuine
        # coordinate polynomials regardless of what the caller's own output
        # array looks like, so they keep the printer-level default.
        for arr_name, matrix in internal:
            codetext, ops = self._array(arr_name, matrix, allocate=True,
                                        ignore_symbols=[arr_name] + ignore, coords=coords)
            code += codetext
            opscount += ops

        codetext, ops = self._array(LHS, RHS, operator=operator, atomic=atomic, coords=coords, horner=horner)
        code += codetext
        opscount += ops
        return code, opscount

    def _generate_header(self, name, LHS, RHS, inputs):
        logger.debug(f"Generating headerfile for LHS = {str(LHS)}")
        # FMMGEN_RESTRICT tells the compiler none of these array arguments
        # ever overlap. True for every call site: the driver always passes
        # DISTINCT cells' M/L arrays and a separate F, never the same buffer
        # twice, so this is a free vectorisation hint rather than a behaviour
        # change. It is a macro rather than a bare keyword because the spelling
        # depends on the compiler mode, not on the `language` option: C99 has
        # `restrict`, C++ only the `__restrict` extension, and pre-C99 C has
        # nothing. Generated C is routinely compiled as C++ (pyximport with
        # CC=g++, or a C++ project including the header), so a fixed spelling
        # breaks one of the two. The macro is defined in the generated header.
        restrict = "FMMGEN_RESTRICT"
        ptr_type = f"{self.precision} * {restrict}"
        types = []
        for arg in map(type, inputs):
            if arg == sp.MatrixSymbol:
                types.append(ptr_type)
            else:
                types.append(self.precision)

        inputs.append(LHS)
        types.append(ptr_type)

        combined_inputs = ", ".join([str(x) + " " + str(y) for x, y in zip(types, inputs)])

        if self.gpu:
            return "__device__ void {}({})".format(name, combined_inputs)
        else:
            return "void {}({})".format(name, combined_inputs)

    # ------------------------------------------------------------------
    # Fortran 90 output
    #
    # Kept apart from the C path above so that path is untouched. The shape
    # differs enough to justify it: declarations must precede statements (so
    # every CSE temporary has to be collected and declared first), arrays are
    # 1-based assumed-size dummies, and aliasing between dummies is forbidden
    # by the language, which gives the compiler for free what FMMGEN_RESTRICT
    # asks for in C.
    # ------------------------------------------------------------------
    def _fortran_reduce(self, name, matrix, coords, horner):
        """Return (decls, temp_lines, exprs) for one output array.

        decls: names of the local scalars (R, Rinv, CSE temporaries) to declare.
        temp_lines: statements computing them, in dependency order.
        exprs: one reduced expression per entry of matrix.
        """
        use_horner = self.horner if horner is None else horner
        if use_horner:
            matrix = sp.Matrix([sp_horner(e) if e.free_symbols else e for e in matrix])

        decls, lines = [], []
        r_squared = " + ".join(f"{s}*{s}" for s in coords)
        if sp.symbols("R") in matrix.free_symbols:
            decls.append("R")
            lines.append(f"R = sqrt({r_squared})")
        if sp.symbols("Rinv") in matrix.free_symbols:
            decls.append("Rinv")
            lines.append(f"Rinv = 1.0d0 / sqrt({r_squared})")

        opscount = 0
        if not self.debug:
            sub_expressions, rmatrix = cse(matrix, optimizations=opts, symbols=SymbolIterator(name))
            exprs = list(sp.Matrix(rmatrix))
            for var, sub_expr in sub_expressions:
                opscount += count_ops(sub_expr)
                decls.append(str(var))
                lines.append(self.statements.assign(str(var), sub_expr).rstrip("\n"))
        else:
            exprs = list(matrix)
        opscount += count_ops(sp.Matrix(exprs))
        return decls, lines, exprs, opscount

    def _fortran_emit(self, target, exprs, operator):
        """Assignment statements for every entry of an output array."""
        lines = []
        for i, e in enumerate(exprs):
            if e == 0 and operator == "+=":
                continue
            lines.append(self.statements.assign(f"{target}({i + 1})", sp.sympify(e), operator).rstrip("\n"))
        return lines

    @staticmethod
    def _fortran_declarations(args, extra_int=(), local=()):
        """Declaration block for a routine with `args` = [(name, kind)]."""
        out = ["implicit none"]
        scalars = [n for n, k in args if k == "scalar"]
        ins = [n for n, k in args if k == "in"]
        inouts = [n for n, k in args if k == "inout"]
        ints = [n for n, k in args if k == "int"]
        if scalars:
            out.append("real(wp), intent(in) :: " + ", ".join(scalars))
        if ins:
            out.append("real(wp), intent(in) :: " + ", ".join(f"{n}(*)" for n in ins))
        if inouts:
            out.append("real(wp), intent(inout) :: " + ", ".join(f"{n}(*)" for n in inouts))
        if ints:
            out.append("integer, intent(in) :: " + ", ".join(ints))
        if extra_int:
            out.append("integer :: " + ", ".join(extra_int))
        # Several short statements, not one: a declaration list at high order
        # runs to hundreds of CSE temporaries, past the 132-column limit.
        for k in range(0, len(local), 8):
            out.append("real(wp) :: " + ", ".join(local[k:k + 8]))
        return out

    def _fortran_routine(self, name, args, decls, body):
        self.signatures[name] = args
        text = f"subroutine {name}({', '.join(n for n, _ in args)})\n"
        text += "\n".join("  " + d for d in decls) + "\n\n"
        text += "\n".join(body) + "\n"
        text += f"end subroutine {name}\n"
        return text

    def _fortran_generate(self, name, LHS, RHS, inputs, operator, internal, horner):
        coords = tuple(s for s in inputs if type(s) is not sp.MatrixSymbol)
        args = [(str(s), "in" if type(s) is sp.MatrixSymbol else "scalar") for s in inputs]
        args.append((LHS, "inout"))

        local, body, opscount = [], [], 0
        for arr_name, matrix in internal:
            d, tmp_lines, exprs, ops = self._fortran_reduce(arr_name, matrix, coords, None)
            opscount += ops
            local += d + [f"{arr_name}({len(matrix)})"]
            body += tmp_lines + self._fortran_emit(arr_name, exprs, "=")
        d, tmp_lines, exprs, ops = self._fortran_reduce(LHS, RHS, coords, horner)
        opscount += ops
        local += d
        body += tmp_lines + self._fortran_emit(LHS, exprs, operator)

        decls = self._fortran_declarations(args, local=local)
        code = self._fortran_routine(name, args, decls, ["  " + b.replace("\n", "\n  ") for b in body])
        header = f"subroutine {name}({', '.join(n for n, _ in args)})\n"
        return header, code, opscount

    def _fortran_generate_batch(self, name, LHS, RHS, symbols, source_size):
        n_out = len(RHS)
        acc = [f"{LHS.lower()}acc{i}" for i in range(n_out)]
        coords = sp.symbols(" ".join(symbols)) if len(symbols) > 1 else (sp.Symbol(symbols[0]),)
        # `begin`/`end` are C names; `end` is a Fortran keyword. Bounds are
        # 1-based and INCLUSIVE here, the natural Fortran convention, and an
        # empty range (ibeg > iend) is a no-op.
        args = ([(f"t{d}", "scalar") for d in symbols]
                + [(f"s{d}", "in") for d in symbols]
                + [("S", "in"), ("ibeg", "int"), ("iend", "int"), (LHS, "inout")])

        # The printer hook rewrites S(k) as it prints, so it must be active
        # for the CSE temporaries as well as the output entries.
        self.printer.moment_stride = source_size
        try:
            d, tmp_lines, exprs, opscount = self._fortran_reduce(LHS, RHS, coords, None)
            acc_lines = []
            for i, e in enumerate(exprs):
                if e != 0:
                    acc_lines.append(self.statements.assign(acc[i], sp.sympify(e), "+=").rstrip("\n"))
        finally:
            self.printer.moment_stride = None

        decls = self._fortran_declarations(args, extra_int=["u"], local=list(symbols) + acc + d)
        body = ["  " + f"{a} = 0.0d0" for a in acc]
        body.append("  !$omp simd reduction(+:" + ",".join(acc) + ")")
        body.append("  do u = ibeg, iend")
        body += [f"    {c} = t{c} - s{c}(u)" for c in symbols]
        body += ["    " + x.replace("\n", "\n    ") for x in tmp_lines + acc_lines]
        body.append("  end do")
        body += [f"  {LHS}({i + 1}) = {LHS}({i + 1}) + {a}" for i, a in enumerate(acc)]
        code = self._fortran_routine(name, args, decls, body)
        return f"subroutine {name}\n", code, opscount

    def generate_batch(self, name, LHS, RHS, symbols, source_size):
        """Emit a batched kernel: one target against a contiguous run of sources.

        The ordinary `generate` emits a function handling a single interaction.
        That function lives in the generated translation unit while the caller
        lives in another, so every interaction costs a real call -- and a call
        in the innermost loop makes vectorisation impossible in principle.
        Inspecting the object code confirmed it: zero vector instructions in
        the P2P loop and three un-inlined call sites.

        This emits the same expression inside a `#pragma omp simd` reduction
        loop over sources, with source coordinates taken from SoA arrays so the
        loads are unit-stride. The expression itself is untouched, so the kernel
        stays general over source order instead of being hand-specialised for
        the Coulomb monopole.
        """
        if self.language == "fortran":
            return self._fortran_generate_batch(name, LHS, RHS, symbols, source_size)
        n_out = len(RHS)
        acc = ["{}acc{}".format(LHS.lower(), i) for i in range(n_out)]
        pr = self.precision
        restrict = "FMMGEN_RESTRICT"  # see _generate_header

        args = ", ".join(
            ["{} t{}".format(pr, d) for d in symbols]
            + ["const {} * {} s{}".format(pr, restrict, d) for d in symbols]
            + ["const {} * {} S".format(pr, restrict), "size_t begin", "size_t end",
               "{} * {} {}".format(pr, restrict, LHS)]
        )
        header = "void {}({})".format(name, args)

        # symbols here IS the coordinate list (strings): 2-wide for the
        # planar batch kernels, 3-wide otherwise.
        coords = sp.symbols(" ".join(symbols)) if len(symbols) > 1 else (sp.Symbol(symbols[0]),)
        body, opscount = self._array(LHS, RHS, operator="+=", atomic=False, coords=coords)

        # Retarget the emitted body into the loop:
        #   S[k]    -> S[source_size*u + k]  (this source's moments)
        #   F[k] += -> facck +=              (private reduction accumulator)
        body = re.sub(
            r"\bS\[(\d+)\]",
            lambda m: "S[{}*u + {}]".format(source_size, m.group(1)),
            body,
        )
        for i in range(n_out):
            body = body.replace("{}[{}] +=".format(LHS, i), "{} +=".format(acc[i]))

        lines = [header + " {"]
        lines += ["{} {} = 0.0;".format(pr, a) for a in acc]
        lines.append("#pragma omp simd reduction(+:" + ",".join(acc) + ")")
        lines.append("for (size_t u = begin; u < end; u++) {")
        lines += ["{} {} = t{} - s{}[u];".format(pr, d, d, d) for d in symbols]
        lines.append(body.rstrip())
        lines.append("}")
        lines += ["{}[{}] += {};".format(LHS, i, a) for i, a in enumerate(acc)]
        lines.append("}")

        return header + ";\n", "\n".join(lines) + "\n", opscount

    def generate(self, name, LHS, RHS, inputs, operator="=", atomic=False, internal=[], ignore=[], horner=None):
        # Plain-Symbol inputs are coordinates (x, y, [z]); MatrixSymbol ones
        # are arrays (M, L, S, ...). Extracted before _generate_header, which
        # mutates `inputs` by appending LHS. Every operator used to take x, y
        # AND z, so R/Rinv's definition (built from this list, see _array)
        # was always safe to hardcode as 3-wide -- the planar 2-argument
        # P2P/P2P_batch variants are the first functions without z at all.
        if self.language == "fortran":
            return self._fortran_generate(name, LHS, RHS, list(inputs), operator, internal, horner)
        coords = tuple(s for s in inputs if type(s) is not sp.MatrixSymbol)
        header = self._generate_header(name, LHS, RHS, inputs)
        code = header + " {\n"
        codetext, opscount = self._generate_body(LHS, RHS, internal, operator, atomic=atomic,
                                                 ignore=ignore, coords=coords, horner=horner)
        code += codetext
        code += "\n}\n"
        header += ";\n"

        return header, code, opscount
