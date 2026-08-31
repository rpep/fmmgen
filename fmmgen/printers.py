from sympy.printing.c import C99CodePrinter as C99Base
from sympy.printing.cxx import CXX11CodePrinter as CXX11Base
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


language_mapping = {
    "c": CCodePrinter,
    "c++": CXXCodePrinter,
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
        # restrict (C99 keyword; C++ has no standard spelling, but every
        # compiler that matters -- GCC, Clang, MSVC, nvcc -- accepts the
        # `__restrict` extension in C++ mode) tells the compiler none of
        # these array arguments ever overlap. True for every call site: the
        # driver always passes DISTINCT cells' M/L arrays and a separate F,
        # never the same buffer twice, so this is a free vectorisation hint
        # rather than a behaviour change.
        restrict = "restrict" if self.language == "c" else "__restrict"
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
        n_out = len(RHS)
        acc = ["{}acc{}".format(LHS.lower(), i) for i in range(n_out)]
        pr = self.precision
        restrict = "restrict" if self.language == "c" else "__restrict"

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
        coords = tuple(s for s in inputs if type(s) is not sp.MatrixSymbol)
        header = self._generate_header(name, LHS, RHS, inputs)
        code = header + " {\n"
        codetext, opscount = self._generate_body(LHS, RHS, internal, operator, atomic=atomic,
                                                 ignore=ignore, coords=coords, horner=horner)
        code += codetext
        code += "\n}\n"
        header += ";\n"

        return header, code, opscount
