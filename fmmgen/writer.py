#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 08:41:10 2019

@author: ryan
"""
import sympy as sp
from fmmgen.printers import FunctionPrinter
from fmmgen.utils import Nterms
import textwrap
from fmmgen.harmonic import Nkeep
from fmmgen.generator import (
    _keep_sets,
    _keep_sets_planar,
    Nterms_planar,
    Nterms_planar_M,
    generate_mappings,
    generate_S2M_operators,
    generate_M_shift_operators,
    generate_L_operators,
    generate_L_shift_operators,
    generate_L2P_operators,
    generate_M2P_operators,
    generate_P2P_operators,
    generate_derivs,
    generate_S2M_operators_compressed,
    generate_M_shift_operators_compressed,
    generate_L_shift_operators_compressed,
    generate_L2P_operators_compressed,
    generate_M2P_operators_compressed,
    generate_M2L_compressed,
    generate_S2M_operators_planar,
    generate_M_shift_operators_planar,
    generate_L_shift_operators_planar,
    generate_L2P_operators_planar,
    generate_M2P_operators_planar,
    generate_M2L_planar,
)

import logging

logger = logging.getLogger(name="fmmgen")


q, x, y, z, R = sp.symbols("q x y z R")
symbols = (x, y, z)


def generate_code(
    order,
    name,
    precision="double",
    cython=False,
    CSE=False,
    harmonic_derivs=False,
    compress=False,
    planar=False,
    include_dir=None,
    src_dir=None,
    potential=True,
    field=True,
    source_order=0,
    atomic=False,
    gpu=False,
    minpow=0,
    language="c",
    save_opscounts=None,
):
    """
    Inputs:

    order, int:
        Expansion order of the FMM operators.

    name, str:
        Name of generated files.

        Note: Cython compiles *.pyx files to the *.c files with the same name,
        which means that we have to use a different name for the cython file.
        We therefore append "wrap" to the end of the name for the cython file,
        and this is therefore the name of the Python module which must be
        imported if using pyximport.

    cython_wrapper, bool:
        Enable generation of a Cython wrapper for the C files.

    CSE, bool:
        Enable common subexpression elimination, to reduce the op count in
        the generated code.

    harmonic_derivs, bool:
        The harmonicity of the Laplace means that at order p, there are only
        2p - 1 independent derivatives. Enabling this option therefore computes
        some derivatives as combinations of others.

    compress, bool:
        Store multipole and local expansions in the trace-free basis, which has
        (p+1)^2 coefficients rather than Nterms(p) = C(p+3,3). See
        fmmgen.harmonic. This shrinks M2L from C(p+6,6) multiplies to
        sum_{a+b<=p} (2a+1)(2b+1) -- measured 1.85x at p=5 and 2.61x at p=8 in
        wall clock on the kernel.

        Orthogonal to harmonic_derivs, which cheapens the derivative array
        instead; all four combinations are valid.

        The compressed operators are emitted ALONGSIDE the uncompressed ones in
        the same files, under a 'c' suffix -- M2Lc_5 beside M2L_5, and wrappers
        M2Lc(...) beside M2L(...). One header, one set of FMMGEN_* defines, so
        a caller can hold both and choose at runtime. P2P is not duplicated, as
        it never touches an expansion. The header additionally gains
        FMMGEN_COMPRESSED and FMMGEN_MULTIPOLESIZE/FMMGEN_LOCALSIZE, since the
        compressed arrays are (p+1)^2 rather than Nterms(order).

        With compress=False nothing above is altered, so the output is
        byte-identical to that of a build without this option.

        Supports source_order 0 and 1. The multipole and local arrays then have
        different retained sets -- s <= |n| <= p against |m| <= p - s -- which
        is handled, but at source_order >= 2 the source shell itself compresses
        (a quadrupole's 6 monomials to 5 retained, since the trace of a
        quadrupole moment contributes nothing to a Laplace field) and S2M would
        have to project the caller's moment array. Rejected rather than guessed.

    planar, bool:
        Sources (and, for the local array's retained set, evaluation points)
        confined to the z=0 plane, with moment vectors that remain full 3D --
        e.g. an atomistic dipole in a monolayer still has a (mux, muy, muz)
        moment. This is a positional degeneracy in the ordinary 1/R kernel,
        not the 2D Laplace log(r) kernel a translationally-invariant-along-z
        source would need.

        Every translation vector is a difference of z=0 coordinates and
        hence identically zero, which the multipole array's construction
        propagates into a hard bound: only n_z <= source_order can ever be
        nonzero (the source's own moment shell is the only place z-content
        enters; every displacement power in M_shift's recursion is dz**k
        with dz=0, killing k >= 1). The local array's bound is independent
        of source_order and instead follows from what L2P asks for: n_z <= 1
        if the field is wanted (Fz needs one extra z-derivative before
        evaluating at dz=0), else n_z <= 0. See fmmgen.generator's "Planar
        operator variants" section for the full derivation.

        Unlike compress, excluded entries are provably zero rather than a
        linear combination of retained ones, so there is no source_order cap
        here -- a quadrupole's own moment shell (n_z up to 2) is exactly what
        the retained set is built to include.

        Emitted ALONGSIDE the uncompressed (and, if also requested,
        compressed) operators under an 'xy' suffix -- M2Lxy_5 beside M2L_5
        (and M2Lc_5) -- guarded by FMMGEN_PLANAR and its own
        FMMGEN_PLANAR_MULTIPOLESIZE/FMMGEN_PLANAR_LOCALSIZE tables. Freely
        combinable with compress: the two are independent linear reductions
        of the same uncompressed operators, not alternatives, and requesting
        both does not produce some fourth hybrid array.

        With planar=False nothing above is altered, so the output is
        byte-identical to that of a build without this option.

    source_order, int:
        If source_order > 0 then we set certain multipole terms to zero
        in the local expansion, and hence they are not used. This is useful if,
        for example, we only have pure dipoles or quadrupoles in the system.

    minpow, int:
        If minpow is set, pow(x, n) expressions are expanded such that
        if n < minpow, the expression in code is printed as multiplications.

        e.g. if a sympy expression is pow(x, 2) + pow(y, 6) and minpow is 5,
        the printed version will be x*x + pow(y, 6)

    save_opcounts, string:
        Filename to save opcounts in
    """
    if save_opscounts:
        f = open(save_opscounts, "w")

    assert language in ["c", "c++"], "Language must be 'c' or 'c++'"
    if language == "c":
        fext = "c"
        hext = "h"
    if language == "c++":
        fext = "cpp"
        hext = "h"

    logger.info(f"Generating FMM operators to order {order}")
    assert precision in ["double", "float"], "Precision must be float or double"
    logger.info(f"Precision = {precision}")

    if compress and source_order > 1:
        raise NotImplementedError(
            "compress=True supports source_order 0 or 1; got "
            f"{source_order}. At degree >= 2 the source shell itself compresses "
            "(6 monomials -> 5 retained), so S2M would have to project the "
            "caller's moment array. See fmmgen.generator."
        )
    if compress:
        logger.info("Harmonic compression enabled")
    if CSE:
        logger.info("CSE Enabled")
        p = FunctionPrinter(precision=precision, debug=False, minpow=minpow)
    else:
        logger.info("CSE Disabled")
        p = FunctionPrinter(precision=precision, debug=True, minpow=minpow)

    header = ""
    body = ""

    x, y, z, q = sp.symbols("x y z q")
    symbols = (x, y, z)

    start = source_order
    if field:
        # No point starting at source_order
        # because no field calculation can be done
        # at this multipole order - the L2P derivative
        # is 0. This is true whether or not we
        # calculate potential, so we always increment
        # by 1 here the start order, such that the
        # field is always calculated.
        #
        # It may be preferable to instead compute
        # the field using a finite-difference approximation,
        # or to enable this as a general option.
        # However, this would be fairly sensitive, would
        # require another user parameter (lengths over which
        # to take the f.d. approximation), a decision on
        # whether to use forward/backward/central differencing,
        # and the order of the approximation. For now, we
        # leave it as a simple symbolic derivative.
        start += 1

    for i in range(start, order):
        print(f"Generating Order {i} operators")
        M_dict, _ = generate_mappings(i, symbols, "grevlex", source_order=source_order)
        L_dict, _ = generate_mappings(i - source_order, symbols, "grevlex", source_order=0)

        source_size = Nterms(source_order) - Nterms(source_order - 1)
        S2M = sp.Matrix(generate_S2M_operators(i, symbols, M_dict, source_order=source_order))

        head, code, S2M_opscount = p.generate(
            f"S2M_{i}",
            "M",
            S2M,
            list(symbols) + [sp.MatrixSymbol("S", source_size, 1)],
            operator="+=",
            atomic=atomic,
        )
        print(f"S2M_{i} opscount = {S2M_opscount}")
        header += head
        body += code + "\n"
        Ms = sp.Matrix(generate_M_shift_operators(i, symbols, M_dict, source_order=source_order))

        head, code, M2M_opscount = p.generate(
            f"M2M_{i}",
            "Ms",
            Ms,
            list(symbols) + [sp.MatrixSymbol("M", Nterms(i), 1)],
            operator="+=",
            atomic=atomic,
        )
        header += head
        body += code + "\n"
        print(f"M2M_{i} opscount = {M2M_opscount}")
        # Two stages here; generate derivs and then the L matrix. Both
        # must be passed to the function printer.
        derivs = sp.Matrix(generate_derivs(i, symbols, M_dict, source_order, harmonic_derivs=harmonic_derivs))
        L = sp.Matrix(generate_L_operators(i, symbols, M_dict, L_dict, source_order=source_order))

        head, code, M2L_opscount = p.generate(
            f"M2L_{i}",
            "L",
            L,
            list(symbols) + [sp.MatrixSymbol("M", Nterms(i), 1)],
            operator="+=",
            atomic=atomic,
            internal=[("D", derivs)],
        )
        header += head
        body += code + "\n"
        print(f"M2L_{i} opscount = {M2L_opscount}")

        Ls = sp.Matrix(generate_L_shift_operators(i, symbols, L_dict, source_order=source_order))
        head, code, L2L_opscount = p.generate(
            f"L2L_{i}",
            "Ls",
            Ls,
            list(symbols) + [sp.MatrixSymbol("L", Nterms(i), 1)],
            operator="+=",
            atomic=atomic,
        )

        header += head
        body += code + "\n"
        print(f"L2L_{i} opscount = {L2L_opscount}")
        L2P = generate_L2P_operators(i, symbols, L_dict, potential=potential, field=field)

        Fs = sp.Matrix(L2P)
        head, code, L2P_opscount = p.generate(
            f"L2P_{i}",
            "F",
            Fs,
            list(symbols) + [sp.MatrixSymbol("L", Nterms(i), 1)],
            operator="+=",
            atomic=atomic,
        )

        header += head
        body += code + "\n"
        print(f"L2P_{i} opscount = {L2P_opscount}")
        M2P = generate_M2P_operators(
            i,
            symbols,
            M_dict,
            potential=potential,
            field=field,
            source_order=source_order,
            harmonic_derivs=harmonic_derivs,
        )
        Fs = sp.Matrix(M2P)
        head, code, M2P_opscount = p.generate(
            f"M2P_{i}",
            "F",
            Fs,
            list(symbols) + [sp.MatrixSymbol("M", Nterms(i), 1)],
            operator="+=",
            atomic=atomic,
        )
        header += head
        body += code + "\n"
        print(f"M2P_{i} opscount = {M2P_opscount}")

        # Compressed variant, emitted ALONGSIDE the above rather than instead
        # of it. One header, one set of FMMGEN_* defines, nothing to include
        # twice and no redefinition to manage -- and every line above is left
        # untouched, so compress=False is byte-identical by construction rather
        # than by inspection.
        #
        # P2P is deliberately not duplicated: it never touches an expansion, so
        # the two variants would emit identical arithmetic.
        #
        # The uncompressed operator lists built above are handed straight in
        # rather than rebuilt: compression is a linear map applied to them, so
        # regenerating from sympy a second time is pure waste. Substitution
        # uses xreplace, not subs -- these are exact atom-for-atom maps and
        # need none of subs' pattern matching.
        if compress:
            # Separate retained sets: the multipole spans s <= |n| <= i while
            # the local spans |m| <= i - s, so these are different index sets
            # of different sizes whenever source_order > 0.
            keep_M, keep_L, dec_M, dec_L = _keep_sets(i, symbols, M_dict, L_dict, source_order)
            Mc = sp.MatrixSymbol("M", len(keep_M), 1)
            Lc = sp.MatrixSymbol("L", len(keep_L), 1)

            for nm, lhs, ops, inputs in (
                ("S2Mc", "M",
                 generate_S2M_operators_compressed(i, symbols, M_dict, keep_M,
                                                   source_order=source_order, dec_M=dec_M,
                                                   ops=list(S2M)),
                 list(symbols) + [sp.MatrixSymbol("S", source_size, 1)]),
                ("M2Mc", "Ms",
                 generate_M_shift_operators_compressed(i, symbols, M_dict, keep_M,
                                                       source_order=source_order, dec_M=dec_M,
                                                       ops=list(Ms)),
                 list(symbols) + [Mc]),
                ("L2Lc", "Ls",
                 generate_L_shift_operators_compressed(i, symbols, L_dict, keep_L,
                                                       source_order=source_order, dec_L=dec_L,
                                                       ops=list(Ls)),
                 list(symbols) + [Lc]),
                ("L2Pc", "F",
                 generate_L2P_operators_compressed(i, symbols, L_dict, keep_L,
                                                   potential=potential, field=field, dec_L=dec_L,
                                                   ops=list(L2P)),
                 list(symbols) + [Lc]),
                ("M2Pc", "F",
                 generate_M2P_operators_compressed(i, symbols, M_dict, keep_M,
                                                   potential=potential, field=field,
                                                   source_order=source_order,
                                                   harmonic_derivs=harmonic_derivs,
                                                   ops=list(M2P)),
                 list(symbols) + [Mc]),
            ):
                head, code, n = p.generate(
                    f"{nm}_{i}", lhs, sp.Matrix(ops), inputs,
                    operator="+=", atomic=atomic,
                )
                header += head
                body += code + "\n"
                print(f"{nm}_{i} opscount = {n}")

            # M2L returns its derivative array too: the compressed contraction
            # only reaches derivatives with third index <= 2, so D is reindexed
            # onto just those rather than emitted in full.
            dv, L_ops = generate_M2L_compressed(
                i, symbols, M_dict, L_dict, keep_M, keep_L,
                source_order=source_order, harmonic_derivs=harmonic_derivs,
                ops=list(L),
            )
            head, code, n = p.generate(
                f"M2Lc_{i}", "L", sp.Matrix(L_ops), list(symbols) + [Mc],
                operator="+=", atomic=atomic, internal=[("D", sp.Matrix(dv))],
            )
            header += head
            body += code + "\n"
            print(f"M2Lc_{i} opscount = {n}")

        # Planar variant, emitted ALONGSIDE the above -- independent of, and
        # freely combinable with, compress (see the planar docstring). Same
        # "hand in the already-built uncompressed lists" pattern.
        if planar:
            keep_M, keep_L = _keep_sets_planar(i, symbols, source_order, field=field)
            Mxy = sp.MatrixSymbol("M", len(keep_M), 1)
            Lxy = sp.MatrixSymbol("L", len(keep_L), 1)

            for nm, lhs, ops, inputs in (
                ("S2Mxy", "M",
                 generate_S2M_operators_planar(i, symbols, M_dict, keep_M,
                                               source_order=source_order,
                                               ops=list(S2M)),
                 list(symbols) + [sp.MatrixSymbol("S", source_size, 1)]),
                ("M2Mxy", "Ms",
                 generate_M_shift_operators_planar(i, symbols, M_dict, keep_M,
                                                   source_order=source_order,
                                                   ops=list(Ms)),
                 list(symbols) + [Mxy]),
                ("L2Lxy", "Ls",
                 generate_L_shift_operators_planar(i, symbols, L_dict, keep_L,
                                                   source_order=source_order,
                                                   ops=list(Ls)),
                 list(symbols) + [Lxy]),
                ("L2Pxy", "F",
                 generate_L2P_operators_planar(i, symbols, L_dict, keep_L,
                                               potential=potential, field=field,
                                               ops=list(L2P)),
                 list(symbols) + [Lxy]),
                ("M2Pxy", "F",
                 generate_M2P_operators_planar(i, symbols, M_dict, keep_M,
                                               potential=potential, field=field,
                                               source_order=source_order,
                                               harmonic_derivs=harmonic_derivs,
                                               ops=list(M2P)),
                 list(symbols) + [Mxy]),
            ):
                head, code, n = p.generate(
                    f"{nm}_{i}", lhs, sp.Matrix(ops), inputs,
                    operator="+=", atomic=atomic,
                )
                header += head
                body += code + "\n"
                print(f"{nm}_{i} opscount = {n}")

            dv, L_ops = generate_M2L_planar(
                i, symbols, M_dict, L_dict, keep_M, keep_L,
                source_order=source_order, harmonic_derivs=harmonic_derivs,
                ops=list(L),
            )
            head, code, n = p.generate(
                f"M2Lxy_{i}", "L", sp.Matrix(L_ops), list(symbols) + [Mxy],
                operator="+=", atomic=atomic, internal=[("D", sp.Matrix(dv))],
            )
            header += head
            body += code + "\n"
            print(f"M2Lxy_{i} opscount = {n}")

        if i == start:
            P2P = sp.Matrix(
                generate_P2P_operators(
                    symbols,
                    M_dict,
                    potential=potential,
                    field=field,
                    source_order=source_order,
                )
            )

            head, code, P2P_opscount = p.generate(
                "P2P",
                "F",
                P2P,
                list(symbols) + [sp.MatrixSymbol("S", Nterms(i), 1)],
                operator="+=",
                atomic=atomic,
            )
            print(f"P2P opscount = {P2P_opscount}")
            header += head
            body += code + "\n"

            # Batched form of the same operator: one target against a
            # contiguous run of sources, inside an omp simd reduction loop.
            # The per-pair P2P above cannot vectorise -- it is one call per
            # interaction, across translation units -- and P2P is the dominant
            # cost of the whole method at realistic particle counts.
            head, code, _ = p.generate_batch(
                "P2P_batch",
                "F",
                P2P,
                [str(sym) for sym in symbols],
                Nterms(source_order) - Nterms(source_order - 1),
            )
            header += head
            body += code + "\n"

        if save_opscounts:
            if i == start:
                f.write(f"P2P,{P2P_opscount}\n")
            f.write(f"S2M_{i},{S2M_opscount}\n")
            f.write(f"M2M_{i},{M2M_opscount}\n")
            f.write(f"M2L_{i},{M2L_opscount}\n")
            f.write(f"L2P_{i},{L2P_opscount}\n")
            f.write(f"L2L_{i},{L2L_opscount}\n")
            f.write(f"M2P_{i},{M2P_opscount}\n")

    # We now generate wrapper functions that cover all orders generated.
    unique_funcs = []
    func_definitions = header.split(";\n")
    # print(f"func_defs = {func_definitions}")
    for func in func_definitions:
        # Must do it this way in order to avoid breaking
        # for expansions > 10.
        function_name = func.split("(")[0]
        # print(f"Function_name = {function_name}")
        end_string = f"_{start}"
        if end_string == function_name[-len(end_string):]:
            # print("Unique!")
            unique_funcs.append(func)
        else:
            pass
            # print(f"{func} not unique")
            # print(f"  {end_string}  {function_name[-len(end_string)-1:]}")

    wrapper_funcs = [f.replace(")", ", int order)").replace(f"_{start}", "") for f in unique_funcs]

    #  print(wrapper_funcs)

    func_definitions += wrapper_funcs
    # print('\n'.join(func_definitions))

    for wfunc, func in zip(wrapper_funcs, unique_funcs):
        # Add to header file
        header += wfunc + ";\n"
        # Create a switch statement that covers all functions:
        code = wfunc + " {\n"
        code += "switch (order) {\n"
        for i in range(start, order):
            code += "  case {}:\n".format(i)
            # print(func)
            replaced_code = (
                func.replace(f"_{start}", f"_{i}")
                .replace("* ", "")
                .replace("double ", "")
                .replace("float ", "")
                .replace("void ", "")
            )
            # print(f"replaced_code: {replaced_code}")
            code += "    " + replaced_code + ";\n    break;\n"
        code += "  }\n}\n"
        # print(code)
        body += code

    if not include_dir:
        f = open(f"{name}.{hext}", "w")
    else:
        f = open(f"{include_dir.rstrip('/')}/{name}.{hext}", "w")
    f.write("#pragma once\n")
    # P2P_batch takes size_t range bounds.
    f.write("#include <cstddef>\n" if language == "c++" else "#include <stddef.h>\n")
    f.write(f"#define FMMGEN_MINORDER {start}\n")
    f.write(f"#define FMMGEN_MAXORDER {order}\n")
    f.write(f"#define FMMGEN_SOURCEORDER {source_order}\n")
    f.write(f"#define FMMGEN_SOURCESIZE {Nterms(source_order) - Nterms(source_order - 1)}\n")
    if potential and not field:
        osize = 1
    elif field and not potential:
        osize = 3
    elif field and potential:
        osize = 4
    f.write(f"#define FMMGEN_OUTPUTSIZE {osize}\n")
    if compress:
        # Only emitted when compressed. Uncompressed, the sizes are Nterms(order)
        # and callers already derive them that way -- adding these defines
        # unconditionally would change the header for existing builds, and
        # compress=False is required to be byte-identical.
        #
        # Indexed by order, so entry [order] is the size at that order; entries
        # below FMMGEN_MINORDER are present only to keep the indexing direct.
        msizes = ", ".join(str(Nkeep(o, source_order)) for o in range(order + 1))
        lsizes = ", ".join(str(Nkeep(o - source_order)) for o in range(order + 1))
        f.write("#define FMMGEN_COMPRESSED 1\n")
        f.write(
            "/* Array sizes for the *c-suffixed* operators only. The trace-free\n"
            "   basis drops the C(p+3,3) - (p+1)^2 redundant coefficients, so these\n"
            "   are (p+1)^2; the unsuffixed operators still take Nterms(order).\n"
            "   Index by order. */\n"
        )
        f.write(f"static const size_t FMMGEN_MULTIPOLESIZE[] = {{{msizes}}};\n")
        f.write(f"static const size_t FMMGEN_LOCALSIZE[] = {{{lsizes}}};\n")
    if planar:
        # Independent of FMMGEN_COMPRESSED -- both can be defined in the same
        # header, since the two variants never touch each other's arrays.
        # The multipole cutoff is source_order (not a fixed 1: see the planar
        # docstring), so unlike FMMGEN_MULTIPOLESIZE this table is NOT the
        # same at every source_order and must be computed per source_order,
        # not assumed to equal FMMGEN_MULTIPOLESIZE even when compress is
        # also set.
        pmsizes = ", ".join(str(Nterms_planar_M(o, source_order)) for o in range(order + 1))
        plsizes = ", ".join(
            str(Nterms_planar(o - source_order, 1 if field else 0)) for o in range(order + 1)
        )
        f.write("#define FMMGEN_PLANAR 1\n")
        f.write(
            "/* Array sizes for the *xy-suffixed* operators only. Retained set is\n"
            "   n_z <= source_order for the multipole array, n_z <= 1 (field) or\n"
            "   <= 0 (potential only) for the local array. Index by order. */\n"
        )
        f.write(f"static const size_t FMMGEN_PLANAR_MULTIPOLESIZE[] = {{{pmsizes}}};\n")
        f.write(f"static const size_t FMMGEN_PLANAR_LOCALSIZE[] = {{{plsizes}}};\n")
    f.write(header)
    f.close()

    if not src_dir:
        f = open(f"{name}.{fext}", "w")
    else:
        f = open(f"{src_dir.rstrip('/')}/{name}.{fext}", "w")

    f.write(f'#include "{name}.{hext}"\n')
    if language == "c":
        f.write('#include "math.h"\n')
    elif language == "c++":
        f.write("#include<cmath>\n")

    f.write(body)
    f.close()

    if cython and gpu:
        raise Warning("Cannot write a Cython wrapper for GPU code; skipping")

    elif cython:
        logger.info(f"Generating Cython wrapper: {name}_wrap.pyx")
        library = f"{name}"

        f = open(f"{name}_decl.pxd", "w")
        pxdcode = textwrap.dedent(
            """\
        cdef extern from "{}.h":
            cdef int FMMGEN_MINORDER
            cdef int FMMGEN_MAXORDER
            cdef int FMMGEN_SOURCEORDER
            cdef int FMMGEN_SOURCESIZE
            cdef int FMMGEN_OUTPUTSIZE
            {}
        """
        )
        f.write(pxdcode.format(name, "\n    ".join(func_definitions)))

        f.close()

        f = open(f"{name}_wrap.pyx", "w")
        # expose the C functions from the header file.
        pyxcode = textwrap.dedent(
            """\
        # cython: language_level=3
        cimport numpy as np
        cimport {}

        FMMGEN_MINORDER = {}.FMMGEN_MINORDER
        FMMGEN_MAXORDER = {}.FMMGEN_MAXORDER
        FMMGEN_SOURCEORDER = {}.FMMGEN_SOURCEORDER
        FMMGEN_SOURCESIZE = {}.FMMGEN_SOURCESIZE
        FMMGEN_OUTPUTSIZE = {}.FMMGEN_OUTPUTSIZE
        """
        ).format(*[name + "_decl"] * 6)

        # Underscores are stripped from the FUNCTION NAME only (so M2M_1 is
        # exposed as M2M1, as it always has been). Applying that substitution
        # to the whole signature, as this once did, also rewrites parameter
        # types -- "size_t" became "sizet" and the generated .pyx would not
        # compile as soon as any operator took one.
        for funcname in func_definitions:
            if not funcname:
                continue
            head, argstr = funcname.split("(", 1)
            argstr = argstr.rsplit(")", 1)[0]
            cname = head.split()[-1]
            pyxcode += "cpdef " + cname.replace("_", "") + "(" + argstr.replace(" *", "[:]") + "):\n"

            function_name = cname
            args = funcname.split("(")[1][:-1].split(",")
            processed_args = []

            for arg in args:
                if "*" in arg:
                    arg = arg.replace("* ", "&") + "[0]"
                processed_args.append(arg.split(" ")[-1])

            pyxcode += "    " + name + "_decl." + function_name + "(" + ", ".join(processed_args) + ")\n\n"

        f.write(pyxcode)
        f.close()

        f = open(f"{name}_wrap.pyxbld", "w")

        # print(library)

        logger.info(f"Generating Cython buildfile: {name}_wrap.pyxbld")
        bldcode = textwrap.dedent(
            """\
        import numpy as np

        def make_ext(modname, pyxfilename):
            from distutils.extension import Extension
            return Extension(name = modname,
                             sources=[pyxfilename, '{}'],
                             include_dirs=[np.get_include(), '.'],
                             library_dirs=['.'],
                             extra_link_args=[],
                             extra_compile_args=['-O3', '-fopenmp'])
        """
        ).format(library + ".c", library)

        f.write(bldcode)
        f.close()
