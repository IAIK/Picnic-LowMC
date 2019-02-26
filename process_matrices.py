#! /usr/bin/env sage

from __future__ import unicode_literals

from sage.all import GF, matrix, vector, copy
from six.moves import range
from generate_matrices import Instance
import pickle, io
import sys
import math

F = GF(2)


def calc_dot_from_full_rank(mat):
    """
    calculate the matrix \dot{M} from M, where M is M^1*
    in our implementation, M is horizontally and vertically mirrored and transposed
    """
    dot_size = mat.ncols()
    M_dot = matrix(F, dot_size, dot_size)
    current_rank = 0
    idx = 0
    indices = []
    while current_rank < dot_size:
        if idx >= mat.nrows():
            raise RuntimeError("mat should have full rank, something is wrong")
        M_dot[current_rank, :] = mat[idx, :]
        if current_rank == M_dot.rank():
            indices.append(idx)
        current_rank = M_dot.rank()
        idx += 1

    return M_dot, indices


def snip_r_wedge(mat_wedge, indices):
    """
    Calculate the minimal representation of R_wedge, and the corresponding shuffle information
    """
    snipped_rows = mat_wedge.nrows() - mat_wedge.ncols()
    M_snip = mat_wedge[mat_wedge.nrows() - snipped_rows : mat_wedge.nrows(), :]

    idx = 0
    for correction in indices:
        M_snip[idx, :] = mat_wedge[correction, :]
        idx += 1

    return M_snip


def gen_masks_from_cols(cols, num_sboxes, statesize):
    word = 0
    set_bits = 3 * num_sboxes
    for col in cols:
        if col < (statesize - 64):
            raise ValueError("Invalid col")
        word |= 1 << (col % 64)
        set_bits -= 1
    for i in range(set_bits):
        word |= 1 << (63 - i)
    return word


def combine_words(w):
    return long("".join(reversed([str(x) for x in w])), base=2)


def uint_constant_fmtstr(width):
    hexwidth = width / 4
    return "UINT{}_C(0x{{:0{}x}})".format(width, hexwidth)


def print_vector(output, name, typename, m, width=64):
    cols = len(m)
    rcols = cols / width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name} = {{ {{".format(type=typename, name=name))
    tmp = []
    for b in range(0, cols, width):
        w = m[b : b + width]
        tmp.append(formatstr.format(combine_words(w)))
    output.write(", ".join(tmp))
    output.write("} };\n")


def print_matrix(output, name, typename, m, width=64):
    rows = m.nrows()
    cols = m.ncols()
    rcols = cols / width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name}".format(type=typename, name=name))
    output.write("[{}] = {{\n".format(rows))
    for r in m.rows():
        output.write("  { {")
        tmp = []
        for b in range(0, cols, width):
            w = r[b : b + width]
            tmp.append(formatstr.format(combine_words(w)))
        output.write(", ".join(tmp))
        output.write("} },\n")
    output.write("};\n")


def calc_rowstride(rcols, width):
    bound = 128 / width
    if rcols > bound:
        return ((rcols * (width / 8) + 31) & ~31) / (width / 8)
    else:
        return ((rcols * (width / 8) + 15) & ~15) / (width / 8)


def print_row_t(output, entries, formatstr):
    while len(entries) % 4 != 0:
        entries.append(formatstr.format(0))
    strentries = []
    for idx in range(len(entries) / 4):
        strentries.append(
            "{{{{ {} }}}}".format(", ".join(entries[idx * 4 : (idx + 1) * 4]))
        )
    output.write(", ".join(strentries))


def print_matrix_mzd(output, name, typename, m, width=64):
    rows = m.nrows()
    cols = m.ncols()
    rcols = (cols + width - 1) / width
    rowstride = calc_rowstride(rcols, width)
    totalcols = rowstride * width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name}[]".format(type=typename, name=name))
    output.write(" = {\n")

    take_rows = 2 if rowstride == 2 else 1
    for idx in range(0, rows, take_rows):
        tmp = []
        for j in range(idx, min(rows, idx + take_rows)):
            r = m.row(j)
            for b in range(0, totalcols, width):
                w = r[b : b + width] if b < cols else [0] * width
                tmp.append(formatstr.format(combine_words(w)))
        output.write("  ")
        print_row_t(output, tmp, formatstr)
        output.write(",\n")
    output.write("};\n")


def print_vector_mzd(output, name, typename, m, width=64):
    rows = 1
    cols = len(m)
    rcols = (cols + width - 1) / width
    rowstride = calc_rowstride(rcols, width)
    totalcols = rowstride * width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name}[]".format(type=typename, name=name))
    output.write(" = {\n")
    output.write("  ")
    tmp = []
    for b in range(0, totalcols, width):
        w = m[b : b + width] if b < cols else [0] * width
        tmp.append(formatstr.format(combine_words(w)))
    print_row_t(output, tmp, formatstr)
    # output.write(",\n")
    output.write("};\n")


def print_mzd(w, width=64):
    b64 = []
    for t in range(0, len(w), width):
        s = w[t : t + width]
        b4 = []
        for j in range(0, width, 4):
            sj = s[j : j + 4]
            b4.append("".join("1" if x == 1 else " " for x in sj))
        b64.append(":".join(b4))
    print("[" + "|".join(b64) + "]")


def main(blocksize=256, keysize=256, rounds=19, sboxes=10):
    with io.open(
        "matrices_and_constants_{}_{}_{}.pickle".format(blocksize, keysize, rounds),
        "rb",
    ) as matfile:
        inst = pickle.load(matfile)

    if inst.n != blocksize or inst.k != keysize or inst.r != rounds:
        raise ValueError("Unexpected LowMC instance.")
    if blocksize != keysize:
        raise ValueError("Only blocksize == keysize is currently supported!")
    if sboxes not in (10, 1):
        raise ValueError("Pre-computation only implemented for m = 10 and m = 1")

    P = matrix(F, blocksize, blocksize)
    for i in range(blocksize):
        P[blocksize - i - 1, i] = 1

    Ls = [P * matrix(F, L) * P for L in inst.L]
    Ks = [P * matrix(F, K) * P for K in inst.K]
    Cs = [P * vector(F, C) for C in inst.R]

    Lt = [m.transpose() for m in Ls]
    K0t = Ks[0].transpose()

    Li = [m.inverse() for m in Lt]
    LiK = [Ks[i + 1].transpose() * Li[i] for i in range(inst.r)]
    LiC = [Cs[i] * Li[i] for i in range(inst.r)]

    mod_Li = [copy(Li[i]) for i in range(inst.r)]
    for j in range(inst.r):
        mod_Li[j][inst.n - 3 * sboxes :, : inst.n] = matrix(F, 3 * sboxes, inst.n)

    if sboxes == 10:
        precomputed_key_matrix = None
        precomputed_key_matrix_nl = matrix(F, inst.n, (sboxes * 3 + 2) * inst.r)
        precomputed_constant = None
        precomputed_constant_nl = vector(F, (sboxes * 3 + 2) * inst.r)

        for round in range(inst.r):
            tmp = copy(LiK[round])
            tmpC = copy(LiC[round])

            for i in range(round + 1, inst.r):
                x = LiK[i]
                c = LiC[i]
                for j in range(i - 1, round - 1, -1):
                    x = x * mod_Li[j]
                    c = c * mod_Li[j]
                tmp += x
                tmpC += c

            # non-linear part
            idx = round * (3 * sboxes + 2)
            precomputed_key_matrix_nl[:inst.n, idx + 2:idx + 3 * sboxes + 2] = tmp[:inst.n, inst.n - 3*sboxes:]
            precomputed_constant_nl[idx + 2:idx + 3 * sboxes + 2] = tmpC[inst.n - 3 * sboxes:]

            # linear part
            if round == 0:
                tmp[:, inst.n - 3 * sboxes :] = matrix(F, inst.n, 3 * sboxes)
                tmpC[inst.n - 3 * sboxes :] = vector(F, 3 * sboxes)
                precomputed_key_matrix = tmp + K0t
                precomputed_constant = tmpC
    elif sboxes == 1:
        num_64bit_words = (inst.r + 20) / 21
        precomputed_key_matrix = None
        precomputed_key_matrix_nl = matrix(F, inst.n, 64 * num_64bit_words)
        precomputed_constant = None
        precomputed_constant_nl = vector(F, 64 * num_64bit_words)

        for round in range(inst.r):
            tmp = copy(LiK[round])
            tmpC = copy(LiC[round])

            for i in range(round + 1, inst.r):
                x = LiK[i]
                c = LiC[i]
                for j in range(i - 1, round - 1, -1):
                    x = x * mod_Li[j]
                    c = c * mod_Li[j]
                tmp += x
                tmpC += c

            # non-linear part
            idx = 1 + ((round % 21) * 3) + 64 * (round / 21)
            precomputed_key_matrix_nl[:inst.n, idx:idx+3] = tmp[:inst.n, inst.n - 3*sboxes:]
            precomputed_constant_nl[idx:idx + 3] = tmpC[inst.n - 3 * sboxes:]

            # linear part
            if round == 0:
                tmp[:, inst.n - 3 * sboxes :] = matrix(F, inst.n, 3 * sboxes)
                tmpC[inst.n - 3 * sboxes :] = vector(F, 3 * sboxes)
                precomputed_key_matrix = tmp + K0t
                precomputed_constant = tmpC
    # RRKC precomputation done

    Z_i = []
    R_full = []
    R_wedge_snipped = []
    R_dot = []
    R_dot_inv = []
    R_wedge = []
    R_cols = []
    R_masks = []
    T_vee = []

    # i = 1
    R_full.append(Lt[0][:, 0 : inst.n - 3 * sboxes])
    Rdot, colR = calc_dot_from_full_rank(R_full[0])
    R_dot.append(Rdot)
    R_cols.append(colR)
    R_masks.append(gen_masks_from_cols(colR, sboxes, inst.n))
    R_dot_inv.append(R_dot[0].inverse())

    R_wedge.append(R_full[0] * R_dot_inv[0])
    Z_i.append(Lt[0][:, inst.n - 3 * sboxes : inst.n])
    R_wedge_snipped.append(snip_r_wedge(R_wedge[0], colR))

    # i = 2...r-1
    for i in range(1,rounds-1):
        T_vee.append(R_dot[i-1] * Lt[i][0:inst.n-sboxes*3,:])
        R = matrix(F, inst.n, inst.n-sboxes*3)
        R[0:inst.n-3*sboxes, :] = T_vee[-1][0:inst.n-3*sboxes, 0:inst.n-3*sboxes]
        R[inst.n-3*sboxes:inst.n, :] = Lt[i][inst.n-3*sboxes:inst.n, 0:inst.n-3*sboxes]
        R_full.append(R)
        Rdot, colR = calc_dot_from_full_rank(R_full[i])
        R_dot.append(Rdot)
        R_cols.append(colR)
        R_masks.append(gen_masks_from_cols(colR, sboxes, inst.n))
        R_dot_inv.append(R_dot[i].inverse())
        R_wedge.append(R_full[i]*R_dot_inv[i])
        Z_i.append(T_vee[i-1][0:inst.n-3*sboxes, inst.n-3*sboxes:inst.n].transpose().augment(Lt[i][inst.n-sboxes*3:inst.n, inst.n-sboxes*3:inst.n].transpose()).transpose())
        R_wedge_snipped.append(snip_r_wedge(R_wedge[-1],colR))

    # i = r
    T_vee.append(R_dot[rounds-2] * Lt[rounds-1][0:inst.n-sboxes*3,:])
    Z_r = T_vee[-1].transpose().augment(Lt[-1][inst.n-sboxes*3:inst.n, :].transpose()).transpose()
    # RLL precomputation done

    with io.open("lowmc_{}_{}_{}.h".format(inst.n, inst.k, inst.r), "w") as matfile:
        matfile.write('''#ifndef LOWMC_{inst.n}_{inst.k}_{inst.r}_H
#define LOWMC_{inst.n}_{inst.k}_{inst.r}_H

#include "lowmc_pars.h"

#if !defined(MUL_M4RI)
extern const lowmc_t lowmc_{inst.n}_{inst.k}_{inst.r};
#else
extern lowmc_t lowmc_{inst.n}_{inst.k}_{inst.r};
#endif

#endif
'''.format(s=inst.n / 64, inst=inst))

    with io.open('lowmc_{}_{}_{}.c'.format(inst.n, inst.k, inst.r), 'w') as matfile:
        matfile.write('''#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stddef.h>

#include "lowmc_{inst.n}_{inst.k}_{inst.r}.h"

'''.format(inst=inst))

        typename = 'block_t'

        matfile.write('#if !defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)\n')
        for i, L in enumerate(Ls):
            print_matrix_mzd(matfile, 'L_{}'.format(i),
                    typename, L.transpose())
            matfile.write('\n')
        matfile.write('#endif\n')

        matfile.write('#if !defined(REDUCED_ROUND_KEY_COMPUTATION)')
        for i, K in enumerate(Ks):
            matfile.write('\n')
            print_matrix_mzd(matfile, 'K_{}'.format(i),
                    typename, K.transpose())

        for i, C in enumerate(Cs):
            matfile.write('\n')
            print_vector_mzd(matfile, 'C_{}'.format(i),
                    typename, C)
        matfile.write('#endif\n')

        matfile.write('#if defined(REDUCED_ROUND_KEY_COMPUTATION)\n')
        print_matrix_mzd(matfile, 'precomputed_round_key_matrix_linear_part',
                    typename, precomputed_key_matrix)
        matfile.write('\n')

        print_matrix_mzd(matfile, 'precomputed_round_key_matrix_non_linear_part',
                    typename, precomputed_key_matrix_nl)
        matfile.write('\n')

        print_vector_mzd(matfile, 'precomputed_constant_linear_part', typename, precomputed_constant)
        matfile.write('\n')
        print_vector_mzd(matfile, 'precomputed_constant_non_linear_part', typename, precomputed_constant_nl)

        matfile.write('\n')
        matfile.write('#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)\n')
        print_matrix_mzd(matfile, 'Z_r', typename, Z_r)
        for i, Z in enumerate(Z_i):
            matfile.write('\n')
            print_matrix_mzd(matfile, 'Zi_{}'.format(i), typename, Z.transpose())

        for i, R in enumerate(R_wedge_snipped):
            matfile.write('\n')
            print_matrix_mzd(matfile, 'Ri_{}'.format(i), typename, R)

        matfile.write('#endif\n')
        matfile.write('#endif\n\n')

        matfile.write(
'''#if defined(MUL_M4RI)
static lowmc_round_t rounds[{inst.r}] = {{
#else
static const lowmc_round_t rounds[{inst.r}] = {{
#endif
'''.format(inst=inst))
        for i in range(inst.r):
            if i != inst.r-1:
                matfile.write(
'''
  {{
#if !defined(REDUCED_ROUND_KEY_COMPUTATION)
#if defined(MUL_M4RI)
    K_{j}, L_{i}, C_{i}, NULL, NULL
#else
    K_{j}, L_{i}, C_{i}
#endif
#else
#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)
    Zi_{i}, Ri_{i}, {R_mask},
#else
#if defined(MUL_M4RI)
    L_{i}, NULL
#else
    L_{i}
#endif
#endif
#endif
  }},'''.format(inst=inst, i=i, j=i+1, R_mask=uint_constant_fmtstr(64).format(R_masks[i])))
            else:
                matfile.write(
'''
  {{
#if !defined(REDUCED_ROUND_KEY_COMPUTATION)
#if defined(MUL_M4RI)
    K_{j}, L_{i}, C_{i}, NULL, NULL
#else
    K_{j}, L_{i}, C_{i}
#endif
#else
#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)
    NULL, NULL, 0,
#else
#if defined(MUL_M4RI)
    L_{i}, NULL
#else
    L_{i}
#endif
#endif
#endif
  }},'''.format(inst=inst, i=i, j=i+1))

        matfile.write(
'''
}};

#if defined(MUL_M4RI)
lowmc_t lowmc_{inst.n}_{inst.k}_{inst.r} = {{
#else
const lowmc_t lowmc_{inst.n}_{inst.k}_{inst.r} = {{
#endif
  {m}, {inst.n}, {inst.r}, {inst.k},
#if defined(REDUCED_ROUND_KEY_COMPUTATION)
  precomputed_round_key_matrix_linear_part,
#else
  K_0,
#endif
#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)
  Z_r,
#endif
#if defined(MUL_M4RI)
  NULL,
#endif
  rounds,
#if defined(REDUCED_ROUND_KEY_COMPUTATION)
  precomputed_round_key_matrix_non_linear_part,
#if defined(MUL_M4RI)
  NULL,
#endif
  precomputed_constant_linear_part,
  precomputed_constant_non_linear_part,
#endif
}};

'''.format(inst=inst, m=sboxes))

if __name__ == "__main__":
    import sys

    if len(sys.argv) == 5:
        blocksize, keysize, rounds, sboxes = map(int, sys.argv[1:5])
        main(blocksize, keysize, rounds, sboxes)
    else:
        main()
