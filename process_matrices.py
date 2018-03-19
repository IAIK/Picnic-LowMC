from __future__ import unicode_literals

from sage.all import GF, matrix, vector, copy

from generate_matrices import Instance
import pickle, io
import sys
import math


def combine_words(w):
    return long("".join(reversed([str(x) for x in w])), base=2)


def uint_constant_fmtstr(width):
    hexwidth = width / 4
    return "UINT{}_C(0x{{:0{}x}})".format(width, hexwidth)


def print_vector(output, name, typename, m, width=64):
    cols = len(m)
    rcols = cols / width
    formatstr = uint_constant_fmtstr(width)

    output.write('static const {type} {name} = {{ {{'.format(type=typename, name=name))
    tmp = []
    for b in range(0, cols, width):
        w = m[b:b+width]
        tmp.append(formatstr.format(combine_words(w)))
    output.write(', '.join(tmp))
    output.write('} };\n');


def print_matrix(output, name, typename, m, width=64):
    rows = m.nrows()
    cols = m.ncols()
    rcols = cols / width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name}".format(type=typename, name=name))
    output.write("[{}] = {{\n".format(rows));
    for r in m.rows():
        output.write("  { {")
        tmp = []
        for b in range(0, cols, width):
            w = r[b:b+width]
            tmp.append(formatstr.format(combine_words(w)))
        output.write(', '.join(tmp))
        output.write("} },\n")
    output.write("};\n");

def calc_rowstride(rcols, width):
    bound = 128 / width
    if rcols >= bound:
        return ((rcols * (width / 8) + 31) & ~31) / (width / 8);
    else:
        return ((rcols * (width / 8) + 15) & ~15) / (width / 8);

def print_matrix_mzd(output, name, typename, m, width=64):
    rows = m.nrows()
    cols = m.ncols()
    rcols = (cols + width - 1) / width
    rowstride = calc_rowstride(rcols, width)
    totalcols = rowstride * width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name}".format(type=typename, name=name))
    output.write(" = {{ {}, {}, {}, {}, {{ 0 }}, {{\n".format(rows, cols, rcols, rowstride))
    for r in m.rows():
        output.write("  ")
        tmp = []
        for b in range(0, totalcols, width):
            w = r[b:b+width] if b < cols else [0] * width
            tmp.append(formatstr.format(combine_words(w)))
        output.write(', '.join(tmp))
        output.write(",\n")
    output.write("}};\n");

def print_vector_mzd(output, name, typename, m, width=64):
    rows = 1
    cols = len(m)
    rcols = (cols + width - 1) / width
    rowstride = calc_rowstride(rcols, width)
    totalcols = rowstride * width
    formatstr = uint_constant_fmtstr(width)

    output.write("static const {type} {name}".format(type=typename, name=name))
    output.write(" = {{ {}, {}, {}, {}, {{ 0 }}, {{\n".format(rows, cols, rcols, rowstride))
    output.write("  ")
    tmp = []
    for b in range(0, totalcols, width):
        w = m[b:b+width] if b < cols else [0] * width
        tmp.append(formatstr.format(combine_words(w)))
    output.write(', '.join(tmp))
    output.write(",\n")
    output.write("}};\n");


def print_mzd(w, width=64):
    b64 = []
    for t in range(0, len(w), width):
        s = w[t:t+width]
        b4 = []
        for j in range(0, width, 4):
            sj = s[j:j+4]
            b4.append(''.join('1' if x == 1 else ' ' for x in sj))
        b64.append(':'.join(b4))
    print('[' + '|'.join(b64) + ']')


def main(blocksize=256, keysize=256, rounds=19, sboxes=10):
    with io.open('matrices_and_constants_{}_{}_{}.pickle'.format(blocksize, keysize, rounds), 'rb') as matfile:
        inst = pickle.load(matfile)

    if inst.n != blocksize or inst.k != keysize or inst.r != rounds:
        raise ValueError("Unexpected LowMC instance.")
    if blocksize != keysize:
        raise ValueError("Only blocksize == keysize is currently supported!")

    F = GF(2)
    P = matrix(F, blocksize, blocksize)
    for i in range(blocksize):
        P[blocksize - i - 1, i] = 1

    Ls = [P * matrix(F, L) * P for L in inst.L]
    Ks = [P * matrix(F, K) * P for K in inst.K]
    Cs = [P * vector(F, C) for C in inst.R]

    Kt = [m.transpose() for m in Ks]
    Lt = [m.transpose() for m in Ls]
    Li = [m.inverse() for m in Lt]
    LiK = [Kt[i + 1] * Li[i] for i in xrange(inst.r)]

    mod_Li = [copy(Li[i]) for i in range(inst.r)]
    for j in range(inst.r):
        mod_Li[j][inst.n - 3*sboxes:, :inst.n] = matrix(3 * sboxes, inst.n)

    precomputed_values = [matrix(inst.n, inst.n) for i in range(inst.r)]
    precomputed_matrix = matrix(inst.n, (sboxes * 3 + 2) * inst.r)

    for round in range(inst.r):
        tmp = matrix(inst.n, inst.n)
        tmp += copy(LiK[round])

        for i in range(round + 1, inst.r):

            x = copy(LiK[i])
            for j in range(i - 1, round - 1, -1):
                    x = x * mod_Li[j]
            tmp += x

        precomputed_matrix[:inst.n, round * (3 * sboxes + 2): round * (3 * sboxes + 2) + 3 * sboxes] = tmp[:inst.n, inst.n - 3*sboxes:]
        if round:
            tmp[:,:inst.n - 3*sboxes] = matrix(inst.n, inst.n - 3*sboxes)
        else: 
            tmp[:,inst.n - 3*sboxes:] = matrix(inst.n, 3*sboxes)
        precomputed_values[round] = tmp

    # print(precomputed_matrix.str())

    with io.open('lowmc_{}_{}_{}.h'.format(inst.n, inst.k, inst.r), 'w') as matfile:
        matfile.write('''#ifndef LOWMC_{inst.n}_{inst.k}_{inst.r}
#define LOWMC_{inst.n}_{inst.k}_{inst.r}

#include <stdint.h>

#include "mzd_additional.h"


const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_linear_layer(uint32_t r);
const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_round_key(uint32_t r);
const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_round_const(uint32_t r);
const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_precomputed_round_key_matrix_non_linear_part();
const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_precomputed_round_key_matrix_linear_part();

#endif
'''.format(s=inst.n / 64, inst=inst, b = int(math.ceil(inst.r * (3 * sboxes + 2) / 64.0))))

    with io.open('lowmc_{}_{}_{}.c'.format(inst.n, inst.k, inst.r), 'w') as matfile:
        matfile.write('''#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stddef.h>

#include "lowmc_{inst.n}_{inst.k}_{inst.r}.h"

'''.format(inst=inst))

        typename = 'mzd_local_t'

        for i, L in enumerate(Ls):
            print_matrix_mzd(matfile, 'L_{}_{}_{}_{}'.format(inst.n, inst.k, inst.r, i),
                    typename, L.transpose())
            matfile.write('\n')

        for i, K in enumerate(Ks):
            print_matrix_mzd(matfile, 'K_{}_{}_{}_{}'.format(inst.n, inst.k, inst.r, i),
                    typename, K.transpose())
            matfile.write('\n')

        for i, C in enumerate(Cs):
            print_vector_mzd(matfile, 'C_{}_{}_{}_{}'.format(inst.n, inst.k, inst.r, i),
                    typename, C)
            matfile.write('\n')

        print_matrix_mzd(matfile, 'precomputed_round_key_matrix_linear_part_{}_{}_{}'.format(inst.n, inst.k, inst.r),
                    'mzd_local_t', precomputed_values[0])
        matfile.write('\n')

        print_matrix_mzd(matfile, 'precomputed_round_key_matrix_non_linear_part_{}_{}_{}'.format(inst.n, inst.k, inst.r),
                    'mzd_local_t', precomputed_matrix)
        matfile.write('\n')

        matfile.write(
'''const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_linear_layer(uint32_t r) {{
  switch(r) {{
    default:
      return NULL;'''.format(inst=inst, s=inst.n / 64))

        for i in range(inst.r):
            matfile.write('''
    case {i}:
      return &L_{inst.n}_{inst.k}_{inst.r}_{i};'''.format(i=i, inst=inst))

        matfile.write('\n  }\n}\n\n')

        matfile.write(
'''const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_round_key(uint32_t r) {{
  switch(r) {{
    default:
      return NULL;'''.format(inst=inst, s=inst.n / 64))

        for i in range(inst.r + 1):
            matfile.write('''
    case {i}:
      return &K_{inst.n}_{inst.k}_{inst.r}_{i};'''.format(i=i, inst=inst))

        matfile.write('\n  }\n}\n\n')

        matfile.write(
'''const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_round_const(uint32_t r) {{
  switch(r) {{
    default:
      return NULL;'''.format(inst=inst, s=inst.n / 64))

        for i in range(inst.r):
            matfile.write('''
    case {i}:
      return &C_{inst.n}_{inst.k}_{inst.r}_{i};'''.format(i=i, inst=inst))

        matfile.write('\n  }\n}\n')

        matfile.write(
'''const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_precomputed_round_key_matrix_non_linear_part() {{
  '''.format(inst=inst, s=inst.n / 64))
        matfile.write('return &precomputed_round_key_matrix_non_linear_part_{inst.n}_{inst.k}_{inst.r};'.format(inst=inst))
        matfile.write('\n  }\n')


        matfile.write(
'''const mzd_local_t* lowmc_{inst.n}_{inst.k}_{inst.r}_get_precomputed_round_key_matrix_linear_part() {{
  '''.format(inst=inst, s=inst.n / 64))
        matfile.write('return &precomputed_round_key_matrix_linear_part_{inst.n}_{inst.k}_{inst.r};'.format(inst=inst))
        matfile.write('\n  }\n')


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 4:
        blocksize, keysize, rounds = map(int, sys.argv[1:4])
        main(blocksize, keysize, rounds)
    else:
        main()
