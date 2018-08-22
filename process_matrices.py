#! /usr/bin/env sage

from __future__ import unicode_literals

from sage.all import GF, matrix, vector, copy
from six.moves import range

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
    if rcols > bound:
        return ((rcols * (width / 8) + 31) & ~31) / (width / 8)
    else:
        return ((rcols * (width / 8) + 15) & ~15) / (width / 8)

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
    LiK = [Kt[i + 1] * Li[i] for i in range(inst.r)]
    LiC = [Cs[i] * Li[i] for i in range(inst.r)]

    mod_Li = [copy(Li[i]) for i in range(inst.r)]
    for j in range(inst.r):
        mod_Li[j][inst.n - 3*sboxes:, :inst.n] = matrix(F, 3 * sboxes, inst.n)

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
            tmp[:,inst.n - 3*sboxes:] = matrix(F, inst.n, 3*sboxes)
            tmpC[inst.n - 3*sboxes:] = vector(F, 3*sboxes)
            precomputed_key_matrix = tmp
            precomputed_constant = tmpC

    with io.open('lowmc_{}_{}_{}.h'.format(inst.n, inst.k, inst.r), 'w') as matfile:
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

        typename = 'mzd_local_t'

        for i, L in enumerate(Ls):
            print_matrix_mzd(matfile, 'L_{}'.format(i),
                    typename, L.transpose())
            matfile.write('\n')

        matfile.write('#if !defined(REDUCED_LINEAR_LAYER)')
        for i, K in enumerate(Ks):
            matfile.write('\n')
            print_matrix_mzd(matfile, 'K_{}'.format(i),
                    typename, K.transpose())

        for i, C in enumerate(Cs):
            matfile.write('\n')
            print_vector_mzd(matfile, 'C_{}'.format(i),
                    typename, C)
        matfile.write('#endif\n')

        matfile.write('#if defined(REDUCED_LINEAR_LAYER)\n')
        print_matrix_mzd(matfile, 'precomputed_round_key_matrix_linear_part',
                    'mzd_local_t', precomputed_key_matrix + Ks[0].transpose())
        matfile.write('\n')

        print_matrix_mzd(matfile, 'precomputed_round_key_matrix_non_linear_part',
                    'mzd_local_t', precomputed_key_matrix_nl)
        matfile.write('\n')

        print_vector_mzd(matfile, 'precomputed_constant_linear_part', typename, precomputed_constant)
        matfile.write('\n')
        print_vector_mzd(matfile, 'precomputed_constant_non_linear_part', typename, precomputed_constant_nl)
        matfile.write('#endif\n\n')

        matfile.write(
'''#if defined(MUL_M4RI)
static lowmc_round_t rounds[{inst.r}] = {{
#else
static const lowmc_round_t rounds[{inst.r}] = {{
#endif
'''.format(inst=inst))
        for i in range(inst.r):
            matfile.write(
'''
  {{
#if !defined(REDUCED_LINEAR_LAYER)
#if defined(MUL_M4RI)
    &K_{j}, &L_{i}, &C_{i}, NULL, NULL
#else
    &K_{j}, &L_{i}, &C_{i}
#endif
#else
#if defined(MUL_M4RI)
    &L_{i}, NULL
#else
    &L_{i}
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
#if defined(REDUCED_LINEAR_LAYER)
  &precomputed_round_key_matrix_linear_part,
#else
  &K_0,
#endif
#if defined(MUL_M4RI)
  NULL,
#endif
  rounds,
#if defined(REDUCED_LINEAR_LAYER)
  &precomputed_round_key_matrix_non_linear_part,
#if defined(MUL_M4RI)
  NULL,
#endif
  &precomputed_constant_linear_part,
  &precomputed_constant_non_linear_part,
#endif
#if defined(WITH_CUSTOM_INSTANCES)
  {{ NULL, NULL, NULL, NULL }},
  false
#endif
}};
'''.format(inst=inst, m=sboxes))

if __name__ == '__main__':
    import sys

    if len(sys.argv) == 4:
        blocksize, keysize, rounds = map(int, sys.argv[1:4])
        main(blocksize, keysize, rounds)
    else:
        main()
