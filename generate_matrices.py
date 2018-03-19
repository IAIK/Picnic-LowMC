
import pickle
import io


class Instance(object):
    def __init__(self, n, k, r, L, K, R):
        self.n = n
        self.k = k
        self.r = r
        self.L = L
        self.K = K
        self.R = R


def main(blocksize=256, keysize=256, rounds=19):
    ''' Use the global parameters `blocksize`, `keysize` and `rounds`
        to create the set of matrices and constants for the corresponding
        LowMC instance. Save those in a file named
        `matrices_and_constants.dat`.
    '''
    gen = grain_ssg()
    linlayers = []
    for _ in range(rounds):
        linlayers.append(instantiate_matrix(blocksize, blocksize, gen))

    round_constants = []
    for _ in range(rounds):
        constant = [next(gen) for _ in range(blocksize)]
        round_constants.append(constant)

    roundkey_matrices = []
    for _ in range(rounds + 1):
        mat = instantiate_matrix(blocksize, keysize, gen)
        roundkey_matrices.append(mat)

    with open('matrices_and_constants_{}_{}_{}.dat'.format(blocksize, keysize, rounds), 'w') as matfile:
        s = 'LowMC matrices and constants\n'\
            '============================\n'\
            'Block size: ' + str(blocksize) + '\n'\
            'Key size: ' + str(keysize) + '\n'\
            'Rounds: ' + str(rounds) + '\n\n'
        matfile.write(s)
        s = 'Linear layer matrices\n'\
            '---------------------'
        matfile.write(s)
        for r in range(rounds):
            s = '\nLinear layer ' + str(r + 1) + ':\n'
            for row in linlayers[r]:
                s += str(row) + '\n'
            matfile.write(s)

        s = '\nRound constants\n'\
              '---------------------'
        matfile.write(s)
        for r in range(rounds):
            s = '\nRound constant ' + str(r + 1) + ':\n'
            s += str(round_constants[r]) + '\n'
            matfile.write(s)

        s = '\nRound key matrices\n'\
              '---------------------'
        matfile.write(s)
        for r in range(rounds + 1):
            s = '\nRound key matrix ' + str(r) + ':\n'
            for row in roundkey_matrices[r]:
                s += str(row) + '\n'
            matfile.write(s)

    inst = Instance(blocksize, keysize, rounds, linlayers, roundkey_matrices,
            round_constants)
    with io.open('matrices_and_constants_{}_{}_{}.pickle'.format(blocksize, keysize, rounds), 'wb') as matfile:
        pickle.dump(inst, matfile, protocol=2)


def instantiate_matrix(n, m, gen):
    ''' Instantiate a matrix of maximal rank using bits from the
        generatator `gen`.
    '''
    while True:
        mat = []
        for _ in range(n):
            row = []
            for _ in range(m):
                row.append(next(gen))
            mat.append(row)
        if rank(mat) >= min(n, m):
            return mat

def rank(matrix):
    ''' Determine the rank of a binary matrix. '''
    # Copy matrix
    mat = [[x for x in row] for row in matrix]
    
    n = len(matrix)
    m = len(matrix[0])
    for c in range(m):
        if c > n - 1:
            return n
        r = c
        while mat[r][c] != 1:
            r += 1
            if r >= n:
                return c
        mat[c], mat[r] = mat[r], mat[c]
        for r in range(c + 1, n):
            if mat[r][c] == 1:
                for j in range(m):
                    mat[r][j] ^= mat[c][j]
    return m


def grain_ssg():
    ''' A generator for using the Grain LSFR in a self-shrinking generator. '''
    state = [1 for _ in range(80)]
    index = 0
    # Discard first 160 bits
    for _ in range(160):
        state[index] ^= state[(index + 13) % 80] ^ state[(index + 23) % 80]\
                        ^ state[(index + 38) % 80] ^ state[(index + 51) % 80]\
                        ^ state[(index + 62) % 80]
        index += 1
        index %= 80
    choice = False
    while True:
        state[index] ^= state[(index + 13) % 80] ^ state[(index + 23) % 80]\
                        ^ state[(index + 38) % 80] ^ state[(index + 51) % 80]\
                        ^ state[(index + 62) % 80]
        choice = state[index]
        index += 1
        index %= 80
        state[index] ^= state[(index + 13) % 80] ^ state[(index + 23) % 80]\
                        ^ state[(index + 38) % 80] ^ state[(index + 51) % 80]\
                        ^ state[(index + 62) % 80]
        if choice == 1:
            yield state[index]
        index += 1
        index %= 80


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 4:
        blocksize, keysize, rounds = map(int, sys.argv[1:4])
        main(blocksize, keysize, rounds)
    else:
        main()
