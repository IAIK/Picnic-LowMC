from sage.all import matrix, vector, GF, random_matrix, random_vector, copy, timeit

# from process_matrices import print_mzd
from generate_matrices import Instance
from process_matrices import print_vector
import pickle, io

F = GF(2)

def state_to_byte_array(s):
    a = []
    for i in range(0, len(s), 8):
        tmp = 0
        for j in range(8):
            tmp = tmp | int(s[i+j]) << j
        a.insert(0, tmp)
    return a

def print_hex_state(s):
    print "".join([hex(x)[2:].zfill(2).upper() for x in state_to_byte_array(s)])


def S(a, b, c):
    return (a + b*c, a + b + a*c, a + b + c + a*b)

def print_mzd(w, do_reverse=False):
    if do_reverse:
        print("".join(str(x) for x in reversed(w)))
    else:
        print("".join(str(x) for x in w))

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
        assert idx < mat.nrows(), "mat should have full rank, something is wrong"
        M_dot[current_rank,:] = mat[idx,:]
        if current_rank == M_dot.rank():
            indices.append(idx)
        current_rank = M_dot.rank()
        idx += 1

    return M_dot, indices



class LowMC(object):
    def __init__(self, n, k, s, r):
        self.n = n
        self.k = k
        self.s = s
        self.r = r

        with io.open('matrices_and_constants_{}_{}_{}.pickle'.format(n, k, r), 'rb') as matfile:
            inst = pickle.load(matfile)

        if inst.n != n or inst.k != k or inst.r != r:
            raise ValueError("Unexpected LowMC instance.")

        P = matrix(F, n, n)
        for i in range(n):
            P[n - i - 1, i] = 1

        self.L = [P * matrix(F, L) * P for L in inst.L]
        self.K = [P * matrix(F, K) * P for K in inst.K]
        self.C = [P * vector(F, R) for R in inst.R]

        self.Lt = [L.transpose() for L in self.L]
        self.Kt = [K.transpose() for K in self.K]

        self.rrkc_precomputations_done = False
        self.rll_precomputations_done = False

    def rrkc_precomputations(self):
        self.Li = [m.inverse() for m in self.Lt]
        self.LiK = [self.Kt[i + 1] * self.Li[i] for i in range(self.r)]
        self.LiC = [self.C[i] * self.Li[i] for i in range(self.r)]

        mod_Li = [copy(self.Li[i]) for i in range(self.r)]
        for j in range(self.r):
            mod_Li[j][self.n - 3*self.s:, :self.n] = matrix(F, 3 * self.s, self.n)

        self.precomputed_key_matrix = None
        self.precomputed_key_matrix_nl = matrix(F, self.n, (self.s * 3) * self.r)
        self.precomputed_constant = None
        self.precomputed_constant_nl = vector(F, (self.s * 3) * self.r)

        for round in range(self.r):
            tmp = copy(self.LiK[round])
            tmpC = copy(self.LiC[round])

            for i in range(round + 1, self.r):
                x = self.LiK[i]
                c = self.LiC[i]
                for j in range(i - 1, round - 1, -1):
                    x = x * mod_Li[j]
                    c = c * mod_Li[j]
                tmp += x
                tmpC += c

            # non-linear part
            idx = round * (3 * self.s)
            self.precomputed_key_matrix_nl[:self.n, idx:idx + 3 * self.s] = tmp[:self.n, self.n - 3*self.s:]
            self.precomputed_constant_nl[idx:idx + 3 * self.s] = tmpC[self.n - 3 * self.s:]

            # linear part
            if round == 0:
                tmp[:,self.n - 3*self.s:] = matrix(F, self.n, 3*self.s)
                tmpC[self.n - 3*self.s:] = vector(F, 3*self.s)
                self.precomputed_key_matrix = tmp
                self.precomputed_constant = tmpC

        self.rrkc_precomputations_done = True

    def rll_precomputations(self):
        self.R = []
        self.R_dot = []
        self.R_dot_inv = []
        self.R_wedge = []
        self.R_cols = []
        self.T_vee = []
        
        # i = 1
        self.R.append(self.Lt[0][:, 0:self.n-3*self.s])
        Rdot, colR = calc_dot_from_full_rank(self.R[0])
        self.R_dot.append(Rdot)
        self.R_cols.append(colR)
        self.R_dot_inv.append(self.R_dot[0].inverse())

        self.R_wedge.append(self.R[0]*self.R_dot_inv[0])

        # i = 2...r-1
        for i in range(1,self.r-1):
            self.T_vee.append(self.R_dot[i-1] * self.Lt[i][0:self.n-self.s*3,:])
            R = matrix(F, self.n, self.n-self.s*3)
            R[0:self.n-3*self.s, :] = self.T_vee[-1][0:self.n-3*self.s, 0:self.n-3*self.s]
            R[self.n-3*self.s:self.n, :] = self.Lt[i][self.n-3*self.s:self.n, 0:self.n-3*self.s]
            self.R.append(R)
            Rdot, colR = calc_dot_from_full_rank(self.R[i])
            self.R_dot.append(Rdot)
            self.R_cols.append(colR)
            self.R_dot_inv.append(self.R_dot[i].inverse())
            self.R_wedge.append(self.R[i]*self.R_dot_inv[i])

        # i = r
        self.T_vee.append(self.R_dot[self.r-2] * self.Lt[self.r-1][0:self.n-self.s*3,:])

        self.rll_precomputations_done = True


    def keygen(self):
        return random_vector(F, self.k)

    def S(self, s):
        for i in xrange(self.s):
            t = s[self.n - (3*(i+1)) : self.n - 3*(i)]
            s[self.n - (3*(i+1)) : self.n - 3*(i)] = S(t[0], t[1], t[2])
        return s

    def S_nl(self, s):
        for i in xrange(self.s):
            t = s[self.s*3 - (3*(i+1)) : self.s*3 - 3*(i)]
            s[self.s*3 - (3*(i+1)) : self.s*3 - 3*(i)] = S(t[0], t[1], t[2])
        return s

    def enc(self, sk, p):
        s = self.K[0] * sk + p
        for i in xrange(self.r):
            s = self.S(s)
            s = self.L[i] * s
            s = self.C[i] + s
            s = self.K[i + 1] * sk + s
        return s

    def enc_transposed(self, sk, p):
        s = sk * self.Kt[0] + p
        for i in xrange(self.r):
            s = self.S(s)
            s = s * self.Lt[i]
            s = self.C[i] + s
            s = sk * self.Kt[i + 1] + s
        return s

    def enc_rrkc(self, sk, p):
        
        #precomputations
        if not self.rrkc_precomputations_done:
            self.rrkc_precomputations()

        #---------------------------
        #ENCRYPTION
        #---------------------------
        v = sk * self.precomputed_key_matrix_nl + self.precomputed_constant_nl
        s = p + sk * (self.Kt[0] + self.precomputed_key_matrix) + self.precomputed_constant

        #calculate non-linear part
        for i in range(self.r):
            s = self.S(s)
            s[self.n - 3*self.s:] += v[i * (3 * self.s): (i+1) * (3 * self.s)]
            s = s * self.Lt[i]

        return s

    def enc_rrkc_rll(self, sk, p):
        
        #precomputations
        if not self.rrkc_precomputations_done:
            self.rrkc_precomputations()

        if not self.rll_precomputations_done:
            self.rll_precomputations()

        #---------------------------
        #ENCRYPTION
        #---------------------------
        v = sk * self.precomputed_key_matrix_nl + self.precomputed_constant_nl
        x = p + self.precomputed_constant
        x = x + sk * (self.Kt[0] + self.precomputed_key_matrix)

        # round 1
        y = self.S(x)
        y[self.n - 3*self.s:] += v[0 * (3 * self.s): (0+1) * (3 * self.s)]
        x_0 = y * self.Lt[0][:,self.n-3*self.s:self.n]
        z_1 = y * self.R_wedge[0]
        # rounds 2..r-1
        for i in range(1,self.r-1):
            # print x_0
            y_0 = self.S_nl(x_0)
            y_0 += v[i * (3 * self.s): (i+1) * (3 * self.s)]
            x_0 = y_0 * self.Lt[i][self.n-self.s*3:self.n, self.n-self.s*3:self.n] # L00(y_0)
            x_0 = x_0 + z_1 * self.T_vee[i-1][0:self.n-3*self.s, self.n-3*self.s:self.n] 
            z_1 = vector(F, z_1.list() + y_0.list()) * self.R_wedge[i]
        # round r
        y_0 = self.S_nl(x_0)
        y_0 += v[(self.r-1) * (3 * self.s): (self.r) * (3 * self.s)]
        x_0 = y_0 * self.Lt[-1][self.n-self.s*3:self.n, :]
        x_0 = x_0 + z_1 * self.T_vee[-1]
        
        return x_0


def test_1():
    k = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    e = "1010001100111111101000101100011100101000011100101100111101110110111001111011101010001001100010110111011010100010010110100110001001001000011010100111000000010001000001111100000111111001000101100100011111110011111011000110100001110000001111001010110111001000"

    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in k])
    e = vector(F, [int(x) for x in e])

    lowmc = LowMC(256, 256, 10, 38)
    assert lowmc.enc(k, p) == e
    assert lowmc.enc_transposed(k, p) == e
    assert lowmc.enc_rrkc(k, p) == e

def test_2():
    k = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001010101111111111"
    e = "01001001010000111101010111100110110001110010110011111101000100110001011100000010011000010000110101001011000100010101010011001101"

    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in p])
    e = vector(F, [int(x) for x in e])

    lowmc = LowMC(128, 128, 10, 20)
    assert lowmc.enc(k, p) == e
    assert lowmc.enc_transposed(k, p) == e
    assert lowmc.enc_rrkc(k, p) == e

def test_3():
    k = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001111111111010101"
    e = "00001110111100011011000111000100001100010011100011101110111001010100001110101011001001101111100111010000010011100000110001110000"
    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in p])
    e = vector(F, [int(x) for x in e])

    lowmc = LowMC(128, 128, 10, 20)
    assert lowmc.enc(k, p) == e
    assert lowmc.enc_transposed(k, p) == e
    assert lowmc.enc_rrkc(k, p) == e
    assert lowmc.enc_rrkc_rll(k, p) == e

def test_4():
    k = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001111111111010101"
    e = "01000101110101111100111101010011110100110110000100010000110001001000001001111001100000000000011000101110010110110001110100111010"
    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in p])
    e = vector(F, [int(x) for x in e])

    lowmc = LowMC(128, 128, 1, 182)
    assert lowmc.enc(k, p) == e
    assert lowmc.enc_transposed(k, p) == e
    assert lowmc.enc_rrkc(k, p) == e
    assert lowmc.enc_rrkc_rll(k, p) == e

def test_5():
    k = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001111111111010101"
    e = "01000101001111000001100001111101001001010001001111101000001101110110110100000111010100100011011100000111010100001011001000000000"

    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in p])
    e = vector(F, [int(x) for x in e])

    # lowmc = LowMC(128, 128, 1, 192)
    # lowmc.rrkc_precomputations()
    # lowmc.rll_precomputations()
    # with io.open('tmplowmc', 'wb') as matfile:
        # lowmc = pickle.dump(lowmc, matfile)
    with io.open('tmplowmc', 'rb') as matfile:
        lowmc = pickle.load(matfile)
    # exit(1)
    # assert lowmc.enc(k, p) == e
    # assert lowmc.enc_transposed(k, p) == e
    # assert lowmc.enc_rrkc(k, p) == e
    assert lowmc.enc_rrkc_rll(k, p) == e
# test_1()
# test_2()
# test_3()
test_4()
# test_5()
