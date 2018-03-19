from sage.all import matrix, vector, GF, random_matrix, random_vector, copy, timeit

# from process_matrices import print_mzd
from generate_matrices import Instance
from process_matrices import print_vector
import pickle, io

n = 256
k = n
s = 10
r = 38

F = GF(2)



def S(a, b, c):
    return (a + b*c, a + b + a*c, a + b + c + a*b)

def print_mzd(w, do_reverse=False):
    if do_reverse:
        print("".join(str(x) for x in reversed(w)))
    else:
        print("".join(str(x) for x in w))


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
        self.LiK = [self.K[i + 1] * self.L[i].inverse() for i in xrange(self.r)]

        self.Lt = [L.transpose() for L in self.L]
        self.Kt = [K.transpose() for K in self.K]


    def keygen(self):
        return random_vector(F, self.k)

    def S(self, s):
        for i in xrange(self.s):
            t = s[self.n - (3*(i+1)) : self.n - 3*(i)]
            s[self.n - (3*(i+1)) : self.n - 3*(i)] = S(t[0], t[1], t[2])
        return s

    def enc(self, sk, p):
        s = self.K[0] * sk + p
        for i in xrange(self.r):
            s = self.S(s)
            s = self.L[i] * s
            s = self.C[i] + s
            s = self.K[i + 1] * sk + s
        return s

    def enc1(self, sk, p):
        s = self.K[0].transpose() * sk + p
        for i in xrange(self.r):
            s = self.S(s)
            s = self.L[i].transpose() * s
            s = self.C[i] + s
            s = self.K[i + 1].transpose() * sk + s
        return s

    def enc2(self, sk, p):
        s = sk * self.Kt[0] + p
        for i in xrange(self.r):
            s = self.S(s)
            s = s * self.Lt[i]
            s = self.C[i] + s
            s = sk * self.Kt[i + 1] + s
        return s

    def enc3(self, sk, p):
        #precomputations
        
        mod_inv = [matrix(self.n, self.n) for i in range(self.r)]
        
        for i in range(0, self.r):
            mod_inv[i] = self.L[i].inverse()
            mod_inv[i][self.n - 3*self.s:, :self.n] = matrix(3*self.s, self.n)

        precomputed_values = [matrix(self.n, self.n) for i in range(self.r)]

        for round in range(self.r):
            tmp = matrix(self.n, self.n)
            tmp += copy(self.LiK[round])

            for i in range(round + 1, self.r):
                x = copy(self.LiK[i])

                for j in range(i - 1, round - 1, -1):
                     x = x * mod_inv[j]

                tmp += x
            precomputed_values[round] = tmp
            
        #---------------------------
        #ENCRYPTION
        #---------------------------

        s = self.K[0].transpose() * sk  + p

        #calculate linear part
        s[:self.n - 3*self.s] +=  (sk * precomputed_values[0])[:self.n - 3*self.s]

        
        #calculate non-linear part
        for round in range(self.r):
            s = self.S(s)
            s[self.n - 3*self.s:] += (sk * precomputed_values[round][:n, self.n - 3*self.s:])
            s = s * self.L[round]
            s = self.C[round] + s


        return s

    def enc4(self, sk, p):
        
        #precomputations
        
        mod_inv = [matrix(self.n, self.n) for i in range(self.r)]
        
        for i in range(0, self.r):
            mod_inv[i] = self.L[i].inverse()
            mod_inv[i][self.n - 3*self.s:, :self.n] = matrix(3*self.s, self.n)

        precomputed_values = [matrix(self.n, self.n) for i in range(self.r)]
        precomputed_matrix = matrix(self.n, (self.s * 3 + 2) * self.r)

        for round in range(self.r):
            tmp = matrix(self.n, self.n)
            tmp += copy(self.LiK[round])

            for i in range(round + 1, self.r):
                x = copy(self.LiK[i])

                for j in range(i - 1, round - 1, -1):
                     x = x * mod_inv[j]

                tmp += x
            precomputed_values[round] = tmp
            precomputed_matrix[:n, round * (3 * self.s + 2): round * (3 * self.s + 2) + 30] = tmp[:n, self.n - 3*self.s:]
        precomputed_matrix = sk * precomputed_matrix
        #---------------------------
        #ENCRYPTION
        #---------------------------

        s = self.K[0].transpose() * sk  + p

        #calculate linear part
        s[:self.n - 3*self.s] +=  (sk * precomputed_values[0])[:self.n - 3*self.s]

        #calculate non-linear part
        for round in range(self.r):
            s = self.S(s)
            s[self.n - 3*self.s:] += precomputed_matrix[round * (3 * self.s + 2): round * (3 * self.s + 2) + 30]

            s = s * self.L[round]
            s = self.C[round] + s

        return s


def test_1():
    k = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    e = "1010001100111111101000101100011100101000011100101100111101110110111001111011101010001001100010110111011010100010010110100110001001001000011010100111000000010001000001111100000111111001000101100100011111110011111011000110100001110000001111001010110111001000"

    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in k])
    e = vector(F, [int(x) for x in e])

    lowmc = LowMC(256, 256, 10, 38)
    assert lowmc.enc(k, p) == e
    assert lowmc.enc2(k, p) == e

def test_2():
    k = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
    p = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001010101111111111"
    e = "01001001010000111101010111100110110001110010110011111101000100110001011100000010011000010000110101001011000100010101010011001101"

    k = vector(F, [int(x) for x in k])
    p = vector(F, [int(x) for x in p])
    e = vector(F, [int(x) for x in e])

    lowmc = LowMC(128, 128, 10, 20)
    assert lowmc.enc(k, p) == e
    assert lowmc.enc2(k, p) == e

test_1()
test_2()


lowmc = LowMC(n, k, s, r)


sk = vector(F, [1] + [0] * (k-1))
pt = vector(F, [1] + [0] * (n-1))

lowmc.enc4(sk, pt)
print(lowmc.enc1(sk, pt) == lowmc.enc4(sk, pt))
