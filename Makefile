all: lowmc_256_256_38.c lowmc_192_192_30.c lowmc_128_128_20.c

lowmc_256_256_38.c lowmc_192_192_30.c lowmc_128_128_20.c: process_matrices.py
matrices_and_constants_256_256_38.pickle matrices_and_constants_192_192_30.pickle matrices_and_constants_128_128_20.pickle: generate_matrices.py

lowmc_256_256_38.c: matrices_and_constants_256_256_38.pickle
	sage process_matrices.py 256 256 38

lowmc_192_192_30.c: matrices_and_constants_192_192_30.pickle
	sage process_matrices.py 192 192 30

lowmc_128_128_20.c: matrices_and_constants_128_128_20.pickle
	sage process_matrices.py 128 128 20

matrices_and_constants_256_256_38.pickle:
	sage generate_matrices.py 256 256 38

matrices_and_constants_192_192_30.pickle:
	sage generate_matrices.py 192 192 30

matrices_and_constants_128_128_20.pickle:
	sage generate_matrices.py 128 128 20
