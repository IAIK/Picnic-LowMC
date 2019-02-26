LOWMC_C = \
	lowmc_256_256_38.c \
	lowmc_192_192_30.c \
	lowmc_128_128_20.c \
	lowmc_256_256_363.c \
	lowmc_192_192_284.c \
	lowmc_128_128_182.c
LOWMC_PICKLE = $(patsubst lowmc_%,matrices_and_constants_%,$(patsubst %.c,%.pickle,$(LOWMC_C)))
PROCESS_MATRICES=process_matrices.py

all: $(LOWMC_C)

$(LOWMC_C): $(PROCESS_MATRICES)
$(LOWMC_PICKLE): generate_matrices.py

lowmc_256_256_38.c: matrices_and_constants_256_256_38.pickle
	sage $(PROCESS_MATRICES) 256 256 38 10

lowmc_192_192_30.c: matrices_and_constants_192_192_30.pickle
	sage $(PROCESS_MATRICES) 192 192 30 10

lowmc_128_128_20.c: matrices_and_constants_128_128_20.pickle
	sage $(PROCESS_MATRICES) 128 128 20 10

lowmc_256_256_363.c: matrices_and_constants_256_256_363.pickle
	sage $(PROCESS_MATRICES) 256 256 363 1

lowmc_192_192_284.c: matrices_and_constants_192_192_284.pickle
	sage $(PROCESS_MATRICES) 192 192 284 1

lowmc_128_128_182.c: matrices_and_constants_128_128_182.pickle
	sage $(PROCESS_MATRICES) 128 128 182 1

matrices_and_constants_256_256_38.pickle:
	sage generate_matrices.py 256 256 38

matrices_and_constants_192_192_30.pickle:
	sage generate_matrices.py 192 192 30

matrices_and_constants_128_128_20.pickle:
	sage generate_matrices.py 128 128 20

matrices_and_constants_256_256_363.pickle:
	sage generate_matrices.py 256 256 363

matrices_and_constants_192_192_284.pickle:
	sage generate_matrices.py 192 192 284

matrices_and_constants_128_128_182.pickle:
	sage generate_matrices.py 128 128 182
