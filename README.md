Picnic: Post-Quantum Signatures -- LowMC constant generation
============================================================

The [Picnic](https://microsoft.github.io/Picnic/) is a familiy of digital signature schemes secure
against attacks by quantum computers. This repository contains the code to generate the LowMC
constants as used by the [optimized implementation](https://github.com/IAIK/Picnic). The generated
constantes are compatible with the reduced round key computation and optimized linear layer
evaluation as discussed in the [paper](https://eprint.iacr.org/2018/772) introducing an optimal
representation of LowMC.

The LowMC constants are generated using the [script](https://github.com/LowMC/lowmc) from the LowMC
authors.

Generating constants
--------------------

After installing [SageMath](https://www.sagemath.org/), simply run `make`:
```sh
make
```

License
-------

The code is licensed under the MIT license.
