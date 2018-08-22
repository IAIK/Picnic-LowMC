Picnic: Post-Quantum Signatures -- LowMC constant generation
============================================================

The [Picnic](https://microsoft.github.io/Picnic/) digital signature schemes
secure against attacks by quantum computers. This repository contains the code
to generate the LowMC constants as used by the [optimized
implementation](https://github.com/IAIK/Picnic).

The LowMC constants are generated using the
[script](https://github.com/LowMC/lowmc) from the LowMC authors.

Generating constants
--------------------

After installing [SageMath](https://www.sagemath.org/), simply run `make`:
```sh
make
```

License
-------

The code is licensed under the MIT license.
