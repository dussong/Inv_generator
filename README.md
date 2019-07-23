# Invariant generator

This package generates invariants polynomials (primary and secondary) based on invariant theory for permutation groups [see Derksen, H., & Kemper, G. (2015). Computational Invariant Theory. Springer] and generates Julia code for an efficient implementation of these polynomials.

To compute the invariant polynomials, you need an access to the computer algebra system Magma [Bosma, W., Cannon, J., & Playoust, C. (1997). The Magma Algebra System I: The User Language. Journal of Symbolic Computation, 24(3), 235–265].

To compute the invariants, modify the parameter file `parameters.conf` with the right parameters, the most important being the definition of the group.
Then run `main_run.sh`.

The folder `src_nbody` contains the invariant computation for the permutation of particles. `NBody` corresponds to the number of particles.

After the run, the results (invariants and Julia code) should be in the indicated folder (e.g. `data`).
