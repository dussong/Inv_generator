# Inv_generator

This package generates invariants polynomials (primary and secondary) for permutation groups and generates Julia code for an efficient implementation of these polynomials.


To compute the invariant polynomials, you need an access to the computer algebra system Magma. To run the

The folder `src_nbody` contains the 


1 - run generate_invariants.sh with wanted Nbody and degree parameters.
This creates files in data

2- NOT NEEDED ANYMORE
2 - copy the files
a - NB_?_deg_?_irr_invariants.jl
b - NB_?_deg_?_prim_invariants.jl
c - NB_?_deg_?_non_efficient_invariants.jl
d - NB_?_deg_?_relations_invariants.jl
into the data folder

3 NOT NEEDED ANYMORE
3 - do the following manipulations on these files :
In NB_?_deg_?_non_efficient_invariants.jl
a - check that no line starts with a + (otherwise wring computation of the invariants - check in particular pv[17] and pv[18]), and if there is, move the + to the previous line.

4 - run src_jl/generate_text_invariants.jl.
This generates the invariants in:
NB_?_deg_?_invariants.jl

5 - NOT NEEDED ANYMORE
5 - copy-paste the definition of the invariants into the corresponding file in test and use the relations between irreducible secondaries and secondaries to compute all of them.


6 - you can run test_invariants.jl to test them. For the moment, 4 primary invariants are not right.
