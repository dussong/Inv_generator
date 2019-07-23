// Running the magma calculation for the computation of the invariants of a given group.


//Logfile name
logfile := "logNbody_output.txt";

//Attach package from Braams
// Attach("pack_opt_primaries.m");

//Prescribe logfile
SetLogFile(logfile: Overwrite := true);

//Define the field
K := RationalField();
Kb := RationalField();
//Define the symmetry group (line that can vary)
// G:= GROUP_DEF;
// Gb := GROUP_DEF;
// Examples:
// G := PermutationGroup<9| (3,4)(6,7)(8,9) >; //4B
// Gb := PermutationGroup<9| (3,4)(6,7)(8,9) >;
G := PermutationGroup<12| (3,4,5)(7,8,9)(10,11,12) >; //4B
Gb := PermutationGroup<12| (3,4,5)(7,8,9)(10,11,12) >;
// G := PermutationGroup<15| (3,4)(8,9)(12,13), (3,4,5,6)(8,9,10,11)(12,13,14,15) >; //4B
// Gb := PermutationGroup<15| (3,4)(8,9)(12,13), (3,4,5,6)(8,9,10,11)(12,13,14,15) >;
// G := PermutationGroup<9| (1,2,3)(4,5,6)(7,8,9)>; //4B
// Gb := PermutationGroup<9| (1,2,3)(4,5,6)(7,8,9)>;
// G := PermutationGroup<6| (1,2)(3,4)(5,6)>; //4B
// Gb := PermutationGroup<6| (1,2)(3,4)(5,6)>;
// G := PermutationGroup<10 | (2,3)(5,6)(10,9), (2,3,4)(5,6,7)(8,10,9), (1,3,4,2)(5,6,10,9)(7,8)>; //5B
// Gb := PermutationGroup<10 | (2,3)(5,6)(10,9), (2,3,4)(5,6,7)(8,10,9), (1,3,4,2)(5,6,10,9)(7,8)>; //5B

//Define the invariant ring
Rpre := InvariantRing(G, K);
Rpre;

R := InvariantRing(Gb, Kb);
R;

// Compute the primary invariants
R0pre := PrimaryInvariants(Rpre);
Rpre;
R;

//To get the variables names
Q := PolynomialRing(R);

// Enforce to have only pure symmetrized monomials
// R0 := [&+[(LeadingMonomial(R0pre[i]))^g:g in Group(Rpre)]: i in [1..#R0pre]];
// R0 := [&+[(Q.1)^g:g in Group(Rpre)],&+[(Q.4)^g:g in Group(Rpre)],&+[(Q.1*Q.1)^g:g in Group(Rpre)],&+[(Q.1*Q.4)^g:g in Group(Rpre)],&+[(Q.4*Q.4*Q.4)^g:g in Group(Rpre)],&+[(Q.1*Q.1*Q.4)^g:g in Group(Rpre)]];
// R0 := [&+[(Q.1)^g:g in Group(Rpre)],&+[(Q.5)^g:g in Group(Rpre)],&+[(Q.1*Q.1)^g:g in Group(Rpre)],&+[(Q.5*Q.5)^g:g in Group(Rpre)],&+[(Q.1*Q.5)^g:g in Group(Rpre)],&+[(Q.1*Q.1*Q.1)^g:g in Group(Rpre)],&+[(Q.5*Q.5*Q.5)^g:g in Group(Rpre)],&+[(Q.5*Q.6*Q.7)^g:g in Group(Rpre)],&+[(Q.1*Q.1*Q.1*Q.2)^g:g in Group(Rpre)],&+[(Q.5*Q.6*Q.7*Q.7)^g:g in Group(Rpre)]];
// R0 := [&+[(Q.1)^g:g in Group(Rpre)],
//        &+[(Q.1*Q.1)^g:g in Group(Rpre)],
//        &+[(Q.1*Q.1*Q.1)^g:g in Group(Rpre)],
//        &+[(Q.1*Q.1*Q.1*Q.1)^g:g in Group(Rpre)],
//        &+[(Q.5)^g:g in Group(Rpre)],
//       &+[(Q.5*Q.5)^g:g in Group(Rpre)],
//       &+[(Q.5*Q.5*Q.5)^g:g in Group(Rpre)],
//       &+[(Q.5*Q.5*Q.5*Q.5)^g:g in Group(Rpre)],
//       &+[(Q.9)^g:g in Group(Rpre)],
//      &+[(Q.9*Q.9)^g:g in Group(Rpre)],
//      &+[(Q.9*Q.9*Q.9)^g:g in Group(Rpre)],
//      &+[(Q.9*Q.9*Q.9*Q.9)^g:g in Group(Rpre)]];



// R`PrimaryInvariants := R0;
R0 := R0pre;
// R`PrimaryInvariants := R0;

// Construct a list of the group elements
printf "nb_group_elements_begin\n";
#G;
printf "nb_group_elements_end\n";

GList := [Eltseq(g) : g in G];
printf "group_def_begin\n";
printf "const Gr_elts = [\n";
for i:=1 to #G do
  GList[i];
  printf ",";
end for;
printf "] \n";
printf "group_def_end\n";

// Display the nb of primary invariants
printf "nb_primary_invariants_begin\n";
#R0;
printf "nb_primary_invariants_end\n";

// Display the primary invariants
printf "primary_invariants_begin\n";
for i:=1 to #R0 do
  printf "prim[";
  printf IntegerToString(i);
  printf "]=";
  R0[i];
  printf "\n";
end for;
printf "primary_invariants_end\n";

// Display the monomials of the primary invariants
printf "prim_inv_mon_list_begin\n";
for i:=1 to #R0 do
  printf "prim[";
  printf IntegerToString(i);
  printf "]=";
  Exponents(LeadingMonomial(R0[i]));
  printf "\n";
end for;
printf "prim_inv_mon_list_end\n";

// Display the degrees of the primary invariants
printf "prim_inv_tot_degree_begin\n";
[TotalDegree(R0[i]): i in [1..#R0]];
printf "prim_inv_tot_degree_end\n";

// Computation of irreducible secondary invariants
R2 := IrreducibleSecondaryInvariants(R);

// Nb of irreducible secondaries
printf "nb_irr_sec_invariants_begin\n";
#R2;
printf "nb_irr_sec_invariants_end\n";

// Write irreducible secondaries
printf "irreducible_sec_invariants_begin\n";
for i:=1 to #R2 do
  printf "pv[";
  printf IntegerToString(i);
  printf "]=";
  R2[i];
  printf "\n";
end for;
printf "irreducible_sec_invariants_end\n";

// Write monomials for irreducible secondaries
printf "irrsec_mon_list_begin\n";
for i:=1 to #R2 do
  printf "pv[";
  printf IntegerToString(i);
  printf "]=";
  Exponents(LeadingMonomial(R2[i]));
  printf "\n";
end for;
printf "irrsec_mon_list_end\n";

// Write degrees for irreducible secondaries
printf "irr_sec_degrees_begin\n";
[TotalDegree(R2[i]): i in [1..#R2]];
printf "irr_sec_degrees_end\n";


// Compute and write relations between irreducible secondaries and all secondaries
A, Q := Algebra(R);
A;
printf "inv_relations_begin\n";
for i:=1 to #Q do
  printf "v[";
  printf IntegerToString(i);
  printf "] = ";
  Q[i];
end for;
printf "inv_relations_end\n";

// Compute all secondary invariants for their number and degrees (would be possible to do it differently)
R1 := SecondaryInvariants(R);

printf "nb_secondaries_begin\n";
#Q;
printf "nb_secondaries_end\n";

printf "degrees_secondaries_begin\n";
[TotalDegree(R1[i]): i in [1..#R1]];
printf "degrees_secondaries_end\n";


exit;
