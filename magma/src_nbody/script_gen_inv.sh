#!/bin/bash

# Script generating files containing primary and secondary invariants for N-body terms up to a given DEGREE

# # -------------------------------------------
# # Paramters
# # -------------------------------------------
# NBODY=5
# DEGREE=9 #maximal polynomial DEGREE
#
# PREFSEC="SEC" #prefix for the secondary invariants
# PREFIRRSEC="IS" #prefix for the irreducible secondary invariants

NBlengths=$((($NBODY*($NBODY-1))/2))

# printing the parameters
ECHO NBODY order= $NBODY
ECHO Nb of lengths= $NBlengths
ECHO Polynomial DEGREE= $DEGREE

#Define the file names used later and print out their names
mkdir -p "$OUTPUT_DIR"NB_$NBODY""_deg_$DEGREE

filename_log=$OUTPUT_DIR"NB_$NBODY""_deg_$DEGREE/NB_$NBODY""_deg_$DEGREE""_log.txt"
fn_jl_check=$OUTPUT_DIR"NB_$NBODY""_deg_$DEGREE/NB_$NBODY""_deg_$DEGREE""_non_efficient_invariants.jl"
fn_jl_irr_inv=$OUTPUT_DIR"NB_$NBODY""_deg_$DEGREE/NB_$NBODY""_deg_$DEGREE""_irr_invariants.jl"
fn_jl_prim_inv=$OUTPUT_DIR"NB_$NBODY""_deg_$DEGREE/NB_$NBODY""_deg_$DEGREE""_prim_invariants.jl"
fn_jl_sec_rel_inv=$OUTPUT_DIR"NB_$NBODY""_deg_$DEGREE/NB_$NBODY""_deg_$DEGREE""_relations_invariants.jl"

ECHO Output files:

ECHO $filename_log
ECHO $fn_jl_check
ECHO $fn_jl_irr_inv
ECHO $fn_jl_prim_inv
ECHO $fn_jl_sec_rel_inv


# -------------------------------------------
# Magma part
# -------------------------------------------
if [ $MAGMA_RUN -eq 1 ]
then
#put the paramters into the input file
cp magma_gen_inv_NBODY.m NBODY_run.m;

sed -i -e "s/DEGREE/$DEGREE/g" NBODY_run.m;
sed -i -e "s/NBODY/$NBODY/g" NBODY_run.m;

#connect to galois and copy the input files
scp pack_opt_primaries.m $SSH_ADDRESS: ;
scp NBODY_run.m $SSH_ADDRESS: ;

#run the magma computation
ssh $SSH_ADDRESS << EOF
magma NBODY_run.m
EOF

#remove now useless file
rm NBODY_run.m;

#copy the output on the local machine
scp $SSH_ADDRESS:logNbody_output.txt .;

# change the name of the output file
mv logNbody_output.txt $filename_log
fi

# Generate julia file with function computing primary and secondary invariants (not efficient but hopefully correct)
# Pick lines with primaries, irreducible secondaries and secondaries
cp $filename_log $fn_jl_check
# remove things that are not secondaries or primaries
sed -i '' '/v\[1\]/,$!d' $fn_jl_check
sed -i '' '/Total/d' $fn_jl_check

# remove line with number of secondaries
TEMPVAR=$(grep "Nb_secondary_invariants" $fn_jl_check)
NBsecondaries=${TEMPVAR#*=}
#print number of secondaries
ECHO "Nb of secondaries="$NBsecondaries

sed -i '' '/ Nb_secondary_invariants/d' $fn_jl_check

# remove line with number of irreducible secondaries
TEMPVAR=$(grep "Nb_irr_sec_invariants" $fn_jl_check)
NBirr_sec=${TEMPVAR#*=}
ECHO "Nb of irreducible secondaries="$NBirr_sec

sed -i '' '/ Nb_irr_sec_invariants/d' $fn_jl_check

sed -i '' '/ Names_of_variables/,/end_names_of_variables/d' $fn_jl_check

awk '/pv$/ { printf("%s\t", $0); next } 1' $fn_jl_check


# replace variables for the secondaries and the primaries
# and d(ik,il) -> x[i]
count=$NBlengths
for a in `seq $(($NBODY-2)) -1 0`; do
	for b in `seq $((($NBODY-1))) -1 $(($a+1))`; do
		ECHO $a,$b, $count
		OLD="d(i$a,i$b)" ;
		NEW="x\[$count\]" ;
		sed -i '' "s/$OLD/$NEW/g" $fn_jl_check
		count=$(($count-1))
	done
done

cp $fn_jl_check $fn_jl_irr_inv
cp $fn_jl_check $fn_jl_prim_inv
cp $fn_jl_check $fn_jl_sec_rel_inv

sed -i '' '/SYM/d' $fn_jl_check
sed -i '' '/SYYM/d' $fn_jl_check

#write the julia file for the check function
echo "v=zeros($NBsecondaries"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "pv=zeros($NBirr_sec"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "prim=zeros($NBlengths"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "function invariants_check(x)" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check


echo "return prim, v, pv"  >> $fn_jl_check
echo "" >> $fn_jl_check
echo "end" >> $fn_jl_check
# echo "x = rand($NBlengths"")"  >> $fn_jl_check
# echo "display(invariants_Q$NBlengths""_check(x))"  >> $fn_jl_check

# Generate a file with only monomials of irreducible secondaries
# ---------------------------------------------------------
sed -i '' '/SYM/!d' $fn_jl_irr_inv
sed -i '' 's/SYM/ /' $fn_jl_irr_inv

# Generate a file with only monomials of primaries
# ---------------------------------------------------------
sed -i '' '/SYYM/!d' $fn_jl_prim_inv
sed -i '' 's/SYYM/ /' $fn_jl_prim_inv

# Generate a file with relations between irreducible and secondary invariants
# ---------------------------------------------------------
sed -i '' '/ v\[/!d' $fn_jl_sec_rel_inv
sed -i '' "s/pv\[/$PREFIRRSEC/g" $fn_jl_sec_rel_inv
sed -i '' "s/v\[/$PREFSEC/g" $fn_jl_sec_rel_inv
sed -i '' "s/\]/ /g" $fn_jl_sec_rel_inv

# Moving lines starting with +
# ---------------------------------------------------------
gsed -i '$!N;s/\n\s*+/ +/;P;D' $fn_jl_check


#remove useless file
rm NBODY_run.m-e
