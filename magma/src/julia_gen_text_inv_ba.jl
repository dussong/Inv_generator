using Combinatorics, StaticArrays, Calculus

include("invariants_generator_ba.jl")


# Parameters (given in the parameters.conf file)
GROUP_NAME="BA_5B"
prefsec="SEC" #prefix for the secondaries
prefirrsec="IS" #prefix for the irreducible secondaries
prefprim="P" #prefix for the primaries
OUTPUT_DIR="../data/"

#TODO: include degrees of invariants
# --------------

include(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_invariants.jl")
include(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_invariants.jl")
include(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_group_elements.jl")
# -------------------------------------------
#
# Generate irreducible secondaries
#
# -------------------------------------------
filenameirrsec = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text";

# -------------------------------------------
#
# Define data files
#
# -------------------------------------------
filenameirrsecdata = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_invariants.jl";
preword = "# Irreducible secondaries for group "*GROUP_NAME*"\n"

max_exp_irrsec = generate_invariants(filenameirrsecdata,filenameirrsec,preword,prefirrsec,pv,GROUP_NAME)

#Nb of irreducible secondary invariants
NBirrsec = length(pv)
# -------------------------------------------
#
# Generate primary invariants files
#
# -------------------------------------------
filenameprim = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_text";
filenameprimdata = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_invariants.jl";
preword = "# Primary invariants for "*GROUP_NAME*"  \n"

max_exp_prim = generate_invariants(filenameprimdata,filenameprim,preword,prefprim,prim,GROUP_NAME)

#Nb of primary invariants
NBprim = length(prim)
# -------------------------------------------
#
# Secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_relations_invariants.jl";
# -------------------------------------------
#
# Derivatives of secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec_d = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_relations_invariants_derivatives.jl";

open(filenamesec_d, "w") do f
end

fileI = open(filenamesec)
line = readlines(fileI)
NBsec = length(line)
part1,part2 = split(line[1], "=")
repl1 = replace(part1, prefsec, "d"*prefsec)
open(filenamesec_d, "a") do f
    write(f, repl1, " = @SVector zeros($NBprim) \n")
end

# Construct an expression containing all irreducible secondaries as variables
Variables = Symbol[];
for k=1:NBirrsec
    push!(Variables, parse(prefirrsec*"$k"))
end

for i=2:NBsec
    part1,part2 = split(line[i], " = ")
    repl1 = replace(part1, prefsec, "d"*prefsec)
    open(filenamesec_d, "a") do f
        write(f, repl1, " = ")
    end

    ex2 = parse(part2)
    if typeof(ex2) == Symbol
        open(filenamesec_d, "a") do f
            write(f, "d"*part2)
        end
    else
        Der = differentiate(ex2,Variables)
        Ind_nonzero_der = find(Der)
        for (k,ind) in enumerate(Ind_nonzero_der)
            open(filenamesec_d, "a") do f
                write(f, " + d",prefirrsec, "$ind*", "$(Der[ind])")
            end
        end
    end
    open(filenamesec_d, "a") do f
        write(f, "\n")
    end
end

# -------------------------------------------
#
# Generate function with all invariants
#
# -------------------------------------------
file = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_invariants.jl";

open(file, "w") do f
    write(f, "module ", GROUP_NAME, " \n\n")
    write(f, "using NBodyIPs.FastPolys \n")
    write(f, "using StaticArrays \n")
    write(f, "using BenchmarkTools: @btime \n\n")
    write(f, "import NBodyIPs.tdegrees \n\n")

    # #write the definition of the degrees of the invariants
    # degree_inv = open(io->read(io), filename_deg)
    # write(f, degree_inv)
    # write(f, "\n\n")

    # write the definition for simplex permutation
    sim_perm = open(io->read(io), OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_group_elements.jl")
    write(f, sim_perm)

    # write the definition of the constant vectors
    prim1 = open(io->read(io), filenameprim*"1.jl")
    # read(filenameprim1)
    write(f, prim1)
    irrsec1 = open(io->read(io), filenameirrsec*"1.jl")
    # read(filenameirrsec1)
    write(f, irrsec1)

    # write the definitions of the types
    write(f, "\n")
    prim2 = open(io->read(io), filenameprim*"2.jl")
    # read(filenameprim2)
    write(f, prim2)
    irrsec2 = open(io->read(io), filenameirrsec*"2.jl")
    # read(filenameirrsec2)
    write(f, irrsec2)

    write(f, "\n\n")

    # write the name of the function
    write(f, "function invariants(x1::SVector{$NBprim, T}) where {T}\n")

    # write the precomputed powers of x
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end

    # write the primary invariants
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Primaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    prim3 = open(io->read(io), filenameprim*"3.jl")
    # read(filenameprim3)
    write(f, prim3)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec3 = open(io->read(io), filenameirrsec*"3.jl")
    # read(filenameirrsec3)
    write(f, irrsec3)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec)
    # read(filenamesec)
    write(f, sec)

    #write the return part
    write(f, "\n\n")
    write(f, "return (@SVector [")
    for i=1:NBprim
        write(f, prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NBsec
        write(f, prefsec, "$i,")
    end
    write(f, "])\n end")



    # -------------------------------------------
    #
    # Generate derivatives of the invariants
    #
    # -------------------------------------------
    write(f, "\n\n\n\n")
    write(f, "function invariants_d(x1::SVector{$NBprim, T}) where {T}\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end
    write(f, "\n   dx1 = @SVector ones($NBprim)\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   dx$i = $i * x$im \n")
    end

    # write the primary invariants
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Primaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    prim4 = open(io->read(io), filenameprim*"4.jl")
    # read(filenameprim4)
    write(f, prim4)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec3 = open(io->read(io), filenameirrsec*"3.jl")
    # read(filenameirrsec3)
    write(f, irrsec3)

    write(f, "\n\n")
    irrsec4 = open(io->read(io), filenameirrsec*"4.jl")
    # read(filenameirrsec4)
    write(f, irrsec4)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec_d)
    # read(filenamesec_d)
    write(f, sec)

    #write the return part
    write(f, "\n\n")
    write(f, "return (")
    for i=1:NBprim
        write(f, "d", prefprim, "$i,")
    end
    write(f, "), (")
    for i=1:NBsec
        write(f, "d", prefsec, "$i,")
    end
    write(f, ")\n end")

# -------------------------------------------
#
# Generate both invariants and derivatives of the invariants
#
# -------------------------------------------
    write(f, "\n\n\n\n")
    write(f, "function invariants_ed(x1::SVector{$NBprim, T}) where {T}\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end
    write(f, "\n   dx1 = @SVector ones($NBprim)\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   dx$i = $i * x$im \n")
    end

    # write the primary invariants
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Primaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    prim5 = open(io->read(io), filenameprim*"5.jl")
    # read(filenameprim5)
    write(f, prim5)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec5 = open(io->read(io), filenameirrsec*"5.jl")
    # read(filenameirrsec5)
    write(f, irrsec5)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec)
    # read(filenamesec)
    write(f, sec)

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec_d)
    # read(filenamesec_d)
    write(f, sec)

    #write the return part
    write(f, "\n\n")
    write(f, "return (@SVector [")
    for i=1:NBprim
        write(f, prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NBsec
        write(f, prefsec, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NBprim
        write(f, "d", prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NBsec
        write(f, "d", prefsec, "$i,")
    end
    write(f, "])\n end \n\n")
    write(f, "end")
end



#Remove the temporary files
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text1.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text2.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text3.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text4.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text5.jl");

rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_text1.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_text2.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_text3.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_text4.jl");
rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_prim_text5.jl");
# rm(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_text_deg.jl");
