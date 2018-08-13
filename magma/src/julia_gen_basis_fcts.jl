using Combinatorics, StaticArrays, Calculus, NBodyIPs

include("invariants_generator_ba.jl")
include("../src_nbody/inv_monomials.jl")

# # Parameters
GROUP_NAME="BA_5B"
OUTPUT_DIR="../data/"

# polynomial degree
deg = 3;

# degree function
function degree(monomial)
    # return sum(monomial)
    return prod(find(monomial))
end



# name of output files
filename = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_basis_fcts"
preword = " "
prefix = "BF"



# --------------
include(OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_group_elements.jl")
# --------------
# nb of variables
Nbvar = length(G_BA_5B[1])

# Generate monomials
function generate_rep_mon_deg(M,total_deg,partial_deg,degree)
    mon_list_out = SVector{M, Int}[]
    for deg = 1:total_deg
        mon_list_out_deg = SVector{M, Int}[]
        Decomp_deg = multisets(deg,M)
        for i=1:length(Decomp_deg)
            monrep = mon_repr(SVector(Decomp_deg[i]...))
            if !(monrep in mon_list_out_deg)&(degree(monrep)<=partial_deg)
                push!(mon_list_out_deg,monrep)
            end
        end
        append!(mon_list_out,mon_list_out_deg)
    end
    return MonList(mon_list_out)
end

Monomial_list = generate_rep_mon_deg(Nbvar,deg,deg,degree)

 max_exp = generate_invariants(~,filename,preword,prefix,Mon(Monomial_list),GROUP_NAME)


#Nb of basis functions
NBbf = length(Mon(Monomial_list))

# -------------------------------------------
#
# Generate function with all basis functions and derivatives
#
# -------------------------------------------
file = OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_basis_functions.jl";

open(file, "w") do f
    write(f, "module ", GROUP_NAME, " \n\n")
    write(f, "using NBodyIPs.FastPolys \n")
    write(f, "using StaticArrays \n")
    write(f, "using BenchmarkTools: @btime \n\n")
    write(f, "import NBodyIPs.tdegrees \n\n")

    # write the definition for simplex permutation
    sim_perm = open(io->read(io), OUTPUT_DIR*GROUP_NAME*"/"*GROUP_NAME*"_group_elements.jl")
    write(f, sim_perm)

    # write the definition of the constant vectors
    w1 = open(io->read(io), filename*"1.jl")
    # read(filenameprim1)
    write(f, w1)

    # write the definitions of the types
    write(f, "\n")
    w2 = open(io->read(io), filename*"2.jl")
    # read(filenameprim2)
    write(f, w2)

    write(f, "\n\n")

# -------------------------------------------
#
# Generate both invariants and derivatives of the invariants
#
# -------------------------------------------
    write(f, "\n\n\n\n")
    write(f, "function invariants_ed(x1::SVector{$Nbvar, T}) where {T}\n")
    for i=2:max_exp
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end
    write(f, "\n   dx1 = @SVector ones($Nbvar)\n")
    for i=2:max_exp
        im = i-1;
        write(f, "   dx$i = $i * x$im \n")
    end

    # write the basis functions
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Basis functions\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    w5 = open(io->read(io), filename*"5.jl")
    # read(filenameprim5)
    write(f, w5)

    #write the return part
    write(f, "\n\n")
    write(f, "return (@SVector [")
    for i=1:NBbf
        write(f, prefix, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NBbf
        write(f, "d", prefix, "$i,")
    end
    write(f, "])\n end \n\n")
    write(f, "end")
end



#Remove the temporary files
rm(filename*"1.jl");
rm(filename*"2.jl");
rm(filename*"3.jl");
rm(filename*"4.jl");
rm(filename*"5.jl");
