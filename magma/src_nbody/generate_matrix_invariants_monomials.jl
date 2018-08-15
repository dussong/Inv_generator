# This code generates the matrix of basis change between the tuples in terms of the invariants and tuples of monomials.
using JLD
using NBodyIPs

include(homedir() * "/Gits/InvariantsGenerator/magma/src_nbody/invariants_generator.jl")
include(homedir() * "/Gits/InvariantsGenerator/magma/src_nbody/inv_monomials.jl")
include(homedir() * "/Gits/InvariantsGenerator/magma/src_nbody/misc.jl")

NBody = 5;
Deg = 6; #maximal total degree
M = Int(NBody*(NBody-1)/2)

# Loading some useful files
filenameirrsecdata = homedir() * "/Gits/InvariantsGenerator/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_irr_invariants.jl";
filenameprimdata = homedir() * "/Gits/InvariantsGenerator/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_prim_invariants.jl";
filenamesec = homedir() * "/Gits/InvariantsGenerator/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_relations_invariants.jl";

prefsec = "SEC" #prefix for the secondaries
prefirrsec = "IS" #prefix for the irreducible secondaries
prefprim = "P" #prefix for the primaries


# -------------------------------------------------------------
#
# Generate monomials with weights: for primaries, irreducible secondaries,
# and secondaries
#
# -------------------------------------------------------------
#

# -------------------------------------------------------------
# for the irreducible secondaries invariants
Mon_irrsec, coef_list_irrsec = generate_inv_mon(filenameirrsecdata,M,Deg)

# generate list of PolyMon for irreducible secondaries
IrrSecMonPol = CPolyMon{M}[]
for i=1:length(Mon_irrsec)
    push!(IrrSecMonPol,CPolyMon(CMonList([mon_repr(SVector(Mon_irrsec[i]...))]),[coef_list_irrsec[i]]))
end
IrrSecMonPol

# -------------------------------------------------------------
# for the primary invariants
Mon_prim, coef_list_prim,deg_prim = generate_inv_mon(filenameprimdata,M,Deg)

# generate list of PolyMon for primaries
PrimMonPol = CPolyMon{M}[]
for i=1:length(Mon_prim)
    push!(PrimMonPol,CPolyMon(CMonList([mon_repr(SVector(Mon_prim[i]...))]),[coef_list_prim[i]]))
end
PrimMonPol

# -------------------------------------------------------------
#  for all the secondary invariants (need to parse the relations)
SecMonPol = CPolyMon{M}[]

fileI = open(filenamesec)
line = readlines(fileI)
# first line contains sec invariant =1, we remove it
push!(SecMonPol,constpoly(M))
for i=2:length(line)
    part1,part2 = split(line[i], "=")
    Part1 = replace(part1, prefsec, "")
    @assert parse(Int64,Part1) == i
    if contains(line[i], "*")
        part2_1,part2_2 = split(part2, "*")
        Part2_1 = replace(part2_1, prefirrsec, "")
        Part2_2 = replace(part2_2, prefirrsec, "")
        int1 = parse(Int64,Part2_1)
        int2 = parse(Int64,Part2_2)

        IrrSec1 = IrrSecMonPol[int1]
        IrrSec2 = IrrSecMonPol[int2]
        push!(SecMonPol,IrrSec1*IrrSec2)
    else
        Part2 = replace(part2, prefirrsec, "")
        int = parse(Int64,Part2)
        push!(SecMonPol,IrrSecMonPol[int])
    end
end
SecMonPol


# -------------------------------------------------------------
#
# Generate tuples in terms of the invariants
#
# -------------------------------------------------------------
InvTup = gen_tuples(NBody,Deg)
@show length(InvTup)


# -------------------------------------------------------------
#
# Compute monomial coefficients for the invariants and their powers
#
# -------------------------------------------------------------
# max power of the invariants necessary to reach the prescribed degree
maxpower = [maximum(InvTup[j][i] for j=1:length(InvTup)) for i=1:M]
maxpower

# precomputation of the necessary powers of primary invariants
powersPrimInv = []
for i=1:M
    powersi = [PrimMonPol[i]]
    for j=1:(maxpower[i]-1)
        push!(powersi,PrimMonPol[i]*powersi[j])
    end
    push!(powersPrimInv,powersi)
end
powersPrimInv

# Computation of all invariant tuples in terms of monomials
InvMonPoly = CPolyMon{M}[]
for i=1:length(InvTup)
    # initialization
    PMonTup = SecMonPol[InvTup[i][end]+1] #initialize with secondary invariant
    for j=1:M
        if InvTup[i][j] > 0
            PMonTup = PMonTup*powersPrimInv[j][InvTup[i][j]] #multiply with powers of primary invariants
            # power(PrimMonPol[j],InvTup[i][j])
        end
    end
    push!(InvMonPoly,PMonTup)
    print("$i ")
end
InvMonPoly


# -------------------------------------------------------------
#
# Compute partial matrices for the basis change and the inverse
# to get the monomials in terms of the invariants tuples
#
# -------------------------------------------------------------
# Monomial basis and corresponding degrees
MonBasis = Mon(generate_rep_mon(M,Deg))
MonBasisDeg = [sum(MonBasis[i]) for i=1:length(MonBasis)]
@assert length(InvTup) == length(MonBasis)

# Degrees of the invariants tuples
InvMonPolyDeg = [sum(Mon(InvMonPoly[i])[1]) for i=1:length(InvMonPoly)]

# Initialize basis change matrix
BasisChange = Matrix{Float64}[]
InvBasisChange = Matrix{Float64}[]
MonBasisInvTup = []

for deg=1:Deg
    # select monomials and invariant tuples with prescribed degree
    IndMonDeg = find( x->(x == deg), MonBasisDeg)
    IndInvDef = find( x->(x == deg), InvMonPolyDeg)

    MonB = MonBasis[IndMonDeg]
    InvB = InvMonPoly[IndInvDef]
    InvTupSelect = InvTup[IndInvDef]
    @assert length(MonB)==length(InvB)

    MBasisChangeDeg = zeros(Float64,length(MonB),length(MonB))

    for j=1:length(InvTupSelect) #for each of the inv tuple
        InvjMon = Mon(InvB[j]) #get monomials of the inv tuple
        InvjCoef = Coef(InvB[j]) #get coef of the inv tuple
        for (i,InvjMoni) in enumerate(InvjMon)
            indi = find(InvjMoni == MonB[k] for k=1:length(MonB))
            @assert length(indi)==1
            MBasisChangeDeg[j,indi[1]] = InvjCoef[i]
        end
    end
    push!(BasisChange,MBasisChangeDeg)
    push!(InvBasisChange,inv(MBasisChangeDeg))

    for i=1:length(MonB)
        nonzeroind = find( x->(abs(x) >1e-10), inv(MBasisChangeDeg)[i,:])
        coeflist = inv(MBasisChangeDeg)[i,nonzeroind]
        monlist = [Mon(InvB[j]) for j in nonzeroind]

        push!(MonBasisInvTup,[monlist,coeflist])
    end
    print(".")

end

MonBasisInvTup

#TODO: Tests to check the implementation
#TODO: use type NBPoly for MonBasisInvTup



# M_basis_change = zeros(Float64,length(MonBasis),length(MonBasis))
# # express all inv_tuples in terms of monomials
#
# for i=1:length(InvTup)
#     InviMon = Mon(InvMonPoly[i])
#     InviCoef = Coef(InvMonPoly[i])
#     for j=1:length(MonBasis)
#         if (MonBasis[j] in InviMon)
#             indj = find([MonBasis[j] == InviMon[k] for k=1:length(InviMon)])
#             @assert length(indj) == 1
#             M_basis_change[i,j] = InviCoef[indj[1]]
#         end
#     end
# end
# save(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_basis_change.jld", "Mbasischange", M_basis_change, "NBody", NBody, "Deg", Deg)
