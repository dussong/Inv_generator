# Compute rt degrees for 5B-BA:

# Primary invariants
prim = Vector{Vector{Int64}}(10)
prim[1]=[ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

prim[2]=[ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ]

prim[3]=[ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

prim[4]=[ 0, 0, 0, 0, 2, 0, 0, 0, 0, 0 ]

prim[5]=[ 1, 0, 0, 0, 1, 0, 0, 0, 0, 0 ]

prim[6]=[ 3, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

prim[7]=[ 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 ]

prim[8]=[ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0 ]

prim[9]=[ 0, 0, 0, 0, 4, 0, 0, 0, 0, 0 ]

prim[10]=[ 4, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

prim

rtdegrees_prim = []
for i=1:length(prim)
    push!(rtdegrees_prim,(sum(prim[i][1:4]),sum(prim[i][5:10])))
end

@show rtdegrees_prim

# With the hack
# (1, 0), (0, 1), (2, 0), (0, 2), (1, 1), (3, 0), (0, 3), (0, 3), (0, 4), (4, 4)
# the hack appears in the very last number

# For the secondaries
include("BA_monomials_deg_secondaries.jl")

sec_rt


rtdeg_secondary = []


for i=1:length(sec_rt)
    rdeg = maximum([sum(sec_rt[i][j][1:4]) for j=1:length(sec_rt[i])])
    tdeg = maximum([sum(sec_rt[i][1][5:10]) for j=1:length(sec_rt[i])])
    push!(rtdeg_secondary,(rdeg,tdeg))
end

@show rtdeg_secondary
