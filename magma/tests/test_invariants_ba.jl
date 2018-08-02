using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
using BenchmarkTools

# Check that group name and variables ok. Modify accordingly the line 16.
GROUP_NAME="BA_5B"
OUTPUT_DIR="../data/"
const dim=10

include(homedir()*"/.julia/v0.6/NBodyIPs/src/fastpolys.jl")
using FastPolys

include(OUTPUT_DIR*GROUP_NAME*"/"GROUP_NAME*"_non_efficient_invariants.jl")
include(OUTPUT_DIR*GROUP_NAME*"/"GROUP_NAME*"_invariants_hack.jl")
# include(OUTPUT_DIR*GROUP_NAME*"/"GROUP_NAME*"_group_elements.jl")

using BA_5B: invariants, invariants_d, invariants_ed, simplex_permutations

using JuLIP.Potentials: evaluate, evaluate_d

all_invariants(r) = vcat(invariants(r)...)  # [I1; I2]
all_invariants_check(r) = vcat(invariants_check(r)...)  # [I1; I2]
ad_invariants(r) = ForwardDiff.jacobian(all_invariants, r)

println("-------------------------------------------")
println("   Testing implementation of `invariants`")
println("-------------------------------------------")


println("[1] Correctness of gradients")
r = 1.0 + SVector(rand(dim)...)
println("---------------")
println("dim = $dim")
println("---------------")

I = all_invariants(r)
dI1, dI2 = invariants_d(r)
dI = [hcat(dI1...)'; hcat(dI2...)']
dIh = zeros(size(dI))
r0 = Vector(r)
errs = []
for p = 2:9
   h = .1^p
   dIh = zeros(size(dI))
   for j = 1:length(r)
      r0[j] += h
      Ih = all_invariants(SVector(r0...))
      dIh[:, j] = (Ih - I) / h
      r0[j] -= h
   end
   push!(errs, vecnorm(dIh - dI, Inf))
   @printf(" %d | %.2e \n", p, errs[end])
end
println("---------------")
@test minimum(errs) <= 1e-3 * maximum(errs)
println()


println("[2] Symmetry")
for n = 1:3
   r = 1.0 + SVector(rand(dim)...)
   I = all_invariants(r)
   for rπ in simplex_permutations(r)
      @test I ≈ all_invariants(SVector(rπ...))
   end
   print(".")
end
println()


println("[3] invariants_ed")
for n = 1:3
   r = 1.0 + SVector(rand(dim)...)
   I1, I2 = invariants(r)
   dI1, dI2 = invariants_d(r)
   J1, J2, dJ1, dJ2 = invariants_ed(r)
   @test all(i ≈ j for (i,j) in zip(I1, J1))
   @test all(i ≈ j for (i,j) in zip(I2, J2))
   @test all(di ≈ dj for (di,dj) in zip(dI1, dJ1))
   @test all(di ≈ dj for (di,dj) in zip(dI2, dJ2))
   print(".")
end
println()

println("[4] comparison with non-efficient implementation")
for n = 1:3
   r = 1.0 + SVector(rand(dim)...)
   Iprim, Isec = invariants(r)
   Iprimcheck, Iseccheck, ~ = invariants_check(r)
   @test Iprim ≈ SVector(Iprimcheck...)
   @test Isec ≈ SVector(Iseccheck...)
   print(".")
end
println()
