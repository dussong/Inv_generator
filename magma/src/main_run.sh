#!/bin/bash

source parameters.conf

. ./script_gen_inv_ba.sh

cat parameters.conf >> julia_run.jl
cat generate_text_invariants_ba.jl >> julia_run.jl

echo "julia run"

julia
# julia_run.jl > logjl

# rm julia_run.jl
