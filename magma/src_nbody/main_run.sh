#!/bin/bash

source parameters_nbody.conf

. ./script_gen_inv.sh

cat parameters_nbody.conf >> julia_run.jl
cat julia_gen_text_inv_nbody.jl >> julia_run.jl

echo "julia run"

j6 julia_run.jl

# rm julia_run.jl
