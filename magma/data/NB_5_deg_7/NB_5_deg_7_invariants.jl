module NB5I 

using NBodyIPs.FastPolys 
using StaticArrays 
using BenchmarkTools: @btime 

import NBodyIPs.tdegrees 

tdegrees(::Val{5}) =
          (@SVector [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, ]),
    (@SVector [0, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, ])

# Primary invariants for NBody=5and deg=7 
 # : definitions at the beginning of the file 
const P1_1 = (1,2,3,4,5,6,7,8,9,10,) 

const P2_1 = (1,1,1,2,2,3,1,1,1,2,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,) 
const P2_2 = (2,3,4,3,4,4,5,6,7,5,6,7,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 

const P3_1 = (1,2,3,4,5,6,7,8,9,10,) 

const P4_1 = (1,1,1,2,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,) 
const P4_2 = (2,2,3,3,5,5,6,5,5,6,7,6,7,8,8,9,6,8,8,9,) 
const P4_3 = (3,4,4,4,6,7,7,8,9,8,9,10,10,9,10,10,7,9,10,10,) 

const P5_1 = (1,2,3,4,5,6,7,8,9,10,) 

const P6_1 = (1,1,2,3,4,) 
const P6_2 = (2,5,5,6,7,) 
const P6_3 = (3,6,8,8,9,) 
const P6_4 = (4,7,9,10,10,) 

const P7_1 = (1,2,3,4,5,6,7,8,9,10,) 

const P8_1 = (1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,1,1,1,2,2,3,4,3,4,2,3,4,) 
const P8_2 = (2,2,2,2,2,2,2,3,4,2,3,4,3,4,3,4,4,4,5,5,5,5,5,5,5,6,6,5,6,7,) 
const P8_3 = (3,3,3,3,3,3,5,5,5,5,6,7,5,5,6,7,6,7,6,6,6,6,7,6,7,7,7,8,8,8,) 
const P8_4 = (4,4,4,4,4,4,6,6,6,8,8,9,8,8,8,9,8,9,7,7,7,8,8,8,9,8,9,9,9,9,) 
const P8_5 = (5,6,7,8,9,10,7,7,7,9,10,10,9,9,10,10,10,10,8,9,10,9,9,10,10,10,10,10,10,10,) 

const P9_1 = (1,2,3,4,5,6,7,8,9,10,) 

const P10_1 = (1,2,3,4,5,6,7,8,9,10,) 

# Irreducible secondaries for NBody=5and deg=7 
 # : definitions at the beginning of the file 
const IS1_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS1_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS2_1 = (1,1,1,2,2,3,5,5,6,8,) 
const IS2_2 = (2,3,4,3,4,4,6,7,7,9,) 
const IS2_3 = (5,6,7,8,9,10,8,9,10,10,) 

const IS3_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS3_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS4_1 = (1,1,1,2,2,3,1,1,1,2,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS4_2 = (2,3,4,3,4,4,5,6,7,5,6,7,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 

const IS5_1 = (1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,8,9,10,5,5,6,7,6,7,8,9,10,8,9,10,) 
const IS5_2 = (2,3,4,1,1,1,3,4,2,2,4,3,1,1,1,2,2,3,6,7,5,5,7,6,5,5,6,9,8,8,) 
const IS5_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,3,4,4,8,9,8,9,10,10,6,7,7,10,10,9,) 

const IS6_1 = (1,1,1,2,2,3,4,3,4,2,3,4,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,6,7,5,6,7,8,9,8,9,10,10,) 
const IS6_2 = (2,2,3,1,1,1,1,1,1,3,2,2,5,5,6,5,5,6,7,6,7,8,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,5,5,8,8,9,5,5,6,7,6,7,) 
const IS6_3 = (3,4,4,3,4,2,2,4,3,4,4,3,6,7,7,8,9,8,9,10,10,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,7,7,6,9,10,10,9,8,10,10,8,9,) 

const IS7_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,8,9,8,9,10,10,) 
const IS7_2 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,7,6,7,6,5,5,6,7,5,5,7,6,) 
const IS7_3 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,9,8,10,10,8,9,) 

const IS8_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS8_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS9_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS9_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS10_1 = (1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,8,9,10,5,5,6,7,6,7,8,9,10,8,9,10,) 
const IS10_2 = (2,3,4,1,1,1,3,4,2,2,4,3,1,1,1,2,2,3,6,7,5,5,7,6,5,5,6,9,8,8,) 
const IS10_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,3,4,4,8,9,8,9,10,10,6,7,7,10,10,9,) 

const IS11_1 = (1,1,1,2,2,3,1,1,1,2,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS11_2 = (2,3,4,3,4,4,5,6,7,5,6,7,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 
const IS11_3 = (5,6,7,8,9,10,2,3,4,1,1,1,3,4,2,2,4,3,8,9,10,6,7,5,5,7,6,10,9,8,) 

const IS12_1 = (1,1,1,2,2,3,4,3,4,2,3,4,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,6,7,5,6,7,8,9,8,9,10,10,) 
const IS12_2 = (2,2,3,1,1,1,1,1,1,3,2,2,5,5,6,5,5,6,7,6,7,8,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,5,5,8,8,9,5,5,6,7,6,7,) 
const IS12_3 = (3,4,4,3,4,2,2,4,3,4,4,3,6,7,7,8,9,8,9,10,10,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,7,7,6,9,10,10,9,8,10,10,8,9,) 

const IS13_1 = (1,1,1,1,1,1,2,2,3,2,2,3,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS13_2 = (2,2,3,4,3,4,3,4,4,3,4,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 
const IS13_3 = (3,4,2,2,4,3,1,1,1,4,3,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,2,2,3,4,3,4,2,3,4,7,6,5,9,8,10,10,8,9,5,6,7,) 

const IS14_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,8,9,8,9,10,10,) 
const IS14_2 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,7,6,7,6,5,5,6,7,5,5,7,6,) 
const IS14_3 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,9,8,10,10,8,9,) 

const IS15_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS15_2 = (2,2,2,2,3,3,1,1,1,1,1,1,1,1,1,1,1,1,3,3,2,2,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,2,2,3,3,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,6,5,5,5,5,6,7,5,5,7,6,5,5,5,5,6,6,5,5,6,7,6,7,) 
const IS15_3 = (3,4,3,4,4,4,3,4,2,2,4,3,3,4,2,2,4,3,4,4,4,3,4,3,5,5,5,5,6,6,5,5,6,7,6,7,5,5,6,7,6,7,8,8,8,9,8,9,2,2,3,4,3,4,2,2,3,4,3,4,3,4,3,4,4,4,3,4,3,4,4,4,6,7,5,5,7,6,6,7,5,5,7,6,5,5,5,5,6,6,9,8,9,8,8,8,7,7,7,6,7,6,8,8,8,9,8,9,6,7,6,7,7,7,9,8,9,8,8,8,) 
const IS15_4 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,8,9,8,9,10,10,9,9,10,10,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,10,10,10,10,9,9,8,9,8,9,10,10,9,9,10,10,10,10,9,8,10,10,8,9,10,10,10,10,9,9,) 

const IS16_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS16_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS17_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS17_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS18_1 = (1,1,1,2,2,3,1,1,1,2,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS18_2 = (2,3,4,3,4,4,5,6,7,5,6,7,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 

const IS19_1 = (1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,8,9,10,5,5,6,7,6,7,8,9,10,8,9,10,) 
const IS19_2 = (2,3,4,1,1,1,3,4,2,2,4,3,1,1,1,2,2,3,6,7,5,5,7,6,5,5,6,9,8,8,) 
const IS19_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,3,4,4,8,9,8,9,10,10,6,7,7,10,10,9,) 

const IS20_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS20_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 
const IS20_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,1,1,1,3,4,2,2,4,3,2,3,4,1,1,1,3,4,2,2,4,3,8,9,8,9,10,10,6,7,5,5,7,6,6,7,5,5,7,6,10,10,9,8,9,8,) 

const IS21_1 = (1,1,1,2,2,3,5,5,6,8,) 
const IS21_2 = (2,3,4,3,4,4,6,7,7,9,) 
const IS21_3 = (5,6,7,8,9,10,8,9,10,10,) 

const IS22_1 = (1,1,1,2,2,3,4,3,4,2,3,4,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,6,7,5,6,7,8,9,8,9,10,10,) 
const IS22_2 = (2,2,3,1,1,1,1,1,1,3,2,2,5,5,6,5,5,6,7,6,7,8,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,5,5,8,8,9,5,5,6,7,6,7,) 
const IS22_3 = (3,4,4,3,4,2,2,4,3,4,4,3,6,7,7,8,9,8,9,10,10,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,7,7,6,9,10,10,9,8,10,10,8,9,) 

const IS23_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS23_2 = (2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 
const IS23_3 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,7,6,7,6,5,5,9,8,10,10,8,9,9,8,10,10,8,9,5,5,6,7,6,7,) 

const IS24_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,8,9,8,9,10,10,) 
const IS24_2 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,7,6,7,6,5,5,6,7,5,5,7,6,) 
const IS24_3 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,9,8,10,10,8,9,) 

const IS25_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS25_2 = (2,2,2,2,3,3,1,1,1,1,1,1,1,1,1,1,1,1,3,3,2,2,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,2,2,3,3,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,6,5,5,5,5,6,7,5,5,7,6,5,5,5,5,6,6,5,5,6,7,6,7,) 
const IS25_3 = (3,4,3,4,4,4,3,4,2,2,4,3,3,4,2,2,4,3,4,4,4,3,4,3,5,5,5,5,6,6,5,5,6,7,6,7,5,5,6,7,6,7,8,8,8,9,8,9,2,2,3,4,3,4,2,2,3,4,3,4,3,4,3,4,4,4,3,4,3,4,4,4,6,7,5,5,7,6,6,7,5,5,7,6,5,5,5,5,6,6,9,8,9,8,8,8,7,7,7,6,7,6,8,8,8,9,8,9,6,7,6,7,7,7,9,8,9,8,8,8,) 
const IS25_4 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,8,9,8,9,10,10,9,9,10,10,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,10,10,10,10,9,9,8,9,8,9,10,10,9,9,10,10,10,10,9,8,10,10,8,9,10,10,10,10,9,9,) 

const IS26_1 = (1,1,1,1,1,1,2,2,3,2,2,3,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS26_2 = (2,2,3,4,3,4,3,4,4,3,4,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 
const IS26_3 = (3,4,2,2,4,3,1,1,1,4,3,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,2,2,3,4,3,4,2,3,4,7,6,5,6,7,5,5,7,6,5,6,7,) 
const IS26_4 = (5,5,6,7,6,7,8,9,10,8,9,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,10,6,7,5,5,7,6,10,9,8,8,9,10,9,8,10,10,8,9,10,9,8,) 

const IS27_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS27_2 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,8,9,8,9,10,10,6,7,5,5,7,6,6,7,5,5,7,6,10,10,9,8,9,8,8,9,8,9,10,10,6,7,5,5,7,6,6,7,5,5,7,6,10,10,9,8,9,8,) 
const IS27_3 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,7,6,7,6,5,5,9,8,10,10,8,9,9,8,10,10,8,9,5,5,6,7,6,7,) 

const IS28_1 = (1,2,3,4,5,6,7,8,9,10,) 

const IS29_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS29_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS30_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS30_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS31_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS31_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 

const IS32_1 = (1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,8,9,10,5,5,6,7,6,7,8,9,10,8,9,10,) 
const IS32_2 = (2,3,4,1,1,1,3,4,2,2,4,3,1,1,1,2,2,3,6,7,5,5,7,6,5,5,6,9,8,8,) 
const IS32_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,3,4,4,8,9,8,9,10,10,6,7,7,10,10,9,) 

const IS33_1 = (1,1,1,2,3,4,2,2,3,4,3,4,1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,5,6,7,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS33_2 = (2,3,4,1,1,1,3,4,2,2,4,3,5,6,7,5,6,7,8,9,8,9,10,10,1,1,1,2,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 
const IS33_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,1,1,1,3,4,2,2,4,3,2,3,4,1,1,1,3,4,2,2,4,3,8,9,8,9,10,10,6,7,5,5,7,6,6,7,5,5,7,6,10,10,9,8,9,8,) 

const IS34_1 = (1,1,1,2,2,3,1,1,1,2,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS34_2 = (2,3,4,3,4,4,5,6,7,5,6,7,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 
const IS34_3 = (5,6,7,8,9,10,2,3,4,1,1,1,3,4,2,2,4,3,8,9,10,6,7,5,5,7,6,10,9,8,) 

const IS35_1 = (1,1,1,2,3,4,2,2,3,4,3,4,5,6,7,8,9,10,5,5,6,7,6,7,8,9,10,8,9,10,) 
const IS35_2 = (2,3,4,1,1,1,3,4,2,2,4,3,1,1,1,2,2,3,6,7,5,5,7,6,5,5,6,9,8,8,) 
const IS35_3 = (5,6,7,5,6,7,8,9,8,9,10,10,2,3,4,3,4,4,8,9,8,9,10,10,6,7,7,10,10,9,) 

const IS36_1 = (1,1,1,2,2,3,4,3,4,2,3,4,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,6,7,5,6,7,8,9,8,9,10,10,) 
const IS36_2 = (2,2,3,1,1,1,1,1,1,3,2,2,5,5,6,5,5,6,7,6,7,8,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,5,5,8,8,9,5,5,6,7,6,7,) 
const IS36_3 = (3,4,4,3,4,2,2,4,3,4,4,3,6,7,7,8,9,8,9,10,10,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,7,7,6,9,10,10,9,8,10,10,8,9,) 

const IS37_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS37_2 = (2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,) 
const IS37_3 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,7,6,7,6,5,5,9,8,10,10,8,9,9,8,10,10,8,9,5,5,6,7,6,7,) 

const IS38_1 = (1,1,1,1,1,1,2,2,3,2,2,3,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,5,5,6,5,5,6,7,6,7,8,8,9,) 
const IS38_2 = (2,2,3,4,3,4,3,4,4,3,4,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,) 
const IS38_3 = (3,4,2,2,4,3,1,1,1,4,3,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,2,2,3,4,3,4,2,3,4,7,6,5,9,8,10,10,8,9,5,6,7,) 

const IS39_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,8,9,8,9,10,10,) 
const IS39_2 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,7,6,7,6,5,5,6,7,5,5,7,6,) 
const IS39_3 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,9,8,10,10,8,9,) 

const IS40_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,) 
const IS40_2 = (2,2,2,2,3,3,1,1,1,1,1,1,1,1,1,1,1,1,3,3,2,2,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,2,2,3,3,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,6,6,5,5,5,5,6,7,5,5,7,6,5,5,5,5,6,6,5,5,6,7,6,7,) 
const IS40_3 = (3,4,3,4,4,4,3,4,2,2,4,3,3,4,2,2,4,3,4,4,4,3,4,3,5,5,5,5,6,6,5,5,6,7,6,7,5,5,6,7,6,7,8,8,8,9,8,9,2,2,3,4,3,4,2,2,3,4,3,4,3,4,3,4,4,4,3,4,3,4,4,4,6,7,5,5,7,6,6,7,5,5,7,6,5,5,5,5,6,6,9,8,9,8,8,8,7,7,7,6,7,6,8,8,8,9,8,9,6,7,6,7,7,7,9,8,9,8,8,8,) 
const IS40_4 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,8,9,8,9,10,10,9,9,10,10,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,10,10,10,10,9,9,8,9,8,9,10,10,9,9,10,10,10,10,9,8,10,10,8,9,10,10,10,10,9,9,) 


# Primary invariants for NBody=5and deg=7 
 # : definitions of the types at the beginning of the file 
const P1 = Val((P1_1,)) 
const P2 = Val((P2_1,P2_2,)) 
const P3 = Val((P3_1,)) 
const P4 = Val((P4_1,P4_2,P4_3,)) 
const P5 = Val((P5_1,)) 
const P6 = Val((P6_1,P6_2,P6_3,P6_4,)) 
const P7 = Val((P7_1,)) 
const P8 = Val((P8_1,P8_2,P8_3,P8_4,P8_5,)) 
const P9 = Val((P9_1,)) 
const P10 = Val((P10_1,)) 
# Irreducible secondaries for NBody=5and deg=7 
 # : definitions of the types at the beginning of the file 
const IS1 = Val((IS1_1,IS1_2,)) 
const IS2 = Val((IS2_1,IS2_2,IS2_3,)) 
const IS3 = Val((IS3_1,IS3_2,)) 
const IS4 = Val((IS4_1,IS4_2,)) 
const IS5 = Val((IS5_1,IS5_2,IS5_3,)) 
const IS6 = Val((IS6_1,IS6_2,IS6_3,)) 
const IS7 = Val((IS7_1,IS7_2,IS7_3,)) 
const IS8 = Val((IS8_1,IS8_2,)) 
const IS9 = Val((IS9_1,IS9_2,)) 
const IS10 = Val((IS10_1,IS10_2,IS10_3,)) 
const IS11 = Val((IS11_1,IS11_2,IS11_3,)) 
const IS12 = Val((IS12_1,IS12_2,IS12_3,)) 
const IS13 = Val((IS13_1,IS13_2,IS13_3,)) 
const IS14 = Val((IS14_1,IS14_2,IS14_3,)) 
const IS15 = Val((IS15_1,IS15_2,IS15_3,IS15_4,)) 
const IS16 = Val((IS16_1,IS16_2,)) 
const IS17 = Val((IS17_1,IS17_2,)) 
const IS18 = Val((IS18_1,IS18_2,)) 
const IS19 = Val((IS19_1,IS19_2,IS19_3,)) 
const IS20 = Val((IS20_1,IS20_2,IS20_3,)) 
const IS21 = Val((IS21_1,IS21_2,IS21_3,)) 
const IS22 = Val((IS22_1,IS22_2,IS22_3,)) 
const IS23 = Val((IS23_1,IS23_2,IS23_3,)) 
const IS24 = Val((IS24_1,IS24_2,IS24_3,)) 
const IS25 = Val((IS25_1,IS25_2,IS25_3,IS25_4,)) 
const IS26 = Val((IS26_1,IS26_2,IS26_3,IS26_4,)) 
const IS27 = Val((IS27_1,IS27_2,IS27_3,)) 
const IS28 = Val((IS28_1,)) 
const IS29 = Val((IS29_1,IS29_2,)) 
const IS30 = Val((IS30_1,IS30_2,)) 
const IS31 = Val((IS31_1,IS31_2,)) 
const IS32 = Val((IS32_1,IS32_2,IS32_3,)) 
const IS33 = Val((IS33_1,IS33_2,IS33_3,)) 
const IS34 = Val((IS34_1,IS34_2,IS34_3,)) 
const IS35 = Val((IS35_1,IS35_2,IS35_3,)) 
const IS36 = Val((IS36_1,IS36_2,IS36_3,)) 
const IS37 = Val((IS37_1,IS37_2,IS37_3,)) 
const IS38 = Val((IS38_1,IS38_2,IS38_3,)) 
const IS39 = Val((IS39_1,IS39_2,IS39_3,)) 
const IS40 = Val((IS40_1,IS40_2,IS40_3,IS40_4,)) 


function invariants_gen(x1::SVector{10, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 
   x5 = x4.*x1 
   x6 = x5.*x1 
   x7 = x6.*x1 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NBody=5and deg=7 
 # : what goes in the function for the evaluation 
P1 = fpoly((x1,) , NB5I.P1) 
P2 = fpoly((x1,x1,) , NB5I.P2) 
P3 = fpoly((x2,) , NB5I.P3) 
P4 = fpoly((x1,x1,x1,) , NB5I.P4) 
P5 = fpoly((x3,) , NB5I.P5) 
P6 = fpoly((x1,x1,x1,x1,) , NB5I.P6) 
P7 = fpoly((x4,) , NB5I.P7) 
P8 = fpoly((x1,x1,x1,x1,x1,) , NB5I.P8) 
P9 = fpoly((x5,) , NB5I.P9) 
P10 = fpoly((x6,) , NB5I.P10) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for NBody=5and deg=7 
 # : what goes in the function for the evaluation 
IS1 = fpoly((x2,x1,) , NB5I.IS1) 
IS2 = fpoly((x1,x1,x1,) , NB5I.IS2) 
IS3 = fpoly((x3,x1,) , NB5I.IS3) 
IS4 = fpoly((x2,x2,) , NB5I.IS4) 
IS5 = fpoly((x2,x1,x1,) , NB5I.IS5) 
IS6 = fpoly((x2,x1,x1,) , NB5I.IS6) 
IS7 = fpoly((x2,x1,x1,) , NB5I.IS7) 
IS8 = fpoly((x4,x1,) , NB5I.IS8) 
IS9 = fpoly((x3,x2,) , NB5I.IS9) 
IS10 = fpoly((x3,x1,x1,) , NB5I.IS10) 
IS11 = fpoly((x2,x2,x1,) , NB5I.IS11) 
IS12 = fpoly((x3,x1,x1,) , NB5I.IS12) 
IS13 = fpoly((x2,x2,x1,) , NB5I.IS13) 
IS14 = fpoly((x3,x1,x1,) , NB5I.IS14) 
IS15 = fpoly((x2,x1,x1,x1,) , NB5I.IS15) 
IS16 = fpoly((x5,x1,) , NB5I.IS16) 
IS17 = fpoly((x4,x2,) , NB5I.IS17) 
IS18 = fpoly((x3,x3,) , NB5I.IS18) 
IS19 = fpoly((x4,x1,x1,) , NB5I.IS19) 
IS20 = fpoly((x3,x2,x1,) , NB5I.IS20) 
IS21 = fpoly((x2,x2,x2,) , NB5I.IS21) 
IS22 = fpoly((x4,x1,x1,) , NB5I.IS22) 
IS23 = fpoly((x3,x2,x1,) , NB5I.IS23) 
IS24 = fpoly((x4,x1,x1,) , NB5I.IS24) 
IS25 = fpoly((x3,x1,x1,x1,) , NB5I.IS25) 
IS26 = fpoly((x2,x2,x1,x1,) , NB5I.IS26) 
IS27 = fpoly((x3,x2,x1,) , NB5I.IS27) 
IS28 = fpoly((x7,) , NB5I.IS28) 
IS29 = fpoly((x6,x1,) , NB5I.IS29) 
IS30 = fpoly((x5,x2,) , NB5I.IS30) 
IS31 = fpoly((x4,x3,) , NB5I.IS31) 
IS32 = fpoly((x5,x1,x1,) , NB5I.IS32) 
IS33 = fpoly((x4,x2,x1,) , NB5I.IS33) 
IS34 = fpoly((x3,x3,x1,) , NB5I.IS34) 
IS35 = fpoly((x3,x2,x2,) , NB5I.IS35) 
IS36 = fpoly((x5,x1,x1,) , NB5I.IS36) 
IS37 = fpoly((x4,x2,x1,) , NB5I.IS37) 
IS38 = fpoly((x3,x3,x1,) , NB5I.IS38) 
IS39 = fpoly((x5,x1,x1,) , NB5I.IS39) 
IS40 = fpoly((x4,x1,x1,x1,) , NB5I.IS40) 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


 SEC1  = 1
 SEC2  = IS1 
 SEC3  = IS2 
 SEC4  = IS3 
 SEC5  = IS4 
 SEC6  = IS5 
 SEC7  = IS6 
 SEC8  = IS7 
 SEC9  = IS8 
 SEC10  = IS9 
 SEC11  = IS10 
 SEC12  = IS11 
 SEC13  = IS12 
 SEC14  = IS13 
 SEC15  = IS14 
 SEC16  = IS15 
 SEC17  = IS1 *IS1 
 SEC18  = IS1 *IS2 
 SEC19  = IS2 *IS2 
 SEC20  = IS16 
 SEC21  = IS17 
 SEC22  = IS18 
 SEC23  = IS19 
 SEC24  = IS20 
 SEC25  = IS21 
 SEC26  = IS22 
 SEC27  = IS23 
 SEC28  = IS24 
 SEC29  = IS25 
 SEC30  = IS26 
 SEC31  = IS27 
 SEC32  = IS1 *IS3 
 SEC33  = IS2 *IS3 
 SEC34  = IS1 *IS4 
 SEC35  = IS2 *IS4 
 SEC36  = IS1 *IS5 
 SEC37  = IS2 *IS5 
 SEC38  = IS1 *IS6 
 SEC39  = IS2 *IS6 
 SEC40  = IS1 *IS7 
 SEC41  = IS2 *IS7 
 SEC42  = IS28 
 SEC43  = IS29 
 SEC44  = IS30 
 SEC45  = IS31 
 SEC46  = IS32 
 SEC47  = IS33 
 SEC48  = IS34 
 SEC49  = IS35 
 SEC50  = IS36 
 SEC51  = IS37 
 SEC52  = IS38 
 SEC53  = IS39 
 SEC54  = IS40 


return (@SVector [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,SEC13,SEC14,SEC15,SEC16,SEC17,SEC18,SEC19,SEC20,SEC21,SEC22,SEC23,SEC24,SEC25,SEC26,SEC27,SEC28,SEC29,SEC30,SEC31,SEC32,SEC33,SEC34,SEC35,SEC36,SEC37,SEC38,SEC39,SEC40,SEC41,SEC42,SEC43,SEC44,SEC45,SEC46,SEC47,SEC48,SEC49,SEC50,SEC51,SEC52,SEC53,SEC54,])
 end



function invariants_d_gen(x1::SVector{10, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 
   x5 = x4.*x1 
   x6 = x5.*x1 
   x7 = x6.*x1 

   dx1 = @SVector ones(10)
   dx2 = 2 * x1 
   dx3 = 3 * x2 
   dx4 = 4 * x3 
   dx5 = 5 * x4 
   dx6 = 6 * x5 
   dx7 = 7 * x6 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NBody=5and deg=7 
 # : what goes in the function for the derivatives 
dP1 = fpoly_d((x1,),(dx1,) , NB5I.P1) 
dP2 = fpoly_d((x1,x1,),(dx1,dx1,) , NB5I.P2) 
dP3 = fpoly_d((x2,),(dx2,) , NB5I.P3) 
dP4 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , NB5I.P4) 
dP5 = fpoly_d((x3,),(dx3,) , NB5I.P5) 
dP6 = fpoly_d((x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,) , NB5I.P6) 
dP7 = fpoly_d((x4,),(dx4,) , NB5I.P7) 
dP8 = fpoly_d((x1,x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,dx1,) , NB5I.P8) 
dP9 = fpoly_d((x5,),(dx5,) , NB5I.P9) 
dP10 = fpoly_d((x6,),(dx6,) , NB5I.P10) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for NBody=5and deg=7 
 # : what goes in the function for the evaluation 
IS1 = fpoly((x2,x1,) , NB5I.IS1) 
IS2 = fpoly((x1,x1,x1,) , NB5I.IS2) 
IS3 = fpoly((x3,x1,) , NB5I.IS3) 
IS4 = fpoly((x2,x2,) , NB5I.IS4) 
IS5 = fpoly((x2,x1,x1,) , NB5I.IS5) 
IS6 = fpoly((x2,x1,x1,) , NB5I.IS6) 
IS7 = fpoly((x2,x1,x1,) , NB5I.IS7) 
IS8 = fpoly((x4,x1,) , NB5I.IS8) 
IS9 = fpoly((x3,x2,) , NB5I.IS9) 
IS10 = fpoly((x3,x1,x1,) , NB5I.IS10) 
IS11 = fpoly((x2,x2,x1,) , NB5I.IS11) 
IS12 = fpoly((x3,x1,x1,) , NB5I.IS12) 
IS13 = fpoly((x2,x2,x1,) , NB5I.IS13) 
IS14 = fpoly((x3,x1,x1,) , NB5I.IS14) 
IS15 = fpoly((x2,x1,x1,x1,) , NB5I.IS15) 
IS16 = fpoly((x5,x1,) , NB5I.IS16) 
IS17 = fpoly((x4,x2,) , NB5I.IS17) 
IS18 = fpoly((x3,x3,) , NB5I.IS18) 
IS19 = fpoly((x4,x1,x1,) , NB5I.IS19) 
IS20 = fpoly((x3,x2,x1,) , NB5I.IS20) 
IS21 = fpoly((x2,x2,x2,) , NB5I.IS21) 
IS22 = fpoly((x4,x1,x1,) , NB5I.IS22) 
IS23 = fpoly((x3,x2,x1,) , NB5I.IS23) 
IS24 = fpoly((x4,x1,x1,) , NB5I.IS24) 
IS25 = fpoly((x3,x1,x1,x1,) , NB5I.IS25) 
IS26 = fpoly((x2,x2,x1,x1,) , NB5I.IS26) 
IS27 = fpoly((x3,x2,x1,) , NB5I.IS27) 
IS28 = fpoly((x7,) , NB5I.IS28) 
IS29 = fpoly((x6,x1,) , NB5I.IS29) 
IS30 = fpoly((x5,x2,) , NB5I.IS30) 
IS31 = fpoly((x4,x3,) , NB5I.IS31) 
IS32 = fpoly((x5,x1,x1,) , NB5I.IS32) 
IS33 = fpoly((x4,x2,x1,) , NB5I.IS33) 
IS34 = fpoly((x3,x3,x1,) , NB5I.IS34) 
IS35 = fpoly((x3,x2,x2,) , NB5I.IS35) 
IS36 = fpoly((x5,x1,x1,) , NB5I.IS36) 
IS37 = fpoly((x4,x2,x1,) , NB5I.IS37) 
IS38 = fpoly((x3,x3,x1,) , NB5I.IS38) 
IS39 = fpoly((x5,x1,x1,) , NB5I.IS39) 
IS40 = fpoly((x4,x1,x1,x1,) , NB5I.IS40) 


# Irreducible secondaries for NBody=5and deg=7 
 # : what goes in the function for the derivatives 
dIS1 = fpoly_d((x2,x1,),(dx2,dx1,) , NB5I.IS1) 
dIS2 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , NB5I.IS2) 
dIS3 = fpoly_d((x3,x1,),(dx3,dx1,) , NB5I.IS3) 
dIS4 = fpoly_d((x2,x2,),(dx2,dx2,) , NB5I.IS4) 
dIS5 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , NB5I.IS5) 
dIS6 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , NB5I.IS6) 
dIS7 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , NB5I.IS7) 
dIS8 = fpoly_d((x4,x1,),(dx4,dx1,) , NB5I.IS8) 
dIS9 = fpoly_d((x3,x2,),(dx3,dx2,) , NB5I.IS9) 
dIS10 = fpoly_d((x3,x1,x1,),(dx3,dx1,dx1,) , NB5I.IS10) 
dIS11 = fpoly_d((x2,x2,x1,),(dx2,dx2,dx1,) , NB5I.IS11) 
dIS12 = fpoly_d((x3,x1,x1,),(dx3,dx1,dx1,) , NB5I.IS12) 
dIS13 = fpoly_d((x2,x2,x1,),(dx2,dx2,dx1,) , NB5I.IS13) 
dIS14 = fpoly_d((x3,x1,x1,),(dx3,dx1,dx1,) , NB5I.IS14) 
dIS15 = fpoly_d((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , NB5I.IS15) 
dIS16 = fpoly_d((x5,x1,),(dx5,dx1,) , NB5I.IS16) 
dIS17 = fpoly_d((x4,x2,),(dx4,dx2,) , NB5I.IS17) 
dIS18 = fpoly_d((x3,x3,),(dx3,dx3,) , NB5I.IS18) 
dIS19 = fpoly_d((x4,x1,x1,),(dx4,dx1,dx1,) , NB5I.IS19) 
dIS20 = fpoly_d((x3,x2,x1,),(dx3,dx2,dx1,) , NB5I.IS20) 
dIS21 = fpoly_d((x2,x2,x2,),(dx2,dx2,dx2,) , NB5I.IS21) 
dIS22 = fpoly_d((x4,x1,x1,),(dx4,dx1,dx1,) , NB5I.IS22) 
dIS23 = fpoly_d((x3,x2,x1,),(dx3,dx2,dx1,) , NB5I.IS23) 
dIS24 = fpoly_d((x4,x1,x1,),(dx4,dx1,dx1,) , NB5I.IS24) 
dIS25 = fpoly_d((x3,x1,x1,x1,),(dx3,dx1,dx1,dx1,) , NB5I.IS25) 
dIS26 = fpoly_d((x2,x2,x1,x1,),(dx2,dx2,dx1,dx1,) , NB5I.IS26) 
dIS27 = fpoly_d((x3,x2,x1,),(dx3,dx2,dx1,) , NB5I.IS27) 
dIS28 = fpoly_d((x7,),(dx7,) , NB5I.IS28) 
dIS29 = fpoly_d((x6,x1,),(dx6,dx1,) , NB5I.IS29) 
dIS30 = fpoly_d((x5,x2,),(dx5,dx2,) , NB5I.IS30) 
dIS31 = fpoly_d((x4,x3,),(dx4,dx3,) , NB5I.IS31) 
dIS32 = fpoly_d((x5,x1,x1,),(dx5,dx1,dx1,) , NB5I.IS32) 
dIS33 = fpoly_d((x4,x2,x1,),(dx4,dx2,dx1,) , NB5I.IS33) 
dIS34 = fpoly_d((x3,x3,x1,),(dx3,dx3,dx1,) , NB5I.IS34) 
dIS35 = fpoly_d((x3,x2,x2,),(dx3,dx2,dx2,) , NB5I.IS35) 
dIS36 = fpoly_d((x5,x1,x1,),(dx5,dx1,dx1,) , NB5I.IS36) 
dIS37 = fpoly_d((x4,x2,x1,),(dx4,dx2,dx1,) , NB5I.IS37) 
dIS38 = fpoly_d((x3,x3,x1,),(dx3,dx3,dx1,) , NB5I.IS38) 
dIS39 = fpoly_d((x5,x1,x1,),(dx5,dx1,dx1,) , NB5I.IS39) 
dIS40 = fpoly_d((x4,x1,x1,x1,),(dx4,dx1,dx1,dx1,) , NB5I.IS40) 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


 dSEC1   = @SVector zeros(10) 
 dSEC2  = dIS1 
 dSEC3  = dIS2 
 dSEC4  = dIS3 
 dSEC5  = dIS4 
 dSEC6  = dIS5 
 dSEC7  = dIS6 
 dSEC8  = dIS7 
 dSEC9  = dIS8 
 dSEC10  = dIS9 
 dSEC11  = dIS10 
 dSEC12  = dIS11 
 dSEC13  = dIS12 
 dSEC14  = dIS13 
 dSEC15  = dIS14 
 dSEC16  = dIS15 
 dSEC17   =  dIS1 *IS1 + IS1 *dIS1 
 dSEC18   =  dIS1 *IS2 + IS1 *dIS2 
 dSEC19   =  dIS2 *IS2 + IS2 *dIS2 
 dSEC20  = dIS16 
 dSEC21  = dIS17 
 dSEC22  = dIS18 
 dSEC23  = dIS19 
 dSEC24  = dIS20 
 dSEC25  = dIS21 
 dSEC26  = dIS22 
 dSEC27  = dIS23 
 dSEC28  = dIS24 
 dSEC29  = dIS25 
 dSEC30  = dIS26 
 dSEC31  = dIS27 
 dSEC32   =  dIS1 *IS3 + IS1 *dIS3 
 dSEC33   =  dIS2 *IS3 + IS2 *dIS3 
 dSEC34   =  dIS1 *IS4 + IS1 *dIS4 
 dSEC35   =  dIS2 *IS4 + IS2 *dIS4 
 dSEC36   =  dIS1 *IS5 + IS1 *dIS5 
 dSEC37   =  dIS2 *IS5 + IS2 *dIS5 
 dSEC38   =  dIS1 *IS6 + IS1 *dIS6 
 dSEC39   =  dIS2 *IS6 + IS2 *dIS6 
 dSEC40   =  dIS1 *IS7 + IS1 *dIS7 
 dSEC41   =  dIS2 *IS7 + IS2 *dIS7 
 dSEC42  = dIS28 
 dSEC43  = dIS29 
 dSEC44  = dIS30 
 dSEC45  = dIS31 
 dSEC46  = dIS32 
 dSEC47  = dIS33 
 dSEC48  = dIS34 
 dSEC49  = dIS35 
 dSEC50  = dIS36 
 dSEC51  = dIS37 
 dSEC52  = dIS38 
 dSEC53  = dIS39 
 dSEC54  = dIS40 


return (dP1,dP2,dP3,dP4,dP5,dP6,dP7,dP8,dP9,dP10,), (dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,dSEC7,dSEC8,dSEC9,dSEC10,dSEC11,dSEC12,dSEC13,dSEC14,dSEC15,dSEC16,dSEC17,dSEC18,dSEC19,dSEC20,dSEC21,dSEC22,dSEC23,dSEC24,dSEC25,dSEC26,dSEC27,dSEC28,dSEC29,dSEC30,dSEC31,dSEC32,dSEC33,dSEC34,dSEC35,dSEC36,dSEC37,dSEC38,dSEC39,dSEC40,dSEC41,dSEC42,dSEC43,dSEC44,dSEC45,dSEC46,dSEC47,dSEC48,dSEC49,dSEC50,dSEC51,dSEC52,dSEC53,dSEC54,)
 end



function invariants_ed_gen(x1::SVector{10, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 
   x5 = x4.*x1 
   x6 = x5.*x1 
   x7 = x6.*x1 

   dx1 = @SVector ones(10)
   dx2 = 2 * x1 
   dx3 = 3 * x2 
   dx4 = 4 * x3 
   dx5 = 5 * x4 
   dx6 = 6 * x5 
   dx7 = 7 * x6 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NBody=5and deg=7 
 # : what goes in the function for the evaluation and derivatives 
P1, dP1 = fpoly_ed((x1,),(dx1,) , NB5I.P1) 
P2, dP2 = fpoly_ed((x1,x1,),(dx1,dx1,) , NB5I.P2) 
P3, dP3 = fpoly_ed((x2,),(dx2,) , NB5I.P3) 
P4, dP4 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , NB5I.P4) 
P5, dP5 = fpoly_ed((x3,),(dx3,) , NB5I.P5) 
P6, dP6 = fpoly_ed((x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,) , NB5I.P6) 
P7, dP7 = fpoly_ed((x4,),(dx4,) , NB5I.P7) 
P8, dP8 = fpoly_ed((x1,x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,dx1,) , NB5I.P8) 
P9, dP9 = fpoly_ed((x5,),(dx5,) , NB5I.P9) 
P10, dP10 = fpoly_ed((x6,),(dx6,) , NB5I.P10) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for NBody=5and deg=7 
 # : what goes in the function for the evaluation and derivatives 
IS1, dIS1 = fpoly_ed((x2,x1,),(dx2,dx1,) , NB5I.IS1) 
IS2, dIS2 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , NB5I.IS2) 
IS3, dIS3 = fpoly_ed((x3,x1,),(dx3,dx1,) , NB5I.IS3) 
IS4, dIS4 = fpoly_ed((x2,x2,),(dx2,dx2,) , NB5I.IS4) 
IS5, dIS5 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , NB5I.IS5) 
IS6, dIS6 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , NB5I.IS6) 
IS7, dIS7 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , NB5I.IS7) 
IS8, dIS8 = fpoly_ed((x4,x1,),(dx4,dx1,) , NB5I.IS8) 
IS9, dIS9 = fpoly_ed((x3,x2,),(dx3,dx2,) , NB5I.IS9) 
IS10, dIS10 = fpoly_ed((x3,x1,x1,),(dx3,dx1,dx1,) , NB5I.IS10) 
IS11, dIS11 = fpoly_ed((x2,x2,x1,),(dx2,dx2,dx1,) , NB5I.IS11) 
IS12, dIS12 = fpoly_ed((x3,x1,x1,),(dx3,dx1,dx1,) , NB5I.IS12) 
IS13, dIS13 = fpoly_ed((x2,x2,x1,),(dx2,dx2,dx1,) , NB5I.IS13) 
IS14, dIS14 = fpoly_ed((x3,x1,x1,),(dx3,dx1,dx1,) , NB5I.IS14) 
IS15, dIS15 = fpoly_ed((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , NB5I.IS15) 
IS16, dIS16 = fpoly_ed((x5,x1,),(dx5,dx1,) , NB5I.IS16) 
IS17, dIS17 = fpoly_ed((x4,x2,),(dx4,dx2,) , NB5I.IS17) 
IS18, dIS18 = fpoly_ed((x3,x3,),(dx3,dx3,) , NB5I.IS18) 
IS19, dIS19 = fpoly_ed((x4,x1,x1,),(dx4,dx1,dx1,) , NB5I.IS19) 
IS20, dIS20 = fpoly_ed((x3,x2,x1,),(dx3,dx2,dx1,) , NB5I.IS20) 
IS21, dIS21 = fpoly_ed((x2,x2,x2,),(dx2,dx2,dx2,) , NB5I.IS21) 
IS22, dIS22 = fpoly_ed((x4,x1,x1,),(dx4,dx1,dx1,) , NB5I.IS22) 
IS23, dIS23 = fpoly_ed((x3,x2,x1,),(dx3,dx2,dx1,) , NB5I.IS23) 
IS24, dIS24 = fpoly_ed((x4,x1,x1,),(dx4,dx1,dx1,) , NB5I.IS24) 
IS25, dIS25 = fpoly_ed((x3,x1,x1,x1,),(dx3,dx1,dx1,dx1,) , NB5I.IS25) 
IS26, dIS26 = fpoly_ed((x2,x2,x1,x1,),(dx2,dx2,dx1,dx1,) , NB5I.IS26) 
IS27, dIS27 = fpoly_ed((x3,x2,x1,),(dx3,dx2,dx1,) , NB5I.IS27) 
IS28, dIS28 = fpoly_ed((x7,),(dx7,) , NB5I.IS28) 
IS29, dIS29 = fpoly_ed((x6,x1,),(dx6,dx1,) , NB5I.IS29) 
IS30, dIS30 = fpoly_ed((x5,x2,),(dx5,dx2,) , NB5I.IS30) 
IS31, dIS31 = fpoly_ed((x4,x3,),(dx4,dx3,) , NB5I.IS31) 
IS32, dIS32 = fpoly_ed((x5,x1,x1,),(dx5,dx1,dx1,) , NB5I.IS32) 
IS33, dIS33 = fpoly_ed((x4,x2,x1,),(dx4,dx2,dx1,) , NB5I.IS33) 
IS34, dIS34 = fpoly_ed((x3,x3,x1,),(dx3,dx3,dx1,) , NB5I.IS34) 
IS35, dIS35 = fpoly_ed((x3,x2,x2,),(dx3,dx2,dx2,) , NB5I.IS35) 
IS36, dIS36 = fpoly_ed((x5,x1,x1,),(dx5,dx1,dx1,) , NB5I.IS36) 
IS37, dIS37 = fpoly_ed((x4,x2,x1,),(dx4,dx2,dx1,) , NB5I.IS37) 
IS38, dIS38 = fpoly_ed((x3,x3,x1,),(dx3,dx3,dx1,) , NB5I.IS38) 
IS39, dIS39 = fpoly_ed((x5,x1,x1,),(dx5,dx1,dx1,) , NB5I.IS39) 
IS40, dIS40 = fpoly_ed((x4,x1,x1,x1,),(dx4,dx1,dx1,dx1,) , NB5I.IS40) 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


 SEC1  = 1
 SEC2  = IS1 
 SEC3  = IS2 
 SEC4  = IS3 
 SEC5  = IS4 
 SEC6  = IS5 
 SEC7  = IS6 
 SEC8  = IS7 
 SEC9  = IS8 
 SEC10  = IS9 
 SEC11  = IS10 
 SEC12  = IS11 
 SEC13  = IS12 
 SEC14  = IS13 
 SEC15  = IS14 
 SEC16  = IS15 
 SEC17  = IS1 *IS1 
 SEC18  = IS1 *IS2 
 SEC19  = IS2 *IS2 
 SEC20  = IS16 
 SEC21  = IS17 
 SEC22  = IS18 
 SEC23  = IS19 
 SEC24  = IS20 
 SEC25  = IS21 
 SEC26  = IS22 
 SEC27  = IS23 
 SEC28  = IS24 
 SEC29  = IS25 
 SEC30  = IS26 
 SEC31  = IS27 
 SEC32  = IS1 *IS3 
 SEC33  = IS2 *IS3 
 SEC34  = IS1 *IS4 
 SEC35  = IS2 *IS4 
 SEC36  = IS1 *IS5 
 SEC37  = IS2 *IS5 
 SEC38  = IS1 *IS6 
 SEC39  = IS2 *IS6 
 SEC40  = IS1 *IS7 
 SEC41  = IS2 *IS7 
 SEC42  = IS28 
 SEC43  = IS29 
 SEC44  = IS30 
 SEC45  = IS31 
 SEC46  = IS32 
 SEC47  = IS33 
 SEC48  = IS34 
 SEC49  = IS35 
 SEC50  = IS36 
 SEC51  = IS37 
 SEC52  = IS38 
 SEC53  = IS39 
 SEC54  = IS40 


 dSEC1   = @SVector zeros(10) 
 dSEC2  = dIS1 
 dSEC3  = dIS2 
 dSEC4  = dIS3 
 dSEC5  = dIS4 
 dSEC6  = dIS5 
 dSEC7  = dIS6 
 dSEC8  = dIS7 
 dSEC9  = dIS8 
 dSEC10  = dIS9 
 dSEC11  = dIS10 
 dSEC12  = dIS11 
 dSEC13  = dIS12 
 dSEC14  = dIS13 
 dSEC15  = dIS14 
 dSEC16  = dIS15 
 dSEC17   =  dIS1 *IS1 + IS1 *dIS1 
 dSEC18   =  dIS1 *IS2 + IS1 *dIS2 
 dSEC19   =  dIS2 *IS2 + IS2 *dIS2 
 dSEC20  = dIS16 
 dSEC21  = dIS17 
 dSEC22  = dIS18 
 dSEC23  = dIS19 
 dSEC24  = dIS20 
 dSEC25  = dIS21 
 dSEC26  = dIS22 
 dSEC27  = dIS23 
 dSEC28  = dIS24 
 dSEC29  = dIS25 
 dSEC30  = dIS26 
 dSEC31  = dIS27 
 dSEC32   =  dIS1 *IS3 + IS1 *dIS3 
 dSEC33   =  dIS2 *IS3 + IS2 *dIS3 
 dSEC34   =  dIS1 *IS4 + IS1 *dIS4 
 dSEC35   =  dIS2 *IS4 + IS2 *dIS4 
 dSEC36   =  dIS1 *IS5 + IS1 *dIS5 
 dSEC37   =  dIS2 *IS5 + IS2 *dIS5 
 dSEC38   =  dIS1 *IS6 + IS1 *dIS6 
 dSEC39   =  dIS2 *IS6 + IS2 *dIS6 
 dSEC40   =  dIS1 *IS7 + IS1 *dIS7 
 dSEC41   =  dIS2 *IS7 + IS2 *dIS7 
 dSEC42  = dIS28 
 dSEC43  = dIS29 
 dSEC44  = dIS30 
 dSEC45  = dIS31 
 dSEC46  = dIS32 
 dSEC47  = dIS33 
 dSEC48  = dIS34 
 dSEC49  = dIS35 
 dSEC50  = dIS36 
 dSEC51  = dIS37 
 dSEC52  = dIS38 
 dSEC53  = dIS39 
 dSEC54  = dIS40 


return (@SVector [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,SEC13,SEC14,SEC15,SEC16,SEC17,SEC18,SEC19,SEC20,SEC21,SEC22,SEC23,SEC24,SEC25,SEC26,SEC27,SEC28,SEC29,SEC30,SEC31,SEC32,SEC33,SEC34,SEC35,SEC36,SEC37,SEC38,SEC39,SEC40,SEC41,SEC42,SEC43,SEC44,SEC45,SEC46,SEC47,SEC48,SEC49,SEC50,SEC51,SEC52,SEC53,SEC54,]), (@SVector [dP1,dP2,dP3,dP4,dP5,dP6,dP7,dP8,dP9,dP10,]), (@SVector [dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,dSEC7,dSEC8,dSEC9,dSEC10,dSEC11,dSEC12,dSEC13,dSEC14,dSEC15,dSEC16,dSEC17,dSEC18,dSEC19,dSEC20,dSEC21,dSEC22,dSEC23,dSEC24,dSEC25,dSEC26,dSEC27,dSEC28,dSEC29,dSEC30,dSEC31,dSEC32,dSEC33,dSEC34,dSEC35,dSEC36,dSEC37,dSEC38,dSEC39,dSEC40,dSEC41,dSEC42,dSEC43,dSEC44,dSEC45,dSEC46,dSEC47,dSEC48,dSEC49,dSEC50,dSEC51,dSEC52,dSEC53,dSEC54,])
 end 

end