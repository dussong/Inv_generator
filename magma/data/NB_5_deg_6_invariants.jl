using StaticArrays
using BenchmarkTools: @btime

# Primary invariants for NBody=5and deg=6
 # : definitions at the beginning of the file
const P1_1 = (1,2,3,4,5,6,7,8,9,10,)

const P2_1 = (1,1,1,2,2,3,1,1,1,2,3,4,2,2,3,4,3,4,5,5,6,5,5,6,7,6,7,8,8,9,)
const P2_2 = (2,3,4,3,4,4,5,6,7,5,6,7,8,9,8,9,10,10,6,7,7,8,9,8,9,10,10,9,10,10,)

const P3_1 = (1,2,3,4,5,6,7,8,9,10,)

const P4_1 = (1,1,1,2,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,)
const P4_2 = (2,2,3,3,5,5,6,5,5,6,7,6,7,8,8,9,6,8,8,9,)
const P4_3 = (3,4,4,4,6,7,7,8,9,8,9,10,10,9,10,10,7,9,10,10,)

const P5_1 = (1,2,3,4,5,6,7,8,9,10,)

const P6_1 = (4,3,2,1,2,2,3,1,1,1,)
const P6_2 = (5,5,6,8,3,4,4,3,4,2,)
const P6_3 = (6,7,7,9,7,6,5,6,7,5,)
const P6_4 = (8,9,10,10,8,9,10,9,8,10,)

const P7_1 = (1,2,3,4,5,6,7,8,9,10,)

const P8_1 = (1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,)
const P8_2 = (2,2,3,3,2,2,2,2,3,3,2,2,3,3,3,3,3,3,2,2,2,2,3,3,2,2,3,4,3,4,3,4,3,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,5,5,5,5,6,6,5,5,6,7,6,7,)
const P8_3 = (4,3,4,4,3,4,4,3,4,4,3,4,4,4,4,4,4,4,3,4,3,4,4,4,5,5,6,6,5,5,6,7,5,5,7,6,6,6,6,6,6,6,6,7,6,7,7,7,6,7,6,7,7,7,8,8,8,8,8,8,)
const P8_4 = (5,5,5,5,6,6,5,5,6,7,6,7,5,5,6,7,6,7,8,8,8,9,8,9,7,6,7,7,6,7,8,8,8,9,8,9,7,7,7,7,7,7,8,8,8,9,8,9,8,8,8,9,8,9,9,9,9,9,9,9,)
const P8_5 = (6,7,6,7,7,7,8,9,8,9,10,10,8,9,8,9,10,10,9,9,10,10,10,10,8,9,8,9,10,10,9,9,10,10,10,10,8,9,8,9,10,10,9,9,10,10,10,10,9,9,10,10,10,10,10,10,10,10,10,10,)

const P9_1 = (1,2,3,4,5,6,7,8,9,10,)

const P10_1 = (1,2,3,4,5,6,7,8,9,10,)

# Irreducible secondaries for NBody=5and deg=6
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

const IS15_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,2,2,3,3,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,5,5,6,7,6,7,5,5,5,5,6,6,6,7,5,5,7,6,8,9,8,8,9,8,5,5,6,7,6,7,5,5,5,5,6,6,6,7,5,5,7,6,8,9,8,8,9,8,)
const IS15_2 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,2,2,3,4,3,4,2,2,3,4,3,4,3,4,3,4,4,4,3,4,3,4,4,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,10,10,9,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,10,10,9,9,10,10,)
const IS15_3 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,7,6,7,6,5,5,9,8,10,10,8,9,9,8,10,10,8,9,5,5,6,7,6,7,)

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


# Primary invariants for NBody=5and deg=6
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
# Irreducible secondaries for NBody=5and deg=6
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
const IS15 = Val((IS15_1,IS15_2,IS15_3,))
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


function invariants(x1::SVector{10, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   x4 = x3.*x1
   x5 = x4.*x1
   x6 = x5.*x1
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NBody=5and deg=6
 # : what goes in the function for the evaluation
P1 = fpoly((x1,) , Main.P1)
P2 = fpoly((x1,x1,) , Main.P2)
P3 = fpoly((x2,) , Main.P3)
P4 = fpoly((x1,x1,x1,) , Main.P4)
P5 = fpoly((x3,) , Main.P5)
P6 = fpoly((x1,x1,x1,x1,) , Main.P6)
P7 = fpoly((x4,) , Main.P7)
P8 = fpoly((x1,x1,x1,x1,x1,) , Main.P8)
P9 = fpoly((x5,) , Main.P9)
P10 = fpoly((x6,) , Main.P10)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for NBody=5and deg=6
 # : what goes in the function for the evaluation
IS1 = fpoly((x2,x1,) , Main.IS1)
IS2 = fpoly((x1,x1,x1,) , Main.IS2)
IS3 = fpoly((x3,x1,) , Main.IS3)
IS4 = fpoly((x2,x2,) , Main.IS4)
IS5 = fpoly((x2,x1,x1,) , Main.IS5)
IS6 = fpoly((x2,x1,x1,) , Main.IS6)
IS7 = fpoly((x2,x1,x1,) , Main.IS7)
IS8 = fpoly((x4,x1,) , Main.IS8)
IS9 = fpoly((x3,x2,) , Main.IS9)
IS10 = fpoly((x3,x1,x1,) , Main.IS10)
IS11 = fpoly((x2,x2,x1,) , Main.IS11)
IS12 = fpoly((x3,x1,x1,) , Main.IS12)
IS13 = fpoly((x2,x2,x1,) , Main.IS13)
IS14 = fpoly((x3,x1,x1,) , Main.IS14)
IS15 = fpoly((x2,x2,x1,) , Main.IS15)
IS16 = fpoly((x5,x1,) , Main.IS16)
IS17 = fpoly((x4,x2,) , Main.IS17)
IS18 = fpoly((x3,x3,) , Main.IS18)
IS19 = fpoly((x4,x1,x1,) , Main.IS19)
IS20 = fpoly((x3,x2,x1,) , Main.IS20)
IS21 = fpoly((x2,x2,x2,) , Main.IS21)
IS22 = fpoly((x4,x1,x1,) , Main.IS22)
IS23 = fpoly((x3,x2,x1,) , Main.IS23)
IS24 = fpoly((x4,x1,x1,) , Main.IS24)
IS25 = fpoly((x3,x1,x1,x1,) , Main.IS25)
IS26 = fpoly((x2,x2,x1,x1,) , Main.IS26)
IS27 = fpoly((x3,x2,x1,) , Main.IS27)



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


return (@SVector [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,SEC13,SEC14,SEC15,SEC16,SEC17,SEC18,SEC19,SEC20,SEC21,SEC22,SEC23,SEC24,SEC25,SEC26,SEC27,SEC28,SEC29,SEC30,SEC31,])
 end
