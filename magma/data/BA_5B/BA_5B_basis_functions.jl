module BA_5B 

using NBodyIPs.FastPolys 
using StaticArrays 
using BenchmarkTools: @btime 

import NBodyIPs.tdegrees 

const G_BA_5B = [
[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
,[ 3, 1, 4, 2, 6, 10, 8, 7, 5, 9 ]
,[ 2, 4, 1, 3, 9, 5, 8, 7, 10, 6 ]
,[ 4, 1, 2, 3, 7, 9, 10, 5, 6, 8 ]
,[ 1, 3, 4, 2, 6, 7, 5, 10, 8, 9 ]
,[ 3, 4, 2, 1, 10, 8, 6, 9, 7, 5 ]
,[ 2, 1, 3, 4, 5, 8, 9, 6, 7, 10 ]
,[ 4, 2, 3, 1, 9, 10, 7, 8, 5, 6 ]
,[ 1, 4, 2, 3, 7, 5, 6, 9, 10, 8 ]
,[ 3, 2, 1, 4, 8, 6, 10, 5, 9, 7 ]
,[ 2, 3, 4, 1, 8, 9, 5, 10, 6, 7 ]
,[ 4, 3, 1, 2, 10, 7, 9, 6, 8, 5 ]
,[ 1, 3, 2, 4, 6, 5, 7, 8, 10, 9 ]
,[ 3, 4, 1, 2, 10, 6, 8, 7, 9, 5 ]
,[ 2, 1, 4, 3, 5, 9, 8, 7, 6, 10 ]
,[ 4, 2, 1, 3, 9, 7, 10, 5, 8, 6 ]
,[ 1, 4, 3, 2, 7, 6, 5, 10, 9, 8 ]
,[ 3, 2, 4, 1, 8, 10, 6, 9, 5, 7 ]
,[ 2, 3, 1, 4, 8, 5, 9, 6, 10, 7 ]
,[ 4, 3, 2, 1, 10, 9, 7, 8, 6, 5 ]
,[ 1, 2, 4, 3, 5, 7, 6, 9, 8, 10 ]
,[ 3, 1, 2, 4, 6, 8, 10, 5, 7, 9 ]
,[ 2, 4, 3, 1, 9, 8, 5, 10, 7, 6 ]
,[ 4, 1, 3, 2, 7, 10, 9, 6, 5, 8 ]
,] 
simplex_permutations(x::SVector{10}) = [x[G_BA_5B[i]] for i=1:24]
  # : definitions at the beginning of the file 
const BF1_1 = (1,2,3,4,) 

const BF2_1 = (5,9,6,8,7,10,) 

const BF3_1 = (1,2,3,4,) 

const BF4_1 = (1,2,1,2,1,3,) 
const BF4_2 = (2,4,3,3,4,4,) 

const BF5_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF5_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF6_1 = (5,9,6,8,7,10,) 

const BF7_1 = (1,2,3,4,) 

const BF8_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF8_2 = (2,4,1,3,4,3,1,2,3,2,1,4,) 

const BF9_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF9_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF10_1 = (1,1,1,2,) 
const BF10_2 = (2,2,3,3,) 
const BF10_3 = (3,4,4,4,) 

const BF11_1 = (5,9,6,8,7,10,5,9,6,8,7,10,) 
const BF11_2 = (1,2,3,2,1,4,2,4,1,3,4,3,) 

const BF12_1 = (5,9,6,8,7,10,) 

const BF13_1 = (1,2,3,4,) 

const BF14_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF14_2 = (2,4,1,3,4,3,1,2,3,2,1,4,) 

const BF15_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF15_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF16_1 = (1,2,1,2,1,3,) 
const BF16_2 = (2,4,3,3,4,4,) 

const BF17_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF17_2 = (2,1,1,3,2,1,1,2,3,1,1,2,) 
const BF17_3 = (3,4,4,4,4,3,3,3,4,2,2,4,) 

const BF18_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF18_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF19_1 = (5,9,6,8,7,10,5,9,6,8,7,10,) 
const BF19_2 = (1,2,3,2,1,4,2,4,1,3,4,3,) 

const BF20_1 = (5,9,6,8,7,10,) 

const BF21_1 = (1,2,3,4,) 

const BF22_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF22_2 = (2,4,1,3,4,3,1,2,3,2,1,4,) 

const BF23_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF23_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF24_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF24_2 = (2,4,1,3,4,3,1,2,3,2,1,4,) 

const BF25_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF25_2 = (2,1,1,3,2,1,1,2,3,1,1,2,) 
const BF25_3 = (3,4,4,4,4,3,3,3,4,2,2,4,) 

const BF26_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF26_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF27_1 = (1,2,1,2,1,3,2,2,3,1,1,1,) 
const BF27_2 = (2,4,3,3,4,4,4,3,4,3,2,4,) 
const BF27_3 = (3,1,4,4,2,1,3,1,2,2,4,3,) 

const BF28_1 = (5,9,6,8,7,10,5,9,6,8,7,10,) 
const BF28_2 = (1,2,3,2,1,4,2,4,1,3,4,3,) 

const BF29_1 = (5,9,6,8,7,10,5,9,6,8,7,10,) 
const BF29_2 = (1,2,3,2,1,4,2,4,1,3,4,3,) 

const BF30_1 = (5,9,6,8,7,10,) 

const BF31_1 = (1,2,3,4,) 

const BF32_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF32_2 = (2,4,1,3,4,3,1,2,3,2,1,4,) 

const BF33_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF33_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF34_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF34_2 = (2,4,1,3,4,3,1,2,3,2,1,4,) 

const BF35_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF35_2 = (2,1,1,3,2,1,1,2,3,1,1,2,) 
const BF35_3 = (3,4,4,4,4,3,3,3,4,2,2,4,) 

const BF36_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF36_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF37_1 = (1,2,1,2,1,3,) 
const BF37_2 = (2,4,3,3,4,4,) 

const BF38_1 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,) 
const BF38_2 = (2,4,1,3,4,3,1,2,3,2,1,4,3,4,1,2,4,2,1,3,2,3,1,4,) 
const BF38_3 = (3,1,4,4,2,1,3,3,4,1,2,2,2,1,4,4,3,1,2,2,4,1,3,3,) 

const BF39_1 = (1,2,3,2,1,4,2,4,1,3,4,3,) 
const BF39_2 = (5,9,6,8,7,10,5,9,6,8,7,10,) 

const BF40_1 = (1,1,1,2,) 
const BF40_2 = (2,2,3,3,) 
const BF40_3 = (3,4,4,4,) 

const BF41_1 = (5,9,6,8,7,10,5,9,6,8,7,10,) 
const BF41_2 = (1,2,3,2,1,4,2,4,1,3,4,3,) 

const BF42_1 = (5,9,6,8,7,10,5,9,6,8,7,10,) 
const BF42_2 = (1,2,3,2,1,4,2,4,1,3,4,3,) 

const BF43_1 = (5,9,6,8,7,10,) 


  # : definitions of the types at the beginning of the file 
const BF1 = Val((BF1_1,)) 
const BF2 = Val((BF2_1,)) 
const BF3 = Val((BF3_1,)) 
const BF4 = Val((BF4_1,BF4_2,)) 
const BF5 = Val((BF5_1,BF5_2,)) 
const BF6 = Val((BF6_1,)) 
const BF7 = Val((BF7_1,)) 
const BF8 = Val((BF8_1,BF8_2,)) 
const BF9 = Val((BF9_1,BF9_2,)) 
const BF10 = Val((BF10_1,BF10_2,BF10_3,)) 
const BF11 = Val((BF11_1,BF11_2,)) 
const BF12 = Val((BF12_1,)) 
const BF13 = Val((BF13_1,)) 
const BF14 = Val((BF14_1,BF14_2,)) 
const BF15 = Val((BF15_1,BF15_2,)) 
const BF16 = Val((BF16_1,BF16_2,)) 
const BF17 = Val((BF17_1,BF17_2,BF17_3,)) 
const BF18 = Val((BF18_1,BF18_2,)) 
const BF19 = Val((BF19_1,BF19_2,)) 
const BF20 = Val((BF20_1,)) 
const BF21 = Val((BF21_1,)) 
const BF22 = Val((BF22_1,BF22_2,)) 
const BF23 = Val((BF23_1,BF23_2,)) 
const BF24 = Val((BF24_1,BF24_2,)) 
const BF25 = Val((BF25_1,BF25_2,BF25_3,)) 
const BF26 = Val((BF26_1,BF26_2,)) 
const BF27 = Val((BF27_1,BF27_2,BF27_3,)) 
const BF28 = Val((BF28_1,BF28_2,)) 
const BF29 = Val((BF29_1,BF29_2,)) 
const BF30 = Val((BF30_1,)) 
const BF31 = Val((BF31_1,)) 
const BF32 = Val((BF32_1,BF32_2,)) 
const BF33 = Val((BF33_1,BF33_2,)) 
const BF34 = Val((BF34_1,BF34_2,)) 
const BF35 = Val((BF35_1,BF35_2,BF35_3,)) 
const BF36 = Val((BF36_1,BF36_2,)) 
const BF37 = Val((BF37_1,BF37_2,)) 
const BF38 = Val((BF38_1,BF38_2,BF38_3,)) 
const BF39 = Val((BF39_1,BF39_2,)) 
const BF40 = Val((BF40_1,BF40_2,BF40_3,)) 
const BF41 = Val((BF41_1,BF41_2,)) 
const BF42 = Val((BF42_1,BF42_2,)) 
const BF43 = Val((BF43_1,)) 






function invariants_ed(x1::SVector{10, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 
   x5 = x4.*x1 
   x6 = x5.*x1 

   dx1 = @SVector ones(10)
   dx2 = 2 * x1 
   dx3 = 3 * x2 
   dx4 = 4 * x3 
   dx5 = 5 * x4 
   dx6 = 6 * x5 
   #------------------------------------------------
   # Basis functions
   #------------------------------------------------

  # : what goes in the function for the evaluation and derivatives 
BF1, dBF1 = fpoly_ed((x1,),(dx1,) , BA_5B.BF1) 
BF2, dBF2 = fpoly_ed((x1,),(dx1,) , BA_5B.BF2) 
BF3, dBF3 = fpoly_ed((x2,),(dx2,) , BA_5B.BF3) 
BF4, dBF4 = fpoly_ed((x1,x1,),(dx1,dx1,) , BA_5B.BF4) 
BF5, dBF5 = fpoly_ed((x1,x1,),(dx1,dx1,) , BA_5B.BF5) 
BF6, dBF6 = fpoly_ed((x2,),(dx2,) , BA_5B.BF6) 
BF7, dBF7 = fpoly_ed((x3,),(dx3,) , BA_5B.BF7) 
BF8, dBF8 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_5B.BF8) 
BF9, dBF9 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_5B.BF9) 
BF10, dBF10 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.BF10) 
BF11, dBF11 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_5B.BF11) 
BF12, dBF12 = fpoly_ed((x3,),(dx3,) , BA_5B.BF12) 
BF13, dBF13 = fpoly_ed((x4,),(dx4,) , BA_5B.BF13) 
BF14, dBF14 = fpoly_ed((x3,x1,),(dx3,dx1,) , BA_5B.BF14) 
BF15, dBF15 = fpoly_ed((x3,x1,),(dx3,dx1,) , BA_5B.BF15) 
BF16, dBF16 = fpoly_ed((x2,x2,),(dx2,dx2,) , BA_5B.BF16) 
BF17, dBF17 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.BF17) 
BF18, dBF18 = fpoly_ed((x2,x2,),(dx2,dx2,) , BA_5B.BF18) 
BF19, dBF19 = fpoly_ed((x3,x1,),(dx3,dx1,) , BA_5B.BF19) 
BF20, dBF20 = fpoly_ed((x4,),(dx4,) , BA_5B.BF20) 
BF21, dBF21 = fpoly_ed((x5,),(dx5,) , BA_5B.BF21) 
BF22, dBF22 = fpoly_ed((x4,x1,),(dx4,dx1,) , BA_5B.BF22) 
BF23, dBF23 = fpoly_ed((x4,x1,),(dx4,dx1,) , BA_5B.BF23) 
BF24, dBF24 = fpoly_ed((x3,x2,),(dx3,dx2,) , BA_5B.BF24) 
BF25, dBF25 = fpoly_ed((x3,x1,x1,),(dx3,dx1,dx1,) , BA_5B.BF25) 
BF26, dBF26 = fpoly_ed((x3,x2,),(dx3,dx2,) , BA_5B.BF26) 
BF27, dBF27 = fpoly_ed((x2,x2,x1,),(dx2,dx2,dx1,) , BA_5B.BF27) 
BF28, dBF28 = fpoly_ed((x3,x2,),(dx3,dx2,) , BA_5B.BF28) 
BF29, dBF29 = fpoly_ed((x4,x1,),(dx4,dx1,) , BA_5B.BF29) 
BF30, dBF30 = fpoly_ed((x5,),(dx5,) , BA_5B.BF30) 
BF31, dBF31 = fpoly_ed((x6,),(dx6,) , BA_5B.BF31) 
BF32, dBF32 = fpoly_ed((x5,x1,),(dx5,dx1,) , BA_5B.BF32) 
BF33, dBF33 = fpoly_ed((x5,x1,),(dx5,dx1,) , BA_5B.BF33) 
BF34, dBF34 = fpoly_ed((x4,x2,),(dx4,dx2,) , BA_5B.BF34) 
BF35, dBF35 = fpoly_ed((x4,x1,x1,),(dx4,dx1,dx1,) , BA_5B.BF35) 
BF36, dBF36 = fpoly_ed((x4,x2,),(dx4,dx2,) , BA_5B.BF36) 
BF37, dBF37 = fpoly_ed((x3,x3,),(dx3,dx3,) , BA_5B.BF37) 
BF38, dBF38 = fpoly_ed((x3,x2,x1,),(dx3,dx2,dx1,) , BA_5B.BF38) 
BF39, dBF39 = fpoly_ed((x3,x3,),(dx3,dx3,) , BA_5B.BF39) 
BF40, dBF40 = fpoly_ed((x2,x2,x2,),(dx2,dx2,dx2,) , BA_5B.BF40) 
BF41, dBF41 = fpoly_ed((x4,x2,),(dx4,dx2,) , BA_5B.BF41) 
BF42, dBF42 = fpoly_ed((x5,x1,),(dx5,dx1,) , BA_5B.BF42) 
BF43, dBF43 = fpoly_ed((x6,),(dx6,) , BA_5B.BF43) 


return (@SVector [BF1,BF2,BF3,BF4,BF5,BF6,BF7,BF8,BF9,BF10,BF11,BF12,BF13,BF14,BF15,BF16,BF17,BF18,BF19,BF20,BF21,BF22,BF23,BF24,BF25,BF26,BF27,BF28,BF29,BF30,BF31,BF32,BF33,BF34,BF35,BF36,BF37,BF38,BF39,BF40,BF41,BF42,BF43,]), (@SVector [dBF1,dBF2,dBF3,dBF4,dBF5,dBF6,dBF7,dBF8,dBF9,dBF10,dBF11,dBF12,dBF13,dBF14,dBF15,dBF16,dBF17,dBF18,dBF19,dBF20,dBF21,dBF22,dBF23,dBF24,dBF25,dBF26,dBF27,dBF28,dBF29,dBF30,dBF31,dBF32,dBF33,dBF34,dBF35,dBF36,dBF37,dBF38,dBF39,dBF40,dBF41,dBF42,dBF43,])
 end 

end