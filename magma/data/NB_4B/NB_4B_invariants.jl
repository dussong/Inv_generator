module NB_4B 

using NBodyIPs.FastPolys 
using StaticArrays 
using BenchmarkTools: @btime 

import NBodyIPs.tdegrees 

const G_NB_4B = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 2, 3, 1, 5, 6, 4 ]
,[ 2, 4, 6, 5, 1, 3 ]
,[ 3, 1, 2, 6, 4, 5 ]
,[ 1, 6, 5, 4, 3, 2 ]
,[ 3, 5, 4, 6, 2, 1 ]
,[ 2, 1, 3, 5, 4, 6 ]
,[ 3, 2, 1, 6, 5, 4 ]
,[ 4, 2, 6, 1, 5, 3 ]
,[ 1, 3, 2, 4, 6, 5 ]
,[ 6, 1, 5, 3, 4, 2 ]
,[ 5, 3, 4, 2, 6, 1 ]
,[ 5, 4, 3, 2, 1, 6 ]
,[ 6, 5, 1, 3, 2, 4 ]
,[ 1, 5, 6, 4, 2, 3 ]
,[ 4, 6, 2, 1, 3, 5 ]
,[ 3, 4, 5, 6, 1, 2 ]
,[ 2, 6, 4, 5, 3, 1 ]
,[ 4, 5, 3, 1, 2, 6 ]
,[ 5, 6, 1, 2, 3, 4 ]
,[ 5, 1, 6, 2, 4, 3 ]
,[ 6, 4, 2, 3, 1, 5 ]
,[ 4, 3, 5, 1, 6, 2 ]
,[ 6, 2, 4, 3, 5, 1 ]
,] 
simplex_permutations(x::SVector{6}) = [x[G_NB_4B[i]] for i=1:24]
# Primary invariants for NB_4B  
 # : definitions at the beginning of the file 
const P1_1 = (1,3,5,2,6,4,) 

const P2_1 = (1,3,5,2,6,4,) 

const P3_1 = (1,1,1,2,1,5,2,2,4,4,3,3,) 
const P3_2 = (2,3,5,3,6,6,4,6,6,5,5,4,) 

const P4_1 = (1,3,5,2,6,4,) 

const P5_1 = (1,1,2,3,) 
const P5_2 = (2,5,4,4,) 
const P5_3 = (3,6,6,5,) 

const P6_1 = (1,3,5,2,6,4,) 

# Irreducible secondaries for group NB_4B
 # : definitions at the beginning of the file 
const IS1_1 = (1,3,5,2,1,6,2,3,4,1,2,6,5,3,1,4,5,6,4,3,2,5,4,6,) 
const IS1_2 = (2,1,1,3,6,5,1,2,2,3,6,4,4,5,5,3,6,1,5,4,4,3,6,2,) 

const IS2_1 = (1,3,5,2,1,6,2,3,4,1,2,6,5,3,1,4,5,6,4,3,2,5,4,6,) 
const IS2_2 = (2,1,1,3,4,3,1,2,1,3,5,3,2,5,4,1,2,1,1,4,4,2,1,2,) 
const IS2_3 = (4,6,2,5,6,5,5,6,2,4,6,4,4,6,5,3,6,3,5,6,5,3,6,3,) 

const IS3_1 = (1,3,2,2,1,3,3,1,2,3,2,1,) 
const IS3_2 = (4,6,5,5,4,6,6,4,5,6,5,4,) 
const IS3_3 = (2,1,1,3,6,5,2,3,6,4,4,5,) 


# Primary invariants for NB_4B  
 # : definitions of the types at the beginning of the file 
const P1 = Val((P1_1,)) 
const P2 = Val((P2_1,)) 
const P3 = Val((P3_1,P3_2,)) 
const P4 = Val((P4_1,)) 
const P5 = Val((P5_1,P5_2,P5_3,)) 
const P6 = Val((P6_1,)) 
# Irreducible secondaries for group NB_4B
 # : definitions of the types at the beginning of the file 
const IS1 = Val((IS1_1,IS1_2,)) 
const IS2 = Val((IS2_1,IS2_2,IS2_3,)) 
const IS3 = Val((IS3_1,IS3_2,IS3_3,)) 


function invariants(x1::SVector{6, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B  
 # : what goes in the function for the evaluation 
P1 = fpoly((x1,) , NB_4B.P1) 
P2 = fpoly((x2,) , NB_4B.P2) 
P3 = fpoly((x1,x1,) , NB_4B.P3) 
P4 = fpoly((x3,) , NB_4B.P4) 
P5 = fpoly((x1,x1,x1,) , NB_4B.P5) 
P6 = fpoly((x4,) , NB_4B.P6) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B
 # : what goes in the function for the evaluation 
IS1 = fpoly((x2,x1,) , NB_4B.IS1) 
IS2 = fpoly((x2,x1,x1,) , NB_4B.IS2) 
IS3 = fpoly((x2,x2,x1,) , NB_4B.IS3) 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1
SEC3  = IS2
SEC4  = IS3
SEC5  = IS1^2
SEC6  = IS2*IS3


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,])
 end



function invariants_d(x1::SVector{6, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 

   dx1 = @SVector ones(6)
   dx2 = 2 * x1 
   dx3 = 3 * x2 
   dx4 = 4 * x3 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B  
 # : what goes in the function for the derivatives 
dP1 = fpoly_d((x1,),(dx1,) , NB_4B.P1) 
dP2 = fpoly_d((x2,),(dx2,) , NB_4B.P2) 
dP3 = fpoly_d((x1,x1,),(dx1,dx1,) , NB_4B.P3) 
dP4 = fpoly_d((x3,),(dx3,) , NB_4B.P4) 
dP5 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , NB_4B.P5) 
dP6 = fpoly_d((x4,),(dx4,) , NB_4B.P6) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B
 # : what goes in the function for the evaluation 
IS1 = fpoly((x2,x1,) , NB_4B.IS1) 
IS2 = fpoly((x2,x1,x1,) , NB_4B.IS2) 
IS3 = fpoly((x2,x2,x1,) , NB_4B.IS3) 


# Irreducible secondaries for group NB_4B
 # : what goes in the function for the derivatives 
dIS1 = fpoly_d((x2,x1,),(dx2,dx1,) , NB_4B.IS1) 
dIS2 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , NB_4B.IS2) 
dIS3 = fpoly_d((x2,x2,x1,),(dx2,dx2,dx1,) , NB_4B.IS3) 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


dSEC1   = @SVector zeros(6) 
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  =  + dIS1*2IS1
dSEC6  =  + dIS2*IS3 + dIS3*IS2


return (dP1,dP2,dP3,dP4,dP5,dP6,), (dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,)
 end



function invariants_ed(x1::SVector{6, T}) where {T}
   x2 = x1.*x1 
   x3 = x2.*x1 
   x4 = x3.*x1 

   dx1 = @SVector ones(6)
   dx2 = 2 * x1 
   dx3 = 3 * x2 
   dx4 = 4 * x3 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for NB_4B  
 # : what goes in the function for the evaluation and derivatives 
P1, dP1 = fpoly_ed((x1,),(dx1,) , NB_4B.P1) 
P2, dP2 = fpoly_ed((x2,),(dx2,) , NB_4B.P2) 
P3, dP3 = fpoly_ed((x1,x1,),(dx1,dx1,) , NB_4B.P3) 
P4, dP4 = fpoly_ed((x3,),(dx3,) , NB_4B.P4) 
P5, dP5 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , NB_4B.P5) 
P6, dP6 = fpoly_ed((x4,),(dx4,) , NB_4B.P6) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group NB_4B
 # : what goes in the function for the evaluation and derivatives 
IS1, dIS1 = fpoly_ed((x2,x1,),(dx2,dx1,) , NB_4B.IS1) 
IS2, dIS2 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , NB_4B.IS2) 
IS3, dIS3 = fpoly_ed((x2,x2,x1,),(dx2,dx2,dx1,) , NB_4B.IS3) 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1
SEC3  = IS2
SEC4  = IS3
SEC5  = IS1^2
SEC6  = IS2*IS3


dSEC1   = @SVector zeros(6) 
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  =  + dIS1*2IS1
dSEC6  =  + dIS2*IS3 + dIS3*IS2


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,]), (@SVector [dP1,dP2,dP3,dP4,dP5,dP6,]), (@SVector [dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,])
 end 

end