module BA_3B 

using NBodyIPs.FastPolys 
using StaticArrays 
using BenchmarkTools: @btime 

import NBodyIPs.tdegrees 

const G_BA_3B = [
[ 1, 2, 3 ]
,[ 2, 1, 3 ]
,] 
simplex_permutations(x::SVector{3}) = [x[G_BA_3B[i]] for i=1:2]
# Primary invariants for BA_3B  
 # : definitions at the beginning of the file 
const P1_1 = (1,2,) 

const P2_1 = (3,) 

const P3_1 = (1,2,) 

# Irreducible secondaries for group BA_3B
 # : definitions at the beginning of the file 

# Primary invariants for BA_3B  
 # : definitions of the types at the beginning of the file 
const P1 = Val((P1_1,)) 
const P2 = Val((P2_1,)) 
const P3 = Val((P3_1,)) 
# Irreducible secondaries for group BA_3B
 # : definitions of the types at the beginning of the file 


function invariants(x1::SVector{3, T}) where {T}
   x2 = x1.*x1 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_3B  
 # : what goes in the function for the evaluation 
P1 = fpoly((x1,) , BA_3B.P1) 
P2 = fpoly((x1,) , BA_3B.P2) 
P3 = fpoly((x2,) , BA_3B.P3) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_3B
 # : what goes in the function for the evaluation 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1


return (@SVector [P1,P2,P3,]), (@SVector [SEC1,])
 end



function invariants_d(x1::SVector{3, T}) where {T}
   x2 = x1.*x1 

   dx1 = @SVector ones(3)
   dx2 = 2 * x1 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_3B  
 # : what goes in the function for the derivatives 
dP1 = fpoly_d((x1,),(dx1,) , BA_3B.P1) 
dP2 = fpoly_d((x1,),(dx1,) , BA_3B.P2) 
dP3 = fpoly_d((x2,),(dx2,) , BA_3B.P3) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_3B
 # : what goes in the function for the evaluation 


# Irreducible secondaries for group BA_3B
 # : what goes in the function for the derivatives 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


dSEC1   = @SVector zeros(3) 


return (dP1,dP2,dP3,), (dSEC1,)
 end



function invariants_ed(x1::SVector{3, T}) where {T}
   x2 = x1.*x1 

   dx1 = @SVector ones(3)
   dx2 = 2 * x1 
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_3B  
 # : what goes in the function for the evaluation and derivatives 
P1, dP1 = fpoly_ed((x1,),(dx1,) , BA_3B.P1) 
P2, dP2 = fpoly_ed((x1,),(dx1,) , BA_3B.P2) 
P3, dP3 = fpoly_ed((x2,),(dx2,) , BA_3B.P3) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_3B
 # : what goes in the function for the evaluation and derivatives 



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1


dSEC1   = @SVector zeros(3) 


return (@SVector [P1,P2,P3,]), (@SVector [SEC1,]), (@SVector [dP1,dP2,dP3,]), (@SVector [dSEC1,])
 end 

end