function invariants_check(x)
prim=zeros(3,1);
pv=zeros(0,1);

v=zeros(1,1);
prim[1]=x[1] + x[2]

prim[2]=x[3]

prim[3]=x[1]^2 + x[2]^2

v[1] = 1
return prim, v, pv

end
