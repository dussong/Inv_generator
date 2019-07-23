
Kt<t11,t12,t21,t22>:=PolynomialRing(Rationals(),4);
I:=ideal<Kt | [t11^2+t12^2-1,t21^2+t22^2-1,t11*t21+t12*t22]>;
// this defines the orthogonal group
A:=Matrix([[t11,t12],[t21,t22]]);
// this defines the natural action
A:=TensorProduct(MatrixAlgebra(Kt,2)!1,A);
A;
// [t11 t12 0 0 0 0]
// [t21 t22 0 0 0 0]
// [ 0 0 t11 t12 0 0]
// [ 0 0 t21 t22 0 0]
// [ 0 0 0 0 t11 t12]
// [ 0 0 0 0 t21 t22]
G := PermutationGroup<6 | (1,3)(2,4), (1,5)(2,6)>;
// Bt:=Matrix([[0,1],[1,0]]);
// B:=TensorProduct(MatrixAlgebra(Kt,2)!1,Bt);
// this defines the action on 3 points
R:=InvariantRing(I,A: LinearlyReductive:=true);
// This only sets up the data structure
Kx<x11,x12,x21,x22,x31,x32>:=PolynomialRing(R);
