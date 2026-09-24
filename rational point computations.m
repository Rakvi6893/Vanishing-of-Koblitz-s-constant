P<x>:=PolynomialRing(Rationals());

//10.24.1.b.1 and 10.24.1.b.2 (both have isomorphic models)
f:=  x^3 - x^2 - 1033*x - 12438; // model taken from LMFDB
C:=EllipticCurve(f);
Rank(C); // 0 true
#TorsionSubgroup(C);// 2 both are cusps

//14.48.2.e.1
f:= 7*x^5 + 63*x^4 + 91*x^3 + 28*x^2 - 7*x; // model taken from LMFDB
C:=HyperellipticCurve(f);
J:=Jacobian(C);
RankBounds(J); // Rank 0
#Chabauty0(J); // 3 points, all are cusps

//14.48.2.e.2
f:= - 7*x^5 + 63*x^4 - 91*x^3 + 28*x^2 + 7*x; // model taken from LMFDB
C:=HyperellipticCurve(f);
J:=Jacobian(C);
RankBounds(J); // Rank 0
#Chabauty0(J); // 3 points, all are cusps

//15.48.1.a.1 and 15.48.1.a.2 (which have isomorphic models)
f:=  x^3 + x^2 - 5*x + 2; g:= x+1; // model taken from LMFDB
C:=EllipticCurve(f,g);
Rank(C); // 0 true
#TorsionSubgroup(C); //8 points, 4 are cusps, 4 give rise to two non-CM j-invariants
