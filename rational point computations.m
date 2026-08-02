//14.96.2-14.e.2.1

P<x>:=PolynomialRing(Rationals());
f:= - 7*x^5 + 63*x^4 - 91*x^3 + 28*x^2 + 7*x; // model taken from LMFDB
C:=HyperellipticCurve(f);
J:=Jacobian(C);
#Chabauty0(J); // 3 points, all are cusps

//14.96.2-14.e.1.1


f:= 7*x^5 + 63*x^4 + 91*x^3 + 28*x^2 - 7*x; // model taken from LMFDB
C:=HyperellipticCurve(f);
J:=Jacobian(C);
#Chabauty0(J); // 3 points, all are cusps

//15.96.1-15.a.1.5



f:=  x^3 + x^2 - 5*x + 2; 
g:= x+1; // model taken from LMFDB
C:=EllipticCurve(f,g);
Rank(C); // 0 true
#TorsionSubgroup(C); //8 points, 4 are cusps, 4 give rise to two non-CM j-invariants


//15.96.1-15.a.2.5



f:=  x^3 + x^2 - 5*x + 2; 
g:= x+1; // model taken from LMFDB
C:=EllipticCurve(f,g);
Rank(C); // 0 true
#TorsionSubgroup(C); //8 points, 4 are cusps, 4 give rise to two non-CM j-invariants

//30.96.1-15.a.1.7


f:=  x^3 + x^2 - 5*x + 2; 
g:= x+1; // model taken from LMFDB
C:=EllipticCurve(f,g);
Rank(C); // 0 true
#TorsionSubgroup(C); //8 points, 4 are cusps, 4 give rise to two non-CM j-invariants

//30.96.1-15.a.2.3



f:=  x^3 + x^2 - 5*x + 2; 
g:= x+1; // model taken from LMFDB
C:=EllipticCurve(f,g);
Rank(C); // 0 true
#TorsionSubgroup(C); //8 points, 4 are cusps, 4 give rise to two non-CM j-invariants
