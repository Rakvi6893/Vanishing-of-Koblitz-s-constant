//14.96.2-14.e.2.1

P<x>:=PolynomialRing(Rationals());
f:= - 7*x^5 + 63*x^4 - 91*x^3 + 28*x^2 + 7*x; // model taken from LMFDB
C:=HyperellipticCurve(f);
J:=Jacobian(C);
#Chabauty0(J); // 3 points, all are cusps

//14.96.2-14.e.1.1

P<x>:=PolynomialRing(Rationals());
f:= 7*x^5 + 63*x^4 + 91*x^3 + 28*x^2 - 7*x; // model taken from LMFDB
C:=HyperellipticCurve(f);
J:=Jacobian(C);
#Chabauty0(J); // 3 points, all are cusps

