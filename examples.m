/* The following functions PrimitiveDivisionPolynomial, IsIntegrallyDefined, MakeIntegrallyDefined and TorsionField 
are taken from galrep.m which is a part of github code associated to  "Computing images of Galois representations attached to elliptic curves", Forum of Mathematics, Sigma 4 (2016) e4 (79 pages), http://dx.doi.org/10.1017/fms.2015.33
Copyright (c) 2015 by Andrew V. Sutherland */

PrimitiveDivisionPolynomial := function(E,m)
    local f;
    f:=DivisionPolynomial(E,m);
    for d in Divisors(m) do if d gt 1 and d lt m then f := ExactQuotient(f,$$(E,d)); end if; end for;
    return f;
end function;

// Magma really wants number fields to be defined by integral monic polynomials, so we make sure this happens
IsIntegrallyDefined := function(K)
    local f;
    if K eq Rationals() then return true; end if;
    if not IsAbsoluteField(K) then return false; end if;
    f := DefiningPolynomial(K);
    return IsMonic(f) and &and[c in Integers():c in Coefficients(f)];
end function;

// Redefines a number field so that it is defined in terms of the absolute minimal polynomial of a generator that is an algebraic integer
MakeIntegrallyDefined := function(K)
    local g;
    while not IsIntegrallyDefined(K) do
        g := SimpleExtension(K).1;
        f := MinimalPolynomial(g);
        g *:= &*PrimeDivisors(LCM([Denominator(c/LeadingCoefficient(f)):c in Coefficients(f)]));
        K:=NumberField(MinimalPolynomial(g));
    end while;
    return K;
end function;


// Returns a pair [P,Q] of independent generators for E[m] (the points P and Q will necessarily have order m), where E is an elliptic curve y^2=x^3+Ax+B with A,B in Q
// Be warned that this is painfully slow: unless the m-division field Q(E[m]) has very small degree you will need to be patient.
TorsionField := function(E,m)
    local C, K, L, EL, x1, y1, y1s, x2, y2, y2s, Q, P, S, phi, f, b, g;

    C := Coefficients(E);
    assert C[1] eq 0 and C[2] eq 0 and C[3] eq 0 and C[4] in Rationals() and C[5] in Rationals(); // To simplify matters, we require E to be in the form y^2=x^3+Ax+B with A,B in Q
    phi:=PrimitiveDivisionPolynomial(E,m);
    roots := Roots(phi);
    if #roots ne Degree(phi) then
        K:=SplittingField(phi);
        return $$(ChangeRing(E,MakeIntegrallyDefined(K)),m);
    end if;
    K:=BaseRing(E);
    L:=K;
    R<x>:=PolynomialRing(K);
    // Our first basis point P (of order m) will have x-coord equal to the first root of phi
    x1:=roots[1][1];
    f:=x^3+C[4]*x+C[5];
    y1s:=Evaluate(f,x1);
    b,y1:=IsSquare(y1s);  // this step is time-consuming
    // if y1 is not in L, extend L so that it is
    if not b then L := NumberField(x^2-y1s); end if;
    if L ne Rationals() and not IsAbsoluteField(L) then L:=AbsoluteField(L); end if;
    return MakeIntegrallyDefined(L);
end function;

PrimitiveDivisionPolynomial := function(E,m)
    local f;
    f:=DivisionPolynomial(E,m);
    for d in Divisors(m) do if d gt 1 and d lt m then f := ExactQuotient(f,$$(E,d)); end if; end for;
    return f;
end function;

// Magma really wants number fields to be defined by integral monic polynomials, so we make sure this happens
IsIntegrallyDefined := function(K)
    local f;
    if K eq Rationals() then return true; end if;
    if not IsAbsoluteField(K) then return false; end if;
    f := DefiningPolynomial(K);
    return IsMonic(f) and &and[c in Integers():c in Coefficients(f)];
end function;

// Redefines a number field so that it is defined in terms of the absolute minimal polynomial of a generator that is an algebraic integer
MakeIntegrallyDefined := function(K)
    local g;
    while not IsIntegrallyDefined(K) do
        g := SimpleExtension(K).1;
        f := MinimalPolynomial(g);
        g *:= &*PrimeDivisors(LCM([Denominator(c/LeadingCoefficient(f)):c in Coefficients(f)]));
        K:=NumberField(MinimalPolynomial(g));
    end while;
    return K;
end function;


// Returns a pair [P,Q] of independent generators for E[m] (the points P and Q will necessarily have order m), where E is an elliptic curve y^2=x^3+Ax+B with A,B in Q
// Be warned that this is painfully slow: unless the m-division field Q(E[m]) has very small degree you will need to be patient.
TorsionField := function(E,m)
    local C, K, L, EL, x1, y1, y1s, x2, y2, y2s, Q, P, S, phi, f, b, g;

    C := Coefficients(E);
    assert C[1] eq 0 and C[2] eq 0 and C[3] eq 0 and C[4] in Rationals() and C[5] in Rationals(); // To simplify matters, we require E to be in the form y^2=x^3+Ax+B with A,B in Q
    phi:=PrimitiveDivisionPolynomial(E,m);
    roots := Roots(phi);
    if #roots ne Degree(phi) then
        K:=SplittingField(phi);
        return $$(ChangeRing(E,MakeIntegrallyDefined(K)),m);
    end if;
    K:=BaseRing(E);
    L:=K;
    R<x>:=PolynomialRing(K);
    // Our first basis point P (of order m) will have x-coord equal to the first root of phi
    x1:=roots[1][1];
    f:=x^3+C[4]*x+C[5];
    y1s:=Evaluate(f,x1);
    b,y1:=IsSquare(y1s);  // this step is time-consuming
    // if y1 is not in L, extend L so that it is
    if not b then L := NumberField(x^2-y1s); end if;
    if L ne Rationals() and not IsAbsoluteField(L) then L:=AbsoluteField(L); end if;
    return MakeIntegrallyDefined(L);
end function;


/* 450.c1 */
D:=CremonaDatabase();
E:=EllipticCurve(D,"450d4");
E := WeierstrassModel(E);
R<x> := PolynomialRing(Rationals());

K3 := TorsionField(E,3);
K5 := TorsionField(E,5);

// Check that Q(E[3]) contains sqrt(5).
has_root, _ := HasRoot(x^2 - 5, K3);
assert has_root;

// Check that Q(E[5]) contains sqrt(-3).
has_root, _ := HasRoot(x^2 + 3, K5);
assert has_root;

//E(F_7) eq 10

E7:=Curve(Reduction(E,7));
assert #Points(E7) eq 10;
