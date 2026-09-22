// Rule out congruence obstructions for the non-CM j-invariants
// arising from X_0(11), X_{S_4}(13), X_0(17), and X_0(37).

ChangeDirectory("/Users/Jacob/OpenImage-master"); // Change this line as needed.
load "main/FindOpenImage.m";

Q := Rationals();

// Non-CM rational j-invariants on the four modular curves.
JData := [
    <"X0(11)",    Q!(-11^2)>,
    <"X0(11)",    Q!(-11*131^3)>,

    <"XS4(13)",   Q!(2^4*5*13^4*17^3) / 3^13>,
    <"XS4(13)",   Q!(-2^12*5^3*11*13^4) / 3^13>,
    <"XS4(13)",   Q!(2^18*3^3*13^4*127^3*139^3*157^3*283^3*929) /
                       (5^13*61^13)>,

    <"X0(17)",    Q!(-17*373^3) / 2^17>,
    <"X0(17)",    Q!(-17^2*101^3) / 2>,

    <"X0(37)",    Q!(-7*11^3)>,
    <"X0(37)",    Q!(-7*137^3*2083^3)>
];

// Find A in G such that both det(I-A) and det(I+A) are units.
// For a 2 by 2 matrix A,
// det(I-A) = 1 - tr(A) + det(A),
// det(I+A) = 1 + tr(A) + det(A).
function RulingOutMatrix(G)
    A := Random(G);

    while not (IsUnit(1 - Trace(A) + Determinant(A)) and IsUnit(1 + Trace(A) + Determinant(A))) do
        A := Random(G);
    end while;

    return A;
end function;

Results := [* *];

for d in JData do
    name,j := Explode(d);

    E := MinimalModel(EllipticCurveWithjInvariant(j));
    assert jInvariant(E) eq j;

    G := FindOpenImage(E);
    N := #BaseRing(G);

    A := RulingOutMatrix(G);

    assert IsUnit(1 - Trace(A) + Determinant(A));
    assert IsUnit(1 + Trace(A) + Determinant(A));

    Append(~Results,<name,j,E,N,A>);

    printf "\n%o\n",name;
    printf "  j = %o\n",j;
    printf "  E = %o\n",E;
    printf "  N = %o\n",N;
    printf "  A =\n%o\n",A;
end for;



// Output (The matrix A will differ with each run due to randomization)
/*
Loading "RulingOutj.m"
Loading "main/FindOpenImage.m"
Loading "precomputation/ComputeFrobData.m"
Loading "main/GL2GroupTheory.m"
Loading "main/ModularCurves.m"
Loading "data-files/cyclic_invariant_polynomials.m"

X0(11)
  j = -121
  E = Elliptic Curve defined by y^2 + x*y = x^3 - 4661*x + 478402 over Rational
Field
  N = 3784
  A =
[ 464 1891]
[3465  301]

X0(11)
  j = -24729001
  E = Elliptic Curve defined by y^2 + x*y = x^3 - 12740962961911*x +
17505182394367463082 over Rational Field
  N = 57328744
  A =
[ 2124289 16590539]
[37269749 32966376]

XS4(13)
  j = 11225615440/1594323
  E = Elliptic Curve defined by y^2 = x^3 + x^2 - 123812476704729388*x +
14566237418944333871353268 over Rational Field
  N = 156
  A =
[137  11]
[ 65 136]

XS4(13)
  j = -160855552000/1594323
  E = Elliptic Curve defined by y^2 = x^3 + x^2 - 133858561551881833*x +
19010999804442441926361713 over Rational Field
  N = 78
  A =
[30 71]
[67  9]

XS4(13)
  j = 90616364985637924505590372621162077487104/1976504973537020943085705566406\
25
  E = Elliptic Curve defined by y^2 = x^3 -
4726867479322454444344535106724398017409848858810077158596991601303554390512*x +
1250856343289620870267647488715146101994341147279685474086331520166772951551293\
33769534583372313521695249290062734 over Rational Field
  N = 864370
  A =
[122792 469831]
[582371 458685]

X0(17)
  j = -882216989/131072
  E = Elliptic Curve defined by y^2 + x*y = x^3 - 509438683830771638*x +
156894294411435093981033892 over Rational Field
  N = 3776953240
  A =
[2722759850  774168263]
[3670661437 2415990451]

X0(17)
  j = -297756989/2
  E = Elliptic Curve defined by y^2 + x*y = x^3 - 46177215388802138*x +
3819374496993724449276142 over Rational Field
  N = 530003560
  A =
[316632366  50500989]
[ 55999887 323800369]

X0(37)
  j = -9317
  E = Elliptic Curve defined by y^2 + x*y = x^3 - 53597013*x + 164434147142 over
Rational Field
  N = 2678060
  A =
[1609523 1885183]
[1487881  398482]

X0(37)
  j = -162677523113838677
  E = Elliptic Curve defined by y^2 + x*y = x^3 -
13783321107528019849809024528304263*x + 622843482793159382014668369334075963495\
568880424142 over Rational Field
  N = 266635745093276020
  A =
[161931452026746340 164768536687949139]
[ 36698270818365921 206670171626697099]
*/