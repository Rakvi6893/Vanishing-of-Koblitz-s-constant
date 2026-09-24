ChangeDirectory("/Users/Jacob/OpenImage-master"); // Change this line as needed.
load "main/FindOpenImage.m";

PrimObsLevels := [2,3,5,6,7,10,14,15,30];

function PrimObs(E)
    G := FindOpenImage(E);
    N := #BaseRing(G);
    S := {};

    for m in PrimObsLevels do
        if N mod m eq 0 and not exists{n : n in S | m mod n eq 0} then
            Gm := sub< GL(2, Integers(m)) | Generators(G) >;
            C := ConjugacyClasses(Gm);

            if &and[not IsUnit(1 - Trace(c[3]) + Determinant(c[3])) : c in C] then
                Include(~S,m);
            end if;
        end if;
    end for;

    return S;
end function;


// LMFDB: 14400.ef4
E := EllipticCurve([0, 0, 0, 12660, -248560]);
PrimObs(E); // { 30 }

// LMFDB: 66.c3
E := EllipticCurve([1, 0, 0, -45, 81]);
PrimObs(E); // { 2, 5 }