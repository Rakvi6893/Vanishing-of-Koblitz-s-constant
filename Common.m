// Common data and helper routines for the congruence-obstruction computations.
//
// The scripts that load this file should first change to the directory
// containing Modular.spec and this file.

ChangeDirectory("/Users/jacob/Modular-main"); // Change this line as needed.
AttachSpec("Modular.spec");

// Compute F1(G), using Magma's conjugacy classes and their sizes.
// For a 2 by 2 matrix A, det(I-A) = 1 - tr(A) + det(A).
function F1(G)
    C := ConjugacyClasses(G);
    return (&+[c[2] : c in C | not IsUnit(1 - Trace(c[3]) + Determinant(c[3]))]) / #G;
end function;

// Prime level Galois images of level <= 13; see Table 3 
// of Sutherland's "Computing images of Galois representations 
// attached to elliptic curves" article. 
// Each entry is <prime, index, generators>.
Data := AssociativeArray();

Data["2Cs"] := <2, 6, []>;
Data["2B"]  := <2, 3, [[1,1,0,1]]>;
Data["2Cn"] := <2, 2, [[0,1,1,1]]>;

Data["3Cs.1.1"] := <3, 24, [[1,0,0,2]]>;
Data["3Cs"]     := <3, 12, [[2,0,0,2], [1,0,0,2]]>;
Data["3B.1.1"]  := <3, 8,  [[1,0,0,2], [1,1,0,1]]>;
Data["3B.1.2"]  := <3, 8,  [[2,0,0,1], [1,1,0,1]]>;
Data["3Ns"]     := <3, 6,  [[2,0,0,2], [0,2,1,0], [1,0,0,2]]>;
Data["3B"]      := <3, 4,  [[2,0,0,2], [1,0,0,2], [1,1,0,1]]>;
Data["3Nn"]     := <3, 3,  [[1,0,0,2], [2,1,2,2]]>;

Data["5Cs.1.1"] := <5, 120, [[1,0,0,2]]>;
Data["5Cs.1.3"] := <5, 120, [[3,0,0,4]]>;
Data["5Cs.4.1"] := <5, 60,  [[4,0,0,4], [1,0,0,2]]>;
Data["5Ns.2.1"] := <5, 30,  [[2,0,0,3], [0,1,3,0]]>;
Data["5Cs"]     := <5, 30,  [[2,0,0,3], [1,0,0,2]]>;
Data["5B.1.1"]  := <5, 24,  [[1,0,0,2], [1,1,0,1]]>;
Data["5B.1.2"]  := <5, 24,  [[2,0,0,1], [1,1,0,1]]>;
Data["5B.1.4"]  := <5, 24,  [[4,0,0,3], [1,1,0,1]]>;
Data["5B.1.3"]  := <5, 24,  [[3,0,0,4], [1,1,0,1]]>;
Data["5Ns"]     := <5, 15,  [[0,4,1,0], [2,0,0,3], [1,0,0,2]]>;
Data["5B.4.1"]  := <5, 12,  [[4,0,0,4], [1,0,0,2], [1,1,0,1]]>;
Data["5B.4.2"]  := <5, 12,  [[4,0,0,4], [2,0,0,1], [1,1,0,1]]>;
Data["5Nn"]     := <5, 10,  [[1,0,0,4], [2,3,4,2]]>;
Data["5B"]      := <5, 6,   [[2,0,0,3], [1,0,0,2], [1,1,0,1]]>;
Data["5S4"]     := <5, 5,   [[0,3,3,4], [2,0,0,2], [3,0,4,4]]>;

Data["7Ns.2.1"] := <7, 112, [[2,0,0,4], [0,1,4,0]]>;
Data["7Ns.3.1"] := <7, 56,  [[3,0,0,5], [0,1,4,0]]>;
Data["7B.1.1"]  := <7, 48,  [[1,0,0,3], [1,1,0,1]]>;
Data["7B.1.3"]  := <7, 48,  [[3,0,0,1], [1,1,0,1]]>;
Data["7B.1.2"]  := <7, 48,  [[2,0,0,5], [1,1,0,1]]>;
Data["7B.1.5"]  := <7, 48,  [[5,0,0,2], [1,1,0,1]]>;
Data["7B.1.6"]  := <7, 48,  [[6,0,0,4], [1,1,0,1]]>;
Data["7B.1.4"]  := <7, 48,  [[4,0,0,6], [1,1,0,1]]>;
Data["7Ns"]     := <7, 28,  [[0,6,1,0], [3,0,0,5], [1,0,0,3]]>;
Data["7B.6.1"]  := <7, 24,  [[6,0,0,6], [1,0,0,3], [1,1,0,1]]>;
Data["7B.6.3"]  := <7, 24,  [[6,0,0,6], [3,0,0,1], [1,1,0,1]]>;
Data["7B.6.2"]  := <7, 24,  [[6,0,0,6], [2,0,0,5], [1,1,0,1]]>;
Data["7Nn"]     := <7, 21,  [[1,0,0,6], [2,5,4,2]]>;
Data["7B.2.1"]  := <7, 16,  [[2,0,0,4], [1,0,0,3], [1,1,0,1]]>;
Data["7B.2.3"]  := <7, 16,  [[2,0,0,4], [3,0,0,1], [1,1,0,1]]>;
Data["7B"]      := <7, 8,   [[3,0,0,5], [1,0,0,3], [1,1,0,1]]>;

Data["13S4"]    := <13, 91, [[3,0,12,9], [2,0,0,2], [9,5,0,6]]>;
Data["13B.3.1"] := <13, 56, [[3,0,0,9], [1,0,0,2], [1,1,0,1]]>;
Data["13B.3.2"] := <13, 56, [[3,0,0,9], [2,0,0,1], [1,1,0,1]]>;
Data["13B.3.4"] := <13, 56, [[3,0,0,9], [4,0,0,7], [1,1,0,1]]>;
Data["13B.3.7"] := <13, 56, [[3,0,0,9], [7,0,0,4], [1,1,0,1]]>;
Data["13B.5.1"] := <13, 42, [[5,0,0,8], [1,0,0,2], [1,1,0,1]]>;
Data["13B.5.2"] := <13, 42, [[5,0,0,8], [2,0,0,1], [1,1,0,1]]>;
Data["13B.5.4"] := <13, 42, [[5,0,0,8], [4,0,0,7], [1,1,0,1]]>;
Data["13B.4.1"] := <13, 28, [[4,0,0,10], [1,0,0,2], [1,1,0,1]]>;
Data["13B.4.2"] := <13, 28, [[4,0,0,10], [2,0,0,1], [1,1,0,1]]>;
Data["13B"]     := <13, 14, [[2,0,0,7], [1,0,0,2], [1,1,0,1]]>;

// Also include full image.
for p in [2,3,5,7,13] do
    Data[IntegerToString(p) cat "GL"] := <p, 1, [Eltseq(g) : g in Generators(GL(2,Integers(p)))]>;
end for;

Groups := AssociativeArray();
Alpha := AssociativeArray();
InBorel := AssociativeArray();

for label in Keys(Data) do
    p, index, gens := Explode(Data[label]);
    A := GL2Ambient(p);
    H := sub<A | gens>;
    assert GL2Index(H) eq index;
    assert GL2DeterminantIndex(H) eq 1;
    Groups[label] := H;
    Alpha[label] := F1(H);
    // A reducible subgroup of GL2 is conjugate into a Borel.
    InBorel[label] := not IsIrreducible(GModule(ChangeRing(H,GF(p))));
end for;

// The prime-by-prime restrictions established in the article.
function AllowedLabels(m,p)
    ell := Max(PrimeDivisors(m));
    labels := [];

    for label in Sort([a : a in Keys(Data) | Data[a][1] eq p]) do
        full := Data[label][2] eq 1;

        if Alpha[label] eq 1 then
            continue;
        end if;

        // Corollary C:full-prime-support.
        if p eq ell and p ge 5 and full then
            continue;
        end if;

        // Proposition P:13-large.
        if ell eq 13 then
            if p in {2,3} and not full then
                continue;
            end if;
            if p eq 13 and not InBorel[label] then
                continue;
            end if;
        end if;

        // Proposition P:7-large.
        if ell eq 7 and m mod 5 eq 0 then
            if p eq 2 and not full then
                continue;
            end if;
            if p eq 3 and InBorel[label] then
                continue;
            end if;
            if p in {5,7} and (full or InBorel[label]) then
                continue;
            end if;
        end if;

        Append(~labels,label);
    end for;
    return labels;
end function;

// Proposition P:level-210.  Level210.m independently recomputes
// this list and checks that it agrees with these four tuples.
Level210Tuples := {
    ["2GL","3GL","5Nn","7Nn"],
    ["2GL","3GL","5S4","7Nn"],
    ["2GL","3GL","5S4","7Ns"],
    ["2GL","3Nn","5Nn","7Nn"]
};

// Construct the CRT product of prime-level groups with pairwise distinct
// prime levels.  The returned subgroup lies in GL_2(Z/mZ), where m is the
// product of the prime levels attached to the labels.
function CRTProduct(labels)
    primes := [Data[label][1] : label in labels];
    assert #Seqset(primes) eq #labels;

    m := &*primes;
    H := &meet [GL2Lift(Groups[label],m) : label in labels];
    assert #H eq &*[#Groups[label] : label in labels];
    return H;
end function;
