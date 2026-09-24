// Finite group-theoretic search for primitive congruence obstructions.

load "Common.m";

CandidateLevels := [6,10,14,15,21,26,30,35,39,42,70,78,105,210];

// Determine whether the coarse modular curve has finitely many Q-points.
// Compute the genus first; construct a model only in genus 0 or 1.
function CoarseCurveStatus(H)
    A := GL2Ambient(BaseRing(H));
    K := GL2IncludeNegativeOne(H);

    // Zywina's convention requires the transpose of the article's group.
    Kt := sub<A | {Transpose(k) : k in Generators(K)}>;
    N,Kt := GL2Level(Kt);
    g := GL2Genus(Kt : NoGenusData := true);

    if g ge 2 then
        return "finite",N,g;
    end if;

    try
        M := CreateModularCurveRec(Kt);
        point := exists{O : O in M`cusp_orbits | #O eq 1};
        if g eq 0 and point then
            return "infinite",N,g;
        end if;

        M := FindModelOfXG(M);
        if #M`psi eq 0 then
            // The model is the projective line itself.
            assert g eq 0 and M`model_degree eq 1;
            return "infinite",N,g;
        end if;

        R := Parent(M`psi[1]);
        X := Curve(ProjectiveSpace(BaseRing(R),Rank(R)-1),M`psi);

        if g eq 0 then
            if HasRationalPoint(Conic(X)) then
                return "infinite",N,g;
            end if;
            return "finite",N,g;
        end if;

        C := GenusOneModel(X);
        J := MinimalModel(Jacobian(C));
        lo,hi := RankBounds(J);
        if hi eq 0 then
            return "finite",N,g;
        end if;
        if not point and not IsLocallySoluble(C) then
            return "finite",N,g;
        end if;
        if lo gt 0 then
            // Finding a point proves infinitude; failing to find one
            // within this bound does not prove that X(Q) is empty.
            if point or #PointSearch(X,1000) gt 0 then
                return "infinite",N,g;
            end if;
        end if;
    catch e
        printf "  Coarse curve at level %o, genus %o: %o\n",N,g,e;
        assert false;
    end try;
    return "unresolved",N,g;
end function;

// Recursively descend through maximal subgroups, returning the updated
// lists of previously seen groups and recorded obstructions.
function Descend(H,labels,seen,rows)
    m := #BaseRing(H);
    P := Sort(PrimeDivisors(m));
    A := GL2Ambient(m);

    // (1) The exact local images must be unchanged.
    if exists{i : i in [1..#P] | GL2Project(H,P[i]) ne Groups[labels[i]]} then
        return seen,rows;
    end if;

    // (2) The determinant must be surjective.
    if GL2DeterminantIndex(H) ne 1 then
        return seen,rows;
    end if;

    C := ConjugacyClasses(H);

    // (3) Test the article's complex conjugation condition.
    if not exists{c : c in C | c[3]^2 eq Identity(H) and Trace(c[3]) eq 0 and Determinant(c[3]) eq -1} then
        return seen,rows;
    end if;

    // Obstruction tests are constant on conjugacy classes of H.

    // (4) No proper divisor may already give an obstruction.
    // It is enough to test the maximal proper divisors m/p.
    delta := [Integers()!(1 - Trace(c[3]) + Determinant(c[3])) : c in C];
    if exists{p : p in P | &and[GCD(a,m div p) gt 1 : a in delta]} then
        return seen,rows;
    end if;

    // Deduplicate in the full ambient GL2 after the cheap rejection tests.
    sig := Sort([<c[1],c[2]> : c in C]);
    if exists{s : s in seen | #s[1] eq #H and s[2] eq sig and IsConjugate(A,H,s[1])} then
        return seen,rows;
    end if;
    Append(~seen,<H,sig>);

    if &and[GCD(a,m) gt 1 : a in delta] then
        // The obstruction is primitive by (4).
        status,N,g := CoarseCurveStatus(H);
        Append(~rows,<labels,H,N,g,status>);
        printf "  Found candidate of genus %o: %o\n",g,status;

        // Stop only when the coarse curve has finitely many
        // rational points.  Otherwise continue below H.
        if status eq "finite" then
            return seen,rows;
        end if;
    end if;

    // Reverse preserves the previous order of visiting subgroups.
    for S in Reverse(MaximalSubgroups(H)) do
        seen,rows := Descend(S`subgroup,labels,seen,rows);
    end for;
    return seen,rows;
end function;

// At a fixed modulus, discard finite rows contained up to conjugacy
// in a strictly larger finite row.  Keep every other row unchanged.
function KeepMaximalFiniteRows(rows)
    result := [* *];
    for r in rows do
        keep := true;
        if r[5] eq "finite" then
            H := r[2];
            for s in rows do
                K := s[2];
                if s[5] ne "finite" or #K le #H then
                    continue;
                end if;
                contained,_ := GL2IsConjugateSubgroup(K,H);
                if contained then
                    keep := false;
                    break;
                end if;
            end for;
        end if;
        if keep then
            Append(~result,r);
        end if;
    end for;
    return result;
end function;

function SearchLevel(m)
    assert m in CandidateLevels;
    P := Sort(PrimeDivisors(m));
    choices := [AllowedLabels(m,p) : p in P];
    rows := [* *];

    for t in CartesianProduct(choices) do
        labels := [t[i] : i in [1..#P]];

        if m eq 210 and labels notin Level210Tuples then
            continue;
        end if;
        // Lemma L:Sum.  AllowedLabels already enforces Alpha < 1.
        if &+[Alpha[label] : label in labels] lt 1 then
            continue;
        end if;

        print "Starting", labels;
        // Via CRT, this is the product of the local groups.
        H0 := CRTProduct(labels);
        assert #BaseRing(H0) eq m;
        seen := [* *];
        seen,rows := Descend(H0,labels,seen,rows);
    end for;
    return KeepMaximalFiniteRows(rows);
end function;

// Data from LMFDB Beta with relevant modular curve data.
columns := ["label", "name", "level", "index", "genus", "generators"];
data := [
[*"6.24.0-3.a.1.1", 6, 24, 0, [[1, 0, 3, 5], [2, 3, 3, 2]]*],
[*"6.8.0-3.a.1.1", 6, 8, 0, [[2, 5, 3, 4], [5, 1, 0, 5]]*],
[*"6.8.0-3.a.1.2", 6, 8, 0, [[4, 3, 3, 2], [5, 5, 3, 4]]*],
[*"10.120.0-5.a.1.1", 10, 120, 0, [[1, 3, 5, 2], [9, 9, 0, 7]]*],
[*"10.24.0-5.a.2.2", 10, 24, 0, [[3, 8, 3, 5], [4, 7, 5, 8]]*],
[*"10.24.0-5.a.1.1", 10, 24, 0, [[1, 9, 0, 9], [8, 3, 7, 8]]*],
[*"14.48.0-7.a.1.1", 14, 48, 0, [[11, 9, 8, 5], [12, 1, 1, 12]]*],
[*"14.48.0-7.a.2.1", 14, 48, 0, [[6, 3, 7, 2], [13, 11, 5, 12]]*],
[*"10.24.1.b.1", 10, 24, 1, [[0, 7, 7, 3], [8, 9, 5, 6]]*],
[*"10.24.1.b.2", 10, 24, 1, [[5, 6, 3, 3], [7, 3, 6, 5]]*],
[*"14.48.2.e.2", 14, 48, 2, [[2, 7, 3, 6], [11, 1, 10, 9]]*],
[*"14.48.2.e.1", 14, 48, 2, [[4, 9, 1, 5], [9, 13, 4, 11]]*],
[*"15.72.3.e.1", 15, 72, 3, [[9, 7, 5, 6], [9, 10, 8, 6], [13, 6, 6, 13], [13, 6, 12, 2]]*],
[*"15.72.3.e.2", 15, 72, 3, [[0, 14, 11, 3], [3, 14, 2, 0], [9, 7, 5, 6], [14, 0, 0, 14]]*],
[*"15.48.1.a.1", 15, 48, 1, [[1, 5, 0, 8], [1, 11, 0, 14], [4, 2, 0, 7], [14, 4, 0, 2]]*],
[*"15.48.1.a.2", 15, 48, 1, [[1, 14, 0, 4], [7, 13, 0, 14], [8, 7, 0, 1], [14, 12, 0, 4]]*],
[*"21.144.7.b.1", 21, 144, 7, [[7, 6, 18, 2], [12, 1, 20, 0], [18, 2, 7, 6], [20, 9, 0, 1]]*],
[*"21.144.7.b.2", 21, 144, 7, [[0, 2, 20, 3], [10, 3, 18, 16], [11, 18, 3, 19], [12, 10, 20, 9]]*],
[*"21.96.3.a.1", 21, 96, 3, [[1, 4, 0, 10], [1, 15, 0, 8], [13, 4, 0, 5], [20, 5, 0, 13]]*],
[*"21.96.3.a.2", 21, 96, 3, [[8, 12, 0, 13], [10, 10, 0, 8], [11, 2, 0, 20], [16, 3, 0, 13]]*],
[*"30.72.5.bn.2", 30, 72, 5, [[1, 5, 14, 29], [2, 25, 11, 26], [28, 15, 27, 19], [29, 20, 1, 11]]*],
[*"30.72.5.bn.1", 30, 72, 5, [[9, 5, 5, 24], [17, 19, 4, 1], [18, 5, 11, 21], [25, 2, 17, 17]]*],
[*"42.144.10.t.2", 42, 144, 10, [[13, 28, 7, 41], [22, 35, 29, 26], [29, 29, 28, 41], [37, 26, 11, 29]]*],
[*"42.144.10.t.1", 42, 144, 10, [[0, 5, 19, 3], [6, 35, 7, 27], [37, 25, 8, 19], [39, 22, 4, 39]]*]
];
RecFormat := recformat<label,name,level,index,genus,generators,subgroup>;


function create_record(row)
    out := rec<RecFormat|label:=row[1],level:=row[2],index:=row[3],genus:=row[4],generators:=row[5]>;
    subgroup := out`level eq 1 select sub<GL(2,Integers())|> else sub<GL(2,Integers(out`level))|out`generators>;
    out`subgroup := subgroup;
    return out;
end function;

function make_data()
    return [create_record(row) : row in data];
end function;

IdentifyCrv := function(H, data)
    n, G := GL2Level(H);
    Gl := GL(2, Integers(n));
    idx := Index(Gl, G);
    for d in data do
        if d[2] eq n and d[3] eq idx then
            H := sub<Gl | d[5]>;
            if IsConjugate(Gl, G, H) then
                return d[1];
            end if;
        end if;
    end for;
    return "Not found";
end function;

Results := AssociativeArray();
for m in CandidateLevels do
    printf "\nSearching m = %o\n",m;
    time Results[m] := SearchLevel(m);

    for status in ["finite","infinite","unresolved"] do
        printf "  %o: %o recorded groups\n",status, #[r : r in Results[m] | r[5] eq status];
    end for;

    assert #[r : r in Results[m] | r[5] eq "unresolved"] eq 0;

    printf "\nModular curves for m = %o:\n",m;
    for r in Results[m] do
        labels := r[1];
        H := r[2];
        status := r[5];

        if status eq "infinite" then
            // Identify the fine modular curve X_H.
            printf "  %o: infinite, fine curve %o\n", labels, IdentifyCrv(H, data);
        elif status eq "finite" then
            // Identify the coarse modular curve X_<H,-I>.
            K := GL2IncludeNegativeOne(H);
            printf "  %o: finite, coarse curve %o\n", labels, IdentifyCrv(K, data);
        end if;
    end for;
end for;

// Output
/*
Searching m = 6
Starting [ 2Cn, 3B ]
Starting [ 2Cn, 3Cs ]
Starting [ 2GL, 3B ]
  Found candidate of genus 0: infinite
  Found candidate of genus 0: infinite
Starting [ 2GL, 3Cs ]
  Found candidate of genus 0: infinite
Starting [ 2GL, 3GL ]
Starting [ 2GL, 3Ns ]
Time: 0.630
  finite: 0 recorded groups
  infinite: 3 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 6:
  [ 2GL, 3B ]: infinite, fine curve 6.8.0-3.a.1.2
  [ 2GL, 3B ]: infinite, fine curve 6.8.0-3.a.1.1
  [ 2GL, 3Cs ]: infinite, fine curve 6.24.0-3.a.1.1

Searching m = 10
Starting [ 2GL, 5B ]
Starting [ 2GL, 5B.1.3 ]
  Found candidate of genus 1: finite
Starting [ 2GL, 5B.1.4 ]
  Found candidate of genus 1: finite
Starting [ 2GL, 5B.4.1 ]
  Found candidate of genus 0: infinite
Starting [ 2GL, 5B.4.2 ]
  Found candidate of genus 0: infinite
Starting [ 2GL, 5Cs ]
Starting [ 2GL, 5Cs.1.3 ]
  Found candidate of genus 5: finite
Starting [ 2GL, 5Cs.4.1 ]
  Found candidate of genus 0: infinite
Starting [ 2GL, 5Ns ]
Time: 0.240
  finite: 2 recorded groups
  infinite: 3 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 10:
  [ 2GL, 5B.1.3 ]: finite, coarse curve 10.24.1.b.2
  [ 2GL, 5B.1.4 ]: finite, coarse curve 10.24.1.b.1
  [ 2GL, 5B.4.1 ]: infinite, fine curve 10.24.0-5.a.1.1
  [ 2GL, 5B.4.2 ]: infinite, fine curve 10.24.0-5.a.2.2
  [ 2GL, 5Cs.4.1 ]: infinite, fine curve 10.120.0-5.a.1.1

Searching m = 14
Starting [ 2Cn, 7B.1.4 ]
Starting [ 2Cn, 7B.1.6 ]
Starting [ 2GL, 7B.1.2 ]
Starting [ 2GL, 7B.1.4 ]
  Found candidate of genus 2: finite
Starting [ 2GL, 7B.1.5 ]
Starting [ 2GL, 7B.1.6 ]
  Found candidate of genus 2: finite
Starting [ 2GL, 7B.2.1 ]
Starting [ 2GL, 7B.2.3 ]
Starting [ 2GL, 7B.6.1 ]
  Found candidate of genus 0: infinite
Starting [ 2GL, 7B.6.3 ]
  Found candidate of genus 0: infinite
Starting [ 2GL, 7Ns.2.1 ]
Time: 0.100
  finite: 2 recorded groups
  infinite: 2 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 14:
  [ 2GL, 7B.1.4 ]: finite, coarse curve 14.48.2.e.2
  [ 2GL, 7B.1.6 ]: finite, coarse curve 14.48.2.e.1
  [ 2GL, 7B.6.1 ]: infinite, fine curve 14.48.0-7.a.1.1
  [ 2GL, 7B.6.3 ]: infinite, fine curve 14.48.0-7.a.2.1

Searching m = 15
Starting [ 3B, 5B ]
Starting [ 3B, 5B.1.3 ]
Starting [ 3B, 5B.1.4 ]
Starting [ 3B, 5B.4.1 ]
  Found candidate of genus 1: finite
  Found candidate of genus 1: finite
  Found candidate of genus 1: finite
Starting [ 3B, 5B.4.2 ]
  Found candidate of genus 1: finite
  Found candidate of genus 1: finite
  Found candidate of genus 1: finite
Starting [ 3B, 5Cs ]
Starting [ 3B, 5Cs.1.3 ]
Starting [ 3B, 5Cs.4.1 ]
  Found candidate of genus 17: finite
  Found candidate of genus 17: finite
  Found candidate of genus 9: finite
Starting [ 3B, 5Ns ]
Starting [ 3Cs, 5B ]
Starting [ 3Cs, 5B.1.3 ]
Starting [ 3Cs, 5B.1.4 ]
Starting [ 3Cs, 5B.4.1 ]
  Found candidate of genus 9: finite
  Found candidate of genus 5: finite
Starting [ 3Cs, 5B.4.2 ]
  Found candidate of genus 9: finite
  Found candidate of genus 5: finite
Starting [ 3Cs, 5Cs ]
Starting [ 3Cs, 5Cs.1.3 ]
Starting [ 3Cs, 5Cs.4.1 ]
  Found candidate of genus 73: finite
  Found candidate of genus 37: finite
Starting [ 3Cs, 5Ns ]
Starting [ 3GL, 5B.4.1 ]
Starting [ 3GL, 5B.4.2 ]
Starting [ 3GL, 5Cs.4.1 ]
Starting [ 3Ns, 5B ]
Starting [ 3Ns, 5B.1.3 ]
Starting [ 3Ns, 5B.1.4 ]
Starting [ 3Ns, 5B.4.1 ]
  Found candidate of genus 5: finite
  Found candidate of genus 3: finite
Starting [ 3Ns, 5B.4.2 ]
  Found candidate of genus 5: finite
  Found candidate of genus 3: finite
Starting [ 3Ns, 5Cs ]
Starting [ 3Ns, 5Cs.1.3 ]
Starting [ 3Ns, 5Cs.4.1 ]
  Found candidate of genus 19: finite
  Found candidate of genus 37: finite
Time: 2.980
  finite: 4 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 15:
  [ 3B, 5B.4.1 ]: finite, coarse curve 15.48.1.a.1
  [ 3B, 5B.4.2 ]: finite, coarse curve 15.48.1.a.2
  [ 3Ns, 5B.4.1 ]: finite, coarse curve 15.72.3.e.1
  [ 3Ns, 5B.4.2 ]: finite, coarse curve 15.72.3.e.2

Searching m = 21
Starting [ 3B, 7B ]
Starting [ 3B, 7B.1.2 ]
Starting [ 3B, 7B.1.4 ]
Starting [ 3B, 7B.1.5 ]
Starting [ 3B, 7B.1.6 ]
Starting [ 3B, 7B.2.1 ]
Starting [ 3B, 7B.2.3 ]
Starting [ 3B, 7B.6.1 ]
  Found candidate of genus 5: finite
  Found candidate of genus 3: finite
  Found candidate of genus 5: finite
Starting [ 3B, 7B.6.2 ]
Starting [ 3B, 7B.6.3 ]
  Found candidate of genus 3: finite
  Found candidate of genus 5: finite
  Found candidate of genus 5: finite
Starting [ 3B, 7Ns.2.1 ]
Starting [ 3B, 7Ns.3.1 ]
Starting [ 3Cs, 7B ]
Starting [ 3Cs, 7B.1.2 ]
Starting [ 3Cs, 7B.1.4 ]
Starting [ 3Cs, 7B.1.5 ]
Starting [ 3Cs, 7B.1.6 ]
Starting [ 3Cs, 7B.2.1 ]
Starting [ 3Cs, 7B.2.3 ]
Starting [ 3Cs, 7B.6.1 ]
  Found candidate of genus 13: finite
  Found candidate of genus 25: finite
Starting [ 3Cs, 7B.6.2 ]
Starting [ 3Cs, 7B.6.3 ]
  Found candidate of genus 25: finite
  Found candidate of genus 13: finite
Starting [ 3Cs, 7Ns.2.1 ]
Starting [ 3Cs, 7Ns.3.1 ]
Starting [ 3GL, 7B.1.4 ]
Starting [ 3GL, 7B.1.6 ]
Starting [ 3GL, 7B.6.1 ]
Starting [ 3GL, 7B.6.3 ]
Starting [ 3Ns, 7B.1.4 ]
Starting [ 3Ns, 7B.1.6 ]
Starting [ 3Ns, 7B.2.1 ]
Starting [ 3Ns, 7B.2.3 ]
Starting [ 3Ns, 7B.6.1 ]
  Found candidate of genus 13: finite
  Found candidate of genus 7: finite
Starting [ 3Ns, 7B.6.3 ]
  Found candidate of genus 13: finite
  Found candidate of genus 7: finite
Starting [ 3Ns, 7Ns.2.1 ]
Time: 2.320
  finite: 4 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 21:
  [ 3B, 7B.6.1 ]: finite, coarse curve 21.96.3.a.1
  [ 3B, 7B.6.3 ]: finite, coarse curve 21.96.3.a.2
  [ 3Ns, 7B.6.1 ]: finite, coarse curve 21.144.7.b.1
  [ 3Ns, 7B.6.3 ]: finite, coarse curve 21.144.7.b.2

Searching m = 26
Starting [ 2GL, 13B.3.1 ]
Starting [ 2GL, 13B.3.2 ]
Time: 0.040
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 26:

Searching m = 30
Starting [ 2Cn, 3B, 5B ]
Starting [ 2Cn, 3B, 5B.1.3 ]
Starting [ 2Cn, 3B, 5B.1.4 ]
Starting [ 2Cn, 3B, 5B.4.1 ]
Starting [ 2Cn, 3B, 5B.4.2 ]
Starting [ 2Cn, 3B, 5Cs ]
Starting [ 2Cn, 3B, 5Cs.1.3 ]
Starting [ 2Cn, 3B, 5Cs.4.1 ]
Starting [ 2Cn, 3B, 5Nn ]
Starting [ 2Cn, 3B, 5Ns ]
Starting [ 2Cn, 3B, 5Ns.2.1 ]
Starting [ 2Cn, 3B, 5S4 ]
Starting [ 2Cn, 3Cs, 5B ]
Starting [ 2Cn, 3Cs, 5B.1.3 ]
Starting [ 2Cn, 3Cs, 5B.1.4 ]
Starting [ 2Cn, 3Cs, 5B.4.1 ]
Starting [ 2Cn, 3Cs, 5B.4.2 ]
Starting [ 2Cn, 3Cs, 5Cs ]
Starting [ 2Cn, 3Cs, 5Cs.1.3 ]
Starting [ 2Cn, 3Cs, 5Cs.4.1 ]
Starting [ 2Cn, 3Cs, 5Nn ]
Starting [ 2Cn, 3Cs, 5Ns ]
Starting [ 2Cn, 3Cs, 5Ns.2.1 ]
Starting [ 2Cn, 3Cs, 5S4 ]
Starting [ 2Cn, 3GL, 5B ]
Starting [ 2Cn, 3GL, 5B.1.3 ]
Starting [ 2Cn, 3GL, 5B.1.4 ]
Starting [ 2Cn, 3GL, 5B.4.1 ]
Starting [ 2Cn, 3GL, 5B.4.2 ]
Starting [ 2Cn, 3GL, 5Cs ]
Starting [ 2Cn, 3GL, 5Cs.1.3 ]
Starting [ 2Cn, 3GL, 5Cs.4.1 ]
Starting [ 2Cn, 3GL, 5Ns ]
Starting [ 2Cn, 3Nn, 5B ]
Starting [ 2Cn, 3Nn, 5B.1.3 ]
Starting [ 2Cn, 3Nn, 5B.1.4 ]
Starting [ 2Cn, 3Nn, 5B.4.1 ]
Starting [ 2Cn, 3Nn, 5B.4.2 ]
Starting [ 2Cn, 3Nn, 5Cs ]
Starting [ 2Cn, 3Nn, 5Cs.1.3 ]
Starting [ 2Cn, 3Nn, 5Cs.4.1 ]
Starting [ 2Cn, 3Ns, 5B ]
Starting [ 2Cn, 3Ns, 5B.1.3 ]
Starting [ 2Cn, 3Ns, 5B.1.4 ]
Starting [ 2Cn, 3Ns, 5B.4.1 ]
Starting [ 2Cn, 3Ns, 5B.4.2 ]
Starting [ 2Cn, 3Ns, 5Cs ]
Starting [ 2Cn, 3Ns, 5Cs.1.3 ]
Starting [ 2Cn, 3Ns, 5Cs.4.1 ]
Starting [ 2Cn, 3Ns, 5Nn ]
Starting [ 2Cn, 3Ns, 5Ns ]
Starting [ 2Cn, 3Ns, 5Ns.2.1 ]
Starting [ 2Cn, 3Ns, 5S4 ]
Starting [ 2GL, 3B, 5B ]
Starting [ 2GL, 3B, 5B.1.3 ]
  Found candidate of genus 5: finite
Starting [ 2GL, 3B, 5B.1.4 ]
  Found candidate of genus 5: finite
Starting [ 2GL, 3B, 5B.4.1 ]
  Found candidate of genus 1: finite
  Found candidate of genus 5: finite
  Found candidate of genus 5: finite
  Found candidate of genus 1: finite
  Found candidate of genus 5: finite
  Found candidate of genus 1: finite
Starting [ 2GL, 3B, 5B.4.2 ]
  Found candidate of genus 5: finite
  Found candidate of genus 5: finite
  Found candidate of genus 1: finite
  Found candidate of genus 1: finite
  Found candidate of genus 5: finite
  Found candidate of genus 1: finite
Starting [ 2GL, 3B, 5Cs ]
Starting [ 2GL, 3B, 5Cs.1.3 ]
  Found candidate of genus 29: finite
Starting [ 2GL, 3B, 5Cs.4.1 ]
  Found candidate of genus 9: finite
  Found candidate of genus 29: finite
  Found candidate of genus 29: finite
  Found candidate of genus 17: finite
  Found candidate of genus 17: finite
  Found candidate of genus 29: finite
Starting [ 2GL, 3B, 5Nn ]
Starting [ 2GL, 3B, 5Ns ]
Starting [ 2GL, 3B, 5Ns.2.1 ]
Starting [ 2GL, 3B, 5S4 ]
Starting [ 2GL, 3Cs, 5B ]
Starting [ 2GL, 3Cs, 5B.1.3 ]
  Found candidate of genus 17: finite
Starting [ 2GL, 3Cs, 5B.1.4 ]
  Found candidate of genus 17: finite
Starting [ 2GL, 3Cs, 5B.4.1 ]
  Found candidate of genus 17: finite
  Found candidate of genus 9: finite
  Found candidate of genus 17: finite
  Found candidate of genus 5: finite
Starting [ 2GL, 3Cs, 5B.4.2 ]
  Found candidate of genus 5: finite
  Found candidate of genus 17: finite
  Found candidate of genus 17: finite
  Found candidate of genus 9: finite
Starting [ 2GL, 3Cs, 5Cs ]
Starting [ 2GL, 3Cs, 5Cs.1.3 ]
  Found candidate of genus 97: finite
Starting [ 2GL, 3Cs, 5Cs.4.1 ]
  Found candidate of genus 97: finite
  Found candidate of genus 97: finite
  Found candidate of genus 73: finite
  Found candidate of genus 37: finite
Starting [ 2GL, 3Cs, 5Nn ]
Starting [ 2GL, 3Cs, 5Ns ]
Starting [ 2GL, 3Cs, 5Ns.2.1 ]
Starting [ 2GL, 3Cs, 5S4 ]
Starting [ 2GL, 3GL, 5B ]
Starting [ 2GL, 3GL, 5B.1.3 ]
Starting [ 2GL, 3GL, 5B.1.4 ]
Starting [ 2GL, 3GL, 5B.4.1 ]
Starting [ 2GL, 3GL, 5B.4.2 ]
Starting [ 2GL, 3GL, 5Cs ]
Starting [ 2GL, 3GL, 5Cs.1.3 ]
Starting [ 2GL, 3GL, 5Cs.4.1 ]
Starting [ 2GL, 3GL, 5Nn ]
Starting [ 2GL, 3GL, 5Ns ]
Starting [ 2GL, 3GL, 5Ns.2.1 ]
Starting [ 2GL, 3GL, 5S4 ]
Starting [ 2GL, 3Nn, 5B ]
Starting [ 2GL, 3Nn, 5B.1.3 ]
Starting [ 2GL, 3Nn, 5B.1.4 ]
Starting [ 2GL, 3Nn, 5B.4.1 ]
  Found candidate of genus 9: finite
  Found candidate of genus 5: finite
  Found candidate of genus 9: finite
  Found candidate of genus 5: finite
Starting [ 2GL, 3Nn, 5B.4.2 ]
  Found candidate of genus 9: finite
  Found candidate of genus 9: finite
  Found candidate of genus 5: finite
  Found candidate of genus 5: finite
Starting [ 2GL, 3Nn, 5Cs ]
Starting [ 2GL, 3Nn, 5Cs.1.3 ]
Starting [ 2GL, 3Nn, 5Cs.4.1 ]
  Found candidate of genus 49: finite
  Found candidate of genus 49: finite
  Found candidate of genus 25: finite
  Found candidate of genus 25: finite
Starting [ 2GL, 3Nn, 5Nn ]
Starting [ 2GL, 3Nn, 5Ns ]
Starting [ 2GL, 3Nn, 5Ns.2.1 ]
Starting [ 2GL, 3Nn, 5S4 ]
Starting [ 2GL, 3Ns, 5B ]
Starting [ 2GL, 3Ns, 5B.1.3 ]
  Found candidate of genus 9: finite
  Found candidate of genus 17: finite
Starting [ 2GL, 3Ns, 5B.1.4 ]
  Found candidate of genus 9: finite
  Found candidate of genus 17: finite
Starting [ 2GL, 3Ns, 5B.4.1 ]
  Found candidate of genus 17: finite
  Found candidate of genus 5: finite
  Found candidate of genus 3: finite
  Found candidate of genus 9: finite
  Found candidate of genus 9: finite
  Found candidate of genus 9: finite
Starting [ 2GL, 3Ns, 5B.4.2 ]
  Found candidate of genus 5: finite
  Found candidate of genus 17: finite
  Found candidate of genus 9: finite
  Found candidate of genus 3: finite
  Found candidate of genus 9: finite
  Found candidate of genus 9: finite
Starting [ 2GL, 3Ns, 5Cs ]
Starting [ 2GL, 3Ns, 5Cs.1.3 ]
  Found candidate of genus 49: finite
  Found candidate of genus 97: finite
Starting [ 2GL, 3Ns, 5Cs.4.1 ]
  Found candidate of genus 49: finite
  Found candidate of genus 97: finite
  Found candidate of genus 37: finite
  Found candidate of genus 19: finite
  Found candidate of genus 49: finite
  Found candidate of genus 49: finite
Starting [ 2GL, 3Ns, 5Nn ]
Starting [ 2GL, 3Ns, 5Ns ]
Starting [ 2GL, 3Ns, 5Ns.2.1 ]
Starting [ 2GL, 3Ns, 5S4 ]
Time: 95.700
  finite: 8 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 30:
  [ 2GL, 3B, 5B.4.1 ]: finite, coarse curve 15.48.1.a.1
  [ 2GL, 3B, 5B.4.2 ]: finite, coarse curve 15.48.1.a.2
  [ 2GL, 3Nn, 5B.4.1 ]: finite, coarse curve 30.72.5.bn.2
  [ 2GL, 3Nn, 5B.4.1 ]: finite, coarse curve 30.72.5.bn.2
  [ 2GL, 3Nn, 5B.4.2 ]: finite, coarse curve 30.72.5.bn.1
  [ 2GL, 3Nn, 5B.4.2 ]: finite, coarse curve 30.72.5.bn.1
  [ 2GL, 3Ns, 5B.4.1 ]: finite, coarse curve 15.72.3.e.1
  [ 2GL, 3Ns, 5B.4.2 ]: finite, coarse curve 15.72.3.e.2

Searching m = 35
Time: 0.000
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 35:

Searching m = 39
Time: 0.000
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 39:

Searching m = 42
Starting [ 2Cn, 3B, 7B ]
Starting [ 2Cn, 3B, 7B.1.2 ]
Starting [ 2Cn, 3B, 7B.1.4 ]
Starting [ 2Cn, 3B, 7B.1.5 ]
Starting [ 2Cn, 3B, 7B.1.6 ]
Starting [ 2Cn, 3B, 7B.2.1 ]
Starting [ 2Cn, 3B, 7B.2.3 ]
Starting [ 2Cn, 3B, 7B.6.1 ]
Starting [ 2Cn, 3B, 7B.6.2 ]
Starting [ 2Cn, 3B, 7B.6.3 ]
Starting [ 2Cn, 3B, 7Nn ]
Starting [ 2Cn, 3B, 7Ns ]
Starting [ 2Cn, 3B, 7Ns.2.1 ]
Starting [ 2Cn, 3B, 7Ns.3.1 ]
Starting [ 2Cn, 3Cs, 7B ]
Starting [ 2Cn, 3Cs, 7B.1.2 ]
Starting [ 2Cn, 3Cs, 7B.1.4 ]
Starting [ 2Cn, 3Cs, 7B.1.5 ]
Starting [ 2Cn, 3Cs, 7B.1.6 ]
Starting [ 2Cn, 3Cs, 7B.2.1 ]
Starting [ 2Cn, 3Cs, 7B.2.3 ]
Starting [ 2Cn, 3Cs, 7B.6.1 ]
Starting [ 2Cn, 3Cs, 7B.6.2 ]
Starting [ 2Cn, 3Cs, 7B.6.3 ]
Starting [ 2Cn, 3Cs, 7Nn ]
Starting [ 2Cn, 3Cs, 7Ns ]
Starting [ 2Cn, 3Cs, 7Ns.2.1 ]
Starting [ 2Cn, 3Cs, 7Ns.3.1 ]
Starting [ 2Cn, 3GL, 7B ]
Starting [ 2Cn, 3GL, 7B.1.2 ]
Starting [ 2Cn, 3GL, 7B.1.4 ]
Starting [ 2Cn, 3GL, 7B.1.5 ]
Starting [ 2Cn, 3GL, 7B.1.6 ]
Starting [ 2Cn, 3GL, 7B.2.1 ]
Starting [ 2Cn, 3GL, 7B.2.3 ]
Starting [ 2Cn, 3GL, 7B.6.1 ]
Starting [ 2Cn, 3GL, 7B.6.2 ]
Starting [ 2Cn, 3GL, 7B.6.3 ]
Starting [ 2Cn, 3GL, 7Ns ]
Starting [ 2Cn, 3GL, 7Ns.2.1 ]
Starting [ 2Cn, 3GL, 7Ns.3.1 ]
Starting [ 2Cn, 3Nn, 7B.1.4 ]
Starting [ 2Cn, 3Nn, 7B.1.6 ]
Starting [ 2Cn, 3Nn, 7B.2.1 ]
Starting [ 2Cn, 3Nn, 7B.2.3 ]
Starting [ 2Cn, 3Nn, 7B.6.1 ]
Starting [ 2Cn, 3Nn, 7B.6.3 ]
Starting [ 2Cn, 3Nn, 7Ns.2.1 ]
Starting [ 2Cn, 3Ns, 7B ]
Starting [ 2Cn, 3Ns, 7B.1.2 ]
Starting [ 2Cn, 3Ns, 7B.1.4 ]
Starting [ 2Cn, 3Ns, 7B.1.5 ]
Starting [ 2Cn, 3Ns, 7B.1.6 ]
Starting [ 2Cn, 3Ns, 7B.2.1 ]
Starting [ 2Cn, 3Ns, 7B.2.3 ]
Starting [ 2Cn, 3Ns, 7B.6.1 ]
Starting [ 2Cn, 3Ns, 7B.6.2 ]
Starting [ 2Cn, 3Ns, 7B.6.3 ]
Starting [ 2Cn, 3Ns, 7Nn ]
Starting [ 2Cn, 3Ns, 7Ns ]
Starting [ 2Cn, 3Ns, 7Ns.2.1 ]
Starting [ 2Cn, 3Ns, 7Ns.3.1 ]
Starting [ 2GL, 3B, 7B ]
Starting [ 2GL, 3B, 7B.1.2 ]
Starting [ 2GL, 3B, 7B.1.4 ]
  Found candidate of genus 11: finite
Starting [ 2GL, 3B, 7B.1.5 ]
Starting [ 2GL, 3B, 7B.1.6 ]
  Found candidate of genus 11: finite
Starting [ 2GL, 3B, 7B.2.1 ]
Starting [ 2GL, 3B, 7B.2.3 ]
Starting [ 2GL, 3B, 7B.6.1 ]
  Found candidate of genus 11: finite
  Found candidate of genus 3: finite
  Found candidate of genus 11: finite
  Found candidate of genus 11: finite
  Found candidate of genus 5: finite
  Found candidate of genus 5: finite
Starting [ 2GL, 3B, 7B.6.2 ]
Starting [ 2GL, 3B, 7B.6.3 ]
  Found candidate of genus 11: finite
  Found candidate of genus 3: finite
  Found candidate of genus 11: finite
  Found candidate of genus 5: finite
  Found candidate of genus 5: finite
  Found candidate of genus 11: finite
Starting [ 2GL, 3B, 7Nn ]
Starting [ 2GL, 3B, 7Ns ]
Starting [ 2GL, 3B, 7Ns.2.1 ]
Starting [ 2GL, 3B, 7Ns.3.1 ]
Starting [ 2GL, 3Cs, 7B ]
Starting [ 2GL, 3Cs, 7B.1.2 ]
Starting [ 2GL, 3Cs, 7B.1.4 ]
  Found candidate of genus 37: finite
Starting [ 2GL, 3Cs, 7B.1.5 ]
Starting [ 2GL, 3Cs, 7B.1.6 ]
  Found candidate of genus 37: finite
Starting [ 2GL, 3Cs, 7B.2.1 ]
Starting [ 2GL, 3Cs, 7B.2.3 ]
Starting [ 2GL, 3Cs, 7B.6.1 ]
  Found candidate of genus 37: finite
  Found candidate of genus 13: finite
  Found candidate of genus 25: finite
  Found candidate of genus 37: finite
Starting [ 2GL, 3Cs, 7B.6.2 ]
Starting [ 2GL, 3Cs, 7B.6.3 ]
  Found candidate of genus 37: finite
  Found candidate of genus 13: finite
  Found candidate of genus 25: finite
  Found candidate of genus 37: finite
Starting [ 2GL, 3Cs, 7Nn ]
Starting [ 2GL, 3Cs, 7Ns ]
Starting [ 2GL, 3Cs, 7Ns.2.1 ]
Starting [ 2GL, 3Cs, 7Ns.3.1 ]
Starting [ 2GL, 3GL, 7B ]
Starting [ 2GL, 3GL, 7B.1.2 ]
Starting [ 2GL, 3GL, 7B.1.4 ]
Starting [ 2GL, 3GL, 7B.1.5 ]
Starting [ 2GL, 3GL, 7B.1.6 ]
Starting [ 2GL, 3GL, 7B.2.1 ]
Starting [ 2GL, 3GL, 7B.2.3 ]
Starting [ 2GL, 3GL, 7B.6.1 ]
Starting [ 2GL, 3GL, 7B.6.2 ]
Starting [ 2GL, 3GL, 7B.6.3 ]
Starting [ 2GL, 3GL, 7Nn ]
Starting [ 2GL, 3GL, 7Ns ]
Starting [ 2GL, 3GL, 7Ns.2.1 ]
Starting [ 2GL, 3GL, 7Ns.3.1 ]
Starting [ 2GL, 3Nn, 7B ]
Starting [ 2GL, 3Nn, 7B.1.2 ]
Starting [ 2GL, 3Nn, 7B.1.4 ]
  Found candidate of genus 19: finite
Starting [ 2GL, 3Nn, 7B.1.5 ]
Starting [ 2GL, 3Nn, 7B.1.6 ]
  Found candidate of genus 19: finite
Starting [ 2GL, 3Nn, 7B.2.1 ]
Starting [ 2GL, 3Nn, 7B.2.3 ]
Starting [ 2GL, 3Nn, 7B.6.1 ]
  Found candidate of genus 19: finite
  Found candidate of genus 10: finite
  Found candidate of genus 10: finite
Starting [ 2GL, 3Nn, 7B.6.2 ]
Starting [ 2GL, 3Nn, 7B.6.3 ]
  Found candidate of genus 19: finite
  Found candidate of genus 10: finite
  Found candidate of genus 10: finite
Starting [ 2GL, 3Nn, 7Nn ]
Starting [ 2GL, 3Nn, 7Ns ]
Starting [ 2GL, 3Nn, 7Ns.2.1 ]
Starting [ 2GL, 3Nn, 7Ns.3.1 ]
Starting [ 2GL, 3Ns, 7B ]
Starting [ 2GL, 3Ns, 7B.1.2 ]
Starting [ 2GL, 3Ns, 7B.1.4 ]
  Found candidate of genus 19: finite
  Found candidate of genus 37: finite
Starting [ 2GL, 3Ns, 7B.1.5 ]
Starting [ 2GL, 3Ns, 7B.1.6 ]
  Found candidate of genus 37: finite
  Found candidate of genus 19: finite
Starting [ 2GL, 3Ns, 7B.2.1 ]
Starting [ 2GL, 3Ns, 7B.2.3 ]
Starting [ 2GL, 3Ns, 7B.6.1 ]
  Found candidate of genus 7: finite
  Found candidate of genus 37: finite
  Found candidate of genus 19: finite
  Found candidate of genus 19: finite
  Found candidate of genus 13: finite
  Found candidate of genus 19: finite
Starting [ 2GL, 3Ns, 7B.6.2 ]
Starting [ 2GL, 3Ns, 7B.6.3 ]
  Found candidate of genus 7: finite
  Found candidate of genus 19: finite
  Found candidate of genus 37: finite
  Found candidate of genus 19: finite
  Found candidate of genus 19: finite
  Found candidate of genus 13: finite
Starting [ 2GL, 3Ns, 7Nn ]
Starting [ 2GL, 3Ns, 7Ns ]
Starting [ 2GL, 3Ns, 7Ns.2.1 ]
Starting [ 2GL, 3Ns, 7Ns.3.1 ]
Time: 233.100
  finite: 8 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 42:
  [ 2GL, 3B, 7B.6.1 ]: finite, coarse curve 21.96.3.a.1
  [ 2GL, 3B, 7B.6.3 ]: finite, coarse curve 21.96.3.a.2
  [ 2GL, 3Nn, 7B.6.1 ]: finite, coarse curve 42.144.10.t.2
  [ 2GL, 3Nn, 7B.6.1 ]: finite, coarse curve 42.144.10.t.2
  [ 2GL, 3Nn, 7B.6.3 ]: finite, coarse curve 42.144.10.t.1
  [ 2GL, 3Nn, 7B.6.3 ]: finite, coarse curve 42.144.10.t.1
  [ 2GL, 3Ns, 7B.6.1 ]: finite, coarse curve 21.144.7.b.1
  [ 2GL, 3Ns, 7B.6.3 ]: finite, coarse curve 21.144.7.b.2

Searching m = 70
Starting [ 2GL, 5Nn, 7Ns ]
Starting [ 2GL, 5Nn, 7Ns.2.1 ]
Starting [ 2GL, 5Nn, 7Ns.3.1 ]
Starting [ 2GL, 5Ns, 7Nn ]
Starting [ 2GL, 5Ns, 7Ns ]
Starting [ 2GL, 5Ns, 7Ns.2.1 ]
Starting [ 2GL, 5Ns, 7Ns.3.1 ]
Starting [ 2GL, 5Ns.2.1, 7Ns ]
Starting [ 2GL, 5Ns.2.1, 7Ns.2.1 ]
Starting [ 2GL, 5Ns.2.1, 7Ns.3.1 ]
Starting [ 2GL, 5S4, 7Ns ]
Starting [ 2GL, 5S4, 7Ns.2.1 ]
Starting [ 2GL, 5S4, 7Ns.3.1 ]
Time: 53.540
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 70:

Searching m = 78
Starting [ 2GL, 3GL, 13B ]
Starting [ 2GL, 3GL, 13B.3.1 ]
Starting [ 2GL, 3GL, 13B.3.2 ]
Starting [ 2GL, 3GL, 13B.3.4 ]
Starting [ 2GL, 3GL, 13B.3.7 ]
Starting [ 2GL, 3GL, 13B.4.1 ]
Starting [ 2GL, 3GL, 13B.4.2 ]
Starting [ 2GL, 3GL, 13B.5.1 ]
Starting [ 2GL, 3GL, 13B.5.2 ]
Starting [ 2GL, 3GL, 13B.5.4 ]
Time: 10.180
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 78:

Searching m = 105
Starting [ 3GL, 5Nn, 7Ns.2.1 ]
Starting [ 3GL, 5Ns, 7Ns ]
Starting [ 3GL, 5Ns, 7Ns.2.1 ]
Starting [ 3GL, 5Ns, 7Ns.3.1 ]
Starting [ 3GL, 5Ns.2.1, 7Ns.2.1 ]
Starting [ 3GL, 5S4, 7Ns.2.1 ]
Starting [ 3Nn, 5Ns, 7Ns.2.1 ]
Starting [ 3Ns, 5Nn, 7Ns ]
Starting [ 3Ns, 5Nn, 7Ns.2.1 ]
Starting [ 3Ns, 5Nn, 7Ns.3.1 ]
Starting [ 3Ns, 5Ns, 7Nn ]
Starting [ 3Ns, 5Ns, 7Ns ]
Starting [ 3Ns, 5Ns, 7Ns.2.1 ]
Starting [ 3Ns, 5Ns, 7Ns.3.1 ]
Starting [ 3Ns, 5Ns.2.1, 7Ns ]
Starting [ 3Ns, 5Ns.2.1, 7Ns.2.1 ]
Starting [ 3Ns, 5Ns.2.1, 7Ns.3.1 ]
Starting [ 3Ns, 5S4, 7Ns ]
Starting [ 3Ns, 5S4, 7Ns.2.1 ]
Starting [ 3Ns, 5S4, 7Ns.3.1 ]
Time: 431.210
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 105:

Searching m = 210
Starting [ 2GL, 3GL, 5Nn, 7Nn ]
Starting [ 2GL, 3GL, 5S4, 7Nn ]
Starting [ 2GL, 3GL, 5S4, 7Ns ]
Starting [ 2GL, 3Nn, 5Nn, 7Nn ]
Time: 6742.860
  finite: 0 recorded groups
  infinite: 0 recorded groups
  unresolved: 0 recorded groups

Modular curves for m = 210:
*/