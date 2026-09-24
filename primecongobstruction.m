load "Common.m";

// Compute F1
for lbl in Keys(Data) do
    ell := Data[lbl][1];
    if ell notin [2,3,5,7] then continue; end if;
    G := sub<GL(2, Integers(ell)) | Data[lbl][3]>;
    assert Index(GL(2, Integers(ell)), G) eq Data[lbl][2];
    if F1(G) eq 1 then lbl; end if;
end for;

/*Output 5B.1.1
5B.1.2
5Cs.1.1
2Cs
7B.1.1
7B.1.3
2B
3B.1.1
3B.1.2
3Cs.1.1 */
