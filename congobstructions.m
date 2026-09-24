ChangeDirectory("/Users/Jacob/OpenImage-master"); // Change this line as needed.
load "main/FindOpenImage.m";
load "Common.m";

function HasPrimitiveCongObstruction(G)
    m:=#BaseRing(G);
    S:=[];
    for d in Divisors(Squarefree(m)) do;
        if d eq 1 or d eq m then continue; end if;
        pi:=hom<GL(2,Integers(m))->GL(2,Integers(d))|[GL(2,Integers(d))!GL(2,Integers(m)).i:i in [1..#Generators(GL(2,Integers(m)))]]>;
        Gd:=pi(G); 
        S:= S cat [F1(Gd)];
    end for;
    if F1(G) eq 1 and 1 notin S then return true;
    else return false; end if;
end function;


D:=CremonaDatabase();
for label in ["450d4","450d3","450d1","450d2","14400dy4","14400dy3","14400dy1","14400dy2"] do;
label;
E:=EllipticCurve(D, label);
G:=FindOpenImage(E);
level:=#BaseRing(G);
for d in Divisors(Squarefree(level)) do;
    if d eq 1 then continue; end if;
    d;
pi:=hom<GL(2,Integers(level))->GL(2,Integers(d))|[GL(2,Integers(d))!GL(2,Integers(level)).i:i in [1..#Generators(GL(2,Integers(level)))]]>;
 Gd:=pi(G);
 HasPrimitiveCongObstruction(Gd);
 
 end for;
 print ".........";
 end for;
 
/*output 450d4
2
false
3
false
5
false
6
false
10
false
15
true
30
false
.........
450d3
2
false
3
false
5
false
6
false
10
false
15
true
30
false
.........
450d1
2
false
3
false
5
false
6
false
10
false
15
true
30
false
.........
450d2
2
false
3
false
5
false
6
false
10
false
15
true
30
false
.........
14400dy4
2
false
3
false
5
false
6
false
10
false
15
false
30
true
.........
14400dy3
2
false
3
false
5
false
6
false
10
false
15
false
30
true
.........
14400dy1
2
false
3
false
5
false
6
false
10
false
15
false
30
true
.........
14400dy2
2
false
3
false
5
false
6
false
10
false
15
false
30
true
......... */
