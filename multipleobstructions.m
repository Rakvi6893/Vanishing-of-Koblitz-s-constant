gps:=[<"6.24.0-3.a.1.1",6,[[1, 0, 3, 5], [2, 3, 3, 2]]>,<"6.8.0-3.a.1.1",6,[[2, 5, 3, 4], [5, 1, 0, 5]]>,<"6.8.0-3.a.1.2",6,[[4, 3, 3, 2], [5, 5, 3, 4]]>,<"10.120.0-5.a.1.1",10,[[1, 3, 5, 2], [9, 9, 0, 7]]>,<"10.24.0-5.a.2.2",10,[[3, 8, 3, 5], [4, 7, 5, 8]]>,<"10.24.0-5.a.1.1",10,[[1, 9, 0, 9], [8, 3, 7, 8]]>,<"14.48.0-7.a.1.1",14,[[11, 9, 8, 5], [12, 1, 1, 12]]>,<"14.48.0-7.a.2.1",14,[[6, 3, 7, 2], [13, 11, 5, 12]]>];

D:=CremonaDatabase();
for label in ["450d4","450d3","450d1","450d2","14400dy4","14400dy3","14400dy1","14400dy2"] do;
label;
E:=EllipticCurve(D, label);
G:=FindOpenImage(E);
level:=#BaseRing(G);
for group in gps do;
    
    gp:=sub<GL(2,Integers(group[2]))|group[3]>;
    pi1:=hom<GL(2,Integers(level*7))->GL(2,Integers(level))|[GL(2,Integers(level))!GL(2,Integers(level*7)).i:i in [1..#Generators(GL(2,Integers(level*7)))]]>;
pi:=hom<GL(2,Integers(level*7))->GL(2,Integers(group[2]))|[GL(2,Integers(group[2]))!GL(2,Integers(level*7)).i:i in [1..#Generators(GL(2,Integers(level*7)))]]>;
 Gd:=pi(G@@pi1);
 if IsConjugate(GL(2,Integers(group[2])),Gd,gp) eq true then;group[1];end if;
 
 end for;
 print ".........";
 end for;
 
// curves that have mod 15 or mod 30 obstruction do not admit obstruction modulo any other level


