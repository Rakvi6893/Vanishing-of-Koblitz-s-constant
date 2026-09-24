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

pi15:=hom<GL(2,Integers(120))->GL(2,Integers(15))|[GL(2,Integers(15))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
pi30:=hom<GL(2,Integers(120))->GL(2,Integers(30))|[GL(2,Integers(30))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;



// testing coarse curves

for label in ["50a1","50a2","50a3","50a4"] do;
    E:=EllipticCurve(D, label);
    G:=FindOpenImage(E);
    G15:=pi15(G); assert #LowIndexSubgroups(G15,2) eq 4;
    G30:=pi30(G); assert #LowIndexSubgroups(G30,2) eq 8;
    S:={g:g in Generators(G15)};
    S:=S join {-Identity(GL(2,Integers(15)))};
    G15coarse:=sub<GL(2,Integers(15))|[g:g in S]>;
    if HasPrimitiveCongObstruction(G15coarse) eq true then;label;end if;
    S:={g:g in Generators(G30)}; 
    S:=S join {-Identity(GL(2,Integers(30)))};
    G30coarse:=sub<GL(2,Integers(30))|[g:g in S]>;
    if HasPrimitiveCongObstruction(G30coarse)eq true then ;label;end if;
    
end for;

for label in ["50a1","50a2","50a3","50a4"] do;
 label;
 E:=EllipticCurve(D, label);


 for d in [1,-2,-3,5,-15,6,-10,30] do;
 
  Ed:=QuadraticTwist(E,d);
  Gd:=FindOpenImage(Ed);
  Gd15:=pi15(Gd);
  Gd30:=pi30(Gd);
  if HasPrimitiveCongObstruction(Gd15) eq true then 
     CremonaReference(Ed);
  end if;
 
  if HasPrimitiveCongObstruction(Gd30) eq true then
     CremonaReference(Ed);
  end if;
 
 
  end for;
  print ".........";
 end for;

 /* 50a1
450d3
14400dy3
.........
50a2
450d4
14400dy4
.........
50a3
450d1
14400dy1
.........
50a4
450d2
14400dy2
.........*/
 

