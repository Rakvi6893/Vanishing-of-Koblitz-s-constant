D:=CremonaDatabase();
H:=sub<GL(2,Integers(15))|[[11, 10, 0, 8], [11, 12, 0, 14], [14, 5, 0, 13]]>;//15.96.1-15.a.1.5
E:=EllipticCurve(D, "50a3");//50.a2
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(15))|[GL(2,Integers(15))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G15:=pi(G);
 Index(GL(2,Integers(15)),G15); //192
 -Identity(GL(2,Integers(15))) in G15; // false
 #LowIndexSubgroups(G15,2);//4
 for d in [-3,5,-15] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G15d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=120) do;
         IsConjugate(GL(2,Integers(15)),T`subgroup,G15d);
     end for;
end for;

//-15 is the only twist whose mod-15 image is contained in H

print "************************************************";

E:=EllipticCurve(D, "50a4");//50.a4
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(15))|[GL(2,Integers(15))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G15:=pi(G);
 Index(GL(2,Integers(15)),G15); //192
 -Identity(GL(2,Integers(15))) in G15; // false
 #LowIndexSubgroups(G15,2);//4
 for d in [-3,5,-15] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G15d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=120) do;
         IsConjugate(GL(2,Integers(15)),T`subgroup,G15d);
     end for;
end for;

//-15 is the only twist whose mod-15 image is contained in H



H:=sub<GL(2,Integers(15))|[[2, 8, 0, 11], [4, 14, 0, 14], [13, 1, 0, 1]]>;//15.96.1-15.a.2.5
E:=EllipticCurve(D, "50a1");//50.a3
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(15))|[GL(2,Integers(15))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G15:=pi(G);
 Index(GL(2,Integers(15)),G15); //192
 -Identity(GL(2,Integers(15))) in G15; // false
 #LowIndexSubgroups(G15,2);//4
 for d in [-3,5,-15] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G15d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=120) do;
         IsConjugate(GL(2,Integers(15)),T`subgroup,G15d);
     end for;
end for;

//-15 is the only twist whose mod-15 image is contained in H

print "************************************************";

E:=EllipticCurve(D, "50a2");//50.a1
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(15))|[GL(2,Integers(15))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G15:=pi(G);
 Index(GL(2,Integers(15)),G15); //192
 -Identity(GL(2,Integers(15))) in G15; // false
 #LowIndexSubgroups(G15,2);//4
 for d in [-3,5,-15] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G15d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=120) do;
         IsConjugate(GL(2,Integers(15)),T`subgroup,G15d);
     end for;
end for;

//-15 is the only twist whose mod-15 image is contained in H


H:=sub<GL(2,Integers(30))|[[19, 0, 12, 17], [19, 15, 27, 26], [29, 20, 18, 7], [29, 20, 27, 29]]>;//30.96.1-15.a.2.3
E:=EllipticCurve(D, "50a1");//50.a3
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(30))|[GL(2,Integers(30))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G30:=pi(G);
 Index(GL(2,Integers(30)),G30); //192
 -Identity(GL(2,Integers(30))) in G30; // false
 #LowIndexSubgroups(G30,2);//8
 for d in [-2,-3,5,-15,6,-10,30] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G30d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=720) do;
         IsConjugate(GL(2,Integers(30)),T`subgroup,G30d);
     end for;
end for;

//30 is the only twist whose mod-15 image is contained in H

print "************************************************";

E:=EllipticCurve(D, "50a2");//50.a1
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(30))|[GL(2,Integers(30))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G30:=pi(G);
 Index(GL(2,Integers(30)),G30); //192
 -Identity(GL(2,Integers(30))) in G30; // false
 #LowIndexSubgroups(G30,2);//8
 for d in [-2,-3,5,-15,6,-10,30] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G30d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=720) do;
         IsConjugate(GL(2,Integers(30)),T`subgroup,G30d);
     end for;
end for;

//30 is the only twist whose mod-15 image is contained in H

print "************************************************";

D:=CremonaDatabase();
H:=sub<GL(2,Integers(30))|[[1, 5, 12, 11], [4, 5, 15, 29], [8, 15, 9, 16], [26, 25, 27, 19]]>;//30.96.1-15.a.1.7
E:=EllipticCurve(D, "50a3");//50.a2
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(30))|[GL(2,Integers(30))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G30:=pi(G);
 Index(GL(2,Integers(30)),G30); //192
 -Identity(GL(2,Integers(30))) in G30; // false
 #LowIndexSubgroups(G30,2);//8
 for d in [-2,-3,5,-15,6,-10,30] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G30d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=720) do;
         IsConjugate(GL(2,Integers(30)),T`subgroup,G30d);
     end for;
end for;

//30 is the only twist whose mod-15 image is contained in H

print "************************************************";

E:=EllipticCurve(D, "50a4");//50.a4
G:=FindOpenImage(E);
pi:=hom<GL(2,Integers(120))->GL(2,Integers(30))|[GL(2,Integers(30))!GL(2,Integers(120)).i:i in [1..#Generators(GL(2,Integers(120)))]]>;
 G30:=pi(G);
 Index(GL(2,Integers(30)),G30); //192
 -Identity(GL(2,Integers(30))) in G30; // false
 #LowIndexSubgroups(G30,2);//8
 for d in [-2,-3,5,-15,6,-10,30] do;
     d;
     Ed:=QuadraticTwist(E,d);
     CremonaReference(Ed);
     G30d:=pi(FindOpenImage(Ed));
     for T in Subgroups(H:OrderEqual:=720) do;
         IsConjugate(GL(2,Integers(30)),T`subgroup,G30d);
     end for;
end for;

//30 is the only twist whose mod-15 image is contained in H

print "************************************************";
