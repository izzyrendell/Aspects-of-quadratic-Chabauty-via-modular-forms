N := 13;
Gr := GL(2, IntegerRing(N));
c := Gr ! [0,12, 11,0];
b := Gr ! [1,8, 2,12];
G := sub<Gr| b,c>;