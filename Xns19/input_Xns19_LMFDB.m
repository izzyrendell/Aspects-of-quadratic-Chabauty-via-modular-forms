N := 19;
Gr := GL(2, IntegerRing(N));
c := Gr ! [7,10, 13,12];
b := Gr ! [10,10, 18,0];
G := sub<Gr| b,c>;