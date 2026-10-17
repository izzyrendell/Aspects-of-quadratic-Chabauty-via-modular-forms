N := 29;
Gr := GL(2, IntegerRing(N));
c := Gr ! [14,11, 26,15];
b := Gr ! [15,12, 17,3];
G := sub<Gr| b,c>;