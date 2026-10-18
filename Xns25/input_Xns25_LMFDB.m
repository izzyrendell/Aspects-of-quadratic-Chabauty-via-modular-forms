N := 25;
Gr := GL(2, IntegerRing(N));
c := Gr ! [5,21, 16,20];
b := Gr ! [11,9, 23,14];
G := sub<Gr| b,c>;