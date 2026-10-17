N := 23;
Gr := GL(2, IntegerRing(N));
c := Gr ! [2,10, 8,21];
b := Gr ! [2,19, 4,6];
G := sub<Gr| b,c>; // These are the generators from LMFDB,