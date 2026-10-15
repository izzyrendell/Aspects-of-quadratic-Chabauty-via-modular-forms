N := 19;
Gr := GL(2, IntegerRing(N));
// LMFDB generators of the normaliser of a nonsplit Cartan subgroup of GL(2,Z/19)
c := Gr ! [7,10, 13,12];
b := Gr ! [10,10, 18,0];
G := sub<Gr| b,c>;
