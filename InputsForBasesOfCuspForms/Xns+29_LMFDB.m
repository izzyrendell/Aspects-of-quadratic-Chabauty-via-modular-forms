N := 29;
Gr := GL(2, IntegerRing(N));
// LMFDB generators of the normaliser of a nonsplit Cartan subgroup of GL(2,Z/29)
c := Gr ! [14, 11, 26, 15];
b := Gr ! [15, 12, 17, 3];
G := sub<Gr| b,c>;
