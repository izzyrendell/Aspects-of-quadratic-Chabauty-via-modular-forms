N := 23;
Gr := GL(2, IntegerRing(N));
// LMFDB generators of the normaliser of a nonsplit Cartan subgroup of GL(2,Z/23)
c := Gr ! [2, 10, 8, 21];
b := Gr ! [2, 19, 4, 6];
G := sub<Gr| b,c>;
