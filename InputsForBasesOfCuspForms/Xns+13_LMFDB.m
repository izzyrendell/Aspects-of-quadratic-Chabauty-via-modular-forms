N := 13;
Gr := GL(2, IntegerRing(N));
// LMFDB generators of the normaliser of a nonsplit Cartan subgroup of GL(2,Z/13)
c := Gr ! [0,12, 11,0];
b := Gr ! [1,8, 2,12];
G := sub<Gr| b,c>; // This is the normalizer of a nonsplit Cartan
