N := 17;
Gr := GL(2, IntegerRing(N));
// LMFDB generators of the normaliser of a nonsplit Cartan subgroup of GL(2,Z/17)
c := Gr ! [7,3, 13,10];
b := Gr ! [7,16, 1,8];
G := sub<Gr| b,c>; // This is the normalizer of a nonsplit Cartan
