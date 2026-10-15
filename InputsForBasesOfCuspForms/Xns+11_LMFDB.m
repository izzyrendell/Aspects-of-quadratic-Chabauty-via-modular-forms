N := 11;
Gr := GL(2, IntegerRing(N));
// LMFDB generators of the normaliser of a nonsplit Cartan subgroup of GL(2,Z/11)
c := Gr ! [2,7, 5,9];
b := Gr ! [6,1, 6,5];
G := sub<Gr| b,c>; // This is the normalizer of a nonsplit Cartan
