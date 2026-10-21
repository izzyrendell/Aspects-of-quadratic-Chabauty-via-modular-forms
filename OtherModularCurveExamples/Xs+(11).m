// this code requires Zywina's Modular code: https://github.com/davidzywina/Modular

// setup: level N, subgroup G in GL(2,Z/N)
N := 11;
Gr := GL(2, IntegerRing(N));
c := Gr ! [0,8, 2,0];
b := Gr ! [9,0, 0,7]; // from beta.lmfdb.org
G := sub<Gr| b,c>; 
// Order(G); // sanity check; has the right order 2*(N-1)^2

// set up to generate spaces of cusp forms
X := CreateModularCurveRec(G);
g := X`genus;
N2 := 61; // =(N^2+1)/2

// weight 2 forms (the holomorphic differentials)
M2 := FindModularForms(2,X);
S2 := FindCuspForms(M2);
S2 := IncreaseModularFormPrecision(S2, [100,100,100,100,100,100]);
cusp_forms := [ S2`F0[i] : i in [1..#S2`F0]]; //expansions of all basis elements (at all cusps) up to O(qw^100)
cusp_forms_at_infinity := [cusp_forms[i][1]: i in [1..#cusp_forms]]; // expansions of both basis elements at infinity

// calculation of denominator (delta) 
M12 := FindModularForms(12,X);
j := 8;
S12superspecial := FindModularFormsWithVanishingConditions(M12,[N2,1,1,1,1,1]);
S12superspecial := IncreaseModularFormPrecision(S12superspecial, [500,12,12,12,12,12]);
delta := -S12superspecial`F0[1][1]/N^6; // expansion of denominator at infinity

// checking our denominator is in fact an eta-quotient
M := 2000;
R<q> := PuiseuxSeriesRing(Rationals(), M);
eta_quot := function(lambda, prec);
A := &*[(1-q^(lambda*n)): n in [1..prec]];
return q^(lambda/24)*A;
end function;
eta_quotient := eta_quot(1, 1000)^12*eta_quot(N^2, 1000)^12;
// print delta/eta_quotient;

// calculation of numerator 
S12special := FindModularFormsWithVanishingConditions(M12,[53,1,1,1,1,1]); // N2-j=53     
S12special := IncreaseModularFormPrecision(S12special,[100,12,12,12,12,12]);
numerator_candidates := [S12special`F0[i][1] : i in [1..#S12special`F0]]; // expansions of these forms at infinity

// construction of basis candidate
M := 100;
R<q> := LaurentSeriesRing(RationalField(),M);
varphi_dR := numerator_candidates[1]/delta;
basis_candidate := cusp_forms_at_infinity cat [varphi_dR*f: f in cusp_forms_at_infinity];

// verifying we have a basis - only need to consider expansions at infinity
CupModForms:=function(f,g); // def [w_g] u [w_f]
A:=Integral(f/q); // send f to integral of omega f
B:=A*g;
return Coefficient(B,0); // taking residue
end function;
CMx:=Transpose(Matrix(#basis_candidate,#basis_candidate,[[CupModForms(basis_candidate[i],basis_candidate[j]):j in [1..#basis_candidate]]:i in [1..#basis_candidate]]));
//print Determinant(CMx); // this is nonzero so we have found a basis!

