// we load our input curve: level N, subgroup G of Gr := GL_2(F_N)
// load "input_Xns11_LMFDB.m";

print "Sanity check that our subgroup of GL(2,F_N) has the correct order for nsplit+ group:", Order(G) eq 2*(N^2 - 1);
// Zywina's Code: Create modular curve X
X := CreateModularCurveRec(G);

// fix a prec for the q-expansions
prec := 25;
nu := (N-1)/2;

// weight 2 forms
M2 := FindModularForms(2,X);

// increase precision
M2 := IncreaseModularFormPrecision(M2,[prec : i in [1..nu]]);

// cusp forms
S2 := FindCuspForms(M2);
cusp_forms := [ S2`F0[i][1] : i in [1..#S2`F0]]; // only the q-exp at the first cusp
cusp_forms_all_cusps := S2`F0;

// weight 12 forms
M12 := FindModularForms(12,X);

// increase precision
M12 := IncreaseModularFormPrecision(M12,[prec : i in [1..nu]]);

// cusp forms
S12 := FindCuspForms(M12);

// save magma session thus far
// save "nsplit11formscomputed";
