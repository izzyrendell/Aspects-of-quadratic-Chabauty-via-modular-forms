// This file has to be placed in the folder Modular
AttachSpec("Modular.spec");

print "Sanity check that our subgroup of GL(2,F_N) has the correct order for nsplit+ group:", Order(G) eq 2*(N^2 - 1); // will be wrong for non-prime N
// Zywina's Code: Create modular curve X
X := CreateModularCurveRec(G);

// weight 2 forms
T := Time();
M2 := FindModularForms(2,X);
printf "Finding weight 2 forms took %o seconds\n", Time(T);

// weight 12 forms
T := Time();
M12 := FindModularForms(12,X);
printf "Finding weight 12 forms took %o seconds\n", Time(T);

// if you want, save the forms for later
// save "formscomputed";