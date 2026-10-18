// fix a prec for the q-expansions
prec := 42;
nu := #X`cusps;
g := X`genus;

// increase precision
T1 := Time();
M2 := IncreaseModularFormPrecision(M2,[prec : i in [1..nu]]);
printf "Increasing the q-expansion precision to %o took %o seconds\n", prec, Time(T1);

// cusp forms
T2 := Time();
S2 := FindCuspForms(M2);
printf "Finding weight 2 cusp forms took %o seconds\n", Time(T2);
cusp_forms := [ S2`F0[i][1] : i in [1..#S2`F0]]; // only the q-exp at the first cusp
cusp_forms_all_cusps := S2`F0;

// increase precision
T3 := Time();
M12 := IncreaseModularFormPrecision(M12,[prec : i in [1..nu]]);
printf "Increasing the q-expansion precision to %o took %o seconds\n", prec, Time(T3);

// cusp forms
T4 := Time();
S12 := FindCuspForms(M12);
printf "Finding weight 12 cusp forms took %o seconds\n", Time(T4);

printf "Total time: %o seconds\n", Time(T1);