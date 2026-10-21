q_prec := 200;
R<q>:=LaurentSeriesRing(Rationals(),q_prec);
level := 37;
k := 12;

// computing weight 2 forms to give 1st half of basis
C2:=CuspForms(level,2); 
B2:=Basis(C2,q_prec); 

// computing $\varphi_{dR}$ from weight 12 forms
Ck:=CuspForms(level,k); 
Bk:=Basis(Ck,q_prec);
varphi_dR := Bk[#Bk-1]/Bk[#Bk]; // naming quotient 

// checking our denominator is a scalar factor of w_N(\Delta)
delta := Delta(q);
w_N_delta := level^6*Evaluate(delta, q^level);
// print Bk[#Bk]/w_N_delta;

// constructing our basis candidate
B2edit:= [varphi_dR*B2[i]: i in [1..#B2]];
B2cat:= B2 cat B2edit; // test if this is a basis for H1dR

// cup product
CupModForms:=function(f,g);
A:=Integral(f/q);
B:=A*g;
return Coefficient(B,0);
end function;

// cup product matrix
CMx:=Matrix(Rationals(),#B2cat,#B2cat,[[CupModForms(B2cat[i],B2cat[j]):j in [1..#B2cat]]:i in [1..#B2cat]]);
// print Determinant(CMx); // this is equal to 1/4 so we have a basis!
