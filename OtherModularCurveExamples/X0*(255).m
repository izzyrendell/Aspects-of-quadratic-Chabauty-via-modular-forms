T := Time();
level := 255; 
prime1 := 3;
prime2 := 5;
prime3 := 17;
genus := 2;
M := 250;
R<q>:=LaurentSeriesRing(Rationals(),M);

// weight 2 cusp forms invariant under wp = prime1, wq = prime2, wr = prime3
C2:=CuspForms(level,2);
B2:=Basis(C2,M); 
A2p:=AtkinLehnerOperator(C2,prime1); // matrix of wp acting on all cusp forms
A2q:=AtkinLehnerOperator(C2,prime2); // matrix of wq acting on all cusp forms
A2r:=AtkinLehnerOperator(C2,prime3); // matrix of wr acting on all cusp forms
K2 := Basis(Kernel(A2p-1) meet Kernel(A2q-1) meet Kernel(A2r-1)); // this has dimension 2 as expected 
M2:=[&+[R!K2[i][j]*B2[j]:j in [1..#B2]]:i in [1..#K2]]; // invariant under W(N)

// same process as above with k=12
k := 12;
Ck:=CuspForms(level,k);
Bk:=Basis(Ck,M); 
Akp:=AtkinLehnerOperator(Ck,prime1); // matrix of wp acting on all cusp forms
Akq:=AtkinLehnerOperator(Ck,prime2); // matrix of wq acting on all cusp forms
Akr:=AtkinLehnerOperator(Ck,prime3); // matrix of wr acting on all cusp forms
Kk := Basis(Kernel(Akp-1) meet Kernel(Akq-1) meet Kernel(Akr-1));
Mk:=[&+[R!Kk[i][j]*Bk[j]:j in [1..#Bk]]:i in [1..#Kk]]; // invariant under W(N)
denom_candidate := Mk[#Mk];

// checking whether this is an eta quotient 
eta_prec := 300;
S<q_N>:= PuiseuxSeriesRing(Rationals(), M); 
eta_quot := function(delta, prec);
A := &*[(1-q_N^(delta*n)): n in [1..prec]];
return q_N^(delta/24)*A;
end function;
eta_quotient := (eta_quot(1, eta_prec)*eta_quot(prime1, eta_prec)*eta_quot(prime2, eta_prec)*eta_quot(prime3, eta_prec)*eta_quot(prime1*prime2, eta_prec)*eta_quot(prime1*prime3, eta_prec)*eta_quot(prime2*prime3, eta_prec)*eta_quot(level, eta_prec))^3; 
// print denom_candidate/eta_quotient;

// computing varphi_dR
for j in [1..#Mk-1] do
  if Valuation((Mk[#Mk-j]/Mk[#Mk])*M2[#M2]) lt 0 then
    j_dR := j;
    varphi_dR := Mk[#Mk-j_dR]/Mk[#Mk];
    break j;
  end if;
end for;

// constructing H1dR basis
M2edit:=[varphi_dR*M2[i]: i in [1..#M2]]; // multiplying wt 2 cusp forms by quotient
M2cat:=M2 cat M2edit; // will test if these are a basis for H1dR

q:=Parent(B2[1]).1;

CupModForms:=function(f,g); // def [w_g] u [w_f]
A:=Integral(f/q); // send f to integral of omega f
B:=A*g;
return Coefficient(B,0); // taking residue
end function;

CMx:=Transpose(Matrix(Rationals(),#M2cat,#M2cat,[[CupModForms(M2cat[i],M2cat[j]):j in [1..#M2cat]]:i in [1..#M2cat]]));

// det below is nonzero so we have found a basis
print Determinant(CMx);
print Time(T);


