// to run after compute_Hecke_all_cusps.m
// so we already have the Hecke matrix T = T_p
// w.r.t. a symplectic basis

Tpmat := T;
genus := g;
CMxsymp := Transpose(BC_M2cat_symp) * C * BC_M2cat_symp; // this is the standard symplectic pairing matrix
Zmat:=(Trace(Tpmat)*ScalarMatrix(2*genus, 1) - 2*genus*Tpmat)*CMxsymp^(-1); // the matrix for the correspondence
ZmatR := ChangeRing(Zmat, R);
printf "The matrix of the correspondence associated with Hecke T_p with p = %o is:\n %o \n", p, Zmat;

// Find basis of log differentials with poles only at cusps
// candidates are h*f/Delta dq/q with h wt 12 valuation N-1 at cusps and f wt 2
S12special_num2 := FindModularFormsWithVanishingConditions(M12,[N-1 : i in [1..nu]]);
logdiff_num := S12special_num2`F0; // candidates for h (except first one, which is Delta)
printf "We have %o candidates for numerators of log differentials.\n", #logdiff_num - 1;
log_diffs := &cat [ [ [logdiff_num[i][c]*f[c]/(qw*denom_all_cusps[c]) : c in [1..nu]] : i in [2..#logdiff_num]] : f in cusp_forms_all_cusps]; // these are our candidates of log differentials
printf "We have %o candidates for log differentials. Their q-exp at 1st cusp are:\n%o \n", #log_diffs, [x[1] : x in log_diffs];
// need to narrow it down to nu-1 of them
// that form a basis together with the holom diff
// as candidates just take the first nu-1 of log_diffs
test_diff := [ [ R !  (symp_basis_all_cusps[i][c]/qw)  : c in [1..nu] ] : i in [1..g]] cat [log_diffs[i] : i in [1..(nu-1)]];
bound := Ceiling((g+nu-1)/nu);
ncoeff := [Coefficient(x[1],-1) : x in test_diff];
ncoeffQ := CoefficientMatrixFromSequence(N,ncoeff);
for n in [0..bound-2] do
    ncoeff_new := [Coefficient(x[1],n) : x in test_diff];
    ncoeffQ_new := CoefficientMatrixFromSequence(N,ncoeff_new);
    ncoeffQ := VerticalJoin(ncoeffQ, ncoeffQ_new);
end for;
printf "The rank of the candidate basis of log diff is: %o. It should be: %o. They are equal: %o\n", Rank(ncoeffQ), g+nu-1, Rank(ncoeffQ) eq g+nu-1;
// so test_diff is a basis of log diff, and the remaining omega are
omega_log_vec := [test_diff[i] : i in [g+1..g+nu-1]];

function ColumnMatrix(seq)
    K := Parent(seq[1]);
    n := #seq;
    return Matrix(K, n, 1, seq);
end function;

// Now Algorithm 5.20 from AWS notes
// (1) we use the local coordinate qw at the cusps
// (2) omega_vec
omega_vec := [ [ R !  (x[c]/qw)  : c in [1..nu] ] : x in symp_basis_all_cusps]; // sympl basis differentials (only dqw missing)
// (3) Omega_vec
Omega_vec := [[Integral(-omega_vec[x][c]) : c in [1..nu]] : x in [1..2*g]]; // primitives

// (4) eta (typo in [AWS (66)], I think)
OmTZdOm := [(Transpose(ColumnMatrix([R ! (x[c]) : x in Omega_vec])) * ZmatR * ColumnMatrix([R ! (-x[c]) : x in omega_vec]))[1][1] : c in [1..nu]];
Res_cond := [Coefficient(OmTZdOm[c], -1) : c in [1..nu]]; // RHS of linear system we need to solve
RightHandSide := CoefficientMatrixFromSequence(N,[Res_cond[1]]); // RHS Q-coeff-wise
for a in [2..#Res_cond] do
    RHS_new := CoefficientMatrixFromSequence(N,[Res_cond[a]]);
    RightHandSide := VerticalJoin(RightHandSide, RHS_new);
end for;

res_omega_log_vec := [ [ Coefficient(omega_log_vec[i][c], -1) : c in [1..nu] ] : i in [1..nu-1]];
res_omega_log_vec_matrix := CoefficientMatrixFromSequence(N, [x[1] : x in res_omega_log_vec]);
for c in [2..nu] do
    res_at_cusp_c := CoefficientMatrixFromSequence(N,[x[c] : x in res_omega_log_vec]);
    res_omega_log_vec_matrix := VerticalJoin(res_omega_log_vec_matrix, res_at_cusp_c);
end for;
alpha, _ := Solution(Transpose(res_omega_log_vec_matrix),Transpose(RightHandSide));
eta := [&+[alpha[1][a] * omega_log_vec[a][c] : a in [1..nu-1]] : c in [1..nu]];
printf "The q-exp at 1st cusp of the differential eta is: %o\n", eta[1];
printf "The residues of eta at the cusps are as desired: %o\n", [Res_cond[c] eq Coefficient(eta[c],-1) : c in [1..nu]];

// (5) g_x
g_fct := [Integral( OmTZdOm[c] - eta[c]) : c in [1..nu]]; // primitives

// (6) gamma_Fil and b_Fil
RRoffset := 0;
done := false;

while not done do // try larger and larger pole orders // I thought RR says this is ok for RRoffset = 0 ...
  porder := Maximum(Floor((2*g-2)/nu) + 1, -Valuation(symp_basis_all_cusps[g+1][1])) + RRoffset; // function gamma_Fil lives in RR space L(porder*cusps) 
  S12special_num3 := FindModularFormsWithVanishingConditions(M12,[N-porder : i in [1..nu]]);
  gammaFil_num := S12special_num3`F0; // candidates for numerator
  printf "gammaFil lives in RR space with pole order at worst %o at cusps. It has dimension %o.\n", porder, #gammaFil_num;
  gammaFil_basis := [ [gammaFil_num[i][c]/(denom_all_cusps[c]) : c in [1..nu]] : i in [1..#gammaFil_num]]; // gammaFil lives in the space spanned by these

// Need to find gammaFil, a linear comb of gammaFil_basis, and constants bFil s.t. gammaFil - sum_j bFil[j] * (integral of omega[j]) has a prescribed principal part

// Similar and partially copied from QCMod/hodge.m
// for RHS, need the OmTZNNTOm term + principal parts of all gx 
  v:=[]; // coefficients of principal parts of all gx 
    cnt:=0;
    for c:=1 to nu do
      for j:=-porder to -1 do
        cnt +:= 1;
        v[cnt]:= Coefficient(g_fct[c],j);
      end for;
    end for;

  Nmatrix := ZeroMatrix(R, 2*g, g);
  for i in [1..g] do
    Nmatrix[g+i][i] := 1;
  end for;
  OmTZNNTOm := [(Transpose(ColumnMatrix([R ! (x[c]) : x in Omega_vec])) * ZmatR * Nmatrix * Transpose(Nmatrix) * ColumnMatrix([R ! (-x[c]) : x in omega_vec]))[1][1] : c in [1..nu]];

// now adjust RHS
  cnt:=0;
    for c:=1 to nu do
      for j:=-porder to -1 do
        cnt +:= 1;
        v[cnt] +:= Coefficient(OmTZNNTOm[c],j);
      end for;
    end for;

// now the matrix whose "entries" are the principal parts of gammaFil_basis and of Omega[g+j]

  rows:=[];

    for i:=1 to g do
      row:=[];
      cnt:=0;
      for j:=1 to nu do
        for k:=-porder to -1 do
          cnt +:= 1;
          row[cnt]:=Coefficient(Omega_vec[i+g][j],k); // coefficients of principal part of Omegax_{i+g} at jth point at infinity
        end for; 
      end for;
      Append(~rows,row);
    end for;

    for i:=1 to #gammaFil_num do
      row:=[];
      cnt:=0;
      for j:=1 to nu do
        for k:=-porder to -1 do
          cnt +:= 1;
          row[cnt]:=Coefficient(gammaFil_basis[i][j],k); // coefficients of principal part of gammaFil_basis at jth point at infinity
        end for; 
      end for;
      Append(~rows,row);
    end for;

  suc,sol := IsConsistent(Matrix(rows),Vector(v));
  if suc then
    done:=true;
  else // if no success, increase the degree in x
    RRoffset:=RRoffset+1;
  end if;
    
end while;

printf "RRoffset is %o \n", RRoffset;

// read off bFil from solution

  bFil:=[];
  for i:=1 to g do
    bFil[i] := K!sol[i];
  end for;

printf "The vector bFil is %o \n", bFil;

// read off gammaFil from solution

gammaFil := [ &+[ sol[i] * gammaFil_basis[i][c] : i in [1..#gammaFil_num]] : c in [1..nu]];

printf "The q-exp at 1st cusp of the function gammaFil is %o \n", gammaFil[1];

// Since I don't choose a base point b, I ignore the normalising condition gammaFil(b) = 0
