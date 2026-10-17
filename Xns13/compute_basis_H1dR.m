// restore "formscomputed";

k := 9;

H1dRbasistime := Time();

// Setup for Laurent series expansions (the q-expansions of the cusp forms)
K<z> := CyclotomicField(N);
w := N; // width
M := prec; // precision for q-expansions
R<qw>:=LaurentSeriesRing(K,M);

// denominator is +/- Delta function
T := Time();
S12special_denom := FindModularFormsWithVanishingConditions(M12,[N : i in [1..nu]]);
denom := S12special_denom`F0[1][1]; // only q-exp at first cusp
denom_all_cusps := S12special_denom`F0[1]; // q-exp at all cusps
printf "Finding the cusp form with highest order of vanishing took %o seconds.\nIts q-expansion at the first cusp is %o \n", Time(T), denom;

// 2 Functions to convert mat or seq with entries in Q(zetaN) to matrices with entries in Q

function CoefficientMatrix(n, M)
    // Create the cyclotomic field
    K<z> := CyclotomicField(n);

    phi := Degree(K);
    Q := Rationals();
    N := Nrows(M);
    P := Ncols(M);

    // Initialize result: same number of rows, φ(n) times more columns
    Mcoeffs := ZeroMatrix(Q, N, P*phi);

    for i in [1..N] do
        for j in [1..P] do
            v := Eltseq(M[i,j]);  // coefficients in basis 1, ζ, ζ^2, ...
            for k in [1..#v] do
                Mcoeffs[i, (j-1)*phi + k] := v[k];
            end for;
        end for;
    end for;

    return Mcoeffs;
end function;

function CoefficientMatrixFromSequence(n, v)
    // Define the cyclotomic field
    K<z> := CyclotomicField(n);

    phi := Degree(K);
    Q := Rationals();
    m := #v;  // number of elements in the sequence

    // Initialize φ(n) × m matrix over Q
    M := ZeroMatrix(Q, phi, m);

    for j in [1..m] do
        coeffs := Eltseq(v[j]);  // coefficients of v_j in basis 1, ζ, ζ², ...
        for k in [1..#coeffs] do
            M[k, j] := coeffs[k];
        end for;
    end for;

    return M;
end function;

// function for residue of f * g d/dq at q = 0, i.e. the constant term of f*g
function ddqres(f,g)
	a := f*g;
	return Coefficient(a,0);
end function;

// now we find numerators F so that, with fdR := F/denom, we have fdR * (wt 2 form) dq/q is a 2nd kind differential
T := Time();
S12special_num := FindModularFormsWithVanishingConditions(M12,[k : i in [1..nu]]); // modular forms that have order at least k at every cusp
printf "Computing wt 12 forms with valuation at least %o took %o seconds.\n", k, Time(T);
wt12valkforms_all_cusps := S12special_num`F0; // expansions of these at all cusps
// fdR_candidates := [R! h / denom : h in wt12valkforms]; // q-exp of candidates for fdR at first cusp
fdR_candidates_all_cusps := [ [R! h[c] / denom_all_cusps[c] : c in [1..nu]] : h in wt12valkforms_all_cusps]; // q-exp of fdR's at all cusps
// res_candidates := [ [ddqres(fdR_candidates[i],cusp_forms[a]) : i in [1..#wt12valkforms]] : a in [1..g]]; // residue of fdR_cand[i] * cusp_form[a] at first cusp
res_candidates_all_cusps := [ [ [ ddqres(fdR_candidates_all_cusps[i][c],cusp_forms_all_cusps[a][c]) : c in [1..nu]] : i in [1..#wt12valkforms_all_cusps]] : a in [1..g]]; // residues at all cusps
res_matrix_list_all_cusps := [ [CoefficientMatrixFromSequence(N, [res_candidates_all_cusps[a][i][c] : i in [1..#wt12valkforms_all_cusps]]) : c in [1..nu]] : a in [1..g]]; // expand residues in Q-basis of Q(zeta_N)
res_matrix := res_matrix_list_all_cusps[1][1];
for a in [1..g] do
    for c in [1..nu] do
        if a+c gt 2 then
            res_matrix := VerticalJoin(res_matrix, res_matrix_list_all_cusps[a][c]);
        end if;
    end for;
end for;
Nspace := Nullspace(Transpose(res_matrix));
d := Dimension(Nspace);
B := Basis(Nspace);
fdR_cand_new_all_cusps := [ [ &+[ B[j][i] * fdR_candidates_all_cusps[i][c] : i in [1..#wt12valkforms_all_cusps] ] : c in [1..nu] ] : j in [2..d]];

printf "The numerator candidates have valuation at least %o at all cusps. There are %o of them.\n", k, #fdR_cand_new_all_cusps;

// list of fdR * (holom diff); first #fdR_cand_new entries are fdR * f_1 etc
seckindcand_all_cusps := &cat [ [ [ fdR_cand_new_all_cusps[i][c] * R!cusp_forms_all_cusps[a][c] : c in [1..nu] ] : i in [1..#fdR_cand_new_all_cusps] ] : a in [1..g]];

// different pole orders seem to increase our chances of finding a basis
seckindcand := [x[1] : x in seckindcand_all_cusps];
mat := Matrix(#seckindcand, g, [[Coefficient(seckindcand[i], j+k-N): j in [1..g]]: i in [1..#seckindcand]]); // entries are coefficients of seckindcand with negative q-exponent; lie in Q(zetaN)
Mcoeff := CoefficientMatrix(N, mat); // matrix with coeff in Q
a, b := EchelonForm(Mcoeff);
nu := Integers() ! nu;
seckindcand_new_all_cusps := [ [ &+ [ b[Integers() ! 1+j*nu][i] * seckindcand_all_cusps[i][c] : i in [1..#seckindcand] ] : c in [1..nu] ] : j in [0..(g-1)]]; // these g functions lie in the vector space spanned by seckindcand and have different, increasing pole orders

// Our candidates for bases of H1dR consist of elements of cusp_forms (of weight 2) and seckindcand_new-functions (as usual, we omit the dq/q)
// To test whether we indeed get a basis, we use cup products again

// function that calculates the cup product contribution (at the 1st cusp)
// we integrate the second argument and multiply with the first argument, then take its residue
CupModForms:=function(g,f); // def [w_g] u [w_f]
    A:=Integral(f/qw); // send f to integral of omega f
    B:=A*g;
    return Coefficient(B,0); // taking residue
end function;

// candidate basis
M2cat_all_cusps := cusp_forms_all_cusps cat seckindcand_new_all_cusps;
// compute cup product matrix
C := Matrix(2*g,2*g,[ [ &+[CupModForms(M2cat_all_cusps[i][c],M2cat_all_cusps[j][c]) : c in [1..nu] ] : j in [1..(2*g)] ] : i in [1..(2*g)]]);
printf "The cup product matrix is:\n%o \n", C;
printf "Its determinant is: %o \n", Determinant(C);
if Determinant(C) eq 0 then
    printf "The following forms DO NOT form a basis of H1dR: %o \n", M2cat_all_cusps;
else
    printf "A basis of H1dR is given by the classes of f(q) dq/q with the q-exp of f at all cusps given by: %o \n", M2cat_all_cusps;
end if;

// the stuff without "_all_cusps" is
fdR_cand_new := [x[1] : x in fdR_cand_new_all_cusps];
seckindcand_new := [ x[1] : x in seckindcand_new_all_cusps];
M2cat := [x[1] : x in M2cat_all_cusps];

printf "Total time for basis of H1dR is %o \n", Time(H1dRbasistime);