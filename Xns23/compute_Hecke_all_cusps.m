// We want the matrix of the Hecke operator T_p w.r.t. the basis of H1dR computed before

TimeforHecke := Time();

p := 3; // want to calculate T_p
printf "We calculate T_p with p = %o\n", p;

// we first need to find a 2x2 matrix alpha with Z-coeff with det = p, whose reduction mod N lies in G, and whose 2nd column is divisible by p
// often, one can do this quickly by hand; to automate the process, need a few function about lifting matrices

// need a function that ensures we find a coprime lift of a tuple of integers mod N
function CoprimeLift(a,b,N)
    if IntegerRing(N) ! b ne IntegerRing(N) ! 0 then
        t := 0;
        while GCD(a+N*t, b) ne 1 do
            t +:= 1;
        end while;

        return a + N*t, b;

    elif IntegerRing(N) ! a ne IntegerRing(N) ! 0 then
        t := 0;
        while GCD(a, b+N*t) ne 1 do
            t +:= 1;
        end while;

        return a, b+N*t;

    end if;

    error "Both numbers divisible by N";
end function;

// a function that lifts a 2x2 matrix A with Z/N coeff with det = p mod N (p any integer not divisible by N) to a 2x2 Z-matrix of det = p whose reduction mod N is A
function LiftDetPmatrix(A, N, pdet)
    if Determinant(A) ne IntegerRing(N) ! pdet then
        return false;
    else
        a := Integers()!A[1,1];
        b := Integers()!A[1,2];
        c := Integers()!A[2,1];
        d := Integers()!A[2,2];

        c, d := CoprimeLift(c,d,N); // ensure gcd(c,d) = 1
        k := (a*d - b*c - pdet) div N;

        e, u, v := XGCD(c, d);
        x := -v*k;
        y := u*k;

        M := Matrix(Integers(),2,2,[a+N*x,b+N*y, c,d]); 
        if Determinant(M) eq pdet then // double check that det M = pdet
            return M;
        end if;
    end if;
end function;

// a function that finds a 2x2 matrix with Z-coeff with det = pdet whose reduction lies in G
function LiftDetP(G, N, pdet)
    det_p_mat := [ g : g in G | Determinant(g) eq pdet ];
    for i in [1..#det_p_mat] do
        A := det_p_mat[i];
        return LiftDetPmatrix(A, N, pdet);
    end for;
end function;

A := LiftDetP(G, N, p);
printf "The 2x2 Z-matrix with det %o that reduces to G mod %o is A = \n%o \n", p, N, A;

// additionally want the second column to be divisible by p; so let's left and right multiply by matrices from Gamma := matrices in SL(2,Z) that reduce to G mod N

Gdet1 := [ g : g in G | Determinant(g) eq 1 ];
Gamma_gens := [LiftDetPmatrix(A, N, 1) : A in Gdet1]; // some elements of Gamma

function FindSecondColumnDivisibleByp(A, gensH, p) // want to find H1 and H2 such that H1*A*H2 has second column divisible by p

    Z := Integers();
    F := GF(p);

    // Reduce A and generators mod p

    Ap := ChangeRing(A, F);
    gensHp := [ ChangeRing(g, F) : g in gensH ];

    // Build H mod p

    Gp := SL(2, F);
    Hp := sub< Gp | gensHp >;

    // Enumerate Hp and store words

    Idp := Identity(Gp);

    Hplist := [ Idp ];
    Words := [ [] ];

    i := 1;

    while i le #Hplist do

        h := Hplist[i];
        w := Words[i];

        for j in [1..#gensHp] do

            newh := h * gensHp[j];

            if newh notin Hplist then
                Append(~Hplist, newh);
                Append(~Words, w cat [j]);
            end if;

        end for;

        i +:= 1;

    end while;

    // Now test all pairs

    for i1 in [1..#Hplist] do
        for i2 in [1..#Hplist] do

            h1 := Hplist[i1];
            h2 := Hplist[i2];

            M := h1 * Ap * h2;

            if M[1,2] eq 0 and M[2,2] eq 0 then

                // Lift to integer matrices

                H1 := Identity(MatrixAlgebra(Z,2));
                H2 := Identity(MatrixAlgebra(Z,2));

                for j in Words[i1] do
                    H1 *:= gensH[j];
                end for;

                for j in Words[i2] do
                    H2 *:= gensH[j];
                end for;

                return true, H1, H2;

            end if;

        end for;
    end for;

    return false, _, _;

end function;

function FindUprightDivisibleByp(A, gensH, p) // similar to above, but finds H1 such that H1*A has upper right entry divisble by p

    Z := Integers();
    F := GF(p);

    // Reduce A and generators mod p

    Ap := ChangeRing(A, F);
    gensHp := [ ChangeRing(g, F) : g in gensH ];

    // Build H mod p

    Gp := SL(2, F);
    Hp := sub< Gp | gensHp >;

    // Enumerate Hp and store words

    Idp := Identity(Gp);

    Hplist := [ Idp ];
    Words := [ [] ];

    i := 1;

    while i le #Hplist do

        h := Hplist[i];
        w := Words[i];

        for j in [1..#gensHp] do

            newh := h * gensHp[j];

            if newh notin Hplist then
                Append(~Hplist, newh);
                Append(~Words, w cat [j]);
            end if;

        end for;

        i +:= 1;

    end while;

    // Now test all pairs

    for i1 in [1..#Hplist] do

            h1 := Hplist[i1];

            M := h1 * Ap;

            if M[1,2] eq 0 then

                // Lift to integer matrices

                H1 := Identity(MatrixAlgebra(Z,2));

                for j in Words[i1] do
                    H1 *:= gensH[j];
                end for;

                return true, H1;

            end if;
    end for;

    return false, _;

end function;

flag, H1, H2 := FindSecondColumnDivisibleByp(A, Gamma_gens, p);
alpha := H1 * A * H2; // in [CKL] language: alpha
printf "The suitable 2x2 Z-matrix with det %o that reduces to G mod %o is alpha = \n%o \n", p, N, alpha;


// The Hecke operator is the double coset operator [Gamma_G alpha Gamma_G] 
// Calculation of coset reps for Gamma_G \ Gamma_G alpha Gamma_G; this will be the list cosetrepsHecke


// the second column of alpha is divisible by p, so really alpha = alpha1 * alpha2: (in paper, these are alpha1' and alpha2')
alpha1 := MatrixAlgebra(Integers(),2) ! [alpha[1,1],Integers() ! alpha[1,2]/p, alpha[2,1],Integers() ! alpha[2,2]/p];
alpha2 := MatrixAlgebra(Integers(),2) ! [1,0, 0,p];
// we thus see that alpha^-1 SL(2,Z) alpha \cap SL(2,Z) = alpha2^-1 SL(2,Z) alpha2 \cap SL(2,Z) = Gamma^0(p)
// coset reps for this group are given by the "usual" [Diamond-Shurman, p.170, §5.2]
cosetrepsforGamma0 := [MatrixAlgebra(Integers(),2) ! [1,j, 0,1] : j in [0..p-1]] cat [MatrixAlgebra(Integers(),2) ! [p,p-1, 1,1]]; // u_j in the paper
printf "Coset reps beta for Gamma^0(p) in SL(2,Z) are: %o \n", cosetrepsforGamma0;

// For each beta in cosetrepsforGamma0, find gamma in Gamma0 such that gamma*beta lies in Gamma_G (i.e. reduces mod N to G)
// gamma*beta has determinant 1, so lies in Gdet1
// in general, gamma*beta could be any matrix in Gdet1; in practice, we'll mostly just get that gamma*beta reduces to id mod N
// ansatz: gamma is [1,p*b, 0,1], then gamma*beta will have shape [1,p*b+j, 0,1] and choose b so this is id; this works easily for the first p matrices
listofgamma := [ MatrixAlgebra(Integers(),2) ! [1,p*Integers() ! (IntegerRing(N) ! -(j-1)/p), 0,1]: j in [1..p]]; 
// for last cosetrepsforGamma0 matrix, the same ansatz may fail

Glastrow := [ g : g in Gdet1 | g[2,1] eq 1 and g[2,2] eq 1];
if Glastrow ne [] then // same ansatz works for gamma works and leads to gamma*beta = [p(b+1),bla, 1,1]
    u := Glastrow[1]; 
    V := MatrixAlgebra(Integers(),2) ! [1,p*Integers() ! (IntegerRing(N) ! u[1,1]/p -1), 0,1];
    listofgamma := listofgamma cat [V];
else // if ansatz fails, need to do something else:
    uinverse := (Gr ! cosetrepsforGamma0[p+1])^-1; // gamma will be a lift of uinverse to SL(2,Z) that is also lower triangular mod p
    V := LiftDetPmatrix(uinverse, N, 1); // first an arbitrary lift of uinverse to SL(2,Z)
    flag, H1 := FindUprightDivisibleByp(V, Gamma_gens, p);
    V := H1 * V;
    listofgamma := listofgamma cat [V];
end if;

printf "A list of elts gamma of Gamma^0(p) such that gamma*beta lies in Gamma_G: %o \n", listofgamma; // gamma_j in the paper
check1 := [Gr ! ChangeRing(listofgamma[i] * cosetrepsforGamma0[i],IntegerRing(N)) in G : i in [1..p+1]]; // if this gives all true, we've done it correctly
printf "gamma*beta does indeed lie in Gamma_G: %o \n", check1;
check2 := [listofgamma[i][1,2] mod p eq 0: i in [1..p+1]];
printf "All elements of the list of gammas have upper right entry divisible by p: %o \n", check2;
check3 := [Determinant(listofgamma[i]) eq 1: i in [1..p+1]];
printf "All elements of the list of gammas have determinant 1: %o \n", check3;

cosetrepsforGammaG := [listofgamma[i] * cosetrepsforGamma0[i] : i in [1..p+1]];
// coset reps for Gamma_G \ Gamma_G alpha Gamma_G (which is what I need for T_p) are then given by
cosetrepsHecke := [alpha * cosetrepsforGammaG[i] : i in [1..p+1]]; // alpha_j in the paper
printf "Matrices for the Hecke operator: %o \n", cosetrepsHecke;
check4 := [Determinant(cosetrepsHecke[i]) eq p : i in [1..p+1]];
printf "All matrices here have determinant %o: %o \n", p, check4;

// For each c = cosetrepsHecke[i](infinity), need to find a matrix gamma in Gamma and cusp sigma[k](infinity) such that c = gamma sigma[k](infinity)

// matrices representing the cusps from Zywina's Modular Code
sigma := [MatrixAlgebra(Integers(),2) ! s : s in X`cusps];

// a function that, on input an integer matrix alpha with det = p, outputs an SL(2,Z)-matrix epsilon s.t. alpha(infty) = epsilon(infty)
SLmatrix := function(alpha)
    a := alpha[1][1];
    b := alpha[1][2];
    c := alpha[2][1];
    d := alpha[2][2];
    g, m, n := ExtendedGreatestCommonDivisor(a,c);
    h := GreatestCommonDivisor(b,d);
    if g eq 1 and h eq p then
        return MatrixAlgebra(Integers(),2) ! [a,b/p, c,d/p];
    elif g eq 1 and h eq 1 then
        return MatrixAlgebra(Integers(),2) ! [a,-n, c,m];
    elif g eq p then
        return MatrixAlgebra(Integers(),2) ! [a/p,alpha[1][2], c/p,alpha[2][2]];
    end if;

    error "matrix does not have determinant p"; // if the gcd is anything but 1 and p, we made a mistake earlier
end function;

// for each j and k, we find an SL(2,Z)-matrix epsilon s.t. cosetrepsHecke[j] * sigma[k] (infty) = epsilon[j][k](infty)
epsilon := [[SLmatrix(cosetrepsHecke[j]*sigma[k]) : k in [1..nu]] : j in [1..p+1]];
// note that, in particular, this means that delta[j][k] := epsilon[j][k]^-1 * cosetrepsHecke[j] * sigma[k] is upper triangular
delta := [[epsilon[j][k]^-1 * cosetrepsHecke[j]*sigma[k] : k in [1..nu]] : j in [1..p+1]];
check5 := [[delta[j][k][2,1] eq 0 : k in [1..nu]] : j in [1..p+1]];
printf "All matrices delta are upper triangular: %o \n", check5;

// a function that, on input j and k, finds gamma[j][k] in Gamma, an index[j][k] and an upper triangular matrix beta = sgn*[1,h,0,1] s.t.
// gamma*sigma[index]*beta = epsilon[j][k]
// the function simply tries all values for index and small values for h
Index_Gamma_UpTriang := function(j,k)
    for eps in [1,-1] do
        for h in [1..100] do
            b := eps * MatrixAlgebra(Integers(),2) ! [1,h, 0,1];
            for i in [1..nu] do
                gamma := epsilon[j][k] * b * (sigma[i]^-1);
                if Gr ! gamma in G then
                    return i, gamma, b^-1;
                end if;
            end for;
        end for;
    end for;

    error "No h found in range [1..100]";
end function;

index_gamma_beta := [ [ <index, gamma, beta> where index, gamma, beta := Index_Gamma_UpTriang(j,k) : k in [1..nu] ] : j in [1..p+1] ];
index := [ [index_gamma_beta[j][k][1] : k in [1..nu]] : j in [1..p+1] ];
// gamma := [ [index_gamma_beta[j][k][2] : k in [1..nu]] : j in [1..p+1] ]; // we do not need the gammas later
beta := [ [index_gamma_beta[j][k][3] : k in [1..nu]] : j in [1..p+1] ];
betadelta := [ [beta[j][k] * delta[j][k] : k in [1..nu]] : j in [1..p+1] ]; // those are the matrices with whom we slash later; beta_{j,k} in the paper



// temporarily extend coefficients to a larger cyclotomic field
L<zetawp> := CyclotomicField(w*p);
phi:= hom< K -> L | zetawp^p >;
S<qwp> := LaurentSeriesRing(L,p*M);
psi := hom< R -> S | qwp^p>;

// a function that, on input a single qw-expansion f and an upper triangular matrix B, outputs the qw-expansion of f | B
slashUpTriang := function(f,B)
    if B[2][1] eq 0 then
        a := B[1][1]; b := B[1][2]; d := B[2][2];
        fslashB := 0;
        if d eq p then
            for n in [Valuation(f)..(AbsolutePrecision(f)-1)] do
                fslashB := fslashB + a/d * Coefficient(f, n) * zetawp^(b*n) * qwp^(a*n);
            end for;
            return fslashB;
        elif d eq 1 then
            for n in [Valuation(f)..(AbsolutePrecision(f)-1)] do
                fslashB := fslashB + a/d * Coefficient(f, n) * z^(b*n) * qw^(a*n);
            end for;
            return psi(fslashB);
        end if;
    end if;

    error "matrix not upper triangular or bottom right entry neither 1 nor p";
end function;

// a function that, on input a (weakly) modular form f given as a sequence of nu qw-expansions, outputs T_p(f) as a sequence of nu qw-expansions
// atm, this uses the matrices and indices computed above
Hecke := function(f)
    Tpf := [&+ [ slashUpTriang(f[index[j][k]], betadelta[j][k]) : j in [1..p+1] ] : k in [1..nu]];
    return Tpf;
end function;

// a function that, on input an expansion in S, returns the expansion of f in R up to no_of_coeffs_in_R coefficients from the (val_qwp(f)/p) in S
ChangeLaurentSeriesFromStoR := function(f);
    coeffs_in_K := [K!Coefficient(f, j): j in [Valuation(f)..(AbsolutePrecision(f)-1)]]; // first making coefficients in L into K
    return &+[coeffs_in_K[1+p*i]*qw^(Valuation(f)/p+i): i in [0..Integers() ! (RelativePrecision(f)/p-1)]]; // reconstructing expansion into R from S
end function;

// To get the Hecke matrix on H1dR from this, use the cup product again
Hecke_symp_all_cusps := [Hecke(symp_basis_all_cusps[i]) : i in [1..2*g]]; // Hecke image of symplectic basis
Hecke_symp_K_all_cusps := [ [ ChangeLaurentSeriesFromStoR(Hecke_symp_all_cusps[i][k]) : k in [1..nu]] : i in [1..2*g] ]; // in qw, coeff in K
Hecke_symp_K_one_cusp := [Hecke_symp_K_all_cusps[i][1] : i in [1..2*g]];
printf "For the symplectic basis b_i of H1dR, the q-exp of T_p(b_i) at the first cusp are: %o \n", Hecke_symp_K_one_cusp;

// matrix with entries T_p(b_i) cup b_j
M_matrix := Matrix(2*g, 2*g,[[ &+[CupModForms(Hecke_symp_K_all_cusps[i][c],symp_basis_all_cusps[j][c]) : c in [1..nu] ] : j in [1..2*g]] : i in [1..2*g]]);
// for Hecke-matrix T[i,j] we have M[i,j] = - T[j+g,i] or M[i,j] = T[j-g,i] depending on whether j is <= g or not
T := Matrix(2*g, 2*g, [[Rationals() ! 0 : i in [1..2*g]] : j in [1..2*g]]);
for i in [1..2*g] do
    for j in [1..g] do
        T[j,i] := M_matrix[i,j+g];
        T[j+g,i] := -M_matrix[i,j];
    end for;
end for;
printf "The matrix for T_p on H1dR w.r.t. the symplectic basis b_i is:\n%o \n", T;
printf "Its characteristic polynomial is:\n%o \n", Factorization(CharacteristicPolynomial(T));

printf "Total time for computing Hecke matrix is %o seconds\n", Time(TimeforHecke);