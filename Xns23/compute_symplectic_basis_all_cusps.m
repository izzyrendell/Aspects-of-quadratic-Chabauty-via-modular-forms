Timeforsympl := Time();

// function that does one step in symplectic Gram-Schmidt
function SymplecticGramSchmidtStep(C)
    g := Integers() ! (Nrows(C)/2);
    test := true;
    k := 1;
    while test do
        if C[1,g+k] ne 0 then
            test := false;
        else
            k := k+1;
        end if;
    end while;

    a := C[1,g+k];
    S := MatrixAlgebra(Rationals(),2*g) ! 0;
    S[1,1] := 1; // first column of S is unit vector with 1 as first entry
    S[g+k,g+1] := 1/a; // g+1-st column of S is unit vector with 1 in (g+k)-entry, divided by a
    for i in [2..g] do
        S[1,i] := -C[i,g+k]/a; // i-th column of S is unit vector at i, modified at the first entry
        S[i,i] := 1;
    end for;
    for i in [(g+2)..2*g] do
        if i-g gt k then 
            S[i,i] := 1; // i-th column of S is unit vector at i, modified at first and g+k-th entry
            S[1,i] := -C[i,g+k]/a;
            S[g+k,i] := C[i,1]/a;
        else
            S[i-1,i] := 1; // i-th column of S is unit vector at i-1, modified at first and g+k-th entry
            S[1,i] := -C[i-1,g+k]/a;
            S[g+k,i] := C[i-1,1]/a;
        end if;
    end for;

    return S;

end function;

// function that removes certain rows and columns from a matrix
RemoveRowsCols := function(M, inds)
    // inds = sequence of row/column indices to remove
    keep := [i : i in [1..Nrows(M)] | not i in inds];
    return Submatrix(M, keep, keep);
end function;

// another function that inserts "unit vector" rows and columns
InsertPadBlocks := function(M, g)

    n := Nrows(M);
    r := Integers() ! (n/2);

    if r gt g then
        error "Need r <= g";
    end if;

    R := BaseRing(M);
    N := ZeroMatrix(R, 2*g, 2*g);

    off := g - r;

    // -------------------------------------------------
    // Insert identity vectors in both halves
    // -------------------------------------------------

    for i in [1..off] do
        N[i,i] := 1;
        N[g + i, g + i] := 1;
    end for;

    // -------------------------------------------------
    // Embed the 2r x 2r block matrix
    // -------------------------------------------------

    for i in [1..r] do
        for j in [1..r] do

            // A block
            N[off + i, off + j] := M[i,j];

            // B block
            N[off + i, g + off + j] := M[i, j + r];

            // C block
            N[g + off + i, off + j] := M[i + r, j];

            // D block
            N[g + off + i, g + off + j] := M[i + r, j + r];

        end for;
    end for;

    return N;

end function;

// executing entire symplecticGramSchmidt:
// if f = [f[i] : i in [1..2g]] is a basis and C denotes the matrix of an alternating pairing
// i.e. f[i] cup f[j] = C[i,j]
// then the output S of this function is the base change matrix to a symplectic basis
// i.e. if b[j] = sum_i S[i,j] f[i]
// then b[i] cup b[j] = (S^T * C * S)[i,j] has the "usual symplectic block-diagonal form" [0,I, -I,0]
function SymplecticBCFromIsotropicHalf(C)
    dimension1 := Nrows(C);
    dimension2 := Ncols(C);
    if dimension1 ne dimension2 then
        return false;
    else
        if dimension1 mod 2 ne 0 then
            return false;
        else
            g := Integers() ! (dimension1/2);
            if &and [ &and [C[i,j] eq 0 : j in [1..g]] : i in [1..g]] then
                Sfull := MatrixAlgebra(Rationals(),2*g) ! 1;
                test := g;
                while test gt 0 do
                    S := SymplecticGramSchmidtStep(C);
                    Stilde := InsertPadBlocks(S, g);

                    Sfull := Sfull * Stilde;
                    Cnew := Transpose(S) * C * S;
                    C := RemoveRowsCols(Cnew, [1,test+1]);
                    test := test - 1;
                end while;

                return Sfull;
            else
                return false;
            end if;
        end if;
    end if;
end function;

BC_M2cat_symp := SymplecticBCFromIsotropicHalf(C);
printf "The base change matrix from our given basis M2cat of H1dR to a symplectic basis is:\n%o \n", BC_M2cat_symp;
symp_basis_all_cusps := [ [ &+ [BC_M2cat_symp[i,j] * M2cat_all_cusps[i][c] : i in [1..2*g]] : c in [1..nu] ] : j in [1..2*g]];
printf "The q-exp at all cusps of a symplectic basis are: %o \n", symp_basis_all_cusps;
symp_basis_one_cusp := [x[1] : x in symp_basis_all_cusps];

printf "Total time for computing symplectic basis is %o\n", Time(Timeforsympl);