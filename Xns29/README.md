Open magma in a folder where you can use Zywina's Modular.
Do the following steps in order:
  - load "input_Xns29_LMFDB.m";
  - load "compute_cusp_forms_prec1.m";

This computes the spaces of weight 2 and weight 12 modular forms of the chosen level.

In the next file, change prec to the required precision. We choose prec = 50.
  - load "increase_prec.m";

In the next file, change k to the required numerator valuation. We choose k=24.
  - load "compute_basis_H1dR_oneholdiff.m";
  - load "compute_symplectic_basis_all_cusps.m";

In the next file, change the prime p whose Hecke operator on H1dR you want to compute. We choose p=3.
  - load "compute_Hecke_all_cusps.m";
  - load "compute_Hodge_filt_Xns29.m";
(Note: In the computation of the Hodge filtration, the log differentials we find are, in fact, not a basis of the full space of log differentials; that's why there is a "false" in the output. However, we do find eta in the subspace we computed.)
