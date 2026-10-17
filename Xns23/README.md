Open magma in a folder where you can use Zywina's Modular.
Do the following steps in order:
  - load "input_Xns23_LMFDB.m";
  - load "compute_cusp_forms_prec1.m";
This computes the spaces of weight 2 and weight 12 modular forms of the chosen level.

In the next file, change prec to the required precision. We choose prec=80.
  - load "increase_prec.m";

In the next file, change k to the required numerator valuation. We choose k=9.
  - load "compute_basis_H1dR.m";
  - load "compute_symplectic_basis_all_cusps.m";

In the next file, change the prime p whose Hecke operator on H1dR you want to compute. We choose p=3.
  - load "compute_Hecke_all_cusps.m";
  - load "compute_Hodge_filt_Xns23.m";
