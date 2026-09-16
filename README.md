# Bayesian Spatial Additive Regression with Lancaster–Šalkauskas Splines and Matérn Gaussian Process Random Effects

Mengyan Jing and Sounak Chakraborty, University of Missouri

Code accompanying the JRSS-C paper. Under construction.

## Repository structure

- `R/` — model code: `ls_basis.R`, `spatial_utils.R`, `ls_interaction.R`, `gibbs_interaction.R`, and `ls_interaction_core.cpp` (optional Rcpp speedups)
- `toy/` — runnable example on simulated data: `run_toy.R` (pure R) and `run_toy_cpp.R` (with the C++ hot loops)
- `scripts/` — application runs for the paper (to follow)
- `data/` — no data is stored; a README explains how to query CDC WONDER and the ACS tables