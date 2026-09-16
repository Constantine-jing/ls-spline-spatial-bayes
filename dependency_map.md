# Dependency map — toy example

Files needed to run `toy/run_toy.R`.

Load order (what `run_toy.R` must `source()`):
1. `ls_basis.R`
2. `spatial_utils.R`
3. `ls_interaction.R`
4. `gibbs_interaction.R` (last)

Optional, before step 4: `Rcpp::sourceCpp("ls_interaction_core.cpp")`

```
ls_basis.R ──────┐
                 ├──> ls_interaction.R ──┐
                 │                       ├──> gibbs_interaction.R
spatial_utils.R ─┼───────────────────────┘
                 │
ls_interaction_core.cpp ──(optional)──> ls_interaction.R
```

---

## ls_basis.R

**does:** Builds the Lancaster–Šalkauskas natural cubic spline design matrix for one or several covariates: chooses equally spaced knots, forms A, C and A⁻¹C, evaluates the Φ/Ψ basis functions, and returns the identified design W = (Φ + Ψ A⁻¹C) T with a closure for building the same design on new x.

**defines:**
- `ls_contrast_T(M)` — sum-to-zero contrast T (M × (M−1))
- `ls_choose_knots_equal(x, M)` — equally spaced knots on [min x, max x]
- `ls_AC(tau)` — natural-spline constraint matrices A, C and A⁻¹C
- `ls_phi_psi(x, tau)` — raw basis matrices Φ, Ψ (n × M)
- `ls_build_one_full(x, M, tau)` — full machinery for one covariate, incl. `design_new()`
- `ls_build_one_train(x, M, tau)` — lightweight train-time wrapper
- `ls_additive_build(X, M_vec, tau_list)` — block design for p covariates
- `ls_additive_design_new(X_new, objs, clip)` — matching design on new data
- `ls_tests()` — self-check (cardinality, partition of unity, contrast)

**packages:** none
**sources:** none
**uses from other files:** none

**used by the toy:** `ls_build_one_full` only. The general-p functions (`ls_build_one_train`, `ls_additive_build`, `ls_additive_design_new`) are used by the application scripts, not the toy.

---

## spatial_utils.R

**does:** Distance and Matérn correlation helpers for the spatial GP, plus a Cholesky-based linear solve.

**defines:**
- `pairdist(coords)` — n × n Euclidean distance matrix
- `pairdist_cross(A, B)` — cross-distance matrix between two location sets
- `matern_cor(D, rho, nu)` — Matérn correlation matrix from distances D
- `matern_cor_cross(D, rho, nu)` — same for a cross-distance matrix
- `solve_chol(U, b)` — solve (UᵀU)x = b given upper Cholesky U

**packages:** none
**sources:** none
**uses from other files:** none

**used by the toy:** `matern_cor` only.

---

## ls_interaction_core.cpp

**does:** Rcpp/Eigen implementations of computational hot spots: Khatri–Rao product for the interaction design, RW2 penalty builders, quadratic forms, Cholesky solves, MVN draws, whitening, block prior precision.

**defines (exported to R):**
- `khatri_rao_cpp(W_u, W_v)` — row-wise Kronecker product
- `rw2_penalty_cpp(d)`, `rw2_2d_penalty_raw_cpp(M_u, M_v)` — RW2 penalties
- `quad_form_cpp`, `chol_solve_cpp`, `sample_mvn_chol_cpp`, `log_quad_lik_cpp`, `whiten_cpp` — linear-algebra helpers
- `build_block_prior_precision_cpp(d_vec, tau2_s, kappa2, eps_ridge)` — prior precision Q0

**packages:** `Rcpp`, `RcppEigen`; needs a C++ compiler (Rtools on Windows, Xcode tools on Mac)
**sources:** none
**uses from other files:** none

**used by the toy:** nothing at present. `ls_interaction.R` calls `khatri_rao_cpp` only when `use_cpp = TRUE`, and `run_interaction_poc` uses the default `use_cpp = FALSE`. The other nine functions are not called from any R file; the sampler does its linear algebra in base R. The toy runs identically with or without this file.

---

## ls_interaction.R

**does:** Builds the tensor-product interaction design W_uv = W_u ⊙ W_v, the identified 2D RW2 penalty (T_u ⊗ T_v)ᵀ (K_u ⊗ I + I ⊗ K_v)(T_u ⊗ T_v), and the optional design-level ANOVA orthogonalization that residualizes W_uv against [1, W_u, W_v] (the Prop. 4.1 fix).

**defines:**
- `build_rw2_penalty_1d(d)` — 1D RW2 penalty K = D2ᵀD2
- `build_rw2_2d_penalty(K_u_raw, K_v_raw, T_u, T_v)` — identified 2D Kronecker-sum penalty
- `khatri_rao_rowwise_R(W_u, W_v)` — pure-R row-wise Kronecker product
- `anova_orthogonalize_W_uv(W_uv, W_u, W_v, ridge)` — residualize W_uv against main effects
- `ls_build_interaction(obj_u, obj_v, use_cpp, orthogonalize)` — main constructor
- `ls_interaction_design_new(X_u_new, X_v_new, recipe, clip)` — matching W_uv on new data
- `ls_build_all_interactions(full_objs, orthogonalize)` — all p(p−1)/2 pairs
- `ls_assemble_full_design(W_main, int_list, n)` — intercept + mains + interactions with column maps
- `ls_interaction_tests(n, M, seed, ...)` — self-check

**packages:** none
**sources:** none

**uses from other files:**
- `ls_basis.R`: inputs `obj_u`, `obj_v` are outputs of `ls_build_one_full`
- `ls_interaction_core.cpp`: `khatri_rao_cpp`, only if `use_cpp = TRUE`

**used by the toy:** `build_rw2_penalty_1d`, `ls_build_interaction`

---

## gibbs_interaction.R

**does:** The p = 2 collapsed Gibbs sampler for y = μ + f₁(x₁) + f₂(x₂) + f₁₂(x₁, x₂) + b(s) + ε with the Matérn GP b marginalized out: MVN block draw for η = (μ, β₁, β₂, β₁₂), MH steps for σ², τ², ρ, conjugate IG draws for the smoothing variances, posterior mean of b. Also contains the toy driver.

**defines:**
- `.has_cpp` — load-time flag; only prints a message
- `build_interaction_prior_precision(...)` — prior precision Q0 for η
- `gibbs_interaction_sampler(y, H, D, nu, col_map_main, K_main_list, col_map_int, K_int_list, ...)` — the sampler
- `run_interaction_poc(n, M, n_iter, n_burn, seed, verbose)` — simulate → build designs → fit → RMSEs (**the toy**)
- `plot_interaction_poc(poc, outfile)` — PDF of fitted vs true curves/surfaces and trace plots

**packages:** none (`Rcpp` is checked but never required)
**sources:** none — the header says to source the other files first, but does not do it

**uses from other files:**
- `ls_basis.R`: `ls_build_one_full`
- `ls_interaction.R`: `build_rw2_penalty_1d`, `ls_build_interaction`
- `spatial_utils.R`: `matern_cor`

---

## Decisions

- **C++ in the toy.** Options: (a) mandatory — `run_toy.R` compiles the C++ and compares R vs C++ designs; a stranger without a compiler cannot run it. (b) automatic — try to compile; if it works, compare and use C++; if not, continue in R with one message. No user switch either way.

## Later (comment cleanup pass)

- `gibbs_interaction.R` and `ls_interaction_core.cpp` refer to `gibbs_stage_c_full.R`, which is not in the public repo. Reword or drop.
- Header of `gibbs_interaction.R` lists `khatri_rao_rowwise_R` from `ls_interaction.R`; actual dependency is `build_rw2_penalty_1d`.
- `gibbs_interaction_sampler` has a "DIAGNOSTIC HOOK (May 2026)" argument block. Decide whether it stays.
- Mixed line endings (CRLF/LF). Fix once with `.gitattributes` containing `* text=auto`; commit with `.gitignore`.
