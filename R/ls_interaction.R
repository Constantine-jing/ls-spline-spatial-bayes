# ============================================================
# ls_interaction.R
#
# Tensor product LS interaction basis for bivariate surface f_{u,v}(x_u, x_v)
#
# ls_build_interaction() has an `orthogonalize = TRUE` option that
# residualizes W_uv against [1, W_u, W_v] at construction time, so that
#
#     1' W_uv = 0,   W_u' W_uv = 0,   W_v' W_uv = 0
#
# holds exactly at the OBSERVATIONS. This is a data-level ANOVA orthogonality
# condition and is strictly stronger than the coefficient-grid ANOVA enforced
# by T_u (x) T_v (Proposition 4.1). It is needed because the cardinal LS basis
# has columns whose data-mean is generally nonzero, so the Khatri-Rao
# product W_u (*) W_v has columns whose inner products with [1, W_u, W_v] are
# not zero -- producing column-space overlap between main-effect and
# interaction blocks and contaminating beta_main posteriors via beta_uv
# posterior uncertainty.
#
# After residualization:
#   * The likelihood y = ... + W_uv beta_uv + ... is preserved exactly under a
#     silent reparametrization where M_uv A_uv beta_uv (the projected piece) is
#     absorbed into the intercept and main-effect coefficients. With the existing
#     flat intercept prior (kappa2=1e6) and adaptive main-effect smoothing, this
#     absorption is free -- those priors easily accommodate a constant shift.
#   * The prior on beta_uv via K_uv stays exactly the same. K_uv penalizes
#     roughness of the underlying theta-surface; the orthogonalization changes
#     only the parameterization of the interaction column space, not the prior
#     measure on it.
#
# Depends on: ls_basis.R   (ls_build_one_full, ls_contrast_T, etc.)
#
# orthogonalize = FALSE (the default) skips the residualization.
# ============================================================


# ============================================================
# build_rw2_penalty_1d(d)
# Standard 1D RW2 penalty, size d x d.
# ============================================================
build_rw2_penalty_1d <- function(d) {
  if (d < 3) return(matrix(0, d, d))
  D2 <- matrix(0, d - 2, d)
  for (i in seq_len(d - 2)) {
    D2[i, i]     <-  1
    D2[i, i + 1] <- -2
    D2[i, i + 2] <-  1
  }
  crossprod(D2)
}


# ============================================================
# build_rw2_2d_penalty(K_u_raw, K_v_raw, T_u, T_v)
#
# 2D Kronecker-sum RW2 penalty on the raw theta_uv grid, then transformed
# to the identified (M_u-1)(M_v-1) parametrization via T_u (x) T_v.
# Unchanged from original.
# ============================================================
build_rw2_2d_penalty <- function(K_u_raw, K_v_raw, T_u, T_v) {
  M_u <- nrow(K_u_raw)
  M_v <- nrow(K_v_raw)

  K_uv_raw <- kronecker(K_u_raw, diag(M_v)) + kronecker(diag(M_u), K_v_raw)

  T_uv <- kronecker(T_u, T_v)
  K_uv_beta <- t(T_uv) %*% K_uv_raw %*% T_uv

  list(K_uv_raw = K_uv_raw, K_uv_beta = K_uv_beta, T_uv = T_uv)
}


# ============================================================
# khatri_rao_rowwise_R(W_u, W_v)
#
# Row-wise Kronecker product:  result[i, ] = W_u[i, ] (x) W_v[i, ]
# Output: n x ((M_u-1)*(M_v-1))
# Unchanged from original.
# ============================================================
khatri_rao_rowwise_R <- function(W_u, W_v) {
  n   <- nrow(W_u)
  d_u <- ncol(W_u)
  d_v <- ncol(W_v)
  stopifnot(nrow(W_v) == n)

  W_u_rep <- W_u[, rep(seq_len(d_u), each  = d_v), drop = FALSE]
  W_v_rep <- W_v[, rep(seq_len(d_v), times = d_u), drop = FALSE]
  W_u_rep * W_v_rep
}


# ============================================================
# anova_orthogonalize_W_uv(W_uv, W_u, W_v, ridge = 1e-10)
#
# Residualize each column of W_uv against [1_n, W_u, W_v] at the observations.
# Returns:
#   W_uv_ortho : n x q  where q = ncol(W_uv)
#   A_uv       : ((1 + d_u + d_v)) x q  projection coefficients
#                (W_uv = M_uv %*% A_uv + W_uv_ortho with W_uv_ortho perp M_uv)
#   M_inv_MtM  : (M_uv' M_uv)^-1   stored for prediction at new X
#
# Algorithm: compute QR of M_uv = [1, W_u, W_v] for stable projection.
# A small ridge is used only to stabilize the (M' M)^-1 storage when M_uv has
# near-rank-deficient columns (it should be full rank in practice with M >= 4).
# ============================================================
anova_orthogonalize_W_uv <- function(W_uv, W_u, W_v, ridge = 1e-10) {
  n <- nrow(W_uv)
  stopifnot(nrow(W_u) == n, nrow(W_v) == n)

  M_uv <- cbind(1, W_u, W_v)   # n x (1 + d_u + d_v)

  # Solve A_uv = (M' M)^-1 M' W_uv via QR, more stable than explicit inverse.
  # We need A_uv for prediction-at-new-X (see ls_interaction_design_new).
  qr_M <- qr(M_uv)
  A_uv <- qr.coef(qr_M, W_uv)        # (1 + d_u + d_v) x q
  W_uv_ortho <- W_uv - M_uv %*% A_uv

  # Also store (M' M + ridge I)^-1 for prediction (alternate path).
  # We don't actually need this if A_uv is stored AND prediction uses the
  # *training-time* M_uv structure -- which is correct, see prediction code.
  # But we keep the QR object so prediction can reuse the same projection.

  list(
    W_uv_ortho = W_uv_ortho,
    A_uv       = A_uv,
    M_uv       = M_uv,             # stored for prediction at new X
    qr_M       = qr_M
  )
}


# ============================================================
# ls_build_interaction(obj_u, obj_v, use_cpp = FALSE, orthogonalize = FALSE)
#
# Build everything needed for the interaction surface f_{u,v}(x_u, x_v):
#
#   W_uv  : n x (M_u-1)*(M_v-1)  identified design matrix
#           If orthogonalize=TRUE, this is the design-level ANOVA-orthogonalized
#           version: W_uv_ortho = (I - P_M) (W_u (*) W_v) where P_M is the
#           projector onto [1, W_u, W_v].
#   K_uv  : (M_u-1)*(M_v-1) x (M_u-1)*(M_v-1)  identified 2D RW2 penalty
#           (unchanged by orthogonalization -- still the prior on beta_uv)
#   T_uv  : M_u*M_v x (M_u-1)*(M_v-1)  combined contrast
#   recipe: stores everything needed for prediction at new (X_u, X_v)
#
# orthogonalize = TRUE  : recommended (Project 1 fix for basis-correlation)
# orthogonalize = FALSE : original behavior (default for backward compatibility)
# ============================================================
ls_build_interaction <- function(obj_u, obj_v, use_cpp = FALSE,
                                 orthogonalize = FALSE) {
  W_u <- obj_u$W
  W_v <- obj_v$W
  stopifnot(!is.null(W_u), !is.null(W_v))
  stopifnot(nrow(W_u) == nrow(W_v))

  M_u <- obj_u$M
  M_v <- obj_v$M
  T_u <- obj_u$T
  T_v <- obj_v$T

  # 1. Khatri-Rao product (raw, before any orthogonalization)
  if (use_cpp && exists("khatri_rao_cpp")) {
    W_uv_raw <- khatri_rao_cpp(W_u, W_v)
  } else {
    W_uv_raw <- khatri_rao_rowwise_R(W_u, W_v)
  }

  # 2. Optional design-level ANOVA orthogonalization
  if (orthogonalize) {
    ortho <- anova_orthogonalize_W_uv(W_uv_raw, W_u, W_v)
    W_uv  <- ortho$W_uv_ortho
    A_uv  <- ortho$A_uv
  } else {
    W_uv  <- W_uv_raw
    A_uv  <- NULL
  }

  # 3. 2D RW2 penalty on identified beta_uv (unchanged by orthogonalization)
  K_u_raw <- build_rw2_penalty_1d(M_u)
  K_v_raw <- build_rw2_penalty_1d(M_v)
  pen <- build_rw2_2d_penalty(K_u_raw, K_v_raw, T_u, T_v)
  K_uv  <- pen$K_uv_beta
  T_uv  <- pen$T_uv

  # 4. Recipe for prediction at new data
  recipe <- list(
    M_u = M_u, M_v = M_v,
    T_u = T_u, T_v = T_v,
    AinvC_u = obj_u$AinvC,
    AinvC_v = obj_v$AinvC,
    tau_u = obj_u$tau,
    tau_v = obj_v$tau,
    design_new_u = obj_u$design_new,
    design_new_v = obj_v$design_new,
    # NEW: orthogonalization recipe
    orthogonalize = orthogonalize,
    A_uv          = A_uv          # NULL if !orthogonalize
  )

  list(
    W_uv          = W_uv,
    K_uv          = K_uv,
    T_uv          = T_uv,
    K_u_raw       = K_u_raw,
    K_v_raw       = K_v_raw,
    M_u = M_u, M_v = M_v,
    recipe        = recipe,
    orthogonalize = orthogonalize
  )
}


# ============================================================
# ls_interaction_design_new(X_u_new, X_v_new, recipe, clip = TRUE)
#
# Build the interaction design matrix at NEW data points, applying the SAME
# orthogonalization that was applied at training time (if any).
#
# Concretely: at training time we set
#   W_uv_train = W_uv_train_raw - M_train %*% A_uv
# where A_uv = (M_train' M_train)^-1 M_train' W_uv_train_raw.
#
# At prediction time, we set
#   W_uv_new = W_uv_new_raw - M_new %*% A_uv
# i.e. we subtract the SAME projection coefficients (frozen at training time)
# from the new Khatri-Rao product, evaluated at the new data points. This is
# the analog of "subtract the training-set mean from new points" for one-sample
# centering, generalized to a multi-column projection.
# ============================================================
ls_interaction_design_new <- function(X_u_new, X_v_new, recipe, clip = TRUE) {
  W_u_new <- recipe$design_new_u(X_u_new, type = "W", clip = clip)
  W_v_new <- recipe$design_new_v(X_v_new, type = "W", clip = clip)
  W_uv_raw <- khatri_rao_rowwise_R(W_u_new, W_v_new)

  if (isTRUE(recipe$orthogonalize) && !is.null(recipe$A_uv)) {
    M_new <- cbind(1, W_u_new, W_v_new)
    W_uv_raw - M_new %*% recipe$A_uv
  } else {
    W_uv_raw
  }
}


# ============================================================
# ls_build_all_interactions(full_objs, orthogonalize = FALSE)
#
# Build interaction designs for all p*(p-1)/2 covariate pairs.
# orthogonalize is forwarded to each call of ls_build_interaction().
# ============================================================
ls_build_all_interactions <- function(full_objs, orthogonalize = FALSE) {
  p <- length(full_objs)
  if (p < 2) stop("Need at least 2 covariates for interactions.")
  pairs <- combn(p, 2, simplify = FALSE)
  int_list <- vector("list", length(pairs))
  names_list <- character(length(pairs))

  for (k in seq_along(pairs)) {
    u <- pairs[[k]][1]
    v <- pairs[[k]][2]
    int_list[[k]] <- ls_build_interaction(full_objs[[u]], full_objs[[v]],
                                          orthogonalize = orthogonalize)
    names_list[k] <- paste0(u, "_", v)
  }
  names(int_list) <- names_list
  int_list
}


# ============================================================
# Assemble full design matrix H for the model with interactions
# Unchanged from original.
# ============================================================
ls_assemble_full_design <- function(W_main, int_list, n) {
  intercept <- matrix(1, n, 1)
  int_blocks <- lapply(int_list, function(x) x$W_uv)

  H <- cbind(intercept, W_main, do.call(cbind, int_blocks))

  p_main <- ncol(W_main)
  p_int_each <- sapply(int_blocks, ncol)

  col_map <- list()
  col_map$intercept <- 0L

  col_map$main <- list()
  start <- 1L
  for (j in seq_len(ncol(W_main) > 0)) {
    col_map$main[[j]] <- start
  }

  col_map$interaction <- vector("list", length(int_list))
  int_start <- 1L + p_main
  for (k in seq_along(int_list)) {
    col_map$interaction[[k]] <- int_start + seq_len(p_int_each[k]) - 1L
    int_start <- int_start + p_int_each[k]
  }

  list(H = H, col_map = col_map, p_main = p_main, p_int_each = p_int_each)
}


# ============================================================
# Sanity checks
# ============================================================
ls_interaction_tests <- function(n = 200, M = 8, seed = 42,
                                 orthogonalize = FALSE) {
  set.seed(seed)
  cat(sprintf("=== ls_interaction sanity checks (orthogonalize=%s) ===\n",
              orthogonalize))

  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)

  bu <- ls_build_one_full(X1, M = M)
  bv <- ls_build_one_full(X2, M = M)

  int <- ls_build_interaction(bu, bv, orthogonalize = orthogonalize)

  expected_cols <- (M - 1)^2
  stopifnot(nrow(int$W_uv) == n)
  stopifnot(ncol(int$W_uv) == expected_cols)
  cat(sprintf("  W_uv: %d x %d  OK\n", nrow(int$W_uv), ncol(int$W_uv)))

  stopifnot(nrow(int$K_uv) == expected_cols)
  cat(sprintf("  K_uv: %d x %d  OK\n", nrow(int$K_uv), ncol(int$K_uv)))

  ev <- eigen(int$K_uv, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(ev > -1e-8))
  cat(sprintf("  K_uv PSD: min eig = %.2e  OK\n", min(ev)))

  # NEW: design-level ANOVA orthogonality check
  W_u <- bu$W; W_v <- bv$W
  cat("\n  -- Design-level ANOVA orthogonality --\n")
  ip_one <- max(abs(colSums(int$W_uv)))
  ip_u   <- max(abs(crossprod(W_u, int$W_uv)))
  ip_v   <- max(abs(crossprod(W_v, int$W_uv)))
  cat(sprintf("    max |1' W_uv|         = %.3e\n", ip_one))
  cat(sprintf("    max |W_u' W_uv|       = %.3e\n", ip_u))
  cat(sprintf("    max |W_v' W_uv|       = %.3e\n", ip_v))
  if (orthogonalize) {
    if (max(ip_one, ip_u, ip_v) < 1e-8) {
      cat("    -> ANOVA-orthogonal at observations (as expected).\n")
    } else {
      warning("Orthogonalization did not produce zero inner products!")
    }
  } else {
    cat("    (without orthogonalization, these can be large -- known issue)\n")
  }

  # Prediction at new data
  X1_new <- runif(50, 0, 1)
  X2_new <- runif(50, 0, 1)
  W_uv_new <- ls_interaction_design_new(X1_new, X2_new, int$recipe)
  stopifnot(nrow(W_uv_new) == 50, ncol(W_uv_new) == expected_cols)
  cat(sprintf("\n  Prediction design at new X: %d x %d  OK\n",
              nrow(W_uv_new), ncol(W_uv_new)))

  cat("=== checks passed ===\n")
  invisible(int)
}
