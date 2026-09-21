# ============================================================
# gibbs_interaction.R
# Collapsed Gibbs sampler for the Bayesian spatial additive model
# with an LS tensor-product interaction surface (p = 2 covariates)
#
# Model:
#   y(s_r) = mu + f_1(X_1r) + f_2(X_2r) + f_{12}(X_1r, X_2r) + b(s_r) + eps_r
#
#   f_j   : LS natural cubic spline with RW2 prior
#   f_{12}: LS tensor-product interaction, W_uv * beta_uv, RW2 prior
#   b     : Matern GP spatial effect, marginalized out (collapsed Gibbs)
#   eps   : iid N(0, tau2)
#
# Priors:
#   mu          ~ N(0, kappa2)
#   beta_j      ~ N(0, tau2_s_j * K_j^-)     j = 1, 2
#   beta_12     ~ N(0, tau2_s_12 * K_12^-)   (interaction, RW2 2D)
#   tau2_s_j    ~ IG(a_s, b_s)   (conjugate)
#   tau2_s_12   ~ IG(a_s, b_s)   (conjugate, separate smoothing variance)
#   sigma2      ~ IG(a_sigma, b_sigma)  via MH
#   tau2        ~ IG(a_tau, b_tau)      via MH
#   rho         ~ logN(log_rho_mu, log_rho_sd^2)  via MH
#
# Gibbs steps:
#   1) eta = (mu, beta_1, beta_2, beta_12) | y, sigma2, tau2, rho, {tau2_s}  ~ MVN
#   2) sigma2 | rest   MH
#   3) tau2   | rest   MH
#   4) rho    | rest   MH
#   5) tau2_s_j | beta_j  ~ IG  (for j=1,2 and j=12)
#   6) b_hat (derived posterior mean, stored)
#
# Depends on:
#   ls_basis.R          (ls_build_one_full)
#   ls_interaction.R    (ls_build_interaction, build_rw2_penalty_1d)
#   spatial_utils.R     (matern_cor)
#   ls_interaction_core.cpp  (optional; see toy/run_toy_cpp.R)
#
# Usage: see toy/run_toy.R for a complete worked example.
# ============================================================

# Report whether the compiled C++ helpers are present (optional; pure-R
# fallbacks are used otherwise).
.has_cpp <- tryCatch({
  requireNamespace("Rcpp", quietly = TRUE) &&
  exists("khatri_rao_cpp")
}, error = function(e) FALSE)

if (.has_cpp) message("Using C++ hot loops.") else
  message("Using pure R fallbacks (compile ls_interaction_core.cpp for speedup).")


# ============================================================
# build_interaction_prior_precision()
#
# Builds the FULL prior precision Q0 for eta = (mu, beta_1, beta_2, beta_12).
#
# col_map_full: list with
#   $main      : list of 0-indexed column ranges for main effect blocks
#   $interaction: list of 0-indexed column ranges for interaction blocks
# K_int_list: list of identified 2D RW2 penalty matrices (one per interaction pair)
# tau2_s: named/ordered vector (tau2_s_1, tau2_s_2, tau2_s_12)
# ============================================================
build_interaction_prior_precision <- function(col_map_full, K_main_list,
                                              K_int_list, tau2_s_main,
                                              tau2_s_int, kappa2 = 1e6,
                                              eps_ridge = 1e-6) {
  n_main <- if (length(K_main_list) > 0) sum(vapply(K_main_list, nrow, integer(1))) else 0L
  n_int  <- if (length(K_int_list)  > 0) sum(vapply(K_int_list,  nrow, integer(1))) else 0L
  p_total <- 1L + n_main + n_int
  Q0 <- matrix(0, p_total, p_total)

  # Intercept
  Q0[1, 1] <- 1 / kappa2

  # Main effect blocks
  for (j in seq_along(K_main_list)) {
    idx <- 1 + col_map_full$main[[j]]   # 0-indexed -> 1-indexed
    Kj  <- K_main_list[[j]]
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s_main[j]) * Kj + eps_ridge * diag(length(idx))
  }

  # Interaction blocks
  for (k in seq_along(K_int_list)) {
    idx <- 1 + col_map_full$interaction[[k]]
    Kk  <- K_int_list[[k]]
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s_int[k]) * Kk + eps_ridge * diag(length(idx))
  }

  Q0
}


# ============================================================
# gibbs_interaction_sampler()
#
# Main sampler.
#
# Inputs:
#   y            : response vector (length n)
#   H            : full design matrix [1 | W_1 | W_2 | W_12]  (n x p)
#   D            : n x n distance matrix between locations
#   nu           : Matern smoothness (fixed)
#   col_map_main : list of 0-indexed column ranges for each main-effect block in H
#   K_main_list  : list of identified 1D RW2 penalty matrices, one per main effect
#   col_map_int  : list of 0-indexed column ranges for each interaction block in H
#   K_int_list   : list of identified 2D RW2 penalty matrices, one per interaction
#
# Optional clamps (for sensitivity checks): any of tau2_s_main_fixed,
# tau2_s_int_fixed, sigma2_fixed, tau2_fixed, rho_fixed may be set to a
# finite value, in which case that component is held fixed and its update
# step is skipped. NA / NULL (the defaults) leave the sampler untouched.
# ============================================================
gibbs_interaction_sampler <- function(
    y, H, D, nu = 1.5,
    col_map_main,    # list: 0-indexed col ranges for main effects in H (after intercept)
    K_main_list,     # list: identified 1D RW2 penalty for each main effect
    col_map_int,     # list: 0-indexed col ranges for each interaction block in H
    K_int_list,      # list: identified 2D RW2 penalty for each interaction
    n_iter  = 5000,
    n_burn  = 1000,
    n_thin  = 1,
    kappa2  = 1e6,
    a_sigma = 2, b_sigma = 1,
    a_tau   = 2, b_tau   = 0.3,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = -1.6, log_rho_sd = 1.0,
    mh_sd_log_sigma2 = 0.3,
    mh_sd_log_tau2   = 0.3,
    mh_sd_log_rho    = 0.2,
    eps_ridge = 1e-6,
    init      = NULL,
    jitter    = 1e-8,
    verbose   = TRUE,
    # Optional clamps on the main-effect smoothing variances: numeric vector
    # of length n_main. A finite entry j holds tau2_s_main[j] fixed at that
    # value and skips its IG draw; NA entries (or NULL) get the usual draw.
    tau2_s_main_fixed = NULL,
    # Same for the interaction smoothing variances (length n_int).
    tau2_s_int_fixed = NULL,
    # Optional scalar clamps on the variance components. When finite, the
    # MH step for that component is skipped and the value held fixed.
    sigma2_fixed = NA_real_,
    tau2_fixed   = NA_real_,
    rho_fixed    = NA_real_,
    # Additive-only switch. When FALSE, Step 5b (interaction smoothing
    # variance draw) is skipped. The caller must then pass empty
    # col_map_int / K_int_list and no interaction columns in H.
    fit_interactions = TRUE
) {
  y <- as.numeric(y)
  H <- as.matrix(H)
  D <- as.matrix(D)
  n <- length(y)
  p <- ncol(H)

  n_main <- length(col_map_main)
  n_int  <- length(col_map_int)

  # Internal consistency between the additive-only switch and the
  # caller-supplied interaction structures.
  if (!fit_interactions && n_int > 0L) {
    stop("fit_interactions = FALSE but col_map_int has length ",
         n_int, "; caller must pass empty col_map_int / K_int_list ",
         "(and not include interaction columns in H).")
  }

  # Pre-build col_map_full for Q0 constructor
  col_map_full <- list(main = col_map_main, interaction = col_map_int)

  # Helper: Matern correlation
  compute_R <- function(rho) matern_cor(D, rho = rho, nu = nu)

  # Log marginal likelihood (quadratic part only; log-det separated)
  log_marg_lik <- function(resid, sigma2, tau2, R) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- tryCatch(chol(Sigma + diag(jitter, n)), error = function(e) NULL)
    if (is.null(L)) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    alpha  <- forwardsolve(t(L), resid)
    -0.5 * (logdet + sum(alpha^2))
  }

  log_ig_prior <- function(x, a, b) -(a + 1) * log(x) - b / x

  log_lognormal_prior <- function(rho, mu, sd)
    dnorm(log(rho), mean = mu, sd = sd, log = TRUE) - log(rho)

  compute_b_postmean <- function(resid, sigma2, tau2, R) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    Sigma_inv <- chol2inv(chol(Sigma + diag(jitter, n)))
    as.vector(sigma2 * R %*% Sigma_inv %*% resid)
  }

  # --- Initialization ---
  if (is.null(init)) {
    sigma2      <- 1.0
    tau2        <- 0.5
    rho         <- 0.2
    eta         <- rep(0, p)
    tau2_s_main <- rep(1.0, n_main)
    tau2_s_int  <- rep(1.0, n_int)
  } else {
    sigma2      <- init$sigma2
    tau2        <- init$tau2
    rho         <- init$rho
    eta         <- init$eta
    tau2_s_main <- init$tau2_s_main
    tau2_s_int  <- init$tau2_s_int
  }
  # Apply scalar clamps if requested (override init / default)
  if (is.finite(sigma2_fixed)) sigma2 <- sigma2_fixed
  if (is.finite(tau2_fixed))   tau2   <- tau2_fixed
  if (is.finite(rho_fixed))    rho    <- rho_fixed
  R <- compute_R(rho)

  # --- Storage ---
  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples      <- matrix(NA, n_keep, p)
  b_samples        <- matrix(NA, n_keep, n)
  sigma2_samples   <- numeric(n_keep)
  tau2_samples     <- numeric(n_keep)
  rho_samples      <- numeric(n_keep)
  tau2_s_main_samp <- matrix(NA, n_keep, n_main)
  tau2_s_int_samp  <- matrix(NA, n_keep, n_int)

  accept <- c(sigma2 = 0, tau2 = 0, rho = 0)
  keep_idx <- 0

  # ============================================================
  # MAIN LOOP
  # ============================================================
  for (iter in seq_len(n_iter)) {

    # ----------------------------------------------------------
    # Step 1: Draw eta | y, sigma2, tau2, rho, {tau2_s}
    #
    # eta = (mu, beta_1, beta_2, beta_12) ~ MVN
    # Q_eta = H^T Sigma^{-1} H + Q0
    # m_eta = Q_eta^{-1} H^T Sigma^{-1} y
    # ----------------------------------------------------------
    Q0 <- build_interaction_prior_precision(
      col_map_full, K_main_list, K_int_list,
      tau2_s_main, tau2_s_int, kappa2, eps_ridge
    )

    Sigma    <- sigma2 * R + tau2 * diag(n)
    L_Sigma  <- chol(Sigma + diag(jitter, n))
    y_w      <- forwardsolve(t(L_Sigma), y)
    H_w      <- forwardsolve(t(L_Sigma), H)

    Q_eta <- crossprod(H_w) + Q0
    U_eta <- chol(Q_eta)
    V_eta <- chol2inv(U_eta)
    m_eta <- as.vector(V_eta %*% crossprod(H_w, y_w))

    z   <- rnorm(p)
    eta <- m_eta + as.vector(backsolve(U_eta, z))

    resid <- as.vector(y - H %*% eta)

    # ----------------------------------------------------------
    # Steps 2-4: MH for sigma2, tau2, rho
    # ----------------------------------------------------------

    # sigma2
    if (!is.finite(sigma2_fixed)) {
      sigma2_prop <- exp(log(sigma2) + rnorm(1, 0, mh_sd_log_sigma2))
      lp_c <- log_marg_lik(resid, sigma2,      tau2, R) + log_ig_prior(sigma2,      a_sigma, b_sigma)
      lp_p <- log_marg_lik(resid, sigma2_prop, tau2, R) + log_ig_prior(sigma2_prop, a_sigma, b_sigma)
      if (log(runif(1)) < (lp_p - lp_c)) { sigma2 <- sigma2_prop; accept["sigma2"] <- accept["sigma2"] + 1 }
    }

    # tau2
    if (!is.finite(tau2_fixed)) {
      tau2_prop <- exp(log(tau2) + rnorm(1, 0, mh_sd_log_tau2))
      lp_c <- log_marg_lik(resid, sigma2, tau2,      R) + log_ig_prior(tau2,      a_tau, b_tau)
      lp_p <- log_marg_lik(resid, sigma2, tau2_prop, R) + log_ig_prior(tau2_prop, a_tau, b_tau)
      if (log(runif(1)) < (lp_p - lp_c)) { tau2 <- tau2_prop; accept["tau2"] <- accept["tau2"] + 1 }
    }

    # rho
    if (!is.finite(rho_fixed)) {
      rho_prop <- exp(log(rho) + rnorm(1, 0, mh_sd_log_rho))
      R_prop   <- compute_R(rho_prop)
      lp_c <- log_marg_lik(resid, sigma2, tau2, R)      + log_lognormal_prior(rho,      log_rho_mu, log_rho_sd)
      lp_p <- log_marg_lik(resid, sigma2, tau2, R_prop) + log_lognormal_prior(rho_prop, log_rho_mu, log_rho_sd)
      if (log(runif(1)) < (lp_p - lp_c)) { rho <- rho_prop; R <- R_prop; accept["rho"] <- accept["rho"] + 1 }
    }

    # ----------------------------------------------------------
    # Step 5a: Draw tau2_s_j | beta_j  (main effects, conjugate IG)
    # ----------------------------------------------------------
    for (j in seq_len(n_main)) {
      # Skip the IG draw if this component is clamped
      if (!is.null(tau2_s_main_fixed) &&
          length(tau2_s_main_fixed) >= j &&
          is.finite(tau2_s_main_fixed[j])) {
        tau2_s_main[j] <- tau2_s_main_fixed[j]
        next
      }
      idx    <- 1 + col_map_main[[j]]
      beta_j <- eta[idx]
      K_j    <- K_main_list[[j]]
      qf     <- as.numeric(t(beta_j) %*% K_j %*% beta_j)
      rank_j <- nrow(K_j) - 2
      if (rank_j < 1) rank_j <- 1
      tau2_s_main[j] <- 1 / rgamma(1, shape = a_smooth + rank_j / 2,
                                       rate  = b_smooth + qf / 2)
    }

    # ----------------------------------------------------------
    # Step 5b: Draw tau2_s_12 | beta_12  (interaction, conjugate IG)
    # Skipped when fit_interactions = FALSE.
    # ----------------------------------------------------------
    if (fit_interactions) {
      for (k in seq_len(n_int)) {
        # Skip the IG draw if this component is clamped
        if (!is.null(tau2_s_int_fixed) &&
            length(tau2_s_int_fixed) >= k &&
            is.finite(tau2_s_int_fixed[k])) {
          tau2_s_int[k] <- tau2_s_int_fixed[k]
          next
        }
        idx    <- 1 + col_map_int[[k]]
        beta_k <- eta[idx]
        K_k    <- K_int_list[[k]]
        qf     <- as.numeric(t(beta_k) %*% K_k %*% beta_k)
        # Rank of 2D RW2 penalty: (M_u-1)*(M_v-1) - 1
        rank_k <- nrow(K_k) - 1
        if (rank_k < 1) rank_k <- 1
        tau2_s_int[k] <- 1 / rgamma(1, shape = a_smooth + rank_k / 2,
                                        rate  = b_smooth + qf / 2)
      }
    }

    # ----------------------------------------------------------
    # Step 6: Posterior mean of b (derived)
    # ----------------------------------------------------------
    b_mean <- compute_b_postmean(resid, sigma2, tau2, R)

    # ----------------------------------------------------------
    # Store
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]          <- eta
      b_samples[keep_idx, ]            <- b_mean
      sigma2_samples[keep_idx]         <- sigma2
      tau2_samples[keep_idx]           <- tau2
      rho_samples[keep_idx]            <- rho
      tau2_s_main_samp[keep_idx, ]     <- tau2_s_main
      tau2_s_int_samp[keep_idx, ]      <- tau2_s_int
    }

    if (verbose && (iter %% 500 == 0)) {
      cat(sprintf("  iter %d/%d  rho=%.3f  sigma2=%.3f  tau2=%.3f  tau2_s_int=%.4f\n",
                  iter, n_iter, rho, sigma2, tau2, tau2_s_int[1]))
    }
  }

  cat(sprintf("  Accept rates: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
              accept["sigma2"] / n_iter, accept["tau2"] / n_iter, accept["rho"] / n_iter))

  list(
    eta_samples      = eta_samples,
    b_samples        = b_samples,
    sigma2_samples   = sigma2_samples,
    tau2_samples     = tau2_samples,
    rho_samples      = rho_samples,
    tau2_s_main_samp = tau2_s_main_samp,
    tau2_s_int_samp  = tau2_s_int_samp,
    n_burn = n_burn, n_thin = n_thin, n_iter = n_iter,
    accept_rate = accept / n_iter
  )
}
