# ============================================================
# ls_basis.R
# Lancaster-Salkauskas (LS) natural cubic spline basis.
#
# For one covariate x with knots tau_1 < ... < tau_M, the LS basis is
# cardinal: f(x) = sum_m theta_m B_m(x) with theta_m = f(tau_m).
# The design is Z = Phi + Psi %*% A^{-1} C (n x M), and the identified
# design used in the regression is W = Z %*% T (n x (M-1)), where T
# enforces sum(theta) = 0.
#
# Main entry points: ls_build_one_full(), ls_additive_build().
# ============================================================

# ------------------------------------------------------------
# ls_contrast_T(M): contrast matrix T (M x (M-1)) so that
# theta = T %*% beta satisfies sum(theta) = 0.
# ------------------------------------------------------------
ls_contrast_T <- function(M) {
  stopifnot(M >= 2)
  T <- matrix(0, nrow = M, ncol = M - 1)
  T[1, ] <- -1
  for (m in 2:M) T[m, m - 1] <- 1
  T
}

# ------------------------------------------------------------
# ls_choose_knots_equal(x, M): M equally spaced knots on [min(x), max(x)].
# ------------------------------------------------------------
ls_choose_knots_equal <- function(x, M) {
  x <- as.numeric(x)
  if (M < 4) stop("Need M >= 4 knots for cubic spline.")
  xmin <- min(x); xmax <- max(x)
  if (!is.finite(xmin) || !is.finite(xmax)) stop("x contains non-finite values.")
  if (xmin == xmax) stop("x has zero range.")
  seq(xmin, xmax, length.out = as.integer(M))
}

# ------------------------------------------------------------
# ls_AC(tau): natural-spline constraint matrices A and C (M x M) for the
# knot vector tau, and the operator A^{-1} C that maps knot ordinates to
# knot slopes. Returns list(A, C, AinvC, tau).
# ------------------------------------------------------------
ls_AC <- function(tau) {
  tau <- sort(unique(as.numeric(tau)))
  M <- length(tau)
  if (M < 4) stop("Need at least 4 knots.")
  if (any(diff(tau) <= 0)) stop("Knots must be strictly increasing.")
  
  omega     <- rep(NA_real_, M)
  omega_bar <- rep(NA_real_, M)
  
  for (m in 2:(M - 1)) {
    h_m   <- tau[m]   - tau[m - 1]
    h_mp1 <- tau[m+1] - tau[m]
    omega[m]     <- h_m / (h_m + h_mp1)
    omega_bar[m] <- 1 - omega[m]
  }
  
  # A
  A <- matrix(0, M, M)
  diag(A) <- 2
  A[1, 2] <- 1
  A[M, M-1] <- 1
  for (m in 2:(M - 1)) {
    A[m, m-1] <- omega[m]
    A[m, m+1] <- omega_bar[m]
  }
  
  # C (then *3)
  C <- matrix(0, M, M)
  h2 <- tau[2] - tau[1]
  hM <- tau[M] - tau[M-1]
  C[1, 1]   <- -1 / h2
  C[1, 2]   <-  1 / h2
  C[M, M-1] <- -1 / hM
  C[M, M]   <-  1 / hM
  
  for (m in 2:(M - 1)) {
    h_m   <- tau[m]   - tau[m - 1]   # h_m
    h_mp1 <- tau[m+1] - tau[m]       # h_{m+1}
    
    C[m, m-1] <- - omega[m]     / h_m
    C[m, m+1] <-   omega_bar[m] / h_mp1
    C[m, m]   <-   omega[m]     / h_m - omega_bar[m] / h_mp1
  }
  C <- 3 * C
  
  AinvC <- solve(A, C)
  
  list(A = A, C = C, AinvC = AinvC, tau = tau)
}

# ------------------------------------------------------------
# ls_phi_psi(x, tau): evaluate the raw Hermite-type basis functions
# Phi_m(x) (value) and Psi_m(x) (slope) at each x. Both are n x M.
# ------------------------------------------------------------
ls_phi_psi <- function(x, tau) {
  # Note: we use "<=" only on the LAST interval so x==tau_M is handled.
  x <- as.numeric(x)
  tau <- sort(unique(as.numeric(tau)))
  M <- length(tau)
  n <- length(x)
  
  Phi <- matrix(0, n, M)
  Psi <- matrix(0, n, M)
  
  # m = 1
  h2 <- tau[2] - tau[1]
  idx <- (x >= tau[1]) & (x < tau[2])
  Phi[idx, 1] <- (2 / h2^3) * (x[idx] - tau[2])^2 * (x[idx] - tau[1] + 0.5 * h2)
  Psi[idx, 1] <- (1 / h2^2) * (x[idx] - tau[2])^2 * (x[idx] - tau[1])
  
  # m = M
  hM <- tau[M] - tau[M-1]
  idx <- (x >= tau[M-1]) & (x <= tau[M])
  Phi[idx, M] <- -(2 / hM^3) * (x[idx] - tau[M-1])^2 * (x[idx] - tau[M] - 0.5 * hM)
  Psi[idx, M] <-  (1 / hM^2) * (x[idx] - tau[M-1])^2 * (x[idx] - tau[M])
  
  # interior m = 2..M-1
  for (m in 2:(M - 1)) {
    h_m   <- tau[m]   - tau[m - 1]
    h_mp1 <- tau[m+1] - tau[m]
    
    idx1 <- (x >= tau[m-1]) & (x <  tau[m])
    idx2 <- (x >= tau[m])   & (x <  tau[m+1])
    
    Phi[idx1, m] <- -(2 / h_m^3)   * (x[idx1] - tau[m-1])^2 * (x[idx1] - tau[m] - 0.5 * h_m)
    Phi[idx2, m] <-  (2 / h_mp1^3) * (x[idx2] - tau[m+1])^2 * (x[idx2] - tau[m] + 0.5 * h_mp1)
    
    Psi[idx1, m] <- (1 / h_m^2)    * (x[idx1] - tau[m-1])^2 * (x[idx1] - tau[m])
    Psi[idx2, m] <- (1 / h_mp1^2)  * (x[idx2] - tau[m+1])^2 * (x[idx2] - tau[m])
  }
  
  list(Phi = Phi, Psi = Psi)
}

# ------------------------------------------------------------
# ls_build_one_full(x, M = NULL, tau = NULL)
# ------------------------------------------------------------
# PURPOSE:
#   Build the COMPLETE LS-basis machinery for ONE covariate vector x.
#   This is the "debug/inspect everything" constructor.
#
# WHAT IT DOES (conceptually):
#   1) Choose knots tau on [min(x), max(x)] if tau not provided (needs M).
#   2) Build A and C matrices (natural spline constraints) and compute A^{-1} C.
#   3) Evaluate LS basis functions Phi(x), Psi(x) at each observation.
#   4) Form the LS design:
#        Z = Phi + Psi %*% (A^{-1} C)        # n x M
#      so that f(x) = Z %*% theta, where theta are knot ordinates.
#   5) Build the identifiability contrast:
#        T  enforces sum(theta)=0 via theta = T %*% beta
#        W  = Z %*% T                        # n x (M-1)
#      so f(x) = W %*% beta with identified parameters.
#   6) Returns a closure design_new() so you can build Z/W for new x later
#      using the SAME tau and SAME A^{-1}C and SAME T.
#
# RETURNS (key items):
#   - tau      : knot locations used
#   - AinvC    : A^{-1}C operator (encodes natural-spline constraints)
#   - Phi, Psi : raw basis matrices (n x M) (useful for sanity checks)
#   - Z        : full LS design (n x M), NOT identified
#   - T        : contrast matrix (M x (M-1)) enforcing sum(theta)=0
#   - W        : identified LS design (n x (M-1)) = Z %*% T
#   - design_new(x_new, type=c("Z","W"), clip=FALSE):
#       builds matching design for new x (prediction/test time)
#
# WHEN TO USE:
#   Use this when developing / verifying correctness / inspecting Phi/Psi/Z/W.

ls_build_one_full <- function(x, M = NULL, tau = NULL) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("x contains non-finite values.")
  
  if (is.null(tau)) {
    if (is.null(M)) stop("Provide M when tau is NULL.")
    tau <- ls_choose_knots_equal(x, M)
  } else {
    tau <- sort(as.numeric(tau))
    if (length(tau) < 4) stop("Need at least 4 knots.")
    if (any(diff(tau) <= 0)) stop("tau must be strictly increasing.")
    M <- length(tau)
  }
  
  ac <- ls_AC(tau)
  bp <- ls_phi_psi(x, ac$tau)
  Phi <- bp$Phi
  Psi <- bp$Psi
  
  Z <- Phi + Psi %*% ac$AinvC
  colnames(Z) <- paste0("B", 1:ncol(Z))
  
  T <- ls_contrast_T(M)
  W <- Z %*% T
  colnames(W) <- paste0("W", 1:ncol(W))
  
  design_new <- function(x_new, type = c("Z", "W"), clip = FALSE) {
    type <- match.arg(type)
    x_new <- as.numeric(x_new)
    if (clip) x_new <- pmin(pmax(x_new, ac$tau[1]), ac$tau[length(ac$tau)])
    
    bp_new <- ls_phi_psi(x_new, ac$tau)
    Z_new <- bp_new$Phi + bp_new$Psi %*% ac$AinvC
    if (type == "Z") return(Z_new)
    Z_new %*% T
  }
  
  list(M = M, tau = ac$tau, AinvC = ac$AinvC, T = T,
       Phi = Phi, Psi = Psi, Z = Z, W = W,
       design_new = design_new)
}

# ------------------------------------------------------------
# ls_build_one_train(x, M, tau = NULL)
# ------------------------------------------------------------
# PURPOSE:
#   Build the LS design for ONE covariate for TRAINING, but return a
#   LIGHTWEIGHT object that stores only what you need for later prediction.
#
# WHAT IT DOES:
#   - Calls ls_build_one_full() internally
#   - Keeps training matrices (W and Z), and keeps a compact "obj"
#     (the reusable recipe) for building designs on new x.
#
# RETURNS:
#   - W   : identified design (n x (M-1))  [typically used in the regression]
#   - Z   : full design       (n x M)      [sometimes useful for debugging]
#   - obj : list containing:
#       * M, tau, AinvC, T
#       * design_new() closure (build Z/W for new x using same settings)
#
# WHEN TO USE:
#   This is the standard "train-time" constructor for one covariate.
ls_build_one_train <- function(x, M, tau = NULL) {
  full <- ls_build_one_full(x, M = M, tau = tau)
  obj <- list(M = full$M, tau = full$tau, AinvC = full$AinvC, T = full$T,
              design_new = full$design_new)
  list(W = full$W, Z = full$Z, obj = obj)
}

# ------------------------------------------------------------
# ls_additive_build(X, M_vec, tau_list = NULL)
# ------------------------------------------------------------
# PURPOSE:
#   Build LS-basis designs for MULTIPLE covariates (additive model),
#   and concatenate them into one block design matrix.
#
# INPUTS:
#   - X      : n x p matrix of covariates (each column is one covariate)
#   - M_vec  : either length-1 (same M for all) or length-p (per covariate)
#   - tau_list (optional): list of length p, where tau_list[[j]] gives
#     knots for covariate j (useful to fix knots across replications)
#
# WHAT IT DOES:
#   For each covariate column j:
#     - calls ls_build_one_train(X[,j], M=M_vec[j], tau=tau_list[[j]])
#     - gets W_j (n x (M_j-1)) and an obj_j (recipe for new data)
#   Then:
#     - W_block = cbind(W_1, W_2, ..., W_p)
#     - records where each block sits via col_map (for interpretation)
#
# RETURNS:
#   - W       : block design matrix (n x sum_j(M_j-1)) for fitting
#   - W_list  : list of each covariate's W_j
#   - objs    : list of per-covariate "recipes" (each has tau, AinvC, T, design_new)
#   - col_map : list mapping covariate j -> column indices in the big W
#
# WHEN TO USE:
#   Use this right before fitting the additive model:
#     y ≈ intercept/linear terms + W %*% beta + (spatial/random/noise part)
ls_additive_build <- function(X, M_vec, tau_list = NULL) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  
  if (length(M_vec) == 1) M_vec <- rep(M_vec, p)
  stopifnot(length(M_vec) == p)
  
  objs <- vector("list", p)
  W_list <- vector("list", p)
  col_map <- vector("list", p)
  col_start <- 1
  
  for (j in 1:p) {
    xj <- X[, j]
    tau_j <- if (!is.null(tau_list)) tau_list[[j]] else NULL
    
    tmp <- ls_build_one_train(xj, M = M_vec[j], tau = tau_j)
    objs[[j]]   <- tmp$obj
    W_list[[j]] <- tmp$W
    
    kj <- ncol(tmp$W)
    col_map[[j]] <- col_start:(col_start + kj - 1)
    col_start <- col_start + kj
  }
  
  W_block <- do.call(cbind, W_list)
  colnames(W_block) <- paste0("W", seq_len(ncol(W_block)))
  list(W = W_block, W_list = W_list, objs = objs, col_map = col_map)
}

# ------------------------------------------------------------
# ls_additive_design_new(X_new, objs, clip = TRUE)
# ------------------------------------------------------------
# PURPOSE:
#   Build the matching additive-model LS design matrix for NEW covariate data.
#   This is the "prediction/test-time" constructor.
#
# INPUTS:
#   - X_new : n_new x p matrix of new covariates
#   - objs  : the list returned by ls_additive_build() (per-covariate recipes)
#   - clip  : if TRUE, clamp new x values into [tau_1, tau_M] before building
#             designs (avoids out-of-range behavior with compact support)
#
# WHAT IT DOES:
#   For each covariate j:
#     - uses objs[[j]]$design_new(X_new[,j], type="W", clip=clip)
#   Then:
#     - W_new = cbind(W_new_1, ..., W_new_p)
#
# RETURNS:
#   - W_new : n_new x sum_j(M_j-1) block design, consistent with training W
#
# WHEN TO USE:
#   Use after fitting, when you want predictions/evaluations on test/new data.

ls_additive_design_new <- function(X_new, objs, clip = TRUE) {
  X_new <- as.matrix(X_new)
  p <- ncol(X_new)
  stopifnot(length(objs) == p)
  
  W_list_new <- vector("list", p)
  for (j in 1:p) {
    W_list_new[[j]] <- objs[[j]]$design_new(X_new[, j], type = "W", clip = clip)
  }
  W_new <- do.call(cbind, W_list_new)
  colnames(W_new) <- paste0("W", seq_len(ncol(W_new)))
  W_new
}

# quick LS sanity tests
ls_tests <- function() {
  M <- 6
  tau <- seq(0, 1, length.out = M)
  x_knots <- tau
  bp <- ls_phi_psi(x_knots, tau)
  stopifnot(max(abs(bp$Phi - diag(M))) < 1e-12)
  stopifnot(max(abs(bp$Psi)) < 1e-12)
  
  set.seed(1)
  x <- runif(2000, 0, 1)
  bp2 <- ls_phi_psi(x, tau)
  stopifnot(max(abs(rowSums(bp2$Phi) - 1)) < 1e-12)
  
  T <- ls_contrast_T(M)
  beta <- rnorm(M-1)
  theta <- as.vector(T %*% beta)
  stopifnot(abs(sum(theta)) < 1e-12)
  
  out <- ls_build_one_full(x, tau = tau)
  stopifnot(max(abs(out$Z %*% theta - out$W %*% beta)) < 1e-10)
  TRUE
}

