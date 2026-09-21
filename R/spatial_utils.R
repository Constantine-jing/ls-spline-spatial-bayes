# ============================================================
# spatial_utils.R
# Distance matrices and the Matern correlation function for the
# spatial Gaussian process. Correlation is set to 1 at distance 0.
# ============================================================

# pairdist(coords): n x n Euclidean distance matrix from an n x 2 coordinate matrix
pairdist <- function(coords) {
  as.matrix(dist(as.matrix(coords)))
}

# pairdist_cross(A, B): distances between rows of A (n_A x 2) and rows of B (n_B x 2)
pairdist_cross <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  dx <- outer(A[,1], B[,1], "-")
  dy <- outer(A[,2], B[,2], "-")
  sqrt(dx^2 + dy^2)
}

# matern_cor(D, rho, nu): Matern correlation matrix with range rho and
# smoothness nu, evaluated at the distance matrix D
matern_cor <- function(D, rho = 0.25, nu = 1.5) {
  D <- as.matrix(D)
  if (rho <= 0) stop("rho must be > 0")
  if (nu  <= 0) stop("nu must be > 0")
  
  R <- matrix(0, nrow(D), ncol(D))
  diag(R) <- 1
  
  idx <- (D > 0)
  if (any(idx)) {
    Z <- D[idx] / rho
    R[idx] <- (2^(1 - nu) / gamma(nu)) * (Z^nu) * besselK(Z, nu)
  }
  R
}

# matern_cor_cross(D, rho, nu): same for a rectangular cross-distance matrix
matern_cor_cross <- function(D, rho = 0.25, nu = 1.5) {
  D <- as.matrix(D)
  if (rho <= 0) stop("rho must be > 0")
  if (nu  <= 0) stop("nu must be > 0")
  
  out <- matrix(0, nrow(D), ncol(D))
  idx <- (D > 0)
  
  if (any(idx)) {
    Z <- D[idx] / rho
    out[idx] <- (2^(1 - nu) / gamma(nu)) * (Z^nu) * besselK(Z, nu)
  }
  
  # distance==0 => correlation=1
  out[!idx] <- 1
  out
}

# solve_chol(U, b): solve (U'U) x = b given the upper Cholesky factor U
solve_chol <- function(U, b) {
  forwardsolve(t(U), backsolve(U, b))
}
