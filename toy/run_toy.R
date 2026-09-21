# ============================================================
# toy/run_toy.R
#
# Minimal worked example on simulated data. Pure R; no compiler needed.
#
# Model fitted:
#   y(s) = mu + f1(x1) + f2(x2) + f12(x1, x2) + b(s) + eps
#   f1, f2  : LS natural cubic splines with RW2 priors
#   f12     : LS tensor-product interaction, orthogonalized against
#             [1, W1, W2] at the observations (design-level ANOVA)
#   b(s)    : Matern Gaussian process, integrated out (collapsed Gibbs)
#   eps     : iid N(0, tau2)
#
# Steps:
#   1. simulate n = 200 spatial observations with known f1, f2, f12
#   2. build the LS design matrices W1, W2 and the interaction design W12
#   3. stack H = [1 | W1 | W2 | W12] and record which columns are which
#   4. run the collapsed Gibbs sampler
#   5. recover f1, f2, f12 from the posterior mean and compare to the truth
#   6. write trace plots and fitted-vs-true plots to toy/toy_plots.pdf
#
# Run from the repository root:
#   source("toy/run_toy.R")          # inside R / RStudio
#   Rscript toy/run_toy.R            # from a terminal
# Takes about a minute on a laptop.
# ============================================================

# use_cpp is set to TRUE by toy/run_toy_cpp.R; otherwise pure R
if (!exists("use_cpp")) use_cpp <- FALSE

root <- if (file.exists("R/ls_basis.R")) "." else ".."
source(file.path(root, "R", "ls_basis.R"))
source(file.path(root, "R", "spatial_utils.R"))
source(file.path(root, "R", "ls_interaction.R"))
source(file.path(root, "R", "gibbs_interaction.R"))

# ---- Settings (the paper's simulations use n = 1000, M = 20, n_iter >= 5000) ----
n      <- 200      # observations
M      <- 8        # knots per covariate
n_iter <- 1500     # MCMC iterations
n_burn <- 500      # burn-in
seed   <- 1

# ---- 1. Simulate data ----
set.seed(seed)
locs <- matrix(runif(2 * n), n, 2)          # locations on the unit square
D    <- as.matrix(dist(locs))               # pairwise distances
X1   <- runif(n)
X2   <- runif(n)

# True functions. Each has mean zero, and f12 has zero marginal means,
# so the three components are separately identified (ANOVA decomposition).
true_f1  <- sin(2 * pi * X1)
true_f2  <- cos(2 * pi * X2)
true_f12 <- 1.5 * sin(2 * pi * X1) * sin(2 * pi * X2)

true_mu     <- 1.0
true_sigma2 <- 0.5      # GP partial sill
true_rho    <- 0.3      # GP range
true_nu     <- 1.5      # Matern smoothness (fixed in the fit)
true_tau2   <- 0.25     # nugget / noise variance

R_true <- matern_cor(D, rho = true_rho, nu = true_nu)
b_true <- as.vector(t(chol(R_true + diag(1e-8, n))) %*% rnorm(n)) * sqrt(true_sigma2)
eps    <- rnorm(n, 0, sqrt(true_tau2))
y      <- true_mu + true_f1 + true_f2 + true_f12 + b_true + eps

cat(sprintf("Simulated n = %d, M = %d.  y: mean = %.2f, sd = %.2f\n", n, M, mean(y), sd(y)))

# ---- 2. Design matrices and penalties ----
obj1 <- ls_build_one_full(X1, M = M)        # LS basis for x1
obj2 <- ls_build_one_full(X2, M = M)        # LS basis for x2
W1 <- obj1$W                                # n x (M-1), identified
W2 <- obj2$W

# 1D RW2 penalties on the identified coefficients
K1 <- t(obj1$T) %*% build_rw2_penalty_1d(M) %*% obj1$T
K2 <- t(obj2$T) %*% build_rw2_penalty_1d(M) %*% obj2$T

# Tensor-product interaction design and its 2D RW2 penalty
int12 <- ls_build_interaction(obj1, obj2, use_cpp = use_cpp, orthogonalize = TRUE)
W12   <- int12$W_uv                         # n x (M-1)^2
K12   <- int12$K_uv

cat(sprintf("Designs: W1 %dx%d, W2 %dx%d, W12 %dx%d\n",
            nrow(W1), ncol(W1), nrow(W2), ncol(W2), nrow(W12), ncol(W12)))

# ---- 3. Stack H and record column blocks ----
H  <- cbind(1, W1, W2, W12)                 # intercept in column 1
d1 <- ncol(W1); d2 <- ncol(W2); d12 <- ncol(W12)

# Column ranges of each block, 0-indexed after the intercept
# (the sampler uses 1 + idx internally)
col_map_main <- list(seq(1, d1), seq(d1 + 1, d1 + d2))
col_map_int  <- list(seq(d1 + d2 + 1, d1 + d2 + d12))

# ---- 4. Run the collapsed Gibbs sampler ----
t0 <- proc.time()
gs <- gibbs_interaction_sampler(
  y = y, H = H, D = D, nu = true_nu,
  col_map_main = col_map_main, K_main_list = list(K1, K2),
  col_map_int  = col_map_int,  K_int_list  = list(K12),
  n_iter = n_iter, n_burn = n_burn, verbose = TRUE
)
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Sampler: %d iterations in %.0f sec\n", n_iter, elapsed))

# ---- 5. Posterior summaries ----
eta_hat <- colMeans(gs$eta_samples)
f1_hat  <- as.vector(W1  %*% eta_hat[1 + col_map_main[[1]]])
f2_hat  <- as.vector(W2  %*% eta_hat[1 + col_map_main[[2]]])
f12_hat <- as.vector(W12 %*% eta_hat[1 + col_map_int[[1]]])

# Each component is identified up to a constant (absorbed by the
# intercept), so compare after centering.
ctr  <- function(v) v - mean(v)
rmse <- function(a, b) sqrt(mean((ctr(a) - ctr(b))^2))

cat(sprintf("RMSE (centered)  f1 = %.3f   f2 = %.3f   f12 = %.3f\n",
            rmse(f1_hat, true_f1), rmse(f2_hat, true_f2), rmse(f12_hat, true_f12)))
cat(sprintf("sigma2  posterior mean = %.3f   (true %.2f)\n", mean(gs$sigma2_samples), true_sigma2))
cat(sprintf("tau2    posterior mean = %.3f   (true %.2f)\n", mean(gs$tau2_samples),   true_tau2))
cat(sprintf("rho     posterior mean = %.3f   (true %.2f)\n", mean(gs$rho_samples),    true_rho))

# ---- 6. Plots ----
outfile <- file.path(root, "toy", "toy_plots.pdf")
pdf(outfile, width = 10, height = 8)

par(mfrow = c(2, 2))
plot(gs$sigma2_samples, type = "l", main = "sigma2 trace", ylab = "sigma2")
abline(h = true_sigma2, col = "red", lty = 2)
plot(gs$tau2_samples, type = "l", main = "tau2 trace", ylab = "tau2")
abline(h = true_tau2, col = "red", lty = 2)
plot(gs$rho_samples, type = "l", main = "rho trace", ylab = "rho")
abline(h = true_rho, col = "red", lty = 2)
plot(gs$tau2_s_int_samp[, 1], type = "l", main = "interaction smoothing variance trace",
     ylab = "tau2_s_12")

par(mfrow = c(1, 2))
o <- order(X1)
plot(X1[o], ctr(true_f1)[o], type = "l", lwd = 2, col = "red",
     main = "f1(x1): true vs posterior mean", xlab = "x1", ylab = "f1")
lines(X1[o], ctr(f1_hat)[o], lwd = 2, lty = 2, col = "steelblue")
legend("topright", c("true", "estimate"), col = c("red", "steelblue"), lwd = 2, lty = 1:2)
o <- order(X2)
plot(X2[o], ctr(true_f2)[o], type = "l", lwd = 2, col = "red",
     main = "f2(x2): true vs posterior mean", xlab = "x2", ylab = "f2")
lines(X2[o], ctr(f2_hat)[o], lwd = 2, lty = 2, col = "steelblue")
legend("topright", c("true", "estimate"), col = c("red", "steelblue"), lwd = 2, lty = 1:2)

par(mfrow = c(1, 2))
pal  <- colorRampPalette(c("blue", "white", "red"))(100)
brks <- seq(-1.6, 1.6, length.out = 101)
plot(X1, X2, pch = 16, cex = 0.8, col = pal[cut(ctr(true_f12), brks, labels = FALSE)],
     main = "f12(x1, x2): true", xlab = "x1", ylab = "x2")
plot(X1, X2, pch = 16, cex = 0.8, col = pal[cut(ctr(f12_hat), brks, labels = FALSE)],
     main = "f12(x1, x2): posterior mean", xlab = "x1", ylab = "x2")

invisible(dev.off())
cat("Plots written to toy/toy_plots.pdf\n")
