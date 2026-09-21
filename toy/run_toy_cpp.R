# ============================================================
# toy/run_toy_cpp.R
#
# Same example as run_toy.R, with the C++ helpers compiled.
# Needs Rcpp, RcppEigen, and a C++ compiler (Rtools on Windows,
# Xcode command-line tools on macOS, gcc/g++ on Linux).
#
# What it adds over run_toy.R:
#   - compiles R/ls_interaction_core.cpp
#   - checks that the C++ and R versions of the interaction design agree
#   - reports how much faster the C++ version is
#   - then runs the full example with the C++ path enabled
#
# Run from the repository root:
#   source("toy/run_toy_cpp.R")
#   Rscript toy/run_toy_cpp.R
# ============================================================

root <- if (file.exists("R/ls_basis.R")) "." else ".."

# ---- Compile ----
Rcpp::sourceCpp(file.path(root, "R", "ls_interaction_core.cpp"))
cat("Compiled R/ls_interaction_core.cpp\n")

# ---- Check R and C++ interaction designs agree ----
source(file.path(root, "R", "ls_basis.R"))
source(file.path(root, "R", "ls_interaction.R"))

set.seed(1)
n <- 1000; M <- 20
obj_u <- ls_build_one_full(runif(n), M = M)
obj_v <- ls_build_one_full(runif(n), M = M)

t_R   <- system.time(W_R   <- khatri_rao_rowwise_R(obj_u$W, obj_v$W))["elapsed"]
t_cpp <- system.time(W_cpp <- khatri_rao_cpp(obj_u$W, obj_v$W))["elapsed"]

cat(sprintf("Interaction design %d x %d:  R %.3f s,  C++ %.3f s,  max |diff| = %.1e\n",
            nrow(W_R), ncol(W_R), t_R, t_cpp, max(abs(W_R - W_cpp))))
stopifnot(max(abs(W_R - W_cpp)) < 1e-10)

# ---- Run the example with the C++ path ----
use_cpp <- TRUE
source(file.path(root, "toy", "run_toy.R"))
