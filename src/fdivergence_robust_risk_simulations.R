############################################################
## F-divergence Toolkit – Robust Risk Assessment Examples
##
## This script implements:
##  - Preliminary functions for F-divergence worst-case analysis
##  - Experiment 1: Inventory pooling (lognormal case)
##  - Experiment 2: Proportional reinsurance under dependence
##  - Additional Example 1: Stability test on standardizers
##  - Additional Example 2: Generalized power curve
##  - Additional Example 3: Nonparametric reference model
##
## Requirements:
##   install.packages("nloptr")
##   install.packages("numDeriv")
##   install.packages("rpart")
############################################################

# Load required packages (stop with a clear message if not installed)
if (!requireNamespace("nloptr", quietly = TRUE)) {
  stop("Package 'nloptr' is required but not installed.")
}
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Package 'numDeriv' is required but not installed.")
}
if (!requireNamespace("rpart", quietly = TRUE)) {
  stop("Package 'rpart' is required but not installed.")
}

library(nloptr)   # For non-linear optimization
library(numDeriv) # For numerical differentiation
library(rpart)    # For classification trees

################################################################################
################################################################################
#####                         Preliminary Functions                        #####
################################################################################
################################################################################

# Estimate alpha = (alpha1, alpha2) for worst-case density
alpha_est <- function(X, kappa, F, T) {
  # Initial values for optimization
  a2 <- sqrt(2 * kappa / var(X))
  a1 <- -2 * mean(X)
  
  # Objective function to be minimized
  eval_f0 <- function(a) {
    val_F <- mean(F(T(a[1] + a[2] * X)))
    val_T <- mean(T(a[1] + a[2] * X))
    # Penalize deviation from F-constraint and normalization of T
    abs(val_F - kappa) * 10 + (val_T - 1)^2
  }
  
  res <- nloptr::nloptr(
    x0   = c(a1, a2),
    eval_f = eval_f0,
    lb   = c(-Inf, -Inf),
    ub   = c(Inf, Inf),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel  = 1e-12
    )
  )
  alpha <- res$solution
  return(alpha)
}

# Worst-case expectation of X given alpha and T
WC_est <- function(X, alpha, T) {
  mean(X * T(alpha[1] + alpha[2] * X))
}

# Asymptotic variance of the worst-case expectation estimator
v_est <- function(U, X, F, T) {
  TpU  <- numDeriv::grad(T, U)
  FpTU <- numDeriv::grad(F, T(U))
  
  B11 <- mean(TpU)
  B12 <- mean(FpTU * TpU)
  B21 <- mean(TpU * X)
  B22 <- mean(FpTU * TpU * X)
  
  w1 <- mean(X * TpU)
  w2 <- mean((X^2) * TpU)
  
  B <- rbind(
    c(B11, B12, 0),
    c(B21, B22, 0),
    c(0,   0,   1)
  )
  w <- rbind(w1, w2, -1)
  
  Sig <- cov(cbind(T(U), F(T(U)), X * T(U)))
  
  v2 <- t(w) %*% t(solve(B)) %*% Sig %*% solve(B) %*% w
  return(as.numeric(v2))
}

# Construct F from H by numerical integration + spline interpolation
F_spline <- function(H) {
  # Partition of input range [0, 10]
  par   <- seq(0, 10, length.out = 10^5)
  delta <- abs(par[1] - par[length(par)]) / length(par)
  
  # Remove the point 0 to avoid log(0)
  par <- par[-1]
  n   <- length(par)
  
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- H(log(par[i]))
  }
  
  # Integration split into 2 parts (for numerical stability / shape)
  split_idx <- round(n / 10)
  out1      <- out[1:split_idx]
  out2      <- out[(split_idx + 1):n]
  
  F1_val <- (cumsum(out1) - sum(out1)) * delta
  F2_val <- cumsum(out2) * delta
  
  F_val <- c(F1_val, F2_val)
  return(splinefun(par, F_val, method = "fmm"))
}

# Toolkit for Lognormal case
Lognormal_Toolkit <- function(r, sig, theta) {
  
  H_log <- function(y) {
    if (y >= 0) {
      return((exp((r * y * (theta * sig)^r + 1)^(1 / r) - 1) - 1) / (theta * sig)^r)
    }
    if (y < 0) {
      return(y)
    }
  }
  H_log <- Vectorize(H_log)
  
  H_inv_log <- function(x) {
    if (x >= 0) {
      return(((log((theta * sig)^r * x + 1) + 1)^r - 1) / (r * (theta * sig)^r))
    }
    if (x < 0) {
      return(x)
    }
  }
  H_inv_log <- Vectorize(H_inv_log)
  
  F_log <- F_spline(H_log)
  
  T_log <- function(x) {
    exp(H_inv_log(x))
  }
  
  list(H = H_log, H_inv = H_inv_log, T = T_log, F = F_log)
}

# Toolkit for Weibull case
Weibull_Toolkit <- function(k, theta) {
  
  H_w <- function(y) {
    if (y >= 0) {
      return(((y + 1)^(theta / k) - 1) * k / theta)
    }
    if (y < 0) {
      return(y)
    }
  }
  H_w <- Vectorize(H_w)
  
  H_inv_w <- function(x) {
    if (x >= 0) {
      return((theta * x / k + 1)^(k / theta) - 1)
    }
    if (x < 0) {
      return(x)
    }
  }
  H_inv_w <- Vectorize(H_inv_w)
  
  F_w <- F_spline(H_w)
  
  T_w <- function(x) {
    exp(H_inv_w(x))
  }
  
  list(H = H_w, H_inv = H_inv_w, T = T_w, F = F_w)
}

################################################################################
################################################################################
#####              Simulation Experiment 1: Inventory Pooling              #####
################################################################################
################################################################################

Toolkit <- Lognormal_Toolkit(r = 2, sig = sqrt(2), theta = 2)

# Cost functions
Cd <- function(q, X) {
  h * sum(pmax(q - X, 0)) + p * sum(pmax(X - q, 0))
}
Cc <- function(Q, X) {
  h * max(Q - sum(X), 0) + p * max(sum(X) - Q, 0)
}

nsim  <- 2^13
n_max <- 40
a     <- 100
b     <- sqrt(2)
h     <- 1
p     <- 1

# Simulate demands: D[trial, location]
D <- array(NA_real_, dim = c(nsim, n_max))
for (isim in seq_len(nsim)) {
  set.seed(isim)
  D[isim, ] <- a * exp(b * rnorm(n_max))
}

sig_wc    <- sig_bc    <- numeric(n_max - 1)
EX_nominal <- EX_wc    <- EX_bc <- numeric(n_max - 1)

for (n in 2:n_max) {
  
  # Cost for decentralized case
  Cd_obj <- function(q) {
    out <- 0
    for (isim in seq_len(nsim)) {
      out <- out + Cd(q, D[isim, 1:n])
    }
    out / nsim
  }
  
  # Optimal inventory per location (decentralized)
  qv <- round(optimize(Cd_obj, c(0, 40000))$minimum)
  
  # Cost for centralized case
  Cc_obj <- function(Q) {
    out <- 0
    for (isim in seq_len(nsim)) {
      out <- out + Cc(Q, D[isim, 1:n])
    }
    out / nsim
  }
  
  # Optimal total inventory (centralized)
  Qv <- round(optimize(Cc_obj, c(0, 40000))$minimum)
  
  # Nominal gains per location from pooling
  X <- numeric(nsim)
  for (isim in seq_len(nsim)) {
    X[isim] <- (Cd(qv, D[isim, 1:n]) - Cc(Qv, D[isim, 1:n])) / n
  }
  
  EX_nominal[n - 1] <- mean(X)
  
  # Standardize for numerical stability
  Z <- (X - mean(X)) / sd(X)
  
  # Worst-case mean
  alpha_wc        <- alpha_est(Z, kappa = 0.1, F = Toolkit$F, T = Toolkit$T)
  EX_wc[n - 1]    <- WC_est(Z, alpha_wc, T = Toolkit$T) * sd(X) + mean(X)
  U_wc            <- alpha_wc[1] + alpha_wc[2] * Z
  v_wc            <- v_est(U_wc, Z, F = Toolkit$F, T = Toolkit$T) * var(X)
  sig_wc[n - 1]   <- sqrt(v_wc / 1000)
  
  # Best-case mean
  alpha_bc        <- alpha_est(-Z, kappa = 0.1, F = Toolkit$F, T = Toolkit$T)
  EX_bc[n - 1]    <- -WC_est(-Z, alpha_bc, T = Toolkit$T) * sd(X) + mean(X)
  U_bc            <- alpha_bc[1] + alpha_bc[2] * (-Z)
  v_bc            <- v_est(U_bc, -Z, F = Toolkit$F, T = Toolkit$T) * var(X)
  sig_bc[n - 1]   <- sqrt(v_bc / 1000)
}

plot(
  2:n_max, EX_wc,
  type = "l",
  ylim = c(0, 200),
  col = "darkred",
  lwd = 2,
  ylab = "Gains from inventory pooling",
  xlab = "n",
  xaxt = "none",
  yaxt = "none",
  main = "Reproduction of Experiment 1: Inventory Pooling"
)
axis(1, seq(0, 40, 5))
axis(2, seq(0, 200, 20))

lines(2:n_max, EX_wc + 2 * sig_wc, col = "darkred", lwd = 2, lty = 3)
lines(2:n_max, EX_bc,              col = "darkblue", lwd = 2)
lines(2:n_max, EX_bc - 2 * sig_bc, col = "darkblue", lwd = 2, lty = 3)
lines(2:n_max, EX_nominal,         col = "black",   lwd = 2, lty = 2)

legend(
  "bottomright",
  c("Worst Case Mean", "Nominal Mean", "Best Case Mean"),
  bg  = "white",
  col = c("darkred", "black", "darkblue"),
  lty = c(1, 2, 1),
  lwd = c(2, 2, 2),
  cex = 1.2
)

################################################################################
################################################################################
### Experiment 2: Proportional Reinsurance under different Claim Dependences ###
################################################################################
################################################################################

Toolkit <- Weibull_Toolkit(k = 0.3, theta = 2)

# Model parameters
rho <- 0.3
n   <- 5
k   <- seq(0.3, 0.7, length.out = n)
lambda <- 1.1^(-1 / k)

# Copula covariance matrix
v     <- 2
Sig_c <- rbind(
  c(1,   -0.25, 0.25, 0,    0.25),
  c(-0.25, 1,   0,    0.25, 0.25),
  c(0.25, 0,    1,    0,    0.25),
  c(0,    0.25, 0,    1,    0.5),
  c(0.25, 0.25, 0.25, 0.5,  1)
)

nsim <- 80000

L_iid <- L_G <- L_t <- M_t <- array(NA_real_, dim = c(nsim, n))

set.seed(4102)

# Step 1: iid uniforms
U_iid <- matrix(runif(nsim * n), ncol = n)

# Step 2: iid Normal -> multivariate Normal via Cholesky
M_G <- qnorm(U_iid) %*% chol(Sig_c)

# Step 3: multivariate t-distribution
for (isim in seq_len(nsim)) {
  M_t[isim, ] <- M_G[isim, ] * sqrt(v / rchisq(1, v))
}

# Step 4: copula samples
U_G <- pnorm(M_G)
U_t <- pt(M_t, df = v)

# Step 5: Weibull marginals
for (i in seq_len(n)) {
  L_iid[, i] <- qweibull(U_iid[, i], k[i], lambda[i])
  L_G[, i]   <- qweibull(U_G[, i],   k[i], lambda[i])
  L_t[, i]   <- qweibull(U_t[, i],   k[i], lambda[i])
}

R_iid <- rho * rowSums(L_iid)
R_G   <- rho * rowSums(L_G)
R_t   <- rho * rowSums(L_t)

Z_iid <- (R_iid - mean(R_iid)) / sd(R_iid)
Z_G   <- (R_G   - mean(R_G))   / sd(R_G)
Z_t   <- (R_t   - mean(R_t))   / sd(R_t)

kappa  <- seq(0, 0.2, length.out = 21)
EX_iid <- EX_G <- EX_t <- numeric(length(kappa))

for (i in seq_along(kappa)) {
  
  alpha_iid <- alpha_est(Z_iid, kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_iid[i] <- WC_est(Z_iid, alpha_iid, T = Toolkit$T) * sd(R_iid) + mean(R_iid)
  
  alpha_G   <- alpha_est(Z_G,   kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_G[i]   <- WC_est(Z_G, alpha_G, T = Toolkit$T)   * sd(R_G)   + mean(R_G)
  
  alpha_t   <- alpha_est(Z_t,   kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_t[i]   <- WC_est(Z_t, alpha_t, T = Toolkit$T)   * sd(R_t)   + mean(R_t)
  
  if (i %% 3 == 0) cat(i, "-th trial >> \n")
}

plot(
  kappa, EX_iid,
  type = "l",
  lwd  = 2,
  ylim = c(3.5, 7.5),
  xaxt = "none",
  yaxt = "none",
  xlab = expression(kappa),
  ylab = "Worst-case expected value",
  main = "Reproduction of Experiment 2"
)
axis(1, seq(0, 0.2, 0.02))
axis(2, seq(3.5, 7.5, 0.5))

lines(kappa, EX_G, lwd = 2, lty = 2, col = "darkred")
lines(kappa, EX_t, lwd = 2, lty = 4, col = "darkblue")
abline(h = 4, lty = 2)

legend(
  "bottomright",
  c("IID", "Gaussian Copula", "t - copula"),
  bg  = "white",
  col = c("black", "darkred", "darkblue"),
  lty = c(1, 2, 4),
  lwd = c(2, 2, 2),
  cex = 1.2
)

################################################################################
################################################################################
###           Additional Example 1: Stability Test on Standardizers          ###
################################################################################
################################################################################

Toolkit <- Weibull_Toolkit(k = 1, theta = 2)

nsim <- 2^12
sig_A <- sig_B <- numeric(nsim)

for (isim in seq_len(nsim)) {
  set.seed(isim)
  X <- rpois(500, lambda = pi)
  sig_A[isim] <- sqrt(mean(X))
  sig_B[isim] <- sqrt(var(X))
}

Z_A <- (sig_A - mean(sig_A)) / sd(sig_A)
Z_B <- (sig_B - mean(sig_B)) / sd(sig_B)

kappa <- seq(0, 0.5, length.out = 21)
EX_A  <- EX_B <- numeric(length(kappa))

for (i in seq_along(kappa)) {
  alpha_A <- alpha_est(Z_A, kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_A[i] <- WC_est(Z_A, alpha_A, T = Toolkit$T) * sd(sig_A) + mean(sig_A)
  
  alpha_B <- alpha_est(Z_B, kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_B[i] <- WC_est(Z_B, alpha_B, T = Toolkit$T) * sd(sig_B) + mean(sig_B)
  
  if (i %% 3 == 0) cat(i, "-th trial >> \n")
}

plot(
  kappa, EX_B,
  col  = "darkred",
  type = "l",
  main = "Comparison of Standardizers",
  xlab = expression(kappa),
  ylab = "E(X)",
  xaxt = "none",
  yaxt = "none",
  lwd  = 2
)
axis(1, seq(0, 0.5, 0.05))
axis(2, seq(1.77, 1.85, 0.01))

lines(kappa, EX_A, col = "darkblue", type = "l", lwd = 2)

legend(
  "bottomright",
  c("A: Parametric", "B: Non-parametric"),
  bg  = "white",
  col = c("darkblue", "darkred"),
  lty = c(1, 1),
  lwd = c(2, 2),
  cex = 1.2
)

################################################################################
################################################################################
###               Additional Example 2: Generalized Power Curve              ###
################################################################################
################################################################################

Toolkit <- Weibull_Toolkit(k = 0.5, theta = 2)

ks.test0 <- function(x, F0) {
  n <- length(x)
  i <- seq_len(n)
  u <- F0(sort(x))
  Kpos <- max(i / n - u)
  Kneg <- max(u - (i - 1) / n)
  K <- sqrt(n) * max(Kpos, Kneg)
  K
}

ad.test0 <- function(x, F0) {
  n <- length(x)
  u <- F0(sort(x))
  i <- seq_len(n)
  A <- -n - (sum((2 * i - 1) * (log(u) + log(1 - rev(u)))))/n
  A
}

nsim <- 2000
n    <- 200

# Critical values for n = 200, alpha = 10%
q_ks <- 1.208
q_ad <- 1.923

psi_A <- psi_B <- numeric(nsim)

F0 <- function(x) 0.3 * pnorm(x, 0, 1) + 0.7 * pnorm(x, 1, 2)

for (isim in seq_len(nsim)) {
  set.seed(isim)
  I <- rbinom(n, 1, 0.3)
  X <- I * rnorm(n, 0, 1) + (1 - I) * rnorm(n, 1, 2)
  psi_A[isim] <- as.numeric(ks.test0(X, F0) > q_ks)
  psi_B[isim] <- as.numeric(ad.test0(X, F0) > q_ad)
}

kappa  <- seq(0, 1, length.out = 21)
EX_pA  <- EX_pB <- numeric(length(kappa))

for (i in seq_along(kappa)) {
  alpha_pA <- alpha_est(psi_A, kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_pA[i] <- WC_est(psi_A, alpha_pA, T = Toolkit$T)
  
  alpha_pB <- alpha_est(psi_B, kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_pB[i] <- WC_est(psi_B, alpha_pB, T = Toolkit$T)
  
  if (i %% 3 == 0) cat(i, "-th trial >> \n")
}

plot(
  kappa, EX_pB,
  col  = "darkred",
  type = "l",
  main = "Comparison of Tests",
  xlab = expression(kappa),
  ylab = "Power",
  xaxt = "none",
  yaxt = "none",
  lwd  = 2
)
axis(1, seq(0, 1, 0.1))
axis(2, seq(0, 1, 0.1))

lines(kappa, EX_pA, col = "darkblue", type = "l", lwd = 2)

legend(
  "bottomright",
  c("A: KS-Test", "B: AD-Test"),
  bg  = "white",
  col = c("darkblue", "darkred"),
  lty = c(1, 1),
  lwd = c(2, 2),
  cex = 1.2
)

################################################################################
################################################################################
###           Additional Example 3: Nonparametric Reference Model            ###
################################################################################
################################################################################

## NOTE:
## Original large COVID dataset is from:
## https://www.kaggle.com/datasets/gpreda/coronavirus-2019ncov
## Here we assume a cleaned series has been saved as:
##   data/COVIDHKdata.csv
## containing a single column "x" (daily new cases).

# X: cleaned daily new cases (Hong Kong)
X <- read.csv(file.path("data", "COVIDHKdata.csv"))$x
n <- length(X)  # number of observations (fixes missing n in original code)

Toolkit <- Weibull_Toolkit(k = 2, theta = 3)

nsim  <- 2^12
X_bar <- numeric(nsim)

# Bootstrap replicate of the sample mean
for (isim in seq_len(nsim)) {
  set.seed(isim)
  X_b <- X[sample(seq_len(n), n, replace = TRUE)]
  X_bar[isim] <- mean(X_b)
}

Z_bar <- (X_bar - mean(X_bar)) / sd(X_bar)

kappa   <- seq(0, 1, length.out = 21)
EX_bar  <- sig_bar <- numeric(length(kappa))

for (i in seq_along(kappa)) {
  alpha_bar   <- alpha_est(Z_bar, kappa = kappa[i], F = Toolkit$F, T = Toolkit$T)
  EX_bar[i]   <- WC_est(Z_bar, alpha_bar, T = Toolkit$T) * sd(X_bar) + mean(X_bar)
  
  if (i %% 3 == 0) cat(i, "-th trial >> \n")
  
  U_bar       <- alpha_bar[1] + alpha_bar[2] * Z_bar
  v_bar       <- v_est(U_bar, Z_bar, F = Toolkit$F, T = Toolkit$T) * var(X_bar)
  sig_bar[i]  <- sqrt(v_bar)
}

Upper <- EX_bar + 0.5 * sig_bar
Lower <- EX_bar - 0.5 * sig_bar

plot(
  kappa, Upper,
  col  = "darkred",
  lty  = 2,
  lwd  = 3,
  type = "l",
  xlab = expression(kappa),
  ylab = "sup E(X)",
  ylim = c(23, 32),
  main = "Hong Kong COVID data application"
)
polygon(
  c(kappa, rev(kappa)),
  c(Lower, rev(Upper)),
  col = rgb(1, 0, 0, 0.5)
)
lines(kappa, Lower, col = "darkred", lty = 2, lwd = 3)
lines(kappa, EX_bar, col = "darkred", lwd = 3, type = "l")

legend(
  "bottomright",
  c("Worst Case Expectation", "Confidence Bound"),
  bg  = "white",
  col = c("darkred", "darkred"),
  lty = c(1, 2),
  lwd = c(2, 2),
  cex = 1.2
)
