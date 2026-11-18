# =============================================================================
# Sanity checks: does the code do what the design says?
# Run from the repository root:  Rscript checks/sanity_checks.R
# Every check stops with an error if it fails.
# =============================================================================

source("R/functions.R")
set.seed(1)

# 1. Data generation: with a very large trial, the simulated data should show
#    the parameters we put in (treatment effect at the last visit = delta,
#    correlation between neighbouring visits = rho).
n_per_arm_saved <- n_per_arm
n_per_arm <- 20000
big <- simulate_trial(delta = 0.4)
diff_last <- mean(big$Y[big$arm == 1, n_visits]) - mean(big$Y[big$arm == 0, n_visits])
resid <- big$Y - b_base * big$base - outer(big$arm, (1:n_visits) / n_visits) * 0.4
cor12 <- cor(resid[, 1], resid[, 2])
cat(sprintf("Check 1: effect at last visit = %.3f (target 0.4), corr(visit 1, 2) = %.3f (target %.1f)\n",
            diff_last, cor12, rho))
stopifnot(abs(diff_last - 0.4) < 0.03, abs(cor12 - rho) < 0.02)

# 2. Dropout is monotone (once missing, always missing) and visit 1 is always observed.
for (mech in c("MCAR", "MAR", "MNAR")) {
  Y <- apply_dropout(big, mech)
  stopifnot(!anyNA(Y[, 1]))
  for (j in 2:n_visits) stopifnot(all(is.na(Y[is.na(Y[, j - 1]), j])))
}
cat("Check 2: dropout is monotone, visit 1 always observed\n")

# 3. MCAR dropout rate matches theory: P(missing at last visit) = 1 - (1 - p)^(visits - 1)
p <- plogis(alpha_drop)
theory <- 1 - (1 - p)^(n_visits - 1)
observed <- mean(is.na(apply_dropout(big, "MCAR")[, n_visits]))
cat(sprintf("Check 3: MCAR dropout at last visit = %.3f (theory %.3f)\n", observed, theory))
stopifnot(abs(observed - theory) < 0.01)

# 4. MAR dropout really depends on the observed outcome: patients who dropped out
#    after visit 1 had a lower visit-1 value than those who stayed.
Y <- apply_dropout(big, "MAR")
gone <- is.na(Y[, 2])
cat(sprintf("Check 4: mean visit-1 value, dropped = %.2f, stayed = %.2f\n",
            mean(Y[gone, 1]), mean(Y[!gone, 1])))
stopifnot(mean(Y[gone, 1]) < mean(Y[!gone, 1]))
n_per_arm <- n_per_arm_saved

# 5. With no missing data, all three methods must give the same estimate:
#    CC and LOCF use the same complete data, and MMRM with an unstructured
#    covariance and visit-specific effects reduces to ANCOVA at the last visit.
trial <- simulate_trial(delta = 0.4)
ests <- sapply(methods, function(m) m(trial$Y, trial$arm, trial$base)[["est"]])
cat("Check 5: estimates without missing data:", paste(names(ests), round(ests, 6), collapse = ", "), "\n")
stopifnot(max(ests) - min(ests) < 1e-6)

# 6. Reproducibility: the same seed gives exactly the same replicate.
set.seed(42); tr <- simulate_trial(0.4); Y <- apply_dropout(tr, "MAR")
set.seed(42); tr2 <- simulate_trial(0.4); Y2 <- apply_dropout(tr2, "MAR")
stopifnot(identical(Y, Y2))
cat("Check 6: same seed gives identical data\n")

cat("\nAll checks passed.\n")
