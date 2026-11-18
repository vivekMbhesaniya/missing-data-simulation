# =============================================================================
# Functions for the simulation study: how missing data affects the treatment effect estimate
# in a randomized two-arm trial with a longitudinal endpoint
#
# Methods compared (all estimate the treatment difference at the last visit):
#   1. CC    : complete-case ANCOVA (only patients who finished the study)
#   2. LOCF  : last observation carried forward, then ANCOVA
#   3. MMRM  : mixed model for repeated measures, unstructured covariance
#              (fitted with nlme::gls, which ships with R)
#
# Missing-data mechanisms: none, MCAR, MAR, MNAR
# Performance measures: bias, empirical SE, model SE, 95% CI coverage,
# rejection rate (type I error when delta = 0, power when delta > 0),
# each with its Monte Carlo standard error (MCSE).
#
# Planned with the ADEMP structure from Morris, White & Crowther (2019),
# "Using simulation studies to evaluate statistical methods", Stat Med.
# =============================================================================


# This file only defines settings and functions. Nothing runs when it is sourced.
# run_simulation.R runs the study; checks/sanity_checks.R tests the functions.

library(nlme)

# ---- Settings ---------------------------------------------------------------
n_per_arm <- 100        # patients per arm
n_visits  <- 4          # post-baseline visits (e.g. weeks 4, 8, 12, 16)
rho       <- 0.6        # correlation between neighbouring visits (AR(1) structure)
sd_y      <- 1          # residual SD per visit (outcome on a standardized scale)
b_base    <- 0.5        # effect of baseline value on later visits
n_sim     <- 500      # replicates per scenario. Test with 50 first (about 3 min);
                        # 500 takes roughly 30 to 45 min on a laptop
alpha_drop <- -2.2      # sets the dropout rate: about 27 % missing at the last visit under MCAR

# ---- D: Data-generating mechanism -------------------------------------------
# Higher y = more improvement. The treatment effect grows linearly over visits
# and equals delta at the last visit.
simulate_trial <- function(delta) {
  n   <- 2 * n_per_arm
  arm <- rep(0:1, each = n_per_arm)
  base <- rnorm(n)

  # AR(1) covariance between visits
  Sigma <- sd_y^2 * rho^abs(outer(1:n_visits, 1:n_visits, "-"))
  E <- matrix(rnorm(n * n_visits), n) %*% chol(Sigma)

  trend <- outer(arm, (1:n_visits) / n_visits) * delta
  Y <- b_base * base + trend + E
  list(arm = arm, base = base, Y = Y)
}

# Monotone dropout: once a patient drops out, all later visits are missing.
# Visit 1 is always observed. alpha sets the overall dropout rate,
# beta sets how strongly dropout depends on the outcome.
apply_dropout <- function(trial, mechanism, alpha = alpha_drop, beta = 1.0) {
  Y <- trial$Y
  n <- nrow(Y)
  if (mechanism == "none") return(Y)
  dropped <- rep(FALSE, n)
  for (j in 2:n_visits) {
    p <- switch(mechanism,
      MCAR = rep(plogis(alpha), n),                 # unrelated to anything
      MAR  = plogis(alpha - beta * Y[, j - 1]),     # poor observed improvement -> drop out
      MNAR = plogis(alpha - beta * Y[, j])          # depends on the value we will never see
    )
    dropped <- dropped | (runif(n) < p)
    Y[dropped, j] <- NA
  }
  Y
}

# ---- M: Methods ---------------------------------------------------------------
# Each method returns: estimate, SE, lower and upper 95% CI, p-value
fit_ancova <- function(y_last, arm, base) {
  fit <- lm(y_last ~ arm + base)
  s <- summary(fit)$coefficients["arm", ]
  ci <- confint(fit)["arm", ]
  c(est = s[["Estimate"]], se = s[["Std. Error"]],
    lo = ci[[1]], hi = ci[[2]], p = s[["Pr(>|t|)"]])
}

method_cc <- function(Y, arm, base) {
  keep <- !is.na(Y[, n_visits])
  fit_ancova(Y[keep, n_visits], arm[keep], base[keep])
}

method_locf <- function(Y, arm, base) {
  y_last <- apply(Y, 1, function(r) tail(r[!is.na(r)], 1))
  fit_ancova(y_last, arm, base)
}

method_mmrm <- function(Y, arm, base) {
  n <- nrow(Y)
  long <- data.frame(
    id    = rep(seq_len(n), times = n_visits),
    visit = factor(rep(seq_len(n_visits), each = n)),
    arm   = rep(arm, times = n_visits),
    base  = rep(base, times = n_visits),
    y     = as.vector(Y)
  )
  long <- long[!is.na(long$y), ]
  long$vnum <- as.integer(long$visit)

  # Visit-specific treatment and baseline effects, unstructured covariance:
  # corSymm = free correlations, varIdent = free variance per visit.
  fit <- gls(y ~ visit + visit:arm + visit:base,
             data = long,
             correlation = corSymm(form = ~ vnum | id),
             weights = varIdent(form = ~ 1 | visit),
             method = "REML", na.action = na.omit)

  term <- paste0("visit", n_visits, ":arm")
  est <- coef(fit)[[term]]
  se  <- sqrt(vcov(fit)[term, term])
  # gls uses residual df (no Satterthwaite); fine for n = 200, note it as a limitation
  df  <- fit$dims$N - fit$dims$p
  q   <- qt(0.975, df)
  c(est = est, se = se, lo = est - q * se, hi = est + q * se,
    p = 2 * pt(-abs(est / se), df))
}

methods <- list(CC = method_cc, LOCF = method_locf, MMRM = method_mmrm)

# ---- Run one scenario ---------------------------------------------------------
run_scenario <- function(delta, mechanism, scen_id) {
  res <- vector("list", n_sim)
  drop_rate <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    set.seed(base_seed + 100000 * scen_id + i)   # each replicate reproducible on its own
    trial <- simulate_trial(delta)
    Y <- apply_dropout(trial, mechanism)
    drop_rate[i] <- mean(is.na(Y[, n_visits]))
    res[[i]] <- do.call(rbind, lapply(names(methods), function(m) {
      out <- tryCatch(methods[[m]](Y, trial$arm, trial$base),
                      error = function(e) c(est = NA, se = NA, lo = NA, hi = NA, p = NA))
      data.frame(rep = i, method = m, t(out))
    }))
  }
  out <- do.call(rbind, res)
  out$delta <- delta
  out$mechanism <- mechanism
  attr(out, "drop_rate") <- mean(drop_rate)
  out
}

# ---- E: Estimand and P: Performance measures -----------------------------------
# Estimand: difference in mean outcome at the last visit between arms,
# had all patients stayed in the study (hypothetical strategy, ICH E9(R1)).
summarise_perf <- function(d) {
  d <- d[!is.na(d$est), ]
  true <- d$delta[1]
  k <- nrow(d)
  bias  <- mean(d$est) - true
  empse <- sd(d$est)
  cover <- mean(d$lo <= true & d$hi >= true)
  rej   <- mean(d$p < 0.05)
  data.frame(
    n_ok       = k,
    bias       = bias,
    bias_mcse  = empse / sqrt(k),
    emp_se     = empse,
    model_se   = mean(d$se),
    coverage   = cover,
    cover_mcse = sqrt(cover * (1 - cover) / k),
    reject     = rej,
    reject_mcse = sqrt(rej * (1 - rej) / k)
  )
}

