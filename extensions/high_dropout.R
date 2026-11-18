# =============================================================================
# Extension: what happens when about half of the patients drop out?
# Same design as the main study, true effect 0.4, but a higher dropout rate.
# Run from the repository root:  source("extensions/high_dropout.R")
# =============================================================================

source("R/functions.R")
alpha_drop <- -1.5      # higher dropout (main study: -2.2)
n_sim <- 500

mechs <- c("MCAR", "MAR", "MNAR")
res <- list(); drop <- c()
for (s in seq_along(mechs)) {
  cat("Running", mechs[s], "...\n")
  r <- run_scenario(delta = 0.4, mechanism = mechs[s], scen_id = 100 + s)  # new seeds, not the main study's
  drop[s] <- attr(r, "drop_rate")
  res[[s]] <- r
}
raw <- do.call(rbind, res)

perf <- do.call(rbind, lapply(split(raw, list(raw$mechanism, raw$method)), function(d)
  cbind(mechanism = d$mechanism[1], method = d$method[1], summarise_perf(d))))
perf <- perf[order(match(perf$mechanism, mechs), perf$method), ]
rownames(perf) <- NULL

cat("\nShare missing at the last visit:", paste(mechs, round(drop, 2), collapse = ", "), "\n\n")
print(perf[, c("mechanism", "method", "bias", "bias_mcse", "coverage", "reject")], digits = 3)
write.csv(perf, "results/high_dropout.csv", row.names = FALSE)
cat("\nSaved to results/high_dropout.csv\n")
