# =============================================================================
# Runs the full simulation study and writes tables and figures to results/
# Run from the repository root:  Rscript run_simulation.R
# =============================================================================

library(ggplot2)
source("R/functions.R")

# Optional: number of replicates from the command line, e.g. Rscript run_simulation.R 50
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) n_sim <- as.integer(args[1])
cat("Replicates per scenario:", n_sim, "\n")

# ---- Run all scenarios ----------------------------------------------------------
scenarios <- expand.grid(delta = c(0, 0.4),
                         mechanism = c("none", "MCAR", "MAR", "MNAR"),
                         stringsAsFactors = FALSE)

t0 <- Sys.time()
all_res <- list(); drop_info <- c()
for (s in seq_len(nrow(scenarios))) {
  cat(sprintf("Scenario %d/%d: delta = %.1f, %s\n",
              s, nrow(scenarios), scenarios$delta[s], scenarios$mechanism[s]))
  r <- run_scenario(scenarios$delta[s], scenarios$mechanism[s], s)
  drop_info <- c(drop_info, attr(r, "drop_rate"))
  all_res[[s]] <- r
}
raw <- do.call(rbind, all_res)
cat("Run time:", format(Sys.time() - t0), "\n")
scenarios$dropout_last_visit <- round(drop_info, 3)
print(scenarios)

perf <- do.call(rbind, lapply(split(raw, list(raw$delta, raw$mechanism, raw$method)),
  function(d) cbind(delta = d$delta[1], mechanism = d$mechanism[1],
                    method = d$method[1], summarise_perf(d))))
perf$mechanism <- factor(perf$mechanism, levels = c("none", "MCAR", "MAR", "MNAR"))
perf <- perf[order(perf$delta, perf$mechanism, perf$method), ]
rownames(perf) <- NULL
print(perf, digits = 3)

dir.create("results", showWarnings = FALSE)
write.csv(raw,  "results/raw_estimates.csv", row.names = FALSE)
write.csv(perf, "results/performance.csv",   row.names = FALSE)

# ---- Figures ---------------------------------------------------------------------
p_bias <- ggplot(subset(perf, delta == 0.4),
                 aes(x = mechanism, y = bias, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = bias - 1.96 * bias_mcse, ymax = bias + 1.96 * bias_mcse),
                  position = position_dodge(width = 0.5)) +
  labs(title = "Bias of the treatment effect at the last visit (true effect 0.4)",
       x = "Missing-data mechanism", y = "Bias (estimate - truth)", colour = "Method") +
  theme_minimal(base_size = 13)

p_cov <- ggplot(subset(perf, delta == 0.4),
                aes(x = mechanism, y = coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_pointrange(aes(ymin = coverage - 1.96 * cover_mcse,
                      ymax = coverage + 1.96 * cover_mcse),
                  position = position_dodge(width = 0.5)) +
  labs(title = "Coverage of the 95% confidence interval",
       x = "Missing-data mechanism", y = "Coverage", colour = "Method") +
  theme_minimal(base_size = 13)

p_t1e <- ggplot(subset(perf, delta == 0),
                aes(x = mechanism, y = reject, colour = method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_pointrange(aes(ymin = reject - 1.96 * reject_mcse,
                      ymax = reject + 1.96 * reject_mcse),
                  position = position_dodge(width = 0.5)) +
  labs(title = "Type I error (true effect 0)",
       x = "Missing-data mechanism", y = "Rejection rate", colour = "Method") +
  theme_minimal(base_size = 13)

ggsave("results/bias.png",     p_bias, width = 8, height = 4.5, dpi = 150)
ggsave("results/coverage.png", p_cov,  width = 8, height = 4.5, dpi = 150)
ggsave("results/type1.png",    p_t1e,  width = 8, height = 4.5, dpi = 150)
cat("Done. Results in the 'results' folder.\n")
