# Missing Data in a Randomized Clinical Trial: A Simulation Study in R

How much does patient dropout distort the estimated treatment effect, and which analysis method holds up? This project compares three common methods across four missing-data mechanisms with a Monte Carlo simulation study.

## Design (ADEMP, Morris, White & Crowther 2019)

**Aim.** Compare bias, precision, confidence-interval coverage, type I error and power of three analysis methods when patients drop out of a two-arm trial.

**Data-generating mechanism.**
- Two arms, 100 patients each, one baseline value and 4 post-baseline visits
- Outcome on a standardized scale (higher = more improvement), correlated within patients (AR(1), rho = 0.6)
- True treatment effect grows linearly and reaches delta at the last visit: delta = 0 (null) and delta = 0.4
- Monotone dropout after visit 1, about 30 % missing at the last visit:
  - MCAR: dropout unrelated to the outcome
  - MAR: dropout more likely after a poor *observed* value (lack of efficacy)
  - MNAR: dropout depends on the current, *unobserved* value

**Estimand.** Difference in mean outcome at the last visit between treatment and control, had all patients stayed in the study (hypothetical strategy for dropout, ICH E9(R1)).

**Methods.**
1. Complete-case ANCOVA at the last visit, adjusted for baseline
2. LOCF (last observation carried forward), then ANCOVA
3. MMRM with unstructured covariance and visit-specific treatment and baseline effects, fitted with `nlme::gls`

**Performance measures.** Bias, empirical SE, average model SE, coverage of the 95 % CI, rejection rate at alpha = 0.05, each with its Monte Carlo standard error. 500 replicates per scenario, 8 scenarios.

## How to run

```r
# from the repository root, R 4.x with ggplot2 installed
Rscript checks/sanity_checks.R      # about 10 seconds
Rscript run_simulation.R 50         # quick test, a few minutes
Rscript run_simulation.R            # full run, 500 replicates per scenario
```

Results go to `results/` (performance table, raw estimates, three figures).

## How I checked the code

`checks/sanity_checks.R` verifies that:
- a very large simulated trial reproduces the effect size and correlation that were put in
- dropout is monotone and visit 1 is always observed
- the MCAR dropout rate matches the theoretical value 1 - (1 - p)^3
- under MAR, patients who drop out had lower observed values than those who stay
- without missing data, CC, LOCF and MMRM give identical estimates (as theory says)
- the same seed reproduces the same data exactly

Every replicate has its own seed, so any single replicate can be rerun on its own.

## Results

**Validation.** Without dropout, all three methods give identical estimates and the power is 81.4 %, close to the theoretical 80 % for an effect of 0.4 SD with 100 patients per arm.

**Type I error.** With no true effect, the rejection rate stays between 3.2 % and 5.6 % for every method and mechanism, which is within Monte Carlo error of the nominal 5 %.

**Bias and coverage.** LOCF underestimates the treatment effect in every dropout scenario, by 14 % under MCAR, 11 % under MAR and 18 % under MNAR, because it carries forward values from earlier visits where the treatment effect was still smaller. Under MAR its 95 % CI covers the true effect in only 91.6 % of trials. MMRM is unbiased under MCAR and MAR (bias within 1 MCSE) and keeps more power than complete-case analysis under MAR (68.6 % vs 60.4 %). Under MNAR all methods are biased; MMRM is the least biased (about 8 % of the effect), which shows why sensitivity analyses for MNAR are needed.

**A misleading result.** Under MCAR, LOCF has the highest power (71.0 %) despite its bias, because it treats carried-forward values as real observations and so understates its uncertainty.

![Bias](results/bias.png)
![Coverage](results/coverage.png)
![Type I error](results/type1.png)

## Extension: about half of the patients drop out

`extensions/high_dropout.R` raises the dropout rate (`alpha_drop = -1.5`) so that about 47 % of patients are missing at the last visit, with true effect 0.4 and 500 new replicates per mechanism. LOCF's bias roughly doubles to about 24 % of the true effect under all three mechanisms, and under MCAR its 95 % CI covers the truth in only 87.8 % of trials. Complete-case analysis, unbiased within Monte Carlo error at 30 % dropout, becomes biased under MAR (-6 %, 2.7 MCSE). MMRM stays unbiased under MAR (bias -0.0004) and keeps the highest power under MCAR and MAR. Results: `results/high_dropout.csv`.

## Limitations and next steps

- MMRM degrees of freedom use the residual df from `gls`; the `mmrm` package offers Kenward-Roger and Satterthwaite
- Only one MNAR mechanism; reference-based imputation (e.g. jump-to-reference with `rbmi`) and tipping-point sensitivity analyses are the natural next step
- Fixed sample size and effect size; a larger grid would show how results scale
- Cross-check of the MMRM estimates against SAS PROC MIXED is planned

## References

- ICH E9(R1) Addendum on estimands and sensitivity analysis in clinical trials (2019)
- Morris TP, White IR, Crowther MJ (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine* 38:2074-2102
- National Research Council (2010). *The Prevention and Treatment of Missing Data in Clinical Trials*
