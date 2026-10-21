 
library(sandwich)
library(dplyr)
library(ggplot2)
 
set.seed(20261019)
B <- 2000
stopifnot(B >= 1)
 
folder <- "Expanded_Modified_Poisson_Results"
cache <- file.path(folder, paste0("cache_v1_B", B))
figdir <- file.path(folder, paste0("Figures_B", B))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
 
# Simulation scenarios
 
scenarios <- expand.grid(n = c(400, 1200, 3000), baseline = c(0.04, 0.12, 0.16, 0.20), female_share = c(0.20, 0.50, 0.70), true_rrr = c(1.0, 1.5), conf_rr_women = c(1.25, 1.50, 2.00), KEEP.OUT.ATTRS = FALSE)
scenarios$scenario <- seq_len(nrow(scenarios))
scenarios$conf_label <- factor(scenarios$conf_rr_women, levels = c(1.25, 1.50, 2.00), labels = c("Same (RR 1.25)", "Moderate (RR 1.50)", "Strong (RR 2.00)"))
scenarios$true_label <- ifelse(scenarios$true_rrr == 1, "No true effect modification", "True effect modification (RRR 1.5)")
scenarios$maximum_risk <- scenarios$baseline * 1.5 * scenarios$true_rrr * scenarios$conf_rr_women * 1.10
stopifnot(nrow(scenarios) == 216, all(scenarios$maximum_risk < 1))
write.csv(scenarios, file.path(folder, "Table_1_simulation_design_216.csv"), row.names = FALSE)
 
# Generate data
 
make_data <- function(n, baseline, female_share, true_rrr, conf_rr_women) {
  sex <- rbinom(n, 1, female_share)
  confounder <- rbinom(n, 1, 0.40)
  exposure <- rbinom(n, 1, plogis(-1 + 0.8 * confounder + 0.2 * sex))
  conf_rr <- ifelse(sex == 1, conf_rr_women, 1.25)
  risk <- baseline * 1.5^exposure * true_rrr^(exposure * sex) * conf_rr^confounder * 1.1^sex
  if (any(!is.finite(risk)) || any(risk < 0 | risk >= 1)) stop("Invalid event probability")
  data.frame(outcome = rbinom(n, 1, risk), exposure, sex, confounder)
}
 
# Modified Poisson regression
 
fit_mp <- function(formula, data) {
 fit <- suppressWarnings(glm(formula, family = poisson(link = "log"), data = data, control = glm.control(maxit = 50)))
  if (!isTRUE(fit$converged) || fit$rank != length(coef(fit))) stop("Model did not converge")
  if (any(!is.finite(coef(fit)))) stop("Non-finite coefficient")
  V <- sandwich::vcovHC(fit, type = "HC0")
  if (any(!is.finite(V))) stop("Non-finite robust covariance")
  list(beta = coef(fit), V = V)
}
 
extract_term <- function(model, term) {
  if (!(term %in% names(model$beta))) stop("Missing term: ", term)
  c(beta = unname(model$beta[term]), variance = unname(model$V[term, term]))
}
 
wald_ci <- function(beta, variance) {
  if (!is.finite(beta) || !is.finite(variance) || variance <= 0) stop("Invalid Wald variance")
  se <- sqrt(unname(variance))
  beta <- unname(beta)
  c(beta = beta, se = se, low = beta - qnorm(.975) * se, high = beta + qnorm(.975) * se)
}
 
# Compare the three models
 
compare_three <- function(dat) {
  events <- with(dat, tapply(outcome, list(factor(sex, levels = 0:1), factor(exposure, levels = 0:1)), sum))
  if (anyNA(events) || any(events == 0)) stop("No events in a sex/exposure cell")
 
  men <- fit_mp(outcome ~ exposure + confounder, subset(dat, sex == 0))
  women <- fit_mp(outcome ~ exposure + confounder, subset(dat, sex == 1))
  pooled <- fit_mp(outcome ~ exposure * sex + confounder, dat)
  full <- fit_mp(outcome ~ (exposure + confounder) * sex, dat)
 
  m <- extract_term(men, "exposure")
  w <- extract_term(women, "exposure")
  b <- extract_term(pooled, "exposure:sex")
  c <- extract_term(full, "exposure:sex")
 
  A <- wald_ci(w["beta"] - m["beta"], w["variance"] + m["variance"])
  Bmodel <- wald_ci(b["beta"], b["variance"])
  Cmodel <- wald_ci(c["beta"], c["variance"])
 
  data.frame(beta_A = unname(A["beta"]), se_A = unname(A["se"]), low_A = unname(A["low"]), high_A = unname(A["high"]),
             beta_B = unname(Bmodel["beta"]), se_B = unname(Bmodel["se"]), low_B = unname(Bmodel["low"]), high_B = unname(Bmodel["high"]),
             beta_C = unname(Cmodel["beta"]), se_C = unname(Cmodel["se"]), low_C = unname(Cmodel["low"]), high_C = unname(Cmodel["high"]))
}
 
# Run the simulations
 
for (i in seq_len(nrow(scenarios))) {
  file <- file.path(cache, sprintf("scenario_%03d.rds", i))
  if (file.exists(file)) {
    cat(sprintf("Scenario %3d/%d already completed\n", i, nrow(scenarios)))
    next
  }
  set.seed(20261019 + i)
  s <- scenarios[i, ]
  replications <- vector("list", B)
  cat(sprintf("Scenario %3d/%d | N=%d, risk=%.2f, women=%.2f, RRR=%.1f, confRR=%.2f\n", i, nrow(scenarios), s$n, s$baseline, s$female_share, s$true_rrr, s$conf_rr_women))
 
  for (r in seq_len(B)) {
    result <- tryCatch({
      dat <- make_data(s$n, s$baseline, s$female_share, s$true_rrr, s$conf_rr_women)
      compare_three(dat)
    }, error = function(e) structure(list(message = conditionMessage(e)), class = "simulation_error"))
 
    if (inherits(result, "simulation_error")) {
      replications[[r]] <- data.frame(scenario = i, replication = r, failed = 1L, failure_reason = result$message)
    } else {
      replications[[r]] <- data.frame(scenario = i, replication = r, failed = 0L, failure_reason = NA_character_, result)
    }
  }
  saveRDS(dplyr::bind_rows(replications), file)
}
 
# Combine the results
 
files <- file.path(cache, sprintf("scenario_%03d.rds", scenarios$scenario))
if (!all(file.exists(files))) stop("Some scenarios are missing from the cache")
runs <- bind_rows(lapply(files, readRDS)) %>% left_join(scenarios, by = "scenario")
write.csv(runs, file.path(folder, "All_simulation_replicates.csv"), row.names = FALSE)
 
valid <- runs %>% filter(failed == 0) %>% mutate(
  diff_BA = beta_B - beta_A, abs_diff_BA = abs(diff_BA),
  diff_CA = beta_C - beta_A, se_diff_CA = se_C - se_A,
  width_A = high_A - low_A, width_B = high_B - low_B, width_C = high_C - low_C,
  endpoint_shift_BA = (abs(low_B - low_A) + abs(high_B - high_A)) / 2,
  overlap_AB = pmax(low_A, low_B) <= pmin(high_A, high_B),
  cover_A = low_A <= log(true_rrr) & high_A >= log(true_rrr),
  cover_B = low_B <= log(true_rrr) & high_B >= log(true_rrr),
  cover_C = low_C <= log(true_rrr) & high_C >= log(true_rrr),
  reject_A = low_A > 0 | high_A < 0, reject_B = low_B > 0 | high_B < 0, reject_C = low_C > 0 | high_C < 0,
  decision_disagree_AB = reject_A != reject_B
)
 
if (nrow(valid) == 0) stop("No successful datasets")
write.csv(valid, file.path(folder, "Successful_fits_all_models.csv"), row.names = FALSE)
write.csv(runs %>% filter(failed == 1), file.path(folder, "Model_failures.csv"), row.names = FALSE)
 
# Summarise results
 
fail_table <- runs %>% group_by(scenario) %>% summarise(
  attempted = n(), successful = sum(failed == 0), failure_pct = 100 * mean(failed), .groups = "drop"
)
 
results_table <- valid %>% group_by(scenario) %>% summarise(
  successful_runs = n(),
  mean_beta_A = mean(beta_A), mean_beta_B = mean(beta_B), mean_beta_C = mean(beta_C),
  bias_A = mean(beta_A - log(first(true_rrr))),
  bias_B = mean(beta_B - log(first(true_rrr))),
  bias_C = mean(beta_C - log(first(true_rrr))),
  mean_se_A = mean(se_A), mean_se_B = mean(se_B), mean_se_C = mean(se_C),
  width_A = mean(width_A), width_B = mean(width_B), width_C = mean(width_C),
  coefficient_B_minus_A = mean(diff_BA), mean_abs_diff_BA = mean(abs_diff_BA),
  mean_endpoint_shift_BA = mean(endpoint_shift_BA),
  overlap_AB_pct = 100 * mean(overlap_AB),
  coverage_A_pct = 100 * mean(cover_A), coverage_B_pct = 100 * mean(cover_B), coverage_C_pct = 100 * mean(cover_C),
  rejection_A_pct = 100 * mean(reject_A), rejection_B_pct = 100 * mean(reject_B), rejection_C_pct = 100 * mean(reject_C),
  decision_disagreement_AB_pct = 100 * mean(decision_disagree_AB),
  max_abs_coefficient_A_C = max(abs(diff_CA)), max_abs_SE_A_C = max(abs(se_diff_CA)),
  .groups = "drop"
)
 
summary_table <- scenarios %>% left_join(fail_table, by = "scenario") %>% left_join(results_table, by = "scenario") %>%
  mutate(coverage_A_mcse_pp = sqrt(coverage_A_pct * (100 - coverage_A_pct) / successful_runs),
         coverage_B_mcse_pp = sqrt(coverage_B_pct * (100 - coverage_B_pct) / successful_runs),
         coverage_C_mcse_pp = sqrt(coverage_C_pct * (100 - coverage_C_pct) / successful_runs))
 
write.csv(summary_table, file.path(folder, "Table_2_all_216_scenarios.csv"), row.names = FALSE)
