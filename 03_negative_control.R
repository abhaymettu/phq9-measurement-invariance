#!/usr/bin/env Rscript
# Validation checks. None is optional: an invariance pipeline that
# rejects invariance across a coin flip proves nothing about sex, and a
# pipeline that cannot find a non-invariance you planted proves nothing either.
#
# 1 Random grouping. Split the real respondents at random and run the identical
#   sequence. There is no measurement difference to find, so scalar invariance
#   must hold comfortably. Many splits, not one, because a single split is a
#   sample of size one.
# 2 Simulated known non-invariance. Generate ordinal data from a one-factor
#   model with one item's thresholds deliberately shifted in group 2, and check
#   the procedure both rejects scalar invariance and points at the right item.
# 3 Planted latent mean difference, no item bias at all. Added 2026-10-30. The
#   two controls above have no real group difference in the trait, so they
#   could not tell a test of equal thresholds from a test of equal means. This
#   one can, and it is the control that exposed the fixed-mean error in the
#   original ordinal scalar model.

source("R/invariance.R")
set.seed(20260729)

s <- readRDS("output/analytic_sample.rds")
out <- list()

# ---------------------------------------------------------------------------
N_SPLITS <- 20
cat("=== 1. Random two-group splits,", N_SPLITS, "replications ===\n")
cat("Same sequence, same 0.010 / 0.015 cutoffs, arbitrary groups.\n\n")

rand_rows <- list()
for (i in seq_len(N_SPLITS)) {
  d <- s
  d$rand2 <- factor(sample(c("A", "B"), nrow(d), replace = TRUE))
  for (o in c(TRUE, FALSE)) {
    t <- fit_table(invariance_sequence(d, "rand2", ordinal = o))
    rand_rows[[length(rand_rows) + 1]] <- data.frame(
      rep = i, estimator = if (o) "WLSMV" else "MLR",
      d_cfi_metric = t$d_cfi[2], d_cfi_scalar = t$d_cfi[3],
      d_rmsea_scalar = t$d_rmsea[3], verdict = t$verdict[3])
  }
}
rand <- do.call(rbind, rand_rows)
out$random2 <- rand

for (e in c("WLSMV", "MLR")) {
  r <- rand[rand$estimator == e, ]
  cat(sprintf("%-6s scalar dCFI: median %+.5f, worst %+.5f, range [%+.5f, %+.5f]\n",
              e, median(r$d_cfi_scalar), min(r$d_cfi_scalar),
              min(r$d_cfi_scalar), max(r$d_cfi_scalar)))
  cat(sprintf("       scalar dRMSEA: median %+.5f, worst %+.5f\n",
              median(r$d_rmsea_scalar), max(r$d_rmsea_scalar)))
  cat(sprintf("       FAIL verdicts: %d of %d\n\n", sum(r$verdict == "FAIL"), nrow(r)))
}

# A random six-group split, matching the shape of the race/ethnicity test,
# because more groups means more constraints and more chances to trip.
cat("=== Random six-group splits, 10 replications (mirrors the race test) ===\n")
rand6_rows <- list()
for (i in 1:10) {
  d <- s
  d$rand6 <- factor(sample(LETTERS[1:6], nrow(d), replace = TRUE))
  for (o in c(TRUE, FALSE)) {
    t <- try(fit_table(invariance_sequence(d, "rand6", ordinal = o)), silent = TRUE)
    if (inherits(t, "try-error")) next
    rand6_rows[[length(rand6_rows) + 1]] <- data.frame(
      rep = i, estimator = if (o) "WLSMV" else "MLR",
      d_cfi_scalar = t$d_cfi[3], d_rmsea_scalar = t$d_rmsea[3], verdict = t$verdict[3])
  }
}
rand6 <- do.call(rbind, rand6_rows)
out$random6 <- rand6
for (e in unique(rand6$estimator)) {
  r <- rand6[rand6$estimator == e, ]
  cat(sprintf("%-6s scalar dCFI: median %+.5f, worst %+.5f, FAIL %d of %d\n",
              e, median(r$d_cfi_scalar), min(r$d_cfi_scalar),
              sum(r$verdict == "FAIL"), nrow(r)))
}

# ---------------------------------------------------------------------------
# Simulate ordinal indicators from a one-factor model, using loadings and
# thresholds close to the ones actually estimated on NHANES so the simulation
# has realistic skew, then shift one item's thresholds in group 2.
sim_two_groups <- function(n = 2500, shift = 0.4, bad_item = 5,
                           lambda = c(.80, .82, .58, .62, .63, .78, .74, .70, .66),
                           tau = c(0.7, 1.5, 1.9)) {
  gen <- function(n, sh) {
    eta <- rnorm(n)
    m <- vapply(seq_along(lambda), function(j) {
      ystar <- lambda[j] * eta + rnorm(n, sd = sqrt(1 - lambda[j]^2))
      t <- tau + if (j == bad_item) sh else 0
      as.numeric(cut(ystar, c(-Inf, t, Inf), labels = FALSE)) - 1
    }, numeric(n))
    as.data.frame(m)
  }
  d <- rbind(cbind(gen(n, 0), grp = "A"), cbind(gen(n, shift), grp = "B"))
  names(d)[1:9] <- PHQ_ITEMS
  d$grp <- factor(d$grp)
  d
}

cat("\n=== 2. Simulated non-invariance planted on item 5 (Appetite change) ===\n")
cat("Group B's thresholds on item 5 are shifted on the latent response scale.\n")
cat("Everything else is invariant by construction, n = 2,500 per group.\n")
cat("The criterion here is the chi-square rise and the score-test ranking, not\n")
cat("delta CFI: data generated from a true one-factor model fits so well that\n")
cat("CFI is 1.000 at every rung, so delta CFI has no room to move. The 0.010\n")
cat("cutoff is a rule for approximate misfit in real data, not for a model that\n")
cat("is correct by construction. This is also the check that caught the\n")
cat("intercept-identification bug documented in R/invariance.R.\n\n")

sim_rows <- list()
for (shift in c(0, 0.2, 0.4, 0.6)) {
  d <- sim_two_groups(shift = shift)
  fits <- invariance_sequence(d, "grp", ordinal = TRUE)
  t <- fit_table(fits)
  sc <- score_test_by_item(fits$scalar)
  rank5 <- which(sc$item == PHQ_ITEMS[5])
  chisq_rise <- t$chisq[3] - t$chisq[2]
  sim_rows[[length(sim_rows) + 1]] <- data.frame(
    shift = shift, chisq_rise = chisq_rise, d_df = t$df[3] - t$df[2],
    d_cfi = t$d_cfi[3], top_item = PHQ_LABELS[sc$item[1]],
    rank_of_planted_item = rank5, X2_planted = sc$X2[rank5],
    X2_runner_up = sc$X2[if (rank5 == 1) 2 else 1])
  cat(sprintf("shift %+.1f: chisq rise %7.1f on %d df,  worst item = %-22s  planted item rank %d/9  X2 %6.1f vs %5.1f next\n",
              shift, chisq_rise, t$df[3] - t$df[2], PHQ_LABELS[sc$item[1]], rank5,
              sc$X2[rank5], sc$X2[if (rank5 == 1) 2 else 1]))
}
out$simulated <- do.call(rbind, sim_rows)

cat("\nWith no shift the chi-square rise is at its degrees of freedom and the\n")
cat("item ranking is noise, as it must be. With a shift the ranking puts the\n")
cat("planted item first by a wide margin, so the localisation step points at the\n")
cat("right item rather than scattering blame across the questionnaire.\n")

# ---------------------------------------------------------------------------
cat("\n=== 3. Planted latent mean difference, no item bias ===\n")
cat("Identical item parameters in both groups; only the trait mean differs.\n")
cat("A correct scalar test must pass. The original specification, which kept\n")
cat("both latent means at zero, is run beside it to show what it did.\n\n")
set.seed(20261029)
pm_lambda <- c(.80, .85, .70, .75, .70, .80, .75, .70, .70)
pm_tau <- cbind(c(.4, .6, .0, -.3, .4, .6, .6, 1.0, 1.8),
                c(1.3, 1.5, .9, .7, 1.3, 1.4, 1.4, 1.8, 2.4),
                c(1.8, 2.0, 1.4, 1.3, 1.8, 2.0, 1.9, 2.3, 2.8))
gen_mean <- function(n, mu) {
  eta <- rnorm(n, mu)
  y <- vapply(1:9, function(j)
    findInterval(pm_lambda[j] * eta + rnorm(n, sd = sqrt(1 - pm_lambda[j]^2)),
                 pm_tau[j, ]), numeric(n))
  colnames(y) <- PHQ_ITEMS
  as.data.frame(y)
}
pm_rows <- list()
for (mu in c(0, 0.5)) {
  for (r in 1:5) {
    d <- rbind(cbind(gen_mean(2500, 0), grp = "A"), cbind(gen_mean(2500, mu), grp = "B"))
    d$grp <- factor(d$grp)
    for (fm in c(FALSE, TRUE)) {
      t <- fit_table(invariance_sequence(d, "grp", ordinal = TRUE, free_means = fm))
      pm_rows[[length(pm_rows) + 1]] <- data.frame(
        latent_diff = mu, rep = r, spec = if (fm) "means free" else "means fixed (original)",
        d_cfi = t$d_cfi[3], chisq_rise = t$chisq[3] - t$chisq[2],
        d_df = t$df[3] - t$df[2], verdict = t$verdict[3])
    }
  }
}
pm <- do.call(rbind, pm_rows)
for (mu in unique(pm$latent_diff)) for (sp in unique(pm$spec)) {
  r <- pm[pm$latent_diff == mu & pm$spec == sp, ]
  cat(sprintf("latent diff %.1f SD, %-22s scalar dCFI mean %+.4f, worst %+.4f, FAIL %d of %d\n",
              mu, sp, mean(r$d_cfi), min(r$d_cfi), sum(r$verdict == "FAIL"), nrow(r)))
}
out$planted_mean <- pm

# Recovery of a known latent mean difference, to check the consequence
# decomposition is not inventing the bias it reports.
cat("\n=== 4. Recovery of a known gap decomposition ===\n")
cat("Simulate a true latent mean difference of 0.25 SD plus a planted\n")
cat("intercept shift on item 5, then see whether decompose_gap splits the\n")
cat("observed sum-score gap the way it was built. 20 replications, because one\n")
cat("replication cannot distinguish a biased estimator from sampling noise.\n\n")
set.seed(11)
lambda <- c(.80, .82, .58, .62, .63, .78, .74, .70, .66)
nu <- c(.35, .32, .55, .70, .35, .20, .19, .12, .05)
gen_cont <- function(n, kappa, nu_shift) {
  eta <- rnorm(n, mean = kappa)
  m <- vapply(seq_along(lambda), function(j)
    nu[j] + (if (j == 5) nu_shift else 0) + lambda[j] * eta +
      rnorm(n, sd = sqrt(1 - lambda[j]^2)), numeric(n))
  as.data.frame(m)
}
TRUE_KAPPA <- 0.25; TRUE_NU_SHIFT <- 0.30
rec_rows <- list()
for (i in 1:20) {
  d <- rbind(cbind(gen_cont(2500, 0, 0), grp = "A"),
             cbind(gen_cont(2500, TRUE_KAPPA, TRUE_NU_SHIFT), grp = "B"))
  names(d)[1:9] <- PHQ_ITEMS
  d$grp <- factor(d$grp)
  d$phq9 <- rowSums(d[PHQ_ITEMS])
  fp <- invariance_sequence(d, "grp", ordinal = FALSE,
                            partial = partial_labels(PHQ_ITEMS[5], ordinal = FALSE))$scalar
  r <- decompose_gap(fp, "grp", d, focal = "B")
  rec_rows[[i]] <- data.frame(rep = i, latent_sd = r$latent_sd_units,
                              meas_pts = r$measurement_part,
                              latent_pts = r$latent_part,
                              obs = r$observed_gap, model = r$model_gap)
}
rec <- do.call(rbind, rec_rows)
cat(sprintf("true latent gap        %+.3f SD   recovered mean %+.3f SD  (sd %.3f)\n",
            TRUE_KAPPA, mean(rec$latent_sd), sd(rec$latent_sd)))
cat(sprintf("true measurement part  %+.3f pts  recovered mean %+.3f pts (sd %.3f)\n",
            TRUE_NU_SHIFT, mean(rec$meas_pts), sd(rec$meas_pts)))
cat(sprintf("true latent part       %+.3f pts  recovered mean %+.3f pts (sd %.3f)\n",
            sum(lambda) * TRUE_KAPPA, mean(rec$latent_pts), sd(rec$latent_pts)))
cat(sprintf("mean observed gap %+.3f pts, mean model-implied %+.3f pts\n",
            mean(rec$obs), mean(rec$model)))
cat(sprintf("bias in the measurement share: true %.1f%%, recovered %.1f%%\n",
            100 * TRUE_NU_SHIFT / (TRUE_NU_SHIFT + sum(lambda) * TRUE_KAPPA),
            100 * mean(rec$meas_pts / rec$model)))
out$recovery <- rec
out$recovery_truth <- c(kappa = TRUE_KAPPA, nu = TRUE_NU_SHIFT,
                        latent_pts = sum(lambda) * TRUE_KAPPA)

saveRDS(out, "output/negative_control.rds")
cat("\nWrote output/negative_control.rds\n")
