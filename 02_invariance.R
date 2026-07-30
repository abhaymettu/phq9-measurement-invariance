#!/usr/bin/env Rscript
# The measurement invariance sequence for the PHQ-9, across sex, age band and
# race/ethnicity, under both estimators. Then partial invariance to localise the
# offending items, then the consequence for a naive sum-score comparison.

source("R/invariance.R")
suppressMessages(library(survey))
options(survey.lonely.psu = "adjust")

s <- readRDS("output/analytic_sample.rds")
res <- list()

# ---------------------------------------------------------------------------
cat("=== Baseline one-factor CFA, single group, n =", nrow(s), "===\n")
for (o in c(TRUE, FALSE)) {
  f <- fit_cfa(s, ordinal = o)
  cat(sprintf("%-6s CFI %.4f  RMSEA %.4f  SRMR %.4f  chisq(%d) %.1f\n",
              if (o) "WLSMV" else "MLR",
              fitMeasures(f, "cfi.scaled"), fitMeasures(f, "rmsea.scaled"),
              fitMeasures(f, "srmr"), fitMeasures(f, "df"),
              fitMeasures(f, "chisq.scaled")))
  res$baseline[[if (o) "WLSMV" else "MLR"]] <-
    standardizedSolution(f)[1:9, c("rhs", "est.std")]
}
cat("\nStandardised loadings (WLSMV):\n")
b <- res$baseline$WLSMV
print(data.frame(item = PHQ_LABELS[b$rhs], loading = round(b$est.std, 3)), row.names = FALSE)

# ---------------------------------------------------------------------------
# Invariance sequences. Both estimators for every grouping variable, because
# the whole point is that the verdict can depend on that choice.
groups <- c(sex = "sex", age = "age_grp", race = "race")
for (gname in names(groups)) {
  for (o in c(TRUE, FALSE)) {
    est <- if (o) "WLSMV" else "MLR"
    fits <- invariance_sequence(s, groups[[gname]], ordinal = o)
    tab <- fit_table(fits)
    tab$group <- gname; tab$estimator <- est
    res$seq[[paste(gname, est)]] <- tab
    cat(sprintf("\n=== %s, %s (%d groups) ===\n", gname, est,
                length(lavInspect(fits$configural, "group.label"))))
    print(tab[, c("model", "chisq", "df", "cfi", "rmsea", "srmr", "d_cfi", "d_rmsea", "verdict")],
          row.names = FALSE, digits = 4)
    res$scores[[paste(gname, est)]] <- score_test_by_item(fits$scalar)
  }
}

cat("\nDecision rule applied: scalar invariance is rejected if CFI drops by more\n")
cat("than 0.010 or RMSEA rises by more than 0.015 relative to the metric model\n")
cat("(Cheung and Rensvold 2002; Chen 2007). The scaled chi-square difference\n")
cat("test is not used as the rule: at n = 5,068 it is significant for every\n")
cat("model here, including the random-split negative control, so it carries no\n")
cat("information about practical non-invariance.\n")

cat("\n=== Chi-square difference tests, shown to make that point concrete ===\n")
fits_sex <- invariance_sequence(s, "sex", ordinal = FALSE)
print(anova(fits_sex$configural, fits_sex$metric, fits_sex$scalar))

# ---------------------------------------------------------------------------
# Sex is the case that fails under both estimators, so it is the focus.
freed_by_est <- list()
for (o in c(TRUE, FALSE)) {
  est <- if (o) "WLSMV" else "MLR"
  cat(sprintf("\n=== Localising non-invariance across sex (%s scalar model) ===\n", est))
  sc <- res$scores[[paste("sex", est)]]
  sc$item_label <- PHQ_LABELS[sc$item]
  sc$somatic <- sc$item %in% SOMATIC_ITEMS
  print(sc[, c("item_label", "X2", "somatic")], row.names = FALSE, digits = 4)

  cat("\nReleasing constraints one item at a time, worst score statistic first:\n")
  release <- character(0); path <- list()
  for (it in sc$item) {
    release <- c(release, it)
    t <- fit_table(invariance_sequence(s, "sex", ordinal = o,
                                       partial = partial_labels(release, ordinal = o)))
    path[[length(path) + 1]] <-
      data.frame(estimator = est, k = length(release),
                 freed = paste(PHQ_LABELS[release], collapse = " + "),
                 d_cfi = t$d_cfi[3], d_rmsea = t$d_rmsea[3], verdict = t$verdict[3])
    cat(sprintf("  free %-2d item(s), latest [%s]: dCFI %+.4f  dRMSEA %+.4f  %s\n",
                length(release), PHQ_LABELS[it], t$d_cfi[3], t$d_rmsea[3], t$verdict[3]))
    if (t$verdict[3] == "pass") break
  }
  res$partial_path[[est]] <- do.call(rbind, path)
  freed_by_est[[est]] <- release
  cat("Partial scalar invariance reached after freeing:",
      paste(PHQ_LABELS[release], collapse = ", "), "\n")
}
freed <- freed_by_est$MLR
res$freed <- freed_by_est
cat("\nItems flagged by both estimators:",
    paste(PHQ_LABELS[intersect(freed_by_est$WLSMV, freed_by_est$MLR)], collapse = ", "), "\n")

# ---------------------------------------------------------------------------
cat("\n=== Consequence for a naive female minus male comparison ===\n")
fit_partial <- invariance_sequence(s, "sex", ordinal = FALSE,
                                   partial = partial_labels(freed, ordinal = FALSE))$scalar
fit_full <- fits_sex$scalar
dp <- decompose_gap(fit_partial, "sex", s, focal = "Female")
df <- decompose_gap(fit_full, "sex", s, focal = "Female")
res$gap <- list(partial = dp, full = df)

cat(sprintf("observed gap in mean PHQ-9 sum, unweighted:      %+.3f points\n", dp$observed_gap))
cat(sprintf("model-implied gap, partial scalar model:          %+.3f points\n", dp$model_gap))
cat(sprintf("  attributable to a real latent difference:       %+.3f points\n", dp$latent_part))
cat(sprintf("  attributable to items behaving differently:     %+.3f points\n", dp$measurement_part))
cat(sprintf("  share of the model-implied gap that is measurement: %.1f%%\n",
            100 * dp$measurement_part / dp$model_gap))
cat(sprintf("latent mean difference, partial model:  %+.3f SD\n", dp$latent_sd_units))
cat(sprintf("latent mean difference, forcing full scalar invariance: %+.3f SD\n",
            df$latent_sd_units))
cat(sprintf("forcing invariance overstates the latent sex difference by %.1f%%\n",
            100 * (df$latent_sd_units / dp$latent_sd_units - 1)))

# Robustness: with std.lv the reference group's factor variance is fixed at 1,
# so the split between the latent and measurement parts is not exactly
# reference-group invariant. Refit with the rows reversed, which flips which
# group lavaan calls group 1, and report both.
s_rev <- s[order(s$sex, decreasing = TRUE), ]
dp_rev <- decompose_gap(
  invariance_sequence(s_rev, "sex", ordinal = FALSE,
                      partial = partial_labels(freed, ordinal = FALSE))$scalar,
  "sex", s_rev, focal = "Female")
cat(sprintf("reference-group check: measurement part %+.3f with males as reference,\n",
            dp$measurement_part))
cat(sprintf("  %+.3f with females as reference (share %.1f%% vs %.1f%%)\n",
            dp_rev$measurement_part,
            100 * dp$measurement_part / dp$model_gap,
            100 * dp_rev$measurement_part / dp_rev$model_gap))
res$gap$reversed <- dp_rev

# ---------------------------------------------------------------------------
# Survey-weighted descriptives. The CFA above is unweighted; this is the
# population-level version of the naive comparison that an investigator would
# actually publish, and it is the quantity the measurement bias contaminates.
cat("\n=== Survey-weighted naive comparison (survey package) ===\n")
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
                 data = s, nest = TRUE)
cat("represents", format(round(sum(weights(des))), big.mark = ","), "US adults\n")
wm <- svyby(~phq9, ~sex, des, svymean)
print(wm)
wt_gap <- as.numeric(wm$phq9[wm$sex == "Female"] - wm$phq9[wm$sex == "Male"])
ttest <- svyttest(phq9 ~ sex, des)
cat(sprintf("weighted female minus male gap: %+.3f points, 95%% CI [%.3f, %.3f], p = %.2g\n",
            wt_gap, confint(ttest)[1], confint(ttest)[2], ttest$p.value))
res$weighted <- list(means = wm, gap = wt_gap, ci = confint(ttest), p = ttest$p.value,
                     n = nrow(s), pop = sum(weights(des)))

cat(sprintf("\nApplying the same %.1f%% measurement share to the weighted gap:\n",
            100 * dp$measurement_part / dp$model_gap))
cat(sprintf("  %+.3f of the %+.3f weighted point gap is measurement, not symptoms\n",
            wt_gap * dp$measurement_part / dp$model_gap, wt_gap))

saveRDS(res, "output/invariance_results.rds")
cat("\nWrote output/invariance_results.rds\n")
