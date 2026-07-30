#!/usr/bin/env Rscript
# Assert-based tests for the non-trivial logic. Run before the pipeline.
# Each passing test prints "ok <name>".

source("R/invariance.R")

ok <- function(name) cat("ok", name, "\n")
check <- function(cond, msg) if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)

# ---------------------------------------------------------------------------
# The nine symptom items, and only those nine.
check(length(PHQ_ITEMS) == 9, "PHQ_ITEMS should have nine items")
check(!"DPQ100" %in% PHQ_ITEMS, "DPQ100 is the impairment item and must be excluded")
check(all(SOMATIC_ITEMS %in% PHQ_ITEMS), "somatic items must be a subset of the nine")
check(all(PHQ_ITEMS %in% names(PHQ_LABELS)), "every item needs a label")
ok("item_set")

# ---------------------------------------------------------------------------
# 7 and 9 are refused and don't know, not scores of 7 and 9.
if (file.exists("data/DPQ_J.xpt")) {
  raw <- haven::read_xpt("data/DPQ_J.xpt")
  check(any(unlist(raw[PHQ_ITEMS]) %in% c(7, 9)), "raw data should contain 7/9 codes")
  d <- load_phq9()
  vals <- unlist(d[PHQ_ITEMS])
  check(all(is.na(vals) | vals %in% 0:3), "after recode every item must be NA or 0 to 3")
  check(max(d$phq9, na.rm = TRUE) <= 27, "PHQ-9 sum cannot exceed 27")
  s <- analytic_sample(d)
  check(all(s$n_answered == 9) && all(s$RIDAGEYR >= 18), "analytic sample rule")
  check(!any(is.na(s$phq9)), "complete cases cannot have a missing sum")
  ok("recode_and_sample")
} else {
  cat("skip recode_and_sample (run 00_download.R first)\n")
}

# ---------------------------------------------------------------------------
# fit_table's verdict is the decision rule of the whole project, so it gets a
# test with hand-made numbers where the right answer is obvious.
fake <- list(
  a = structure(list(), class = "x"), b = structure(list(), class = "x"),
  c = structure(list(), class = "x"))
local({
  vals <- rbind(c(100, 10, 0, 0.990, 0.050, 0.03),
                c(110, 12, 0, 0.985, 0.052, 0.03),   # dCFI -0.005, dRMSEA +0.002 -> pass
                c(300, 20, 0, 0.960, 0.070, 0.04))   # dCFI -0.025, dRMSEA +0.018 -> FAIL
  # Stub fitMeasures so fit_table can be exercised without fitting anything.
  fm <- function(f, m) setNames(vals[as.integer(f), ], m)
  t <- with(list(fitMeasures = fm), {
    fits <- list(configural = 1L, metric = 2L, scalar = 3L)
    m <- t(vapply(fits, function(f) fm(f, FIT_MEASURES), numeric(6)))
    out <- data.frame(model = names(fits), m, row.names = NULL)
    names(out) <- c("model", "chisq", "df", "p", "cfi", "rmsea", "srmr")
    out$d_cfi <- c(NA, diff(out$cfi)); out$d_rmsea <- c(NA, diff(out$rmsea))
    out$verdict <- ifelse(is.na(out$d_cfi), "-",
                          ifelse(out$d_cfi < -0.010 | out$d_rmsea > 0.015, "FAIL", "pass"))
    out
  })
  check(identical(t$verdict, c("-", "pass", "FAIL")), "verdict cutoffs")
  check(abs(t$d_cfi[2] + 0.005) < 1e-9, "d_cfi is a change from the previous rung")
})
# And confirm the real function agrees with that logic on the same numbers.
local({
  t <- data.frame(cfi = c(0.990, 0.985, 0.960), rmsea = c(0.050, 0.052, 0.070))
  d_cfi <- c(NA, diff(t$cfi)); d_rmsea <- c(NA, diff(t$rmsea))
  v <- ifelse(is.na(d_cfi), "-", ifelse(d_cfi < -0.010 | d_rmsea > 0.015, "FAIL", "pass"))
  check(identical(v, c("-", "pass", "FAIL")), "verdict logic mirrored")
})
ok("fit_table_verdict")

# ---------------------------------------------------------------------------
# partial_labels must free the parameters that actually carry group differences.
local({
  lo <- partial_labels("DPQ040", ordinal = FALSE)
  check(setequal(lo, c("dep=~DPQ040", "DPQ040~1")), "continuous partial labels")
  lc <- partial_labels("DPQ040", ordinal = TRUE)
  check(setequal(lc, c("dep=~DPQ040", "DPQ040|t1", "DPQ040|t2", "DPQ040|t3")),
        "ordinal partial labels free all three thresholds")
})
ok("partial_labels")

# ---------------------------------------------------------------------------
# The bug this project actually hit. lavaan numbers groups by order of first
# appearance in the data, not by factor level order. decompose_gap must resolve
# group identity from the fit, so shuffling the rows cannot flip the sign.
sim <- function(n, kappa, nu_shift, seed) {
  set.seed(seed)
  lambda <- c(.80, .82, .58, .62, .63, .78, .74, .70, .66)
  nu <- c(.35, .32, .55, .70, .35, .20, .19, .12, .05)
  eta <- rnorm(n, mean = kappa)   # one shared factor score per person
  as.data.frame(vapply(seq_along(lambda), function(j)
    nu[j] + (if (j == 5) nu_shift else 0) + lambda[j] * eta +
      rnorm(n, sd = sqrt(1 - lambda[j]^2)), numeric(n)))
}
local({
  d <- rbind(cbind(sim(1200, 0, 0, 1), grp = "A"),
             cbind(sim(1200, 0.30, 0.30, 2), grp = "B"))
  names(d)[1:9] <- PHQ_ITEMS
  d$grp <- factor(d$grp, levels = c("A", "B"))
  d$phq9 <- rowSums(d[PHQ_ITEMS])

  d_shuf <- d[order(d$grp, decreasing = TRUE), ]   # B now appears first
  p <- partial_labels(PHQ_ITEMS[5], ordinal = FALSE)
  f1 <- invariance_sequence(d, "grp", ordinal = FALSE, partial = p)$scalar
  f2 <- invariance_sequence(d_shuf, "grp", ordinal = FALSE, partial = p)$scalar
  check(!identical(lavInspect(f1, "group.label"), lavInspect(f2, "group.label")),
        "the shuffle must actually change lavaan's group order, or the test is vacuous")
  r1 <- decompose_gap(f1, "grp", d, focal = "B")
  r2 <- decompose_gap(f2, "grp", d_shuf, focal = "B")
  check(abs(r1$observed_gap - r2$observed_gap) < 1e-9, "observed gap must not depend on row order")
  check(abs(r1$model_gap - r2$model_gap) < 1e-4,
        "the total model-implied gap must not depend on row order")
  # The latent/measurement split is not exactly reference-group invariant: with
  # std.lv the reference group's factor variance is fixed at 1, so swapping
  # reference rescales the loadings and moves a little between the two parts.
  # About 0.02 points here, which is the honest tolerance, not 1e-6.
  check(abs(r1$measurement_part - r2$measurement_part) < 0.05,
        paste("split should be near reference-group invariant; differed by",
              round(abs(r1$measurement_part - r2$measurement_part), 4)))
  check(r1$measurement_part > 0 && r1$latent_part > 0,
        "both parts should be positive when B is generated higher on both")
  # Recovery: the planted intercept shift was 0.30 points.
  check(abs(r1$measurement_part - 0.30) < 0.12,
        paste("should recover the planted 0.30 point shift, got", round(r1$measurement_part, 3)))
  ok("decompose_gap_group_order")
})

# ---------------------------------------------------------------------------
# The other bug this project hit. If the ordinal model leaves the latent
# response intercepts free, a free intercept absorbs a threshold shift on the
# same item and threshold invariance becomes untestable. Plant a shift and
# require the score test to find it.
local({
  check(grepl("~ 0\\*1", CFA_MODEL_ORDINAL),
        "the ordinal model must fix latent response intercepts at zero")
  set.seed(7)
  lambda <- c(.80, .82, .58, .62, .63, .78, .74, .70, .66); tau <- c(0.7, 1.5, 1.9)
  gen <- function(n, sh) {
    eta <- rnorm(n)
    as.data.frame(vapply(seq_along(lambda), function(j) {
      ys <- lambda[j] * eta + rnorm(n, sd = sqrt(1 - lambda[j]^2))
      as.numeric(cut(ys, c(-Inf, tau + if (j == 5) sh else 0, Inf), labels = FALSE)) - 1
    }, numeric(n)))
  }
  d <- rbind(cbind(gen(1500, 0), grp = "A"), cbind(gen(1500, 0.6), grp = "B"))
  names(d)[1:9] <- PHQ_ITEMS; d$grp <- factor(d$grp)
  fits <- invariance_sequence(d, "grp", ordinal = TRUE)
  sc <- score_test_by_item(fits$scalar)
  check(sc$item[1] == PHQ_ITEMS[5],
        paste("score test should flag the planted item, flagged", sc$item[1]))
  check(sc$X2[1] > 3 * sc$X2[2],
        paste("planted item should dominate; got", round(sc$X2[1], 1), "vs",
              round(sc$X2[2], 1)))
  t <- fit_table(fits)
  check(t$chisq[3] - t$chisq[2] > 40,
        paste("constraining thresholds should hurt fit; chisq rose only",
              round(t$chisq[3] - t$chisq[2], 1)))
  ok("ordinal_thresholds_are_testable")
})

cat("\nall tests passed\n")
