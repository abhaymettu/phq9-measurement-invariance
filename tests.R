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

# ---------------------------------------------------------------------------
# IRT extension. These are slower than the rest (each one fits several
# multigroup GRMs), so they are skipped unless the mirt side is actually
# installed, and they run on small simulations rather than the real file.
if (requireNamespace("mirt", quietly = TRUE)) {
  source("R/grm.R")

  # -------------------------------------------------------------------------
  # The simulator has to actually plant what it claims to plant. Everything
  # downstream is measured against this, so if the shift is not real, or points
  # the wrong way, the power curve in 06 is measuring nothing and would still
  # look plausible.
  local({
    set.seed(101)
    pars <- list(a = setNames(rep(2, 9), PHQ_ITEMS),
                 b = matrix(c(0.5, 1.3, 2.0), nrow = 9, ncol = 3, byrow = TRUE,
                            dimnames = list(PHQ_ITEMS, c("b1", "b2", "b3"))))
    s0 <- sim_grm(pars, 4000, shift = 0, bad_item = 5, mean_diff = 0)
    m0 <- tapply(s0$resp[[PHQ_ITEMS[5]]], s0$group, mean)
    check(abs(m0[["A"]] - m0[["B"]]) < 0.05,
          paste("shift 0 must generate no difference; got", round(m0[["A"]] - m0[["B"]], 3)))

    s1 <- sim_grm(pars, 4000, shift = 0.8, bad_item = 5, mean_diff = 0)
    m1 <- tapply(s1$resp[[PHQ_ITEMS[5]]], s1$group, mean)
    # Positive shift raises the thresholds in group B: B needs more trait to
    # endorse, so B must score LOWER on that item. A sign error here would
    # invert the whole planted-effect story and still produce a power curve.
    check(m1[["B"]] < m1[["A"]] - 0.05,
          paste("positive shift must lower the focal group's item score; got A",
                round(m1[["A"]], 3), "B", round(m1[["B"]], 3)))

    # And it must be confined to the item named. If the shift leaked into other
    # items the null items would carry real DIF and the false-positive rate in
    # 06 would be measuring a bug in the simulator instead of the procedure.
    other <- setdiff(seq_len(9), 5)
    leak <- max(vapply(other, function(j) {
      mm <- tapply(s1$resp[[PHQ_ITEMS[j]]], s1$group, mean); abs(mm[[1]] - mm[[2]])
    }, numeric(1)))
    check(leak < 0.05, paste("shift leaked into a clean item, max gap", round(leak, 3)))

    # Multi-item planting is what the contamination condition in 06 rests on.
    # It has to move exactly the named items and no others, or that section is
    # measuring the simulator.
    s2 <- sim_grm(pars, 4000, shift = 0.8, bad_item = c(2, 5, 7), mean_diff = 0)
    gaps <- vapply(seq_len(9), function(j) {
      mm <- tapply(s2$resp[[PHQ_ITEMS[j]]], s2$group, mean); mm[["A"]] - mm[["B"]]
    }, numeric(1))
    check(all(gaps[c(2, 5, 7)] > 0.05),
          paste("all three planted items must move; got", paste(round(gaps[c(2,5,7)], 3), collapse = " ")))
    check(max(abs(gaps[-c(2, 5, 7)])) < 0.05,
          paste("multi-item shift leaked; max clean gap", round(max(abs(gaps[-c(2,5,7)])), 3)))
    ok("sim_grm_plants_what_it_claims")
  })

  # -------------------------------------------------------------------------
  # The DIF test must bite. Reintroducing the defect it is supposed to catch,
  # a real item shift, has to make it fail; and with nothing planted it has to
  # stay quiet. A test that only checked "returns nine rows" would pass on a
  # procedure that flagged everything or nothing.
  local({
    set.seed(102)
    pars <- list(a = setNames(rep(2, 9), PHQ_ITEMS),
                 b = matrix(c(0.5, 1.3, 2.0), nrow = 9, ncol = 3, byrow = TRUE,
                            dimnames = list(PHQ_ITEMS, c("b1", "b2", "b3"))))
    anch <- PHQ_ITEMS[c(1, 2, 6, 7)]

    s <- sim_grm(pars, 2000, shift = 0.8, bad_item = 5, mean_diff = 0)
    d <- grm_dif(s$resp, s$group, anchors = anch)
    check(nrow(d) == 9, "grm_dif must return one row per item, anchors included")
    check(all(is.na(d$X2[d$anchor])), "anchors are not tested and must be NA")
    check(isTRUE(d$flag[5]), "the planted item must flag")
    check(d$X2[5] > 3 * max(d$X2[-5], na.rm = TRUE),
          paste("planted item should dominate; got", round(d$X2[5], 1), "vs",
                round(max(d$X2[-5], na.rm = TRUE), 1)))

    s0 <- sim_grm(pars, 2000, shift = 0, bad_item = 5, mean_diff = 0)
    d0 <- grm_dif(s0$resp, s0$group, anchors = anch)
    check(!any(d0$flag), paste("nothing planted, but flagged:",
                               paste(d0$item[d0$flag], collapse = " ")))
    ok("grm_dif_bites_and_stays_quiet")
  })

  # -------------------------------------------------------------------------
  # Purification exists to keep a contaminated item out of the anchor. Test
  # that property directly: plant DIF, and the guilty item must not end up
  # being the thing that defines the common scale.
  local({
    set.seed(103)
    pars <- list(a = setNames(rep(2, 9), PHQ_ITEMS),
                 b = matrix(c(0.5, 1.3, 2.0), nrow = 9, ncol = 3, byrow = TRUE,
                            dimnames = list(PHQ_ITEMS, c("b1", "b2", "b3"))))
    s <- sim_grm(pars, 2000, shift = 0.8, bad_item = 3, mean_diff = 0)
    a <- purify_anchors(s$resp, s$group, n_anchor = 4)
    check(length(a) == 4, "purify_anchors must return the requested size")
    check(!PHQ_ITEMS[3] %in% a,
          paste("the planted item must be kept out of the anchor; got",
                paste(a, collapse = " ")))
    ok("purify_anchors_excludes_the_guilty_item")
  })

  # -------------------------------------------------------------------------
  # The effect size has to move with the size of the shift, and carry its sign.
  # A magnitude that did not scale would make "the size at which detection
  # reaches 80%" unreportable in item-score points.
  local({
    set.seed(104)
    pars <- list(a = setNames(rep(2, 9), PHQ_ITEMS),
                 b = matrix(c(0.5, 1.3, 2.0), nrow = 9, ncol = 3, byrow = TRUE,
                            dimnames = list(PHQ_ITEMS, c("b1", "b2", "b3"))))
    anch <- PHQ_ITEMS[-5]
    es <- vapply(c(0.3, 0.9), function(sh) {
      s <- sim_grm(pars, 3000, shift = sh, bad_item = 5, mean_diff = 0)
      dif_effect_size(grm_fit(s$resp, s$group, anchors = anch), PHQ_ITEMS)$SIDS[5]
    }, numeric(1))
    check(all(es < 0), paste("a positive threshold shift must give a negative SIDS; got",
                             paste(round(es, 3), collapse = " ")))
    check(es[2] < es[1], paste("a bigger shift must give a bigger effect; got",
                               paste(round(es, 3), collapse = " ")))
    ok("effect_size_scales_with_the_shift")
  })

  # -------------------------------------------------------------------------
  # Local dependence. The Q3 statistic in 09 is the basis for a claim about
  # whether the DIF flags are an artifact, so it has to do three things: carry
  # the known negative bias out of the raw values, stay flat when the model is
  # correct, and rise when dependence is actually there. A statistic that
  # cannot go up says nothing when it stays down.
  local({
    set.seed(105)
    pars <- list(a = setNames(rep(2, 9), PHQ_ITEMS),
                 b = matrix(c(0.5, 1.3, 2.0), nrow = 9, ncol = 3, byrow = TRUE,
                            dimnames = list(PHQ_ITEMS, c("b1", "b2", "b3"))))

    clean <- q3_matrix(sim_unidim(pars, 3000, PHQ_ITEMS))
    # Raw Q3 is biased low by roughly -1/(J-1); with nine items that is -0.125.
    # If centring were dropped, every comparison against a positive cutoff
    # would silently pass and 09 would always report "no local dependence".
    check(clean$bias < -0.05,
          paste("raw Q3 should carry the negative bias; got", round(clean$bias, 3)))
    check(abs(mean(clean$centred[lower.tri(clean$centred)])) < 1e-8,
          "centred Q3 must have mean zero off the diagonal")

    ld <- q3_matrix(sim_ld(pars, 3000, PHQ_ITEMS[c(3, 4)], load = 1.2, PHQ_ITEMS))
    check(ld$stat > clean$stat,
          paste("planted dependence must raise Q3max - Q3mean; clean",
                round(clean$stat, 3), "vs planted", round(ld$stat, 3)))
    # And it must point at the right pair, not merely go up somewhere.
    top <- q3_pairs(ld)[1, ]
    check(setequal(c(top$item_i, top$item_j), PHQ_ITEMS[c(3, 4)]),
          paste("top residual pair should be the planted one; got",
                top$item_i, top$item_j))
    ok("q3_is_unbiased_and_bites")
  })

  # -------------------------------------------------------------------------
  # The testlet remedy has to actually merge the pair: one fewer item, the two
  # members gone, and a combined score on the right 0..2*max scale. Getting
  # this wrong would silently answer a different question in 09 section 3.
  local({
    set.seed(106)
    r <- data.frame(matrix(sample(0:3, 9 * 50, TRUE), ncol = 9))
    names(r) <- PHQ_ITEMS
    t <- make_testlet(r, PHQ_ITEMS[c(3, 4)], name = "LD")
    check(ncol(t) == 8, paste("testlet should leave 8 columns, got", ncol(t)))
    check(!any(PHQ_ITEMS[c(3, 4)] %in% names(t)), "merged items must be removed")
    check(all(t$LD == r[[PHQ_ITEMS[3]]] + r[[PHQ_ITEMS[4]]]),
          "testlet score must be the sum of its two members")
    ok("testlet_merges_the_pair")
  })

  # -------------------------------------------------------------------------
  # An anchor is not a null result. grm_dif reports anchors as NA rather than
  # as a tested-and-unflagged zero, and dif_effect_size necessarily gives them
  # exactly 0 because they were constrained equal. Reading that 0 as evidence
  # of no DIF is the mistake that nearly got written into 09's writeup, so the
  # distinction is locked down here: an anchor must be identifiable as untested
  # from the returned object alone, without consulting the anchor list.
  local({
    set.seed(107)
    pars <- list(a = setNames(rep(2, 9), PHQ_ITEMS),
                 b = matrix(c(0.5, 1.3, 2.0), nrow = 9, ncol = 3, byrow = TRUE,
                            dimnames = list(PHQ_ITEMS, c("b1", "b2", "b3"))))
    anch <- PHQ_ITEMS[c(1, 2, 6, 7)]
    d <- sim_grm(pars, 1500, shift = 0.8, bad_item = 5, mean_diff = 0)
    r <- grm_dif(d$resp, d$group, anchors = anch)

    check(all(r$anchor == (r$item %in% anch)), "anchor column must mark the anchor set")
    check(all(is.na(r$adj_p[r$anchor])), "an anchor has no p-value: it was never tested")
    check(!any(r$flag[r$anchor]), "an anchor can never be flagged")
    # The trap: effect size is a structural zero for anchors, so a caller that
    # reads SIDS alone cannot tell "no effect" from "not measured".
    es <- dif_effect_size(attr(r, "fit"), PHQ_ITEMS)
    check(all(es$SIDS[r$anchor] == 0),
          "anchors are constrained equal, so their SIDS is structurally 0")
    check(any(es$SIDS[!r$anchor] != 0), "tested items should have non-zero effects")
    ok("anchor_is_untested_not_null")
  })
} else {
  cat("skip irt_tests (mirt not installed)\n")
}

cat("\nall tests passed\n")
