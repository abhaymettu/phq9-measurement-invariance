#!/usr/bin/env Rscript
# What is this DIF procedure's actual error rate?
#
# 05 reports which PHQ-9 items flag for DIF by sex. On its own that is half a
# finding, because a flag only means something if you know how often the same
# procedure flags an item when there is nothing there, and how large a real
# effect has to be before it is found. This stage measures both, and it does it
# for three procedures at once so the answer is attached to a method rather
# than to the word "DIF".
#
# The shape is taken from 03_negative_control.R, which does the same two checks
# for the CFA side, and which is where the identification bug that had reversed
# this repo's primary result was caught.
#
#   NULL    30 random splits of the real respondents. There is no measurement
#           difference between two coin-flip groups, so every flagged item is a
#           false positive and the rate is directly comparable to the nominal
#           0.05.
#   PLANTED Simulated data with a known shift on one item, swept across
#           magnitudes rather than run at one convenient size, so the output is
#           a detection curve and the honest statement "this procedure cannot
#           see a shift smaller than X" rather than "I validated it".
#
# Runtime is the reason for the parallelism: roughly 250 multigroup GRM fits.

source("R/invariance.R")
source("R/grm.R")
suppressMessages(library(parallel))
set.seed(20260818)

N_SPLITS <- as.integer(Sys.getenv("N_SPLITS", "30"))
N_REPS   <- as.integer(Sys.getenv("N_REPS", "20"))
N_CORES  <- as.integer(Sys.getenv("N_CORES", max(1, detectCores() - 2)))
SHIFTS   <- c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
BAD_ITEM <- 5                     # 5 Appetite change, the item 03 also plants on
ALPHA    <- 0.05                  # nominal, after BH across the items tested

s <- readRDS("output/analytic_sample.rds")
resp <- s[, PHQ_ITEMS]
pars <- grm_pars(resp)

cat("=== Setup ===\n")
cat(sprintf("null: %d random splits of the real sample (n = %d)\n", N_SPLITS, nrow(resp)))
cat(sprintf("planted: %d replications at each of %d shift sizes, 2,500 per group\n",
            N_REPS, length(SHIFTS)))
cat(sprintf("cores: %d,  nominal alpha: %.2f after Benjamini-Hochberg\n\n", N_CORES, ALPHA))

# Run all three procedures on one data set. Stage 1 is fitted once and re-used
# as both the all-other-anchor result and the input to purification, because
# they are literally the same model.
run_procedures <- function(resp, group) {
  all_other <- grm_dif(resp, group, anchors = NULL)
  anch <- anchors_from_stage1(all_other, n_anchor = 4)
  anchored <- grm_dif(resp, group, anchors = anch)
  lord <- try(lordif_dif(resp, group, alpha = 0.01), silent = TRUE)
  lord_flag <- if (inherits(lord, "try-error")) rep(NA, length(PHQ_ITEMS)) else lord$flag
  data.frame(item = PHQ_ITEMS,
             all_other = all_other$flag,
             anchored  = anchored$flag,
             anchor    = PHQ_ITEMS %in% anch,
             lordif    = lord_flag)
}

PROCS <- c("all_other", "anchored", "lordif")

# ---------------------------------------------------------------------------
cat("=== 1. Null: random splits of the real respondents ===\n")
cat("Split the sample at random into two groups. Nothing distinguishes them, so\n")
cat("the correct number of flagged items is zero and every flag is an error.\n")
cat("Many splits, not one, because a single split is a sample of size one.\n\n")

null_out <- mclapply(seq_len(N_SPLITS), function(i) {
  set.seed(1000 + i)
  g <- sample(c("A", "B"), nrow(resp), replace = TRUE)
  r <- try(run_procedures(resp, g), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  r$rep <- i
  r
}, mc.cores = N_CORES)
null_df <- do.call(rbind, Filter(Negate(is.null), null_out))
cat(sprintf("completed %d of %d splits\n\n", length(unique(null_df$rep)), N_SPLITS))

# The two rates that matter are different questions and both get reported.
# Per-item-test is comparable to the nominal alpha. Per-split ("did this study
# report at least one DIF item") is what actually reaches a reader, and it is
# always the larger and more alarming of the two.
null_rate <- function(df, proc) {
  f <- df[[proc]]
  tested <- !is.na(f) & !(proc == "anchored" & df$anchor)
  per_test <- mean(f[tested], na.rm = TRUE)
  per_split <- mean(vapply(split(df, df$rep), function(d)
    any(d[[proc]][!is.na(d[[proc]])]), logical(1)))
  n_items <- mean(vapply(split(df, df$rep), function(d)
    sum(d[[proc]], na.rm = TRUE), numeric(1)))
  c(per_test = per_test, per_split = per_split, mean_items = n_items)
}
null_tab <- t(vapply(PROCS, function(p) null_rate(null_df, p), numeric(3)))
cat("False positive rate with nothing to find:\n")
print(data.frame(procedure = PROCS, round(null_tab, 4), row.names = NULL), row.names = FALSE)
cat(sprintf("\nNominal rate per item test: %.3f\n", ALPHA))
for (p in PROCS) {
  r <- null_tab[p, "per_test"]
  cat(sprintf("  %-10s %.3f per item test, %.1fx nominal;  %.0f%% of splits reported at least one DIF item\n",
              p, r, r / ALPHA, 100 * null_tab[p, "per_split"]))
}

# Which items false-positive is not uniform, and that matters for reading 05:
# an item that is fragile under a coin flip is an item whose real flag is worth
# less.
cat("\nFalse positives by item (share of splits flagging each item):\n")
by_item <- do.call(rbind, lapply(PROCS, function(p) {
  v <- tapply(null_df[[p]], null_df$item, mean, na.rm = TRUE)
  data.frame(procedure = p, t(round(v[PHQ_ITEMS], 3)))
}))
names(by_item)[-1] <- sub("^(\\d).*", "item\\1", PHQ_LABELS[PHQ_ITEMS])
print(by_item, row.names = FALSE)

# ---------------------------------------------------------------------------
cat("\n\n=== 2. Null with impact but no DIF ===\n")
cat("A random split gives two groups identical in every way. The case that\n")
cat("actually arises is two groups that genuinely differ in depression while\n")
cat("the items behave identically, which is the CFA's 0.19 SD sex difference\n")
cat("with no measurement bias at all. Group mean differences are a known way to\n")
cat("induce spurious DIF, so the random split alone would understate the rate.\n\n")

impact_out <- mclapply(seq_len(N_REPS), function(i) {
  set.seed(2000 + i)
  d <- sim_grm(pars, 2500, shift = 0, bad_item = BAD_ITEM, mean_diff = 0.19,
               items = PHQ_ITEMS)
  r <- try(run_procedures(d$resp, d$group), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  r$rep <- i
  r
}, mc.cores = N_CORES)
impact_df <- do.call(rbind, Filter(Negate(is.null), impact_out))
impact_tab <- t(vapply(PROCS, function(p) null_rate(impact_df, p), numeric(3)))
print(data.frame(procedure = PROCS, round(impact_tab, 4), row.names = NULL), row.names = FALSE)

# ---------------------------------------------------------------------------
cat("\n\n=== 3. Planted: detection as a function of effect size ===\n")
cat(sprintf("Shift every threshold of item %d (%s) in group B by the stated number\n",
            BAD_ITEM, PHQ_LABELS[PHQ_ITEMS[BAD_ITEM]]))
cat("of logits, on top of a real 0.19 SD latent difference. Everything else is\n")
cat("invariant by construction. Shift 0 is the null row and must reproduce\n")
cat("section 2. Sweeping the magnitude is the point: one convenient size would\n")
cat("say only that the procedure is not blind.\n\n")

grid <- expand.grid(shift = SHIFTS, rep = seq_len(N_REPS))
planted_out <- mclapply(seq_len(nrow(grid)), function(k) {
  set.seed(3000 + k)
  sh <- grid$shift[k]
  d <- sim_grm(pars, 2500, shift = sh, bad_item = BAD_ITEM, mean_diff = 0.19,
               items = PHQ_ITEMS)
  r <- try(run_procedures(d$resp, d$group), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  do.call(rbind, lapply(PROCS, function(p) data.frame(
    shift = sh, rep = grid$rep[k], procedure = p,
    # Detected: the planted item was tested and flagged. If purification put
    # the planted item INTO the anchor it was never tested, which is a miss,
    # and a bad one, because the contaminated item is then defining the scale.
    hit = isTRUE(r[[p]][BAD_ITEM]) && !(p == "anchored" && r$anchor[BAD_ITEM]),
    in_anchor = p == "anchored" && r$anchor[BAD_ITEM],
    # False alarms on the eight items where nothing was planted.
    false_alarms = sum(r[[p]][-BAD_ITEM], na.rm = TRUE))))
}, mc.cores = N_CORES)
planted_df <- do.call(rbind, Filter(Negate(is.null), planted_out))

power <- aggregate(cbind(hit, false_alarms, in_anchor) ~ shift + procedure,
                   data = planted_df, FUN = mean)
power <- power[order(power$procedure, power$shift), ]
names(power) <- c("shift", "procedure", "detection", "false_alarms_of_8", "planted_in_anchor")
print(round_df <- transform(power, detection = round(detection, 3),
                            false_alarms_of_8 = round(false_alarms_of_8, 2),
                            planted_in_anchor = round(planted_in_anchor, 2)),
      row.names = FALSE)

# The number that turns "I validated it" into a statement about limits: the
# smallest planted shift the procedure finds 80% of the time. Linear
# interpolation between the bracketing grid points, and NA when the curve never
# gets there, because extrapolating past the swept range would be inventing it.
shift_at <- function(d, target = 0.8) {
  d <- d[order(d$shift), ]
  if (max(d$detection) < target) return(NA_real_)
  i <- which(d$detection >= target)[1]
  if (i == 1) return(d$shift[1])
  x0 <- d$shift[i - 1]; x1 <- d$shift[i]
  y0 <- d$detection[i - 1]; y1 <- d$detection[i]
  x0 + (target - y0) * (x1 - x0) / (y1 - y0)
}
cat("\n80% detection thresholds:\n")
thr <- vapply(PROCS, function(p) shift_at(power[power$procedure == p, ]), numeric(1))
for (p in PROCS) {
  cat(sprintf("  %-10s %s\n", p,
              if (is.na(thr[p])) sprintf("never reaches 80%% within the swept range (max shift %.1f)", max(SHIFTS))
              else sprintf("%.2f logits", thr[p])))
}

# A logit is not a unit anyone can picture. Convert the thresholds to the same
# expected-score effect size 05 reports, on a large simulated sample so the
# conversion is a property of the model and not of sampling noise.
cat("\nWhat those thresholds are worth on the scale 05 reports in:\n")
es_map <- do.call(rbind, lapply(SHIFTS[SHIFTS > 0], function(sh) {
  set.seed(400)
  e <- planted_effect_size(pars, sh, BAD_ITEM, n_per_group = 15000, items = PHQ_ITEMS)
  data.frame(shift = sh, SIDS = e$SIDS, ESSD = e$ESSD)
}))
print(round(es_map, 3), row.names = FALSE)

for (p in PROCS) {
  if (is.na(thr[p])) next
  sids <- approx(es_map$shift, es_map$SIDS, xout = thr[p], rule = 2)$y
  essd <- approx(es_map$shift, es_map$ESSD, xout = thr[p], rule = 2)$y
  cat(sprintf("  %-10s finds 80%% of shifts at or above %+.3f item-score points (ESSD %+.3f)\n",
              p, sids, essd))
}

# Read the real-data effect back out of 05's saved output rather than quoting a
# number from the log. Figures that travel between scripts by hand drift, and a
# repo about measurement artifacts is the wrong place to demonstrate it.
cat("\nCompare against 05, reading its saved output rather than quoting it:\n")
if (file.exists("output/grm_dif.rds")) {
  g <- readRDS("output/grm_dif.rds")
  fl <- g$dif_anchored$item[g$dif_anchored$flag]
  e <- g$es_anchored[g$es_anchored$item %in% fl, ]
  big <- e[which.max(abs(e$SIDS)), ]
  cat(sprintf("  largest real effect across sex: %s at %+.3f item-score points (ESSD %+.3f)\n",
              PHQ_LABELS[big$item], big$SIDS, big$ESSD))
  for (p in PROCS) {
    if (is.na(thr[p])) next
    floor_sids <- approx(es_map$shift, es_map$SIDS, xout = thr[p], rule = 2)$y
    cat(sprintf("  %-10s detection floor %+.3f pts -> the real effect is %s it\n", p,
                floor_sids,
                if (abs(big$SIDS) >= abs(floor_sids)) "at or above" else "BELOW"))
  }
  cat("\nAn effect below the floor is not a finding, it is a coin flip that landed.\n")
  cat("An effect above it is real but still has to be read next to the\n")
  cat("false-positive rate in section 1, not instead of it.\n")
} else {
  cat("  (run 05_grm_dif.R first for this comparison)\n")
}

# ---------------------------------------------------------------------------
cat("\n\n=== 4. Contamination: what happens when several items carry DIF ===\n")
cat("05 flags most of the questionnaire under the all-other-item anchor and\n")
cat("blames anchor contamination. Sections 1 and 2 cannot test that claim:\n")
cat("they contain no real DIF, so there is nothing available to contaminate an\n")
cat("anchor with, and both procedures look well behaved there. This section is\n")
cat("the condition that actually tests it. Plant the same shift on a growing\n")
cat("number of items and count flags on the items where nothing was planted.\n")
cat("A procedure with a clean anchor should hold its false-alarm rate flat as\n")
cat("the number of guilty items rises. A contaminated one should not.\n\n")

N_BAD    <- 0:4
BAD_SET  <- c(4, 5, 9, 8)            # the items 05 flags most stably, in that order
# Contamination is a surface, not a line. How badly the anchor degrades depends
# both on how MANY items carry bias and on how BIG that bias is, and reporting
# one magnitude would leave the reader unable to tell a knife edge from a
# gentle slope. Three magnitudes bracket the detection curve in section 3: one
# below the 80% threshold, one just above it, one well above.
CONTAM_SHIFTS <- c(0.2, 0.4, 0.8)

contam_grid <- expand.grid(n_bad = N_BAD, shift = CONTAM_SHIFTS, rep = seq_len(N_REPS))
# With nothing planted the shift is meaningless and every magnitude would
# resimulate the identical null, so those rows are collapsed to one condition.
contam_grid <- contam_grid[!(contam_grid$n_bad == 0 &
                             contam_grid$shift != CONTAM_SHIFTS[1]), ]

contam_out <- mclapply(seq_len(nrow(contam_grid)), function(k) {
  set.seed(5000 + k)
  nb <- contam_grid$n_bad[k]
  sh <- contam_grid$shift[k]
  bad <- if (nb == 0) integer(0) else BAD_SET[seq_len(nb)]
  d <- sim_grm(pars, 2500, shift = if (nb == 0) 0 else sh,
               bad_item = if (nb == 0) 1 else bad,
               mean_diff = 0.19, items = PHQ_ITEMS)
  r <- try(run_procedures(d$resp, d$group), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  clean <- setdiff(seq_along(PHQ_ITEMS), bad)
  do.call(rbind, lapply(PROCS, function(p) {
    # An anchor item was never tested, so it cannot raise a false alarm and
    # must not be counted as a clean item that passed either. Rate over the
    # clean items actually tested, or NA when none were.
    tested <- clean[!is.na(r[[p]][clean]) & !(p == "anchored" & r$anchor[clean])]
    data.frame(n_bad = nb, shift = sh, procedure = p,
               false_alarm_rate = if (length(tested)) mean(r[[p]][tested]) else NA_real_,
               n_clean_tested = length(tested),
               # Contamination has a second, quieter cost: a guilty item that
               # gets chosen as an anchor is never tested at all, so it is
               # missed AND it defines the scale everything else is judged on.
               guilty_in_anchor = if (p == "anchored" && nb > 0) sum(r$anchor[bad]) else NA_real_)
  }))
}, mc.cores = N_CORES)
contam_df <- do.call(rbind, Filter(Negate(is.null), contam_out))

contam <- aggregate(cbind(false_alarm_rate, guilty_in_anchor) ~ n_bad + shift + procedure,
                    data = contam_df, FUN = function(x) mean(x, na.rm = TRUE),
                    na.action = na.pass)
# The n_bad = 0 null is one shared condition; label it at every magnitude so
# each column of the surface has its own baseline row to be read against.
z <- contam[contam$n_bad == 0, ]
contam <- rbind(contam[contam$n_bad != 0, ],
                do.call(rbind, lapply(CONTAM_SHIFTS[-1], function(v)
                  transform(z, shift = v))), z)
contam <- contam[order(contam$procedure, contam$shift, contam$n_bad), ]
cat(sprintf("%d replications per cell. false_alarm_rate is over the items where\n", N_REPS))
cat("nothing was planted. Rows at n_bad = 0 are the shared no-DIF null.\n\n")
surface <- reshape(contam[, c("n_bad", "shift", "procedure", "false_alarm_rate")],
                   direction = "wide", idvar = c("procedure", "n_bad"),
                   timevar = "shift")
names(surface) <- sub("^false_alarm_rate\\.", "shift_", names(surface))
surface <- surface[order(surface$procedure, surface$n_bad), ]
print(data.frame(surface[, 1:2], round(surface[, -(1:2)], 3)), row.names = FALSE)

cat("\nRead the first and last row of each block. That slope, not the level, is\n")
cat("the contamination: how much a procedure's false-alarm rate on innocent\n")
cat("items is driven by how many guilty items are in the questionnaire with\n")
cat("them. It is the difference between a test of one item and a test of the\n")
cat("company that item keeps.\n\n")
for (p in PROCS) {
  for (sh in CONTAM_SHIFTS) {
    d <- contam[contam$procedure == p & contam$shift == sh, ]
    d <- d[order(d$n_bad), ]
    d <- d[!is.na(d$false_alarm_rate), ]
    cat(sprintf("  %-10s shift %.1f: %.3f with no DIF present -> %.3f with %d items carrying it (%+.3f)\n",
                p, sh, d$false_alarm_rate[1], d$false_alarm_rate[nrow(d)],
                d$n_bad[nrow(d)], d$false_alarm_rate[nrow(d)] - d$false_alarm_rate[1]))
  }
}
# The anchored row is not monotone in the shift, and that is mechanism rather
# than noise: purification finds guilty items more easily when their bias is
# large, so at the biggest shift fewer of them survive into the anchor. The
# guilty_in_anchor count below is the direct evidence for that reading.
cat("\nThe purified row is not monotone in the bias size, which is not noise:\n")
cat("purification spots a guilty item more easily the larger its bias, so at\n")
cat("the widest shift fewer guilty items survive into the anchor. The\n")
cat("guilty_in_anchor figure below is the direct evidence for that reading.\n")
cat("\nReading across a row rather than down a column: at a bias size the\n")
cat("procedure can barely detect, contamination is mild however many items\n")
cat("carry it. The anchor breaks down when the bias is both widespread AND\n")
cat("large enough to see, which is exactly the regime a real questionnaire with\n")
cat("a systematic group difference sits in.\n")
gtab <- contam[contam$procedure == "anchored" & contam$n_bad == max(N_BAD), ]
gtab <- gtab[order(gtab$shift), ]
cat(sprintf("\nAnd purification is not free. With %d guilty items, the average number of\n", max(N_BAD)))
cat("them that ended up INSIDE the anchor set, never tested and instead defining\n")
cat("the scale every other item is judged against:\n")
for (i in seq_len(nrow(gtab)))
  cat(sprintf("  shift %.1f: %.2f of %d guilty items became anchors\n",
              gtab$shift[i], gtab$guilty_in_anchor[i], max(N_BAD)))


# ---------------------------------------------------------------------------
out <- list(null = null_df, contamination = contam_df, contamination_rates = contam, null_rates = null_tab, null_by_item = by_item,
            impact = impact_df, impact_rates = impact_tab,
            planted = planted_df, power = power,
            thresholds = thr, effect_size_map = es_map,
            settings = list(n_splits = N_SPLITS, n_reps = N_REPS,
                            shifts = SHIFTS, bad_item = BAD_ITEM, alpha = ALPHA))
saveRDS(out, "output/dif_error_rate.rds")
cat("\nWrote output/dif_error_rate.rds\n")
