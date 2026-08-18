#!/usr/bin/env Rscript
# Is the DIF in 05 real, or is it local dependence?
#
# 05 checks the unidimensionality assumption with the limited-information M2
# statistic. That was the wrong check for this purpose and this stage exists to
# replace it. M2 is a global fit index: it asks whether the one-factor model
# reproduces the whole covariance structure. Two items can share variance
# beyond the factor while global fit stays excellent, and a locally dependent
# pair inflates DIF flags on that pair, because the extra shared variance has
# nowhere else to go.
#
# The reason this matters here rather than in general: the PHQ-9 somatic items
# (sleep, energy, appetite) are the pair most often reported as locally
# dependent in this literature, and item 4 (low energy) is the largest effect
# in 05. So the leading alternative explanation for the headline finding was,
# until this stage, untested.
#
# Three questions, in order:
#   1 Is there local dependence, judged against a null simulated from this
#     instrument rather than against a cutoff quoted from a paper?
#   2 Would the statistic have found it if it were there? (Rule 1 of this repo:
#     prove the check bites before believing a negative.)
#   3 Do the DIF flags in 05 survive the standard remedies for it?

source("R/invariance.R")
source("R/grm.R")
suppressMessages(library(parallel))
set.seed(20260818)

N_NULL  <- as.integer(Sys.getenv("N_NULL", "200"))
N_CORES <- as.integer(Sys.getenv("N_CORES", max(1, detectCores() - 2)))

s <- readRDS("output/analytic_sample.rds")
resp <- s[, PHQ_ITEMS]
grp <- as.character(s$sex)
prev <- readRDS("output/grm_dif.rds")
anchors0 <- prev$anchors
flags0 <- prev$dif_anchored$item[prev$dif_anchored$flag]

lab <- function(x) PHQ_LABELS[x]
# Item numbers for compact display. Variants introduce a testlet item that is
# not one of the nine, so fall back to its own name rather than emitting NA.
short <- function(x) ifelse(x %in% names(PHQ_LABELS),
                            sub("^(\\d+).*", "\\1", PHQ_LABELS[x]), x)

# ---------------------------------------------------------------------------
cat("=== 1. Residual dependence in the observed data ===\n")
q3 <- q3_matrix(resp)
cat(sprintf("Raw Q3 carries a mechanical negative bias of about -1/(J-1) = %.3f.\n", -1/8))
cat(sprintf("Observed mean off-diagonal Q3 here: %.4f. Everything below is centred\n", q3$bias))
cat("on that mean, which is the standard correction; reading raw Q3 against a\n")
cat("'greater than 0.2' rule would find nothing on a nine-item scale, ever.\n\n")

pairs <- q3_pairs(q3)
top <- head(pairs, 6)
cat("Strongest residual correlations, centred:\n")
print(data.frame(pair = paste(lab(top$item_i), "+", lab(top$item_j)),
                 raw = round(top$raw, 3), centred = round(top$centred, 3)),
      row.names = FALSE)
cat(sprintf("\nQ3max - Q3mean = %.3f\n", q3$stat))

# The general statistic is not the one that matters for this repo's question.
# Dependence between two anchors, or between two untested items, does not
# inflate anything. Dependence between a FLAGGED item and an ANCHOR item does,
# and that is the cell to look at.
cat("\n--- The pairing that actually inflates DIF ---\n")
cat("A flagged item sharing residual variance with an anchor item is the\n")
cat("mechanism: the surplus variance sits inside the very items defining the\n")
cat("scale that the flagged item is being judged against.\n\n")
cat("Anchor set from 05:", paste(lab(anchors0), collapse = " | "), "\n\n")
fa <- do.call(rbind, lapply(flags0, function(it) {
  z <- q3$centred[it, anchors0]
  data.frame(flagged = lab(it), worst_anchor = lab(anchors0[which.max(z)]),
             centred_q3 = max(z))
}))
fa <- fa[order(-fa$centred_q3), ]
print(transform(fa, centred_q3 = round(centred_q3, 3)), row.names = FALSE)

# ---------------------------------------------------------------------------
cat("\n\n=== 2. What that number means, calibrated on this instrument ===\n")
cat(sprintf("Simulating %d samples from a CORRECT unidimensional GRM with this\n", N_NULL))
cat("questionnaire's own parameters and sample size. No local dependence exists\n")
cat("in any of them, so this is what Q3 looks like when the assumption holds.\n")
cat("A cutoff quoted from a simulation study run on other items at another n\n")
cat("would be a weaker basis for the verdict than this.\n\n")

pars <- grm_pars(resp)
null_stats <- unlist(mclapply(seq_len(N_NULL), function(i) {
  set.seed(7000 + i)
  r <- try(q3_matrix(sim_unidim(pars, nrow(resp), PHQ_ITEMS)), silent = TRUE)
  if (inherits(r, "try-error")) return(NA_real_)
  r$stat
}, mc.cores = N_CORES))
null_stats <- null_stats[!is.na(null_stats)]

p_emp <- mean(null_stats >= q3$stat)
cat(sprintf("null Q3max - Q3mean over %d replications: median %.3f, 95th pct %.3f, max %.3f\n",
            length(null_stats), median(null_stats), quantile(null_stats, 0.95), max(null_stats)))
cat(sprintf("observed: %.3f    empirical p = %.3f\n", q3$stat, p_emp))
cat(sprintf("VERDICT: %s\n", if (p_emp < 0.05)
  "the observed dependence is larger than a correct unidimensional model produces" else
  "the observed dependence is within what a correct unidimensional model produces"))

# ---------------------------------------------------------------------------
cat("\n=== Does that statistic bite? ===\n")
cat("A negative result is worth nothing from a check that cannot go positive.\n")
cat("Planting local dependence on the sleep/energy pair by adding a second\n")
cat("factor loading on those two items only, then re-running the identical\n")
cat("statistic. If this does not move, section 2 says nothing.\n\n")

LD_PAIR <- PHQ_ITEMS[c(3, 4)]
bite <- do.call(rbind, lapply(c(0.5, 1.0), function(load) {
  st <- unlist(mclapply(1:20, function(i) {
    set.seed(8000 + i)
    r <- try(q3_matrix(sim_ld(pars, nrow(resp), LD_PAIR, load, PHQ_ITEMS)), silent = TRUE)
    if (inherits(r, "try-error")) NA_real_ else r$stat
  }, mc.cores = N_CORES))
  st <- st[!is.na(st)]
  data.frame(planted_loading = load, median_stat = median(st),
             pct_above_null_95 = mean(st > quantile(null_stats, 0.95)))
}))
print(transform(bite, median_stat = round(median_stat, 3),
                pct_above_null_95 = round(pct_above_null_95, 2)), row.names = FALSE)
cat(sprintf("\nFor reference the null median is %.3f and the observed value is %.3f.\n",
            median(null_stats), q3$stat))

# ---------------------------------------------------------------------------
cat("\n\n=== 3. Do the DIF flags survive the remedies? ===\n")
cat("Whatever section 2 concluded, the flags get re-run under the standard\n")
cat("treatments for local dependence. A verdict of 'not significant' is not a\n")
cat("licence to skip the robustness check, because the pair with the largest\n")
cat("residual correlation is also the pair carrying the largest DIF effect,\n")
cat("and that coincidence is the thing a reviewer will ask about.\n\n")

worst <- pairs[1, ]
wp <- c(worst$item_i, worst$item_j)
cat(sprintf("Treating the strongest pair: %s + %s (centred Q3 %+.3f)\n",
            lab(wp[1]), lab(wp[2]), worst$centred))
second <- pairs[2, ]
sp <- c(second$item_i, second$item_j)
cat(sprintf("and the runner up:          %s + %s (centred Q3 %+.3f)\n\n",
            lab(sp[1]), lab(sp[2]), second$centred))

run_variant <- function(name, r, note, forbid = character(0)) {
  items <- names(r)
  st1 <- grm_dif(r, grp, anchors = NULL, items = items)
  elig <- st1[!st1$item %in% forbid, ]
  a <- anchors_from_stage1(elig, n_anchor = min(4, nrow(elig)))
  d <- grm_dif(r, grp, anchors = a, items = items)
  es <- dif_effect_size(attr(d, "fit"), items)
  list(name = name, note = note, anchors = a, dif = d, es = es,
       flagged = d$item[d$flag])
}

variants <- list()
variants$baseline <- list(name = "baseline (05)", note = "all nine items, purified anchor",
                          anchors = anchors0, dif = prev$dif_anchored,
                          es = prev$es_anchored, flagged = flags0)

# Remedy A. Drop the anchor-side member of the pair. This is the one that
# answers "does item 4 still flag", because item 4 itself stays in the model.
dropped <- if (wp[1] %in% anchors0) wp[1] else wp[2]
kept <- setdiff(wp, dropped)
variants$dropA <- run_variant(
  sprintf("drop %s", short(dropped)),
  resp[, setdiff(PHQ_ITEMS, dropped), drop = FALSE],
  sprintf("removes the anchor-side member of the top pair, so %s is tested with its dependent partner gone",
          short(kept)))

TESTLET <- sprintf("LD%s+%s", short(wp[1]), short(wp[2]))
# Remedy B. Collapse the pair into a testlet. Principled, but it costs the
# ability to test either member on its own, so it answers a different question:
# is there DIF in that content area at all.
variants$testletB <- run_variant(
  sprintf("testlet %s+%s", short(wp[1]), short(wp[2])),
  make_testlet(resp, wp, name = TESTLET),
  "scores the dependent pair as one item instead of two")

# Remedy C. Keep all nine items but forbid the dependent partners from ever
# being anchors. The least destructive remedy and the most direct test of the
# mechanism, since the mechanism is specifically anchor-side contamination.
forbid <- unique(c(wp, sp))
variants$forbidC <- run_variant(
  "LD partners barred from anchor",
  resp,
  "all nine items kept; the four items in the top two dependent pairs cannot serve as anchors",
  forbid = forbid)

cat("Flag sets under each treatment (item numbers):\n\n")
for (v in variants) {
  cat(sprintf("  %-32s anchors %-16s flags %s\n", v$name,
              paste(short(v$anchors), collapse = ","),
              paste(short(v$flagged), collapse = ",")))
  if (!is.null(v$note)) cat(sprintf("  %-32s   (%s)\n", "", v$note))
}

# The specific survival question, in effect-size terms, for the items that can
# still be tested under each treatment.
cat("\nEffect size on the items of interest, where still testable:\n")
watch <- unique(c(kept, flags0, TESTLET))
surv <- do.call(rbind, lapply(variants, function(v) {
  row <- data.frame(treatment = v$name)
  for (it in watch) {
    # An anchor has a SIDS of exactly 0 by construction because it was never
    # tested. Printing that as "+0.000" is indistinguishable from a real
    # near-zero estimate and would read as "tested and cleared", which is a
    # materially different claim. Anchors are labelled as untested instead.
    val <- if (it %in% v$es$item) v$es$SIDS[v$es$item == it] else NA_real_
    fl <- it %in% v$flagged
    row[[short(it)]] <- if (!it %in% v$es$item) "-"
                        else if (it %in% v$anchors) "anchor"
                        else sprintf("%+.3f%s", val, if (fl) "*" else "")
  }
  row
}))
print(surv, row.names = FALSE)
cat("* = flagged at adjusted p < 0.05.  '-' = not in the model under that treatment.\n")
cat("'anchor' = the item served as an anchor and was therefore NOT tested. That is\n")
cat("not the same as being tested and cleared, and the difference matters: an item\n")
cat("that becomes an anchor under a remedy has had its baseline flag neither\n")
cat("confirmed nor contradicted.\n")

base_set <- sort(short(flags0))
for (v in variants[-1]) {
  now <- sort(short(v$flagged))
  cat(sprintf("\n%s: %s\n", v$name,
              if (identical(now, base_set)) "identical flag set to baseline" else
                sprintf("flag set moves from {%s} to {%s}",
                        paste(base_set, collapse = ","), paste(now, collapse = ","))))
  lost <- setdiff(flags0, v$flagged)
  if (length(lost)) {
    became <- lost[lost %in% v$anchors]
    tested <- setdiff(lost, c(became, setdiff(PHQ_ITEMS, v$es$item)))
    if (length(tested))
      cat(sprintf("    tested and no longer flagged: %s\n",
                  paste(short(tested), collapse = ",")))
    if (length(became))
      cat(sprintf("    became an anchor, so NOT tested: %s\n",
                  paste(short(became), collapse = ",")))
  }
}

# Summed across the remedies, the distinction that matters for the writeup.
cat("\nAcross all three remedies:\n")
for (it in flags0) {
  tested_in <- vapply(variants[-1], function(v)
    it %in% v$es$item && !it %in% v$anchors, logical(1))
  flagged_in <- vapply(variants[-1], function(v) it %in% v$flagged, logical(1))
  cat(sprintf("  %-24s tested in %d of %d remedies, flagged in %d\n",
              lab(it), sum(tested_in), length(variants) - 1, sum(flagged_in)))
}

# ---------------------------------------------------------------------------
saveRDS(list(q3 = q3[c("raw", "centred", "bias", "stat")], pairs = pairs,
             flagged_vs_anchor = fa, null_stats = null_stats, p_emp = p_emp,
             bite = bite, variants = lapply(variants, function(v)
               v[c("name", "note", "anchors", "dif", "es", "flagged")]),
             survival = surv, top_pair = wp, second_pair = sp),
        "output/local_dependence.rds")
write.csv(transform(pairs, item_i = lab(item_i), item_j = lab(item_j)),
          "output/table7_local_dependence.csv", row.names = FALSE)
cat("\nWrote output/local_dependence.rds and output/table7_local_dependence.csv\n")
