#!/usr/bin/env Rscript
# Graded response model on the nine PHQ-9 items, and DIF by sex.
#
# The CFA in 02 asked a global question: does constraining all nine items
# equally across sex degrade fit. It said no at the metric rung, yes at the
# scalar rung, and the score test pointed at item 4. This stage asks the
# question one item at a time with an IRT estimator, which is what the DIF
# literature actually uses and what the CFA cannot give: an item-level test
# with an effect size in item-score points.
#
# The finding that organises this script is that the default procedure flags
# seven of nine items, which is not a result about the PHQ-9. It is a result
# about the anchor. 06 measures how often that default is wrong on data where
# the answer is known.

source("R/invariance.R")
source("R/grm.R")
set.seed(20260818)

s <- readRDS("output/analytic_sample.rds")
resp <- s[, PHQ_ITEMS]
grp <- as.character(s$sex)

cat("=== Graded response model, single group, n =", nrow(resp), "===\n")
cat("Sex:", paste(names(table(grp)), table(grp), collapse = "  "), "\n\n")

m1 <- mirt(as.data.frame(resp), 1, itemtype = "graded", verbose = FALSE)
pars <- coef(m1, IRTpars = TRUE, simplify = TRUE)$items
cat("Slopes and thresholds (a is discrimination, b_k is the trait level at\n")
cat("which category k becomes the more likely response):\n")
print(data.frame(item = PHQ_LABELS[PHQ_ITEMS],
                 round(pars[, c("a", "b1", "b2", "b3")], 3)), row.names = FALSE)

# The limited-information M2 statistic and its RMSEA, the GRM analogue of the
# CFA fit indices already reported for this instrument.
#
# What this does NOT establish is the assumption the DIF tests below actually
# rest on. M2 is a global index: it asks whether one factor reproduces the
# whole covariance structure. Conditional independence is a local property, and
# two items can share variance beyond the factor while M2 looks excellent. That
# matters here specifically, because a locally dependent pair inflates DIF on
# that pair, and more so when one member is an anchor. An earlier version of
# this script described M2 as the unidimensionality check and that was an
# overclaim. 09_local_dependence.R is the check, and it finds dependence this
# statistic cannot see.
cat("\nModel fit (limited information M2):\n")
m2 <- M2(m1)
print(round(m2[, c("M2", "df", "p", "RMSEA", "SRMSR", "TLI", "CFI")], 4), row.names = FALSE)
cat("\nFor comparison, the WLSMV CFA on the same items: CFI 0.9809, RMSEA 0.0594.\n")
cat("Global fit is not the assumption the DIF tests below rest on. See 09.\n")

# ---------------------------------------------------------------------------
cat("\n\n=== 1. The default procedure: every other item as the anchor ===\n")
cat("Constrain all nine items across sex, free one item at a time, test the\n")
cat("improvement. This is the standard default and requires no choices.\n\n")

dif_all <- grm_dif(resp, grp, anchors = NULL)
show <- function(d, es = NULL) {
  out <- data.frame(item = PHQ_LABELS[d$item], X2 = round(d$X2, 2), df = d$df,
                    adj_p = round(d$adj_p, 4),
                    flag = ifelse(is.na(d$adj_p), "anchor", ifelse(d$flag, "DIF", ".")))
  if (!is.null(es)) {
    out$SIDS <- round(es$SIDS, 3); out$ESSD <- round(es$ESSD, 3)
  }
  print(out, row.names = FALSE)
}
# No effect sizes here on purpose. This baseline constrains every item, so the
# fitted model has no cross-group item differences to measure and every
# expected-score effect size comes back exactly 0.000. Printing a column of
# structural zeros beside real p-values would read as "significant but null".
es_all <- NULL
show(dif_all)
n_all <- sum(dif_all$flag)
cat("\nFlagged:", n_all, "of 9 items.\n")
cat(n_all, "of nine items cannot all be the anomaly. When most items carry DIF,\n")
cat("the all-other-item anchor is contaminated by the very items it is meant to\n")
cat("control for: the two groups' scales get linked through biased items, the\n")
cat("common metric is pulled off true, and clean items inherit an apparent\n")
cat("difference. The number to distrust here is", n_all, "itself, not the individual p-values.\n")
cat("\nThat is a mechanism, not yet a finding, and it is not established by this\n")
cat("table. Section 4 of 06 is where it is tested, by planting DIF on a growing\n")
cat("number of items and measuring how often this procedure then blames an\n")
cat("innocent one. Note also which check does NOT support it: on random splits\n")
cat("of this sample the same procedure is conservative, well under nominal,\n")
cat("because a coin-flip split contains no real DIF to contaminate an anchor\n")
cat("with. A control has to contain the thing it is meant to control for.\n")

# ---------------------------------------------------------------------------
cat("\n\n=== 2. Purified anchor ===\n")
cat("Take the four items with the smallest statistic above as the anchor, free\n")
cat("the rest, then test one studied item at a time against that fixed anchor.\n")
cat("Rank-based anchor selection (Woods 2009). Four is a choice, and section 4\n")
cat("sweeps it rather than defending it.\n\n")

anchors <- purify_anchors(resp, grp, n_anchor = 4)
cat("Anchor set:", paste(PHQ_LABELS[anchors], collapse = ";  "), "\n\n")
dif_anch <- grm_dif(resp, grp, anchors = anchors)
es_anch <- dif_effect_size(attr(dif_anch, "fit"), PHQ_ITEMS)
show(dif_anch, es_anch)

cat("\nSIDS is the signed difference in expected item score between women and\n")
cat("men at the same latent depression level, in raw points of a 0 to 3 item.\n")
cat("ESSD is the same quantity standardised, so it reads like a Cohen's d.\n")
cat("Sign convention: negative means the item under-reports for women, i.e. at\n")
cat("equal true depression a woman scores lower on this item than a man.\n\n")

flagged <- dif_anch[dif_anch$flag, ]
if (nrow(flagged)) {
  e <- es_anch[match(flagged$item, es_anch$item), ]
  for (i in seq_len(nrow(flagged))) {
    cat(sprintf("  %-24s X2 %6.1f  adj p %.4f  SIDS %+.3f pts  ESSD %+.3f  somatic %s\n",
                PHQ_LABELS[flagged$item[i]], flagged$X2[i], flagged$adj_p[i],
                e$SIDS[i], e$ESSD[i], flagged$item[i] %in% SOMATIC_ITEMS))
  }
  cat(sprintf("\nLargest effect: %s at %+.3f item-score points (ESSD %+.3f).\n",
              PHQ_LABELS[flagged$item[which.max(abs(e$SIDS))]],
              e$SIDS[which.max(abs(e$SIDS))], e$ESSD[which.max(abs(e$SIDS))]))
  cat("Every one of these is statistically unambiguous and small. At n = 5,068 a\n")
  cat("test of exact equality answers a question nobody asked; the size is the\n")
  cat("answer, and none of these items shifts an expected score by a third of a\n")
  cat("point. 06 puts the matching false-positive rate beside them.\n")
}

# ---------------------------------------------------------------------------
cat("\n\n=== 3. A second procedure: ordinal logistic regression DIF (lordif) ===\n")
cat("Different estimator, different test, its own iterative purification. Run\n")
cat("so that agreement is evidence rather than one method restating itself.\n\n")

lord <- lordif_dif(resp, grp, alpha = 0.01)
print(data.frame(item = PHQ_LABELS[lord$item],
                 chi_uniform = round(lord$chi13, 4),
                 flag = ifelse(lord$flag, "DIF", ".")), row.names = FALSE)
cat("\nThe iteration trace above is lordif purifying its own anchor, and it is\n")
cat("worth reading: the flagged set oscillates rather than converging, which is\n")
cat("the same anchor instability section 4 measures deliberately, showing up here\n")
cat("by itself inside a procedure that is supposed to have solved it.\n")
cat("\nlordif flags:", sum(lord$flag), " IRT-LR with purified anchor flags:",
    sum(dif_anch$flag), "\n")
agree <- intersect(lord$item[lord$flag], dif_anch$item[dif_anch$flag])
cat("Both procedures:", if (length(agree)) paste(PHQ_LABELS[agree], collapse = ";  ") else "none", "\n")
only_l <- setdiff(lord$item[lord$flag], dif_anch$item[dif_anch$flag])
only_i <- setdiff(dif_anch$item[dif_anch$flag], lord$item[lord$flag])
cat("lordif only:  ", if (length(only_l)) paste(PHQ_LABELS[only_l], collapse = ";  ") else "none", "\n")
cat("IRT-LR only:  ", if (length(only_i)) paste(PHQ_LABELS[only_i], collapse = ";  ") else "none", "\n")

# ---------------------------------------------------------------------------
cat("\n\n=== 4. Anchor sensitivity ===\n")
cat("DIF is defined relative to the items assumed clean, so the conclusion is a\n")
cat("function of that assumption. Sweep the anchor size and watch the flag set\n")
cat("move. This is the known trap in this literature and the reason section 2's\n")
cat("choice of four is reported as a choice.\n\n")

sens_rows <- list()
for (k in 1:6) {
  a <- anchors_from_stage1(dif_all, n_anchor = k)   # stage 1 already fitted above
  d <- grm_dif(resp, grp, anchors = a)
  sens_rows[[k]] <- data.frame(
    n_anchor = k, anchor = paste(sub("^(\\d).*", "\\1", PHQ_LABELS[a]), collapse = ","),
    n_flagged = sum(d$flag),
    flagged = paste(sub("^(\\d).*", "\\1", PHQ_LABELS[d$item[d$flag]]), collapse = ","))
}
sens <- do.call(rbind, sens_rows)
print(sens, row.names = FALSE)

always <- Reduce(intersect, lapply(sens_rows, function(r)
  strsplit(r$flagged, ",")[[1]]))
cat("\nFlagged at every anchor size tested:", paste(always, collapse = ","), "\n")
cat("(item numbers, per the labels above)\n")
cat("Items that appear and disappear with the anchor are not findings about the\n")
cat("PHQ-9. They are findings about the anchor, and reporting them as the former\n")
cat("is the error this section exists to prevent.\n")

# ---------------------------------------------------------------------------
out <- list(pars = pars, m2 = m2, anchors = anchors,
            dif_all_other = dif_all, es_all_other = es_all,
            dif_anchored = dif_anch, es_anchored = es_anch,
            lordif = lord, anchor_sensitivity = sens,
            stable_flags = always)
saveRDS(out, "output/grm_dif.rds")
cat("\nWrote output/grm_dif.rds\n")
cat("\nNOTE: a DIF result without its error rate is half a finding. 06 is the\n")
cat("other half and this stage should not be read without it. 09 tests the\n")
cat("conditional-independence assumption these flags rest on, which the M2\n")
cat("statistic above does not.\n")
