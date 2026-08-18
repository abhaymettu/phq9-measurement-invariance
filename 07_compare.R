#!/usr/bin/env Rscript
# Does the IRT DIF test flag the same item the CFA called non-invariant?
#
# This is the optional stage and the interesting one. Two estimators, two
# parameterisations, one instrument, one grouping variable. If they agree, the
# non-invariance is a property of the PHQ-9 rather than of lavaan. If they
# disagree, the disagreement is the finding, because everything published about
# sex bias in this questionnaire rests on one of these two methods and almost
# never on both.
#
# Nothing is refit here. Everything is read from the saved output of 02, 05 and
# 06, which is also the rule that keeps the numbers in the writeup identical to
# the numbers in the logs.

source("R/invariance.R")
source("R/grm.R")

need <- c("output/invariance_results.rds", "output/grm_dif.rds")
missing <- need[!file.exists(need)]
if (length(missing)) stop("run the earlier stages first, missing: ",
                          paste(missing, collapse = ", "))

inv <- readRDS("output/invariance_results.rds")
grm <- readRDS("output/grm_dif.rds")

cat("=== What each method was asked ===\n")
cat("CFA (02): does constraining this item's thresholds and loading equal across\n")
cat("  sex degrade global fit, measured by a score test on the constraint.\n")
cat("IRT (05): holding the latent trait fixed, do women and men differ in the\n")
cat("  probability of endorsing this item, measured by a likelihood ratio test\n")
cat("  against a fixed anchor set.\n")
cat("These are related but not the same question, and the ranking below is the\n")
cat("only place in this repo where that difference is visible.\n\n")

# The CFA side. 02 stored the WLSMV score test by item; recompute the ranking
# from the stored object rather than retyping the log.
cfa_rank <- inv$scores[["sex WLSMV"]]
if (is.null(cfa_rank)) {
  stop("expected invariance_results.rds$scores to carry \"sex WLSMV\"; found: ",
       paste(names(inv$scores), collapse = ", "))
}
cfa_rank <- cfa_rank[order(-cfa_rank$X2), ]
cfa_rank$cfa_rank <- seq_len(nrow(cfa_rank))

irt <- grm$dif_anchored
irt$irt_rank <- rank(-replace(irt$X2, is.na(irt$X2), -Inf), ties.method = "min")
irt$irt_rank[is.na(irt$X2)] <- NA          # anchors were never tested
es <- grm$es_anchored

cmp <- data.frame(
  item = PHQ_LABELS[PHQ_ITEMS],
  somatic = PHQ_ITEMS %in% SOMATIC_ITEMS,
  cfa_X2 = round(cfa_rank$X2[match(PHQ_ITEMS, cfa_rank$item)], 1),
  cfa_rank = cfa_rank$cfa_rank[match(PHQ_ITEMS, cfa_rank$item)],
  irt_X2 = round(irt$X2[match(PHQ_ITEMS, irt$item)], 1),
  irt_rank = irt$irt_rank[match(PHQ_ITEMS, irt$item)],
  irt_flag = ifelse(irt$anchor[match(PHQ_ITEMS, irt$item)], "anchor",
                    ifelse(irt$flag[match(PHQ_ITEMS, irt$item)], "DIF", ".")),
  SIDS = round(es$SIDS[match(PHQ_ITEMS, es$item)], 3))

cat("=== Item by item ===\n")
print(cmp, row.names = FALSE)

# ---------------------------------------------------------------------------
cat("\n=== Agreement ===\n")
# The CFA released exactly one item to reach partial scalar invariance; that is
# its verdict, not the whole ranking.
cfa_released <- inv$freed$WLSMV
if (!is.null(cfa_released)) {
  cat("CFA released to reach partial invariance:",
      paste(PHQ_LABELS[cfa_released], collapse = ";  "), "\n")
}
irt_flag <- irt$item[which(irt$flag)]
cat("IRT flagged (purified anchor):     ", paste(PHQ_LABELS[irt_flag], collapse = ";  "), "\n")
cat("IRT flagged at every anchor size:  ", paste(grm$stable_flags, collapse = ","),
    " (item numbers)\n")

both <- intersect(cfa_released %||% character(0), irt_flag)
cat("\nBoth methods:", if (length(both)) paste(PHQ_LABELS[both], collapse = ";  ") else "none", "\n")

# Rank correlation over the items the IRT side actually tested. Spearman on
# five items is a weak statistic and is labelled as such rather than dressed up.
usable <- !is.na(cmp$irt_rank)
rho <- suppressWarnings(cor(cmp$cfa_X2[usable], cmp$irt_X2[usable], method = "spearman"))
cat(sprintf("\nSpearman correlation of the two rankings, over the %d items the IRT\n",
            sum(usable)))
cat(sprintf("side tested: %+.2f. With %d points this is descriptive, not a test;\n",
            rho, sum(usable)))
cat("it is here to show the direction of agreement, not to certify it.\n")

# ---------------------------------------------------------------------------
cat("\n=== The somatic reading ===\n")
cat("The standing claim in this literature is that somatic items (sleep,\n")
cat("energy, appetite, psychomotor) are what make the PHQ-9 behave differently\n")
cat("across sex. Both methods get asked whether they support it.\n\n")
som <- cmp$somatic
cat(sprintf("CFA: mean score statistic, somatic %.1f vs non-somatic %.1f\n",
            mean(cmp$cfa_X2[som]), mean(cmp$cfa_X2[!som])))
cat(sprintf("IRT: %d of %d flagged items are somatic (%d of 9 items are)\n",
            sum(PHQ_ITEMS %in% irt_flag & som), length(irt_flag), sum(som)))
neg <- cmp[!is.na(cmp$SIDS) & cmp$SIDS < 0 & cmp$irt_flag == "DIF", ]
cat(sprintf("IRT: %d flagged items have a negative SIDS, meaning they under-report\n", nrow(neg)))
if (nrow(neg)) {
  cat("     for women at equal depression:\n")
  for (i in seq_len(nrow(neg)))
    cat(sprintf("       %-24s %+.3f pts  somatic %s\n", neg$item[i], neg$SIDS[i], neg$somatic[i]))
}

# ---------------------------------------------------------------------------
cat("\n=== Reading this honestly ===\n")
if (file.exists("output/dif_error_rate.rds")) {
  err <- readRDS("output/dif_error_rate.rds")
  r <- err$null_rates["anchored", "per_test"]
  ps <- err$null_rates["anchored", "per_split"]
  cat(sprintf("On random splits of this same sample the anchored IRT procedure\n"))
  cat(sprintf("flags an item %.1f%% of the time with nothing there, and %.0f%% of\n",
              100 * r, 100 * ps))
  cat(sprintf("splits produce at least one flag. The agreement above has to clear\n"))
  cat(sprintf("that bar before it means the two methods found the same real thing\n"))
  cat(sprintf("rather than the same artifact.\n"))
} else {
  cat("Run 06_dif_error_rate.R: agreement between two methods means nothing\n")
  cat("until the false-positive rate of each is on the page beside it.\n")
}

saveRDS(list(comparison = cmp, spearman = rho, both = both),
        "output/irt_cfa_comparison.rds")
cat("\nWrote output/irt_cfa_comparison.rds\n")
