#!/usr/bin/env Rscript
# Audit the raw files before modelling anything. Row counts, duplicate keys,
# the 7/9 codes, item-level missingness, category sparsity, and how the
# complete-case rule moves the sample.

source("R/invariance.R")
suppressMessages(library(psych))

dpq  <- read_xpt("data/DPQ_J.xpt")
demo <- read_xpt("data/DEMO_J.xpt")

cat("=== Raw files ===\n")
cat(sprintf("DPQ_J   %d rows, %d cols\n", nrow(dpq), ncol(dpq)))
cat(sprintf("DEMO_J  %d rows, %d cols\n", nrow(demo), ncol(demo)))
cat("DPQ columns:", paste(names(dpq), collapse = " "), "\n")
cat(sprintf("duplicate SEQN: DPQ %d, DEMO %d\n",
            sum(duplicated(dpq$SEQN)), sum(duplicated(demo$SEQN))))
cat(sprintf("DPQ100 present in file: %s (excluded from the factor, it is the\n",
            "DPQ100" %in% names(dpq)))
cat("  functional-impairment item, not a symptom item)\n\n")

cat("=== Special codes on the nine items, raw ===\n")
codes <- t(vapply(PHQ_ITEMS, function(v) {
  x <- dpq[[v]]
  c(refused_7 = sum(x == 7, na.rm = TRUE), dk_9 = sum(x == 9, na.rm = TRUE),
    na = sum(is.na(x)), max_valid = suppressWarnings(max(x[x <= 3], na.rm = TRUE)))
}, numeric(4)))
print(codes)
cat("\nAny value outside 0-3 after recode:",
    sum(!is.na(unlist(dpq[PHQ_ITEMS])) & !unlist(dpq[PHQ_ITEMS]) %in% c(0:3, 7, 9)), "\n\n")

dat <- load_phq9()
cat("=== Sample flow ===\n")
cat(sprintf("DEMO inner-join DPQ (DPQ is asked of 12+ only): %d\n", nrow(dat)))
a <- dat[dat$RIDAGEYR >= 18, ]
cat(sprintf("age 18+:                                       %d\n", nrow(a)))
s <- analytic_sample(dat)
cat(sprintf("all nine items answered (complete case):        %d  (%.1f%% of 18+)\n",
            nrow(s), 100 * nrow(s) / nrow(a)))
cat(sprintf("dropped for partial or missing PHQ-9:           %d\n", nrow(a) - nrow(s)))
cat(sprintf("of those, answered zero items:                  %d\n",
            sum(a$n_answered == 0)))
cat(sprintf("of those, answered 1 to 8 items:                %d\n",
            sum(a$n_answered > 0 & a$n_answered < 9)))
cat(sprintf("zero MEC weight in the analytic sample:         %d\n", sum(s$WTMEC2YR == 0)))

cat("\n=== Is missingness related to the outcome? ===\n")
# Partial responders cannot be scored, so check whether they look different on
# what they did answer. If they are much more symptomatic, complete-case is
# dropping the sick tail and the sample is not missing at random.
part <- a[a$n_answered > 0 & a$n_answered < 9, ]
cat(sprintf("mean per-item score, partial responders: %.3f (n = %d)\n",
            mean(rowSums(part[PHQ_ITEMS], na.rm = TRUE) / part$n_answered), nrow(part)))
cat(sprintf("mean per-item score, complete responders: %.3f (n = %d)\n",
            mean(s$phq9 / 9), nrow(s)))

cat("\n=== Item distributions in the analytic sample (proportions) ===\n")
tab <- t(vapply(PHQ_ITEMS, function(v) prop.table(table(factor(s[[v]], 0:3))), numeric(4)))
rownames(tab) <- PHQ_LABELS[rownames(tab)]
print(round(tab, 4))
cat("\nSparsest cell across all items and categories:", sprintf("%.4f", min(tab)), "\n")
cat("Smallest cell count:", min(round(tab * nrow(s))), "\n")
cat("Categories with fewer than 10 respondents anywhere:",
    sum(round(tab * nrow(s)) < 10), "\n")
cat("This sparsity, concentrated in category 3 of item 9, is the reason WLSMV\n")
cat("with a polychoric basis is preferred over treating the items as continuous.\n")

cat("\n=== Group sizes for the invariance tests ===\n")
print(table(s$sex))
print(table(s$age_grp))
print(table(s$race))

cat("\n=== Reliability and descriptives ===\n")
al <- psych::alpha(s[PHQ_ITEMS], warnings = FALSE)
cat(sprintf("Cronbach alpha: %.3f\n", al$total$raw_alpha))
om <- suppressMessages(psych::omega(s[PHQ_ITEMS], nfactors = 1, plot = FALSE, flags = FALSE))
cat(sprintf("Omega total:    %.3f\n", om$omega.tot))
cat(sprintf("PHQ-9 sum: mean %.2f, sd %.2f, median %d, max %d\n",
            mean(s$phq9), sd(s$phq9), median(s$phq9), max(s$phq9)))
cat(sprintf("share scoring >= 10 (moderate or worse), unweighted: %.1f%%\n",
            100 * mean(s$phq9 >= 10)))

saveRDS(s, "output/analytic_sample.rds")
cat("\nWrote output/analytic_sample.rds\n")
