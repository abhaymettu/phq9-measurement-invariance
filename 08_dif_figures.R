#!/usr/bin/env Rscript
# Two figures for the IRT stage. Both exist to be read together with 05: a DIF
# result and its error rate on the same page, which is the pairing the whole
# extension is for.

source("R/invariance.R")
suppressMessages({library(ggplot2)})

err <- readRDS("output/dif_error_rate.rds")
grm <- readRDS("output/grm_dif.rds")

PROC_LAB <- c(all_other = "All-other-item anchor\n(the common default)",
              anchored  = "Purified 4-item anchor",
              lordif    = "Ordinal logistic (lordif)")
PAL <- c(all_other = "#b4462f", anchored = "#2c3e50", lordif = "#9aa6b2")
base_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "top")

# ---------------------------------------------------------------------------
# Figure 4. Detection curve, with the real data's largest observed effect drawn
# on the same axis. The point of putting them together is that the reader can
# see the real finding sitting near the floor of what the procedure can see,
# which no table of p-values communicates.
pw <- err$power
pw$procedure <- factor(pw$procedure, levels = names(PROC_LAB))

# Convert the x axis from logits to the effect size 05 reports in, so both
# panels of the story are in item-score points.
map <- rbind(data.frame(shift = 0, SIDS = 0, ESSD = 0), err$effect_size_map)
pw$sids <- abs(approx(map$shift, map$SIDS, xout = pw$shift, rule = 2)$y)

fl <- grm$dif_anchored$item[grm$dif_anchored$flag]
e_real <- grm$es_anchored[grm$es_anchored$item %in% fl, ]
biggest <- max(abs(e_real$SIDS))
big_item <- e_real$item[which.max(abs(e_real$SIDS))]

p4 <- ggplot(pw, aes(sids, detection, colour = procedure)) +
  annotate("rect", xmin = 0, xmax = biggest, ymin = -Inf, ymax = Inf,
           fill = "#b4462f", alpha = 0.07) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey30") +
  geom_vline(xintercept = biggest, linetype = "dotted", colour = "#b4462f") +
  annotate("text", x = biggest, y = 0.06, hjust = -0.05, size = 2.9, colour = "#b4462f",
           label = sprintf("largest real sex effect\nin this data: %.3f pts (%s)",
                           biggest, sub("^(\\d).*", "item \\1", PHQ_LABELS[big_item]))) +
  annotate("text", x = max(pw$sids), y = 0.8, hjust = 1, vjust = -0.6,
           size = 2.9, colour = "grey25", label = "80% detection") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_colour_manual(values = PAL, labels = PROC_LAB, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "How large does item bias have to be before the test finds it?",
       subtitle = sprintf("%d replications per point, 2,500 per group, on top of a real 0.19 SD group difference",
                          err$settings$n_reps),
       x = "Planted bias, absolute item-score points (SIDS)",
       y = "Share of replications detecting it",
       caption = paste("The shaded strip is the range where the real PHQ-9 sex effects live.",
                       "Every one of them sits at or near the floor of what these\nprocedures can",
                       "reliably see, which is the reason the flags in 05 are reported with this",
                       "figure and not on their own.")) +
  base_theme
ggsave("output/fig4_detection_curve.png", p4, width = 8, height = 5.2, dpi = 200)

# ---------------------------------------------------------------------------
# Figure 5. The contamination result: false alarms on innocent items as a
# function of how many items around them are guilty. This is the figure that
# justifies not using the default procedure.
cm <- err$contamination_rates
cm$procedure <- factor(cm$procedure, levels = names(PROC_LAB))

cm$panel <- factor(sprintf("planted bias %.1f logits", cm$shift),
                   levels = sprintf("planted bias %.1f logits", sort(unique(cm$shift))))

p5 <- ggplot(cm, aes(n_bad, false_alarm_rate, colour = procedure)) +
  facet_wrap(~ panel) +
  geom_hline(yintercept = err$settings$alpha, linetype = "dashed", colour = "grey30") +
  annotate("text", x = max(cm$n_bad), y = err$settings$alpha, hjust = 1, vjust = 1.7,
           size = 2.9, colour = "grey25", label = "nominal 0.05") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_colour_manual(values = PAL, labels = PROC_LAB, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "A DIF test is not a test of one item. It is a test of the company it keeps.",
       subtitle = paste("False positives on items where nothing was planted, as more of the other items carry real bias.",
                        "Swept across three magnitudes: how badly an anchor degrades depends on the size of the bias,",
                        "not only on how many items carry it.", sep = "\n"),
       x = "Number of items carrying planted bias (of 9)",
       y = "False positive rate on the innocent items",
       caption = paste("With nothing planted every procedure is conservative, which is why the",
                       "random-split control alone looks reassuring\nand cannot settle this.",
                       "The default anchor degrades as soon as there is real bias for it to be",
                       "contaminated by.")) +
  base_theme
ggsave("output/fig5_anchor_contamination.png", p5, width = 8, height = 5.2, dpi = 200)

# ---------------------------------------------------------------------------
# Figure 6. The local-dependence verdict. The conventional cutoff and the
# instrument's own null disagree here, and that disagreement is the argument
# for simulating rather than citing, so both belong on the same axis.
ld <- readRDS("output/local_dependence.rds")
nulls <- data.frame(stat = ld$null_stats)

p6 <- ggplot(nulls, aes(stat)) +
  geom_histogram(bins = 30, fill = "#9aa6b2", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = quantile(ld$null_stats, 0.95), linetype = "dashed",
             colour = "grey30") +
  geom_vline(xintercept = ld$q3$stat, colour = "#b4462f", linewidth = 0.8) +
  geom_vline(xintercept = 0.20, linetype = "dotted", colour = "#2c3e50", linewidth = 0.6) +
  annotate("text", x = quantile(ld$null_stats, 0.95), y = Inf, vjust = 1.8, hjust = -0.05,
           size = 2.9, colour = "grey25",
           label = sprintf("null 95th pct %.3f", quantile(ld$null_stats, 0.95))) +
  annotate("text", x = ld$q3$stat, y = Inf, vjust = 3.4, hjust = -0.05,
           size = 3.1, colour = "#b4462f",
           label = sprintf("observed %.3f\n(sleep + energy)", ld$q3$stat)) +
  annotate("text", x = 0.20, y = Inf, vjust = 1.8, hjust = 1.05,
           size = 2.9, colour = "#2c3e50",
           label = "conventional cutoff 0.20:\nwould report no dependence") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "The usual local-dependence cutoff would have cleared this instrument",
       subtitle = sprintf("Q3max - Q3mean across %d samples simulated from a CORRECT unidimensional model with this questionnaire's
own parameters and sample size. The observed value exceeds all %d.",
                          length(ld$null_stats), length(ld$null_stats)),
       x = "Q3max - Q3mean", y = "simulated datasets with no local dependence",
       caption = paste("Raw Q3 is centred to remove its mechanical -1/(J-1) bias.",
                       "A cutoff derived on other items at another sample size\nis a weaker basis for the verdict",
                       "than the instrument's own null, and here it is also the wrong one.")) +
  base_theme
ggsave("output/fig6_local_dependence_null.png", p6, width = 8, height = 5, dpi = 200)

# ---------------------------------------------------------------------------
# Figure 7. Whether the flags survive the remedies. The direction is the point:
# treating the dependence makes the headline item larger, not smaller.
# An item serving as an anchor has a SIDS of exactly 0 by construction: it was
# never tested. Plotting that as an estimate would imply the effect moved when
# in fact it was not measured, so anchor points are drawn hollow and labelled
# as untested rather than silently read as "no effect".
SHORT <- c("baseline (05)" = "baseline", "testlet 4+3" = "testlet 4+3",
           "LD partners barred from anchor" = "no-LD anchor")
surv <- do.call(rbind, lapply(ld$variants, function(v) {
  keep <- v$es$item[v$es$item %in% names(PHQ_LABELS)]
  data.frame(treatment = v$name, item = keep,
             SIDS = v$es$SIDS[match(keep, v$es$item)],
             anchor = keep %in% v$anchors,
             flagged = keep %in% v$flagged)
}))
watch <- unique(unlist(lapply(ld$variants, function(v) v$flagged)))
surv <- surv[surv$item %in% watch, ]
# Short axis labels; any variant name not in SHORT (the drop remedy names the
# item it dropped, which varies with the data) keeps its own name.
shorten <- function(x) { m <- SHORT[x]; ifelse(is.na(m), x, m) }
lev <- vapply(ld$variants, function(v) v$name, "")
surv$treatment <- factor(shorten(surv$treatment), levels = unname(shorten(lev)))
surv$label <- sub("^(\\d+).*", "item \\1", PHQ_LABELS[surv$item])
surv$status <- ifelse(surv$anchor, "anchor (not tested)",
                      ifelse(surv$flagged, "flagged", "tested, not flagged"))
surv$status <- factor(surv$status,
                      levels = c("flagged", "tested, not flagged", "anchor (not tested)"))

p7 <- ggplot(surv, aes(treatment, SIDS, group = item)) +
  facet_wrap(~ label, nrow = 1) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(data = surv[!surv$anchor, ], colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = status, shape = status), size = 2.3, stroke = 0.9) +
  scale_colour_manual(values = c("flagged" = "#b4462f",
                                 "tested, not flagged" = "#9aa6b2",
                                 "anchor (not tested)" = "grey45"), name = NULL, drop = FALSE) +
  scale_shape_manual(values = c("flagged" = 16, "tested, not flagged" = 16,
                                "anchor (not tested)" = 1), name = NULL, drop = FALSE) +
  labs(title = "Treating the local dependence makes the headline effect larger, not smaller",
       subtitle = "Signed expected-score difference by sex, under each standard remedy for local dependence",
       x = NULL, y = "SIDS, item-score points",
       caption = paste("The contaminated anchor was partly absorbing the difference it was meant to measure against.",
                       "Item 7 is the one\nflag that does not survive. Items 4, 5 and 9 survive every remedy,",
                       "the same three that survive every anchor size.\nHollow points are items serving as anchors under that treatment: not tested, so not an estimate of zero.")) +
  base_theme +
  theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1))
ggsave("output/fig7_ld_remedies.png", p7, width = 9, height = 4.6, dpi = 200)

# ---------------------------------------------------------------------------
# Committed tables, same convention as 04.
# Table 4 carries its own power context. A flagged item and the magnitude the
# design can actually resolve were previously in different tables, which left
# the comparison for the reader to notice; every flag here sits close enough to
# the detection floor that not noticing it would change the conclusion.
# `power_at_own_effect` is the share of replications in 06 that detected a
# planted effect of this item's observed size, read off the anchored procedure's
# own detection curve. It is the honest answer to "would this have been found
# again", and for most of these items the answer is no.
pw_a <- err$power[err$power$procedure == "anchored", ]
pw_a <- pw_a[order(pw_a$shift), ]
sids_grid <- abs(approx(map$shift, map$SIDS, xout = pw_a$shift, rule = 2)$y)
floor80 <- abs(approx(map$shift, map$SIDS,
                      xout = err$thresholds[["anchored"]], rule = 2)$y)
own <- abs(grm$es_anchored$SIDS)
power_at_own <- approx(sids_grid, pw_a$detection, xout = own, rule = 2)$y
tested <- !grm$dif_anchored$anchor

t4 <- data.frame(
  item = PHQ_LABELS[grm$dif_anchored$item],
  grm$dif_anchored[, c("X2", "df", "adj_p", "anchor", "flag")],
  SIDS = grm$es_anchored$SIDS, ESSD = grm$es_anchored$ESSD,
  detect_floor_80 = ifelse(tested, round(floor80, 4), NA),
  effect_vs_floor = ifelse(tested, round(own / floor80, 2), NA),
  power_at_own_effect = ifelse(tested, round(power_at_own, 2), NA))
write.csv(t4, "output/table4_grm_dif_by_sex.csv", row.names = FALSE)

cat("\nTable 4, each flag beside the effect size this design can resolve:\n")
cat(sprintf("80%% detection floor for the anchored procedure: %.3f item-score points\n\n",
            floor80))
print(data.frame(item = t4$item,
                 SIDS = round(t4$SIDS, 3),
                 flag = ifelse(t4$anchor, "anchor", ifelse(t4$flag, "DIF", ".")),
                 vs_floor = t4$effect_vs_floor,
                 power = t4$power_at_own_effect), row.names = FALSE)
cat("\nvs_floor is the observed effect divided by the smallest effect the design\n")
cat("resolves 80% of the time. power is the chance a real effect of that exact\n")
cat("size would have been detected at all. Read together: only the largest flag\n")
cat("is near the resolving power of this study, and the rest would not reliably\n")
cat("reappear in a replication of the same size.\n")
write.csv(err$power, "output/table5_dif_detection.csv", row.names = FALSE)
write.csv(err$contamination_rates, "output/table6_anchor_contamination.csv", row.names = FALSE)

cat("Wrote output/fig4_detection_curve.png, output/fig5_anchor_contamination.png\n")
cat("Wrote output/fig6_local_dependence_null.png, output/fig7_ld_remedies.png\n")
cat("Wrote output/table4_grm_dif_by_sex.csv, table5_dif_detection.csv, table6_anchor_contamination.csv\n")
