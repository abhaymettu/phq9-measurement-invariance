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

p5 <- ggplot(cm, aes(n_bad, false_alarm_rate, colour = procedure)) +
  geom_hline(yintercept = err$settings$alpha, linetype = "dashed", colour = "grey30") +
  annotate("text", x = max(cm$n_bad), y = err$settings$alpha, hjust = 1, vjust = -0.7,
           size = 2.9, colour = "grey25", label = "nominal 0.05") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_colour_manual(values = PAL, labels = PROC_LAB, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "A DIF test is not a test of one item. It is a test of the company it keeps.",
       subtitle = "False positives on items where nothing was planted, as more of the other items carry real bias",
       x = "Number of items carrying planted bias (of 9)",
       y = "False positive rate on the innocent items",
       caption = paste("With nothing planted every procedure is conservative, which is why the",
                       "random-split control alone looks reassuring\nand cannot settle this.",
                       "The default anchor degrades as soon as there is real bias for it to be",
                       "contaminated by.")) +
  base_theme
ggsave("output/fig5_anchor_contamination.png", p5, width = 8, height = 5.2, dpi = 200)

# ---------------------------------------------------------------------------
# Committed tables, same convention as 04.
write.csv(data.frame(item = PHQ_LABELS[grm$dif_anchored$item],
                     grm$dif_anchored[, c("X2", "df", "adj_p", "anchor", "flag")],
                     SIDS = grm$es_anchored$SIDS, ESSD = grm$es_anchored$ESSD),
          "output/table4_grm_dif_by_sex.csv", row.names = FALSE)
write.csv(err$power, "output/table5_dif_detection.csv", row.names = FALSE)
write.csv(err$contamination_rates, "output/table6_anchor_contamination.csv", row.names = FALSE)

cat("Wrote output/fig4_detection_curve.png, output/fig5_anchor_contamination.png\n")
cat("Wrote output/table4_grm_dif_by_sex.csv, table5_dif_detection.csv, table6_anchor_contamination.csv\n")
