#!/usr/bin/env Rscript
# Figures and the committed results tables.

source("R/invariance.R")
suppressMessages({library(ggplot2); library(tidyr)})

res <- readRDS("output/invariance_results.rds")
nc  <- readRDS("output/negative_control.rds")
s   <- readRDS("output/analytic_sample.rds")

PAL <- c(WLSMV = "#2c3e50", MLR = "#9aa6b2")
base_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "top")

# ---------------------------------------------------------------------------
# Figure 1. The change in CFI at the scalar step for every grouping variable
# and both estimators, against the 0.010 cutoff, with the random-split null
# band drawn from the negative control so the reader can see what "nothing to
# find" actually looks like on the same axis.
seq_all <- do.call(rbind, res$seq)
sc <- seq_all[seq_all$model == "scalar", ]
sc$group <- factor(sc$group, levels = c("sex", "race", "age"),
                   labels = c("Sex\n(2 groups)", "Race/ethnicity\n(6 groups)",
                              "Age band\n(2 groups)"))
sc$estimator <- factor(sc$estimator, levels = c("WLSMV", "MLR"))
null_band <- range(c(nc$random2$d_cfi_scalar, nc$random6$d_cfi_scalar))
n_splits <- nrow(nc$random2) / 2 + nrow(nc$random6) / 2

p1 <- ggplot(sc, aes(x = group, y = d_cfi, fill = estimator)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = null_band[1], ymax = null_band[2],
           fill = "#7fb3a0", alpha = 0.28) +
  geom_text(data = data.frame(group = factor(levels(sc$group)[1], levels(sc$group))),
            aes(y = null_band[2],
                label = sprintf("full range across %d random splits: nothing to find", n_splits)),
            inherit.aes = FALSE, x = 1, hjust = 0, vjust = -0.9,
            size = 2.9, colour = "#2f6b56") +
  geom_hline(yintercept = -0.010, linetype = "dashed", colour = "grey30") +
  geom_text(data = data.frame(group = factor(levels(sc$group)[3], levels(sc$group))),
            aes(y = -0.0104, label = "non-invariance cutoff, -0.010"),
            inherit.aes = FALSE, x = 3, hjust = 0.5, vjust = 1.6,
            size = 2.9, colour = "grey25") +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%+.4f", d_cfi)),
            position = position_dodge(width = 0.72), vjust = -0.7, size = 3,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL) +
  scale_y_continuous(labels = function(x) sprintf("%+.3f", x),
                     expand = expansion(mult = c(0.04, 0.16))) +
  labs(title = "The PHQ-9 fails scalar invariance across sex, under both estimators",
       subtitle = "Change in CFI when item thresholds or intercepts are forced equal across groups.\nBars below the dashed line mean group mean comparisons are not on a common scale.",
       x = NULL, y = "Change in CFI, scalar vs metric", fill = "Estimator",
       caption = paste0("NHANES 2017-2018 depression screener, adults 18+ with all nine items answered, n = ",
                        nrow(s), ". Unweighted CFA.\nWLSMV treats items as ordered with a polychoric basis; MLR treats them as continuous. Cutoff from Cheung and Rensvold (2002).")) +
  base_theme + theme(panel.grid.major.x = element_blank())
ggsave("output/fig1_invariance_verdict.png", p1, width = 8.2, height = 5.2, dpi = 150)

# ---------------------------------------------------------------------------
# Figure 2. Which items misbehave across sex, somatic ones marked.
sco <- do.call(rbind, lapply(names(res$scores), function(k) {
  if (!grepl("^sex", k)) return(NULL)
  d <- res$scores[[k]]; d$estimator <- sub("^sex ", "", k); d
}))
sco$label <- PHQ_LABELS[sco$item]
sco$cluster <- ifelse(sco$item %in% SOMATIC_ITEMS, "Somatic", "Cognitive-affective")
ord <- aggregate(X2 ~ label, sco, sum)
sco$label <- factor(sco$label, levels = ord$label[order(ord$X2)])

p2 <- ggplot(sco, aes(x = X2, y = label, fill = cluster)) +
  geom_col(width = 0.68) +
  facet_wrap(~estimator, scales = "free_x") +
  scale_fill_manual(values = c(Somatic = "#b4553f", `Cognitive-affective` = "#8d99a6")) +
  labs(title = "Low energy is the item that behaves most differently for women and men",
       subtitle = "Score (Lagrange multiplier) statistic for each item's cross-group equality constraints\nin the scalar model. Larger means the constraint fits worse. Both estimators put a somatic\nitem first, matching the published pattern of sex differential item functioning on the PHQ-9.",
       x = "Score statistic on the item's equality constraints", y = NULL, fill = NULL,
       caption = paste0("NHANES 2017-2018, n = ", nrow(s),
                        ". WLSMV constrains three thresholds and one loading per item; MLR constrains one intercept and one loading,\nso the two panels are not on a comparable scale and are shown with free axes.")) +
  base_theme + theme(panel.grid.major.y = element_blank())
ggsave("output/fig2_offending_items.png", p2, width = 8.6, height = 4.6, dpi = 150)

# ---------------------------------------------------------------------------
# Figure 3. The consequence: how the observed sex gap in sum scores splits.
g <- res$gap$partial
bars <- data.frame(
  part = factor(c("Real latent difference", "Item non-invariance"),
                levels = c("Real latent difference", "Item non-invariance")),
  pts = c(g$latent_part, g$measurement_part))
# ggplot stacks in reverse factor-level order, so the label midpoints must be
# computed in that same order or they land on the wrong segment.
o <- order(bars$part, decreasing = TRUE)
bars$mid[o] <- cumsum(bars$pts[o]) - bars$pts[o] / 2
p3 <- ggplot(bars, aes(x = 1, y = pts, fill = part)) +
  geom_col(width = 0.34) +
  geom_hline(yintercept = g$observed_gap, linetype = "dotted") +
  annotate("text", x = 0.76, y = g$observed_gap, hjust = 1.05, vjust = 0.5, size = 3,
           colour = "grey25",
           label = sprintf("observed gap %+.2f points", g$observed_gap)) +
  geom_text(aes(x = 1.24, y = mid,
                label = sprintf("%+.2f pts\n(%.0f%%)", pts, 100 * pts / g$model_gap)),
            size = 3.3, fontface = "bold", colour = "grey15", lineheight = 0.95) +
  scale_fill_manual(values = c("Real latent difference" = "#2c3e50",
                               "Item non-invariance" = "#b4553f")) +
  scale_x_continuous(limits = c(0.70, 1.40)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  coord_flip() +
  labs(title = sprintf("About %.0f%% of the female-male PHQ-9 gap is measurement, not symptoms",
                       100 * g$measurement_part / g$model_gap),
       subtitle = sprintf("Decomposition of the %+.2f point model-implied gap in mean PHQ-9 sum score under the partially\ninvariant model. The latent difference is %+.2f SD; forcing full invariance would report %+.2f SD.",
                          g$model_gap, g$latent_sd_units, res$gap$full$latent_sd_units),
       x = NULL, y = "Points of the 0 to 27 PHQ-9 sum score", fill = NULL,
       caption = paste0("NHANES 2017-2018, n = ", nrow(s),
                        ", unweighted CFA in the MLR parameterisation where the sum score is linear in the parameters.\nSurvey-weighted gap is ",
                        sprintf("%+.2f", res$weighted$gap), " points, 95% CI [",
                        sprintf("%.2f, %.2f", res$weighted$ci[1], res$weighted$ci[2]),
                        "]. The latent/measurement split is unchanged to three decimals when females are the reference group.")) +
  base_theme + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
ggsave("output/fig3_gap_decomposition.png", p3, width = 8.2, height = 3.6, dpi = 150)

# ---------------------------------------------------------------------------
# Committed tables, so a reviewer can read the numbers without running R.
write.csv(seq_all[, c("group", "estimator", "model", "chisq", "df", "cfi", "rmsea",
                      "srmr", "d_cfi", "d_rmsea", "verdict")],
          "output/table1_invariance_sequence.csv", row.names = FALSE)
sco$somatic <- sco$cluster == "Somatic"
write.csv(sco[order(sco$estimator, -sco$X2), c("estimator", "item", "label", "X2", "somatic")],
          "output/table2_item_score_tests.csv", row.names = FALSE)
write.csv(nc$simulated, "output/table3_simulation_recovery.csv", row.names = FALSE)

cat("Wrote:\n"); cat(paste0("  ", list.files("output"), collapse = "\n"), "\n")
