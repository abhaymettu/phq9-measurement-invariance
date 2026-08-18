# IRT extension: graded response model and DIF by sex, plus the machinery the
# error-rate stage needs. Companion to R/invariance.R, same rule: the numbered
# scripts are thin drivers, everything testable lives here.
#
# The CFA side of this repo asks whether a set of cross-group equality
# constraints degrades global fit. The IRT side asks a narrower question one
# item at a time: holding the latent trait fixed, do the two groups have
# different probabilities of endorsing this item. The two can disagree, and
# 07_compare.R is where that is checked rather than assumed.

suppressMessages({
  library(mirt)
})

# Item parameters that vary across groups in a graded response model: the
# discrimination and the three category intercepts. Testing all four is a test
# of both uniform and non-uniform DIF. Testing only the intercepts would be a
# uniform-only test and would miss an item whose slope differs.
GRM_PARS <- c("a1", "d1", "d2", "d3")

# lordif needs 1-based category codes; mirt is happy with 0-based.
as_lordif_resp <- function(resp) as.data.frame(lapply(resp, function(x) as.integer(x) + 1L))

# ---------------------------------------------------------------------------
# Fitting

#' Fit a two-group graded response model.
#'
#' `anchors` are the items whose parameters are held equal across groups. They
#' are what puts the two groups' latent scales on a common metric; the group
#' mean and variance are freed so the anchors, not an arbitrary constraint, do
#' the linking. With `anchors = NULL` every item is constrained, which is the
#' all-other-item baseline used by the drop-one test.
grm_fit <- function(resp, group, anchors = NULL, items = names(resp)) {
  inv <- c(if (is.null(anchors)) items else anchors, "free_means", "free_var")
  multipleGroup(as.data.frame(resp), 1, group = as.character(group),
                itemtype = "graded", invariance = inv, verbose = FALSE)
}

#' Single-group GRM parameters, in the a / b (slope-threshold) metric.
#' Used to make the simulation in 06 resemble the real instrument rather than
#' a textbook one.
grm_pars <- function(resp) {
  m <- mirt(as.data.frame(resp), 1, itemtype = "graded", verbose = FALSE)
  p <- coef(m, IRTpars = TRUE, simplify = TRUE)$items
  list(a = p[, "a"], b = p[, c("b1", "b2", "b3"), drop = FALSE])
}

# ---------------------------------------------------------------------------
# DIF

#' Item-level DIF likelihood-ratio tests.
#'
#' Two schemes, and the difference between them is the whole anchor problem:
#'
#' "all_other" constrains every item, then frees one item at a time and tests
#'   the improvement. Every other item is the anchor. If several items really
#'   do have DIF, the anchor is contaminated by them, the common metric is
#'   pulled off true, and clean items inherit an apparent difference. This is
#'   the standard default and it is the one that flags seven of nine here.
#'
#' "anchored" constrains only `anchors`, leaves the rest free, then adds a
#'   constraint to one studied item at a time. The anchor set is chosen up
#'   front by purify_anchors(), so a contaminated item cannot be doing the
#'   linking.
#'
#' Returns one row per tested item. `p` is unadjusted; `adj_p` is
#' Benjamini-Hochberg across the items tested in this call, because nine tests
#' at alpha .05 is a false-positive machine and 06 measures exactly that.
grm_dif <- function(resp, group, anchors = NULL, items = names(resp),
                    p_adjust = "BH") {
  fit <- grm_fit(resp, group, anchors = anchors, items = items)
  tested <- if (is.null(anchors)) items else setdiff(items, anchors)
  d <- DIF(fit, which.par = GRM_PARS,
           scheme = if (is.null(anchors)) "drop" else "add",
           items2test = tested, p.adjust = p_adjust)
  out <- data.frame(item = rownames(d), X2 = d$X2, df = d$df,
                    p = d$p, adj_p = d$adj_p, row.names = NULL)
  # Anchors were never tested. Recording them as NA rather than dropping them
  # keeps the row count at nine, so a caller cannot silently compare a
  # nine-item run against a five-item one.
  miss <- setdiff(items, out$item)
  if (length(miss)) {
    out <- rbind(out, data.frame(item = miss, X2 = NA_real_, df = NA_real_,
                                 p = NA_real_, adj_p = NA_real_))
  }
  out <- out[match(items, out$item), ]
  out$anchor <- out$item %in% (anchors %||% character(0))
  out$flag <- !is.na(out$adj_p) & out$adj_p < 0.05
  rownames(out) <- NULL
  attr(out, "fit") <- fit
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Pick an anchor set by a two-stage purification.
#'
#' Stage 1 is the contaminated all-other-item run. Stage 2 takes the
#' `n_anchor` items with the smallest stage-1 statistic as the anchor. This is
#' the rank-based anchor (Woods 2009): a fixed-size anchor rather than "every
#' item that failed to flag", because when seven of nine flag, the survivors
#' are two items and a two-item anchor is thinner than the contamination it is
#' meant to fix.
#'
#' The choice of n_anchor is not obviously right, which is why 05 sweeps it
#' instead of asserting one number.
purify_anchors <- function(resp, group, n_anchor = 4, items = names(resp)) {
  anchors_from_stage1(grm_dif(resp, group, anchors = NULL, items = items), n_anchor)
}

#' The stage-2 half on its own, for callers that already have a stage-1 run in
#' hand. Sweeping the anchor size re-uses one stage 1 for every size, and each
#' replication in 06 re-uses its own. Refitting the identical baseline once per
#' anchor size would multiply the error-rate stage's runtime for nothing.
anchors_from_stage1 <- function(stage1, n_anchor = 4) {
  stage1$item[order(stage1$X2, decreasing = FALSE)][seq_len(n_anchor)]
}

#' Ordinal logistic regression DIF (Crane / Choi), a genuinely different
#' procedure: it matches on an IRT trait estimate but tests with a regression,
#' and it purifies its own anchor set iteratively. Run alongside the IRT-LR
#' tests so that agreement is evidence rather than one method restating itself.
lordif_dif <- function(resp, group, alpha = 0.01, items = names(resp)) {
  l <- lordif::lordif(as_lordif_resp(resp), as.character(group),
                      criterion = "Chisqr", alpha = alpha,
                      minCell = 5, maxIter = 10)
  data.frame(item = items, chi12 = l$stats$chi12, chi13 = l$stats$chi13,
             chi23 = l$stats$chi23, flag = as.logical(l$flag), row.names = NULL)
}

# ---------------------------------------------------------------------------
# Effect size
#
# A p-value at n = 5,068 answers "is it exactly zero", which is never the
# interesting question. These are expected-score effect sizes: how many raw
# points of this 0-3 item do the two groups differ by at the same latent
# depression level.
#   SIDS  signed difference in expected item score, averaged over the focal
#         group's trait distribution. Reads in item-score points.
#   ESSD  the same difference standardised, so it reads like a Cohen's d.
#   max.D the largest signed gap at any point on the trait scale.

dif_effect_size <- function(fit, items = NULL) {
  es <- as.data.frame(empirical_ES(fit, DIF = TRUE))
  out <- data.frame(item = rownames(es), SIDS = es$SIDS, UIDS = es$UIDS,
                    ESSD = es$ESSD, max_D = es$max.D, row.names = NULL)
  if (!is.null(items)) out <- out[match(items, out$item), ]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Simulation, for the planted-effect half of the error-rate check

#' Simulate two groups from a graded response model.
#'
#' `pars` comes from grm_pars() on the real data, so the simulated instrument
#' has this questionnaire's skew and discrimination rather than round numbers.
#' `shift` moves every threshold of `bad_item` in group B by that many logits:
#' positive means group B needs more of the trait to endorse, i.e. the item
#' under-reports for them. `shift = 0` generates a true null.
#'
#' `bad_item` may name more than one item. That is not a convenience: the
#' anchor-contamination claim is a claim about what happens when SEVERAL items
#' carry DIF at once, and it cannot be tested by planting one.
#'
#' `mean_diff` puts a real latent difference between the groups. It defaults to
#' 0.19 SD, the value the CFA estimated for sex, because a DIF procedure's
#' error rate under a group difference in true trait level is the case that
#' actually arises, and impact is a known way to induce spurious DIF.
sim_grm <- function(pars, n_per_group = 2500, shift = 0, bad_item = 5,
                    mean_diff = 0.19, items = names(pars$a)) {
  stopifnot(all(bad_item >= 1), all(bad_item <= length(pars$a)))
  d_from_b <- function(b, a) -a * b
  make_d <- function(sh) {
    b <- pars$b
    b[bad_item, ] <- b[bad_item, ] + sh
    t(vapply(seq_len(nrow(b)), function(j) d_from_b(b[j, ], pars$a[j]),
             numeric(ncol(b))))
  }
  a <- matrix(pars$a, ncol = 1)
  A <- simdata(a, make_d(0), n_per_group, itemtype = rep("graded", length(pars$a)))
  B <- simdata(a, make_d(shift), n_per_group, itemtype = rep("graded", length(pars$a)),
               mu = mean_diff)
  resp <- as.data.frame(rbind(A, B))
  names(resp) <- items
  list(resp = resp, group = rep(c("A", "B"), each = n_per_group))
}

#' The truth a planted shift corresponds to on the effect-size scale, so the
#' power curve can be read in item-score points and not only in logits.
planted_effect_size <- function(pars, shift, bad_item = 5, n_per_group = 20000,
                                items = names(pars$a)) {
  s <- sim_grm(pars, n_per_group, shift, bad_item, mean_diff = 0, items = items)
  fit <- grm_fit(s$resp, s$group, anchors = items[-bad_item])
  es <- dif_effect_size(fit, items)
  es[bad_item, ]
}
