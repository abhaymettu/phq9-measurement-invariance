# Core logic for the PHQ-9 measurement invariance project.
# The numbered scripts are thin drivers; everything testable lives here.

suppressMessages({
  library(haven); library(dplyr); library(lavaan)
})

# The nine PHQ-9 symptom items, in questionnaire order.
# DPQ100 is the "how difficult" functional-impairment question and is NOT one of
# the nine symptom items. Including it would put a non-symptom indicator on the
# factor and inflate every loading estimate.
PHQ_ITEMS <- paste0("DPQ0", c(10, 20, 30, 40, 50, 60, 70, 80, 90))

# Short labels for figures and tables, matching the standard PHQ-9 wording.
PHQ_LABELS <- c(
  DPQ010 = "1 Little interest",
  DPQ020 = "2 Depressed mood",
  DPQ030 = "3 Sleep problems",
  DPQ040 = "4 Low energy",
  DPQ050 = "5 Appetite change",
  DPQ060 = "6 Feeling bad about self",
  DPQ070 = "7 Trouble concentrating",
  DPQ080 = "8 Psychomotor change",
  DPQ090 = "9 Suicidal ideation"
)

# Somatic vs cognitive-affective split, the standard two-cluster reading of the
# PHQ-9 used in the sex-differential-item-functioning literature.
SOMATIC_ITEMS <- c("DPQ030", "DPQ040", "DPQ050", "DPQ080")

# ---------------------------------------------------------------------------
# Data

#' Merge DPQ_J with DEMO_J, recode the PHQ-9 items, build grouping variables.
#'
#' NHANES codes 7 = refused and 9 = don't know on every DPQ item. Left as-is
#' they would be read as scores of 7 and 9 on a 0 to 3 scale, which would wreck
#' every covariance in the model. Recode first, always.
load_phq9 <- function(dir = "data") {
  dpq  <- read_xpt(file.path(dir, "DPQ_J.xpt"))
  demo <- read_xpt(file.path(dir, "DEMO_J.xpt"))

  dpq <- dpq %>%
    mutate(across(all_of(PHQ_ITEMS), ~ ifelse(.x %in% c(7, 9), NA_real_, as.numeric(.x))))

  demo %>%
    select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, WTMEC2YR, SDMVPSU, SDMVSTRA) %>%
    inner_join(dpq %>% select(SEQN, all_of(PHQ_ITEMS)), by = "SEQN") %>%
    mutate(
      n_answered = rowSums(!is.na(across(all_of(PHQ_ITEMS)))),
      phq9 = ifelse(n_answered == 9, rowSums(across(all_of(PHQ_ITEMS))), NA_real_),
      sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
      # Two bands, split at the median-ish 50, to keep group sizes large enough
      # for a nine-indicator ordinal CFA. More bands means thinner cells and
      # more empty response categories on item 9.
      age_grp = factor(ifelse(RIDAGEYR < 50, "18-49", "50+"),
                       levels = c("18-49", "50+")),
      race = factor(RIDRETH3,
                    levels = c(1, 2, 3, 4, 6, 7),
                    labels = c("Mexican American", "Other Hispanic", "White",
                               "Black", "Asian", "Other/Multi"))
    ) %>%
    as.data.frame()
}

#' Complete-case analytic sample: adults 18+ with all nine items answered.
analytic_sample <- function(dat) {
  dat[dat$RIDAGEYR >= 18 & dat$n_answered == 9, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Models

CFA_MODEL <- paste("dep =~", paste(PHQ_ITEMS, collapse = " + "))

# For ordered indicators the latent response intercepts must be fixed at zero.
# This is not cosmetic. Left at lavaan's multigroup default they are freed in
# every group but the first the moment thresholds are constrained, and a free
# intercept absorbs a threshold shift on the same item exactly. Threshold
# invariance then becomes untestable: the model reports it holds no matter what.
# The simulated-non-invariance control in 03_negative_control.R is what caught
# this. Before the fix, a planted +0.6 threshold shift on one item produced a
# chi-square rise of 9 on 9 df and a score statistic of 5 for the guilty item,
# indistinguishable from noise. After it, the same shift gives a rise of 133 on
# 18 df and the guilty item's score statistic is 134 against 9 for the runner up.
# Location must be carried by the thresholds alone (Wu and Estabrook 2016).
CFA_MODEL_ORDINAL <- paste(c(CFA_MODEL, paste0(PHQ_ITEMS, " ~ 0*1")), collapse = "\n")

#' Fit the one-factor CFA, either ordinal (WLSMV) or continuous (MLR).
#'
#' `equal` is passed to lavaan's group.equal. For the ordinal parameterisation
#' the invariance rungs are loadings, then loadings plus thresholds (Wu and
#' Estabrook); for the continuous one they are loadings, then loadings plus
#' intercepts.
fit_cfa <- function(dat, group = NULL, equal = NULL, ordinal = TRUE,
                    partial = NULL) {
  args <- list(model = if (ordinal) CFA_MODEL_ORDINAL else CFA_MODEL,
               data = dat, group = group,
               group.equal = equal, group.partial = partial,
               std.lv = TRUE)
  if (ordinal) {
    args$ordered <- PHQ_ITEMS
    args$estimator <- "WLSMV"
  } else {
    args$estimator <- "MLR"
  }
  args <- args[!vapply(args, is.null, logical(1))]
  suppressWarnings(do.call(lavaan::cfa, args))
}

FIT_MEASURES <- c("chisq.scaled", "df", "pvalue.scaled",
                  "cfi.scaled", "rmsea.scaled", "srmr")

#' One row of fit indices per model, with change-from-previous-rung columns.
#'
#' Deliberately does NOT report a chi-square difference test as the decision
#' rule. At n of roughly 5,000 the scaled chi-square difference is significant
#' for trivially small misfit, so it rejects invariance essentially always.
#' The conventional practical cutoffs are used instead: a drop in CFI greater
#' than 0.010 or a rise in RMSEA greater than 0.015 flags non-invariance
#' (Cheung and Rensvold 2002; Chen 2007).
fit_table <- function(fits) {
  m <- t(vapply(fits, function(f) fitMeasures(f, FIT_MEASURES), numeric(6)))
  out <- data.frame(model = names(fits), m, row.names = NULL)
  names(out) <- c("model", "chisq", "df", "p", "cfi", "rmsea", "srmr")
  out$d_cfi   <- c(NA, diff(out$cfi))
  out$d_rmsea <- c(NA, diff(out$rmsea))
  out$verdict <- ifelse(is.na(out$d_cfi), "-",
                        ifelse(out$d_cfi < -0.010 | out$d_rmsea > 0.015,
                               "FAIL", "pass"))
  out
}

#' Configural, metric, scalar for a given grouping variable.
invariance_sequence <- function(dat, group, ordinal = TRUE, partial = NULL) {
  rungs <- if (ordinal) {
    list(configural = NULL,
         metric     = c("loadings"),
         scalar     = c("loadings", "thresholds"))
  } else {
    list(configural = NULL,
         metric     = c("loadings"),
         scalar     = c("loadings", "intercepts"))
  }
  lapply(rungs, function(eq)
    fit_cfa(dat, group = group, equal = eq, ordinal = ordinal, partial = partial))
}

# ---------------------------------------------------------------------------
# Partial invariance

#' Rank items by how badly their equality constraints fit, using the
#' univariate Lagrange-multiplier (score) tests on the scalar model's
#' constraints. Nothing is refit, so this is cheap and it is the standard way
#' to localise non-invariance rather than releasing constraints by hand.
#'
#' Returns one row per item: the summed score statistic over all of that item's
#' cross-group equality constraints.
score_test_by_item <- function(fit_scalar, items = PHQ_ITEMS) {
  st <- lavTestScore(fit_scalar, univariate = TRUE)$uni
  pt <- parTable(fit_scalar)
  # Constraint labels .p1. etc. index rows of the parameter table.
  idx <- function(lab) as.integer(sub("^\\.p(\\d+)\\.$", "\\1", lab))
  st$item <- pt$lhs[idx(st$lhs)]
  st$item <- ifelse(st$item %in% items, st$item, pt$rhs[idx(st$lhs)])
  agg <- aggregate(X2 ~ item, data = st[st$item %in% items, ], FUN = sum)
  agg[order(-agg$X2), ]
}

#' Free every cross-group constraint belonging to `items`.
#' For ordinal models that means the item's thresholds and its loading.
partial_labels <- function(items, ordinal = TRUE) {
  loads <- paste0("dep=~", items)
  if (!ordinal) return(c(loads, paste0(items, "~1")))
  c(loads, as.vector(outer(items, 1:3, function(i, t) paste0(i, "|t", t))))
}

# ---------------------------------------------------------------------------
# Consequence for a naive group comparison

#' Decompose the observed group gap in PHQ-9 sum score into the part the model
#' attributes to a real latent difference and the part it attributes to items
#' behaving differently across groups.
#'
#' Done in the continuous (MLR) parameterisation on purpose, because there the
#' model-implied sum score is linear in the parameters and the decomposition is
#' exact and reads in sum-score points:
#'
#'   E[sum | g] = sum_j nu_jg + (sum_j lambda_jg) * kappa_g
#'
#' With group 1 as reference (kappa = 0) the gap splits into
#'   latent part    = (sum lambda) * kappa_2
#'   measurement part = sum_j (nu_j2 - nu_j1)     (zero for constrained items)
#' The measurement part is the bias a naive t-test on sum scores would absorb.
#'
#' `focal` names the group the gap is computed FOR: the returned gap is
#' focal minus the other group.
#'
#' Group order matters and is a live trap. lavaan numbers groups by order of
#' first appearance in the data frame, NOT by factor level order, so group 1 is
#' whichever respondent happens to sit in row 1 of the file. Reading group "1"
#' as the first factor level silently flips the sign of the whole
#' decomposition. Group identity is therefore resolved through
#' lavInspect(fit, "group.label"), never assumed. A test locks this down.
decompose_gap <- function(fit_partial, group_var, dat, focal) {
  labs <- lavInspect(fit_partial, "group.label")
  if (length(labs) != 2) stop("decompose_gap expects exactly two groups")
  if (!focal %in% labs) stop("focal group '", focal, "' is not in the fit")
  f <- as.character(which(labs == focal))       # lavaan group number, focal
  r <- as.character(which(labs != focal))       # reference

  est <- parameterEstimates(fit_partial)
  sum_by_group <- function(rows) {
    e <- est[rows, ]
    vapply(split(e$est, factor(e$group, levels = c(1, 2))), sum, numeric(1))
  }
  sum_lam <- sum_by_group(est$op == "=~")
  sum_nu  <- sum_by_group(est$op == "~1" & est$lhs %in% PHQ_ITEMS)
  kappa   <- sum_by_group(est$op == "~1" & est$lhs == "dep")

  latent_part <- sum_lam[[f]] * kappa[[f]] - sum_lam[[r]] * kappa[[r]]
  meas_part   <- sum_nu[[f]] - sum_nu[[r]]

  obs <- tapply(dat$phq9, as.character(dat[[group_var]]), mean, na.rm = TRUE)
  list(focal = focal, reference = labs[labs != focal],
       observed_gap = as.numeric(obs[[focal]] - obs[[labs[labs != focal]]]),
       model_gap = latent_part + meas_part,
       latent_part = latent_part,
       measurement_part = meas_part,
       latent_sd_units = kappa[[f]] - kappa[[r]])
}
