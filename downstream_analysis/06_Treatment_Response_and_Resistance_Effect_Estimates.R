# =============================================================================
# Figure 5 | Clinical response and resistance effect estimates
# =============================================================================
# Scientific scope: Score distributions, stratified odds ratios and adjusted logistic models.
# Usage: Rscript 06_Treatment_Response_and_Resistance_Effect_Estimates.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: Clinical CSV with MIRROR score/call, response, resistance and covariates.
# Outputs: Panel figures and numerical estimates in MIRROR_OUTPUT_DIR.
# Dependencies: ggplot2; patchwork optional.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# Missing calls remain missing; a threshold is reported only when separable.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))


CLIN_PATH <- Sys.getenv("MIRROR_CLINICAL_CSV", unset = "<PATH_TO_CLINICAL_CSV>")
OUT_DIR   <- Sys.getenv("MIRROR_OUTPUT_DIR", unset = "<PATH_TO_OUTPUT_DIRECTORY>")
BASE_SIZE <- 8

assert_input_file <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || grepl("<[^>]+>", path) || !file.exists(path) || dir.exists(path))
    stop(label, " must be set to an existing input file: ", path, call. = FALSE)
  if (file.access(path, 4L) != 0L)
    stop(label, " is not readable: ", path, call. = FALSE)
  invisible(path)
}
prepare_output_dir <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || grepl("<[^>]+>", path))
    stop(label, " must be set to an output directory.", call. = FALSE)
  if (file.exists(path) && !dir.exists(path))
    stop(label, " points to a file: ", path, call. = FALSE)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path) || file.access(path, 2L) != 0L)
    stop(label, " is not writable: ", path, call. = FALSE)
  invisible(path)
}
assert_input_file(CLIN_PATH, "CLIN_PATH")
prepare_output_dir(OUT_DIR, "OUT_DIR")

PAL_RESP <- c(CR = "#2C5F8A", PR = "#9DC3DF", SD = "#F5C98A", PD = "#A33C2A")
COL_MAIN <- "#31577F"

theme_fig <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text  = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size),
      axis.line  = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      plot.title = element_text(size = base_size, hjust = 0, face = "plain",
                                lineheight = 1.1),
      legend.position = "none",
      plot.margin = margin(4, 6, 4, 4)
    )
}
fmt_p <- function(p) {
  if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
}

#   Multi_PRE_HRD         →  MIRROR call(1 = HRD+, 0 = HRD-)

required_columns <- c("sample_id", "HRD_prob", "Multi_PRE_HRD", "therapy_response",
                      "platinum_resistant", "dataset", "age", "stage_group",
                      "grade_group", "residual_disease", "cohort")
read_clinical <- function(path) {
  for (encoding in c("UTF-8-BOM", "UTF-8", "GB18030")) {
    candidate <- tryCatch(read.csv(path, stringsAsFactors = FALSE,
                                   check.names = FALSE, fileEncoding = encoding),
                          error = function(e) NULL)
    if (!is.null(candidate) && all(required_columns %in% names(candidate))) {
      message("Clinical CSV encoding: ", encoding)
      return(candidate)
    }
  }
  stop("Unable to read the required clinical fields with UTF-8 or GB18030 encoding.")
}
raw <- read_clinical(CLIN_PATH)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) stop("Missing clinical fields: ", paste(missing_columns, collapse = ", "))
if (!nrow(raw)) stop("The clinical CSV has no observations.")
calls <- suppressWarnings(as.numeric(raw$Multi_PRE_HRD))
if (any(!is.na(calls) & !calls %in% 0:1)) stop("Multi_PRE_HRD must contain 0, 1 or NA.")
if (any(is.na(calls) != is.na(raw$Multi_PRE_HRD))) stop("Multi_PRE_HRD contains nonnumeric values.")
resistance <- suppressWarnings(as.numeric(raw$platinum_resistant))
if (any(!is.na(resistance) & !resistance %in% 0:1) ||
    any(is.na(resistance) != is.na(raw$platinum_resistant)))
  stop("platinum_resistant must contain 0, 1 or NA.")
message(sprintf("Loaded %d rows x %d columns.", nrow(raw), ncol(raw)))

map_resp <- c("Complete Remission/Response" = "CR",
              "Partial Remission/Response"  = "PR",
              "Stable Disease"              = "SD",
              "Progressive Disease"         = "PD")
map_set  <- c(Train = "TCGA train", Val = "TCGA validation",
              Internal = "TCGA internal test", CHCAMS = "CHCAMS",
              PLCO = "PLCO", HMUCH = "HMUCH")
map_res  <- c("No Macroscopic disease" = "0-10 mm", "1-10 mm" = "0-10 mm",
              "<=10 mm" = "0-10 mm", "≤10 mm" = "0-10 mm",
              "11-20 mm" = ">10 mm", ">20 mm" = ">10 mm", ">10 mm" = ">10 mm")

clin <- data.frame(
  patient_id   = raw$sample_id,
  mirror_score = as.numeric(raw$HRD_prob),
  mirror_call  = factor(ifelse(is.na(calls), NA_character_,
                                ifelse(calls == 1, "HRD+", "HRD-")),
                        levels = c("HRD-", "HRD+")),
  best_response = factor(unname(map_resp[raw$therapy_response]),
                         levels = c("CR", "PR", "SD", "PD")),
  resistant     = resistance,
  analysis_set  = factor(unname(map_set[raw$dataset]),
                         levels = c("TCGA train", "TCGA validation",
                                    "TCGA internal test", "CHCAMS",
                                    "PLCO", "HMUCH")),
  age        = as.numeric(raw$age),
  figo_stage = factor(ifelse(raw$stage_group %in% c("III", "IV"), "III-IV",
                             ifelse(raw$stage_group %in% c("I", "II"), "I-II", NA)),
                      levels = c("I-II", "III-IV")),
  grade      = factor(ifelse(raw$grade_group %in% "High grade", "G3-G4",
                             ifelse(raw$grade_group %in% "Low grade", "G1-G2", NA)),
                      levels = c("G1-G2", "G3-G4")),
  residual   = factor(unname(map_res[raw$residual_disease]),
                      levels = c("0-10 mm", ">10 mm")),
  cohort     = factor(ifelse(raw$cohort == "TCGA", "TCGA", "External"),
                      levels = c("External", "TCGA")),
  stringsAsFactors = FALSE
)


negative_scores <- clin$mirror_score[clin$mirror_call %in% "HRD-" & is.finite(clin$mirror_score)]
positive_scores <- clin$mirror_score[clin$mirror_call %in% "HRD+" & is.finite(clin$mirror_score)]
if (!length(negative_scores) || !length(positive_scores)) stop("Both MIRROR call groups need finite scores.")
MIRROR_THRESHOLD <- if (max(negative_scores) < min(positive_scores))
  (max(negative_scores) + min(positive_scores)) / 2 else NA_real_
if (!is.finite(MIRROR_THRESHOLD)) warning("Calls are not separable; the score threshold is omitted.")
if (is.finite(MIRROR_THRESHOLD))
  message(sprintf("MIRROR operating threshold = %.4f", MIRROR_THRESHOLD))

db <- clin[!is.na(clin$best_response) & !is.na(clin$mirror_score), ]
if (nrow(db) < 5L || nlevels(droplevels(db$best_response)) < 2L)
  stop("Panel B needs at least five valid cases across two response categories.")
db$xi <- as.numeric(db$best_response)
set.seed(1)
db$xjit <- db$xi - 0.26 + runif(nrow(db), -0.075, 0.075)

half_violin_df <- function(d, grp_col, val_col, width = 0.32, offset = 0.16) {
  lv <- levels(d[[grp_col]])
  out <- lapply(seq_along(lv), function(i) {
    v <- d[[val_col]][d[[grp_col]] == lv[i]]
    v <- v[is.finite(v)]
    if (length(v) < 3 || diff(range(v)) == 0) return(NULL)
    dn   <- stats::density(v, n = 512)
    keep <- dn$x >= min(v) & dn$x <= max(v)
    xx   <- dn$x[keep]
    yy   <- dn$y[keep] / max(dn$y[keep]) * width
    data.frame(grp  = lv[i],
               xpos = c(i + offset + yy, rep(i + offset, length(yy))),
               yval = c(xx, rev(xx)), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  if (is.null(res)) return(data.frame(grp = factor(levels = lv), xpos = numeric(), yval = numeric()))
  res$grp <- factor(res$grp, levels = lv)
  res
}
viol_b <- half_violin_df(db, "best_response", "mirror_score")

kw_b     <- kruskal.test(mirror_score ~ best_response, data = db)
rank_map <- c(PD = 1, SD = 2, PR = 3, CR = 4)     # PD < SD < PR < CR
sp_b <- suppressWarnings(cor.test(db$mirror_score,
                                  rank_map[as.character(db$best_response)],
                                  method = "spearman"))
lab_b <- sprintf("%s\n(n = %d)", levels(db$best_response),
                 as.integer(table(db$best_response)))

p_b <- ggplot() +
  geom_polygon(data = viol_b, aes(x = xpos, y = yval, group = grp, fill = grp),
               colour = NA, alpha = 0.85) +
  geom_point(data = db, aes(x = xjit, y = mirror_score),
             colour = "grey35", size = 0.35, alpha = 0.55) +
  geom_boxplot(data = db, aes(x = xi, y = mirror_score, group = best_response),
               width = 0.13, outlier.shape = NA, fill = NA,
               colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = PAL_RESP) +
  scale_x_continuous(breaks = seq_along(lab_b), labels = lab_b, limits = c(0.4, 4.75)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "MIRROR score",
       title = sprintf("MIRROR score by best response\nKruskal-Wallis %s; Spearman rho = %.3f, %s",
                       fmt_p(kw_b$p.value), unname(sp_b$estimate), fmt_p(sp_b$p.value))) +
  theme_fig()
if (is.finite(MIRROR_THRESHOLD))
  p_b <- p_b + geom_hline(yintercept = MIRROR_THRESHOLD, linetype = "dashed",
                          colour = "grey40", linewidth = 0.3) +
    annotate("text", x = 4.72, y = MIRROR_THRESHOLD, label = "threshold",
             hjust = 1, vjust = -0.5, size = BASE_SIZE / 3.2, colour = "grey30")

dd <- clin[!is.na(clin$resistant) & !is.na(clin$mirror_call), ]
if (!nrow(dd)) stop("Panel D has no evaluable resistance observations.")
dd$res_f <- factor(ifelse(dd$resistant == 1, "yes", "no"), levels = c("yes", "no"))
dd$analysis_set <- droplevels(dd$analysis_set)


make_tab <- function(d) {
  table(factor(as.character(d$mirror_call), levels = c("HRD+", "HRD-")), d$res_f)
}
or_from_table <- function(tab) {
  ft <- try(fisher.test(tab), silent = TRUE)
  ok <- !inherits(ft, "try-error") && is.finite(ft$estimate) &&
    ft$estimate > 0 && all(is.finite(ft$conf.int))
  if (ok) return(c(or = unname(ft$estimate), lo = ft$conf.int[1],
                   hi = ft$conf.int[2], p = ft$p.value))
  a  <- tab[1, 1] + 0.5; b <- tab[1, 2] + 0.5
  cc <- tab[2, 1] + 0.5; dv <- tab[2, 2] + 0.5
  or <- (a * dv) / (b * cc)
  se <- sqrt(1 / a + 1 / b + 1 / cc + 1 / dv)
  c(or = or, lo = or * exp(-1.96 * se), hi = or * exp(1.96 * se),
    p = if (!inherits(ft, "try-error")) ft$p.value else NA_real_)
}

rows <- list(); arr_l <- list()
for (s in levels(dd$analysis_set)) {
  ds  <- dd[dd$analysis_set == s, ]
  tab <- make_tab(ds)

  if (any(rowSums(tab) == 0) || sum(tab[, "yes"]) == 0) {
    message(sprintf("Panel d: skipped uninformative cohort %s (%d cases, %d resistant).",
                    s, sum(tab), sum(tab[, "yes"])))
    next
  }
  est <- or_from_table(tab)
  rows[[s]] <- data.frame(
    label = sprintf("%s\n%.1f%% vs %.1f%%", s,
                    100 * tab[1, 1] / sum(tab[1, ]),
                    100 * tab[2, 1] / sum(tab[2, ])),
    or = est["or"], lo = est["lo"], hi = est["hi"], type = "set",
    stringsAsFactors = FALSE)
  arr_l[[s]] <- matrix(as.numeric(tab), 2, 2)
}
if (!length(arr_l)) stop("Panel D has no informative analysis sets for a pooled odds ratio.")
arr <- array(unlist(arr_l), dim = c(2, 2, length(arr_l)))

mh <- try(mantelhaen.test(arr, correct = FALSE), silent = TRUE)
if (!inherits(mh, "try-error")) {
  mh_or <- unname(mh$estimate); mh_lo <- mh$conf.int[1]
  mh_hi <- mh$conf.int[2];      mh_p  <- mh$p.value
} else mh_or <- mh_lo <- mh_hi <- mh_p <- NA_real_


breslow_day <- function(x, or_hat) {
  K <- dim(x)[3]; X2 <- 0
  if (K < 2L || !is.finite(or_hat) || or_hat <= 0)
    return(list(stat = NA_real_, p = NA_real_))
  for (j in seq_len(K)) {
    mj <- rowSums(x[, , j]); nj <- colSums(x[, , j])
    cf <- c(-mj[1] * nj[1] * or_hat,
            nj[2] - mj[1] + or_hat * (nj[1] + mj[1]),
            1 - or_hat)
    rt <- Re(polyroot(cf))
    lo <- max(0, nj[1] - mj[2]); hi <- min(nj[1], mj[1])
    ok <- rt >= lo - 1e-8 & rt <= hi + 1e-8
    if (!any(ok)) return(list(stat = NA_real_, p = NA_real_))
    aa <- rt[ok][1]
    Va <- 1 / (1 / aa + 1 / (nj[1] - aa) + 1 / (mj[1] - aa) + 1 / (mj[2] - nj[1] + aa))
    X2 <- X2 + (x[1, 1, j] - aa)^2 / Va
  }
  list(stat = X2, p = pchisq(X2, df = K - 1, lower.tail = FALSE))
}
bd <- if (is.na(mh_or)) list(p = NA_real_) else breslow_day(arr, mh_or)

tab_all <- make_tab(dd)
rows[["pooled"]] <- data.frame(
  label = sprintf("Pooled (MH)\n%.1f%% vs %.1f%%",
                  100 * tab_all[1, 1] / sum(tab_all[1, ]),
                  100 * tab_all[2, 1] / sum(tab_all[2, ])),
  or = mh_or, lo = mh_lo, hi = mh_hi, type = "pooled", stringsAsFactors = FALSE)

fd <- do.call(rbind, rows)
fd$label <- factor(fd$label, levels = rev(fd$label))
xr_d <- c(0.05, 4)
fd$lo_c <- pmax(fd$lo, xr_d[1]); fd$hi_c <- pmin(fd$hi, xr_d[2])

p_d <- ggplot(fd, aes(x = or, y = label)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo_c, xmax = hi_c), linewidth = 0.35) +
  geom_point(aes(shape = type), fill = COL_MAIN, colour = COL_MAIN, size = 1.9) +
  scale_shape_manual(values = c(set = 15, pooled = 23)) +
  scale_x_log10(limits = xr_d, breaks = c(0.05, 0.25, 1, 4),
                labels = c("0.05", "0.25", "1", "4")) +
  labs(x = "Odds ratio (HRD+ vs HRD-)", y = NULL,
       title = sprintf("Platinum-resistant recurrence (n = %d)\nMH OR %.2f (95%% CI %.2f-%.2f), %s\nBreslow-Day heterogeneity %s",
                       nrow(dd), mh_or, mh_lo, mh_hi, fmt_p(mh_p), fmt_p(bd$p))) +
  theme_fig() +
  theme(axis.text.y = element_text(hjust = 1, lineheight = 0.95))

clin$cr <- ifelse(is.na(clin$best_response), NA_integer_,
                  as.integer(clin$best_response == "CR"))

tidy_or <- function(fit, keep, labels) {
  s  <- summary(fit)$coefficients
  s  <- s[rownames(s) != "(Intercept)", , drop = FALSE]
  i  <- match(keep, rownames(s))
  data.frame(term = labels,
             or = exp(s[i, "Estimate"]),
             lo = exp(s[i, "Estimate"] - 1.96 * s[i, "Std. Error"]),
             hi = exp(s[i, "Estimate"] + 1.96 * s[i, "Std. Error"]),
             p  = s[i, "Pr(>|z|)"], stringsAsFactors = FALSE)
}

d_cr <- clin[stats::complete.cases(clin[, c("cr", "mirror_call", "age",
                                            "figo_stage", "grade", "residual")]), ]
if (nrow(d_cr) < 20L || length(unique(d_cr$cr)) < 2L ||
    any(vapply(d_cr[c("mirror_call", "figo_stage", "grade", "residual")],
               function(x) nlevels(droplevels(x)) < 2L, logical(1))))
  stop("The complete-response model lacks sufficient complete cases, outcomes or factor levels.")
fit_cr <- glm(cr ~ mirror_call + age + figo_stage + grade + residual,
              data = d_cr, family = binomial())
e_left <- tidy_or(fit_cr,
                  c("mirror_callHRD+", "age", "figo_stageIII-IV",
                    "gradeG3-G4", "residual>10 mm"),
                  c("MIRROR HRD+", "Age (per year)", "FIGO stage III-IV",
                    "Grade G3-G4", "Residual > 10 mm"))

d_rr <- clin[stats::complete.cases(clin[, c("resistant", "mirror_call", "age",
                                            "figo_stage", "grade", "cohort")]), ]
if (nrow(d_rr) < 20L || length(unique(d_rr$resistant)) < 2L ||
    any(vapply(d_rr[c("mirror_call", "figo_stage", "grade", "cohort")],
               function(x) nlevels(droplevels(x)) < 2L, logical(1))))
  stop("The resistance model lacks sufficient complete cases, outcomes or factor levels.")
fit_rr <- glm(resistant ~ mirror_call + age + figo_stage + grade + cohort,
              data = d_rr, family = binomial())
e_right <- tidy_or(fit_rr,
                   c("mirror_callHRD+", "age", "figo_stageIII-IV",
                     "gradeG3-G4", "cohortTCGA"),
                   c("MIRROR HRD+", "Age (per year)", "FIGO stage III-IV",
                     "Grade G3-G4", "Cohort: TCGA"))

forest_e <- function(d, ttl, xr = c(0.15, 6)) {
  d$term <- factor(d$term, levels = rev(d$term))
  d$lo_c <- pmax(d$lo, xr[1]); d$hi_c <- pmin(d$hi, xr[2])
  d$cut_lo <- d$lo < xr[1];    d$cut_hi <- d$hi > xr[2]
  g <- ggplot(d, aes(x = or, y = term)) +
    geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40", linewidth = 0.3) +
    geom_linerange(aes(xmin = lo_c, xmax = hi_c), linewidth = 0.35) +
    geom_point(shape = 15, colour = COL_MAIN, size = 1.9) +
    scale_x_log10(limits = xr, breaks = c(0.25, 1, 4), labels = c("0.25", "1", "4")) +
    labs(x = "Adjusted odds ratio (95% CI)", y = NULL, title = ttl) +
    theme_fig()

  if (any(d$cut_lo)) g <- g + geom_segment(
    data = d[d$cut_lo, ], aes(x = or, xend = xr[1], y = term, yend = term),
    linewidth = 0.35, arrow = arrow(length = unit(1.2, "mm"), type = "closed"))
  if (any(d$cut_hi)) g <- g + geom_segment(
    data = d[d$cut_hi, ], aes(x = or, xend = xr[2], y = term, yend = term),
    linewidth = 0.35, arrow = arrow(length = unit(1.2, "mm"), type = "closed"))
  g
}
p_e1 <- forest_e(e_left,  "Complete response")
p_e2 <- forest_e(e_right, "Resistant recurrence")
title_e <- sprintf("Multivariable logistic models\nCR: n = %d, %d events | Resistance: n = %d, %d events",
                   nrow(d_cr), sum(d_cr$cr), nrow(d_rr), sum(d_rr$resistant))

ggsave(file.path(OUT_DIR, "fig5b_response_distribution.pdf"), p_b,
       width = 85, height = 72, units = "mm")
ggsave(file.path(OUT_DIR, "fig5d_resistance_forest.pdf"), p_d,
       width = 92, height = 72, units = "mm")

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_e <- patchwork::wrap_plots(p_e1, p_e2, nrow = 1) +
    patchwork::plot_annotation(title = title_e,
                               theme = theme(plot.title = element_text(size = BASE_SIZE, hjust = 0, lineheight = 1.1)))
  ggsave(file.path(OUT_DIR, "fig5e_multivariable_forest.pdf"), p_e,
         width = 130, height = 68, units = "mm")
  p_all <- patchwork::wrap_plots(p_b, p_d, p_e1, p_e2, nrow = 2) +
    patchwork::plot_annotation(tag_levels = list(c("b", "d", "e", "")))
  ggsave(file.path(OUT_DIR, "fig5_bde.pdf"), p_all, width = 185, height = 150, units = "mm")
} else {
  message("patchwork is unavailable; exporting the two panel E components separately.")
  ggsave(file.path(OUT_DIR, "fig5e_left_complete_response.pdf"), p_e1,
         width = 68, height = 62, units = "mm")
  ggsave(file.path(OUT_DIR, "fig5e_right_resistant_recurrence.pdf"), p_e2,
         width = 68, height = 62, units = "mm")
}

write.csv(data.frame(kruskal_p = kw_b$p.value, spearman_rho = unname(sp_b$estimate),
                     spearman_p = sp_b$p.value, n = nrow(db)),
          file.path(OUT_DIR, "fig5b_stats.csv"), row.names = FALSE)
write.csv(fd[, c("label", "or", "lo", "hi", "type")],
          file.path(OUT_DIR, "fig5d_stats.csv"), row.names = FALSE)
write.csv(rbind(cbind(model = "Complete response", e_left),
                cbind(model = "Resistant recurrence", e_right)),
          file.path(OUT_DIR, "fig5e_stats.csv"), row.names = FALSE)


cat("\n================ Review ================\n")
cat(sprintf("panel b: n = %d  (CR/PR/SD/PD = %s)\n", nrow(db),
            paste(as.integer(table(db$best_response)), collapse = "/")))
cat(sprintf("panel d: n = %d, MH OR = %.2f (%.2f-%.2f), %s, Breslow-Day %s\n",
            nrow(dd), mh_or, mh_lo, mh_hi, fmt_p(mh_p), fmt_p(bd$p)))
cat(sprintf("panel e: CR n = %d (%d events) | Resistance n = %d (%d events)\n",
            nrow(d_cr), sum(d_cr$cr), nrow(d_rr), sum(d_rr$resistant)))
cat("Output directory: ", normalizePath(OUT_DIR, winslash = "/"), "\n", sep = "")
