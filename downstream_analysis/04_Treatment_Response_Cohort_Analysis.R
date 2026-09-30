# =============================================================================
# Figures 4-5 | Cohort-level treatment response
# =============================================================================
# Scientific scope: RECIST composition and treatment response gradients.
# Usage: Rscript 04_Treatment_Response_Cohort_Analysis.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: Clinical CSV with cohort, split, MIRROR scores/calls and therapy_response.
# Outputs: Publication PDF and statistical annotations in MIRROR_OUTPUT_DIR.
# Dependencies: ggplot2, patchwork, scales.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# Historical sample size is audited without preventing new cohorts from running.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260803)

## ============================================================ user settings ==
DATA_FILE <- Sys.getenv("MIRROR_CLINICAL_CSV", unset = "<PATH_TO_CLINICAL_CSV>")
OUT_DIR   <- Sys.getenv("MIRROR_OUTPUT_DIR", unset = "<PATH_TO_OUTPUT_DIRECTORY>")
CSV_ENC   <- "GB18030"

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
assert_input_file(DATA_FILE, "DATA_FILE")
prepare_output_dir(OUT_DIR, "OUT_DIR")

## -------------------------------------------------------------- packages ----
need <- c("ggplot2", "patchwork", "scales")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))
library(ggplot2); library(patchwork); library(scales)

## ---------------------------------------------------- geometry (mm) ---------
BASE         <- 8
SZ           <- BASE / .pt
SHEET_MARGIN <- 6
SHEET_W      <- 175 + 2 * SHEET_MARGIN
ROW_H        <- 60

## ------------------------------------------------- global theme --------------
theme_pub <- function() {
  theme_classic(base_size = BASE, base_family = "sans") +
    theme(
      text              = element_text(size = BASE, family = "sans", colour = "black"),
      line              = element_line(linewidth = 0.25, colour = "black"),
      axis.line         = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = unit(1.2, "pt"),
      axis.title.x      = element_text(size = BASE, margin = margin(t = 1.5)),
      axis.title.y      = element_text(size = BASE, margin = margin(r = 1.5)),
      axis.text.x       = element_text(size = BASE, colour = "black", margin = margin(t = 0.8)),
      axis.text.y       = element_text(size = BASE, colour = "black", margin = margin(r = 0.8)),
      plot.subtitle     = element_text(size = BASE, hjust = 0, lineheight = 1.1,
                                       margin = margin(b = 2)),
      plot.tag          = element_text(size = BASE + 1, face = "bold"),
      plot.tag.position = c(0, 1),
      legend.title      = element_blank(),
      legend.text       = element_text(size = BASE),
      legend.key.height = unit(6, "pt"),
      legend.key.width  = unit(7, "pt"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(-4, 0, -2, 0),
      legend.background = element_blank(),
      panel.spacing     = unit(3, "pt"),
      plot.margin       = margin(4, 4, 2, 2)
    )
}
theme_set(theme_pub())
update_geom_defaults("text", list(size = SZ, family = "sans", colour = "black"))

## ---------------------------------------------------------- colour palette ---
MINUS      <- "\u2212"
LAB_NEG    <- paste0("pHRD", MINUS)
LAB_POS    <- "pHRD+"
COL_POS    <- "#416993"
COL_RECIST <- c(CR = "#2E5E8E", PR = "#8FB3CE", SD = "#EFC48C", PD = "#B94E2A")
TXT_RECIST <- c(CR = "white",   PR = "black",   SD = "black",   PD = "white")

## ------------------------------------------------------------- helpers -------
fmt_p <- function(p, prefix = "p") {
  ifelse(is.na(p), paste0(prefix, " = NA"),
         ifelse(p < 0.001, paste0(prefix, " < 0.001"),
                paste0(prefix, " = ", formatC(p, format = "f", digits = 3))))
}

wilson <- function(k, n) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  phat <- k / n; z <- 1.959964; den <- 1 + z^2 / n
  ctr  <- (phat + z^2 / (2 * n)) / den
  hw   <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - hw), min(1, ctr + hw))
}

half_violin <- function(values, groups, positions, width = 0.34, adjust = 1) {
  out <- list()
  for (g in names(positions)) {
    v <- values[groups == g]; v <- v[is.finite(v)]
    if (length(v) < 3) next
    dd <- density(v, adjust = adjust, from = min(v), to = max(v), n = 256)
    w  <- dd$y / max(dd$y) * width
    out[[g]] <- data.frame(
      x   = c(positions[[g]], positions[[g]] + w, positions[[g]]),
      y   = c(dd$x[1], dd$x, dd$x[length(dd$x)]),
      grp = g
    )
  }
  if (!length(out)) return(data.frame(x = numeric(), y = numeric(), grp = character()))
  do.call(rbind, out)
}

dev_pdf <- function(file, w, h) {
  if (isTRUE(capabilities("cairo")))
    grDevices::cairo_pdf(file, width = w / 25.4, height = h / 25.4, onefile = FALSE)
  else
    grDevices::pdf(file, width = w / 25.4, height = h / 25.4,
                   family = "sans", useDingbats = FALSE, onefile = FALSE)
}

## ============================================================== load data ====
read_any <- function(path, enc) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  encs <- if (nzchar(enc)) c(enc, "UTF-8-BOM", "UTF-8", "GBK") else
    c("UTF-8-BOM", "UTF-8", "GB18030")
  for (e in encs) {
    x <- try(read.csv(path, fileEncoding = e, check.names = FALSE,
                      stringsAsFactors = FALSE,
                      na.strings = c("", " ", "NA", "NaN", "NULL")), silent = TRUE)
    if (!inherits(x, "try-error") && ncol(x) > 1) { message("CSV encoding: ", e); return(x) }
  }
  stop("Cannot parse the CSV.", call. = FALSE)
}
raw <- read_any(DATA_FILE, CSV_ENC)
required_columns <- c("cohort", "split", "chemotherapy", "HRD_prob",
                      "Multi_PRE_HRD", "therapy_response")
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns))
  stop("Missing treatment-response fields: ", paste(missing_columns, collapse = ", "))
if (!nrow(raw)) stop("The clinical CSV has no observations.")
if (!is.numeric(raw$HRD_prob)) stop("HRD_prob must be numeric.")
if (any(!is.na(raw$Multi_PRE_HRD) & !raw$Multi_PRE_HRD %in% 0:1))
  stop("Multi_PRE_HRD must contain only 0, 1 or NA.")

## ------------------------------------------------- TCGA 547 population -------
## Train + Val (all 350) + Internal with chemotherapy == 1 (197)
tcga_all <- raw[raw$cohort == "TCGA", ]
d <- tcga_all[
  tcga_all$split %in% c("Train", "Val") |
    (tcga_all$split == "Internal" &
       !is.na(tcga_all$chemotherapy) & tcga_all$chemotherapy == 1),
]
cat(sprintf("TCGA analysis population: n = %d\n", nrow(d)))
if (!nrow(d)) stop("The selected TCGA population is empty.")
if (nrow(d) != 547L) warning("Observed n = ", nrow(d), "; historical figure used n = 547.")

## ------------------------------------------------ derived variables ----------
d$call <- factor(d$Multi_PRE_HRD, levels = c(0, 1), labels = c(LAB_NEG, LAB_POS))

negative_scores <- d$HRD_prob[d$Multi_PRE_HRD %in% 0 & is.finite(d$HRD_prob)]
positive_scores <- d$HRD_prob[d$Multi_PRE_HRD %in% 1 & is.finite(d$HRD_prob)]
if (!length(negative_scores) || !length(positive_scores)) stop("Both call groups need finite scores.")
if (max(negative_scores) >= min(positive_scores)) stop("Calls are not separable by one threshold.")
THR <- (max(negative_scores) + min(positive_scores)) / 2

RESP_LV  <- c("Complete Remission/Response", "Partial Remission/Response",
              "Stable Disease", "Progressive Disease")
RESP_ORD <- c("CR", "PR", "SD", "PD")

poolB <- d[!is.na(d$therapy_response) & d$therapy_response %in% RESP_LV, ]
poolB <- poolB[is.finite(poolB$HRD_prob) & !is.na(poolB$Multi_PRE_HRD), , drop = FALSE]
if (nrow(poolB) < 8L || length(unique(poolB$Multi_PRE_HRD)) < 2L ||
    length(unique(poolB$therapy_response)) < 2L)
  stop("Response panels need at least eight finite-score patients, both HRD calls and two response levels.")
poolB$resp <- factor(poolB$therapy_response, levels = RESP_LV, labels = RESP_ORD)
poolB$CR   <- as.integer(poolB$resp == "CR")

cat(sprintf("RECIST pool: n = %d | threshold = %.4f\n", nrow(poolB), THR))
print(table(poolB$resp))
cat("\nHRD call in RECIST pool:\n"); print(table(poolB$call))

## Panel A  |  RECIST stacked bars + chi-square CR vs non-CR
tabA      <- as.data.frame(table(call = poolB$call, resp = poolB$resp))
tabA$N    <- as.numeric(tapply(tabA$Freq, tabA$call, sum)[as.character(tabA$call)])
tabA$prop <- tabA$Freq / tabA$N
tabA$resp <- factor(tabA$resp, levels = RESP_ORD)
tabA$lab  <- ifelse(tabA$prop >= 0.04, sprintf("%.1f%%", 100 * tabA$prop), "")

nA   <- table(poolB$call)
xA_l <- sprintf("%s\n(n = %d)", levels(poolB$call), as.integer(nA[levels(poolB$call)]))

## Chi-square: CR vs non-CR
crtab  <- table(
  factor(poolB$call, levels = c(LAB_NEG, LAB_POS)),
  factor(ifelse(poolB$CR == 1, "CR", "non-CR"), levels = c("CR", "non-CR"))
)
chi_cr  <- chisq.test(crtab, correct = FALSE)
p_chi   <- chi_cr$p.value
chi_lab <- paste0("CR versus non\u2013CR\nChi-squared test, ", fmt_p(p_chi))

pa <- ggplot(tabA, aes(x = call, y = prop, fill = resp)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.56,
           colour = "white", linewidth = 0.2) +
  geom_text(aes(label = lab, colour = resp),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = SZ, lineheight = 1) +
  ## bracket
  annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05,
           linewidth = 0.35, colour = "black") +
  annotate("segment", x = 1, xend = 1, y = 1.01, yend = 1.05,
           linewidth = 0.35, colour = "black") +
  annotate("segment", x = 2, xend = 2, y = 1.01, yend = 1.05,
           linewidth = 0.35, colour = "black") +
  annotate("text", x = 1.5, y = 1.08, label = chi_lab,
           size = SZ, hjust = 0.5, vjust = 0, lineheight = 1.15) +
  scale_fill_manual(values = COL_RECIST, breaks = RESP_ORD) +
  scale_colour_manual(values = TXT_RECIST, guide = "none") +
  scale_x_discrete(labels = xA_l) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.28),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Proportion of patients", tag = "A") +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")

## Panel B  |  MIRROR score by best response
posB        <- setNames(seq_along(RESP_ORD), RESP_ORD)
poolB$xbase <- posB[as.character(poolB$resp)]
poolB$xjit  <- poolB$xbase - 0.26 + runif(nrow(poolB), -0.08, 0.08)
vioB        <- half_violin(poolB$HRD_prob, as.character(poolB$resp), posB, width = 0.33)
vioB$resp   <- factor(vioB$grp, levels = RESP_ORD)

kwB <- kruskal.test(HRD_prob ~ resp, data = poolB)
nB  <- table(poolB$resp)

pb <- ggplot() +
  geom_polygon(data = vioB,
               aes(x = x, y = y, group = resp, fill = resp),
               colour = "grey35", linewidth = 0.2, alpha = 0.9) +
  geom_point(data = poolB,
             aes(x = xjit, y = HRD_prob, colour = resp),
             size = 0.7, alpha = 0.45, stroke = 0) +
  geom_boxplot(data = poolB,
               aes(x = xbase - 0.02, y = HRD_prob, group = resp),
               width = 0.12, outlier.shape = NA, linewidth = 0.25,
               fill = "white", colour = "grey15") +
  geom_hline(yintercept = THR, linetype = "dashed",
             linewidth = 0.25, colour = "grey30") +
  scale_fill_manual(values = COL_RECIST, guide = "none") +
  scale_colour_manual(values = COL_RECIST, guide = "none") +
  scale_x_continuous(
    breaks = posB,
    limits = c(0.5, length(posB) + 0.5),
    labels = sprintf("%s\n(n = %d)", RESP_ORD, as.integer(nB[RESP_ORD])),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = NULL, y = "MIRROR score",
       subtitle = paste0("TCGA only\nKruskal\u2013Wallis test, ", fmt_p(kwB$p.value)),
       tag = "B")

## Panel C  |  CR rate across MIRROR-score quartiles (Cochran-Armitage trend)
qb      <- quantile(poolB$HRD_prob, probs = seq(0, 1, 0.25), na.rm = TRUE)
if (anyDuplicated(qb)) stop("Quartile boundaries are tied; four score bins are unavailable.")
poolB$q <- cut(poolB$HRD_prob, breaks = qb, include.lowest = TRUE,
               labels = c("Q1", "Q2", "Q3", "Q4"))

tabC <- data.frame(
  q = factor(c("Q1", "Q2", "Q3", "Q4"), levels = c("Q1", "Q2", "Q3", "Q4")),
  k = as.numeric(tapply(poolB$CR, poolB$q, sum,  na.rm = TRUE)),
  n = as.numeric(tapply(poolB$CR, poolB$q, length))
)
tabC$rate <- tabC$k / tabC$n
ciC       <- t(mapply(wilson, tabC$k, tabC$n))
tabC$lo   <- ciC[, 1]; tabC$hi <- ciC[, 2]
tabC$xlab <- sprintf("%s\n%.2f\u2013%.2f", tabC$q, qb[-5], qb[-1])
trendC    <- prop.trend.test(tabC$k, tabC$n)

pc <- ggplot(tabC, aes(x = q, y = rate)) +
  geom_col(width = 0.58, fill = COL_POS, alpha = 0.9) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.25) +
  geom_text(aes(y = hi + 0.03, label = sprintf("%.1f%%", 100 * rate)),
            size = SZ, vjust = 0) +
  scale_x_discrete(labels = tabC$xlab) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "MIRROR score quartile", y = "Complete response rate",
       subtitle = paste0("TCGA only\nCochran\u2013Armitage trend test, ",
                         fmt_p(trendC$p.value)),
       tag = "C")

## Assemble and export
fig <- pa + pb + pc +
  plot_layout(widths = c(1, 1.4, 1))

out_pdf <- file.path(OUT_DIR, "Fig_treatment_response_n547.pdf")
dev_pdf(out_pdf, w = SHEET_W, h = ROW_H + 2 * SHEET_MARGIN)
print(fig)
dev.off()
message("Saved: ", out_pdf)
