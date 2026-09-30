# =============================================================================
# Figure 2 | Genomic scar and MIRROR concordance (legacy)
# =============================================================================
# Scientific scope: Paired assay comparisons with classification and provenance audits.
# Usage: Rscript 02_Genomic_Scar_Concordance_Legacy.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: Clinical CSV with genomic scar, MIRROR calls and reference labels.
# Outputs: Figures, numerical tables and audit records in MIRROR_OUTPUT_DIR.
# Dependencies: ggplot2, dplyr, patchwork; ragg and svglite optional.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# Restored from the supplied .bak version; review historical assumptions.
# Optional positional arguments: clinical.csv output_dir n_boot scar_cut.
# =============================================================================

if (.Platform$OS.type == "windows") invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)

CSV_PATH <- if (length(args) >= 1L) args[1] else Sys.getenv("MIRROR_CLINICAL_CSV", unset = "<PATH_TO_CLINICAL_CSV>")
OUT_DIR <- if (length(args) >= 2L) args[2] else Sys.getenv("MIRROR_OUTPUT_DIR", unset = "<PATH_TO_OUTPUT_DIRECTORY>")

N_BOOT <- if (length(args) >= 3L) as.integer(args[3]) else 1000L
SCAR_CUT <- if (length(args) >= 4L) as.numeric(args[4]) else 42
if (!is.finite(N_BOOT) || N_BOOT < 100L) stop("n_boot must be at least 100.")
if (!is.finite(SCAR_CUT)) stop("scar_cut must be finite.")
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
assert_input_file(CSV_PATH, "CSV_PATH")
prepare_output_dir(OUT_DIR, "OUT_DIR")
SEED <- 2024L
set.seed(SEED)

# Preserve the historical HRD-call and assay palettes for reproducibility.
COL_NEG <- "#E0995E"
COL_POS <- "#5B87AE"
COL_HRD <- c(neg = COL_NEG, pos = COL_POS)
ASSAYS <- c("MIRROR", "Genomic HRD score")
COL_ASSAY <- c("MIRROR" = "#9DB9D2", "Genomic HRD score" = "#E3C39A")
COL_ASSAY_LINE <- c("MIRROR" = "#3E6E9B", "Genomic HRD score" = "#B57A38")
SPLITS <- c("TCGA-train", "TCGA-val")
ALL_SPLITS <- c(SPLITS, "TCGA-internal")
FONT <- "sans"
BASE <- 8
TS <- BASE / (72.27 / 25.4)
SMALL <- 6.5 / (72.27 / 25.4)
theme_set(theme_classic(base_size = BASE, base_family = FONT) +
            theme(
              text = element_text(colour = "black"),
              axis.text = element_text(size = 7, colour = "black"),
              axis.title = element_text(size = 8),
              axis.line = element_line(linewidth = 0.25),
              axis.ticks = element_line(linewidth = 0.25),
              axis.ticks.length = grid::unit(0.8, "mm"),
              axis.title.x = element_text(margin = margin(t = 2)),
              axis.title.y = element_text(margin = margin(r = 2)),
              panel.grid = element_blank(),
              panel.spacing = grid::unit(2.3, "mm"),
              strip.background = element_blank(),
              strip.text = element_text(size = 8, margin = margin(b = 2)),
              plot.title = element_text(size = 10, face = "bold", hjust = 0,
                                        margin = margin(b = 2)),
              plot.margin = margin(1.4, 1.4, 1.4, 1.4, "mm"),
              legend.title = element_blank(), legend.text = element_text(size = 7),
              legend.key.size = grid::unit(2.8, "mm"),
              legend.margin = margin(0, 0, 0, 0),
              legend.box.spacing = grid::unit(0.7, "mm")
            ))

# ------------------------------ Statistical helpers -------------------------
bind_rows0 <- function(xs) dplyr::bind_rows(xs)
fmt_p <- function(p) vapply(p, function(x) {
  if (!is.finite(x)) return("P = NA")
  if (x < 0.001) return("P < 0.001")
  sprintf("P = %.3f", x)
}, character(1))
stars <- function(p) ifelse(!is.finite(p), "", ifelse(p < .001, "***",
                                                      ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
quant_ci <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(.025, .975)))
}
safe_median <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
safe_wilcox <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  if (length(unique(c(a, b))) == 1L) return(1)
  suppressWarnings(wilcox.test(a, b, alternative = "two.sided",
                               exact = FALSE, correct = TRUE)$p.value)
}
binary_ok <- function(x) is.finite(x) & x %in% c(0, 1)
auc_fun <- function(score, label) {
  ok <- is.finite(score) & binary_ok(label)
  score <- score[ok]; label <- label[ok]
  m <- sum(label == 1); n <- sum(label == 0)
  if (m == 0 || n == 0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[label == 1]) - m * (m + 1) / 2) / (m * n)
}
# Each unique threshold is processed as one block. Tied positives/negatives
# therefore cannot gain an artificial ordering in the ROC curve.
roc_fun <- function(score, label) {
  ok <- is.finite(score) & binary_ok(label)
  score <- score[ok]; label <- label[ok]
  m <- sum(label == 1); n <- sum(label == 0)
  if (!m || !n) return(data.frame(fpr = numeric(), tpr = numeric()))
  o <- order(score, decreasing = TRUE); s <- score[o]; y <- label[o]
  end <- cumsum(rle(s)$lengths)
  data.frame(fpr = c(0, cumsum(y == 0)[end] / n),
             tpr = c(0, cumsum(y == 1)[end] / m))
}
# Linear interpolation through threshold blocks matches the tie-aware AUC.
tpr_at <- function(r, grid) {
  if (!nrow(r)) return(rep(NA_real_, length(grid)))
  # findInterval chooses the last point of an exact vertical segment while
  # interpolation approaches its first point from the left.
  lo <- pmax(1L, pmin(nrow(r), findInterval(grid, r$fpr)))
  hi <- pmin(lo + 1L, nrow(r))
  dx <- r$fpr[hi] - r$fpr[lo]
  fraction <- ifelse(dx > 0, (grid - r$fpr[lo]) / dx, 0)
  pmin(1, pmax(0, r$tpr[lo] + fraction * (r$tpr[hi] - r$tpr[lo])))
}
wilson <- function(k, n) {
  if (!is.finite(n) || n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(.975); p <- k / n; den <- 1 + z^2 / n
  mid <- (p + z^2 / (2 * n)) / den
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, mid - h), min(1, mid + h))
}
mcnemar_exact <- function(a, b) {
  ok <- binary_ok(a) & binary_ok(b); a <- a[ok]; b <- b[ok]
  a_only <- sum(a == 1 & b == 0); b_only <- sum(a == 0 & b == 1)
  p <- if (!length(a)) NA_real_ else if (a_only + b_only == 0L) 1 else
    binom.test(a_only, a_only + b_only, p = .5, alternative = "two.sided")$p.value
  list(p = p, mirror_only = a_only, scar_only = b_only, n = length(a))
}
# Paired, two-sided DeLong comparison; ties contribute one half.
delong_test <- function(s1, s2, label) {
  ok <- is.finite(s1) & is.finite(s2) & binary_ok(label)
  s1 <- s1[ok]; s2 <- s2[ok]; label <- label[ok]
  m <- sum(label == 1); n <- sum(label == 0)
  a1 <- auc_fun(s1, label); a2 <- auc_fun(s2, label); dif <- a1 - a2
  ans <- list(auc_mirror = a1, auc_genomic = a2, difference = dif,
              se = NA_real_, z = NA_real_, p = NA_real_, lo = NA_real_,
              hi = NA_real_, n_pos = m, n_neg = n)
  if (m < 2L || n < 2L) return(ans)
  psi <- function(x, y) outer(x, y, function(a, b) (a > b) + .5 * (a == b))
  p1 <- psi(s1[label == 1], s1[label == 0])
  p2 <- psi(s2[label == 1], s2[label == 0])
  v10 <- cbind(rowMeans(p1), rowMeans(p2))
  v01 <- cbind(colMeans(p1), colMeans(p2))
  cv <- cov(v10) / m + cov(v01) / n
  variance <- max(0, cv[1, 1] + cv[2, 2] - 2 * cv[1, 2])
  ans$se <- sqrt(variance)
  if (variance <= .Machine$double.eps) {
    ans$z <- if (abs(dif) < 1e-12) 0 else sign(dif) * Inf
    ans$p <- if (abs(dif) < 1e-12) 1 else 0
  } else {
    ans$z <- dif / ans$se
    ans$p <- 2 * pnorm(-abs(ans$z))
  }
  ans$lo <- dif - qnorm(.975) * ans$se
  ans$hi <- dif + qnorm(.975) * ans$se
  ans
}
# Jonckheere-Terpstra, two-sided normal approximation, with score-tie correction.
jt_test <- function(x, grp) {
  ok <- is.finite(x) & !is.na(grp); x <- x[ok]; grp <- droplevels(factor(grp[ok], ordered = TRUE))
  gs <- split(x, grp); gs <- gs[lengths(gs) > 0L]
  ns <- lengths(gs); N <- sum(ns); k <- length(gs)
  ans <- list(J = NA_real_, z = NA_real_, p = NA_real_, variance = NA_real_, n = N)
  if (k < 2L || N < 3L) return(ans)
  J <- 0
  for (i in seq_len(k - 1L)) for (j in seq.int(i + 1L, k)) {
    J <- J + sum(outer(gs[[i]], gs[[j]], "<")) + .5 * sum(outer(gs[[i]], gs[[j]], "=="))
  }
  ts <- as.numeric(table(x))
  variance <- (N * (N - 1) * (2 * N + 5) -
                 sum(ns * (ns - 1) * (2 * ns + 5)) -
                 sum(ts * (ts - 1) * (2 * ts + 5))) / 72 +
    sum(ns * (ns - 1) * (ns - 2)) * sum(ts * (ts - 1) * (ts - 2)) /
    (36 * N * (N - 1) * (N - 2)) +
    sum(ns * (ns - 1)) * sum(ts * (ts - 1)) / (8 * N * (N - 1))
  mu <- (N^2 - sum(ns^2)) / 4
  ans$J <- J; ans$variance <- variance
  if (variance > 0) {
    ans$z <- (J - mu) / sqrt(variance); ans$p <- 2 * pnorm(-abs(ans$z))
  } else if (abs(J - mu) < 1e-12) { ans$z <- 0; ans$p <- 1 }
  ans
}
# Focused invariants: tie-aware ROC agrees with rank AUC; identical paired
# scores have zero AUC difference; empty estimates return NA.
tied <- roc_fun(c(.8, .8, .3, .3), c(1, 0, 1, 0))
trap_auc <- sum(diff(tied$fpr) * (head(tied$tpr, -1) + tail(tied$tpr, -1)) / 2)
stopifnot(abs(trap_auc - .5) < 1e-12, abs(auc_fun(c(.8,.8,.3,.3), c(1,0,1,0)) - .5) < 1e-12,
          is.na(auc_fun(1, 1)), all(is.na(wilson(0, 0))),
          delong_test(c(.1,.2,.8,.9), c(.1,.2,.8,.9), c(0,0,1,1))$p == 1)

vertical_check <- roc_fun(c(.9, .8, .7, .6), c(0, 1, 1, 0))
stopifnot(isTRUE(all.equal(tpr_at(vertical_check, c(.25, .5, .75)), c(0, 1, 1))))
# ---------------------------------- Data ------------------------------------
# Source export uses GB18030 (a superset of GBK). A UTF-8 BOM is also
# tolerated. Other clinical columns are not altered or re-exported.
raw <- read.csv(CSV_PATH, check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "GB18030", na.strings = c("", "NA", "NaN"))
names(raw)[1] <- sub("^\ufeff", "", names(raw)[1])
required <- c("sample_id", "cohort", "dataset", "HRD_prob", "Multi_PRE_HRD",
              "HRD-score", "ai1", "lst1", "hrd-loh", "HRD_label", "has_HRD_label", "HRD_gene")
if (length(setdiff(required, names(raw)))) stop("Missing input columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
to_num <- function(x) suppressWarnings(as.numeric(x))
d <- raw %>% filter(cohort == "TCGA") %>%
  transmute(sample_id = as.character(sample_id), dataset = as.character(dataset),
            mirror_p = to_num(HRD_prob), mirror_call = to_num(Multi_PRE_HRD),
            scar = to_num(.data[["HRD-score"]]), TAI = to_num(ai1), LST = to_num(lst1),
            LOH = to_num(.data[["hrd-loh"]]), label = to_num(HRD_label),
            has_label = to_num(has_HRD_label), gene = ifelse(is.na(HRD_gene), "", HRD_gene))
if (anyDuplicated(d$sample_id)) stop("Duplicate TCGA sample_id values: resolve before patient-level analysis.")
if (any(is.finite(d$mirror_call) & !d$mirror_call %in% c(0, 1))) stop("Multi_PRE_HRD contains values other than 0/1.")
if (any(is.finite(d$label) & !d$label %in% c(0, 1))) stop("HRD_label contains values other than 0/1.")
audit_input <- d %>% group_by(dataset) %>%
  summarise(n_source = n(), n_prob_missing = sum(!is.finite(mirror_p)),
            n_call_missing = sum(!binary_ok(mirror_call)), n_scar_missing = sum(!is.finite(scar)),
            n_reference_0 = sum(has_label == 1 & label == 0, na.rm = TRUE),
            n_reference_1 = sum(has_label == 1 & label == 1, na.rm = TRUE), .groups = "drop")
d <- d %>% filter(is.finite(mirror_p), binary_ok(mirror_call), dataset %in% c("Train", "Val", "Internal")) %>%
  mutate(scar_call = ifelse(is.finite(scar), as.integer(scar >= SCAR_CUT), NA_integer_),
         split = factor(dataset, levels = c("Train", "Val", "Internal"), labels = ALL_SPLITS))
if (any(d$mirror_p < 0 | d$mirror_p > 1)) stop("HRD_prob outside [0, 1].")
main <- d %>% filter(split %in% SPLITS) %>% mutate(split = factor(split, levels = SPLITS))
internal <- d %>% filter(split == "TCGA-internal") %>% mutate(split = droplevels(split))
lab <- main %>% filter(has_label == 1, binary_ok(label))
paired <- lab %>% filter(is.finite(scar))
cor_set <- main %>% filter(is.finite(scar))
if (!nrow(lab)) stop("No labelled Train/Val patients with valid MIRROR predictions.")
negative_max <- if (any(d$mirror_call == 0)) max(d$mirror_p[d$mirror_call == 0]) else NA_real_
positive_min <- if (any(d$mirror_call == 1)) min(d$mirror_p[d$mirror_call == 1]) else NA_real_
separable <- is.finite(negative_max) && is.finite(positive_min) && negative_max < positive_min
THR <- if (separable) (negative_max + positive_min) / 2 else NA_real_
threshold_audit <- data.frame(scar_cut = SCAR_CUT, mirror_negative_max = negative_max,
                              mirror_positive_min = positive_min, mirror_plot_cut = THR, calls_separable = separable,
                              call_source = "Multi_PRE_HRD preserved; inferred MIRROR cut is for plotting only")
if (!separable) warning("Existing MIRROR calls are not separable by one threshold; MIRROR threshold lines are omitted.")
cat("Labelled Train/Val:", nrow(lab), "; HR intact:", sum(lab$label == 0),
    "; HR deficient:", sum(lab$label == 1), "; paired scar:", nrow(paired), "\n")
cat("Genomic HRD cut:", SCAR_CUT, "; inferred MIRROR plotting cut:", THR, "\n")

# Event regexes and mechanism precedence are retained from the supplied script.
# Annotations are NOT reinterpreted as new ground-truth labels.
EV <- c("BRCA1 methylation" = "hypermethelation",
        "BRCA1/2 germline" = "germline", "BRCA1 somatic" = "BRCA1 somatic",
        "BRCA2 somatic" = "BRCA2 somatic", "EMSY amplification" = "EMSY",
        "PTEN loss" = "PTEN", "FANC family" = "FANC|FNAC|C19ORF40|PALB2",
        "RAD50/51 family" = "RAD5", "ATM/ATR/CHEK" = "ATM|ATR|CHEK")
ev_mat <- matrix(vapply(EV, function(p) grepl(p, lab$gene), logical(nrow(lab))), nrow = nrow(lab), ncol = length(EV), dimnames = list(NULL, names(EV)))
lab$n_ev_raw <- rowSums(ev_mat)
lab$n_ev <- pmin(lab$n_ev_raw, 3L)
# Preserve the source analysis: reference HR-intact cases are group 0.
# Any contradictory event annotations are separately audited.
lab$n_ev[lab$label == 0] <- 0
lab$n_ev_f <- factor(lab$n_ev, levels = 0:3, labels = c("0", "1", "2", "3+"))
mech <- rep(NA_character_, nrow(lab)); intact <- lab$label == 0
mech[!intact & grepl("germline", lab$gene)] <- "Germline"
mech[is.na(mech) & !intact & grepl("somatic mutation|Mutation of|ATR mutation", lab$gene)] <- "Somatic"
mech[is.na(mech) & !intact & grepl("hypermethelation", lab$gene)] <- "Epigenetic"
mech[is.na(mech) & !intact & grepl("HOMDEL|Amplification", lab$gene)] <- "Copy number"
mech[intact] <- "HR intact"
lab$mech <- factor(mech, levels = c("HR intact", "Germline", "Somatic", "Epigenetic", "Copy number"))
mechanism_audit <- lab %>% filter(is.na(mech) | (label == 0 & n_ev_raw > 0)) %>%
  select(sample_id, label, gene, mech, n_ev_raw)
if (nrow(mechanism_audit)) warning("Some mechanism/event annotations need review; see audit_mechanism.csv.")
write.csv(data.frame(event = names(EV), pattern = unname(EV), n = colSums(ev_mat)),
          file.path(OUT_DIR, "audit_event_definitions.csv"), row.names = FALSE)

# ------------------------------- Plot helpers -------------------------------
half_violin <- function(y, grp, lv, width = .29, off = .075) {
  out <- list()
  for (i in seq_along(lv)) {
    v <- y[!is.na(grp) & grp == lv[i]]; v <- v[is.finite(v)]
    if (length(v) < 5L || diff(range(v)) == 0) next
    de <- density(v, n = 128, from = min(v), to = max(v))
    w <- de$y / max(de$y) * width
    out[[length(out) + 1L]] <- data.frame(g = lv[i],
                                          x = c(i + off + w, rep(i + off, length(w))),
                                          y = c(de$x, rev(de$x)))
  }
  if (!length(out)) return(data.frame(g = character(), x = numeric(), y = numeric()))
  bind_rows0(out)
}
add_mirror_hline <- function(p) {
  if (is.finite(THR)) p + geom_hline(yintercept = THR, linetype = 2,
                                     linewidth = .25, colour = "grey55") else p
}
raincloud <- function(dat, gvar, pos_lv, xlab, panel) {
  lv <- levels(dat[[gvar]])
  dd <- data.frame(g = dat[[gvar]], y = dat$mirror_p)
  dd <- dd[!is.na(dd$g) & is.finite(dd$y), ]
  dd$xi <- as.integer(dd$g); dd$xj <- dd$xi - .23 + runif(nrow(dd), -.07, .07)
  dd$cls <- ifelse(as.character(dd$g) %in% pos_lv, "pos", "neg")
  hv <- half_violin(dd$y, as.character(dd$g), lv)
  hv$cls <- ifelse(hv$g %in% pos_lv, "pos", "neg")
  p <- ggplot() +
    geom_point(data = dd, aes(xj, y, colour = cls), size = .32, alpha = .35) +
    geom_polygon(data = hv, aes(x, y, group = g, fill = cls, colour = cls),
                 alpha = .28, linewidth = .25) +
    geom_boxplot(data = dd, aes(xi, y, group = xi, colour = cls),
                 width = .16, outlier.shape = NA, linewidth = .3, fill = "white") +
    scale_colour_manual(values = COL_HRD) + scale_fill_manual(values = COL_HRD) +
    scale_x_continuous(breaks = seq_along(lv), labels = lv,
                       limits = c(.48, length(lv) + .48), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, .25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1.15), clip = "off") +
    labs(x = xlab, y = "MIRROR HRD score", title = panel) +
    theme(legend.position = "none")
  add_mirror_hline(p)
}

# A: two-sided Wilcoxon rank-sum tests, BH correction across four mechanisms.
a_stat <- bind_rows0(lapply(levels(lab$mech), function(g) {
  x <- lab$mirror_p[!is.na(lab$mech) & lab$mech == g]
  data.frame(mech = g, n = length(x), median = safe_median(x),
             p = if (g == "HR intact") NA_real_ else
               safe_wilcox(x, lab$mirror_p[!is.na(lab$mech) & lab$mech == "HR intact"]))
}))
a_stat$q <- p.adjust(a_stat$p, method = "BH")
a_stat$mark <- ifelse(a_stat$mech == "HR intact", "ref", stars(a_stat$q))
a_stat$xi <- match(a_stat$mech, levels(lab$mech))
pa <- raincloud(lab, "mech", levels(lab$mech)[-1], "Dominant lesion mechanism", "A") +
  geom_text(data = a_stat, aes(xi, 1.11, label = paste0("n = ", n)), size = SMALL, colour = "grey35") +
  geom_text(data = a_stat, aes(xi, 1.035, label = mark), size = TS) +
  theme(axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1))
# B: ordered event-class count, capped at three or more.
jt <- jt_test(lab$mirror_p, lab$n_ev_f)
b_stat <- bind_rows0(lapply(levels(lab$n_ev_f), function(g) {
  x <- lab$mirror_p[lab$n_ev_f == g]
  data.frame(event_count = g, n = length(x), median = safe_median(x),
             xi = match(g, levels(lab$n_ev_f)))
}))
pb <- raincloud(lab, "n_ev_f", c("1", "2", "3+"), "HR pathway event classes per tumour", "B") +
  geom_line(data = b_stat, aes(xi, median), linewidth = .3, colour = "grey35") +
  geom_point(data = b_stat, aes(xi, median), size = .8, colour = "grey20") +
  geom_text(data = b_stat, aes(xi, 1.11, label = paste0("n = ", n)), size = SMALL, colour = "grey35") +
  annotate("text", x = .52, y = 1.035, hjust = 0,
           label = paste("Jonckheere-Terpstra", fmt_p(jt$p)), size = SMALL)

# C: paired Spearman correlation, two-sided P; percentile bootstrap 95% CI.
cor_stats <- function(dat) {
  bind_rows0(lapply(levels(droplevels(dat$split)), function(sp) {
    g <- dat[dat$split == sp & is.finite(dat$mirror_p) & is.finite(dat$scar), ]
    rho <- p <- NA_real_; ci <- c(NA_real_, NA_real_)
    if (nrow(g) >= 3L && length(unique(g$mirror_p)) > 1L && length(unique(g$scar)) > 1L) {
      z <- suppressWarnings(cor.test(g$mirror_p, g$scar, method = "spearman", exact = FALSE))
      rho <- unname(z$estimate); p <- z$p.value
      boot <- replicate(N_BOOT, {
        i <- sample.int(nrow(g), nrow(g), replace = TRUE)
        suppressWarnings(cor(g$mirror_p[i], g$scar[i], method = "spearman"))
      })
      ci <- quant_ci(boot)
    }
    data.frame(split = sp, n = nrow(g), rho = rho, lo = ci[1], hi = ci[2], p = p)
  }))
}
plot_c <- function(dat, stat, panel = "C") {
  dat <- dat[is.finite(dat$scar) & is.finite(dat$mirror_p), ]
  if (!nrow(dat)) return(ggplot() + theme_void() + labs(title = panel) +
                           annotate("text", x = 0, y = 0, label = "No complete genomic HRD scores"))
  ymax <- max(100, SCAR_CUT, max(dat$scar)) * 1.24
  ss <- stat
  ss$line1 <- sprintf("rs = %.2f (95%% CI %.2f-%.2f)", ss$rho, ss$lo, ss$hi)
  ss$line2 <- sprintf("%s; n = %d", fmt_p(ss$p), ss$n)
  p <- ggplot(dat, aes(mirror_p, scar)) +
    geom_hline(yintercept = SCAR_CUT, linetype = 2, linewidth = .25, colour = "grey55") +
    geom_point(aes(colour = ifelse(mirror_call == 1, "pos", "neg")),
               size = .55, alpha = .65, show.legend = FALSE) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                linewidth = .4, fill = "grey70", alpha = .25) +
    geom_text(data = ss, aes(x = .02, y = ymax * .98, label = line1),
              hjust = 0, vjust = 1, inherit.aes = FALSE, size = SMALL) +
    geom_text(data = ss, aes(x = .02, y = ymax * .885, label = line2),
              hjust = 0, vjust = 1, inherit.aes = FALSE, size = SMALL) +
    annotate("text", x = .985, y = SCAR_CUT + 2, label = paste0("Cut-off ", SCAR_CUT),
             hjust = 1, vjust = 0, size = 5.8 / (72.27 / 25.4), colour = "grey40") +
    scale_colour_manual(values = COL_HRD) +
    scale_x_continuous(breaks = c(0, .5, 1), expand = expansion(mult = .02)) +
    scale_y_continuous(breaks = seq(0, floor(ymax / 25) * 25, 25), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, ymax)) +
    facet_wrap(~ split, nrow = 1, drop = TRUE) +
    labs(x = "MIRROR HRD score", y = "Genomic HRD score\n(TAI + LST + LOH)", title = panel)
  if (is.finite(THR)) p <- p + geom_vline(xintercept = THR, linetype = 2,
                                          linewidth = .25, colour = "grey55")
  p
}
c_stat <- cor_stats(cor_set)
pc <- plot_c(cor_set, c_stat)

# D: four genomic-scar metrics by the PRESERVED MIRROR calls, per cohort.
# Each cell uses its own complete metric data (no unnecessary scar-score filter).
MET <- c(scar = "HRD score", TAI = "TAI", LST = "LST", LOH = "LOH")
make_d <- function(dat, panel = "D") {
  slev <- levels(droplevels(dat$split))
  dl <- bind_rows0(lapply(names(MET), function(v) {
    x <- dat[is.finite(dat[[v]]), ]
    data.frame(split = as.character(x$split), metric = rep(unname(MET[v]), nrow(x)),
               value = x[[v]], grp = ifelse(x$mirror_call == 1, "HRD+", "HRD-"))
  }))
  dl$split <- factor(dl$split, levels = slev); dl$metric <- factor(dl$metric, levels = unname(MET))
  dl$grp <- factor(dl$grp, levels = c("HRD-", "HRD+"))
  dl$cls <- ifelse(dl$grp == "HRD+", "pos", "neg")
  dl$xi <- as.integer(dl$grp); dl$xj <- dl$xi - .23 + runif(nrow(dl), -.065, .065)
  ds <- bind_rows0(lapply(unname(MET), function(m) bind_rows0(lapply(slev, function(sp) {
    ss <- dl[dl$metric == m & dl$split == sp, ]
    a <- ss$value[ss$grp == "HRD+"]; b <- ss$value[ss$grp == "HRD-"]
    data.frame(metric = m, split = sp, n_pos = length(a), n_neg = length(b),
               median_pos = safe_median(a), median_neg = safe_median(b), p = safe_wilcox(a, b))
  }))))
  ds$q <- p.adjust(ds$p, "BH") # audited; primary cell labels retain original raw P convention.
  ds$label <- fmt_p(ds$p)
  ds$metric <- factor(ds$metric, levels = unname(MET)); ds$split <- factor(ds$split, levels = slev)
  mx <- vapply(unname(MET), function(m) max(c(1, dl$value[dl$metric == m]), na.rm = TRUE), numeric(1))
  ds$ytop <- mx[as.character(ds$metric)] * 1.28
  ds$ytxt <- mx[as.character(ds$metric)] * 1.13
  hv <- bind_rows0(lapply(unname(MET), function(m) bind_rows0(lapply(slev, function(sp) {
    ss <- dl[dl$metric == m & dl$split == sp, ]
    h <- half_violin(ss$value, as.character(ss$grp), c("HRD-", "HRD+"), width = .26)
    if (!nrow(h)) return(NULL)
    h$metric <- m; h$split <- sp; h
  }))))
  if (!nrow(hv)) hv <- data.frame(g = character(), x = numeric(), y = numeric(),
                                  metric = character(), split = character())
  hv$metric <- factor(hv$metric, levels = unname(MET)); hv$split <- factor(hv$split, levels = slev)
  hv$cls <- ifelse(hv$g == "HRD+", "pos", "neg")
  p <- ggplot() +
    geom_point(data = dl, aes(xj, value, colour = cls), size = .23, alpha = .30) +
    geom_polygon(data = hv, aes(x, y, group = g, colour = cls, fill = cls),
                 alpha = .3, linewidth = .22) +
    geom_boxplot(data = dl, aes(xi, value, group = xi, colour = cls),
                 width = .16, outlier.shape = NA, linewidth = .28, fill = "white") +
    geom_blank(data = ds, aes(x = 1.5, y = ytop)) +
    geom_text(data = ds, aes(x = 1.5, y = ytxt, label = label), size = SMALL) +
    scale_colour_manual(values = COL_HRD) + scale_fill_manual(values = COL_HRD) +
    scale_x_continuous(breaks = 1:2, labels = c("HRD-", "HRD+"),
                       limits = c(.48, 2.48), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(.025, .025)),
                       breaks = function(x) pretty(x, n = 3)) +
    facet_grid(metric ~ split, scales = "free_y", drop = FALSE) +
    labs(x = "MIRROR classification", y = NULL, title = panel) +
    theme(legend.position = "none", panel.spacing.x = grid::unit(1.5, "mm"),
          panel.spacing.y = grid::unit(1.5, "mm"),
          strip.text.x = element_text(size = 7.3), strip.text.y = element_text(size = 7.3),
          axis.text.x = element_text(size = 6.8))
  list(plot = p, stats = ds, data = dl)
}
d_result <- make_d(main)
pd <- d_result$plot

# E: both ROC curves, bands and paired DeLong tests use exactly the same patients.
GRID <- seq(0, 1, length.out = 101)
e_auc <- e_diff <- e_curve <- e_band <- list()
for (sp in SPLITS) {
  g <- paired[paired$split == sp, ]; y <- g$label
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  dl <- delong_test(g$mirror_p, g$scar, y)
  e_diff[[sp]] <- cbind(data.frame(split = sp, n = nrow(g)), as.data.frame(dl))
  if (!n1 || !n0) {
    warning("ROC skipped for ", sp, ": both reference classes are required.")
    next
  }
  scores <- list("MIRROR" = g$mirror_p, "Genomic HRD score" = g$scar)
  ba <- matrix(NA_real_, N_BOOT, 2)
  bt <- array(NA_real_, c(N_BOOT, length(GRID), 2))
  ip <- which(y == 1); ineg <- which(y == 0)
  # Paired stratified bootstrap: reuse each patient resample for both assays.
  for (b in seq_len(N_BOOT)) {
    ii <- c(ip[sample.int(length(ip), length(ip), replace = TRUE)], ineg[sample.int(length(ineg), length(ineg), replace = TRUE)])
    for (j in seq_along(scores)) {
      ba[b, j] <- auc_fun(scores[[j]][ii], y[ii])
      bt[b, , j] <- tpr_at(roc_fun(scores[[j]][ii], y[ii]), GRID)
    }
  }
  for (j in seq_along(scores)) {
    nm <- names(scores)[j]; ci <- quant_ci(ba[, j]); a <- auc_fun(scores[[j]], y)
    e_auc[[length(e_auc) + 1L]] <- data.frame(split = sp, method = nm, n = nrow(g),
                                              n_pos = n1, n_neg = n0, auc = a, lo = ci[1], hi = ci[2],
                                              label = sprintf("%s: %.3f (%.3f-%.3f)", if (j == 1) "MIRROR" else "Genomic HRD", a, ci[1], ci[2]),
                                              text_y = if (j == 1) .22 else .12)
    e_curve[[length(e_curve) + 1L]] <- cbind(roc_fun(scores[[j]], y), method = nm, split = sp)
    e_band[[length(e_band) + 1L]] <- data.frame(fpr = GRID, split = sp, method = nm,
                                                lo = apply(bt[, , j], 2, function(v) quant_ci(v)[1]),
                                                hi = apply(bt[, , j], 2, function(v) quant_ci(v)[2]))
  }
  dif_ci <- quant_ci(ba[, 1] - ba[, 2])
  e_diff[[sp]]$bootstrap_difference_lo <- dif_ci[1]
  e_diff[[sp]]$bootstrap_difference_hi <- dif_ci[2]
}
e_auc <- bind_rows0(e_auc); e_diff <- bind_rows0(e_diff)
e_curve <- bind_rows0(e_curve); e_band <- bind_rows0(e_band)
e_diff$label <- paste("Paired DeLong", fmt_p(e_diff$p))
e_diff$split <- factor(e_diff$split, levels = SPLITS)
if (nrow(e_curve)) {
  e_curve$split <- factor(e_curve$split, levels = SPLITS)
  e_band$split <- factor(e_band$split, levels = SPLITS)
  e_auc$split <- factor(e_auc$split, levels = SPLITS)
  pe <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = .25) +
    geom_ribbon(data = e_band, aes(fpr, ymin = lo, ymax = hi, fill = method), alpha = .14) +
    geom_path(data = e_curve, aes(fpr, tpr, colour = method), linewidth = .45) +
    annotate("text", x = .03, y = .31, label = "AUC (95% CI)", hjust = 0, size = 5.8 / (72.27 / 25.4)) +
    geom_text(data = e_auc, aes(.03, text_y, label = label, colour = method),
              hjust = 0, size = 5.8 / (72.27 / 25.4)) +
    geom_text(data = e_diff, aes(.03, .035, label = label), hjust = 0,
              size = 5.8 / (72.27 / 25.4)) +
    scale_colour_manual(values = COL_ASSAY_LINE) + scale_fill_manual(values = COL_ASSAY_LINE) +
    scale_x_continuous(breaks = c(0, .5, 1), expand = expansion(mult = .015)) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = expansion(mult = .015)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    facet_wrap(~ split, nrow = 1, drop = FALSE) +
    labs(x = "1 - specificity", y = "Sensitivity", title = "E") +
    theme(legend.position = "none", aspect.ratio = 1, panel.spacing.x = grid::unit(2, "mm"))
} else {
  pe <- ggplot() + theme_void() + labs(title = "E") +
    annotate("text", x = 0, y = 0, label = "Insufficient paired data for ROC analysis")
}

# F: detection rates in alteration groups; HR-intact row shows false-positive rate.
f_rows <- c(names(EV), "HR intact (false positives)")
subset_event <- function(dat, nm) {
  if (nm == "HR intact (false positives)") dat[dat$label == 0, ] else dat[grepl(EV[[nm]], dat$gene), ]
}
detection_rows <- function(use_paired = TRUE) {
  out <- lapply(f_rows, function(nm) {
    g <- subset_event(lab, nm)
    gp <- g[is.finite(g$scar) & binary_ok(g$mirror_call), ]
    gm <- if (use_paired) gp else g
    rbind(data.frame(event = nm, method = "MIRROR", k = sum(gm$mirror_call == 1), n = nrow(gm)),
          data.frame(event = nm, method = "Genomic HRD score", k = sum(gp$scar_call == 1), n = nrow(gp)))
  })
  z <- bind_rows0(out)
  z$rate <- ifelse(z$n > 0, z$k / z$n, NA_real_)
  z$lo <- mapply(function(k, n) wilson(k, n)[1], z$k, z$n)
  z$hi <- mapply(function(k, n) wilson(k, n)[2], z$k, z$n)
  z
}
f_full <- detection_rows(FALSE)
f_dat <- detection_rows(TRUE)
f_sig <- bind_rows0(lapply(f_rows, function(nm) {
  g <- subset_event(lab, nm); g <- g[is.finite(g$scar), ]
  mc <- mcnemar_exact(g$mirror_call, g$scar_call)
  data.frame(event = nm, n_pair = mc$n, mirror_only = mc$mirror_only,
             genomic_only = mc$scar_only, p = mc$p)
}))
f_sig$q <- p.adjust(f_sig$p, method = "BH"); f_sig$mark <- stars(f_sig$q)
for (nm in f_rows) stopifnot(length(unique(f_dat$n[f_dat$event == nm])) == 1L,
                             unique(f_dat$n[f_dat$event == nm]) == f_sig$n_pair[f_sig$event == nm])
f_dat$text <- ifelse(f_dat$n > 0, sprintf("%.1f%% (%d/%d)", 100 * f_dat$rate, f_dat$k, f_dat$n), "NA (0/0)")
ord <- f_dat %>% filter(method == "MIRROR", event != "HR intact (false positives)") %>%
  arrange(rate, event) %>% pull(event)
f_levels <- c("HR intact (false positives)", ord)
f_dat$event <- factor(f_dat$event, levels = f_levels); f_sig$event <- factor(f_sig$event, levels = f_levels)
f_dat$method <- factor(f_dat$method, levels = ASSAYS)
DG <- position_dodge(width = .90)
pf <- ggplot(f_dat, aes(x = rate, y = event, fill = method)) +
  geom_col(position = DG, width = .80, colour = "grey35", linewidth = .18, na.rm = TRUE) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", position = DG,
                width = .24, linewidth = .28, colour = "black", na.rm = TRUE) +
  geom_text(aes(x = 1.035, label = text), position = DG, hjust = 0,
            size = 6.3 / (72.27 / 25.4), show.legend = FALSE) +
  geom_text(data = f_sig, aes(x = 1.63, y = event, label = mark),
            inherit.aes = FALSE, size = SMALL) +
  geom_vline(xintercept = 1, linetype = 3, linewidth = .22, colour = "grey65") +
  scale_fill_manual(values = COL_ASSAY, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1.70), breaks = c(0, .5, 1),
                     labels = c("0", "50", "100"), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) ifelse(x == "HR intact (false positives)",
                                               "HR intact\n(false positives)", x), expand = expansion(add = .65)) +
  labs(x = "Detection rate (%)", y = NULL, title = "F") +
  theme(legend.position = "bottom", legend.justification = "left",
        legend.text = element_text(size = 6.5), legend.key.width = grid::unit(2.5, "mm"),
        axis.ticks.y = element_blank(), axis.text.y = element_text(size = 6.5),
        plot.margin = margin(1.4, 2, 1.4, 1.4, "mm"))

# ----------------------------- Assembly/export ------------------------------
# One patchwork level avoids nested-layout flattening. Left upper block is 68%;
# ROC block is 48% at the bottom. D is 4 metrics x 2 cohorts, not three cohorts.
design <- c(area(1, 1, 1, 2), area(1, 3, 1, 4), area(2, 1, 2, 4),
            area(1, 5, 2, 5), area(3, 1, 3, 3), area(3, 4, 3, 5))
fig <- pa + pb + pc + pd + pe + pf +
  plot_layout(design = design, widths = c(17, 17, 14, 20, 32), heights = c(48, 46, 70))
FIG_W <- 190
FIG_H <- 164
save_figure <- function(p, stem, w_mm, h_mm, tiff = FALSE) {
  dest <- file.path(OUT_DIR, stem)
  grDevices::pdf(paste0(dest, ".pdf"), width = w_mm / 25.4, height = h_mm / 25.4,
                 family = "sans", useDingbats = FALSE, onefile = FALSE)
  tryCatch(print(p), finally = grDevices::dev.off())
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(dest, "_preview.png"), width = w_mm, height = h_mm,
                  units = "mm", res = 200, background = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    if (tiff) {
      ragg::agg_tiff(paste0(dest, ".tiff"), width = w_mm, height = h_mm,
                     units = "mm", res = 600, compression = "lzw", background = "white")
      tryCatch(print(p), finally = grDevices::dev.off())
    }
  } else {
    ggsave(paste0(dest, "_preview.png"), p, width = w_mm, height = h_mm,
           units = "mm", dpi = 200, bg = "white")
    if (tiff) ggsave(paste0(dest, ".tiff"), p, width = w_mm, height = h_mm,
                     units = "mm", dpi = 600, compression = "lzw", bg = "white")
  }
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(dest, ".svg"), width = w_mm / 25.4, height = h_mm / 25.4,
                     bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
}
save_figure(fig, "Figure3_optimized", FIG_W, FIG_H, tiff = TRUE)
save_figure(pd, "Figure3_panel_D", 68, 100)
save_figure(pf, "Figure3_panel_F", 120, 88)
internal_c_stat <- internal_d_stat <- data.frame()
if (nrow(internal)) {
  ic <- internal %>% filter(is.finite(scar))
  internal_c_stat <- cor_stats(ic)
  ip_c <- plot_c(ic, internal_c_stat, "A")
  id <- make_d(internal, "B"); internal_d_stat <- id$stats
  internal_fig <- ip_c + id$plot + plot_layout(widths = c(1.45, 1)) +
    plot_annotation(title = "TCGA internal cohort: additional genomic-scar analyses",
                    theme = theme(plot.title = element_text(size = 10, face = "plain")))
  save_figure(internal_fig, "Additional_TCGA_internal", 150, 100)
}
stats_tables <- list(
  audit_input_counts = audit_input, audit_thresholds = threshold_audit,
  audit_mechanism = mechanism_audit,
  stat_A_mechanism = a_stat, stat_B_event_count = b_stat,
  stat_B_Jonckheere_Terpstra = as.data.frame(jt),
  stat_C_correlation = c_stat, stat_D_genomic_metrics = d_result$stats,
  stat_E_AUC = e_auc, stat_E_paired_DeLong = e_diff,
  stat_F_detection_paired = f_dat, stat_F_paired_McNemar = f_sig,
  audit_F_original_method_specific_denominators = f_full,
  additional_internal_C_correlation = internal_c_stat,
  additional_internal_D_genomic_metrics = internal_d_stat)
for (nm in names(stats_tables)) write.csv(stats_tables[[nm]], file.path(OUT_DIR, paste0(nm, ".csv")),
                                          row.names = FALSE, na = "")
write.csv(paired %>% select(sample_id, split, label, mirror_p, mirror_call, scar, scar_call),
          file.path(OUT_DIR, "audit_paired_analysis_patients.csv"), row.names = FALSE)
notes <- c(
  "File: Figure 2_optimized.R; manuscript panel set: Figure 3 A-F.",
  "Input labels (HRD_label) and MIRROR classifications (Multi_PRE_HRD) are preserved.",
  "HRD+ = blue #5B87AE; HRD- = orange #E0995E. Assay colours remain MIRROR blue / genomic HRD orange.",
  "Main panels include only TCGA-train and TCGA-val. Internal cohort analyses are exported separately.",
  "A/B/F pool labelled training and validation patients; this is not an independent validation analysis.",
  "A: two-sided Wilcoxon rank-sum tests with continuity correction, BH correction across four comparisons.",
  "B: two-sided Jonckheere-Terpstra normal approximation with ties correction; event groups preserve source regexes.",
  "C: Spearman correlation, two-sided approximate P; paired percentile-bootstrap 95% CI.",
  "C regression line: ordinary least squares; grey ribbon: ggplot2 pointwise 95% confidence band for mean fit.",
  "D: two-sided Wilcoxon rank-sum tests with continuity correction; raw P values shown; BH values also in CSV.",
  "A/B/D: box = median and interquartile range; whiskers = most extreme observations within 1.5 x IQR.",
  "E: tie-aware ROC threshold blocks joined linearly; AUC is rank-based with half credit for ties.",
  "E: paired stratified bootstrap for AUC 95% CIs and pointwise ROC bands; paired two-sided DeLong test.",
  "F: both bars, Wilson 95% CIs, and exact two-sided McNemar tests use the same complete paired patients.",
  "F: BH correction across ten comparisons; *** q < 0.001, ** q < 0.01, * q < 0.05, ns q >= 0.05.",
  "F: molecular alteration groups can overlap; HR-intact row reports false-positive rate.",
  "Original method-specific denominators are retained only in a separate audit CSV.",
  paste0("Genomic HRD cut-off = ", SCAR_CUT, ". MIRROR plotting cut-off inferred from existing calls = ", signif(THR, 7), "."),
  paste0("Bootstrap replicates = ", N_BOOT, "; random seed = ", SEED, "."),
  paste0("Input path = ", normalizePath(CSV_PATH, winslash = "/", mustWork = TRUE)),
  paste0("Input MD5 = ", unname(tools::md5sum(CSV_PATH)))
)
writeLines(notes, file.path(OUT_DIR, "README_analysis_notes.txt"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
cat("Completed. Output directory:", normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE), "\n")
print(c_stat)
print(e_diff[, intersect(c("split", "n", "auc_mirror", "auc_genomic", "difference", "p"), names(e_diff))])
print(f_sig[, c("event", "n_pair", "mirror_only", "genomic_only", "p", "q")])
