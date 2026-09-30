# =============================================================================
# Figure 1 | Multimodal HRD model validation
# =============================================================================
# Scientific scope: Discrimination, calibration, clinical utility and subgroup AUC.
# Usage: Rscript 01_Model_Discrimination_and_Calibration.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: Clinical CSV with HRD_prob, Multi_PRE_HRD, HRD_label, has_HRD_label, dataset.
# Outputs: Panel PDFs/PNGs and source_data/ beneath MIRROR_OUTPUT_DIR.
# Dependencies: ggplot2, patchwork, scales.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# =============================================================================


options(stringsAsFactors = FALSE, warn = 1)

.need <- c("ggplot2", "patchwork", "scales")
.miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss)) stop("Install required packages: ", paste(.miss, collapse = ", "))
invisible(lapply(.need, library, character.only = TRUE))

# ======================== Paths ==============================================
DATA_FILE <- Sys.getenv("MIRROR_CLINICAL_CSV", unset = "<PATH_TO_CLINICAL_CSV>")
OUT_DIR   <- Sys.getenv("MIRROR_OUTPUT_DIR", unset = "<PATH_TO_OUTPUT_DIRECTORY>")
SRC_DIR   <- file.path(OUT_DIR, "source_data")
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
prepare_output_dir(SRC_DIR, "SRC_DIR")

# ======================== Global switches ====================================
SCORE_COL       <- "HRD_prob"          # continuous multimodal score
PRED_COL        <- "Multi_PRE_HRD"     # binary call shipped with the CSV
SCORE_WSI       <- "HRD_prob_WSI"
SCORE_RNA       <- "HRD_prob_RNASeq"
USE_GIVEN_CALL  <- TRUE                # FALSE -> re-derive a Youden threshold on Train
N_BOOT          <- 2000                # bootstrap replicates (CIs and ROC bands)
AGE_CUT         <- 58                  # age dichotomy used in panel H
N_MIN_SUB       <- 5                   # a subgroup needs this many cases to be plotted
N_MIN_FLAG      <- 25                  # below this the AUC is flagged with '*'
HIST_BIN_W      <- 0.05                # histogram bin width  (panel E, lower)
HL_GROUPS       <- 10                  # Hosmer-Lemeshow groups (df = HL_GROUPS - 2)
DCA_MAX         <- 0.80                # decision-curve threshold range
FIG_W           <- 190                 # composite width  (mm)
FIG_H           <- 230                 # composite height (mm)

## ---- panel B: which version is embedded in the composite Figure 1 ----------
## "pooled" = Train + Val paired subset (the original behaviour)
## "train"  = TCGA-train paired subset only
## "val"    = TCGA-val   paired subset only
## The other two are still computed and exported as standalone files.
FIG1_PANEL_B <- "pooled"

## ---- risk scale used by panels E and G (see the long comment below) --------
RECALIBRATE        <- TRUE       # FALSE -> both panels use the raw HRD_prob
RECAL_METHOD       <- "platt"    # "platt" (logistic) or "isotonic"
RECAL_FIT_ON_TRAIN <- TRUE       # TRUE -> fit on Train, apply frozen weights to Val

## ---- calibration curve presentation (panel E, upper) ----------------------
CAL_BINNING <- "quantile"        # "quantile" (equal-frequency) or "equal" (equal-width)
CAL_N_BIN   <- 10
CAL_SMOOTH  <- TRUE              # overlay a loess calibration curve
if (!is.numeric(N_BOOT) || length(N_BOOT) != 1L || !is.finite(N_BOOT) ||
    N_BOOT < 100L || N_BOOT %% 1 != 0)
  stop("N_BOOT must be an integer of at least 100.")
if (!FIG1_PANEL_B %in% c("pooled", "train", "val"))
  stop("FIG1_PANEL_B must be pooled, train or val.")
if (!RECAL_METHOD %in% c("platt", "isotonic"))
  stop("RECAL_METHOD must be platt or isotonic.")
if (!CAL_BINNING %in% c("quantile", "equal") || length(CAL_N_BIN) != 1L ||
    !is.finite(CAL_N_BIN) || CAL_N_BIN < 3L)
  stop("Specify a supported CAL_BINNING and at least three calibration bins.")

set.seed(20260805)

# ======================== Palette / theme ====================================
COL_TRAIN <- "#416993"   # TCGA-train, and HRD+ in panels C / E / F
COL_VAL   <- "#AE6534"   # TCGA-val,   and HRD- in panels C / E / F
COL_WSI   <- "#3F7F7A"   # WSI-only    (panel B)
COL_INK   <- "#1A1A1A"; COL_GREY <- "#6B7280"; COL_RULE <- "#B8BEC6"
STRIPE    <- "#F3F5F7"

lighten <- function(col, f = 0.65) {
  m <- grDevices::col2rgb(col)
  grDevices::rgb(t(255 - (255 - m) * (1 - f)), maxColorValue = 255)
}
FILL_TRAIN <- lighten(COL_TRAIN, 0.55)
FILL_VAL   <- lighten(COL_VAL,   0.55)

FONT      <- "sans"          # set to "Arial" if the font is installed
BASE_SIZE <- 8
TXT_MM    <- BASE_SIZE / 2.845276     # geom_text size for BASE_SIZE pt
SPLITS    <- c("Train", "Val")
LBL       <- c(Train = "TCGA-train", Val = "TCGA-val")
COL_SET   <- c(Train = COL_TRAIN, Val = COL_VAL)
COL_TRUTH <- c(`HRD-` = COL_VAL, `HRD+` = COL_TRAIN)

theme_fig <- function(base_size = BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size, base_family = FONT) +
    ggplot2::theme(
      text              = element_text(colour = COL_INK),
      panel.grid        = element_blank(),
      axis.line         = element_line(linewidth = 0.35, colour = COL_INK, lineend = "square"),
      axis.ticks        = element_line(linewidth = 0.35, colour = COL_INK),
      axis.ticks.length = unit(0.9, "mm"),
      axis.text         = element_text(size = base_size, colour = COL_INK),
      axis.title.x      = element_text(size = base_size, margin = margin(t = 1.0, unit = "mm")),
      axis.title.y      = element_text(size = base_size, margin = margin(r = 1.0, unit = "mm")),
      plot.title        = element_text(size = base_size, face = "bold", colour = COL_INK,
                                       hjust = 0, margin = margin(b = 0.8, unit = "mm")),
      plot.subtitle     = element_text(size = base_size, colour = COL_INK, hjust = 0.5,
                                       margin = margin(b = 1.0, unit = "mm")),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, colour = COL_INK),
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size, colour = COL_INK),
      legend.key.size   = unit(2.8, "mm"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      plot.margin       = margin(1.2, 1.6, 1.2, 1.2, "mm")
    )
}

head_panel <- function(p, letter = NULL, sub = NULL) {
  p + labs(title = if (is.null(letter)) " " else letter, subtitle = sub)
}

## legend placed inside the panel; ggplot2 >= 3.5 renamed the argument
legend_inside <- function(x, y, just = c(0, 1)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(x, y),
                   legend.justification = just)
  else
    ggplot2::theme(legend.position = c(x, y), legend.justification = just)
}

save_plot <- function(p, file, width = 90, height = 70, dir = OUT_DIR) {
  ggplot2::ggsave(file.path(dir, paste0(file, ".pdf")), plot = p,
                  width = width, height = height, units = "mm", device = "pdf")
  ggplot2::ggsave(file.path(dir, paste0(file, ".png")), plot = p,
                  width = width, height = height, units = "mm", dpi = 600)
  invisible(TRUE)
}

# =============================================================================
# Statistical primitives
# =============================================================================
auc_fast <- function(y, p) {
  ok <- is.finite(p) & !is.na(y); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 < 1 || n0 < 1) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_curve <- function(y, p) {
  ok <- is.finite(p) & !is.na(y); y <- y[ok]; p <- p[ok]
  o <- order(p, decreasing = TRUE); y <- y[o]; p <- p[o]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  keep <- c(which(diff(p) != 0), length(p))
  data.frame(fpr = c(0, fp[keep] / n0), tpr = c(0, tp[keep] / n1))
}
roc_interp <- function(rc, grid) {
  stats::approx(rc$fpr, rc$tpr, xout = grid, method = "constant",
                f = 1, ties = "max", rule = 2)$y
}

boot_idx <- function(y) {
  i1 <- which(y == 1); i0 <- which(y == 0)
  c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE))
}

auc_ci <- function(y, p, B = N_BOOT) {
  b <- vapply(seq_len(B), function(i) { ib <- boot_idx(y); auc_fast(y[ib], p[ib]) }, 0)
  c(auc = auc_fast(y, p), unname(stats::quantile(b, c(0.025, 0.975), na.rm = TRUE)))
}

roc_band <- function(y, p, grid = seq(0, 1, by = 0.005), B = N_BOOT) {
  M <- matrix(NA_real_, nrow = B, ncol = length(grid))
  for (b in seq_len(B)) {
    ib <- boot_idx(y)
    M[b, ] <- roc_interp(roc_curve(y[ib], p[ib]), grid)
  }
  data.frame(fpr = grid,
             lo = apply(M, 2, stats::quantile, 0.025, na.rm = TRUE),
             hi = apply(M, 2, stats::quantile, 0.975, na.rm = TRUE))
}

## Wilson score interval -- stays finite at 0% and 100%, unlike the Wald CI
wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(est = NA_real_, lo = NA_real_, hi = NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  ph <- x / n; d <- 1 + z^2 / n
  ctr <- (ph + z^2 / (2 * n)) / d
  hw  <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / d
  c(est = ph, lo = max(0, ctr - hw), hi = min(1, ctr + hw))
}

## ---- DeLong covariance for k classifiers on the SAME samples ---------------
delong <- function(y, plist) {
  if (!is.list(plist) || length(plist) < 2L)
    stop("plist must be a list of at least two score vectors")
  k  <- length(plist)
  nm <- names(plist); if (is.null(nm)) nm <- paste0("model", seq_len(k))
  
  fail <- list(auc = stats::setNames(rep(NA_real_, k), nm),
               S   = matrix(NA_real_, k, k), m = 0L, n = 0L, ok = FALSE)
  
  y <- suppressWarnings(as.integer(y))
  if (length(unique(vapply(plist, length, 0L))) != 1L ||
      length(y) != length(plist[[1]])) {
    warning("delong(): y and the score vectors have different lengths.")
    return(fail)
  }
  ok <- !is.na(y) & Reduce(`&`, lapply(plist, is.finite))
  y  <- y[ok]; plist <- lapply(plist, function(p) p[ok])
  
  m <- sum(y == 1L); n <- sum(y == 0L)
  if (m < 2L || n < 2L) {
    warning(sprintf("delong(): paired subset too small (HRD+ %d, HRD- %d).", m, n))
    fail$m <- m; fail$n <- n
    return(fail)
  }
  X <- lapply(plist, function(p) p[y == 1L])
  Y <- lapply(plist, function(p) p[y == 0L])
  V10 <- matrix(NA_real_, m, k); V01 <- matrix(NA_real_, n, k)
  for (r in seq_len(k)) {
    psi <- outer(X[[r]], Y[[r]], function(a, b) (a > b) + 0.5 * (a == b))
    V10[, r] <- rowMeans(psi); V01[, r] <- colMeans(psi)
  }
  auc <- stats::setNames(colMeans(V10), nm)
  S <- tryCatch(stats::cov(V10) / m + stats::cov(V01) / n,
                error = function(e) matrix(NA_real_, k, k))
  if (is.null(dim(S)) || any(dim(S) != k)) S <- matrix(NA_real_, k, k)
  list(auc = auc, S = S, m = m, n = n, ok = TRUE)
}

delong_test <- function(dl, i, j) {
  out <- c(z = NA_real_, p = NA_real_)
  if (!is.list(dl) || is.null(dl$auc) || is.null(dl$S)) return(out)
  k <- length(dl$auc)
  if (!isTRUE(dl$ok) || k < 2L || max(i, j) > k) return(out)
  if (is.null(dim(dl$S)) || any(dim(dl$S) != k)) return(out)
  L <- rep(0, k); L[i] <- 1; L[j] <- -1
  v <- tryCatch(as.numeric(t(L) %*% dl$S %*% L), error = function(e) NA_real_)
  if (length(v) != 1L || !is.finite(v) || v <= 0) return(out)
  z <- as.numeric(dl$auc[i] - dl$auc[j]) / sqrt(v)
  if (length(z) != 1L || !is.finite(z)) return(out)
  c(z = z, p = 2 * stats::pnorm(-abs(z)))
}

boot_auc_diff_p <- function(y, p1, p2, B = N_BOOT) {
  ok <- !is.na(y) & is.finite(p1) & is.finite(p2)
  y <- y[ok]; p1 <- p1[ok]; p2 <- p2[ok]
  if (sum(y == 1) < 2 || sum(y == 0) < 2) return(NA_real_)
  d <- vapply(seq_len(B), function(b) {
    ib <- boot_idx(y); auc_fast(y[ib], p1[ib]) - auc_fast(y[ib], p2[ib])
  }, 0)
  d <- d[is.finite(d)]
  if (!length(d)) return(NA_real_)
  min(1, 2 * min(mean(d <= 0), mean(d >= 0)))
}

## Hosmer-Lemeshow goodness of fit (groups by quantiles of predicted risk)
hl_test <- function(y, p, g = HL_GROUPS) {
  br <- unique(stats::quantile(p, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE))
  if (length(br) < 3) return(list(chisq = NA_real_, df = NA_integer_, p = NA_real_))
  grp <- cut(p, breaks = br, include.lowest = TRUE)
  o1 <- tapply(y, grp, sum); e1 <- tapply(p, grp, sum); nn <- tapply(y, grp, length)
  keep <- is.finite(e1) & e1 > 0 & (nn - e1) > 0
  chi <- sum((o1[keep] - e1[keep])^2 / (e1[keep] * (1 - e1[keep] / nn[keep])))
  df  <- sum(keep) - 2
  list(chisq = chi, df = df, p = stats::pchisq(chi, df, lower.tail = FALSE))
}

## Calibration-in-the-large (intercept) and calibration slope.
## Perfect calibration is intercept 0 / slope 1. A slope below 1 means the
## predictions are too extreme, above 1 that they are too compressed -- which
## is the direction this model sits in before recalibration.
cal_int_slope <- function(y, p) {
  lp <- stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  s <- tryCatch(unname(stats::coef(stats::glm(y ~ lp, family = stats::binomial()))[2]),
                error = function(e) NA_real_)
  i <- tryCatch(unname(stats::coef(stats::glm(y ~ 1 + offset(lp),
                                              family = stats::binomial()))[1]),
                error = function(e) NA_real_)
  c(intercept = i, slope = s)
}

metric_set <- function(y, yh) {
  tp <- sum(y == 1 & yh == 1); fn <- sum(y == 1 & yh == 0)
  fp <- sum(y == 0 & yh == 1); tn <- sum(y == 0 & yh == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if (is.finite(prec) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_
  c(Sensitivity = sens, Specificity = spec, Precision = prec,
    `F1-score` = f1, Accuracy = (tp + tn) / length(y))
}

fmt_p_exact <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x <= 0)   return("< 1e-300")
    if (x < 1e-4) return(sprintf("%.1e", x))
    sprintf("%.3f", x)
  }, character(1), USE.NAMES = FALSE)
}
fmt_p_fig <- function(p) ifelse(is.na(p), "NA",
                                ifelse(p < 1e-4, "P < 0.0001", paste0("P = ", sprintf("%.3f", p))))
fmt_ci <- function(e, lo, hi, d = 3)
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), e, lo, hi)

# =============================================================================
# Drawing primitives -- hand-built boxplot, half violin, dodged histogram
# =============================================================================
bstat <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  q <- unname(stats::quantile(v, c(0.25, 0.5, 0.75))); iqr <- q[3] - q[1]
  data.frame(ymin = min(v[v >= q[1] - 1.5 * iqr]), lower = q[1], middle = q[2],
             upper = q[3], ymax = max(v[v <= q[3] + 1.5 * iqr]), n = length(v))
}
box_layers <- function(bs, w = 0.11, cap = 0.05, lw = 0.4) {
  list(
    geom_segment(data = bs, aes(x = x, xend = x, y = ymin, yend = ymax, colour = col),
                 linewidth = lw * 0.85, lineend = "butt"),
    geom_segment(data = bs, aes(x = x - cap, xend = x + cap, y = ymin, yend = ymin, colour = col),
                 linewidth = lw * 0.85),
    geom_segment(data = bs, aes(x = x - cap, xend = x + cap, y = ymax, yend = ymax, colour = col),
                 linewidth = lw * 0.85),
    geom_rect(data = bs, aes(xmin = x - w, xmax = x + w, ymin = lower, ymax = upper,
                             colour = col), fill = "white", linewidth = lw),
    geom_segment(data = bs, aes(x = x - w, xend = x + w, y = middle, yend = middle,
                                colour = col), linewidth = lw * 1.5)
  )
}
DENS_W <- 0.26; DENS_OFF <- 0.07; PT_OFF <- -0.24; PT_JIT <- 0.075
half_violin <- function(v, xbase, adjust = 1.1, n = 512) {
  v <- v[is.finite(v)]
  if (length(v) < 5) return(NULL)
  d <- stats::density(v, adjust = adjust, n = n, from = min(v), to = max(v))
  w <- DENS_W * d$y / max(d$y)
  data.frame(x = c(xbase + DENS_OFF + w, rep(xbase + DENS_OFF, length(w))),
             y = c(d$x, rev(d$x)))
}
jit <- function(n, w) stats::runif(n, -w, w)

## ---- dodged histogram with hand-placed bars --------------------------------
## position_dodge() resolves bar placement at draw time from the width
## aesthetic, so the bars can drift depending on which bins are populated.
## Here each bin is split explicitly: HRD- takes the left half, HRD+ the right
## half, with a fixed gap. The geometry is fully determined by the data.
hist_bars <- function(v, g, binw = HIST_BIN_W, gap = 0.08) {
  br  <- seq(0, 1, by = binw)
  vv  <- pmin(pmax(v, 0), 1)
  bn  <- cut(vv, breaks = br, include.lowest = TRUE)
  lev <- levels(g)
  tb  <- as.data.frame(table(bin = bn, truth = factor(g, levels = lev)),
                       stringsAsFactors = FALSE)
  tb  <- tb[tb$Freq > 0, , drop = FALSE]
  if (!nrow(tb)) return(tb)
  left <- br[match(tb$bin, levels(bn))]
  k    <- match(tb$truth, lev)              # 1 = HRD-, 2 = HRD+
  half <- binw / 2
  tb$xmin  <- left + (k - 1) * half
  tb$xmax  <- tb$xmin + half * (1 - gap)
  tb$truth <- factor(tb$truth, levels = lev)
  tb
}

# =============================================================================
# Data ingestion
# -----------------------------------------------------------------------------
# The delimiter is detected from the header line: the export is sometimes
# comma- and sometimes tab-separated despite the .csv extension.
# =============================================================================
read_scores <- function(path = DATA_FILE) {
  if (!file.exists(path)) stop("Data file not found: ", path)
  hdr <- readLines(path, n = 1L, warn = FALSE)
  n_tab <- length(gregexpr("\t", hdr, fixed = TRUE)[[1]])
  n_com <- length(gregexpr(",",  hdr, fixed = TRUE)[[1]])
  if (!grepl("\t", hdr, fixed = TRUE)) n_tab <- 0
  if (!grepl(",",  hdr, fixed = TRUE)) n_com <- 0
  sep <- if (n_tab > n_com) "\t" else ","
  rd <- function(enc) utils::read.table(path, sep = sep, header = TRUE, quote = "\"",
                                        comment.char = "", check.names = TRUE,
                                        fileEncoding = enc, na.strings = c("NA", ""))
  df <- try(rd("UTF-8-BOM"), silent = TRUE)
  if (inherits(df, "try-error") || ncol(df) < 5) df <- try(rd("UTF-8"), silent = TRUE)
  if (inherits(df, "try-error") || ncol(df) < 5) df <- rd("GB18030")
  names(df) <- make.names(names(df))
  message(sprintf("Read %d rows x %d columns (separator: %s).",
                  nrow(df), ncol(df), if (sep == "\t") "TAB" else "comma"))
  df
}

num <- function(x) suppressWarnings(as.numeric(x))
dat <- read_scores()

for (cc in c(SCORE_COL, PRED_COL, "HRD_label", "has_HRD_label", "dataset"))
  if (!cc %in% names(dat)) stop("Required column missing: ", cc)

dat$score     <- num(dat[[SCORE_COL]])
dat$pred      <- as.integer(round(num(dat[[PRED_COL]])))
dat$HRD_label <- as.integer(num(dat$HRD_label))
dat$dataset   <- factor(trimws(dat$dataset), levels = c("Train", "Val", "Internal",
                                                        "CHCAMS", "HMUCH", "PLCO"))

lab <- dat[!is.na(dat$has_HRD_label) & num(dat$has_HRD_label) == 1 &
             !is.na(dat$HRD_label) & is.finite(dat$score) & !is.na(dat$pred) &
             dat$dataset %in% SPLITS, ]
lab$dataset <- factor(as.character(lab$dataset), levels = SPLITS)
if (!nrow(lab)) stop("No labelled cases with valid scores and calls.")
for (cohort in SPLITS) {
  y <- lab$HRD_label[lab$dataset == cohort]
  if (length(unique(y)) < 2L) stop("Both HRD classes are required in ", cohort, ".")
}
lab$hr      <- factor(ifelse(lab$HRD_label == 1, "HRD+", "HRD-"),
                      levels = c("HRD-", "HRD+"))

if (!USE_GIVEN_CALL) {
  tr <- lab[lab$dataset == "Train", ]
  cand <- sort(unique(tr$score))
  j <- vapply(cand, function(t) {
    yh <- as.integer(tr$score >= t)
    mm <- metric_set(tr$HRD_label, yh); unname(mm["Sensitivity"] + mm["Specificity"] - 1)
  }, 0)
  THR <- cand[which.max(j)]
  lab$pred <- as.integer(lab$score >= THR)
  message(sprintf("Youden threshold re-derived on Train: %.4f", THR))
}

message(sprintf("Labelled cases: Train n = %d (HRD+ %d) | Val n = %d (HRD+ %d)",
                sum(lab$dataset == "Train"), sum(lab$dataset == "Train" & lab$HRD_label == 1),
                sum(lab$dataset == "Val"),   sum(lab$dataset == "Val"   & lab$HRD_label == 1)))

split_df <- function(ds) lab[lab$dataset == ds, ]

# =============================================================================
# Recalibration -- the risk scale shared by panels E and G
# -----------------------------------------------------------------------------
# The raw model output is a discrimination score, not a probability: before
# recalibration the predictions sit almost entirely between 0.2 and 0.8, which
# is what drives the Brier scores and the Hosmer-Lemeshow P values in the
# uncorrected panel E. Rescaling them onto a probability axis is a standard
# and publishable step -- but only under two conditions.
#
#   1. FIT ON THE TRAINING SPLIT, APPLY FROZEN TO VALIDATION.
#      RECAL_FIT_ON_TRAIN = TRUE does this. Fitting the recalibration inside
#      the validation split and then reporting the resulting curve as
#      validation calibration is circular: the same 106 patients would be used
#      to estimate the correction and to certify it. That number could not be
#      defended if a reviewer asked for the source data, and this script would
#      be the evidence trail.
#
#   2. STATE IT IN THE METHODS AND IN THE FIGURE LEGEND.
#      Something like: "Predicted probabilities were recalibrated by logistic
#      (Platt) rescaling fitted in the training split; the same rescaled risks
#      are shown in panel e and thresholded in panel g." A recalibration that
#      is not disclosed is not a presentational choice.
#
# Two further points that belong in the legend rather than being left implicit:
#   * The TRAINING calibration curve after an in-sample fit lies on the
#     diagonal by construction. It is not evidence of anything and must be
#     labelled as the fitted, not the validated, curve.
#   * Platt and isotonic rescaling are both monotone, so AUC, the ROC curves,
#     the confusion matrices and the subgroup forest are all unchanged. Only
#     panels E and G move. That is the tell that this is a calibration fix and
#     not a change in the model.
#
# The console prints Brier, calibration intercept, calibration slope and the
# Hosmer-Lemeshow test both before and after, so the improvement can be
# reported as a number rather than asserted from the look of the curve.
# =============================================================================
fit_recal <- function(d, method = RECAL_METHOD) {
  if (method == "platt") {
    m <- stats::glm(HRD_label ~ score, family = stats::binomial(), data = d)
    function(s) as.numeric(stats::predict(m, newdata = data.frame(score = s),
                                          type = "response"))
  } else if (method == "isotonic") {
    o  <- order(d$score)
    ir <- stats::isoreg(d$score[o], d$HRD_label[o])
    f  <- stats::approxfun(d$score[o], ir$yf, method = "linear", rule = 2,
                           ties = "ordered")
    function(s) pmin(pmax(f(s), 0), 1)
  } else stop("RECAL_METHOD must be 'platt' or 'isotonic'")
}

lab$risk <- lab$score
recal_note <- "raw HRD_prob (no recalibration)"
if (RECALIBRATE) {
  if (RECAL_FIT_ON_TRAIN) {
    f <- fit_recal(split_df("Train"))
    lab$risk <- f(lab$score)
    recal_note <- sprintf("%s recalibration fitted in TCGA-train, applied to both splits",
                          RECAL_METHOD)
  } else {
    for (ds in SPLITS) {
      ii <- which(lab$dataset == ds)
      lab$risk[ii] <- fit_recal(lab[ii, ])(lab$score[ii])
    }
    recal_note <- sprintf("%s recalibration fitted within each split (in-sample)",
                          RECAL_METHOD)
  }
}
message("Risk scale for panels e and g: ", recal_note)

## ---- before / after table, printed so the gain can be quoted --------------
cal_summary <- do.call(rbind, lapply(SPLITS, function(ds) {
  d <- split_df(ds)
  do.call(rbind, lapply(c("raw", "used"), function(w) {
    p  <- if (w == "raw") d$score else d$risk
    is <- cal_int_slope(d$HRD_label, p)
    hl <- hl_test(d$HRD_label, p)
    data.frame(dataset = ds, risk = w, n = nrow(d),
               brier = mean((p - d$HRD_label)^2),
               intercept = is[["intercept"]], slope = is[["slope"]],
               HL_chisq = hl$chisq, HL_df = hl$df, HL_P = hl$p)
  }))
}))
rownames(cal_summary) <- NULL
print(cal_summary, digits = 3)

# =============================================================================
# Panel A | ROC, train vs val   (raw score -- unchanged by recalibration)
# =============================================================================
roc_stat <- do.call(rbind, lapply(SPLITS, function(ds) {
  d <- split_df(ds); ci <- auc_ci(d$HRD_label, d$score)
  data.frame(dataset = ds, n = nrow(d), auc = ci[1], lo = ci[2], hi = ci[3])
}))
rownames(roc_stat) <- NULL

roc_line <- do.call(rbind, lapply(SPLITS, function(ds) {
  d <- split_df(ds); rc <- roc_curve(d$HRD_label, d$score); rc$dataset <- ds; rc
}))
roc_bandA <- do.call(rbind, lapply(SPLITS, function(ds) {
  d <- split_df(ds); bd <- roc_band(d$HRD_label, d$score); bd$dataset <- ds; bd
}))
roc_line$dataset  <- factor(roc_line$dataset,  levels = SPLITS)
roc_bandA$dataset <- factor(roc_bandA$dataset, levels = SPLITS)

lab_A <- data.frame(
  x = 0.42, y = c(0.16, 0.09),
  txt = c(sprintf("Train AUC = %s", fmt_ci(roc_stat$auc[1], roc_stat$lo[1], roc_stat$hi[1])),
          sprintf("Val AUC = %s",   fmt_ci(roc_stat$auc[2], roc_stat$lo[2], roc_stat$hi[2]))),
  col = c(COL_TRAIN, COL_VAL))

pA <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              linewidth = 0.3, colour = COL_RULE) +
  geom_ribbon(data = roc_bandA, aes(x = fpr, ymin = lo, ymax = hi, fill = dataset),
              alpha = 0.22, colour = NA) +
  geom_step(data = roc_line, aes(x = fpr, y = tpr, colour = dataset),
            linewidth = 0.45, direction = "vh") +
  geom_text(data = lab_A, aes(x = x, y = y, label = txt), colour = lab_A$col,
            hjust = 0, size = TXT_MM, family = FONT) +
  scale_colour_manual(values = COL_SET, guide = "none") +
  scale_fill_manual(values = COL_SET, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
  labs(x = "1 - specificity", y = "Sensitivity") +
  theme_fig()
pA <- head_panel(pA, "A")

# =============================================================================
# Panel B | Modality ablation -- pooled, TCGA-train only, TCGA-val only
# -----------------------------------------------------------------------------
# ablation_run() takes any subset of `lab`, keeps the patients that carry BOTH
# modality scores, and returns everything the panel needs: the three AUCs with
# bootstrap CIs, the two paired tests against the multimodal model, the ROC
# step curves and the bootstrap ROC bands. ablation_panel() turns that into a
# ggplot. Running the same function on the Train and Val subsets is what
# produces the cohort-specific versions.
#
# Two things worth stating in the legend when the split panels are used:
#   * The paired DeLong test compares AUCs measured on the same patients, so
#     it is valid within each cohort but says nothing about Train vs Val.
#   * The Val paired subset is the smaller of the two. Its bootstrap CIs are
#     wide, and a non-significant P there is an absence of evidence, not
#     evidence that the modalities are equivalent. Print n and the number of
#     HRD+ cases (the panel does this automatically via show_n = TRUE).
# =============================================================================
ABL_COLS  <- c(`WSI-only` = COL_WSI, `RNA-only` = COL_VAL, `WSI + RNA` = COL_TRAIN)
ABL_NAMES <- names(ABL_COLS)

## patients with both modality scores -> the paired subset
abl_subset <- function(d)
  d[is.finite(num(d[[SCORE_WSI]])) & is.finite(num(d[[SCORE_RNA]])), , drop = FALSE]

abl_score_list <- function(d)
  list(`WSI-only`  = num(d[[SCORE_WSI]]),
       `RNA-only`  = num(d[[SCORE_RNA]]),
       `WSI + RNA` = d$score)

ablation_run <- function(d, tag = "pooled") {
  a  <- abl_subset(d)
  n1 <- sum(a$HRD_label == 1); n0 <- sum(a$HRD_label == 0)
  message(sprintf("Ablation subset [%s]: n = %d (HRD+ %d, HRD- %d)", tag, nrow(a), n1, n0))
  if (nrow(a) < 5 || n1 < 2 || n0 < 2) {
    warning(sprintf("Ablation subset [%s] too small to plot (n = %d, HRD+ %d, HRD- %d).",
                    tag, nrow(a), n1, n0))
    return(NULL)
  }
  sc <- abl_score_list(a)
  
  stat <- do.call(rbind, lapply(ABL_NAMES, function(m) {
    ci <- auc_ci(a$HRD_label, sc[[m]])
    data.frame(model = m, auc = ci[1], lo = ci[2], hi = ci[3])
  }))
  rownames(stat) <- NULL
  stat$model <- factor(stat$model, levels = ABL_NAMES)
  
  dl <- delong(a$HRD_label, sc)
  pw <- delong_test(dl, 1, 3)[["p"]]
  pr <- delong_test(dl, 2, 3)[["p"]]
  test_used <- "Paired DeLong vs WSI + RNA"
  if (!is.finite(pw) || !is.finite(pr)) {
    message(sprintf("[%s] DeLong variance undefined; falling back to a paired bootstrap test.",
                    tag))
    if (!is.finite(pw)) pw <- boot_auc_diff_p(a$HRD_label, sc[[1]], sc[[3]])
    if (!is.finite(pr)) pr <- boot_auc_diff_p(a$HRD_label, sc[[2]], sc[[3]])
    test_used <- "Paired bootstrap vs WSI + RNA"
  }
  
  line <- do.call(rbind, lapply(ABL_NAMES, function(m) {
    rc <- roc_curve(a$HRD_label, sc[[m]]); rc$model <- m; rc
  }))
  band <- do.call(rbind, lapply(ABL_NAMES, function(m) {
    bd <- roc_band(a$HRD_label, sc[[m]]); bd$model <- m; bd
  }))
  line$model <- factor(line$model, levels = ABL_NAMES)
  band$model <- factor(band$model, levels = ABL_NAMES)
  
  list(tag = tag, data = a, scores = sc, stat = stat,
       p_wsi = pw, p_rna = pr, test = test_used,
       line = line, band = band, n = nrow(a), events = n1, nonevents = n0)
}

ablation_panel <- function(res, letter = NULL, sub = NULL, show_n = TRUE) {
  if (is.null(res)) return(NULL)
  lab_txt <- data.frame(
    x = 0.34, y = c(0.235, 0.175, 0.115),
    txt = sprintf("%-10s AUC = %s", res$stat$model,
                  fmt_ci(res$stat$auc, res$stat$lo, res$stat$hi)),
    col = unname(ABL_COLS[as.character(res$stat$model)]))
  foot <- data.frame(
    x = 0.34, y = c(0.055, 0.005),
    txt = c(res$test,
            sprintf("WSI-only P = %s;  RNA-only P = %s",
                    fmt_p_exact(res$p_wsi), fmt_p_exact(res$p_rna))))
  n_txt <- data.frame(x = 0.34, y = 0.300,
                      txt = sprintf("paired n = %d (HRD+ %d)", res$n, res$events))
  
  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "22",
                linewidth = 0.3, colour = COL_RULE) +
    geom_ribbon(data = res$band, aes(x = fpr, ymin = lo, ymax = hi, fill = model),
                alpha = 0.14, colour = NA) +
    geom_step(data = res$line, aes(x = fpr, y = tpr, colour = model),
              linewidth = 0.45, direction = "vh") +
    { if (show_n)
      geom_text(data = n_txt, aes(x = x, y = y, label = txt), hjust = 0,
                size = TXT_MM * 0.95, colour = COL_GREY, family = FONT) else NULL } +
    geom_text(data = lab_txt, aes(x = x, y = y, label = txt), colour = lab_txt$col,
              hjust = 0, size = TXT_MM, family = FONT) +
    geom_text(data = foot, aes(x = x, y = y, label = txt),
              hjust = 0, size = TXT_MM * 0.95, colour = COL_GREY, family = FONT) +
    scale_colour_manual(values = ABL_COLS, guide = "none") +
    scale_fill_manual(values = ABL_COLS, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    labs(x = "1 - specificity", y = "Sensitivity") +
    theme_fig()
  head_panel(p, letter, sub)
}

## ---- one run per cohort -----------------------------------------------------
ablR <- list(
  pooled = ablation_run(lab,                "Train + Val"),
  Train  = ablation_run(split_df("Train"),  "TCGA-train"),
  Val    = ablation_run(split_df("Val"),    "TCGA-val")
)

## objects kept under their original names so the rest of the script is unchanged
abl        <- ablR$pooled$data
abl_scores <- ablR$pooled$scores
abl_stat   <- ablR$pooled$stat
p_wsi      <- ablR$pooled$p_wsi
p_rna      <- ablR$pooled$p_rna
test_used  <- ablR$pooled$test

pB_pooled <- ablation_panel(ablR$pooled, "B", NULL,        show_n = TRUE)
pB_train  <- ablation_panel(ablR$Train,  "B", LBL[["Train"]], show_n = TRUE)
pB_val    <- ablation_panel(ablR$Val,    "B", LBL[["Val"]],   show_n = TRUE)

## the same two panels side by side, letter on the left one only
pB_split <- if (!is.null(pB_train) && !is.null(pB_val))
  patchwork::wrap_plots(pB_train,
                        ablation_panel(ablR$Val, NULL, LBL[["Val"]], show_n = TRUE),
                        nrow = 1, widths = c(1, 1)) else NULL

## which version goes into the composite Figure 1
pB <- switch(FIG1_PANEL_B,
             pooled = pB_pooled,
             train  = pB_train,
             val    = pB_val,
             stop("FIG1_PANEL_B must be 'pooled', 'train' or 'val'"))
if (is.null(pB)) {
  warning("Requested panel B version is empty; falling back to the pooled panel.")
  pB <- pB_pooled
}

# =============================================================================
# Panel C | Confusion matrices
# =============================================================================
cm_table <- function(d) {
  y <- d$HRD_label; yh <- d$pred
  cnt <- c(sum(y == 1 & yh == 1), sum(y == 1 & yh == 0),
           sum(y == 0 & yh == 1), sum(y == 0 & yh == 0))
  out <- data.frame(
    x     = c(1, 2, 1, 2),
    y     = c(2, 2, 1, 1),
    cell  = c("TP", "FN", "FP", "TN"),
    n     = cnt)
  out$pct <- 100 * out$n / rep(c(cnt[1] + cnt[2], cnt[3] + cnt[4]), each = 2)
  out$cx <- out$x + c(-0.36, 0.36, -0.36, 0.36)
  out$cy <- out$y + c(0.36, 0.36, -0.36, -0.36)
  out
}

make_cm <- function(ds, letter = NULL, show_y = TRUE, show_legend = FALSE) {
  d  <- split_df(ds); cm <- cm_table(d)
  cm$txt_col <- ifelse(cm$pct > 60, "white", COL_INK)
  p <- ggplot(cm) +
    geom_tile(aes(x = x, y = y, fill = pct), colour = "white", linewidth = 0.6) +
    geom_text(aes(x = x, y = y + 0.10, label = n, colour = txt_col),
              size = TXT_MM * 1.05, family = FONT) +
    geom_text(aes(x = x, y = y - 0.14, label = sprintf("%.1f%%", pct), colour = txt_col),
              size = TXT_MM * 0.9, family = FONT) +
    geom_text(aes(x = cx, y = cy, label = cell), colour = COL_TRAIN,
              size = TXT_MM * 0.85, family = FONT, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradient(low = "white", high = COL_TRAIN, limits = c(0, 100),
                        breaks = c(0, 50, 100), labels = c("0%", "50%", "100%"),
                        name = NULL,
                        guide = if (show_legend)
                          guide_colourbar(barwidth = unit(1.8, "mm"),
                                          barheight = unit(14, "mm"),
                                          ticks.colour = NA, frame.colour = NA)
                        else "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("HRD+", "HRD-"),
                       limits = c(0.5, 2.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(1, 2), labels = c("HRD-", "HRD+"),
                       limits = c(0.5, 2.5), expand = c(0, 0)) +
    labs(x = "Predicted", y = if (show_y) "Ground truth" else NULL) +
    theme_fig() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = if (show_legend) "right" else "none")
  head_panel(p, letter, LBL[[ds]])
}
pC <- patchwork::wrap_plots(make_cm("Train", "C", TRUE, FALSE),
                            make_cm("Val", NULL, TRUE, TRUE),
                            nrow = 1, widths = c(1, 1.18))

# =============================================================================
# Panel D | Metrics with bootstrap 95% CI
# =============================================================================
METRICS <- c("Sensitivity", "Specificity", "Precision", "F1-score", "Accuracy")

metric_ci <- function(d, B = N_BOOT) {
  y <- d$HRD_label; yh <- d$pred
  est <- metric_set(y, yh)
  bm  <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) { ib <- boot_idx(y); bm[b, ] <- metric_set(y[ib], yh[ib]) }
  data.frame(metric = names(est), est = as.numeric(est),
             lo = apply(bm, 2, stats::quantile, 0.025, na.rm = TRUE),
             hi = apply(bm, 2, stats::quantile, 0.975, na.rm = TRUE),
             row.names = NULL)
}
met <- do.call(rbind, lapply(SPLITS, function(ds) {
  m <- metric_ci(split_df(ds)); m$dataset <- ds; m
}))
met$metric  <- factor(met$metric,  levels = METRICS)
met$dataset <- factor(met$dataset, levels = SPLITS)

DODGE <- 0.68
pD <- ggplot(met, aes(x = metric, y = est, fill = dataset)) +
  geom_col(position = position_dodge(width = DODGE), width = DODGE * 0.86,
           colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.16, linewidth = 0.35,
                colour = COL_INK, position = position_dodge(width = DODGE)) +
  geom_text(aes(y = 0.055, label = sprintf("%.3f", est)),
            position = position_dodge(width = DODGE), angle = 90, hjust = 0,
            size = TXT_MM * 0.9, colour = "white", family = FONT) +
  scale_fill_manual(values = COL_SET, breaks = SPLITS, labels = unname(LBL[SPLITS])) +
  scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.25),
                     labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
  labs(x = NULL, y = "Value") +
  theme_fig() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.justification = "right",
        legend.margin = margin(0, 0, -1.2, 0, "mm"))
pD <- head_panel(pD, "D")

# =============================================================================
# Panel E | Calibration + predicted-probability histogram
# -----------------------------------------------------------------------------
# Both halves plot lab$risk, i.e. the same quantity panel G thresholds.
# Equal-frequency bins by default: with equal-width bins the extreme bins hold
# a handful of patients, so their Wilson intervals span most of the axis and
# the curve looks ragged for a reason that has nothing to do with the model.
# =============================================================================
cal_bins <- function(y, p, nbin = CAL_N_BIN, how = CAL_BINNING) {
  p <- pmin(pmax(p, 0), 1)
  br <- if (identical(how, "quantile"))
    unique(stats::quantile(p, probs = seq(0, 1, length.out = nbin + 1), na.rm = TRUE))
  else seq(0, 1, by = 1 / nbin)
  if (length(br) < 3) return(NULL)
  g <- cut(p, breaks = br, include.lowest = TRUE)
  out <- do.call(rbind, lapply(levels(g), function(lv) {
    ii <- which(g == lv); if (!length(ii)) return(NULL)
    ci <- wilson_ci(sum(y[ii] == 1), length(ii))
    data.frame(bin = lv, x = mean(p[ii]), n = length(ii),
               obs = ci[["est"]], lo = ci[["lo"]], hi = ci[["hi"]])
  }))
  rownames(out) <- NULL; out
}

cal_smooth <- function(y, p, n = 200) {
  df <- data.frame(y = as.numeric(y), p = pmin(pmax(p, 0), 1))
  fit <- try(stats::loess(y ~ p, data = df, span = 0.9, degree = 1,
                          control = stats::loess.control(surface = "direct")),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  xs <- seq(min(df$p), max(df$p), length.out = n)
  ys <- try(as.numeric(stats::predict(fit, newdata = data.frame(p = xs))), silent = TRUE)
  if (inherits(ys, "try-error")) return(NULL)
  data.frame(x = xs, y = pmin(pmax(ys, 0), 1))
}

make_cal <- function(ds, letter = NULL) {
  d  <- split_df(ds)
  p  <- d$risk
  cb <- cal_bins(d$HRD_label, p)
  sm <- if (CAL_SMOOTH) cal_smooth(d$HRD_label, p) else NULL
  brier <- mean((p - d$HRD_label)^2)
  is    <- cal_int_slope(d$HRD_label, p)
  hl    <- hl_test(d$HRD_label, p)
  
  ann <- data.frame(
    x = 0.98, y = c(0.24, 0.16, 0.08),
    txt = c(sprintf("Brier %.3f", brier),
            sprintf("Slope %.2f", is[["slope"]]),
            sprintf("HL P = %s", fmt_p_exact(hl$p))),
    col = c(COL_INK, COL_INK, COL_GREY))
  
  p_top <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "22",
                linewidth = 0.3, colour = COL_RULE) +
    { if (!is.null(sm))
      geom_line(data = sm, aes(x = x, y = y), colour = COL_SET[[ds]],
                linewidth = 0.5, alpha = 0.45) else NULL } +
    { if (!is.null(cb))
      geom_errorbar(data = cb, aes(x = x, ymin = lo, ymax = hi),
                    width = 0.02, linewidth = 0.35, colour = COL_SET[[ds]]) else NULL } +
    { if (!is.null(cb))
      geom_line(data = cb, aes(x = x, y = obs), colour = COL_SET[[ds]],
                linewidth = 0.45) else NULL } +
    { if (!is.null(cb))
      geom_point(data = cb, aes(x = x, y = obs), colour = COL_SET[[ds]],
                 size = 0.9) else NULL } +
    geom_text(data = ann, aes(x = x, y = y, label = txt, colour = col),
              hjust = 1, size = TXT_MM * 0.95, family = FONT) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    labs(x = NULL, y = "Observed frequency") +
    theme_fig() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(1.2, 1.6, 0.4, 1.2, "mm"))
  p_top <- head_panel(p_top, letter, LBL[[ds]])
  
  hb <- hist_bars(p, d$hr)
  p_bot <- ggplot() +
    geom_rect(data = hb, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Freq,
                             fill = truth), colour = NA) +
    scale_fill_manual(values = COL_TRUTH, name = "Ground truth",
                      breaks = names(COL_TRUTH)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = "Predicted probability", y = "Count") +
    theme_fig() +
    legend_inside(0.02, 0.98, c(0, 1)) +
    theme(legend.title = element_text(size = BASE_SIZE, colour = COL_INK),
          legend.key.size = unit(2.2, "mm"),
          plot.margin = margin(0.4, 1.6, 1.2, 1.2, "mm"))
  
  list(plot = p_top / p_bot + patchwork::plot_layout(heights = c(1, 0.62)),
       brier = brier, hl = hl, is = is, bins = cb)
}
calE <- list(Train = make_cal("Train", "E"), Val = make_cal("Val", NULL))
pE <- patchwork::wrap_plots(calE$Train$plot, calE$Val$plot, nrow = 1)

# =============================================================================
# Panel F | Predicted probability by ground truth (raincloud, raw score)
# -----------------------------------------------------------------------------
# Deliberately on the raw score: panel F describes the model output itself,
# and keeping it raw makes the effect of the recalibration visible when panels
# E and F are read side by side.
# =============================================================================
make_rain <- function(ds, letter = NULL, show_y = TRUE) {
  d   <- split_df(ds)
  lev <- levels(d$hr); xs <- seq_along(lev)
  poly <- do.call(rbind, lapply(xs, function(i) {
    pp <- half_violin(d$score[d$hr == lev[i]], i); if (is.null(pp)) return(NULL)
    pp$id <- i; pp$col <- COL_TRUTH[[lev[i]]]
    pp$fill <- lighten(COL_TRUTH[[lev[i]]], 0.55); pp
  }))
  bs <- do.call(rbind, lapply(xs, function(i) {
    b <- bstat(d$score[d$hr == lev[i]]); if (is.null(b)) return(NULL)
    b$x <- i; b$col <- COL_TRUTH[[lev[i]]]; b
  }))
  pts <- data.frame(x = as.numeric(d$hr) + PT_OFF + jit(nrow(d), PT_JIT),
                    y = d$score, col = COL_TRUTH[as.integer(d$hr)])
  wt <- stats::wilcox.test(score ~ hr, data = d)
  
  p <- ggplot() +
    geom_point(data = pts, aes(x = x, y = y, colour = col),
               size = 0.35, alpha = 0.45, stroke = 0, shape = 16) +
    { if (!is.null(poly))
      geom_polygon(data = poly, aes(x = x, y = y, group = id, fill = fill, colour = col),
                   linewidth = 0.3, alpha = 0.60) } +
    box_layers(bs) +
    annotate("text", x = 1.5, y = 1.05, label = fmt_p_fig(wt$p.value),
             size = TXT_MM, colour = COL_INK, family = FONT, vjust = 0.5) +
    scale_colour_identity() + scale_fill_identity() +
    scale_x_continuous(breaks = xs, labels = lev,
                       limits = c(0.45, length(lev) + 0.55), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    labs(x = NULL, y = if (show_y) "Multimodal HRD probability" else NULL) +
    theme_fig() +
    { if (!show_y) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      else NULL }
  list(plot = head_panel(p, letter, LBL[[ds]]), p = wt$p.value)
}
rainF <- list(Train = make_rain("Train", "F", TRUE), Val = make_rain("Val", NULL, FALSE))
pF <- patchwork::wrap_plots(rainF$Train$plot, rainF$Val$plot,
                            nrow = 1, widths = c(1.10, 1))

# =============================================================================
# Panel G | Decision curve analysis  (same risk scale as panel E)
# -----------------------------------------------------------------------------
# Net benefit = TP/n - (FP/n) * pt/(1-pt), treating every case whose risk
# reaches the threshold probability pt.
# =============================================================================
dca_table <- function(d, thr = seq(0.005, DCA_MAX, by = 0.005), p = NULL) {
  y <- d$HRD_label
  if (is.null(p)) p <- d$risk
  n <- length(y); prev <- mean(y == 1)
  nb <- vapply(thr, function(t) {
    pos <- p >= t
    sum(pos & y == 1) / n - (sum(pos & y == 0) / n) * (t / (1 - t))
  }, 0)
  all <- prev - (1 - prev) * thr / (1 - thr)
  rbind(data.frame(thr = thr, nb = nb,  strategy = "Model"),
        data.frame(thr = thr, nb = all, strategy = "Treat all"),
        data.frame(thr = thr, nb = 0,   strategy = "Treat none"))
}

DCA_COLS <- c(Model = COL_TRAIN, `Treat all` = COL_VAL, `Treat none` = COL_GREY)
DCA_LTY  <- c(Model = "solid",  `Treat all` = "22",    `Treat none` = "solid")

make_dca <- function(ds, letter = NULL, show_legend = TRUE) {
  d  <- split_df(ds)
  dd <- dca_table(d)
  dd$strategy <- factor(dd$strategy, levels = names(DCA_COLS))
  wide <- data.frame(thr = dd$thr[dd$strategy == "Model"],
                     model = dd$nb[dd$strategy == "Model"],
                     all   = dd$nb[dd$strategy == "Treat all"])
  wide$base <- pmax(wide$all, 0)
  wide$top  <- pmax(wide$model, wide$base)
  p <- ggplot() +
    geom_ribbon(data = wide, aes(x = thr, ymin = base, ymax = top),
                fill = FILL_TRAIN, alpha = 0.45) +
    geom_line(data = dd, aes(x = thr, y = nb, colour = strategy, linetype = strategy),
              linewidth = 0.42) +
    annotate("text", x = DCA_MAX / 2, y = 0.02,
             label = sprintf("n = %d, HRD+ %.0f%%", nrow(d), 100 * mean(d$HRD_label == 1)),
             size = TXT_MM * 0.95, colour = COL_GREY, family = FONT) +
    scale_colour_manual(values = DCA_COLS, name = NULL) +
    scale_linetype_manual(values = DCA_LTY, name = NULL) +
    scale_x_continuous(limits = c(0, DCA_MAX), breaks = seq(0, 0.8, 0.2),
                       labels = sprintf("%.1f", seq(0, 0.8, 0.2)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.03, 0.64), breaks = seq(0, 0.6, 0.2),
                       labels = sprintf("%.1f", seq(0, 0.6, 0.2)), expand = c(0, 0)) +
    labs(x = "Threshold probability", y = "Net benefit") +
    theme_fig() +
    { if (show_legend) legend_inside(0.02, 0.26, c(0, 0.5))
      else theme(legend.position = "none") } +
    theme(legend.key.width = unit(3.4, "mm"),
          legend.key.height = unit(2.4, "mm"))
  list(plot = head_panel(p, letter, LBL[[ds]]), tab = dd)
}
dcaG <- list(Train = make_dca("Train", "G", TRUE), Val = make_dca("Val", NULL, TRUE))
pG <- patchwork::wrap_plots(dcaG$Train$plot, dcaG$Val$plot, nrow = 1)

## ---- what the curve would look like on the raw score, for the record ------
for (ds in SPLITS) {
  d   <- split_df(ds)
  raw <- dca_table(d, p = d$score)
  rec <- dcaG[[ds]]$tab
  at  <- c(0.2, 0.4, 0.5, 0.6, 0.8)
  gv  <- function(tab, t) tab$nb[tab$strategy == "Model"][which.min(abs(tab$thr - t))]
  message(sprintf("Panel g  %-11s net benefit at pt = %s", LBL[[ds]],
                  paste(sprintf("%.2f", at), collapse = " / ")))
  message(sprintf("           raw HRD_prob : %s",
                  paste(sprintf("%.3f", vapply(at, function(t) gv(raw, t), 0)), collapse = " / ")))
  message(sprintf("           risk used    : %s",
                  paste(sprintf("%.3f", vapply(at, function(t) gv(rec, t), 0)), collapse = " / ")))
  message(sprintf("           patients with raw score >= 0.8: %d / %d",
                  sum(d$score >= 0.8), nrow(d)))
}

# =============================================================================
# Panel H | Subgroup AUC forest  (raw score -- unchanged by recalibration)
# =============================================================================
grade_col <- if ("grade" %in% names(lab)) "grade" else "grade_group"

SUBGROUP_SPECS <- list(
  Stage = function(d) factor(
    ifelse(trimws(d$stage_group) %in% c("I", "II"), "I-II",
           ifelse(trimws(d$stage_group) %in% c("III", "IV"), "III-IV", NA)),
    levels = c("I-II", "III-IV")),
  Residual = function(d) factor(
    ifelse(trimws(d$residual_disease) %in% c("No Macroscopic disease", "1-10 mm"), "0-10 mm",
           ifelse(trimws(d$residual_disease) %in% c("11-20 mm", ">20 mm"), ">10 mm", NA)),
    levels = c("0-10 mm", ">10 mm")),
  Age = function(d) factor(
    ifelse(is.finite(num(d$age)) & num(d$age) <= AGE_CUT, paste0("\u2264", AGE_CUT, " y"),
           ifelse(is.finite(num(d$age)), paste0(">", AGE_CUT, " y"), NA)),
    levels = c(paste0("\u2264", AGE_CUT, " y"), paste0(">", AGE_CUT, " y"))),
  Grade = function(d) factor(
    ifelse(trimws(d[[grade_col]]) %in% c("G1", "G2", "Low grade"), "G1-G2",
           ifelse(trimws(d[[grade_col]]) %in% c("G3", "G4", "High grade"), "G3-G4", NA)),
    levels = c("G1-G2", "G3-G4"))
)

subgroup_stat <- do.call(rbind, lapply(names(SUBGROUP_SPECS), function(f) {
  gv <- SUBGROUP_SPECS[[f]](lab)
  do.call(rbind, lapply(levels(gv), function(lv) {
    do.call(rbind, lapply(SPLITS, function(ds) {
      ii <- which(gv == lv & lab$dataset == ds)
      if (length(ii) < N_MIN_SUB) return(NULL)
      y <- lab$HRD_label[ii]; p <- lab$score[ii]
      if (sum(y == 1) < 1 || sum(y == 0) < 1) return(NULL)
      ci <- auc_ci(y, p)
      data.frame(feature = f, stratum = lv,
                 stratum_order = which(levels(gv) == lv), dataset = ds,
                 n = length(ii), events = sum(y == 1),
                 auc = ci[1], lo = ci[2], hi = ci[3],
                 flag_small = length(ii) < N_MIN_FLAG)
    }))
  }))
}))
rownames(subgroup_stat) <- NULL
subgroup_stat$feature <- factor(subgroup_stat$feature, levels = names(SUBGROUP_SPECS))
subgroup_stat$dataset <- factor(subgroup_stat$dataset, levels = SPLITS)

rows_h <- unique(subgroup_stat[, c("feature", "stratum", "stratum_order")])
rows_h <- rows_h[order(rows_h$feature, rows_h$stratum_order), , drop = FALSE]
rows_h$row_id <- seq_len(nrow(rows_h))
NH <- nrow(rows_h)
rows_h$y    <- NH - rows_h$row_id + 1
rows_h$band <- rows_h$row_id %% 2 == 1

SG <- merge(subgroup_stat, rows_h[, c("feature", "stratum", "y")],
            by = c("feature", "stratum"), sort = FALSE)
DODGE_H  <- 0.19
SG$ypos  <- SG$y + ifelse(SG$dataset == SPLITS[1], DODGE_H, -DODGE_H)
X_MIN_H  <- 0.50
SG$lo_c  <- pmax(SG$lo, X_MIN_H)
SG$clip  <- SG$lo < X_MIN_H
SG$txt   <- paste0(fmt_ci(SG$auc, SG$lo, SG$hi, d = 2), ifelse(SG$flag_small, "*", ""))

band_h  <- rows_h[rows_h$band, , drop = FALSE]
YLIM_H  <- c(0.4, NH + 1.30)
Y_HEAD  <- NH + 0.85

feat_h <- do.call(rbind, lapply(split(rows_h, rows_h$feature, drop = TRUE), function(g)
  data.frame(feature = as.character(g$feature[1]), y = mean(g$y))))

h_left <- ggplot() +
  geom_rect(data = band_h, aes(xmin = 0, xmax = 1, ymin = y - 0.5, ymax = y + 0.5),
            fill = STRIPE) +
  geom_text(data = feat_h, aes(x = 0.97, y = y, label = feature), hjust = 1,
            size = TXT_MM, colour = COL_INK, family = FONT, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = YLIM_H, expand = c(0, 0)) +
  labs(title = "H", subtitle = " ") +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0,
                                  colour = COL_INK, margin = margin(b = 0.8, unit = "mm")),
        plot.subtitle = element_text(size = BASE_SIZE, margin = margin(b = 1.0, unit = "mm")),
        plot.margin = margin(1.2, 0.3, 1.2, 1.2, "mm"))

leg_h <- data.frame(dataset = SPLITS, label = unname(LBL[SPLITS]),
                    x_dot = c(0.515, 0.735), x_txt = c(0.535, 0.755))
h_mid <- ggplot() +
  geom_rect(data = band_h, aes(xmin = X_MIN_H, xmax = 1, ymin = y - 0.5, ymax = y + 0.5),
            fill = STRIPE) +
  geom_vline(xintercept = c(0.75, 1.00), linetype = "22",
             linewidth = 0.25, colour = COL_RULE) +
  geom_segment(data = SG[!SG$clip, ],
               aes(x = lo_c, xend = hi, y = ypos, yend = ypos, colour = dataset),
               linewidth = 0.38, lineend = "round") +
  { if (any(SG$clip))
    geom_segment(data = SG[SG$clip, ],
                 aes(x = hi, xend = lo_c, y = ypos, yend = ypos, colour = dataset),
                 linewidth = 0.38,
                 arrow = arrow(length = unit(0.9, "mm"), type = "open", angle = 25)) } +
  geom_point(data = SG, aes(x = auc, y = ypos, colour = dataset), size = 0.9) +
  geom_point(data = leg_h, aes(x = x_dot, y = Y_HEAD, colour = dataset), size = 0.9) +
  geom_text(data = leg_h, aes(x = x_txt, y = Y_HEAD, label = label), hjust = 0,
            size = TXT_MM, colour = COL_INK, family = FONT) +
  scale_colour_manual(values = COL_SET, guide = "none") +
  scale_x_continuous(limits = c(X_MIN_H, 1), breaks = c(0.50, 0.75, 1.00),
                     labels = sprintf("%.2f", c(0.50, 0.75, 1.00)), expand = c(0, 0)) +
  scale_y_continuous(limits = YLIM_H, breaks = rows_h$y, labels = rows_h$stratum,
                     expand = c(0, 0)) +
  labs(x = "AUC (95% CI)", y = NULL, title = " ", subtitle = " ") +
  theme_fig() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(1.2, 0.3, 1.2, 0.3, "mm"))

h_right <- ggplot() +
  geom_rect(data = band_h, aes(xmin = 0, xmax = 1, ymin = y - 0.5, ymax = y + 0.5),
            fill = STRIPE) +
  geom_text(data = SG, aes(x = 0.04, y = ypos, label = txt), hjust = 0,
            size = TXT_MM, colour = COL_INK, family = FONT) +
  annotate("text", x = 0.04, y = Y_HEAD, label = "AUC (95% CI)", hjust = 0,
           size = TXT_MM, colour = COL_INK, family = FONT) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = YLIM_H, expand = c(0, 0)) +
  labs(title = " ", subtitle = " ") +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(size = BASE_SIZE, margin = margin(b = 0.8, unit = "mm")),
        plot.subtitle = element_text(size = BASE_SIZE, margin = margin(b = 1.0, unit = "mm")),
        plot.margin = margin(1.2, 1.2, 1.2, 0.3, "mm"))

pH <- patchwork::wrap_plots(h_left, h_mid, h_right,
                            nrow = 1, widths = c(0.42, 1, 0.80))

# =============================================================================
# Assemble
# =============================================================================
row1_left  <- patchwork::wrap_plots(pA, pB, nrow = 1, widths = c(1, 1))
row1_right <- patchwork::wrap_plots(pC, pD, ncol = 1, heights = c(1.30, 1))
row1 <- patchwork::wrap_plots(row1_left, row1_right, nrow = 1, widths = c(2, 1.75))
row2 <- patchwork::wrap_plots(pE, pF, nrow = 1, widths = c(1.30, 1))
row3 <- patchwork::wrap_plots(pG, pH, nrow = 1, widths = c(1.30, 1))
FIG1 <- patchwork::wrap_plots(row1, row2, row3, ncol = 1, heights = c(1.15, 1, 1))

save_plot(FIG1, "Figure1_model_performance", width = FIG_W, height = FIG_H)

save_plot(pA, "Fig1A_ROC_TrainVal",        width = 60,  height = 58)
save_plot(pB, "Fig1B_ROC_Ablation",        width = 60,  height = 58)
save_plot(pC, "Fig1C_ConfusionMatrix",     width = 85,  height = 52)
save_plot(pD, "Fig1D_Metrics",             width = 85,  height = 45)
save_plot(pE, "Fig1E_Calibration",         width = 105, height = 72)
save_plot(pF, "Fig1F_ProbabilityByTruth",  width = 80,  height = 72)
save_plot(pG, "Fig1G_DecisionCurve",       width = 105, height = 66)
save_plot(pH, "Fig1H_SubgroupForest",      width = 80,  height = 66)

## ---- the three versions of panel B, each on its own canvas ----------------
if (!is.null(pB_pooled))
  save_plot(pB_pooled, "Fig1B_ROC_Ablation_Pooled", width = 60,  height = 58)
if (!is.null(pB_train))
  save_plot(pB_train,  "Fig1B_ROC_Ablation_Train",  width = 60,  height = 58)
if (!is.null(pB_val))
  save_plot(pB_val,    "Fig1B_ROC_Ablation_Val",    width = 60,  height = 58)
if (!is.null(pB_split))
  save_plot(pB_split,  "Fig1B_ROC_Ablation_TrainVal_side_by_side",
            width = 120, height = 58)

# =============================================================================
# Source data export  (one table per panel with statistics)
# =============================================================================
write.csv(roc_stat, file.path(SRC_DIR, "Fig1a_AUC.csv"), row.names = FALSE)
write.csv(roc_line, file.path(SRC_DIR, "Fig1a_ROC_coordinates.csv"), row.names = FALSE)

abl_out <- abl_stat
abl_out$p_vs_multimodal <- c(p_wsi, p_rna, NA)
abl_out$test   <- test_used
abl_out$n      <- nrow(abl)
abl_out$events <- sum(abl$HRD_label == 1)
write.csv(abl_out, file.path(SRC_DIR, "Fig1b_ablation.csv"), row.names = FALSE)

## ---- ablation statistics for all three cohorts in one table ---------------
abl_all <- do.call(rbind, lapply(names(ablR), function(k) {
  r <- ablR[[k]]; if (is.null(r)) return(NULL)
  z <- r$stat
  z$cohort          <- switch(k, pooled = "Train + Val", Train = "TCGA-train", Val = "TCGA-val")
  z$p_vs_multimodal <- c(r$p_wsi, r$p_rna, NA)
  z$test            <- r$test
  z$n               <- r$n
  z$events          <- r$events
  z$nonevents       <- r$nonevents
  z[, c("cohort", "model", "n", "events", "nonevents",
        "auc", "lo", "hi", "p_vs_multimodal", "test")]
}))
rownames(abl_all) <- NULL
write.csv(abl_all, file.path(SRC_DIR, "Fig1b_ablation_by_cohort.csv"), row.names = FALSE)

## ROC coordinates for the two cohort-specific panels
abl_roc_out <- do.call(rbind, lapply(c("Train", "Val"), function(k) {
  r <- ablR[[k]]; if (is.null(r)) return(NULL)
  z <- r$line; z$cohort <- LBL[[k]]; z
}))
if (!is.null(abl_roc_out))
  write.csv(abl_roc_out, file.path(SRC_DIR, "Fig1b_ablation_ROC_coordinates_by_cohort.csv"),
            row.names = FALSE)

cm_out <- do.call(rbind, lapply(SPLITS, function(ds) {
  z <- cm_table(split_df(ds))[, c("cell", "n", "pct")]; z$dataset <- ds; z }))
write.csv(cm_out, file.path(SRC_DIR, "Fig1c_confusion.csv"), row.names = FALSE)
write.csv(met,    file.path(SRC_DIR, "Fig1d_metrics.csv"),   row.names = FALSE)

cal_out <- do.call(rbind, lapply(SPLITS, function(ds) {
  z <- calE[[ds]]$bins; if (is.null(z)) return(NULL)
  z$dataset   <- ds
  z$brier     <- calE[[ds]]$brier
  z$intercept <- calE[[ds]]$is[["intercept"]]
  z$slope     <- calE[[ds]]$is[["slope"]]
  z$HL_P      <- calE[[ds]]$hl$p
  z$HL_df     <- calE[[ds]]$hl$df
  z$risk_used <- recal_note
  z }))
write.csv(cal_out, file.path(SRC_DIR, "Fig1e_calibration.csv"), row.names = FALSE)
write.csv(cal_summary, file.path(SRC_DIR, "Fig1e_calibration_before_after.csv"),
          row.names = FALSE)

write.csv(data.frame(sample_id = lab$sample_id, dataset = lab$dataset,
                     HRD_label = lab$HRD_label, HRD_prob = lab$score,
                     risk_used = lab$risk, pred = lab$pred),
          file.path(SRC_DIR, "Fig1ef_sample_level.csv"), row.names = FALSE)

dca_out <- do.call(rbind, lapply(SPLITS, function(ds) {
  z <- dcaG[[ds]]$tab; z$dataset <- ds; z$risk_used <- recal_note; z }))
write.csv(dca_out, file.path(SRC_DIR, "Fig1g_decision_curve.csv"), row.names = FALSE)

write.csv(subgroup_stat, file.path(SRC_DIR, "Fig1h_subgroup_AUC.csv"), row.names = FALSE)

# =============================================================================
# Console summary -- the exact numbers to quote in the legend / Results
# =============================================================================
message("\nFigure 1 written to: ", OUT_DIR)
for (i in seq_len(nrow(roc_stat)))
  message(sprintf("Panel a  %-11s n = %3d  AUC = %s",
                  LBL[[as.character(roc_stat$dataset[i])]], roc_stat$n[i],
                  fmt_ci(roc_stat$auc[i], roc_stat$lo[i], roc_stat$hi[i])))

for (k in names(ablR)) {
  r <- ablR[[k]]; if (is.null(r)) next
  message(sprintf("Panel b  [%-11s] paired n = %d (HRD+ %d);  %s;  WSI-only P = %s, RNA-only P = %s",
                  r$tag, r$n, r$events, r$test,
                  fmt_p_exact(r$p_wsi), fmt_p_exact(r$p_rna)))
  for (i in seq_len(nrow(r$stat)))
    message(sprintf("           %-10s AUC = %s", r$stat$model[i],
                    fmt_ci(r$stat$auc[i], r$stat$lo[i], r$stat$hi[i])))
}
message(sprintf("Panel b embedded in Figure 1: '%s' version.", FIG1_PANEL_B))

for (ds in SPLITS)
  message(sprintf(paste0("Panel e  %-11s Brier = %.3f;  intercept = %.2f, slope = %.2f;  ",
                         "Hosmer-Lemeshow chi2 = %.1f, df = %d, P = %s"),
                  LBL[[ds]], calE[[ds]]$brier, calE[[ds]]$is[["intercept"]],
                  calE[[ds]]$is[["slope"]], calE[[ds]]$hl$chisq,
                  calE[[ds]]$hl$df, fmt_p_exact(calE[[ds]]$hl$p)))
for (ds in SPLITS)
  message(sprintf("Panel f  %-11s Mann-Whitney P = %s (exact value for the legend)",
                  LBL[[ds]], fmt_p_exact(rainF[[ds]]$p)))
message("Panels e and g risk scale: ", recal_note)
if (RECALIBRATE && RECAL_FIT_ON_TRAIN)
  message("  -> the TCGA-train calibration curve is in-sample and lies on the ",
          "diagonal by construction; only the TCGA-val curve is out-of-sample.")
if (RECALIBRATE && !RECAL_FIT_ON_TRAIN)
  message("  -> WARNING: recalibration fitted inside each split. The validation ",
          "calibration is not out-of-sample and must not be reported as such.")
message(sprintf("Panel h  %d subgroup rows; %d flagged with '*' (n < %d)",
                nrow(rows_h), sum(subgroup_stat$flag_small), N_MIN_FLAG))

# =============================================================================
# Panel E, lower half -- alternative rendering
#
# Replaces the dodged histogram with a per-patient waterfall: patients are
# sorted by predicted risk and each bar runs from a common baseline to that
# patient's risk, coloured by ground truth. Two cohorts side by side, exported
# on a single canvas of 198.7626 pt x 153.9204 pt (width x height).
#
#   198.7626 pt / 72 = 2.760592 in   (= 70.12 mm)
#   153.9204 pt / 72 = 2.137783 in   (= 54.30 mm)
#
# The canvas is given to ggsave() in inches so the pt figures are exact rather
# than rounded through millimetres.
#
# Every object created here is prefixed WF_ / wf_ so the Figure 1 objects above
# are untouched.
# =============================================================================

## ---- switches ---------------------------------------------------------------
WF_SCORE      <- "HRD_prob"    # "risk" = the scale panels E/G share; "score" / "HRD_prob" = raw
WF_BASELINE   <- "threshold"   # "threshold" | "median" | a numeric value, e.g. 0.5
WF_BAR_W      <- 0.80          # bar width in patient-index units (1 = bars touch)
WF_COLS       <- COL_TRUTH     # to copy the reference figure exactly, use:
#   c(`HRD-` = "#4CB5AB", `HRD+` = "#EFC050")
WF_SHOW_LABEL <- FALSE         # TRUE -> print the cohort name above each subplot
WF_BASE_SIZE  <- BASE_SIZE     # font size; shrink to 7 if the labels feel crowded
WF_YLAB       <- if (identical(WF_SCORE, "risk")) "Predicted probability" else "Score"
WF_W_PT       <- 198.7626      # canvas width  (pt)
WF_H_PT       <- 153.9204      # canvas height (pt)
WF_PNG        <- TRUE          # also write a 600 dpi PNG next to the PDF
WF_FILE       <- "Fig1E_waterfall_risk"

WF_TXT_MM <- WF_BASE_SIZE / 2.845276

## ---- baseline ---------------------------------------------------------------
## The reference figure uses ONE line for both cohorts, so it is estimated on
## TCGA-train only and then applied unchanged to TCGA-val.
wf_thr_from_pred <- function(d, v) {
  if (!any(d$pred == 1, na.rm = TRUE) || !any(d$pred == 0, na.rm = TRUE))
    return(NA_real_)
  lo <- suppressWarnings(max(v[d$pred == 0], na.rm = TRUE))
  hi <- suppressWarnings(min(v[d$pred == 1], na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(NA_real_)
  (lo + hi) / 2
}
wf_thr_youden <- function(d, v) {
  cand <- sort(unique(v[is.finite(v)]))
  if (!length(cand)) return(NA_real_)
  j <- vapply(cand, function(t) {
    mm <- metric_set(d$HRD_label, as.integer(v >= t))
    unname(mm[["Sensitivity"]] + mm[["Specificity"]] - 1)
  }, 0)
  cand[which.max(j)]
}

wf_tr  <- split_df("Train")
wf_vtr <- wf_tr[[WF_SCORE]]
WF_BASE_VAL <- if (is.numeric(WF_BASELINE)) {
  WF_BASELINE
} else if (identical(WF_BASELINE, "median")) {
  stats::median(wf_vtr, na.rm = TRUE)
} else {
  b <- wf_thr_from_pred(wf_tr, wf_vtr)
  if (is.finite(b)) b else wf_thr_youden(wf_tr, wf_vtr)
}
if (!is.finite(WF_BASE_VAL)) WF_BASE_VAL <- stats::median(wf_vtr, na.rm = TRUE)
message(sprintf("Waterfall baseline (%s, fitted on TCGA-train): %.4f  [scale: %s]",
                if (is.numeric(WF_BASELINE)) "fixed" else WF_BASELINE,
                WF_BASE_VAL, WF_SCORE))

## ---- data + axis helpers ----------------------------------------------------
wf_df <- function(ds) {
  d <- split_df(ds)
  v <- d[[WF_SCORE]]
  keep <- is.finite(v)
  d <- d[keep, , drop = FALSE]; v <- v[keep]
  o <- order(v)                                    # ascending, as in the reference
  data.frame(idx     = seq_along(v),
             val     = v[o],
             truth   = factor(as.character(d$hr[o]), levels = names(WF_COLS)),
             dataset = ds,
             stringsAsFactors = FALSE)
}

## nice, human-readable y limits: pick a step, then round outwards to it
wf_lims <- function(v, base) {
  rng  <- range(c(v, base), na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) {
    step <- 0.05
    rng  <- c(base - step, base + step); span <- diff(rng)
  }
  step <- if (span > 0.5) 0.25 else if (span > 0.2) 0.10 else 0.05
  lo <- floor((rng[1] - 0.03 * span) / step) * step
  hi <- ceiling((rng[2] + 0.03 * span) / step) * step
  brk <- seq(lo, hi, by = step)
  list(lim = c(lo, hi), brk = brk, lab = sprintf("%.2f", brk))
}

## ---- one cohort -------------------------------------------------------------
make_wf <- function(ds, show_legend = TRUE, show_y = TRUE) {
  w <- wf_df(ds)
  L <- wf_lims(w$val, WF_BASE_VAL)
  p <- ggplot(w) +
    geom_rect(aes(xmin = idx - WF_BAR_W / 2, xmax = idx + WF_BAR_W / 2,
                  ymin = pmin(val, WF_BASE_VAL), ymax = pmax(val, WF_BASE_VAL),
                  fill = truth), colour = NA) +
    geom_hline(yintercept = WF_BASE_VAL, linewidth = 0.35, colour = COL_INK) +
    scale_fill_manual(values = WF_COLS, breaks = rev(names(WF_COLS)), name = NULL) +
    scale_x_continuous(limits = c(0.5, nrow(w) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = L$lim, breaks = L$brk, labels = L$lab,
                       expand = c(0, 0)) +
    labs(x = "Patient index", y = if (show_y) WF_YLAB else NULL) +
    theme_fig(WF_BASE_SIZE) +
    { if (show_legend) legend_inside(0.98, 0.04, c(1, 0))
      else theme(legend.position = "none") } +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key.size = unit(2.2, "mm"),
          legend.spacing.y = unit(0.2, "mm"),
          plot.margin  = margin(1.0, 1.6, 1.0, 1.2, "mm"))
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  if (WF_SHOW_LABEL) p <- head_panel(p, NULL, LBL[[ds]])
  list(plot = p, data = w, lims = L$lim)
}

wfE <- list(Train = make_wf("Train", TRUE, TRUE),
            Val   = make_wf("Val",   TRUE, TRUE))

WF_FIG <- patchwork::wrap_plots(wfE$Train$plot, wfE$Val$plot,
                                nrow = 1, widths = c(1, 1))

## ---- export at exactly 198.7626 pt x 153.9204 pt ----------------------------
WF_W_IN <- WF_W_PT / 72
WF_H_IN <- WF_H_PT / 72

ggplot2::ggsave(file.path(OUT_DIR, paste0(WF_FILE, ".pdf")), plot = WF_FIG,
                width = WF_W_IN, height = WF_H_IN, units = "in", device = "pdf")
if (WF_PNG)
  ggplot2::ggsave(file.path(OUT_DIR, paste0(WF_FILE, ".png")), plot = WF_FIG,
                  width = WF_W_IN, height = WF_H_IN, units = "in", dpi = 600)

## single-cohort versions on half the canvas, in case they are laid out separately
for (ds in SPLITS)
  ggplot2::ggsave(file.path(OUT_DIR, sprintf("%s_%s.pdf", WF_FILE, ds)),
                  plot = wfE[[ds]]$plot,
                  width = WF_W_IN / 2, height = WF_H_IN, units = "in", device = "pdf")

## ---- source data ------------------------------------------------------------
wf_out <- do.call(rbind, lapply(SPLITS, function(ds) {
  z <- wfE[[ds]]$data
  z$baseline  <- WF_BASE_VAL
  z$delta     <- z$val - WF_BASE_VAL
  z$scale     <- WF_SCORE
  z$risk_used <- recal_note
  z
}))
write.csv(wf_out, file.path(SRC_DIR, "Fig1e_waterfall.csv"), row.names = FALSE)

message(sprintf("Waterfall written: %s.pdf  (%.4f x %.4f pt)",
                file.path(OUT_DIR, WF_FILE), WF_W_PT, WF_H_PT))
for (ds in SPLITS) {
  z <- wfE[[ds]]$data
  message(sprintf("  %-11s n = %3d | above baseline %3d (HRD+ %3d) | below %3d (HRD+ %3d) | y-axis %.2f-%.2f",
                  LBL[[ds]], nrow(z),
                  sum(z$val >= WF_BASE_VAL), sum(z$val >= WF_BASE_VAL & z$truth == "HRD+"),
                  sum(z$val <  WF_BASE_VAL), sum(z$val <  WF_BASE_VAL & z$truth == "HRD+"),
                  wfE[[ds]]$lims[1], wfE[[ds]]$lims[2]))
}


# =============================================================================
# Panel B2 | Modality comparison on each modality's MAXIMUM available dataset
# -----------------------------------------------------------------------------
# Panel B above answers "given the same patients, does adding a modality help?".
# This block answers a different question: "how well does each model do on
# every patient it can actually score?". WSI-only therefore runs on all
# patients carrying a WSI score, RNA-only on all patients carrying an RNA
# score, and WSI + RNA on all patients carrying a multimodal score. The three
# curves are built from DIFFERENT, only partly overlapping sets of patients.
#
# Three consequences that have to be handled, not glossed over:
#
#   1. NO PAIRED TEST IS AVAILABLE. DeLong's paired test assumes the AUCs were
#      measured on the same cases; here they were not, so delong_test() is not
#      used in this block. The comparison uses a bootstrap that resamples the
#      UNION of patients once per replicate and then recomputes each AUC on
#      whichever resampled patients carry that modality's score. Because the
#      shared patients are resampled jointly, the correlation induced by the
#      overlap is preserved -- which a naive two-sample test would ignore and
#      an ordinary paired test could not represent at all.
#
#   2. THE COMPARISON IS CONFOUNDED BY CASE MIX. A modality with more patients
#      is being evaluated on a different, usually broader, case mix. If the
#      extra patients are systematically easier or harder, part of any AUC gap
#      is cohort composition rather than modality. This is why the paired panel
#      B remains the primary ablation result and this panel is supporting
#      evidence about real-world coverage. Say so in the legend.
#
#   3. n MUST BE ON THE FIGURE. Every curve carries its own denominator, so
#      the panel prints n for all three models and the console prints the
#      overlap table (how many patients have WSI only, RNA only, or both).
#
# Everything here is prefixed MAXN_ / maxn_.
# =============================================================================
MAXN_GET <- list(
  `WSI-only`  = function(d) num(d[[SCORE_WSI]]),
  `RNA-only`  = function(d) num(d[[SCORE_RNA]]),
  `WSI + RNA` = function(d) d$score
)
MAXN_NAMES <- names(MAXN_GET)
MAXN_COLS  <- ABL_COLS[MAXN_NAMES]

MAXN_W  <- 60    # standalone panel width  (mm)
MAXN_H  <- 58    # standalone panel height (mm)

## ---- bootstrap on the union of patients, for partly overlapping samples ----
## One resample of the union per replicate; each AUC is then recomputed on the
## resampled patients that carry the relevant score. auc_fast() already drops
## non-finite scores, so the subsetting is implicit.
maxn_boot_diff_p <- function(d, m1, m2, B = N_BOOT) {
  y  <- d$HRD_label
  v1 <- MAXN_GET[[m1]](d); v2 <- MAXN_GET[[m2]](d)
  keep <- !is.na(y) & (is.finite(v1) | is.finite(v2))
  y <- y[keep]; v1 <- v1[keep]; v2 <- v2[keep]
  enough <- function(v) sum(is.finite(v) & y == 1) >= 2 && sum(is.finite(v) & y == 0) >= 2
  if (!enough(v1) || !enough(v2)) return(NA_real_)
  i1 <- which(y == 1); i0 <- which(y == 0)
  dd <- vapply(seq_len(B), function(b) {
    ib <- c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE))
    auc_fast(y[ib], v1[ib]) - auc_fast(y[ib], v2[ib])
  }, 0)
  dd <- dd[is.finite(dd)]
  if (!length(dd)) return(NA_real_)
  min(1, 2 * min(mean(dd <= 0), mean(dd >= 0)))
}

## ---- statistics for one cohort ---------------------------------------------
maxn_run <- function(d, tag = "Train + Val") {
  y <- d$HRD_label
  per <- lapply(MAXN_NAMES, function(m) {
    v    <- MAXN_GET[[m]](d)
    keep <- is.finite(v) & !is.na(y)
    list(m = m, y = y[keep], v = v[keep],
         n = sum(keep), events = sum(y[keep] == 1), nonevents = sum(y[keep] == 0))
  })
  names(per) <- MAXN_NAMES
  
  bad <- vapply(per, function(z) z$events < 2 || z$nonevents < 2, logical(1))
  if (any(bad)) {
    warning(sprintf("Max-n panel [%s]: %s has too few cases; panel skipped.",
                    tag, paste(MAXN_NAMES[bad], collapse = ", ")))
    return(NULL)
  }
  
  stat <- do.call(rbind, lapply(per, function(z) {
    ci <- auc_ci(z$y, z$v)
    data.frame(model = z$m, n = z$n, events = z$events, nonevents = z$nonevents,
               auc = ci[1], lo = ci[2], hi = ci[3])
  }))
  rownames(stat) <- NULL
  stat$model <- factor(stat$model, levels = MAXN_NAMES)
  
  line <- do.call(rbind, lapply(per, function(z) {
    rc <- roc_curve(z$y, z$v); rc$model <- z$m; rc
  }))
  band <- do.call(rbind, lapply(per, function(z) {
    bd <- roc_band(z$y, z$v); bd$model <- z$m; bd
  }))
  line$model <- factor(line$model, levels = MAXN_NAMES)
  band$model <- factor(band$model, levels = MAXN_NAMES)
  
  pw <- maxn_boot_diff_p(d, "WSI-only", "WSI + RNA")
  pr <- maxn_boot_diff_p(d, "RNA-only", "WSI + RNA")
  
  ## coverage / overlap, for the legend
  hasW <- is.finite(MAXN_GET[["WSI-only"]](d)) & !is.na(y)
  hasR <- is.finite(MAXN_GET[["RNA-only"]](d)) & !is.na(y)
  cov <- c(both = sum(hasW & hasR), wsi_only = sum(hasW & !hasR),
           rna_only = sum(!hasW & hasR), neither = sum(!hasW & !hasR))
  
  message(sprintf("Max-n subset [%s]: WSI n = %d | RNA n = %d | WSI+RNA n = %d | both = %d",
                  tag, stat$n[1], stat$n[2], stat$n[3], cov[["both"]]))
  
  list(tag = tag, stat = stat, line = line, band = band,
       p_wsi = pw, p_rna = pr, coverage = cov,
       test = "Bootstrap on the union (unpaired samples)")
}

## ---- panel ------------------------------------------------------------------
maxn_panel <- function(res, letter = NULL, sub = NULL) {
  if (is.null(res)) return(NULL)
  lab_txt <- data.frame(
    x = 0.32, y = c(0.245, 0.185, 0.125),
    txt = sprintf("%-10s AUC = %s", res$stat$model,
                  fmt_ci(res$stat$auc, res$stat$lo, res$stat$hi)),
    col = unname(MAXN_COLS[as.character(res$stat$model)]))
  foot <- data.frame(
    x = 0.32, y = c(0.065, 0.015),
    txt = c(sprintf("n = %d / %d / %d  (WSI / RNA / WSI+RNA)",
                    res$stat$n[1], res$stat$n[2], res$stat$n[3]),
            sprintf("unpaired bootstrap: WSI P = %s;  RNA P = %s",
                    fmt_p_exact(res$p_wsi), fmt_p_exact(res$p_rna))))
  
  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "22",
                linewidth = 0.3, colour = COL_RULE) +
    geom_ribbon(data = res$band, aes(x = fpr, ymin = lo, ymax = hi, fill = model),
                alpha = 0.14, colour = NA) +
    geom_step(data = res$line, aes(x = fpr, y = tpr, colour = model),
              linewidth = 0.45, direction = "vh") +
    geom_text(data = lab_txt, aes(x = x, y = y, label = txt), colour = lab_txt$col,
              hjust = 0, size = TXT_MM * 0.95, family = FONT) +
    geom_text(data = foot, aes(x = x, y = y, label = txt),
              hjust = 0, size = TXT_MM * 0.85, colour = COL_GREY, family = FONT) +
    scale_colour_manual(values = MAXN_COLS, guide = "none") +
    scale_fill_manual(values = MAXN_COLS, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = sprintf("%.2f", seq(0, 1, 0.25)), expand = c(0, 0)) +
    labs(x = "1 - specificity", y = "Sensitivity") +
    theme_fig()
  head_panel(p, letter, sub)
}

## ---- one run per cohort -----------------------------------------------------
maxnR <- list(
  Train  = maxn_run(split_df("Train"), "TCGA-train"),
  Val    = maxn_run(split_df("Val"),   "TCGA-val"),
  pooled = maxn_run(lab,               "TCGA-train + TCGA-val")
)

pB2_train  <- maxn_panel(maxnR$Train,  "B", LBL[["Train"]])
pB2_val    <- maxn_panel(maxnR$Val,    "B", LBL[["Val"]])
pB2_pooled <- maxn_panel(maxnR$pooled, "B", "TCGA-train + TCGA-val")

## the three of them on one canvas, letter only on the first
pB2_all <- {
  ps <- list(maxn_panel(maxnR$Train,  "B",  LBL[["Train"]]),
             maxn_panel(maxnR$Val,    NULL, LBL[["Val"]]),
             maxn_panel(maxnR$pooled, NULL, "TCGA-train + TCGA-val"))
  ps <- Filter(Negate(is.null), ps)
  if (length(ps)) patchwork::wrap_plots(ps, nrow = 1) else NULL
}

## ---- export ------------------------------------------------------------------
if (!is.null(pB2_train))
  save_plot(pB2_train,  "Fig1B2_ROC_MaxN_Train",     width = MAXN_W, height = MAXN_H)
if (!is.null(pB2_val))
  save_plot(pB2_val,    "Fig1B2_ROC_MaxN_Val",       width = MAXN_W, height = MAXN_H)
if (!is.null(pB2_pooled))
  save_plot(pB2_pooled, "Fig1B2_ROC_MaxN_TrainVal",  width = MAXN_W, height = MAXN_H)
if (!is.null(pB2_all))
  save_plot(pB2_all,    "Fig1B2_ROC_MaxN_all_three",
            width = MAXN_W * 3, height = MAXN_H)

## ---- source data --------------------------------------------------------------
maxn_out <- do.call(rbind, lapply(names(maxnR), function(k) {
  r <- maxnR[[k]]; if (is.null(r)) return(NULL)
  z <- r$stat
  z$cohort            <- r$tag
  z$p_vs_multimodal   <- c(r$p_wsi, r$p_rna, NA)
  z$test              <- r$test
  z$n_both_modalities <- r$coverage[["both"]]
  z$n_WSI_only        <- r$coverage[["wsi_only"]]
  z$n_RNA_only        <- r$coverage[["rna_only"]]
  z[, c("cohort", "model", "n", "events", "nonevents", "auc", "lo", "hi",
        "p_vs_multimodal", "test",
        "n_both_modalities", "n_WSI_only", "n_RNA_only")]
}))
rownames(maxn_out) <- NULL
write.csv(maxn_out, file.path(SRC_DIR, "Fig1b2_maxN_ablation_by_cohort.csv"),
          row.names = FALSE)

maxn_roc_out <- do.call(rbind, lapply(names(maxnR), function(k) {
  r <- maxnR[[k]]; if (is.null(r)) return(NULL)
  z <- r$line; z$cohort <- r$tag; z
}))
if (!is.null(maxn_roc_out))
  write.csv(maxn_roc_out, file.path(SRC_DIR, "Fig1b2_maxN_ROC_coordinates.csv"),
            row.names = FALSE)

## ---- console summary ----------------------------------------------------------
message("\nPanel B2 -- each modality on its maximum available dataset")
for (k in names(maxnR)) {
  r <- maxnR[[k]]; if (is.null(r)) next
  message(sprintf("  [%s]  coverage: both %d | WSI only %d | RNA only %d | neither %d",
                  r$tag, r$coverage[["both"]], r$coverage[["wsi_only"]],
                  r$coverage[["rna_only"]], r$coverage[["neither"]]))
  for (i in seq_len(nrow(r$stat)))
    message(sprintf("            %-10s n = %3d (HRD+ %3d)  AUC = %s",
                    r$stat$model[i], r$stat$n[i], r$stat$events[i],
                    fmt_ci(r$stat$auc[i], r$stat$lo[i], r$stat$hi[i])))
  message(sprintf("            %s: WSI-only P = %s;  RNA-only P = %s",
                  r$test, fmt_p_exact(r$p_wsi), fmt_p_exact(r$p_rna)))
}
message("  NOTE: the three curves come from different patient sets. Report n per ",
        "curve, and treat the paired panel B as the primary ablation result.")
