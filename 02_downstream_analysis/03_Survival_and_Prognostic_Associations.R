# =============================================================================
# Figure 3 | Survival stratification and prognostic modelling
# =============================================================================
# Scientific scope: Kaplan-Meier curves, Cox regression and landmark progression.
# Usage: Rscript 03_Survival_and_Prognostic_Associations.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: MIRROR_CLINICAL_CSV (KM/landmark) and MIRROR_COX_CSV (forest).
# Outputs: survival_panels/ and cox_forest/ under MIRROR_OUTPUT_DIR.
# Dependencies: dplyr, survival, survminer, ggplot2, forestploter.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# The Cox source used a distinct clinical export in the supplied script.
# =============================================================================

KM_DATA_FILE <- Sys.getenv("MIRROR_CLINICAL_CSV", unset = "<PATH_TO_CLINICAL_CSV>")
COX_DATA_FILE <- Sys.getenv("MIRROR_COX_CSV", unset = "<PATH_TO_COX_CLINICAL_CSV>")
OUT_DIR <- Sys.getenv("MIRROR_OUTPUT_DIR", unset = "<PATH_TO_OUTPUT_DIRECTORY>")
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
assert_input_file(KM_DATA_FILE, "KM_DATA_FILE")
assert_input_file(COX_DATA_FILE, "COX_DATA_FILE")
prepare_output_dir(OUT_DIR, "OUT_DIR")
KM_OUT_DIR <- file.path(OUT_DIR, "survival_panels")
COX_OUT_DIR <- file.path(OUT_DIR, "cox_forest")
prepare_output_dir(KM_OUT_DIR, "KM_OUT_DIR")
prepare_output_dir(COX_OUT_DIR, "COX_OUT_DIR")



# Preserve analysis state across sections.

required_km <- c("dplyr", "survival", "survminer", "ggplot2")
missing_km <- required_km[!vapply(required_km, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_km)) stop("Install required packages: ", paste(missing_km, collapse = ", "))
suppressPackageStartupMessages(invisible(lapply(required_km, library, character.only = TRUE)))


KaplanMeier_Plot <- function(PD, Y_labels, titles, type = "l") {
  valid <- is.finite(suppressWarnings(as.numeric(PD$OS.time))) &
    PD$OS %in% c(0, 1) & !is.na(PD$Pre_Label)
  PD <- PD[valid, , drop = FALSE]
  if (nrow(PD) < 5L || nlevels(droplevels(PD$Pre_Label)) < 2L ||
      sum(PD$OS == 1) < 1L)
    stop("Cannot estimate survival for ", titles,
         ": at least five valid observations, two HRD groups and an event are required.")
  height_ratio <- 0.75
  base_font_size <- 8
  line_size <- 0.5
  table_height_ratio <- 0.15
  text_mm <- base_font_size / 2.83464567
  censor_mm <- 0.7 * text_mm

  scale_factor <- 0.75
  curve_pt <- 1.5 / scale_factor
  curve_mm <- curve_pt * 0.3527778  # pt -> mm

  axis_pt_final <- 0.5
  axis_pt_export <- axis_pt_final / scale_factor
  axis_mm <- axis_pt_export * 0.3527778            # pt -> mm
  

  text_mm <- base_font_size / 2.83464567
  fit <- survfit(Surv(OS.time, OS) ~ Pre_Label, data = PD)
  P0 <- ggsurvplot(
    fit = fit,
    title = titles,
    risk.table = TRUE,
    ggtheme = theme_classic(base_size = base_font_size),
    risk.table.col = "strata",
    palette = c( "#416993","#AE6534"),
    pval = TRUE,
    pval.method = TRUE,
    

    pval.size = text_mm,
    pval.method.size = text_mm,
    
    ylab = Y_labels,
    xlab = "Time (years)",
    tables.height = 0.25,
    font.x = c(base_font_size, "plain", "black"),
    font.y = c(base_font_size, "plain", "black"),
    font.tickslab = c(base_font_size, "plain", "black"),
    ylim = c(0, 1),
    size = curve_mm,
    legend.labs = c("HRD-", "HRD+"),
    risk.table.title = "Number at risk",
    conf.int = TRUE
  )
  

  P0$plot <- P0$plot +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = base_font_size),
      axis.title = element_text(size = base_font_size),
      axis.text  = element_text(size = base_font_size),
      legend.title = element_text(size = base_font_size),
      legend.text  = element_text(size = base_font_size)
    )
  

  res_cox <- tryCatch(coxph(Surv(OS.time, OS) ~ Pre_Label, data = PD),
                      error = function(e) NULL)
  if (!is.null(res_cox)) {
    ci <- summary(res_cox)$conf.int
    if (nrow(ci) && all(is.finite(ci[1, c(1, 3, 4)])))
      P0$plot <- P0$plot + annotate("text", x = 4.5, y = 0.12,
        label = sprintf("HR = %.3f (%.3f-%.3f)", ci[1, 1], ci[1, 3], ci[1, 4]),
        size = text_mm, hjust = 0)
  }
  

  P0$table <- P0$table +
    theme(
      plot.title = element_text(size = base_font_size),
      axis.text  = element_text(size = base_font_size),
      axis.title = element_text(size = base_font_size),
      text       = element_text(size = base_font_size)
    ) +
    theme_classic(base_size = base_font_size * 0.8, base_line_size = 0.5) +
    theme(
      axis.line  = element_line(linewidth = axis_mm),
      axis.ticks = element_line(linewidth = axis_mm)
    )
  

  for (i in seq_along(P0$table$layers)) {
    if (inherits(P0$table$layers[[i]]$geom, "GeomText")) {
      P0$table$layers[[i]]$aes_params$size <- text_mm
    }
  }
  for (i in seq_along(P0$plot$layers)) {
    if (inherits(P0$plot$layers[[i]]$geom, "GeomPoint")) {
      P0$plot$layers[[i]]$aes_params$size <- censor_mm* (8/6)
    }
  }
  return(P0)
}
################Allcohort#############
csv_path <- KM_DATA_FILE
raw <- read.csv(
  csv_path,
  fileEncoding = "GB18030",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
km_required <- c("dataset", "cohort", "chemotherapy", "Multi_PRE_HRD", "OS.time", "OS",
                 "has_PFS", "PFS", "PFS.time", "stage_group")
km_missing <- setdiff(km_required, names(raw))
if (length(km_missing)) stop("KM source is missing fields: ", paste(km_missing, collapse = ", "))
if (!nrow(raw)) stop("KM source contains no observations.")
if (any(!is.na(raw$Multi_PRE_HRD) & !raw$Multi_PRE_HRD %in% 0:1))
  stop("Multi_PRE_HRD must contain only 0, 1 or NA.")
raw <- raw[raw$chemotherapy == 1 & !is.na(raw$chemotherapy), ]
####TCGA####
###TCGA-Train
Train<-raw[raw$dataset=='Train',]
Train$Pre_Label <- ifelse(Train$Multi_PRE_HRD ==1, "HRD+", "HRD-")
Train$Pre_Label <- factor(Train$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- Train
P0 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "TCGA-train cohort")
###TCGA-Val
Val<-raw[raw$dataset=='Val',]
Val$Pre_Label <- ifelse(Val$Multi_PRE_HRD ==1, "HRD+", "HRD-")
Val$Pre_Label <- factor(Val$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- Val
P1 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "TCGA-val cohort")
###TCGA-Internal
Internal<-raw[raw$dataset=='Internal',]
Internal$Pre_Label <- ifelse(Internal$Multi_PRE_HRD ==1, "HRD+", "HRD-")
Internal$Pre_Label <- factor(Internal$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- Internal
P2 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "TCGA-Internal cohort")
###TCGA-All
Train<-raw[raw$cohort=='TCGA',]
Train$Pre_Label <- ifelse(Train$Multi_PRE_HRD ==1, "HRD+", "HRD-")
Train$Pre_Label <- factor(Train$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- Train
P3 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "TCGA-all cohort")


####HMUCH####
HMUCH<-raw[raw$dataset=='HMUCH',]
HMUCH$Pre_Label <- ifelse(HMUCH$Multi_PRE_HRD ==1, "HRD+", "HRD-")
HMUCH$Pre_Label <- factor(HMUCH$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- HMUCH
P4 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "HMUCH cohort")


####CHCAMS####
CHCAMS<-raw[raw$dataset=='CHCAMS',]
CHCAMS$Pre_Label <- ifelse(CHCAMS$Multi_PRE_HRD ==1, "HRD+", "HRD-")
CHCAMS$Pre_Label <- factor(CHCAMS$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- CHCAMS
P5 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "CHCAMS cohort")

####PLCO####
PLCO<-raw[raw$dataset=='PLCO',]
PLCO$Pre_Label <- ifelse(PLCO$Multi_PRE_HRD ==1, "HRD+", "HRD-")
PLCO$Pre_Label <- factor(PLCO$Pre_Label, levels = c("HRD-", "HRD+"))
PD <- PLCO
P6 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "PLCO cohort")


splots <- list(P0, P1,P2,P3,P4,P5,P6)
res <- arrange_ggsurvplots(splots, print = TRUE, ncol = 5, nrow = 2)
ggsave(file.path(KM_OUT_DIR, "Overall_Survival_Kaplan_Meier.pdf"), plot = res,
       width = 24, height = 15, units = "cm")

##  Univariable + multivariable Cox forest plot (all covariates binary)
##
##  v3 - bug fixes after inspecting HRD_score_normalized_20260901.csv
##  ---------------------------------------------------------------------------
##  BUG 1  Residual disease vanished from every cohort.
##         The CSV stores residual_disease as "\u226410 mm" / ">10 mm" (already
##         binned, with a Unicode <= sign). The old case_when looked for
##         "No Macroscopic disease" / "1-10 mm" / "11-20 mm" / ">20 mm", so
##         EVERY row fell through to NA and the covariate was silently dropped
##         in all six cohorts. Fixed by normalising the string (\u2264 -> <=,
##         collapse whitespace, lower-case) and matching on patterns, and by
##         REPORTING any value that fails to map instead of silently NA-ing it.
##
##  BUG 2  PLCO lost almost everything.
##         Two independent reasons:
##           - residual_disease is 100% missing in PLCO (164/164) -> not
##             estimable, and that is a property of the data, not of the code;
##           - the 58-year cut is degenerate in PLCO: 160 subjects are >58 and
##             only 4 are <=58 (3 events), which fails MIN_PER_LEVEL = 5.
##             PLCO is a screening cohort, so almost nobody is under 58.
##         Both are now stated explicitly in the console and in a diagnostics
##         CSV rather than causing a row to disappear.
##
##  DESIGN CHANGE
##         Every cohort block now always prints the SAME four covariate rows.
##         A covariate that cannot be estimated shows "\u2014" and its reason is
##         logged, instead of the row being deleted. That is why variables
##         appeared to "go missing": var_usable() == FALSE removed the row
##         entirely, so a data problem looked like a plotting problem.
##         Set SHOW_ALL_ROWS <- FALSE to restore the old drop-the-row behaviour.
##
##  Also kept from v2: Grade removed; Age dichotomised at AGE_CUT = 58.

# Preserve analysis state across sections.

## ---- 0. Packages ------------------------------------------------------------
pkgs <- c("survival", "dplyr", "forestploter", "grid")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(pkgs, library, character.only = TRUE))


## ---- 1. User settings -------------------------------------------------------
data_file <- COX_DATA_FILE
out_dir <- COX_OUT_DIR

endpoint <- "OS"          # "OS" or "PFS"
if (!endpoint %in% c("OS", "PFS")) stop("endpoint must be OS or PFS.")

## ---- age dichotomisation ----------------------------------------------------
AGE_CUT      <- 58
AGE_LAB_LOW  <- sprintf("<=%g", AGE_CUT)   # reference,  Age <= 58
AGE_LAB_HIGH <- sprintf(">%g",  AGE_CUT)   # comparison, Age >  58

## ---- residual disease labels (match how the CSV is already binned) ----------
RES_LAB_LOW  <- "<=10 mm"                  # reference
RES_LAB_HIGH <- ">10 mm"

## ---- behaviour switches -----------------------------------------------------
SHOW_ALL_ROWS <- TRUE     # TRUE  = always print all four covariates per cohort,
#         with "-" where not estimable (recommended)
# FALSE = old behaviour, delete the row

## ---- cohorts ----------------------------------------------------------------
cohort_order <- c(
  "Train"    = "TCGA-train cohort",
  "Val"      = "TCGA-val cohort",
  "Internal" = "TCGA-internal cohort",
  "HMUCH"    = "HMUCH cohort",
  "CHCAMS"   = "CHCAMS cohort",
  "PLCO"     = "PLCO cohort"
)
## Dropped to stay consistent with the rest of Figure 3.
## Set to character(0) to put the internal split back in.
# DROP_COHORTS <- c("Internal")
DROP_COHORTS <- character(0)
if (length(DROP_COHORTS)) {
  hit <- intersect(DROP_COHORTS, names(cohort_order))
  if (length(hit))
    message("Cohorts excluded by DROP_COHORTS: ", paste(hit, collapse = ", "))
  cohort_order <- cohort_order[setdiff(names(cohort_order), DROP_COHORTS)]
}

## ---- reference levels (Grade removed) ---------------------------------------
REF <- list(
  PRISM    = "HRD-",          # other level: "HRD+"
  Age      = AGE_LAB_LOW,     # other level: AGE_LAB_HIGH
  Stage    = "I-II",          # other level: "III-IV"
  Residual = RES_LAB_LOW      # other level: RES_LAB_HIGH
)

## Minimum data requirements before a covariate is entered into a model.
## PLCO's age split is 4 vs 160; lowering MIN_PER_LEVEL below 5 would let that
## through, but the resulting HR would rest on 3 events and is not reportable.
MIN_PER_LEVEL <- 5    # subjects in each level of a factor
MIN_EVENTS    <- 5    # events in the analysable subset
MIN_N         <- 20   # analysable subjects

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


## ---- 2. Read data -----------------------------------------------------------
read_any <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  for (enc in c("GBK", "UTF-8-BOM", "UTF-8")) {
    x <- try(read.csv(path, fileEncoding = enc,
                      stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = c("", " ", "NA", "NaN", "NULL")),
             silent = TRUE)
    if (!inherits(x, "try-error") && ncol(x) > 1) {
      message("CSV encoding: ", enc); return(x)
    }
  }
  stop("Cannot parse the CSV with GBK / UTF-8-BOM / UTF-8.", call. = FALSE)
}
raw <- read_any(data_file)

time_col   <- if (endpoint == "OS") "OS.time" else "PFS.time"
status_col <- if (endpoint == "OS") "OS"      else "PFS"

required_cols <- c(
  "dataset", time_col, status_col, "Multi_PRE_HRD",
  "age", "stage_group", "residual_disease"
)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols))
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))


## ---- 3. Cleaning ------------------------------------------------------------
## Normalise free text before matching: this is what BUG 1 was about.
norm_txt <- function(x) {
  x <- as.character(x)
  x <- gsub("\u2264", "<=", x, fixed = TRUE)   # <=
  x <- gsub("\u2265", ">=", x, fixed = TRUE)   # >=
  x <- gsub("\u00A0", " ", x, fixed = TRUE)    # non-breaking space
  x <- gsub("[[:space:]]+", " ", x)
  tolower(trimws(x))
}
## loud reporting: any raw value that does not map is printed, not swallowed
report_unmapped <- function(raw_vals, mapped, var_name) {
  bad <- sort(unique(raw_vals[!is.na(raw_vals) & is.na(mapped)]))
  if (length(bad))
    message("!! ", var_name, ": these raw values did not match any rule and ",
            "became NA -> ", paste(sprintf("'%s'", bad), collapse = ", "))
  else
    message("   ", var_name, ": all non-missing values mapped successfully.")
}
make_fac <- function(x, ref, other) factor(x, levels = c(ref, other))

res_n   <- norm_txt(raw$residual_disease)
stage_n <- norm_txt(raw$stage_group)

work <- raw %>%
  dplyr::mutate(
    cohort_id = trimws(as.character(dataset)),
    time      = suppressWarnings(as.numeric(.data[[time_col]])),
    status    = suppressWarnings(as.numeric(.data[[status_col]])),
    age_num   = suppressWarnings(as.numeric(age)),
    
    PRISM = dplyr::case_when(
      Multi_PRE_HRD == 1 ~ "HRD+",
      Multi_PRE_HRD == 0 ~ "HRD-",
      TRUE               ~ NA_character_
    ),
    ## age dichotomised at AGE_CUT
    Age = dplyr::case_when(
      is.na(age_num)     ~ NA_character_,
      age_num <= AGE_CUT ~ AGE_LAB_LOW,
      age_num >  AGE_CUT ~ AGE_LAB_HIGH,
      TRUE               ~ NA_character_
    ),
    ## FIGO stage: accepts "I"/"II"/"III"/"IV" and "stage iii" style strings
    Stage = dplyr::case_when(
      grepl("^(stage )?(i|ii)$",   stage_n) ~ "I-II",
      grepl("^(stage )?(iii|iv)$", stage_n) ~ "III-IV",
      TRUE                                  ~ NA_character_
    ),
    ## Residual disease: covers BOTH the already-binned labels used in this CSV
    ## ("<=10 mm", ">10 mm") and the original TCGA categories.
    Residual = dplyr::case_when(
      grepl("no macroscopic|no gross|none|^0 ?mm$|^0$", res_n) ~ RES_LAB_LOW,
      grepl("<=ated?10", res_n)                                ~ RES_LAB_LOW,  # never matches; placeholder kept out
      grepl("<= ?10|0 ?- ?10|1 ?- ?10", res_n)                 ~ RES_LAB_LOW,
      grepl("> ?10|11 ?- ?20|> ?20|>= ?11", res_n)             ~ RES_LAB_HIGH,
      TRUE                                                     ~ NA_character_
    )
  )

cat("\n-- value-mapping check --------------------------------------------\n")
report_unmapped(raw$residual_disease, work$Residual, "residual_disease")
report_unmapped(raw$stage_group,      work$Stage,    "stage_group")
report_unmapped(raw$Multi_PRE_HRD,    work$PRISM,    "Multi_PRE_HRD")
cat("residual_disease raw -> mapped:\n")
print(table(raw = raw$residual_disease, mapped = work$Residual, useNA = "ifany"))

dat <- work %>%
  dplyr::select(cohort_id, time, status, age_num, PRISM, Age, Stage, Residual) %>%
  dplyr::filter(
    cohort_id %in% names(cohort_order),
    is.finite(time), time > 0,
    status %in% c(0, 1)
  ) %>%
  dplyr::mutate(
    PRISM    = make_fac(PRISM,    REF$PRISM,   "HRD+"),
    Age      = make_fac(Age,      AGE_LAB_LOW, AGE_LAB_HIGH),
    Stage    = make_fac(Stage,    REF$Stage,   "III-IV"),
    Residual = make_fac(Residual, RES_LAB_LOW, RES_LAB_HIGH)
  )

cohort_present <- names(cohort_order)[names(cohort_order) %in% unique(dat$cohort_id)]
if (!length(cohort_present))
  stop("No cohort survived cleaning. Values in `dataset`: ",
       paste(sort(unique(raw$dataset)), collapse = ", "))
if (length(setdiff(names(cohort_order), cohort_present)))
  message("Cohorts dropped (no usable ", endpoint, " data): ",
          paste(setdiff(names(cohort_order), cohort_present), collapse = ", "))

cohort_list <- lapply(cohort_present, function(x)
  droplevels(dat[dat$cohort_id == x, , drop = FALSE]))
names(cohort_list) <- unname(cohort_order[cohort_present])

cat("\n-- age split at ", AGE_CUT, " -------------------------------------\n", sep = "")
print(table(cohort = dat$cohort_id, age = dat$Age, useNA = "ifany"))
cat("\n-- residual disease by cohort -------------------------------------\n")
print(table(cohort = dat$cohort_id, residual = dat$Residual, useNA = "ifany"))


## ---- 4. Variable definitions ------------------------------------------------
## `term` must match the coefficient name produced by coxph, i.e.
## paste0(<variable>, <non-reference level>). All covariates are two-level
## factors now, so this holds for Age ("Age>58") and Residual ("Residual>10 mm").
var_defs <- list(
  list(var = "PRISM",    label = "PRISM",
       comparison = paste0("HRD+ vs. ", REF$PRISM),
       term = paste0("PRISM", "HRD+")),
  list(var = "Age",      label = "Age (years)",
       comparison = paste0(AGE_LAB_HIGH, " vs. ", AGE_LAB_LOW),
       term = paste0("Age", AGE_LAB_HIGH)),
  list(var = "Stage",    label = "FIGO stage",
       comparison = paste0("III-IV vs. ", REF$Stage),
       term = paste0("Stage", "III-IV")),
  list(var = "Residual", label = "Residual disease",
       comparison = paste0(RES_LAB_HIGH, " vs. ", RES_LAB_LOW),
       term = paste0("Residual", RES_LAB_HIGH))
)


## ---- 5. Cox helpers ---------------------------------------------------------
fmt_hr <- function(hr, lo, hi) {
  if (!all(is.finite(c(hr, lo, hi)))) return("\u2014")
  sprintf("%.2f (%.2f\u2013%.2f)", hr, lo, hi)
}
fmt_p <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
}
empty_result <- function()
  c(HR = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_)

## Returns ok = TRUE/FALSE AND the reason, so nothing fails silently.
var_check <- function(d, v) {
  if (!v %in% names(d))
    return(list(ok = FALSE, reason = "column absent", detail = ""))
  keep <- !is.na(d[[v]]) & !is.na(d$time) & !is.na(d$status)
  n <- sum(keep)
  if (n == 0)
    return(list(ok = FALSE, reason = "100% missing in this cohort", detail = ""))
  x  <- d[[v]][keep]
  ev <- d$status[keep]
  if (is.factor(x)) x <- droplevels(x) else x <- factor(x)
  tb  <- table(x)
  evb <- tapply(ev, x, function(z) sum(z == 1))
  detail <- paste(sprintf("%s: n=%d/ev=%d", names(tb), as.integer(tb),
                          as.integer(evb)), collapse = "; ")
  if (n < MIN_N)
    return(list(ok = FALSE,
                reason = sprintf("n = %d < MIN_N = %d", n, MIN_N), detail = detail))
  if (sum(ev == 1) < MIN_EVENTS)
    return(list(ok = FALSE,
                reason = sprintf("%d events < MIN_EVENTS = %d", sum(ev == 1),
                                 MIN_EVENTS), detail = detail))
  if (nlevels(x) < 2)
    return(list(ok = FALSE,
                reason = sprintf("only one level present ('%s')", levels(x)[1]),
                detail = detail))
  if (any(tb < MIN_PER_LEVEL)) {
    j <- which.min(tb)
    return(list(ok = FALSE,
                reason = sprintf("level '%s' has only %d subjects (< %d)",
                                 names(tb)[j], as.integer(tb)[j], MIN_PER_LEVEL),
                detail = detail))
  }
  if (any(evb < 1)) {
    j <- which.min(evb)
    return(list(ok = FALSE,
                reason = sprintf("level '%s' has 0 events", names(tb)[j]),
                detail = detail))
  }
  list(ok = TRUE, reason = "ok", detail = detail)
}

safe_cox_fit <- function(data, formula, model_label) {
  withCallingHandlers(
    tryCatch(
      survival::coxph(formula, data = data, ties = "efron",
                      na.action = na.omit,
                      control = survival::coxph.control(iter.max = 50)),
      error = function(e) {
        message("Cox error [", model_label, "]: ", conditionMessage(e)); NULL
      }
    ),
    warning = function(w) {
      message("Cox warning [", model_label, "]: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

extract_cox_row <- function(fit, term) {
  if (is.null(fit)) return(empty_result())
  s <- summary(fit)
  if (!(term %in% rownames(s$conf.int))) return(empty_result())
  c(HR    = unname(s$conf.int[term, "exp(coef)"]),
    lower = unname(s$conf.int[term, "lower .95"]),
    upper = unname(s$conf.int[term, "upper .95"]),
    p     = unname(s$coefficients[term, "Pr(>|z|)"]))
}


## ---- 6. Per-cohort model building -------------------------------------------
diag_log <- list()

build_cohort_block <- function(cohort_name, d) {
  
  ## (a) usability of every covariate, with reasons
  chk <- lapply(var_defs, function(v) var_check(d, v$var))
  names(chk) <- vapply(var_defs, function(v) v$var, character(1))
  avail <- vapply(chk, function(z) z$ok, logical(1))
  
  diag_log[[cohort_name]] <<- data.frame(
    cohort   = cohort_name,
    variable = vapply(var_defs, function(v) v$label, character(1)),
    usable   = avail,
    reason   = vapply(chk, function(z) z$reason, character(1)),
    levels   = vapply(chk, function(z) z$detail, character(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  cat("\n[", cohort_name, "]\n", sep = "")
  for (i in seq_along(var_defs))
    cat(sprintf("   %-18s %-4s %-42s %s\n", var_defs[[i]]$label,
                if (avail[i]) "OK" else "SKIP", chk[[i]]$reason, chk[[i]]$detail))
  
  if (!any(avail)) {
    message("Cohort skipped (no usable covariate): ", cohort_name)
    if (!SHOW_ALL_ROWS) return(NULL)
  }
  
  ## (b) multivariable set: drop covariates that collapse after complete-case
  mv_vars <- names(avail)[avail]
  if (length(mv_vars)) {
    repeat {
      dd  <- droplevels(d[stats::complete.cases(
        d[, c("time", "status", mv_vars), drop = FALSE]), , drop = FALSE])
      bad <- mv_vars[!vapply(mv_vars, function(v) var_check(dd, v)$ok, logical(1))]
      if (!length(bad) || length(mv_vars) <= 1) break
      mv_vars <- setdiff(mv_vars, bad)
    }
    dd <- droplevels(d[stats::complete.cases(
      d[, c("time", "status", mv_vars), drop = FALSE]), , drop = FALSE])
  } else dd <- d[0, , drop = FALSE]
  
  if (length(mv_vars) >= 2 && sum(dd$status == 1) < 5 * length(mv_vars))
    message("Note [", cohort_name, "]: only ", sum(dd$status == 1),
            " events for ", length(mv_vars),
            " covariates - multivariable estimates are unstable.")
  if (length(mv_vars) < 2)
    message("Note [", cohort_name, "]: fewer than two covariates available, ",
            "no multivariable model was fitted.")
  
  multi_fit <- if (length(mv_vars) >= 2) {
    safe_cox_fit(
      dd,
      stats::as.formula(paste("Surv(time, status) ~",
                              paste(sprintf("`%s`", mv_vars), collapse = " + "))),
      paste(cohort_name, "multivariable")
    )
  } else NULL
  
  n_tot <- nrow(d); n_ev <- sum(d$status == 1, na.rm = TRUE)
  
  ## (c) header row
  rows <- list(data.frame(
    is_header  = TRUE,
    Variable   = sprintf("%s (n = %d, events = %d)", cohort_name, n_tot, n_ev),
    Comparison = "",
    uni_HR = NA_real_, uni_lo = NA_real_, uni_hi = NA_real_,
    uni_text = "", uni_p_num = NA_real_, uni_p = "",
    multi_HR = NA_real_, multi_lo = NA_real_, multi_hi = NA_real_,
    multi_text = "", multi_p_num = NA_real_, multi_p = "",
    stringsAsFactors = FALSE
  ))
  
  ## (d) one row per covariate - ALL of them when SHOW_ALL_ROWS is TRUE
  keep_idx <- if (SHOW_ALL_ROWS) seq_along(var_defs) else which(avail)
  k <- 1L
  for (i in keep_idx) {
    v <- var_defs[[i]]
    if (avail[i]) {
      uni_fit <- safe_cox_fit(
        d,
        stats::as.formula(paste("Surv(time, status) ~", sprintf("`%s`", v$var))),
        paste(cohort_name, v$label, "univariable")
      )
      u <- extract_cox_row(uni_fit, v$term)
    } else u <- empty_result()
    m <- if (v$var %in% mv_vars) extract_cox_row(multi_fit, v$term) else empty_result()
    
    k <- k + 1L
    rows[[k]] <- data.frame(
      is_header  = FALSE,
      Variable   = v$label,
      Comparison = v$comparison,
      uni_HR = u["HR"], uni_lo = u["lower"], uni_hi = u["upper"],
      uni_text = fmt_hr(u["HR"], u["lower"], u["upper"]),
      uni_p_num = u["p"], uni_p = fmt_p(u["p"]),
      multi_HR = m["HR"], multi_lo = m["lower"], multi_hi = m["upper"],
      multi_text = fmt_hr(m["HR"], m["lower"], m["upper"]),
      multi_p_num = m["p"], multi_p = fmt_p(m["p"]),
      stringsAsFactors = FALSE
    )
  }
  out <- dplyr::bind_rows(rows)
  out$cohort <- cohort_name
  out
}

cat("\n-- covariate availability per cohort ------------------------------\n")
all_rows <- dplyr::bind_rows(
  lapply(names(cohort_list), function(nm) build_cohort_block(nm, cohort_list[[nm]]))
)
rownames(all_rows) <- NULL

if (!nrow(all_rows) || !any(is.finite(all_rows$uni_HR)))
  stop("No finite Cox estimates were produced. Review the messages above.")


## ---- 7. Data audit + numeric export -----------------------------------------
audit <- dplyr::bind_rows(lapply(names(cohort_list), function(nm) {
  d <- cohort_list[[nm]]
  data.frame(
    cohort      = nm,
    n           = nrow(d),
    events      = sum(d$status == 1, na.rm = TRUE),
    age_median  = round(stats::median(d$age_num, na.rm = TRUE), 1),
    age_min     = suppressWarnings(min(d$age_num, na.rm = TRUE)),
    age_max     = suppressWarnings(max(d$age_num, na.rm = TRUE)),
    n_age_low   = sum(d$Age == AGE_LAB_LOW,  na.rm = TRUE),
    n_age_high  = sum(d$Age == AGE_LAB_HIGH, na.rm = TRUE),
    n_res_low   = sum(d$Residual == RES_LAB_LOW,  na.rm = TRUE),
    n_res_high  = sum(d$Residual == RES_LAB_HIGH, na.rm = TRUE),
    PRISM_na    = sum(is.na(d$PRISM)),
    Age_na      = sum(is.na(d$Age)),
    Stage_na    = sum(is.na(d$Stage)),
    Residual_na = sum(is.na(d$Residual))
  )
}))
cat("\n-- per-cohort audit -----------------------------------------------\n")
print(audit, row.names = FALSE)

diag_tab <- dplyr::bind_rows(diag_log); rownames(diag_tab) <- NULL

result_table <- all_rows %>%
  dplyr::filter(!is_header) %>%
  dplyr::select(cohort, Variable, Comparison,
                uni_HR, uni_lo, uni_hi, uni_p_num,
                multi_HR, multi_lo, multi_hi, multi_p_num)

csv_file   <- file.path(out_dir, sprintf("PRISM_Cox_%s_results.csv", endpoint))
audit_file <- file.path(out_dir, sprintf("PRISM_Cox_%s_audit.csv", endpoint))
diag_file  <- file.path(out_dir, sprintf("PRISM_Cox_%s_diagnostics.csv", endpoint))
write.csv(result_table, csv_file,   row.names = FALSE, fileEncoding = "UTF-8")
write.csv(audit,        audit_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(diag_tab,     diag_file,  row.names = FALSE, fileEncoding = "UTF-8")


## ---- 8. Forest-plot table ---------------------------------------------------
plot_df <- data.frame(
  Variable     = all_rows$Variable,
  Comparison   = all_rows$Comparison,
  uni_forest   = strrep(" ", 24),
  uni_text     = all_rows$uni_text,
  uni_p        = all_rows$uni_p,
  multi_forest = strrep(" ", 24),
  multi_text   = all_rows$multi_text,
  multi_p      = all_rows$multi_p,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(plot_df) <- c(
  "Variable", "Comparison",
  "Univariable Cox analysis", "HR (95% CI)", "P value",
  "Multivariable Cox analysis", "HR (95% CI)\u00A0", "P value\u00A0"
)

cohort_idx <- which(all_rows$is_header)
data_idx   <- which(!all_rows$is_header)

## log-scale axis, clipped to a readable window
tick_pool <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
ci_lower  <- c(all_rows$uni_lo, all_rows$multi_lo)
ci_upper  <- c(all_rows$uni_hi, all_rows$multi_hi)
ci_lower  <- ci_lower[is.finite(ci_lower) & ci_lower > 0]
ci_upper  <- ci_upper[is.finite(ci_upper) & ci_upper > 0]

if (!length(ci_lower) || !length(ci_upper)) {
  xmin <- 0.2; xmax <- 5
} else {
  lo_obs <- max(min(ci_lower), 0.05)
  hi_obs <- min(max(ci_upper), 20)
  xmin <- max(tick_pool[tick_pool <= lo_obs], min(tick_pool))
  xmax <- min(tick_pool[tick_pool >= hi_obs], max(tick_pool))
  if (xmax <= xmin) { xmin <- 0.2; xmax <- 5 }
}
ticks <- tick_pool[tick_pool >= xmin & tick_pool <= xmax]
if (!(1 %in% ticks)) ticks <- sort(unique(c(ticks, 1)))


## ---- 9. Theme ---------------------------------------------------------------
NAVY <- "#1F3A5F"; GREY <- "#6B7280"; INK <- "#111827"

tm <- forest_theme(
  base_size = 8, base_family = "sans",
  ci_pch = 15, ci_col = NAVY, ci_fill = NAVY,
  ci_lwd = 1, ci_lty = 1, ci_Theight = 0.12,
  refline_gp  = gpar(lwd = 0.5, lty = "dashed", col = "#9CA3AF"),
  xaxis_gp    = gpar(col = INK, lwd = 0.5, fontsize = 8),
  arrow_type  = "closed",
  arrow_gp    = gpar(fontsize = 8, fontface = "italic", col = GREY, lwd = 0.5),
  footnote_gp = gpar(fontsize = 7.5, fontface = "italic", col = GREY),
  footnote_parse = FALSE,
  core = list(
    fg_params = list(fontface = "plain", col = INK, hjust = 0, x = 0.02),
    bg_params = list(fill = "white"),
    padding   = unit(c(0.6, 2.5), "mm")
  ),
  colhead = list(
    fg_params = list(fontface = "bold", col = INK, hjust = 0.5, x = 0.5),
    bg_params = list(fill = "white"),
    padding   = unit(c(1.8, 2.5), "mm")
  )
)


## ---- 10. Draw ---------------------------------------------------------------
foot <- paste0(
  "HR, hazard ratio; CI, confidence interval. P values are two-sided. ",
  "Endpoint: ", endpoint, ". Age was dichotomised at ", AGE_CUT, " years. ",
  "Reference groups: ",
  paste(sprintf("%s = %s", names(REF), unlist(REF)), collapse = ", "), ". ",
  "\u2014 indicates a covariate that was unavailable or not estimable in that ",
  "cohort; see the diagnostics file for the reason."
)

p <- forest(
  plot_df,
  est       = list(all_rows$uni_HR, all_rows$multi_HR),
  lower     = list(all_rows$uni_lo, all_rows$multi_lo),
  upper     = list(all_rows$uni_hi, all_rows$multi_hi),
  ci_column = c(3, 6),
  ref_line  = 1,
  x_trans   = "log",
  xlim      = c(xmin, xmax),
  ticks_at  = ticks,
  arrow_lab = c("Lower hazard", "Higher hazard"),
  theme     = tm,
  footnote  = foot
)

## cohort header rows
p <- edit_plot(p, row = cohort_idx, col = 1, gp = gpar(fontface = "bold", col = INK))
p <- edit_plot(p, row = cohort_idx, which = "background",
               gp = gpar(fill = "#F3F4F6", col = NA))

## body rows
p <- edit_plot(p, row = data_idx, col = 1, gp = gpar(fontface = "bold", col = INK))
p <- edit_plot(p, row = data_idx, col = 2, gp = gpar(col = GREY, fontface = "italic"))

## grey out rows that carry no estimate at all
blank_rows <- which(!all_rows$is_header &
                      !is.finite(all_rows$uni_HR) & !is.finite(all_rows$multi_HR))
if (length(blank_rows)) {
  p <- edit_plot(p, row = blank_rows, col = 1, gp = gpar(fontface = "plain", col = GREY))
  p <- edit_plot(p, row = blank_rows, col = 2, gp = gpar(col = GREY, fontface = "italic"))
}

## bold significant P values
sig_uni   <- which(is.finite(all_rows$uni_p_num)   & all_rows$uni_p_num   < 0.05)
sig_multi <- which(is.finite(all_rows$multi_p_num) & all_rows$multi_p_num < 0.05)
if (length(sig_uni))
  p <- edit_plot(p, row = sig_uni,   col = 5, gp = gpar(fontface = "bold", col = NAVY))
if (length(sig_multi))
  p <- edit_plot(p, row = sig_multi, col = 8, gp = gpar(fontface = "bold", col = NAVY))

## borders
p <- add_border(p, part = "header", row = 1, where = "top",    gp = gpar(lwd = 1.2, col = INK))
p <- add_border(p, part = "header", row = 1, where = "bottom", gp = gpar(lwd = 0.8, col = INK))
p <- add_border(p, part = "body", row = nrow(plot_df), where = "bottom",
                gp = gpar(lwd = 1, col = INK))
if (length(cohort_idx) > 1) {
  for (r in cohort_idx[-1])
    p <- add_border(p, part = "body", row = r, where = "top",
                    gp = gpar(lwd = 0.3, col = "#D1D5DB"))
}


## ---- 11. Export -------------------------------------------------------------
win_safe <- function(x) {
  z <- try(iconv(x, from = "UTF-8", to = "CP1252"), silent = TRUE)
  if (inherits(z, "try-error")) return(FALSE)
  !any(is.na(z))
}
all_text  <- c(unlist(plot_df, use.names = FALSE), colnames(plot_df), foot)
use_cairo <- !win_safe(all_text) && isTRUE(capabilities("cairo"))
if (!win_safe(all_text) && !use_cairo)
  warning("Non-WinAnsi glyphs present and cairo is unavailable; ",
          "some characters may not render.")

fig_w_in <- 180 / 25.4 * 1.15
fig_h_in <- 0.17 * nrow(all_rows) + 1.5
pdf_file <- file.path(out_dir, sprintf("PRISM_forest_%s.pdf", endpoint))

if (use_cairo) {
  grDevices::cairo_pdf(pdf_file, width = fig_w_in, height = fig_h_in, onefile = TRUE)
} else {
  grDevices::pdf(pdf_file, width = fig_w_in, height = fig_h_in,
                 useDingbats = FALSE, onefile = TRUE)
}
plot(p)
dev.off()

plot(p)
cat("\nSaved:\n", pdf_file, "\n", csv_file, "\n", audit_file, "\n", diag_file, "\n",
    sep = "")


##  Two-year progression status by MIRROR-predicted HRD group,
##  overall and within FIGO stage III-IV (cohort-specific endpoints PFI/DFS/RFS)
##


library(dplyr)
library(ggplot2)
library(survival)


CSV <- KM_DATA_FILE
OUT_PDF <- file.path(KM_OUT_DIR, "Two_Year_Progression_Stacked.pdf")
GROUP_COL <- "Multi_PRE_HRD"
LANDMARK  <- 2
P_METHOD  <- "logrank"
if (!P_METHOD %in% c("logrank", "fisher"))
  stop("P_METHOD must be logrank or fisher.")

W <- 164.2635; H <- 347.9281
FS   <- 8
FONT <- if (.Platform$OS.type == "windows") "Arial" else "Helvetica"


X_L    <- 1.5
X0     <- 27
BARW   <- 112
X_N    <- X0 + BARW + 3
X_R    <- W - 2
BARH   <- 8.5
PITCH  <- 11
H_HEAD <- 12
H_SUB  <- 11
G_STR  <- 6
G_COH  <- 7
TOP    <- 3
D_AXIS <- 2.0
D_TICK <- 2.5
D_TLAB <- 9.0
D_XTTL <- 19.5
D_LEG  <- 30.0
MIN_W  <- 13

COL_EVENT <- "#C0703A"; COL_FREE <- "#4A7CA5"; COL_CENS <- "#DDDDDD"
INK <- "grey15"; SUBINK <- "grey30"; FAINT <- "grey50"

sz <- FS / .pt                 # pt -> ggplot size


progression_required <- c("has_PFS", "PFS", "PFS.time", "stage_group", GROUP_COL,
                          "chemotherapy", "cohort")
progression_missing <- setdiff(progression_required, names(raw))
if (length(progression_missing))
  stop("Landmark source is missing fields: ", paste(progression_missing, collapse = ", "))
dat <- read.csv(CSV, stringsAsFactors = FALSE, fileEncoding = "GBK") %>%
  filter(has_PFS == 1, !is.na(PFS), !is.na(PFS.time)) %>%
  dplyr::mutate(
    grp = factor(ifelse(.data[[GROUP_COL]] == 1, "HRD+", "HRD\u2212"),
                 levels = c("HRD\u2212", "HRD+")),
    status2 = factor(dplyr::case_when(
      PFS == 1 & PFS.time <= LANDMARK ~ "Event",
      PFS.time >= LANDMARK            ~ "Free",
      TRUE                            ~ "Cens"),
      levels = c("Event", "Free", "Cens")),
    stage_bin = ifelse(stage_group %in% c("III", "IV"), "III\u2013IV", "I\u2013II"))

dat <- dat[dat$chemotherapy == 1 & !is.na(dat$chemotherapy), ]


cohorts <- data.frame(cohort   = c("TCGA", "CHCAMS", "HMUCH"),
                      endpoint = c("PFI",  "DFS",    "RFS"),
                      stringsAsFactors = FALSE)
strata  <- data.frame(key   = c("All", "III\u2013IV"),
                      label = c("All patients", "FIGO stage III\u2013IV"),
                      stringsAsFactors = FALSE)

pfun <- function(d) {
  if (dplyr::n_distinct(d$grp) < 2) return(NA_real_)
  if (P_METHOD == "logrank") {
    s <- survival::survdiff(Surv(PFS.time, PFS) ~ grp, data = d)
    1 - pchisq(s$chisq, length(s$n) - 1)
  } else {
    fisher.test(table(d$grp, factor(d$status2, levels = c("Event", "Free"))))$p.value
  }
}
plab <- function(p) if (is.na(p)) "" else
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)


R <- list(); T <- list(); ir <- 0; it <- 0
addR <- function(xmin, xmax, ymin, ymax, fill) {
  ir <<- ir + 1
  R[[ir]] <<- data.frame(xmin, xmax, ymin, ymax, fill, stringsAsFactors = FALSE)
}
addT <- function(x, y, lab, hjust, col, face = "plain") {
  it <<- it + 1
  T[[it]] <<- data.frame(x, y, lab, hjust, col, face, stringsAsFactors = FALSE)
}

y <- H - TOP
for (k in seq_len(nrow(cohorts))) {
  co <- cohorts$cohort[k]; ep <- cohorts$endpoint[k]
  d_co <- filter(dat, cohort == co)
  

  y <- y - H_HEAD
  addT(X_L, y, co, 0, INK, "bold")
  addT(X_R, y, sprintf("%s \u00B7 n = %d", ep, nrow(d_co)), 1, FAINT)
  
  for (s in seq_len(nrow(strata))) {
    d_st <- if (strata$key[s] == "All") d_co else filter(d_co, stage_bin == "III\u2013IV")
    

    y <- y - H_SUB
    addT(X_L, y, strata$label[s], 0, SUBINK)
    addT(X_R, y, plab(pfun(d_st)), 1, SUBINK)
    

    for (g in levels(dat$grp)) {
      dg <- filter(d_st, grp == g); n_g <- nrow(dg)
      if (n_g == 0) next
      y <- y - PITCH
      pr <- as.numeric(table(dg$status2)) / n_g * 100
      names(pr) <- levels(dat$status2)
      xs <- X0 + c(0, cumsum(pr)) / 100 * BARW
      
      for (m in seq_along(pr)) {
        if (pr[m] <= 0) next
        addR(xs[m], xs[m + 1], y - BARH / 2, y + BARH / 2, names(pr)[m])
        if (xs[m + 1] - xs[m] >= MIN_W)
          addT((xs[m] + xs[m + 1]) / 2, y, sprintf("%.0f", pr[m]), 0.5,
               if (names(pr)[m] == "Cens") "grey25" else "white")
      }
      addT(X0 - 3, y, g,               1, INK)
      addT(X_N,    y, as.character(n_g), 0, FAINT)
    }
    y <- y - G_STR
  }
  y <- y - G_COH
}


y_axis <- y + G_COH - D_AXIS
tick_x <- X0 + c(0, 50, 100) / 100 * BARW

addT(tick_x[1], y_axis - D_TLAB, "0",   0.5, INK)
addT(tick_x[2], y_axis - D_TLAB, "50",  0.5, INK)
addT(tick_x[3], y_axis - D_TLAB, "100", 0.5, INK)
addT(X0 + BARW / 2, y_axis - D_XTTL, "Patients (%)", 0.5, INK)

leg_txt <- c(sprintf("Progression/recurrence \u2264%d y", LANDMARK),
             sprintf("Progression-free at %d y", LANDMARK),
             sprintf("Censored <%d y", LANDMARK))
leg_fil <- c("Event", "Free", "Cens")
for (i in 1:3) {
  yy <- y_axis - D_LEG - (i - 1) * PITCH
  addR(X_L, X_L + 7, yy - 3.5, yy + 3.5, leg_fil[i])
  addT(X_L + 10, yy, leg_txt[i], 0, INK)
}

rect_df <- bind_rows(R)
text_df <- bind_rows(T)


bottom <- min(rect_df$ymin, text_df$y - FS / 2)
message(sprintf("bottom margin = %.2f pt (must be positive)", bottom))
stopifnot(bottom > 0)


p_G <- ggplot() +
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  annotate("segment", x = X0, xend = X0 + BARW, y = y_axis, yend = y_axis,
           linewidth = 0.3, colour = INK) +
  annotate("segment", x = tick_x, xend = tick_x,
           y = y_axis, yend = y_axis - D_TICK, linewidth = 0.3, colour = INK) +
  geom_text(data = text_df,
            aes(x = x, y = y, label = lab, hjust = hjust,
                colour = col, fontface = face),
            size = sz, family = FONT, vjust = 0.5) +
  scale_fill_manual(values = c(Event = COL_EVENT, Free = COL_FREE, Cens = COL_CENS)) +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") +
  theme_void(base_family = FONT) +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, unit = "pt"))


ggsave(OUT_PDF, p_G, width = W / 72, height = H / 72, units = "in",
       device = cairo_pdf)


## ggsave(sub("\\.pdf$", ".png", OUT_PDF), p_G,
##        width = W/72, height = H/72, units = "in", dpi = 600)
