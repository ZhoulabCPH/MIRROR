# =============================================================================
# Figures 4-5 | Survival in unlabeled and BRCA wild-type subsets
# =============================================================================
# Scientific scope: Kaplan-Meier curves and univariable/multivariable Cox estimates.
# Usage: Rscript 05_Unlabeled_and_BRCA_Wildtype_Survival.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: BRCA-annotated clinical CSV; wild-type source value is represented as Unicode.
# Outputs: Survival PDF, Cox CSV and Word table in MIRROR_OUTPUT_DIR.
# Dependencies: dplyr, survival, survminer, ggplot2, ggpubr;
# flextable and officer are optional for Word export.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# Historical dimensions are audited rather than imposed on updated data.
# =============================================================================

required_packages <- c("dplyr", "survival", "survminer", "ggplot2", "ggpubr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(invisible(lapply(required_packages, library, character.only = TRUE)))

KaplanMeier_Plot <- function(PD, Y_labels, titles, type = "l") {
  valid <- is.finite(suppressWarnings(as.numeric(PD$OS.time))) &
    PD$OS %in% c(0, 1) & !is.na(PD$Pre_Label)
  PD <- PD[valid, , drop = FALSE]
  if (nrow(PD) < 5L || nlevels(droplevels(PD$Pre_Label)) < 2L ||
      sum(PD$OS == 1) < 1L)
    stop("Cannot estimate survival for ", titles,
         ": at least five valid observations, two groups and an event are required.")

  height_ratio <- 0.75
  base_font_size <- 8
  line_size <- 0.5
  table_height_ratio<-0.15
  
  
  fit <- survfit(Surv(OS.time, OS) ~ Pre_Label, data = PD)
  P0 <- ggsurvplot(
    fit = fit,
    title = titles,
    risk.table = TRUE,
    ggtheme = theme_classic(base_size = base_font_size),
    risk.table.col = "strata",
    ##76A992  #666292
    palette = c("#416993", "#AE6534"),
    pval = TRUE,
    pval.method = TRUE,
    ylab = Y_labels, 
    xlab = "Time (years)",
    # tables.theme = theme_cleantable(base_size = base_font_size),
    tables.height = 0.25,
    font.x = c(base_font_size, "plain", "black"),
    font.y = c(base_font_size, "plain", "black"),
    font.tickslab = c(base_font_size, "plain", "black"),
    ylim = c(0, 1),
    size = 1.5,
    censor.size = base_font_size*0.8,
    legend.labs = c("High risk", "Low risk"),
    risk.table.title="Number at risk",
    conf.int = TRUE
  )

  P0$plot <- P0$plot + theme(panel.grid = element_blank())
  

  res_cox <- tryCatch(coxph(Surv(OS.time, OS) ~ Pre_Label, data = PD),
                      error = function(e) NULL)
  if (!is.null(res_cox)) {
    ci <- summary(res_cox)$conf.int
    if (nrow(ci) && all(is.finite(ci[1, c(1, 3, 4)])))
      P0$plot <- P0$plot + annotate("text", x = 4.5, y = 0.12,
        label = sprintf("HR = %.3f (%.3f-%.3f)", ci[1, 1], ci[1, 3], ci[1, 4]))
  }

  P0$table <- P0$table + theme(
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8)
  ) + theme_classic(base_size = base_font_size*0.8, base_line_size = 0.5)
  
  return(P0)
}

DATA_FILE <- Sys.getenv("MIRROR_BRCA_CLINICAL_CSV", unset = "<PATH_TO_BRCA_CLINICAL_CSV>")
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
assert_input_file(DATA_FILE, "DATA_FILE")
prepare_output_dir(OUT_DIR, "OUT_DIR")
x <- read.csv(
  DATA_FILE,
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", " ", "NA", "NaN", "NULL")
)

required_columns <- c("sample_id", "cohort", "dataset", "chemotherapy", "OS.time", "OS",
                      "Multi_PRE_HRD", "BRCA_somatic_status", "age", "stage_group",
                      "residual_disease")
missing_columns <- setdiff(required_columns, names(x))
if (length(missing_columns)) stop("Missing clinical fields: ", paste(missing_columns, collapse = ", "))
if (!nrow(x)) stop("The BRCA clinical CSV has no observations.")
if (any(!is.na(x$Multi_PRE_HRD) & !x$Multi_PRE_HRD %in% 0:1))
  stop("Multi_PRE_HRD must contain only 0, 1 or NA.")
if (nrow(x) != 937L || ncol(x) != 48L)
  warning("Clinical dimensions differ from the historical export (937 x 48): ", nrow(x), " x ", ncol(x))

table(
  x$BRCA_somatic_status[x$cohort == "TCGA"],
  useNA = "ifany"
)

TCGA<-x[x$cohort=='TCGA',]
datasets<-TCGA[TCGA$chemotherapy %in% c(1),]
Discovery<-datasets[datasets$dataset %in% c("Train", "Val"),]
OtherTest<-datasets[datasets$dataset=='Internal',]
#HMUCH
HMUCH<-x[x$cohort=='HMUCH',]
HMUCH<-HMUCH[HMUCH$chemotherapy %in% c(1),]
#CHCAMS
CHCAMS<-x[x$cohort=='CHCAMS',]
CHCAMS<-CHCAMS[CHCAMS$chemotherapy %in% c(1),]

PD<-Discovery
PD$HR <- ifelse(PD$Multi_PRE_HRD ==0, "High_Risk", "Low_Risk")
PD$Pre_Label <- factor(PD$HR, levels = c("High_Risk", "Low_Risk"))
P0 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "Reference-labeled ovarian tumors")


PD<-OtherTest
PD$HR <- ifelse(PD$Multi_PRE_HRD ==0, "High_Risk", "Low_Risk")
PD$Pre_Label <- factor(PD$HR, levels = c("High_Risk", "Low_Risk"))
P1 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "HRD status unavailable cohort")


PD<-x[x$cohort=='TCGA',]
PD<-PD[PD$chemotherapy %in% c(1),]
PD=PD[PD$BRCA_somatic_status %in% c("\u91ce\u751f\u578b"),]
PD$HR <- ifelse(PD$Multi_PRE_HRD ==0, "High_Risk", "Low_Risk")
PD$Pre_Label <- factor(PD$HR, levels = c("High_Risk", "Low_Risk"))
P2 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "BRCA1/2 somatic wild-type")

PD<-HMUCH
PD$HR <- ifelse(PD$Multi_PRE_HRD ==0, "High_Risk", "Low_Risk")
PD$Pre_Label <- factor(PD$HR, levels = c("High_Risk", "Low_Risk"))
P3 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "HMUCH cohort")

PD<-CHCAMS
PD$HR <- ifelse(PD$Multi_PRE_HRD ==0, "High_Risk", "Low_Risk")
PD$Pre_Label <- factor(PD$HR, levels = c("High_Risk", "Low_Risk"))
P4 <- KaplanMeier_Plot(PD, Y_labels = "Overall Survival", titles = "CHCAMS cohort")

splots <- list(P0, P1,P2, P3,P4)
res <- arrange_ggsurvplots(splots, print = TRUE, ncol = 5, nrow = 2)
ggsave(file.path(OUT_DIR, "Unlabeled_and_BRCA_Wildtype_Survival.pdf"),
       plot = res, width = 24, height = 15, units = "cm")




AGE_CUT <- 58


BRCAwt <- x[x$cohort == 'TCGA', ]
BRCAwt <- BRCAwt[BRCAwt$chemotherapy %in% c(1), ]
BRCAwt <- BRCAwt[BRCAwt$BRCA_somatic_status %in% c("\u91ce\u751f\u578b"), ]
cat("BRCA1/2 somatic wild-type cohort size: ", nrow(BRCAwt), "\n")

cox_groups <- list(
  "Reference-labeled ovarian tumors" = Discovery,
  "HRD status unavailable cohort" = OtherTest,
  "BRCA1/2 somatic wild-type" = BRCAwt,
  "HMUCH cohort" = HMUCH,
  "CHCAMS cohort" = CHCAMS
)

required <- c("OS.time", "OS", "Multi_PRE_HRD", "age", "stage_group", "residual_disease")
for (group_name in names(cox_groups)) {
  missing <- setdiff(required, names(cox_groups[[group_name]]))
  if (length(missing)) stop(group_name, " is missing fields: ", paste(missing, collapse = ", "))
}


make_cox_data <- function(z) {
  age <- suppressWarnings(as.numeric(as.character(z$age)))
  rd <- trimws(as.character(z$residual_disease))
  rd_big <- !is.na(rd) &
    (grepl(">\\s*(10|20)", rd) | grepl("^11\\s*-\\s*20", rd))
  rd_small <- !is.na(rd) & !rd_big &
    (grepl("10\\s*mm", rd) | grepl("no macroscopic", rd, ignore.case = TRUE))
  data.frame(
    time = suppressWarnings(as.numeric(as.character(z$OS.time))),
    status = suppressWarnings(as.numeric(as.character(z$OS))),
    Risk = factor(ifelse(z$Multi_PRE_HRD == 1, "Low risk",
                         ifelse(z$Multi_PRE_HRD == 0, "High risk", NA)),
                  levels = c("High risk", "Low risk")),
    Age = factor(ifelse(is.na(age), NA,
                        ifelse(age <= AGE_CUT, "<=58", ">58")),
                 levels = c("<=58", ">58")),
    Stage = factor(ifelse(is.na(z$stage_group), NA,
                          ifelse(z$stage_group %in% c("III", "IV"), "III-IV", "I-II")),
                   levels = c("I-II", "III-IV")),
    Residual = factor(ifelse(rd_big, ">10 mm",
                             ifelse(rd_small, "<=10 mm", NA)),
                      levels = c("<=10 mm", ">10 mm"))
  )
}

vars <- c("Risk", "Age", "Stage", "Residual")
labels <- c("MIRROR classification", "Age (years)", "FIGO stage", "Residual disease")
comparisons <- c("Low risk vs. High risk", ">58 vs. <=58",
                 "III-IV vs. I-II", ">10 mm vs. <=10 mm")


usable <- function(dd, v) {
  x <- dd[!is.na(dd[[v]]), , drop = FALSE]
  f <- droplevels(x[[v]])
  nrow(x) >= 20 && sum(x$status == 1) >= 5 && nlevels(f) == 2 &&
    all(table(f) >= 5) && all(tapply(x$status == 1, f, sum) >= 1)
}
fit_cox <- function(dd, vs) {
  if (!length(vs)) return(NULL)
  cc <- dd[complete.cases(dd[, vs, drop = FALSE]), , drop = FALSE]
  if (length(vs) > 1 &&
      (nrow(cc) < 20 || sum(cc$status == 1) < 5 ||
       any(!vapply(vs, function(v) usable(cc, v), logical(1))))) return(NULL)
  suppressWarnings(tryCatch(
    survival::coxph(stats::reformulate(vs, response = "survival::Surv(time, status)"),
                    data = cc, ties = "efron"), error = function(e) NULL))
}
get_result <- function(fit, v) {
  if (is.null(fit)) return(rep(NA_real_, 4))
  s <- summary(fit)
  j <- which(startsWith(rownames(s$conf.int), v))
  if (length(j) != 1) return(rep(NA_real_, 4))
  c(s$conf.int[j, "exp(coef)"], s$conf.int[j, "lower .95"],
    s$conf.int[j, "upper .95"], s$coefficients[j, "Pr(>|z|)"])
}
fmt_hr <- function(x) {
  if (!all(is.finite(x[1:3]))) return("—")
  sprintf("%.3f (%.3f–%.3f)", x[1], x[2], x[3])
}
fmt_p <- function(p) {
  if (!is.finite(p)) return("—")
  if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
}

cox_rows <- list()
header_rows <- integer(0)
for (group_name in names(cox_groups)) {
  z <- cox_groups[[group_name]]
  dd <- make_cox_data(z)
  dd <- dd[is.finite(dd$time) & dd$time > 0 & dd$status %in% c(0, 1), , drop = FALSE]
  header_rows <- c(header_rows, length(cox_rows) + 1L)
  cox_rows[[length(cox_rows) + 1L]] <- c(
    sprintf("%s (n = %d; %d deaths)", group_name, nrow(z), sum(dd$status == 1)),
    "", "", "", "", "")
  ok <- vapply(vars, function(v) usable(dd, v), logical(1))
  mv_vars <- vars[ok]
  mv_fit <- if (length(mv_vars) >= 2) fit_cox(dd, mv_vars) else NULL
  for (j in seq_along(vars)) {
    uni <- get_result(if (ok[j]) fit_cox(dd, vars[j]) else NULL, vars[j])
    multi <- get_result(if (vars[j] %in% mv_vars) mv_fit else NULL, vars[j])
    cox_rows[[length(cox_rows) + 1L]] <- c(
      labels[j], comparisons[j], fmt_hr(uni), fmt_p(uni[4]),
      fmt_hr(multi), fmt_p(multi[4]))
  }
}
cox_table <- as.data.frame(do.call(rbind, cox_rows), check.names = FALSE)
names(cox_table) <- c("Variable", "Comparison", "Univariable HR (95% CI)",
                      "P value", "Multivariable HR (95% CI)", "P value ")
write.csv(cox_table, file.path(OUT_DIR, "Table2_Cox_OS_five_cohorts.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
ft <- flextable::flextable(cox_table)
ft <- flextable::theme_booktabs(ft)
ft <- flextable::bold(ft, i = header_rows, j = 1)
ft <- flextable::align(ft, j = 2:6, align = "center", part = "all")
ft <- flextable::fontsize(ft, size = 8, part = "all")
ft <- flextable::font(ft, fontname = "Arial", part = "all")
ft <- flextable::width(ft, j = 1:6, width = c(1.25, 1.27, 1.36, 0.55, 1.36, 0.55))
ft <- flextable::set_table_properties(ft, layout = "fixed")
ft <- flextable::add_footer_lines(ft, paste(
  "HR, hazard ratio; CI, confidence interval. Cohort n is unchanged from the",
  "survival curves. Cox models exclude missing covariates as needed; multivariable",
  "models use complete cases. Age cut-off: 58 years. —, unavailable or not estimable."
))
doc <- officer::read_docx()
doc <- officer::body_add_par(doc,
                             "Table 2. Univariable and multivariable Cox regression for overall survival",
                             style = "heading 2")
doc <- flextable::body_add_flextable(doc, ft)
print(doc, target = file.path(OUT_DIR, "Table2_Cox_OS_five_cohorts_three_line.docx"))
} else {
  warning("Word export skipped: install flextable and officer; the Cox CSV remains available.")
}
cat("Cox analysis saved to: ", OUT_DIR, "\n")
