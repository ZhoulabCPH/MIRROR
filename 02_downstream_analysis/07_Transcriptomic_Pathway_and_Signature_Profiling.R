# =============================================================================
# Figure 6 | Transcriptomic pathway and signature profiling
# =============================================================================
# Scientific scope: Differential expression, enrichment, ssGSEA and all-RNA cohort audits.
# Usage: Rscript 07_Transcriptomic_Pathway_and_Signature_Profiling.R
# Configuration: Replace the <PATH_TO_...> placeholders below or set the
# corresponding MIRROR_* environment variables before running.
# Inputs: RMA and z-score expression CSVs, clinical CSV, Hallmark GMT; Fges GMT optional.
# Outputs: A4 panels, source tables, diagnostic report and session record in All_RNA_A4_8pt/.
# Dependencies: ggplot2, dplyr, tibble, patchwork; analyse also requires Bioconductor packages.
# Reproducibility: Cohort definitions and statistical thresholds follow the
# source analyses; input data are never modified.
# Select RUN_MODE = "analyse" for a full analysis or "replot" for cached results.
# F_FGES may be empty to invoke the documented provisional marker panel.
# =============================================================================

run_mirror_transcriptomics <- function() {
RUN_MODE <- "analyse"
if (!RUN_MODE %in% c("replot", "analyse")) stop("RUN_MODE must be 'replot' or 'analyse'.")
REPLOT_ADJUSTED <- FALSE
SAVE_A4_CACHE <- getOption("mirror.save_cache", FALSE)

BASE_RESULT_DIR <- Sys.getenv("MIRROR_OUTPUT_DIR", unset = "<PATH_TO_OUTPUT_DIRECTORY>")
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


RESULT_DIR <- file.path(BASE_RESULT_DIR, "All_RNA_A4_8pt")
WORKDIR <- RESULT_DIR

F_EXPR     <- Sys.getenv("MIRROR_RMA_CSV", unset = "<PATH_TO_RMA_EXPRESSION_CSV>")
F_EXPR_Z   <- Sys.getenv("MIRROR_ZSCORE_CSV", unset = "<PATH_TO_ZSCORE_EXPRESSION_CSV>")
F_CLIN     <- Sys.getenv("MIRROR_CLINICAL_CSV", unset = "<PATH_TO_CLINICAL_CSV>")
F_HALLMARK <- Sys.getenv("MIRROR_HALLMARK_GMT", unset = "<PATH_TO_HALLMARK_GMT>")

F_FGES <- Sys.getenv("MIRROR_FGES_GMT", unset = "<OPTIONAL_PATH_TO_FGES_GMT_OR_EMPTY>")

## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
HARMONIZE_SYMBOLS <- TRUE
FGES_MIN_COVERAGE <- 0.60
FGES_MIN_GENES    <- 4
FGES_LABEL_COL    <- 2

DDR_MIN_COVERAGE  <- 0.60
DDR_MIN_GENES     <- 4

## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
COHORT_SETS <- list(
  All_RNA  = c("Train", "Val", "Internal"),
  All      = c("Train", "Val", "Internal"),
  Train    = "Train",
  Val      = "Val",
  Internal = "Internal",
  TrainVal = c("Train", "Val")
)

RUN_SETS         <- c("TrainVal")
RUN_ADJUSTED_FOR <- c("TrainVal")
ADJ_DROP_E       <- TRUE
COMPARE_COHORTS  <- TRUE

RUN_TRUTH       <- FALSE
RUN_CONCORDANCE <- FALSE

## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------

LABEL_COL <- "HRD_label"


FDR_CUT <- 0.05
LFC_CUT <- 0.15

MAX_HEATMAP_ROWS  <- 800
N_PATH_EACH       <- 10
ORDER_D           <- "nes"
DROP_DISEASE_KEGG <- TRUE

FDR_FLOOR <- 1e-300


FIG_W_MM    <- 297
FIG_H_MM    <- 210
BASE_SIZE   <- 8
HEAT_LIM    <- 2
SAVE_PANELS <- FALSE

set.seed(20240601)
prepare_output_dir(WORKDIR, "WORKDIR")
assert_input_file(F_CLIN, "F_CLIN")
assert_input_file(F_EXPR_Z, "F_EXPR_Z")
if (RUN_MODE == "analyse") {
  assert_input_file(F_EXPR, "F_EXPR")
  assert_input_file(F_HALLMARK, "F_HALLMARK")
}
if (identical(F_FGES, "<OPTIONAL_PATH_TO_FGES_GMT_OR_EMPTY>")) F_FGES <- ""
if (nzchar(F_FGES)) assert_input_file(F_FGES, "F_FGES")
previous_directory <- getwd()
on.exit(setwd(previous_directory), add = TRUE)
setwd(WORKDIR)

select_all_rna_samples <- function(dat) {
  required <- c("sample_id", "mask_SeqRNA", "dataset")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("All-RNA selection is missing fields: ", paste(missing, collapse=", "))
  dat <- dat[dat$mask_SeqRNA %in% 1, , drop=FALSE]
  if (!nrow(dat)) stop("The clinical table contains no eligible RNA samples.")
  dat
}


need <- c("ggplot2", "dplyr", "tibble", "patchwork")
if (RUN_MODE == "analyse") need <- c(need, "org.Hs.eg.db", "AnnotationDbi", "GSVA", "fgsea")
miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
if (length(miss))
  stop("Missing required packages: ", paste(miss, collapse = ", "),
       "\n  CRAN: install.packages(c('ggplot2','dplyr','tibble','patchwork'))",
       "\n  Bioconductor: BiocManager::install(c('org.Hs.eg.db','GSVA','fgsea'))")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tibble)
  library(patchwork)
  if (RUN_MODE == "analyse") { library(org.Hs.eg.db); library(GSVA) }
})

HAS_MT      <- requireNamespace("matrixTests", quietly = TRUE)
HAS_MSIGDBR <- requireNamespace("msigdbr",     quietly = TRUE)
HAS_NEWSC   <- requireNamespace("ggnewscale",  quietly = TRUE)
message("Optional packages: matrixTests = ", HAS_MT, ", msigdbr = ", HAS_MSIGDBR,
        ", ggnewscale = ", HAS_NEWSC)

LOG <- file.path(WORKDIR, "A4_diagnostics_report.txt")

logmsg <- function(...) {
  txt <- paste0(...)
  message(txt)
  cat(txt, "\n", file = LOG, append = TRUE, sep = "")
}
logtable <- function(x, row.names = FALSE) {

  txt <- if (is.data.frame(x)) capture.output(print(x, row.names = row.names))
  else                  capture.output(print(x))
  cat(paste(txt, collapse = "\n"), "\n", file = LOG, append = TRUE, sep = "")
  cat(paste(txt, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
cat("MIRROR Figure 6 A4 / 8 pt diagnostics — ", format(Sys.time()), "\n",
    file = LOG, sep = "")



COL_NEG   <- "#3C6489"
COL_POS   <- "#A36134"
COL_NEG_L <- "#C9D5E0"
COL_POS_L <- "#E8D3C4"
COL_NS    <- "#D9D9D9"
COL_TEXT  <- "#121727"
COL_PANEL <- "#F0F1F3"

MINUS   <- "\u2212"
LAB_NEG <- paste0("HRD", MINUS)
LAB_POS <- "HRD+"
UP_ARR  <- "\u2191"
DN_ARR  <- "\u2193"

HEAT_COL <- colorRampPalette(c("#9C4DCC", "#FFFFFF", "#FFD702"))(255)

CAT_LEVELS <- c("Immune", "Metabolism", "Signaling", "Proliferation",
                "DNA damage", "Development", "Cellular component")
CAT_COL <- c(Immune               = "#E8C39E",
             Metabolism           = "#EBB9B0",
             Signaling            = "#C9BEDC",
             Proliferation        = "#AFD3E2",
             `DNA damage`         = "#8FA9C4",
             Development          = "#B5D6C4",
             `Cellular component` = "#9CC3A5")

FG_LEVELS <- c("Anti-tumor immune", "Pro-tumor immune",
               "Angiogenesis / fibroblasts", "Tumor state")
FG_COL <- c(`Anti-tumor immune`          = "#1F4E5F",
            `Pro-tumor immune`           = "#C08A2E",
            `Angiogenesis / fibroblasts` = "#8C2F39",
            `Tumor state`                = "#E8A9A9")

theme_mirror <- function(base = BASE_SIZE) {
  theme_bw(base_size = base) +
    theme(text            = element_text(colour = COL_TEXT),
          axis.text       = element_text(colour = COL_TEXT, size = base),
          axis.title      = element_text(colour = COL_TEXT, size = base),
          panel.grid      = element_blank(),
          panel.border    = element_rect(colour = COL_TEXT, linewidth = 0.3),
          axis.ticks      = element_line(colour = COL_TEXT, linewidth = 0.25),
          axis.ticks.length = unit(0.8, "mm"),
          legend.key.size = unit(2.4, "mm"),
          legend.text     = element_text(size = base),
          legend.title    = element_text(size = base),
          legend.margin   = margin(0, 0, 0, 0),
          legend.box.margin = margin(-2, 0, 0, 0),
          plot.title      = element_text(size = base, hjust = 0.5),
          plot.margin     = margin(1, 1, 1, 1, "mm"))
}

theme_axis_x <- function(base = BASE_SIZE) {
  theme_void(base_size = base) +
    theme(text          = element_text(colour = COL_TEXT),
          axis.text.x   = element_text(colour = COL_TEXT, size = base,
                                       margin = margin(t = 0.6, unit = "mm")),
          axis.title.x  = element_text(colour = COL_TEXT, size = base,
                                       margin = margin(t = 0.6, unit = "mm")),
          axis.ticks.x  = element_line(colour = COL_TEXT, linewidth = 0.25),
          axis.ticks.length = unit(0.8, "mm"),
          legend.key.size = unit(2.4, "mm"),
          legend.text   = element_text(size = base),
          legend.title  = element_text(size = base),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin   = margin(1, 1, 1, 1, "mm"))
}

theme_blank <- function(base = BASE_SIZE) {
  theme_void(base_size = base) +
    theme(text          = element_text(colour = COL_TEXT),
          legend.key.size = unit(2.4, "mm"),
          legend.text   = element_text(size = base),
          legend.title  = element_text(size = base),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin   = margin(1, 2, 1, 1, "mm"))
}

guide_zbar <- function(title = "Z-score")
  guide_colourbar(title = title, title.position = "left", title.vjust = 1,
                  barwidth = unit(11, "mm"), barheight = unit(1.6, "mm"),
                  ticks.colour = NA, frame.colour = COL_TEXT,
                  frame.linewidth = 0.15)



norm_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\.", "-", x)
  hit <- grepl("^TCGA-", x)
  x[hit] <- substr(x[hit], 1, 12)
  x
}


merge_dup_rows <- function(m, what = "") {
  if (!any(duplicated(rownames(m)))) return(m)
  n <- sum(duplicated(rownames(m)))
  s <- rowsum(m, rownames(m))
  m <- s / as.vector(table(rownames(m))[rownames(s)])
  message("    ", what, ": merged ", n, " duplicate gene symbols")
  m
}

##


build_symbol_map <- function(sym) {
  sym <- toupper(trimws(sym))
  u   <- unique(sym)
  map <- setNames(u, u)
  
  official <- suppressMessages(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"))
  need <- setdiff(u, official)
  if (!length(need)) return(map)
  
  ali  <- suppressMessages(AnnotationDbi::keys(org.Hs.eg.db, keytype = "ALIAS"))
  need <- intersect(need, ali)
  if (!length(need)) return(map)
  
  tb <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = need, keytype = "ALIAS", columns = "SYMBOL"))
  tb <- tb[!is.na(tb$SYMBOL) & tb$ALIAS != tb$SYMBOL, , drop = FALSE]
  if (!nrow(tb)) return(map)
  

  cnt_a  <- table(tb$ALIAS)
  drop_b <- !(tb$ALIAS %in% names(cnt_a)[cnt_a == 1])
  rej_b  <- tb[drop_b, , drop = FALSE]
  tb     <- tb[!drop_b, , drop = FALSE]
  

  drop_c <- tb$SYMBOL %in% u
  rej_c  <- tb[drop_c, , drop = FALSE]
  tb     <- tb[!drop_c, , drop = FALSE]
  

  cnt_s  <- table(tb$SYMBOL)
  drop_d <- tb$SYMBOL %in% names(cnt_s)[cnt_s > 1]
  rej_d  <- tb[drop_d, , drop = FALSE]
  tb     <- tb[!drop_d, , drop = FALSE]
  
  if (nrow(tb)) map[tb$ALIAS] <- tb$SYMBOL
  
  mk <- function(d, why) if (!nrow(d)) NULL else
    data.frame(old_symbol = d$ALIAS, new_symbol = d$SYMBOL,
               status = why, stringsAsFactors = FALSE)
  attr(map, "renamed")  <- tb
  attr(map, "audit")    <- do.call(rbind, Filter(Negate(is.null), list(
    mk(tb,    "renamed"),
    mk(rej_b, "rejected: alias maps to >1 gene"),
    mk(rej_c, "rejected: target symbol already present"),
    mk(rej_d, "rejected: target symbol claimed by >1 alias"))))
  map
}

apply_symbol_map <- function(m, map, what) {
  key <- toupper(trimws(rownames(m)))
  new <- unname(map[key])
  new[is.na(new)] <- key[is.na(new)]
  n_chg <- sum(new != key)
  rownames(m) <- new
  message("    ", what, ": ", n_chg, " row names were updated to current HGNC symbols")
  merge_dup_rows(m, what)
}

EXAMPLE_RENAMES <- c(IL8 = "CXCL8", INDO = "IDO1", ELA2 = "ELANE",
                     IL8RA = "CXCR1", IL8RB = "CXCR2", NOS2A = "NOS2",
                     SCYA3 = "CCL3", C10ORF54 = "VSIR",
                     MB21D1 = "CGAS", TMEM173 = "STING1")

report_example_renames <- function(map) {
  old <- names(EXAMPLE_RENAMES)
  got <- unname(map[old])
  ok  <- !is.na(got) & got == unname(EXAMPLE_RENAMES)
  logmsg("    Example symbol mapping check (presence and successful renaming):")
  logtable(data.frame(old_symbol = old,
                      expected   = unname(EXAMPLE_RENAMES),
                      applied    = ifelse(ok, "yes", "no"),
                      stringsAsFactors = FALSE))
  if (any(!ok))
    logmsg("    The old symbols marked applied = no have no probes on this platform; ",
           "their current counterparts cannot be evaluated here: ",
           paste(unname(EXAMPLE_RENAMES)[!ok], collapse = ", "))
}


read_matrix <- function(path, what) {
  if (!file.exists(path)) stop(what, " file does not exist: ", path)
  x  <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(x) < 2L || ncol(x) < 2L)
    stop(what, " must have identifiers and at least one expression column across two rows.")
  id <- as.character(x[[1]]); x <- x[, -1, drop = FALSE]
  if (anyNA(id) || any(!nzchar(trimws(id))))
    stop(what, " contains missing or empty row identifiers.")
  by_sample <- mean(grepl("^TCGA", toupper(id))) > 0.5
  m <- if (by_sample) t(as.matrix(x)) else as.matrix(x)
  if (by_sample) colnames(m) <- id else rownames(m) <- id
  storage.mode(m) <- "numeric"
  if (!any(is.finite(m))) stop(what, " contains no finite expression values.")
  colnames(m) <- norm_id(colnames(m))
  dup <- duplicated(colnames(m))
  if (any(dup)) {
    message("    ", what, ": removed ", sum(dup), " duplicate aliquots")
    m <- m[, !dup, drop = FALSE]
  }
  rownames(m) <- toupper(trimws(rownames(m)))
  m <- merge_dup_rows(m, what)
  message(sprintf("    %s: %d genes x %d samples (%s)", what, nrow(m), ncol(m),
                  if (by_sample) "source rows were samples; transposed" else "source rows were genes"))
  m
}


read_gmt_plain <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  if (!length(ln)) stop("The GMT file is empty: ", path)
  out <- lapply(strsplit(ln, "\t"), function(v) {
    g <- v[-c(1, 2)]; g <- g[nzchar(g)]
    if (!length(g)) NULL else
      data.frame(term = v[1], gene = g, stringsAsFactors = FALSE)
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) stop("The GMT file contains no gene sets with members: ", path)
  unique(do.call(rbind, out))
}


read_fges_gmt <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  lst <- lapply(strsplit(ln, "\t"), function(v) {
    g <- v[-c(1, 2)]; g <- unique(toupper(trimws(g[nzchar(g)])))
    list(id = v[1], label = if (length(v) >= 2 && nzchar(v[2])) v[2] else v[1],
         genes = g)
  })
  ids  <- vapply(lst, `[[`, "", "id")
  labs <- vapply(lst, `[[`, "", "label")
  sets <- lapply(lst, `[[`, "genes"); names(sets) <- ids
  list(sets = sets, label = setNames(labs, ids))
}

load_gmt <- function(path, what, export_symbol_gmt = TRUE) {
  if (!file.exists(path)) stop(what, " GMT file does not exist: ", path)
  g <- read_gmt_plain(path)
  if (mean(grepl("^[0-9]+$", g$gene)) > 0.9) {
    message("    ", what, ": identified Entrez IDs; mapping to gene symbols...")
    map <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = unique(g$gene),
      keytype = "ENTREZID", columns = "SYMBOL"))
    map <- map[!is.na(map$SYMBOL), ]
    g <- merge(g, map, by.x = "gene", by.y = "ENTREZID")
    g <- data.frame(term = g$term, gene = g$SYMBOL, stringsAsFactors = FALSE)
    if (export_symbol_gmt) {
      out <- sub("\\.entrez\\.gmt$", ".symbols.CONVERTED.gmt", path)
      if (out == path) out <- paste0(tools::file_path_sans_ext(path), ".symbols.gmt")
      lst <- split(g$gene, g$term)
      writeLines(vapply(names(lst), function(nm)
        paste(c(nm, "converted_from_entrez", unique(lst[[nm]])), collapse = "\t"),
        character(1)), out)
      message("    Exported gene-symbol GMT: ", out)
    }
  }
  g$gene <- toupper(trimws(g$gene))
  g <- unique(g)
  message("    ", what, ": total ", length(unique(g$term)), " gene sets")
  g
}


coverage_table <- function(sets, universe, source_tag, universe_before = NULL) {
  df <- data.frame(
    source  = source_tag,
    id      = names(sets),
    n_total = vapply(sets, length, 1L),
    n_hit   = vapply(sets, function(g) sum(g %in% universe), 1L),
    missing = vapply(sets, function(g)
      paste(setdiff(g, universe), collapse = ","), ""),
    stringsAsFactors = FALSE)
  df$pct <- round(100 * df$n_hit / df$n_total, 1)
  if (!is.null(universe_before)) {
    df$n_hit_before <- vapply(sets, function(g) sum(g %in% universe_before), 1L)
    df$pct_before   <- round(100 * df$n_hit_before / df$n_total, 1)
    df$pct_gain     <- round(df$pct - df$pct_before, 1)
    df <- df[, c("source", "id", "n_total",
                 "n_hit_before", "pct_before", "n_hit", "pct", "pct_gain",
                 "missing")]
  }
  df[order(df$pct), , drop = FALSE]
}


DISEASE_PAT <- paste0(
  "(ALZHEIMER|PARKINSON|HUNTINGTON|PRION|AMYOTROPHIC|NEURODEGENER",
  "|CARDIOMYOPATHY|MYOCARDITIS|ATHEROSCLEROSIS",
  "|DIABETES|THYROID|ASTHMA|ALLOGRAFT|GRAFT_VERSUS|AUTOIMMUNE|LUPUS|ARTHRITIS",
  "|INFECTION|LEISHMANIA|TOXOPLASMOSIS|MALARIA|TUBERCULOSIS|MEASLES|INFLUENZA",
  "|HEPATITIS|HERPES|PAPILLOMAVIRUS|_HIV|EPSTEIN|LEGIONELLOSIS|PERTUSSIS|SHIGELLOSIS",
  "|SALMONELLA|STAPHYLOCOCCUS|AMOEBIASIS|CHAGAS|COVID",
  "|CANCER|CARCINOMA|LEUKEMIA|MELANOMA|GLIOMA|ADENOMA|CARCINOGENESIS",
  "|ADDICTION|ALCOHOLISM|DEPRESSION|EPILEPSY|COCAINE|AMPHETAMINE|MORPHINE|NICOTINE)")

get_pathway_t2g <- function() {
  if (!HAS_MSIGDBR) return(list(t2g = NULL, label = NULL))
  grab <- function(tries, lab) {
    for (a in tries) {
      df <- tryCatch(do.call(msigdbr::msigdbr, a), error = function(e) NULL)
      if (!is.null(df) && nrow(df) > 0) {
        term <- if ("gs_name"     %in% names(df)) df$gs_name     else df$gs_id
        gene <- if ("gene_symbol" %in% names(df)) df$gene_symbol else df$db_gene_symbol
        return(list(t2g = data.frame(term = term, gene = toupper(gene),
                                     stringsAsFactors = FALSE), label = lab))
      }
    }
    NULL
  }
  r <- grab(list(
    list(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
    list(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
    list(species = "Homo sapiens", category   = "C2", subcategory   = "CP:KEGG"),
    list(species = "Homo sapiens", category   = "C2", subcategory   = "CP:KEGG_LEGACY")),
    "KEGG pathway")
  if (is.null(r))
    r <- grab(list(
      list(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
      list(species = "Homo sapiens", category   = "C5", subcategory   = "GO:BP")),
      "GO biological process")
  if (is.null(r)) return(list(t2g = NULL, label = NULL))
  
  n0 <- length(unique(r$t2g$term))
  if (DROP_DISEASE_KEGG && grepl("KEGG", r$label)) {
    keep <- !grepl(DISEASE_PAT, toupper(r$t2g$term))
    r$t2g <- r$t2g[keep, , drop = FALSE]
    message("    Pathway collection: ", r$label, ", ", n0, " sets; after excluding disease terms: ",
            length(unique(r$t2g$term)), " sets")
  } else message("    Pathway collection: ", r$label, ", total ", n0, " sets")
  r
}

pretty_term <- function(x) {
  s <- tolower(gsub("_", " ",
                    sub("^(KEGG_MEDICUS|KEGG|GOBP|REACTOME|HALLMARK)_", "", x)))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}
trunc_lab <- function(x, n = 40)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)


ora_phyper <- function(genes, t2g, universe, min_size = 10, max_size = 500) {
  t2g      <- t2g[t2g$gene %in% universe, , drop = FALSE]
  genes    <- intersect(unique(genes), universe)
  universe <- unique(universe)
  N <- length(universe); n <- length(genes)
  if (n < 10) return(NULL)
  sets <- split(t2g$gene, t2g$term)
  sz   <- vapply(sets, length, 1L)
  sets <- sets[sz >= min_size & sz <= max_size]
  if (!length(sets)) return(NULL)
  res <- lapply(names(sets), function(nm) {
    S <- unique(sets[[nm]]); K <- length(S)
    hit <- intersect(genes, S); k <- length(hit)
    if (k == 0) return(NULL)
    data.frame(Description = nm, GeneRatio = paste0(k, "/", n),
               BgRatio = paste0(K, "/", N), Count = k,
               FoldEnrich = (k / n) / (K / N),
               pvalue = stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               geneID = paste(hit, collapse = "/"), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(res) || !nrow(res)) return(NULL)
  res$p.adjust <- stats::p.adjust(res$pvalue, method = "BH")
  res[order(res$pvalue), ]
}

gsea_fgsea <- function(rank_vec, t2g, min_size = 10, max_size = 500) {
  sets <- split(t2g$gene, t2g$term)
  sets <- lapply(sets, function(g) intersect(unique(g), names(rank_vec)))
  sz   <- vapply(sets, length, 1L)
  sets <- sets[sz >= min_size & sz <= max_size]
  if (!length(sets)) stop("No gene sets passed size filtering; inspect the identifier convention.")
  fg <- suppressWarnings(as.data.frame(
    fgsea::fgsea(pathways = sets, stats = rank_vec,
                 minSize = min_size, maxSize = max_size, eps = 0)))
  data.frame(ID = fg$pathway, Description = fg$pathway, setSize = fg$size,
             NES = fg$NES, pvalue = fg$pval,
             p.adjust = stats::p.adjust(fg$pval, method = "BH"),
             stringsAsFactors = FALSE)
}

raster_df <- function(mat, xs, ys)
  data.frame(x = rep(xs, each = nrow(mat)),
             y = rep(ys, times = ncol(mat)),
             value = as.vector(mat))

de_wilcox <- function(mat, idx_pos, idx_neg) {
  lfc <- rowMeans(mat[, idx_pos, drop = FALSE]) -
    rowMeans(mat[, idx_neg, drop = FALSE])
  p <- if (HAS_MT) {
    matrixTests::row_wilcoxon_twosample(mat[, idx_pos, drop = FALSE],
                                        mat[, idx_neg, drop = FALSE])$pvalue
  } else {
    apply(mat, 1, function(v)
      tryCatch(wilcox.test(v[idx_pos], v[idx_neg])$p.value, error = function(e) NA))
  }
  data.frame(symbol = rownames(mat), log2fc = as.numeric(lfc),
             P = as.numeric(p), stringsAsFactors = FALSE)
}

de_lm <- function(mat, group, covars) {
  X <- model.matrix(~ group + ., data = as.data.frame(covars))
  j <- grep("^group", colnames(X))
  if (length(j) != 1) stop("The grouping coefficient is absent from the design matrix.")
  Y   <- t(mat)
  fit <- lm.fit(x = X, y = Y)
  b   <- fit$coefficients
  dfr <- nrow(X) - fit$rank
  s2  <- colSums(fit$residuals^2) / dfr
  XtXi <- chol2inv(chol(crossprod(X)))
  se  <- sqrt(s2 * XtXi[j, j])
  tt  <- b[j, ] / se
  data.frame(symbol = rownames(mat), log2fc = as.numeric(b[j, ]),
             P = 2 * stats::pt(-abs(tt), dfr), stringsAsFactors = FALSE)
}



FGES_CATMAP <- c(
  MHCI = "Anti-tumor immune", MHCII = "Anti-tumor immune",
  Coactivation_molecules = "Anti-tumor immune", Effector_cells = "Anti-tumor immune",
  T_cell_traffic = "Anti-tumor immune", NK_cells = "Anti-tumor immune",
  T_cells = "Anti-tumor immune", B_cells = "Anti-tumor immune",
  M1_signatures = "Anti-tumor immune", Th1_signature = "Anti-tumor immune",
  Antitumor_cytokines = "Anti-tumor immune",
  Checkpoint_inhibition = "Pro-tumor immune", Treg = "Pro-tumor immune",
  T_reg_traffic = "Pro-tumor immune", Neutrophil_signature = "Pro-tumor immune",
  Granulocyte_traffic = "Pro-tumor immune", MDSC = "Pro-tumor immune",
  MDSC_traffic = "Pro-tumor immune", Macrophages = "Pro-tumor immune",
  Macrophage_DC_traffic = "Pro-tumor immune", Th2_signature = "Pro-tumor immune",
  Protumor_cytokines = "Pro-tumor immune",
  CAF = "Angiogenesis / fibroblasts", Matrix = "Angiogenesis / fibroblasts",
  Matrix_remodeling = "Angiogenesis / fibroblasts",
  Angiogenesis = "Angiogenesis / fibroblasts",
  Endothelium = "Angiogenesis / fibroblasts",
  Proliferation_rate = "Tumor state", EMT_signature = "Tumor state")

BUILTIN_SIG <- list(
  `MHC I` = c("HLA-A","HLA-B","HLA-C","HLA-E","B2M","TAP1","TAP2","TAPBP","NLRC5","PSMB8","PSMB9"),
  `MHC II` = c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DMA","HLA-DMB","CD74","CIITA"),
  `Co-activation molecules` = c("CD27","CD28","CD40","CD40LG","CD80","CD86","ICOS","ICOSLG","TNFRSF4","TNFRSF9","TNFRSF14","CD226"),
  `Effector cells` = c("GZMA","GZMB","GZMK","PRF1","GNLY","NKG7","IFNG","TBX21","KLRD1","FASLG","CD8A","CD8B"),
  `T cells` = c("CD3D","CD3E","CD3G","CD2","CD5","CD6","TRAT1","ITK","LCK","ZAP70","THEMIS"),
  `T cell traffic` = c("CXCL9","CXCL10","CXCL11","CXCR3","CCL5","CCL4","CX3CL1","CXCL16"),
  `NK cells` = c("KLRF1","KLRD1","NCR1","NCR3","KIR2DL3","KIR2DL4","GNLY","NKG7","FGFBP2","SH2D1B"),
  `B cells` = c("CD19","MS4A1","CD79A","CD79B","TNFRSF13C","BLK","FCRL5","POU2AF1","CR2"),
  `M1 signature` = c("NOS2","IL12A","IL12B","IL23A","TNF","IL1B","CXCL9","CXCL10","SOCS3","CD80"),
  `Th1 signature` = c("TBX21","IFNG","IL12RB2","STAT4","CXCR3","IL2"),
  `Checkpoint molecules` = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","VSIR","IDO1"),
  `Treg` = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","CCR8","IL10"),
  `MDSC` = c("ARG1","CD33","ITGAM","IL4R","S100A8","S100A9","IL10","CXCL8","VEGFA"),
  `Neutrophil signature` = c("FCGR3B","CSF3R","CXCR2","S100A8","S100A9","ELANE","MPO","FPR1"),
  `Macrophages` = c("CD68","CD163","MSR1","MRC1","CSF1R","MARCO","SIGLEC1","VSIG4"),
  `Myeloid cells traffic` = c("CCL2","CCL7","CSF1","CXCL8","CXCL5","CCR2","SPP1"),
  `Immune suppression by myeloid cells` = c("ARG1","IDO1","IL10","CD274","VEGFA","PTGS2","NT5E","ENTPD1"),
  `Protumor cytokines` = c("IL10","TGFB1","TGFB2","TGFB3","IL6","VEGFA","CSF1","IL4","IL13"),
  `Th2 signature` = c("GATA3","IL4","IL5","IL13","CCR4","IL10"),
  `Angiogenesis` = c("VEGFA","VEGFB","VEGFC","KDR","FLT1","TEK","ANGPT1","ANGPT2","PGF","ESM1"),
  `Endothelium` = c("PECAM1","CDH5","VWF","CLDN5","TIE1","ENG","ROBO4","EMCN"),
  `Cancer-associated fibroblasts` = c("FAP","PDGFRA","PDGFRB","ACTA2","COL1A1","COL1A2","THY1","TAGLN","POSTN"),
  `Matrix` = c("COL1A1","COL1A2","COL3A1","COL5A1","FN1","LUM","VCAN","BGN","SPARC","ELN"),
  `Matrix remodeling` = c("MMP1","MMP2","MMP9","MMP11","MMP14","TIMP1","PLOD2","LOX","CTSK"),
  `Tumor proliferation rate` = c("MKI67","CCNB1","CDK1","TOP2A","AURKA","BUB1","PLK1","TYMS","RRM2","PCNA"),
  `EMT signature` = c("VIM","CDH2","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","FN1","MMP2","SPARC"))

BUILTIN_CAT <- c(
  `MHC I` = "Anti-tumor immune", `MHC II` = "Anti-tumor immune",
  `Co-activation molecules` = "Anti-tumor immune", `Effector cells` = "Anti-tumor immune",
  `T cells` = "Anti-tumor immune", `T cell traffic` = "Anti-tumor immune",
  `NK cells` = "Anti-tumor immune", `B cells` = "Anti-tumor immune",
  `M1 signature` = "Anti-tumor immune", `Th1 signature` = "Anti-tumor immune",
  `Checkpoint molecules` = "Pro-tumor immune", `Treg` = "Pro-tumor immune",
  `MDSC` = "Pro-tumor immune", `Neutrophil signature` = "Pro-tumor immune",
  `Macrophages` = "Pro-tumor immune", `Myeloid cells traffic` = "Pro-tumor immune",
  `Immune suppression by myeloid cells` = "Pro-tumor immune",
  `Protumor cytokines` = "Pro-tumor immune", `Th2 signature` = "Pro-tumor immune",
  `Angiogenesis` = "Angiogenesis / fibroblasts", `Endothelium` = "Angiogenesis / fibroblasts",
  `Cancer-associated fibroblasts` = "Angiogenesis / fibroblasts",
  `Matrix` = "Angiogenesis / fibroblasts", `Matrix remodeling` = "Angiogenesis / fibroblasts",
  `Tumor proliferation rate` = "Tumor state", `EMT signature` = "Tumor state")

DDR_SETS <- list(
  `Homologous recombination` = c("BRCA1","BRCA2","RAD51","RAD51B","RAD51C","RAD51D",
                                 "PALB2","XRCC2","XRCC3","BARD1","BRIP1","ATM","ATR",
                                 "CHEK1","CHEK2","MRE11","RAD50","NBN","EXO1","BLM","RBBP8"),
  `Fanconi anemia`           = c("FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG",
                                 "FANCI","FANCL","FANCM","UBE2T","SLX4","ERCC4"),
  `Non-homologous end joining` = c("PRKDC","XRCC4","XRCC5","XRCC6","LIG4","NHEJ1","DCLRE1C"),
  `Base excision repair`     = c("OGG1","APEX1","POLB","XRCC1","LIG3","PARP1","PARP2","MUTYH","UNG","NEIL1"),
  `Mismatch repair`          = c("MLH1","MSH2","MSH6","PMS2","MSH3","PMS1","EXO1"),
  `cGAS-STING`               = c("CGAS","STING1","TBK1","IRF3","IKBKE","IFI16","DDX41"))

HR_GENES <- c("BRCA1","BRCA2","RAD51C","RAD51D","PALB2","ATM","CHEK2",
              "FANCA","FANCI","PTEN","EMSY","MKI67")

run_analysis <- function(tag, samples, grp01, adjust = FALSE) {
  logmsg("")
  logmsg("===== Analysis: ", tag, " (grouping label = ", LABEL_COL, ")=====")
  keep <- !is.na(grp01)
  samples <- samples[keep]; grp01 <- grp01[keep]
  
  cl <- clin[match(samples, clin$sample_id), ]
  cl$group <- factor(ifelse(grp01 == 1, "POS", "NEG"), levels = c("NEG", "POS"))
  cl  <- cl[order(cl$group, cl$HRD_prob), ]
  sm  <- cl$sample_id
  
  E_rma <- expr_rma[, sm, drop = FALSE]
  E_z   <- expr_z[,   sm, drop = FALSE]
  grp   <- cl$group
  idx_neg <- which(grp == "NEG"); idx_pos <- which(grp == "POS")
  if (length(idx_neg) < 10 || length(idx_pos) < 10) {
    logmsg("    Insufficient samples in a group (NEG = ", length(idx_neg),
           ", POS = ", length(idx_pos), "); analysis skipped.")
    return(NULL)
  }
  logmsg(sprintf("    Included %d samples (%s %d / %s %d), %d genes%s",
                 length(sm), LAB_NEG, length(idx_neg), LAB_POS, length(idx_pos),
                 nrow(E_rma), if (adjust) "; immune/stromal scores adjusted" else ""))
  
  deg <- if (adjust) {
    de_lm(E_rma, grp, data.frame(ImmuneScore  = cl$ImmuneScore,
                                 StromalScore = cl$StromalScore))
  } else de_wilcox(E_rma, idx_pos, idx_neg)
  deg$FDR <- p.adjust(deg$P, method = "BH")
  if (all(is.na(deg$P))) stop("All differential-expression estimates are missing; inspect the matrix.")

  n_floor <- sum(deg$FDR < FDR_FLOOR, na.rm = TRUE)
  if (n_floor) {
    logmsg("    ", n_floor, " gene FDR values underflowed and were bounded at ", FDR_FLOOR)
    deg$FDR <- pmax(deg$FDR, FDR_FLOOR)
  }
  
  lfc_med <- median(deg$log2fc, na.rm = TRUE)
  logmsg(sprintf("    Median log2FC = %+.4f; %.1f%% of genes > 0; %d genes with FDR < 0.05",
                 lfc_med, 100 * mean(deg$log2fc > 0, na.rm = TRUE),
                 sum(deg$FDR < 0.05, na.rm = TRUE)))
  qq <- round(quantile(abs(deg$log2fc), c(.5, .9, .95, .99, 1), na.rm = TRUE), 3)
  logmsg("    Absolute log2FC quantiles: ", paste(names(qq), qq, sep = "=", collapse = "  "))
  
  for (v in c("ImmuneScore", "StromalScore")) {
    p <- tryCatch(wilcox.test(cl[[v]][idx_pos], cl[[v]][idx_neg])$p.value,
                  error = function(e) NA)
    logmsg(sprintf("    %s: median %s %+.3f, median %s %+.3f, Wilcoxon P = %.3g",
                   v, LAB_POS, median(cl[[v]][idx_pos], na.rm = TRUE),
                   LAB_NEG, median(cl[[v]][idx_neg], na.rm = TRUE), p))
  }
  rp <- suppressWarnings(cor(cl$ImmuneScore, cl$HRD_prob,
                             method = "spearman", use = "complete.obs"))
  logmsg(sprintf("    Spearman rho between ImmuneScore and HRD_prob = %+.3f", rp))
  

  hg   <- intersect(HR_GENES, deg$symbol)
  miss <- setdiff(HR_GENES, deg$symbol)
  if (length(miss))
    logmsg("    These HR/proliferation genes lack probes and were excluded: ",
           paste(miss, collapse = ", "))
  if (length(hg)) {
    sub <- deg[match(hg, deg$symbol), c("symbol", "log2fc", "FDR")]
    sub$log2fc <- round(sub$log2fc, 3); sub$FDR <- signif(sub$FDR, 3)
    logmsg("    HR/proliferation genes: ")
    logtable(sub)
  }
  
  deg$sig <- "no"
  deg$sig[deg$FDR <= FDR_CUT & deg$log2fc >=  LFC_CUT] <- "up"
  deg$sig[deg$FDR <= FDR_CUT & deg$log2fc <= -LFC_CUT] <- "down"
  deg <- deg[order(deg$FDR), ]
  rownames(deg) <- deg$symbol
  write.csv(deg, paste0(tag, "_DEG_table.csv"), row.names = FALSE)
  n_up <- sum(deg$sig == "up"); n_down <- sum(deg$sig == "down")
  logmsg(sprintf("    %s-high genes: %d; %s-high genes: %d", LAB_POS, n_up, LAB_NEG, n_down))
  

  lfc_c  <- deg$log2fc - lfc_med
  n_up_c <- sum(deg$FDR <= FDR_CUT & lfc_c >=  LFC_CUT, na.rm = TRUE)
  n_dn_c <- sum(deg$FDR <= FDR_CUT & lfc_c <= -LFC_CUT, na.rm = TRUE)
  logmsg(sprintf(paste0("    Sensitivity analysis after subtracting global median log2FC %+.4f: ",
                        "%s-high %d (primary %d); %s-high %d (primary %d). ",
                        "Use this only to describe the global shift; retain the primary analysis."),
                 lfc_med, LAB_POS, n_up_c, n_up, LAB_NEG, n_dn_c, n_down))
  
  if (n_up + n_down == 0) {
    logmsg("    No genes passed the threshold; figure generation skipped.")
    return(NULL)
  }
  deg$sig <- factor(deg$sig, levels = c("no", "down", "up"))
  
  ora_up   <- ora_phyper(deg$symbol[deg$sig == "up"],   path_t2g, universe)
  ora_down <- ora_phyper(deg$symbol[deg$sig == "down"], path_t2g, universe)
  if (!is.null(ora_up))   write.csv(ora_up,   paste0(tag, "_ORA_HRDpos.csv"), row.names = FALSE)
  if (!is.null(ora_down)) write.csv(ora_down, paste0(tag, "_ORA_HRDneg.csv"), row.names = FALSE)
  
  rank_vec <- deg$log2fc; names(rank_vec) <- deg$symbol
  rank_vec <- sort(rank_vec[!is.na(rank_vec)], decreasing = TRUE)
  gse_df <- gsea_fgsea(rank_vec, hall_t2g)
  write.csv(gse_df, paste0(tag, "_Hallmark_GSEA.csv"), row.names = FALSE)
  

  if (length(DDR_SETS_KEPT)) {
    ddr_t2g <- do.call(rbind, lapply(names(DDR_SETS_KEPT), function(nm)
      data.frame(term = nm, gene = DDR_SETS_KEPT[[nm]], stringsAsFactors = FALSE)))
    ddr <- tryCatch(gsea_fgsea(rank_vec, ddr_t2g, min_size = DDR_MIN_GENES,
                               max_size = 100),
                    error = function(e) NULL)
    if (!is.null(ddr)) {
      write.csv(ddr, paste0(tag, "_DDR_GSEA.csv"), row.names = FALSE)
      logmsg("    DDR GSEA (positive NES indicates HRD+ enrichment; coverage-filtered pathways):")
      d2 <- ddr[, c("Description", "setSize", "NES", "p.adjust")]
      d2$NES <- round(d2$NES, 2); d2$p.adjust <- signif(d2$p.adjust, 3)
      logtable(d2)
    }
  }

  S <- ss[, sm, drop = FALSE]
  sig_stat <- NULL
  if (adjust && ADJ_DROP_E) {
    logmsg("    Adjusted analysis omits panel E: signature scores lack covariate adjustment; ",
           "they are neither recalculated nor exported for this comparison.")
  } else {
    sig_stat <- data.frame(
      signature = rownames(S),
      P = apply(S, 1, function(v)
        tryCatch(wilcox.test(v[idx_pos], v[idx_neg])$p.value, error = function(e) NA)),
      diff = rowMeans(S[, idx_pos, drop = FALSE]) - rowMeans(S[, idx_neg, drop = FALSE]),
      stringsAsFactors = FALSE)
    sig_stat$FDR <- p.adjust(sig_stat$P, method = "BH")
    sig_stat$dir <- ifelse(is.na(sig_stat$FDR) | sig_stat$FDR >= 0.05, "ns",
                           ifelse(sig_stat$diff > 0, "up_POS", "up_NEG"))
    sig_stat$cat <- factor(SIG_CAT[sig_stat$signature], levels = FG_LEVELS)
    sig_stat$cat[is.na(sig_stat$cat)] <- "Tumor state"
    sig_stat <- sig_stat %>% arrange(cat, desc(diff)) %>%
      mutate(ypos = rev(seq_len(n())))
    out_nm <- if (adjust)
      paste0(tag, "_signature_stats_UNADJUSTED_scores.csv") else
        paste0(tag, "_signature_stats.csv")
    write.csv(sig_stat, out_nm, row.names = FALSE)
    if (adjust)
      logmsg("    Panel E signature scores are unadjusted (UNADJUSTED_scores in filename).")
    logmsg("    Panel E: ", nrow(sig_stat), " signatures; FDR < 0.05: ",
           sum(sig_stat$dir != "ns"), " sets")
  }
  
  list(tag = tag, adjust = adjust, clin = cl, grp = grp,
       idx_neg = idx_neg, idx_pos = idx_pos, samples = sm,
       expr_z = E_z, deg = deg, n_up = n_up, n_down = n_down,
       ora_up = ora_up, ora_down = ora_down, gse = gse_df,
       ss = S, sig_stat = sig_stat)
}

hallmark_cat <- tibble::tribble(
  ~ID, ~cat,
  "HALLMARK_ALLOGRAFT_REJECTION","Immune","HALLMARK_COAGULATION","Immune",
  "HALLMARK_COMPLEMENT","Immune","HALLMARK_INTERFERON_ALPHA_RESPONSE","Immune",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","Immune","HALLMARK_IL6_JAK_STAT3_SIGNALING","Immune",
  "HALLMARK_INFLAMMATORY_RESPONSE","Immune",
  "HALLMARK_BILE_ACID_METABOLISM","Metabolism","HALLMARK_CHOLESTEROL_HOMEOSTASIS","Metabolism",
  "HALLMARK_FATTY_ACID_METABOLISM","Metabolism","HALLMARK_GLYCOLYSIS","Metabolism",
  "HALLMARK_HEME_METABOLISM","Metabolism","HALLMARK_OXIDATIVE_PHOSPHORYLATION","Metabolism",
  "HALLMARK_XENOBIOTIC_METABOLISM","Metabolism",
  "HALLMARK_ANDROGEN_RESPONSE","Signaling","HALLMARK_ESTROGEN_RESPONSE_EARLY","Signaling",
  "HALLMARK_ESTROGEN_RESPONSE_LATE","Signaling","HALLMARK_IL2_STAT5_SIGNALING","Signaling",
  "HALLMARK_KRAS_SIGNALING_UP","Signaling","HALLMARK_KRAS_SIGNALING_DN","Signaling",
  "HALLMARK_MITOTIC_SPINDLE","Signaling","HALLMARK_NOTCH_SIGNALING","Signaling",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING","Signaling","HALLMARK_HEDGEHOG_SIGNALING","Signaling",
  "HALLMARK_TGF_BETA_SIGNALING","Signaling","HALLMARK_TNFA_SIGNALING_VIA_NFKB","Signaling",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING","Signaling","HALLMARK_APOPTOSIS","Signaling",
  "HALLMARK_HYPOXIA","Signaling","HALLMARK_PROTEIN_SECRETION","Signaling",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","Signaling",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","Signaling",
  "HALLMARK_E2F_TARGETS","Proliferation","HALLMARK_G2M_CHECKPOINT","Proliferation",
  "HALLMARK_MYC_TARGETS_V1","Proliferation","HALLMARK_MYC_TARGETS_V2","Proliferation",
  "HALLMARK_P53_PATHWAY","Proliferation","HALLMARK_MTORC1_SIGNALING","Proliferation",
  "HALLMARK_DNA_REPAIR","DNA damage","HALLMARK_UV_RESPONSE_DN","DNA damage",
  "HALLMARK_UV_RESPONSE_UP","DNA damage",
  "HALLMARK_ADIPOGENESIS","Development","HALLMARK_ANGIOGENESIS","Development",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","Development","HALLMARK_MYOGENESIS","Development",
  "HALLMARK_PANCREAS_BETA_CELLS","Development","HALLMARK_SPERMATOGENESIS","Development",
  "HALLMARK_PEROXISOME","Cellular component","HALLMARK_APICAL_JUNCTION","Cellular component",
  "HALLMARK_APICAL_SURFACE","Cellular component")


A4_W <- 297
A4_H <- 210
FONT_PT <- 8
FONT_FAMILY <- "Helvetica"

mm <- function(x) grid::unit(x, "mm")
text_gp <- function(col = COL_TEXT, face = "plain")
  grid::gpar(fontfamily = FONT_FAMILY, fontsize = FONT_PT,
             fontface = face, col = col, lineheight = 1.12)
text_width <- function(label, face = "plain") {
  vapply(label, function(s) grid::convertWidth(
    grid::grobWidth(grid::textGrob(s, gp = text_gp(face = face))),
    "mm", valueOnly = TRUE), numeric(1))
}
txt <- function(label, x, y, col = COL_TEXT, face = "plain",
                hjust = 0, vjust = 0.5, rot = 0) {
  grid::grid.text(label, x = mm(x), y = mm(A4_H - y),
                  hjust = hjust, vjust = vjust, rot = rot,
                  gp = text_gp(col, face))
}
box <- function(x, y, w, h, fill, border = NA, lwd = 0.55) {
  grid::grid.rect(x = mm(x), y = mm(A4_H - y), width = mm(w), height = mm(h),
                  just = c("left", "top"), gp = grid::gpar(fill = fill, col = border, lwd = lwd))
}
line_mm <- function(x0, y0, x1, y1, col = COL_TEXT, lty = 1, lwd = 0.55,
                    arrow = NULL) {
  grid::grid.segments(mm(x0), mm(A4_H-y0), mm(x1), mm(A4_H-y1),
                      arrow = arrow, gp = grid::gpar(col = col, lty = lty, lwd = lwd, fill = col))
}
wrap_mm <- function(label, width, face = "plain") {
  vapply(label, function(s) {
    words <- strsplit(s, " ", fixed = TRUE)[[1]]
    out <- character(); current <- ""
    for (word in words) {
      candidate <- if (nzchar(current)) paste(current, word) else word
      if (nzchar(current) && text_width(candidate, face) > width) {
        out <- c(out, current); current <- word
      } else current <- candidate
    }
    out <- c(out, current)
    if (any(text_width(out, face) > width)) stop("Label is wider than its allocated column: ", s)
    paste(out, collapse = "\n")
  }, character(1))
}
key <- function(label, color, x, y, face = "plain") {
  box(x, y - 1.25, 3.5, 2.5, color)
  txt(label, x + 4.8, y, face = face)
}
zlegend <- function(x, y, w = 32) {
  txt("Z-score", x, y)
  cols <- matrix(colorRampPalette(c("#9C4DCC", "#FFFFFF", "#FFD702"))(256), nrow = 1)
  grid::grid.raster(as.raster(cols), x = mm(x), y = mm(A4_H-y-2.7),
                    width = mm(w), height = mm(2.8), just = c("left", "top"), interpolate = FALSE)
  box(x, y+2.7, w, 2.8, NA, COL_TEXT, 0.3)
  txt(format(-HEAT_LIM), x, y+8.1, hjust = 0.5)
  txt("0", x+w/2, y+8.1, hjust = 0.5)
  txt(format(HEAT_LIM), x+w, y+8.1, hjust = 0.5)
}
heat <- function(m, x, y, w, h) {
  if (!nrow(m) || !ncol(m)) return(invisible(NULL))
  ix <- 1L + round((pmin(HEAT_LIM, pmax(-HEAT_LIM, m)) + HEAT_LIM) /
                     (2*HEAT_LIM) * (length(HEAT_COL)-1L))
  colors <- matrix(HEAT_COL[ix], nrow(m), ncol(m))
  colors[is.na(m)] <- "#EEEEEE"
  grid::grid.raster(as.raster(colors), x = mm(x), y = mm(A4_H-y),
                    width = mm(w), height = mm(h), just = c("left", "top"), interpolate = FALSE)
}
group_strip <- function(R, x, y, w, h = 5, gap = 2) {
  n <- c(length(R$idx_neg), length(R$idx_pos))
  widths <- (w-gap) * n/sum(n)
  starts <- c(x, x+widths[1]+gap)
  for (i in 1:2) {
    box(starts[i], y, widths[i], h, c(COL_NEG, COL_POS)[i])
    txt(c("HRD-", "HRD+")[i], starts[i]+widths[i]/2, y+h/2,
        col = "white", hjust = 0.5)
  }
  list(x = starts, w = widths)
}


draw_deg_heat <- function(R) {
  txt("A", 10, 12, face = "bold")
  deg <- R$deg; rownames(deg) <- deg$symbol
  up <- deg$symbol[deg$sig == "up"]
  down <- deg$symbol[deg$sig == "down"]
  if (length(up)+length(down) > MAX_HEATMAP_ROWS) {
    k <- floor(MAX_HEATMAP_ROWS*length(up)/(length(up)+length(down)))
    up <- head(up, max(1,k)); down <- head(down, max(1,MAX_HEATMAP_ROWS-k))
  }
  up <- intersect(up[order(-deg[up,"log2fc"])], rownames(R$expr_z))
  down <- intersect(down[order(deg[down,"log2fc"])], rownames(R$expr_z))
  if (!length(up) && !length(down)) stop("No differential-expression genes are available to plot.")
  gb <- group_strip(R, 10, 20, 130, 6, 2.8)
  n <- c(length(up),length(down)); gap <- if (all(n>0)) 2.8 else 0
  heights <- (154-gap)*n/sum(n)
  ypos <- c(29,29+heights[1]+gap)
  for (i in 1:2) for (j in 1:2) {
    genes <- list(up,down)[[i]]
    samples <- list(R$idx_neg,R$idx_pos)[[j]]
    heat(R$expr_z[genes,samples,drop=FALSE], gb$x[j], ypos[i], gb$w[j], heights[i])
  }
  zlegend(107, 189, 32)
}


draw_kegg <- function(R) {
  txt("B", 148, 12, face = "bold")
  box(148,20,139,6,COL_PANEL)
  txt(PATH_LAB,217.5,23,hjust=0.5)
  parts <- lapply(list(R$ora_down,R$ora_up), function(d) {
    if (is.null(d) || !nrow(d)) return(data.frame())
    d <- head(d[order(d$p.adjust),,drop=FALSE],N_PATH_EACH)
    d$label <- pretty_term(d$Description)
    d$label <- gsub("\\bRna\\b", "RNA",d$label)
    d$label <- gsub("\\bDna\\b", "DNA",d$label)
    d$label <- gsub("\\bdna\\b", "DNA",d$label)
    d$label <- gsub("\\brna\\b", "RNA",d$label)
    d$label <- gsub("\\btca\\b", "TCA",d$label)
    d$score <- -log10(pmax(d$p.adjust,FDR_FLOOR))
    d
  })
  vals <- unlist(lapply(parts,function(d) d$score))
  xmax <- max(c(vals,-log10(FDR_CUT)),na.rm=TRUE)*1.08
  bx <- 235; bw <- 52; lx <- 157; lw <- bx-lx-2
  xx <- function(v) bx+bw*v/xmax
  yt <- c(29,108); gh <- 71
  for (k in 1:2) {
    d <- parts[[k]]; box(148,yt[k],6,gh,COL_PANEL)
    txt(c("HRD- enriched","HRD+ enriched")[k],151,yt[k]+gh/2,
        hjust=0.5,rot=90)
    if (!nrow(d)) {
      txt("No eligible pathways",157,yt[k]+gh/2); next
    }
    labs <- wrap_mm(d$label,lw)
    counts <- lengths(strsplit(labs,"\n",fixed=TRUE))
    heights <- pmax(1,counts)*3.4
    if(sum(heights)>gh) stop("KEGG labels exceed the A4 height; reduce N_PATH_EACH.")
    heights <- heights+(gh-sum(heights))/nrow(d)
    cy <- yt[k]+cumsum(heights)-heights/2
    box(bx,yt[k],bw,gh,NA,"#B3B6BC",0.4)
    line_mm(xx(-log10(FDR_CUT)),yt[k],xx(-log10(FDR_CUT)),yt[k]+gh,
            col="#60646B",lty=2,lwd=0.55)
    for(i in seq_len(nrow(d))) {
      fill <- if(is.na(d$p.adjust[i]) || d$p.adjust[i]>=FDR_CUT) COL_NS else c(COL_NEG_L,COL_POS_L)[k]
      box(bx,cy[i]-1.75,bw*d$score[i]/xmax,3.5,fill)
      txt(labs[i],lx,cy[i])
    }
  }
  ticks <- pretty(c(0,xmax),n=4);ticks <- ticks[ticks>=0 & ticks<=xmax]
  for(v in ticks) {
    line_mm(xx(v),179,xx(v),180.3)
    txt(format(v,trim=TRUE),xx(v),183.3,hjust=0.5)
  }
  txt("-log10(FDR)",bx+bw/2,190,hjust=0.5)
  txt(paste0("Dashed line: FDR = ",FDR_CUT),157,199)
}


draw_volcano <- function(R) {
  txt("A",10,12,face="bold")
  d <- R$deg
  d <- d[is.finite(d$log2fc) & is.finite(d$FDR),,drop=FALSE]
  if(!nrow(d)) stop("The volcano panel contains no valid data.")
  d$FDR <- pmax(d$FDR,FDR_FLOOR)
  xr <- max(abs(d$log2fc),LFC_CUT)*1.08
  ymax <- max(-log10(d$FDR),-log10(FDR_CUT))*1.40
  x<-22;y<-20;w<-59;h<-87
  xx<-function(v) x+(v+xr)/(2*xr)*w
  yy<-function(v) y+h-v/ymax*h
  box(x,y,w,h,NA,COL_TEXT)
  line_mm(x,yy(-log10(FDR_CUT)),x+w,yy(-log10(FDR_CUT)),"#777777",2)
  for(v in c(-LFC_CUT,LFC_CUT)) line_mm(xx(v),y+20,xx(v),y+h,"#777777",2)
  for(s in c("no","down","up")) {
    v<-d[d$sig==s,,drop=FALSE]
    color<-c(no=COL_NS,down=COL_NEG,up=COL_POS)[s]
    grid::grid.points(mm(xx(v$log2fc)),mm(A4_H-yy(-log10(v$FDR))),pch=21,
                      size=mm(if(s=="no") 0.80 else 1.80),
                      gp=grid::gpar(fill=grDevices::adjustcolor(color,alpha.f=if(s=="no") 0.5 else 0.9),
                                    col=if(s=="no") NA else "white",lwd=0.25))
  }
  txt(paste0("HRD- associated\ngenes (n = ",R$n_down,")"),x+2,y+5,
      col=COL_NEG)
  txt(paste0("HRD+ associated\ngenes (n = ",R$n_up,")"),x+w-2,y+5,
      col=COL_POS,hjust=1)
  ar<-grid::arrow(length=mm(1.5),type="closed")
  line_mm(x+w*.35,y+15,x+4,y+15,COL_NEG,lwd=0.7,arrow=ar)
  line_mm(x+w*.65,y+15,x+w-4,y+15,COL_POS,lwd=0.7,arrow=ar)
  xt<-pretty(c(-xr,xr),n=3);xt<-xt[xt>=-xr & xt<=xr]
  yt<-pretty(c(0,ymax*.9),n=4);yt<-yt[yt>=0 & yt<=ymax*.9]
  for(v in xt) {line_mm(xx(v),y+h,xx(v),y+h+1);txt(format(v,trim=TRUE),xx(v),y+h+4,hjust=.5)}
  for(v in yt) {line_mm(x-1,yy(v),x,yy(v));txt(format(v,trim=TRUE),x-2,yy(v),hjust=1)}
  txt(if(isTRUE(R$adjust)) "Adjusted log2(FC)" else "log2(FC)",x+w/2,118,hjust=.5)
  txt("-log10(FDR)",12,y+h/2,hjust=.5,rot=90)
}


draw_signatures <- function(R) {
  txt("B",92,12,face="bold")
  s<-R$sig_stat
  if(is.null(s)||!nrow(s)) {
    txt("Signature scores are not covariate-adjusted.",100,50)
    txt("This panel is omitted for the adjusted analysis.",100,55)
    return(invisible(NULL))
  }
  s<-s[order(-s$ypos),,drop=FALSE]
  sh<-t(scale(t(R$ss[s$signature,,drop=FALSE])))
  labels<-wrap_mm(s$signature,58)
  nl<-lengths(strsplit(labels,"\n",fixed=TRUE))
  totalh<-84
  rh<-pmax(1,nl)*2.88
  if(sum(rh)>totalh) stop("Signature labels exceed the 8-pt layout; revise the geometry.")
  rh<-rh+(totalh-sum(rh))/nrow(s)
  tops<-24+c(0,head(cumsum(rh),-1));cy<-tops+rh/2
  gb<-group_strip(R,100,17,109,5,1.6)

  for(i in seq_len(nrow(s))) {
    box(94,tops[i],4,rh[i],unname(FG_COL[as.character(s$cat[i])]))
    for(j in 1:2) heat(sh[i,list(R$idx_neg,R$idx_pos)[[j]],drop=FALSE],
                       gb$x[j],tops[i],gb$w[j],rh[i])
    direction<-as.character(s$dir[i])
    col<-if(direction=="up_POS") COL_POS else if(direction=="up_NEG") COL_NEG else "#888888"
    fdr<-if(is.na(s$FDR[i])) "NA" else if(direction=="ns") "ns" else
      if(s$FDR[i]<.001) "< 0.001" else sprintf("%.3f",s$FDR[i])
    txt(fdr,211,cy[i],col=col)
    if(direction %in% c("up_POS","up_NEG")) {
      dy<-if(direction=="up_POS") -0.8 else 0.8
      line_mm(224.3,cy[i]-dy,224.3,cy[i]+dy,col,lwd=.55,
              arrow=grid::arrow(length=mm(.7),type="open"))
    }
    txt(labels[i],228,cy[i])
  }
  txt("FDR",211,19.5)
  key(FG_LEVELS[1],FG_COL[1],94,113)
  key(FG_LEVELS[2],FG_COL[2],94,118.5)
  key(FG_LEVELS[3],FG_COL[3],134,113)
  key(FG_LEVELS[4],FG_COL[4],134,118.5)
  zlegend(250,113,35)
}

hallmark_data <- function(R) {
  d<-merge(as.data.frame(hallmark_cat),R$gse[,c("ID","NES","p.adjust")],by="ID",all=TRUE,sort=FALSE)
  d$cat[is.na(d$cat)]<-"Cellular component"
  d$missing<-is.na(d$NES)|is.na(d$p.adjust)
  d$type<-ifelse(d$missing,"missing",ifelse(d$p.adjust>=FDR_CUT,"ns",ifelse(d$NES>0,"pos","neg")))
  d$cat<-factor(d$cat,levels=CAT_LEVELS)

  d$label<-pretty_term(d$ID)
  acronym<-c("DNA","UV","RNA","TNFA","NFKB","IL2","IL6","JAK","STAT3","STAT5",
             "PI3K","AKT","MTOR","MTORC1","KRAS","MYC","E2F","G2M","EMT","TGF","DN","UP")
  for(a in acronym) d$label<-gsub(paste0("\\b",tolower(a),"\\b"),a,d$label,ignore.case=TRUE)
  d$label<-gsub("\\bv1\\b","V1",d$label,ignore.case=TRUE)
  d$label<-gsub("\\bv2\\b","V2",d$label,ignore.case=TRUE)
  if(ORDER_D=="alpha") d<-d[order(d$cat,d$label),] else d<-d[order(d$cat,-d$NES,na.last=TRUE),]
  d
}


draw_hallmark <- function(R) {
  txt("C",10,123,face="bold")
  d<-hallmark_data(R);n<-nrow(d)
  x<-13;w<-229;y<-126;bottom<-164
  step<-w/n
  if(step < FONT_PT/72*25.4*1.08) stop("Hallmark labels exceed the horizontal 8-pt budget.")
  faces<-ifelse(d$type %in% c("neg","pos"),"bold","plain")
  widths<-mapply(function(s,f) text_width(s,f),d$label,faces)

  label_top<-bottom+5
  if(label_top+max(widths)>201) {
    bottom<-201-max(widths)-5;label_top<-bottom+5
  }
  if(bottom-y<19) stop("Hallmark labels exceed the A4 vertical budget.")
  values<-d$NES[!d$missing]
  hi<-max(c(0,values))*1.08;lo<-min(c(0,values))*1.08
  if(hi==lo) {hi<-1;lo<--1}
  yy<-function(v) bottom-(v-lo)/(hi-lo)*(bottom-y)
  xc<-x+(seq_len(n)-.5)*step
  for(i in seq_len(n)) {
    col<-c(neg=COL_NEG,pos=COL_POS,ns=COL_NS,missing="#FFFFFF")[d$type[i]]
    if(!d$missing[i]) {
      box(xc[i]-step*.42,min(yy(0),yy(d$NES[i])),step*.84,
          abs(yy(d$NES[i])-yy(0)),col)
    } else {
      line_mm(xc[i]-.6,yy(0)-.6,xc[i]+.6,yy(0)+.6,"#888888")
      line_mm(xc[i]-.6,yy(0)+.6,xc[i]+.6,yy(0)-.6,"#888888")
    }
    box(x+(i-1)*step,bottom+2,step,2,unname(CAT_COL[as.character(d$cat[i])]))
    labelcol<-if(d$type[i]=="neg") COL_NEG else if(d$type[i]=="pos") COL_POS else COL_TEXT
    txt(d$label[i],xc[i],label_top,col=labelcol,face=faces[i],hjust=1,vjust=.5,rot=90)
  }
  line_mm(x,yy(0),x+w,yy(0))
  line_mm(x+w,y,x+w,bottom)
  yt<-pretty(c(lo,hi),n=4);yt<-yt[yt>=lo & yt<=hi]
  for(v in yt) {line_mm(x+w,yy(v),x+w+1,yy(v));txt(format(v,trim=TRUE),x+w+1.8,yy(v))}
  txt("NES",x+w+8,(y+bottom)/2,hjust=.5,rot=270)
  for(i in seq_along(CAT_LEVELS)) key(CAT_LEVELS[i],CAT_COL[i],255,128+(i-1)*4.7)
  key("Enriched in HRD-",COL_NEG,255,166)
  key("Enriched in HRD+",COL_POS,255,171)
  key("Non-significant",COL_NS,255,176)
  if(any(d$missing)) txt("x  Not tested",255,181)
  invisible(d)
}

## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
write_volcano_table <- function(R) {
  if (is.null(R)) return(invisible(NULL))
  d <- R$deg
  d <- d[is.finite(d$log2fc) & is.finite(d$FDR), , drop = FALSE]
  if (!nrow(d)) stop("The volcano panel has no valid data to export.")
  fdr_plot <- pmax(d$FDR, FDR_FLOOR)
  sig <- as.character(d$sig)
  out <- data.frame(
    gene               = d$symbol,
    log2FC             = d$log2fc,
    P                  = d$P,
    FDR                = d$FDR,
    neg_log10_FDR      = -log10(fdr_plot),
    regulation         = sig,
    volcano_group      = unname(c(up   = "HRD+ associated",
                                  down = "HRD- associated",
                                  no   = "Not significant")[sig]),
    comparison         = "HRD+ vs HRD-",
    method             = if (isTRUE(R$adjust))
      "Linear model adjusted for ImmuneScore + StromalScore" else
        "Wilcoxon rank-sum test; log2FC = mean(HRD+) - mean(HRD-)",
    FDR_cutoff         = FDR_CUT,
    log2FC_cutoff      = LFC_CUT,
    n_HRDneg           = length(R$idx_neg),
    n_HRDpos           = length(R$idx_pos),
    stringsAsFactors   = FALSE)

  out <- out[order(match(out$regulation, c("up", "down", "no")), out$FDR), , drop = FALSE]

  if (sum(out$regulation == "up") != R$n_up || sum(out$regulation == "down") != R$n_down)
    stop("Volcano export counts disagree with the plotted up/down counts.")
  f <- file.path(WORKDIR, paste0(R$tag, "_Volcano_table.csv"))
  write.csv(out, f, row.names = FALSE)
  logmsg("Exported volcano source table: ", f, "(", nrow(out), " genes; HRD+-associated ",
         R$n_up, ", HRD--associated ", R$n_down, ")")
  invisible(f)
}

build_figure <- function(R, drop_E = FALSE) {
  if(is.null(R)) return(invisible(NULL))
  save_page<-function(file,draw) {
    grDevices::pdf(file=file,width=A4_W/25.4,height=A4_H/25.4,
                   family=FONT_FAMILY,pointsize=FONT_PT,paper="special",useDingbats=FALSE,
                   title="MIRROR transcriptomic features | A4 | 8 pt",colormodel="srgb")
    on.exit(grDevices::dev.off(),add=TRUE)
    grid::grid.newpage()
    draw()
  }
  f1<-file.path(WORKDIR,paste0(R$tag,"_1_Heatmap_KEGG_A4.pdf"))
  f2<-file.path(WORKDIR,paste0(R$tag,"_2_Volcano_Signatures_Hallmark_A4.pdf"))
  save_page(f1,function() {draw_deg_heat(R);draw_kegg(R)})
  save_page(f2,function() {draw_volcano(R);draw_signatures(R);draw_hallmark(R)})
  write_volcano_table(R)

  if(SAVE_A4_CACHE) saveRDS(list(R=R,PATH_LAB=PATH_LAB,LABEL_COL=LABEL_COL,
                                 FDR_CUT=FDR_CUT,LFC_CUT=LFC_CUT),
                            file.path(WORKDIR,paste0(R$tag,"_A4_plot_data.rds")),compress=FALSE)
  sample_columns<-intersect(c("sample_id","dataset","LABEL","group","mask_SeqRNA","mask_wsi"),names(R$clin))
  write.csv(R$clin[,sample_columns,drop=FALSE],
            file.path(WORKDIR,paste0(R$tag,"_plot_sample_order.csv")),row.names=FALSE)
  logmsg("Exported: ",f1)
  logmsg("Exported: ",f2)
  invisible(c(f1,f2))
}


replot_existing <- function() {
  read_result<-function(name,optional=FALSE) {
    f<-file.path(RESULT_DIR,name)
    if(!file.exists(f)) {if(optional) return(NULL); stop("Missing prior result: ",f)}
    read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
  }
  z<-NULL; scores<-NULL; clinical<-NULL
  tags<-paste0("Fig6_",RUN_SETS)
  if(REPLOT_ADJUSTED) tags<-c(tags,paste0("Fig6_",intersect(RUN_SETS,RUN_ADJUSTED_FOR),"_adjusted"))
  for(tag in tags) {
    cache<-file.path(RESULT_DIR,paste0(tag,"_A4_plot_data.rds"))
    if(file.exists(cache)) {
      saved<-readRDS(cache)
      if(!identical(saved$LABEL_COL,LABEL_COL)||saved$FDR_CUT!=FDR_CUT||saved$LFC_CUT!=LFC_CUT)
        stop("Cached labels or thresholds differ; rerun in analyse mode.")
      PATH_LAB<<-saved$PATH_LAB
      build_figure(saved$R);next
    }
    if(is.null(z)) {
      z<-read_matrix(F_EXPR_Z,"analysis z-score matrix")
      audit<-read_result("symbol_harmonization.csv",optional=!HARMONIZE_SYMBOLS)
      if(HARMONIZE_SYMBOLS && !is.null(audit)) {
        a<-audit[audit$status=="renamed",,drop=FALSE]
        names0<-rownames(z);m<-match(names0,a$old_symbol);ok<-!is.na(m)
        names0[ok]<-a$new_symbol[m[ok]];rownames(z)<-names0
        z<-merge_dup_rows(z,"z-score")
      }
      q<-read_result("ssGSEA_scores_all_samples.csv")
      scores<-as.matrix(q[,-1,drop=FALSE]);rownames(scores)<-q[[1]]
      storage.mode(scores)<-"numeric";colnames(scores)<-norm_id(colnames(scores))
      clinical<-read.csv(F_CLIN,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="GBK")
      clinical$sample_id<-norm_id(clinical$sample_id)
      clinical<-select_all_rna_samples(clinical)
      clinical<-clinical[!duplicated(clinical$sample_id),,drop=FALSE]
      clinical$LABEL<-clinical[[LABEL_COL]]
      if(LABEL_COL=="HRD_label" && "has_HRD_label" %in% names(clinical))
        clinical$LABEL[is.na(clinical$has_HRD_label)|clinical$has_HRD_label!=1]<-NA
      if(!all(as.character(na.omit(clinical$LABEL)) %in% c("0","1"))) stop("Grouping labels must be encoded as 0/1.")
      clinical$LABEL<-as.integer(as.character(clinical$LABEL))
    }
    cohort<-sub("_adjusted$","",sub("^Fig6_","",tag))
    cl<-clinical[clinical$dataset %in% COHORT_SETS[[cohort]] & !is.na(clinical$LABEL) &
                   clinical$sample_id %in% intersect(colnames(z),colnames(scores)),,drop=FALSE]
    cl$group<-factor(ifelse(cl$LABEL==1,"POS","NEG"),levels=c("NEG","POS"))
    cl<-cl[order(cl$group,cl$HRD_prob),,drop=FALSE]
    sm<-cl$sample_id
    if(!length(sm)) stop("No common samples are available for replotting.")
    deg<-read_result(paste0(tag,"_DEG_table.csv"));rownames(deg)<-deg$symbol
    expected<-ifelse(deg$FDR<=FDR_CUT & deg$log2fc>=LFC_CUT,"up",
                     ifelse(deg$FDR<=FDR_CUT & deg$log2fc<=-LFC_CUT,"down","no"))
    if(any(expected!=deg$sig,na.rm=TRUE)) stop("Prior DEG calls differ from the current thresholds; rerun analysis.")
    st<-read_result(paste0(tag,"_signature_stats.csv"),optional=grepl("_adjusted$",tag))
    if(!is.null(st)) {
      if(!all(st$signature %in% rownames(scores))) stop("ssGSEA scores lack required signatures.")

      dcheck<-rowMeans(scores[st$signature,sm[cl$LABEL==1],drop=FALSE])-
        rowMeans(scores[st$signature,sm[cl$LABEL==0],drop=FALSE])
      if(!isTRUE(all.equal(unname(dcheck),st$diff,tolerance=1e-8)))
        stop("Signature contrasts differ from saved results; inspect clinical IDs and labels, then rerun.")
    }
    R<-list(tag=tag,adjust=grepl("_adjusted$",tag),clin=cl,grp=cl$group,
            idx_neg=which(cl$LABEL==0),idx_pos=which(cl$LABEL==1),samples=sm,
            expr_z=z[,sm,drop=FALSE],deg=deg,n_up=sum(deg$sig=="up"),n_down=sum(deg$sig=="down"),
            ora_up=read_result(paste0(tag,"_ORA_HRDpos.csv"),TRUE),
            ora_down=read_result(paste0(tag,"_ORA_HRDneg.csv"),TRUE),
            gse=read_result(paste0(tag,"_Hallmark_GSEA.csv")),ss=scores[,sm,drop=FALSE],sig_stat=st)
    PATH_LAB<<-"KEGG pathway"
    logmsg("Replotting ",tag,": HRD- ",length(R$idx_neg)," / HRD+ ",length(R$idx_pos),
           "; differential genes ",R$n_down," / ",R$n_up,"; signature contrasts verified.")
    build_figure(R)
    RES[[cohort]]<<-R
  }
  invisible(RES)
}


RES <- list()
if (RUN_MODE == "replot") {
  replot_existing()
} else {
  ## ===========================================================================

  ## ===========================================================================
  logmsg("")
  logmsg("[1] Reading source matrices...")
  
  expr_rma <- read_matrix(F_EXPR,   "RMA expression matrix")
  expr_z   <- read_matrix(F_EXPR_Z, "z-score matrix")
  

  universe_before <- intersect(rownames(expr_rma), rownames(expr_z))
  
  if (HARMONIZE_SYMBOLS) {
    logmsg("")
    logmsg("[1b] Harmonizing legacy probe symbols to current HGNC names...")
    

    SYMBOL_MAP <- build_symbol_map(union(rownames(expr_rma), rownames(expr_z)))
    ren   <- attr(SYMBOL_MAP, "renamed")
    audit <- attr(SYMBOL_MAP, "audit")
    
    if (is.null(ren) || !nrow(ren)) {
      logmsg("    No gene symbols require renaming.")
    } else {
      logmsg("    ", nrow(ren), " legacy symbols will be mapped to current HGNC symbols.")
      if (!is.null(audit)) {
        rej <- audit[audit$status != "renamed", , drop = FALSE]
        if (nrow(rej)) {
          logmsg("    Additionally, ", nrow(rej), " candidate mappings were rejected by safety criteria:")
          logtable(as.data.frame(table(rej$status)))
          rej_d <- rej[grepl("claimed by", rej$status), , drop = FALSE]
          if (nrow(rej_d))
            logmsg("    Targets claimed by multiple aliases (formerly merged silently): ",
                   paste(unique(rej_d$new_symbol), collapse = ", "))
        }
        write.csv(audit, "symbol_harmonization.csv", row.names = FALSE)
        logmsg("    Mapping and rejection audit exported to symbol_harmonization.csv.")
      }
      expr_rma <- apply_symbol_map(expr_rma, SYMBOL_MAP, "RMA expression matrix")
      expr_z   <- apply_symbol_map(expr_z,   SYMBOL_MAP, "z-score matrix")
    }
    report_example_renames(SYMBOL_MAP)
  }
  
  clin <- read.csv(F_CLIN, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "GBK")
  for (cc in c("sample_id", "mask_SeqRNA", "dataset", LABEL_COL, "HRD_prob"))
    if (!cc %in% names(clin)) stop("Clinical table is missing a required field: ", cc)
  clin$sample_id <- norm_id(clin$sample_id)

  clin$LABEL <- clin[[LABEL_COL]]

  if (LABEL_COL == "HRD_label" && "has_HRD_label" %in% names(clin)) {
    no_lab <- is.na(clin$has_HRD_label) | clin$has_HRD_label != 1
    clin$LABEL[no_lab] <- NA
  }
  u_lab <- unique(as.character(clin$LABEL[!is.na(clin$LABEL)]))
  if (!length(u_lab)) stop("Label field ", LABEL_COL, " is entirely missing; no eligible samples.")
  if (!all(u_lab %in% c("0", "1")))
    stop("Label field ", LABEL_COL, " must contain 0/1; observed: ",
         paste(u_lab, collapse = ", "))
  clin$LABEL <- as.integer(as.character(clin$LABEL))
  logmsg("    Grouping field: ", LABEL_COL, " (labelled ",
         sum(!is.na(clin$LABEL)), " / ", nrow(clin), " samples)")
  
  clin_rna <- select_all_rna_samples(clin)


  dup_cl <- duplicated(clin_rna$sample_id)
  if (any(dup_cl)) {
    logmsg("    Duplicate clinical sample_id values: ", sum(dup_cl), "; retaining the first record for each: ",
           paste(head(unique(clin_rna$sample_id[dup_cl]), 10), collapse = ", "))
    clin_rna <- clin_rna[!dup_cl, , drop = FALSE]
  }
  
  common <- Reduce(intersect, list(clin_rna$sample_id,
                                   colnames(expr_rma), colnames(expr_z)))
  logmsg(sprintf("    RNA clinical %d; RMA %d; z-score %d; common samples %d",
                 nrow(clin_rna), ncol(expr_rma), ncol(expr_z), length(common)))
  if (length(common) < 20) stop("Sample identifiers do not match across the input tables.")
  
  clin <- clin_rna[match(common, clin_rna$sample_id), ]
  
  gene_common <- intersect(rownames(expr_rma), rownames(expr_z))
  expr_rma <- expr_rma[gene_common, clin$sample_id, drop = FALSE]
  expr_z   <- expr_z[  gene_common, clin$sample_id, drop = FALSE]
  universe <- rownames(expr_rma)
  logmsg("    Shared expression genes: ", length(universe), " (enrichment background); ",
         "before symbol harmonization: ", length(universe_before), " genes")
  

  audit_cols <- intersect(c("sample_id", "cohort", "dataset", "LABEL", "has_HRD_label",
                            "mask_SeqRNA", "mask_wsi"), names(clin))
  sample_audit <- clin[, audit_cols, drop=FALSE]
  sample_audit$in_requested_sets <- clin$dataset %in% unique(unlist(COHORT_SETS[RUN_SETS]))
  sample_audit$used_in_comparison <- sample_audit$in_requested_sets & !is.na(clin$LABEL)
  sample_audit$reason <- ifelse(!sample_audit$in_requested_sets, "outside_selected_dataset",
                                ifelse(is.na(clin$LABEL), "missing_HRD_label", "included"))
  write.csv(sample_audit, "All_RNA_sample_audit.csv", row.names=FALSE)
  included <- sample_audit[sample_audit$used_in_comparison, , drop=FALSE]
  write.csv(included, "All_RNA_analysis_sample_list.csv", row.names=FALSE)
  logmsg("    All-RNA common samples: ", nrow(sample_audit), "; included in primary analysis: ", nrow(included),
         " (HRD- ", sum(included$LABEL==0), " / HRD+ ", sum(included$LABEL==1), ").")
  logmsg("    WSI pairing is not required; missing reference HRD labels are not imputed.")
  
  n_nolab <- sum(is.na(clin$LABEL))
  if (n_nolab)
    logmsg("    Common samples with missing labels: ", n_nolab, " missing ", LABEL_COL,
           "; excluded from differential testing.")
  tb <- table(clin$dataset, factor(clin$LABEL, levels = c(0, 1)))
  logmsg("")
  logmsg("    Sample counts by dataset (rows) and HRD 0/1 (columns); label field = ",
         LABEL_COL, "):")
  logtable(tb, row.names = TRUE)
  
  for (nm in RUN_SETS) {
    if (is.null(COHORT_SETS[[nm]]))
      stop("Unknown RUN_SETS entry '", nm, "' in COHORT_SETS.")
    logmsg(sprintf("    Cohort %-9s = {%s}: %d samples (%s labels in %d)",
                   nm, paste(COHORT_SETS[[nm]], collapse = ", "),
                   sum(clin$dataset %in% COHORT_SETS[[nm]]), LABEL_COL,
                   sum(clin$dataset %in% COHORT_SETS[[nm]] & !is.na(clin$LABEL))))
  }
  
  clin$truth <- NA_integer_
  if (all(c("HRD_label", "has_HRD_label") %in% names(clin))) {
    ok <- clin$has_HRD_label == 1 & !is.na(clin$HRD_label)
    clin$truth[ok] <- as.integer(clin$HRD_label[ok])
  } else if ("HRD_label" %in% names(clin)) {
    clin$truth <- suppressWarnings(as.integer(clin$HRD_label))
  }

  ## ===========================================================================

  ## ===========================================================================
  logmsg("")
  logmsg("[2] Loading functional signatures and evaluating coverage...")
  
  if (nzchar(F_FGES) && file.exists(F_FGES)) {
    fg  <- read_fges_gmt(F_FGES)
    raw <- fg$sets
    logmsg("    Read official Fges GMT: ", length(raw), " signatures")
    

    cat_raw <- FGES_CATMAP[names(raw)]
    if (any(is.na(cat_raw))) {
      bad <- names(raw)[is.na(cat_raw)]
      warning("Unrecognized Fges IDs will be categorized as Tumor state; ",
              "check the GMT version: ", paste(bad, collapse = ", "))
      logmsg("    Unrecognized signature IDs: ", paste(bad, collapse = ", "))
      cat_raw[is.na(cat_raw)] <- "Tumor state"
    }
    

    cov <- coverage_table(raw, universe, "Fges", universe_before)
    cov$category <- as.character(cat_raw[cov$id])
    cov$label    <- if (FGES_LABEL_COL >= 2) fg$label[cov$id] else cov$id
    cov$kept     <- cov$pct >= FGES_MIN_COVERAGE * 100 & cov$n_hit >= FGES_MIN_GENES
    write.csv(cov, "Fges_coverage.csv", row.names = FALSE)
    logmsg("    Signature coverage exported to Fges_coverage.csv.")
    logmsg(sprintf("    Mean Fges coverage rose from %.1f%% to %.1f%% after harmonization.",
                   mean(cov$pct_before), mean(cov$pct)))
    logmsg("    Five signatures with lowest coverage:")
    logtable(head(cov[, c("id", "n_total", "pct_before", "pct", "kept")], 5))
    
    if (any(!cov$kept))
      logmsg("    Signatures below minimum coverage ", round(FGES_MIN_COVERAGE * 100),
             "%; excluded: ", paste(cov$id[!cov$kept], collapse = ", "))
    
    keep_id  <- cov$id[cov$kept]
    SIG_SETS <- lapply(raw[keep_id], function(g) intersect(g, universe))
    nm_show  <- if (FGES_LABEL_COL >= 2) unname(fg$label[keep_id]) else keep_id
    names(SIG_SETS) <- nm_show
    SIG_CAT  <- setNames(as.character(cat_raw[keep_id]), nm_show)
    FGES_SRC <- "Bagaev"
    logmsg("    Final signature count: ", length(SIG_SETS), " for ssGSEA.")
  } else {
    SIG_SETS <- lapply(BUILTIN_SIG, function(g) intersect(toupper(g), universe))
    SIG_SETS <- SIG_SETS[vapply(SIG_SETS, length, 1L) >= FGES_MIN_GENES]
    SIG_CAT  <- BUILTIN_CAT
    FGES_SRC <- "builtin"
    logmsg("    F_FGES is absent; using provisional built-in marker panels (",
           length(SIG_SETS), " sets). Replace these with the official Bagaev GMT before publication.")
  }
  if (!length(SIG_SETS)) stop("No signatures passed coverage filtering; inspect gene identifiers.")
  
  logmsg("")
  logmsg("[3] Calculating ssGSEA scores once for all samples...")
  ss <- tryCatch(
    GSVA::gsva(GSVA::ssgseaParam(expr_rma, SIG_SETS, normalize = TRUE), verbose = FALSE),
    error = function(e)
      GSVA::gsva(expr_rma, SIG_SETS, method = "ssgsea",
                 ssgsea.norm = TRUE, verbose = FALSE))
  write.csv(ss, "ssGSEA_scores_all_samples.csv")
  logmsg("    ssGSEA complete: ", nrow(ss), " signatures x ", ncol(ss), " samples.")
  
  ssz <- t(scale(t(ss)))
  imm_sig <- intersect(names(SIG_CAT)[SIG_CAT %in% c("Anti-tumor immune", "Pro-tumor immune")],
                       rownames(ssz))
  str_sig <- intersect(names(SIG_CAT)[SIG_CAT == "Angiogenesis / fibroblasts"],
                       rownames(ssz))
  logmsg("    Immune score averages ", length(imm_sig), " signatures; stromal score averages ",
         length(str_sig), " sets")
  clin$ImmuneScore  <- as.numeric(colMeans(ssz[imm_sig, , drop = FALSE])[clin$sample_id])
  clin$StromalScore <- as.numeric(colMeans(ssz[str_sig, , drop = FALSE])[clin$sample_id])
  write.csv(clin[, c("sample_id", "dataset", "ImmuneScore", "StromalScore")],
            "TME_scores.csv", row.names = FALSE)
  
  
  ## ===========================================================================

  ## ===========================================================================
  logmsg("")
  logmsg("[4] Loading pathway gene sets...")
  pw <- get_pathway_t2g()
  if (is.null(pw$t2g))
    stop("Gene sets are unavailable. Install msigdbr; newer releases may need ",
         "install.packages('msigdbdf', repos='https://igordot.r-universe.dev')")
  path_t2g <- pw$t2g; PATH_LAB <- pw$label
  hall_t2g <- load_gmt(F_HALLMARK, "Hallmark")
  
  cov_all <- rbind(
    coverage_table(split(hall_t2g$gene, hall_t2g$term), universe, "Hallmark", universe_before),
    coverage_table(split(path_t2g$gene, path_t2g$term), universe, PATH_LAB,  universe_before),
    coverage_table(DDR_SETS, universe, "DDR", universe_before))
  write.csv(cov_all, "geneset_coverage.csv", row.names = FALSE)
  logmsg(sprintf("    Mean Hallmark coverage %.1f%% (prior %.1f%%); mean %s coverage %.1f%% (prior %.1f%%)",
                 mean(cov_all$pct[cov_all$source == "Hallmark"]),
                 mean(cov_all$pct_before[cov_all$source == "Hallmark"]),
                 PATH_LAB,
                 mean(cov_all$pct[cov_all$source == PATH_LAB]),
                 mean(cov_all$pct_before[cov_all$source == PATH_LAB])))
  

  ddr_cov <- cov_all[cov_all$source == "DDR", , drop = FALSE]
  ddr_cov$kept <- ddr_cov$pct >= DDR_MIN_COVERAGE * 100 & ddr_cov$n_hit >= DDR_MIN_GENES
  logmsg("    DDR-specific pathway coverage:")
  ddr_show <- ddr_cov[, c("id", "n_total", "n_hit", "pct", "missing", "kept")]
  ddr_show$missing <- trunc_lab(ddr_show$missing, 45)
  logtable(ddr_show)
  if (any(!ddr_cov$kept)) {
    logmsg("    DDR pathways below coverage threshold ", round(DDR_MIN_COVERAGE * 100),
           "%; excluded because null findings would be uninterpretable: ",
           paste(ddr_cov$id[!ddr_cov$kept], collapse = ", "))
    logmsg("       Missing core genes: ",
           paste(ddr_cov$missing[!ddr_cov$kept], collapse = " | "), ")")
  }
  DDR_SETS_KEPT <- DDR_SETS[ddr_cov$id[ddr_cov$kept]]
  if (!length(DDR_SETS_KEPT))
    logmsg("    No DDR pathways passed coverage filtering; DDR GSEA omitted.")

  ## ===========================================================================

  ## ===========================================================================
  
  RES <- list()
  
  for (nm in RUN_SETS) {
    ds  <- COHORT_SETS[[nm]]

    sel <- clin$dataset %in% ds & !is.na(clin$LABEL)
    if (!any(sel)) {
      logmsg("")
      logmsg("Cohort ", nm, " has no samples with ", LABEL_COL, " labels; skipped.")
      next
    }
    
    R <- run_analysis(paste0("Fig6_", nm),
                      clin$sample_id[sel], clin$LABEL[sel], adjust = FALSE)
    build_figure(R)
    if (!is.null(R)) RES[[nm]] <- R
    
    if (nm %in% RUN_ADJUSTED_FOR) {
      Ra <- run_analysis(paste0("Fig6_", nm, "_adjusted"),
                         clin$sample_id[sel], clin$LABEL[sel], adjust = TRUE)
      build_figure(Ra, drop_E = ADJ_DROP_E)
    }
  }
  
  
  ## ===========================================================================

  ## ===========================================================================
  if (COMPARE_COHORTS && length(RES) >= 2) {
    logmsg("")
    logmsg("===== Between-cohort comparisons =====")

    logmsg("    Train, Val and Internal may be splits of the same TCGA-OV cohort; ",
           "these comparisons measure split-sample consistency rather than external replication; ",
           "state this explicitly in the methods.")
    nms <- names(RES)
    
    summ <- do.call(rbind, lapply(nms, function(k) {
      R <- RES[[k]]
      data.frame(cohort = k, n_total = length(R$samples),
                 n_HRDneg = length(R$idx_neg), n_HRDpos = length(R$idx_pos),
                 n_up = R$n_up, n_down = R$n_down,
                 median_log2fc = round(median(R$deg$log2fc, na.rm = TRUE), 4),
                 pct_positive = round(100 * mean(R$deg$log2fc > 0, na.rm = TRUE), 1),
                 n_FDR05 = sum(R$deg$FDR < 0.05, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
    write.csv(summ, "cohort_summary.csv", row.names = FALSE)
    logmsg("    Cohort overview:")
    logtable(summ)
    
    gsh <- Reduce(intersect, lapply(RES, function(R) R$deg$symbol))
    L <- sapply(RES, function(R) R$deg[match(gsh, R$deg$symbol), "log2fc"])
    rownames(L) <- gsh
    
    cmb <- t(combn(nms, 2))
    cmp <- do.call(rbind, lapply(seq_len(nrow(cmb)), function(i) {
      a <- cmb[i, 1]; b <- cmb[i, 2]
      up_a <- RES[[a]]$deg$symbol[RES[[a]]$deg$sig == "up"]
      up_b <- RES[[b]]$deg$symbol[RES[[b]]$deg$sig == "up"]
      dn_a <- RES[[a]]$deg$symbol[RES[[a]]$deg$sig == "down"]
      dn_b <- RES[[b]]$deg$symbol[RES[[b]]$deg$sig == "down"]
      jac <- function(x, y) if (!length(union(x, y))) NA else
        round(length(intersect(x, y)) / length(union(x, y)), 3)
      data.frame(cohort_A = a, cohort_B = b,
                 pearson_log2fc  = round(cor(L[, a], L[, b], use = "complete.obs"), 3),
                 spearman_log2fc = round(cor(L[, a], L[, b], method = "spearman",
                                             use = "complete.obs"), 3),
                 jaccard_up = jac(up_a, up_b), jaccard_down = jac(dn_a, dn_b),
                 stringsAsFactors = FALSE)
    }))
    write.csv(cmp, "cohort_comparison.csv", row.names = FALSE)
    logmsg("    Pairwise comparison (gene log2FC correlation and DEG overlap):")
    logtable(cmp)
    
    M <- cor(L, use = "complete.obs")
    mdf <- expand.grid(A = nms, B = nms, stringsAsFactors = FALSE)
    mdf$r <- as.vector(M[cbind(match(mdf$A, rownames(M)), match(mdf$B, colnames(M)))])
    pC2 <- ggplot(mdf, aes(A, B, fill = r)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.2f", r)), size = 2.2, colour = COL_TEXT) +
      scale_fill_gradientn(colours = HEAT_COL, limits = c(-1, 1),
                           breaks = c(-1, 0, 1), name = "Pearson r") +
      labs(x = NULL, y = NULL,
           title = expression("Cross-cohort correlation of gene-level log"[2]*"(FC)")) +
      theme_mirror(7) + theme(legend.position = "right")
    tryCatch(ggsave("cohort_comparison.pdf", pC2, width = 110, height = 90,
                    units = "mm", device = cairo_pdf),
             error = function(e) ggsave("cohort_comparison.pdf", pC2, width = 110,
                                        height = 90, units = "mm"))
    logmsg("    Saved cohort comparison plots and tables.")
  }
  
  
  ## ===========================================================================

  ## ===========================================================================
  if (RUN_TRUTH && any(!is.na(clin$truth))) {
    sub <- !is.na(clin$truth)
    Rt <- run_analysis("Fig6_truth", clin$sample_id[sub], clin$truth[sub], adjust = FALSE)
    build_figure(Rt)
  }
  
  if (RUN_CONCORDANCE && any(!is.na(clin$truth)) &&
      "Multi_PRE_HRD" %in% names(clin)) {
    logmsg("")
    logmsg("===== Concordance: model predictions versus genomic reference labels =====")
    sub <- which(!is.na(clin$truth))
    sm  <- clin$sample_id[sub]
    E   <- expr_rma[, sm, drop = FALSE]
    lfc_of <- function(g01) {
      ip <- which(g01 == 1); ineg <- which(g01 == 0)
      if (length(ip) < 10 || length(ineg) < 10) return(NULL)
      rowMeans(E[, ip, drop = FALSE]) - rowMeans(E[, ineg, drop = FALSE])
    }
    l_pred <- lfc_of(clin$Multi_PRE_HRD[sub]); l_true <- lfc_of(clin$truth[sub])
    if (!is.null(l_pred) && !is.null(l_true)) {
      r_p <- cor(l_pred, l_true, use = "complete.obs")
      r_s <- cor(l_pred, l_true, method = "spearman", use = "complete.obs")
      acc <- mean(clin$Multi_PRE_HRD[sub] == clin$truth[sub], na.rm = TRUE)
      logmsg(sprintf("    Paired n = %d; label agreement = %.1f%%; Pearson log2FC r = %.3f; Spearman rho = %.3f",
                     length(sm), 100 * acc, r_p, r_s))
      cdf <- data.frame(pred = as.numeric(l_pred), truth = as.numeric(l_true),
                        gene = names(l_pred), stringsAsFactors = FALSE)
      write.csv(cdf, "FigS_concordance_lfc.csv", row.names = FALSE)
      lim <- max(abs(c(cdf$pred, cdf$truth)), na.rm = TRUE) * 1.05
      pS <- ggplot(cdf, aes(truth, pred)) +
        geom_point(size = 0.25, alpha = 0.35, colour = COL_NS, shape = 16) +
        geom_point(data = subset(cdf, abs(pred) >= LFC_CUT & abs(truth) >= LFC_CUT),
                   aes(colour = pred > 0), size = 0.55, alpha = 0.8, shape = 16) +
        scale_colour_manual(values = c(`TRUE` = COL_POS, `FALSE` = COL_NEG), guide = "none") +
        geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25, colour = "grey40") +
        geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
        annotate("text", x = -lim * 0.95, y = lim * 0.92, hjust = 0, size = 2.1,
                 colour = COL_TEXT,
                 label = sprintf("Pearson r = %.3f\nSpearman rho = %.3f\nn = %d samples",
                                 r_p, r_s, length(sm))) +
        coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
        labs(x = expression("log"[2]*"(FC), genomic HRD label"),
             y = expression("log"[2]*"(FC), MIRROR-predicted HRD")) +
        theme_mirror(7)
      tryCatch(ggsave("FigS_concordance.pdf", pS, width = 90, height = 90,
                      units = "mm", device = cairo_pdf),
               error = function(e) ggsave("FigS_concordance.pdf", pS, width = 90,
                                          height = 90, units = "mm"))
    }
  }

}
writeLines(capture.output(sessionInfo()), file.path(WORKDIR, "A4_sessionInfo.txt"))
logmsg("Complete: A4 landscape 297 x 210 mm with 8-pt typography. Output:",WORKDIR)
}
run_mirror_transcriptomics()
