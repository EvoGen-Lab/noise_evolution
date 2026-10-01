###############################################################
### functions_external.R
### Functions for the comparison with published noise data
### (Section 9), sourced by external_validation.R after
### functions.R, whose fit_counts_offset() fits each count-matrix
### source.
###
### Validation of the MIX.SC NB fit (mean, CV^2, burst frequency, burst size)
### against published yeast noise data. Summary-statistic sources (Newman, Keren,
### Stewart-Ornstein) enter through add_burst_terms(); count-matrix sources
### (Gasch, Nadal-Ribelles, Jackson, Jariani) are read, QC'd, mapped to ORFs and
### fitted with the same offset NB model. All sources then share one
### (ORF, Mean, CV2, Fano, BFREQ, BSIZE) layout, so cor_row() compares any pair.
###
### Outline (function name - purpose):
###     cor_row() - Spearman rank correlation between two vectors, with pairwise-complete filtering; returns one row (n, rho, p).
###     add_burst_terms() - Implied Fano factor, burst size, and burst frequency algebraically recovered from a reported mean and CV^2.
###     read_header_line() - Header line of a delimited file as a character vector, with quotes stripped.
###     get_data_header() - Sample/gene names of a file's header with the row-id field removed, whichever header convention the file uses.
###     fread_matrix() - Fast row-named numeric matrix reader (fread based) with column selection.
###     to_numeric_matrix() - Coerces a table to a numeric matrix, reporting and dropping non-numeric columns.
###     map_to_orf() - Resolves gene identifiers to systematic ORF names (NA when unresolved).
###     fit_scrna_source() - One published scRNA-seq source: ORF collapse, cell/gene QC and offset NB fit, returned on the NB.SC layout (ORF, Mean, CV2, Fano, BFREQ, BSIZE).
###     mean_adjusted_noise() - Loess residual of log CV^2 on log mean: noise relative to genes of similar abundance.
###     noise_corr_rows() - Spearman correlation of mean, CV2, burst frequency and burst size between a source and MIX.SC or between two sources.
###     raw_source_table() - the ORF, Mean, CV2, BFREQ and BSIZE columns of one fitted source, in the layout shared.
###     add_cv2_adj() - adds CV2_ADJ, the loess residual of log CV2 on log mean (mean_adjusted_noise()), to one.
###     adj_corr_row() - agreement of source src's mean-adjusted noise (CV2_ADJ) with MIX.SC's (nb_sc) over the.
###############################################################

## Spearman rank correlation between two vectors on their pairwise-finite
## observations; returns a one-row data.frame (n, rho, p). Rank correlation
## suits the external validation because both the NB-fitted and the external
## abundance and noise statistics are heavy-tailed, and it makes no
## assumption about the shape of that tail. Needs at least 3 finite pairs.
cor_row <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  # exact = FALSE selects the asymptotic t approximation for the Spearman
  # p-value, which remains valid when the ranks contain ties.
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(n = sum(keep), rho = unname(ct$estimate), p = ct$p.value)
}

## Implied Fano factor, burst size, and burst frequency recovered from a
## reported mean and CV^2, for a data frame with columns Mean and CV2. The NB
## identity used for MIX.SC runs in reverse: Fano = mean * CV^2, burst size =
## Fano - 1, burst frequency = mean / (Fano - 1). This is an algebraic
## transform of the reported summary statistics, not a re-fit. Fano <= 1
## (sub-Poissonian) is read as measurement noise and gives NA burst terms.
## Interpretive caveats: burst size and frequency are mRNA-level kinetic
## constructs while the protein-level sources add translational and
## degradation noise; and for arbitrary-unit fluorescence the product
## mean * CV^2 carries the instrument's units, so the Fano > 1 gate and the
## burst values are relative quantities. The Fano column is overwritten with
## mean * CV^2 so every source uses one definition.
add_burst_terms <- function(d) {
  fano <- d$Mean * d$CV2
  ok <- fano > 1
  d$Fano  <- fano
  d$BSIZE <- ifelse(ok, fano - 1, NA_real_)
  d$BFREQ <- ifelse(ok, d$Mean / (fano - 1), NA_real_)
  d
}

## read_header_line() parses only the header line of a delimited file, so the
## header's own field count is known independently of fread()'s header
## detection. The scRNA-seq sources write the header one field short of the
## data rows (the row-id column is unnamed), a layout fread() cannot
## reconcile on its own. Surrounding quotes are stripped from each field.
read_header_line <- function(path, sep = "\t", skip = 0) {
  lines  <- readLines(path, n = skip + 1)
  fields <- strsplit(lines[skip + 1], sep, fixed = TRUE)[[1]]
  gsub('^"|"$', "", fields)
}

## get_data_header() returns the sample/gene names of a file's header with the
## row-id field removed. Two header conventions occur in the scRNA-seq
## sources: the header omits the row-id field (Jackson), so it is one field
## shorter than a data row; or it labels the field, explicitly (Gasch) or with
## a blank (Jariani), so header and data rows have equal length. Comparing the
## header's field count with the first data row's tells them apart, so
## fread_matrix() always receives names for the data columns only.
get_data_header <- function(path, sep = "\t", skip = 0) {
  header.raw <- read_header_line(path, sep = sep, skip = skip)
  data.line  <- readLines(path, n = skip + 2)[skip + 2]
  n.data     <- length(strsplit(data.line, sep, fixed = TRUE)[[1]])
  if (length(header.raw) == n.data - 1) {
    header.raw            # header already omits the row-id field
  } else if (length(header.raw) == n.data) {
    header.raw[-1]         # header carries a label (or blank) for the row-id field; drop it
  } else {
    stop("get_data_header(): header has ", length(header.raw),
         " field(s) but the first data row has ", n.data,
         " field(s) in ", path, "; this file matches neither the ",
         "'header omits row-id' nor the 'header labels row-id' convention")
  }
}

## fread_matrix() is a fast equivalent of as.matrix(read.delim(path, row.names
## = 1)) built on data.table::fread(). keep is a logical vector over
## get_data_header()'s columns selecting which to parse (NULL keeps all), so
## only the columns a source needs are read. Data rows are read with header =
## FALSE; each row's first field becomes the row name and column names come
## from get_data_header(), which makes the +1 column shift below correct for
## every source.
fread_matrix <- function(path, keep = NULL, sep = "\t", skip = 0, ...) {
  header <- get_data_header(path, sep = sep, skip = skip)
  if (is.null(keep)) keep <- rep(TRUE, length(header))
  col.idx <- c(1L, which(keep) + 1L)   # +1 shifts past the row-id field, which the header never names
  dt  <- fread(file = path, sep = sep, skip = skip + 1, header = FALSE, select = col.idx, ...)
  mat <- as.matrix(dt[, -1, with = FALSE])
  colnames(mat) <- header[keep]
  rownames(mat) <- dt[[1]]
  mat
}

## to_numeric_matrix() coerces a table to a numeric matrix. One non-numeric
## column (a metadata field, a stray total-count column, a blank trailing
## column) would otherwise turn the whole as.matrix() result into character.
## Columns that do not convert cleanly are reported by name and dropped, so
## rowsum() and the NB fit receive numeric input.
to_numeric_matrix <- function(df, label) {
  mat <- as.matrix(df)
  if (is.character(mat)) {
    num <- suppressWarnings(apply(mat, 2, as.numeric))
    # a column is genuinely non-numeric if it produced an NA somewhere
    # that was not already a blank entry in the original character data
    bad <- colSums(is.na(num) & mat != "", na.rm = TRUE) > 0
    if (any(bad))
      message(label, ": dropping ", sum(bad), " non-numeric column(s): ",
              paste(head(colnames(mat)[bad], 10), collapse = ", "))
    mat <- num[, !bad, drop = FALSE]
    rownames(mat) <- rownames(df)
  }
  storage.mode(mat) <- "numeric"
  mat
}

map_to_orf <- function(ids, label = "") {
  is.orf <- grepl(ORF.PATTERN, ids)
  out <- ids
  need <- unique(ids[!is.orf])
  if (length(need) > 0) {
    # mapIds() throws when NONE of the supplied keys exist under that
    # keytype, rather than just returning all-NA the way it does for a
    # partial mismatch. That total-mismatch case almost always means the
    # identifiers in "need" are not actually common names at all (most
    # often because the systematic-name regex above failed to match the
    # real format this source uses), so it is reported rather than left
    # to crash the whole script, with a sample of the unmatched strings
    # printed so the real format can be diagnosed directly.
    looked.up <- tryCatch(
      suppressMessages(mapIds(org.Sc.sgd.db, keys = need, keytype = "GENENAME",
                               column = "ORF", multiVals = "first")),
      error = function(e) {
        message(label, ": ", length(need), " identifier(s) did not match the systematic-ORF ",
                "pattern and also failed as GENENAME lookups. Sample: ",
                paste(head(need, 10), collapse = ", "))
        setNames(rep(NA_character_, length(need)), need)
      })
    out[!is.orf] <- unname(looked.up[match(ids[!is.orf], need)])
  }
  n.unresolved <- sum(is.na(out))
  if (n.unresolved > 0)
    message(label, ": ", n.unresolved, " of ", length(out), " row(s) did not resolve to an ORF and will be dropped")
  out
}

## fit_scrna_source: takes the raw genes-by-cells count matrix of one published single-cell source to a
## data.frame keyed by ORF, with the Mean/CV2/Fano/BFREQ/BSIZE content NB.SC carries, so every source
## merges against NB.SC the same way.
##   1. Row names raw_ids resolve to systematic ORFs (map_to_orf); rows that never map are dropped and rows
##      sharing an ORF are summed (TSS- or isoform-level rows collapse to one gene-level row).
##   2. QC drops cells whose library size is below min_cell_count, then genes detected in fewer than
##      min_cells_expr of the remaining cells. The thresholds are arguments because the sources differ
##      widely in depth and cell number.
##   3. Exposure is each cell's library size divided by the source's mean library size, so MU is in
##      counts per typical cell (the unit MIX.SC uses) and CV2 = 1/MU + 1/DISP is at this dataset's own
##      depth. BFREQ = DISP and BSIZE = MU/DISP, the definitions documented for NB.SC.
## cl is optional and passed to fit_counts_offset(): NULL fits on one core; a PSOCK cluster from
## parallel::makeCluster() splits genes across workers (the driver reuses one cluster for all sources).
fit_scrna_source <- function(mat, raw_ids, label, min_cell_count, min_cells_expr, cl = NULL) {
  stopifnot(is.numeric(min_cell_count), length(min_cell_count) == 1, min_cell_count >= 0,
            is.numeric(min_cells_expr), length(min_cells_expr) == 1, min_cells_expr >= 0)
  orf <- map_to_orf(raw_ids, label)
  mat <- mat[!is.na(orf), , drop = FALSE]
  orf <- orf[!is.na(orf)]
  mat <- rowsum(mat, group = orf)

  mat <- mat[, colSums(mat) >= min_cell_count, drop = FALSE]
  mat <- mat[rowSums(mat > 0) >= min_cells_expr, , drop = FALSE]

  exposure <- colSums(mat) / mean(colSums(mat))
  fit <- fit_counts_offset(mat, exposure, cl = cl)
  data.frame(
    ORF   = rownames(fit),
    Mean  = fit$MU,
    CV2   = fit$CV^2,
    Fano  = fit$FANO,
    BFREQ = fit$BFREQ,
    BSIZE = fit$BSIZE,
    stringsAsFactors = FALSE
  )
}

## mean_adjusted_noise() expresses each gene's noise relative to genes of
## similar abundance in the same dataset, the idea behind Newman et al.'s DM
## statistic: the residual of log CV^2 from a loess trend on log mean. The
## trend is fitted within each source, so each dataset keeps its own capture
## efficiency and units and the residual isolates gene-specific noise that can
## be compared across datasets. Genes with non-positive or non-finite mean or
## CV^2 return NA.
mean_adjusted_noise <- function(mean, cv2, span = 0.3) {
  lm_ <- log(mean)
  lc  <- log(cv2)
  ok  <- is.finite(lm_) & is.finite(lc)
  res <- rep(NA_real_, length(mean))
  fit <- loess(lc[ok] ~ lm_[ok], span = span)   # smooth abundance trend, as in Newman's DM
  res[ok] <- lc[ok] - predict(fit)              # residual = noise beyond the expected level at that abundance
  res
}

## noise_corr_rows: Spearman correlation of the mean, CV2, burst frequency and burst size between two noise
## tables, one row per statistic. kind "vs_mix": key is a source name and tables[[key]] is that source merged
## with NB.SC on ORF (columns MU/Mean, CV2.x/CV2.y, ...), compared against MIX.SC. kind "pairwise": key is a
## pair of source names and the two tables in tables are merged on ORF. The first column names the source
## or the comparison.
noise_corr_rows <- function(key, tables, kind = c("vs_mix", "pairwise")) {
  kind <- match.arg(kind)
  stats <- c("mean", "CV2", "burst_frequency", "burst_size")
  if (kind == "vs_mix") {
    d <- tables[[key]]
    x <- c("MU", "CV2.x", "BFREQ.x", "BSIZE.x"); y <- c("Mean", "CV2.y", "BFREQ.y", "BSIZE.y")
  } else {
    d <- merge(tables[[key[1]]], tables[[key[2]]], by = "ORF")
    x <- paste0(c("Mean", "CV2", "BFREQ", "BSIZE"), ".x"); y <- paste0(c("Mean", "CV2", "BFREQ", "BSIZE"), ".y")
  }
  out <- do.call(rbind, unname(Map(cor_row, lapply(d[x], log), lapply(d[y], log))))
  if (kind == "vs_mix") out$source <- key else out$comparison <- paste(key[1], "vs", key[2])
  out$statistic <- stats
  out
}

## raw_source_table: the ORF, Mean, CV2, BFREQ and BSIZE columns of one fitted source, in the layout shared
## by all published sources.
raw_source_table <- function(d) {
  data.frame(ORF = d$ORF, Mean = d$Mean, CV2 = d$CV2, BFREQ = d$BFREQ, BSIZE = d$BSIZE,
             stringsAsFactors = FALSE)
}

## add_cv2_adj: adds CV2_ADJ, the loess residual of log CV2 on log mean (mean_adjusted_noise()), to one
## source table d.
add_cv2_adj <- function(d) {
  d$CV2_ADJ <- mean_adjusted_noise(d$Mean, d$CV2)
  d
}

## adj_corr_row: agreement of source src's mean-adjusted noise (CV2_ADJ) with MIX.SC's (nb_sc) over the
## genes both report. The adjusted values are residuals on the log scale already, so no further transform.
adj_corr_row <- function(src, all_raw, nb_sc) {
  d <- merge(nb_sc[, c("ORF", "CV2_ADJ")], all_raw[[src]][, c("ORF", "CV2_ADJ")], by = "ORF")
  out <- cor_row(d$CV2_ADJ.x, d$CV2_ADJ.y)
  out$source <- src
  out
}
