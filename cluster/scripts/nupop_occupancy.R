###############################################################
### nupop_occupancy.R
### Reads the NuPoP inputs, predicts nucleosome occupancy for
### every chromosome that carries a promoter in both species, and
### saves each gene's promoter occupancy track. Mirrors the
### structure of gene_boot.R.
###
### Inputs  : nupop_inputs.rda (NUPOP.INPUTS.SC, NUPOP.INPUTS.SE)
### Output  : nupop_output.rda (NUPOP.OCC.SC, NUPOP.OCC.SE, one
###           occupancy vector per gene, with the promoter coords
###           and any reduced-flank or unscored regions attached)
###
### Each window runs in its own forked process, so the master keeps
### going and prints one progress line per round of windows.
### Requires the Bioconductor package NuPoP on the cluster:
###   BiocManager::install("NuPoP")
###############################################################

library('parallel')
library('NuPoP')

source("functions.R")

## Scores every chromosome in one species' inputs and returns each gene's
## promoter occupancy track. Rounds of windows run until every window is
## scored or reaches min_core_bp. Minus-strand slices are reversed so the
## promoter-proximal end sits at the end of every vector, matching
## tata_box_score() and poly_at_tract() on the same gene's sequence. The
## coords the tracks came from travel with the result, so
## score_promoters_nupop() can confirm they match the current promoters.
nupop_occupancy_cluster <- function(inputs, species = 7, model = 4,
                                    window_bp = 250000, flank = 7000, fallback_flank = 2000,
                                    min_core_bp = 5000, cores = 1) {
  ## Runs a table of core windows (seqid, start, end) with `flank` bp of
  ## context, one forked process per window, and returns one result per
  ## row. A window whose process ends early returns NULL or a try-error,
  ## which the caller reads as "split and rerun".
  nupop_run_windows <- function(chroms, tasks, flank, species, model, cores, work_root) {
    ## Predicts occupancy for one sequence in the current process and returns
    ## positions core_from to core_to. Each call works in its own temporary
    ## folder, because NuPoP writes its prediction file into the working
    ## directory. The folder name carries the process ID, which is unique
    ## among running processes, so forked workers that start from the same
    ## tempfile() state still receive separate folders. The folders sit under
    ## work_root, which nupop_occupancy_cluster() places beside R's session
    ## temp directory rather than inside it, so each worker's files stay
    ## independent of every other worker.
    nupop_predict_window <- function(seg_seq, core_from, core_to, species = 7, model = 4,
                                     work_root = tempdir()) {
      n <- nchar(seg_seq)
      ## NuPoP scores sequences of at least 148 bp, one nucleosome plus one base
      if (n < 148) return(rep(NA_real_, core_to - core_from + 1))
      wd <- tempfile(pattern = sprintf("nupop_pid%d_", Sys.getpid()), tmpdir = work_root)
      dir.create(wd, recursive = TRUE)
      old_wd <- setwd(wd)
      on.exit({ setwd(old_wd); unlink(wd, recursive = TRUE) }, add = TRUE)

      st <- seq(1, n, by = 80)
      writeLines(c(">window", substring(seg_seq, st, pmin(st + 79, n))), "window.fa")
      invisible(utils::capture.output(NuPoP::predNuPoP("window.fa", species = species, model = model)))
      pred <- paste0("window.fa_Prediction", model, ".txt")
      if (!file.exists(pred)) return(NULL)

      ## Columns are Position, P.start, Occup, N/L, Affinity. NULL skips a column.
      tab  <- scan(pred, what = list(0L, NULL, 0, NULL, NULL), skip = 1, quiet = TRUE)
      occ  <- rep(NA_real_, n)
      keep <- tab[[1]] >= 1 & tab[[1]] <= n & tab[[3]] >= 0
      occ[tab[[1]][keep]] <- tab[[3]][keep]
      occ[core_from:core_to]
    }

    parallel::mclapply(seq_len(nrow(tasks)), function(i) {
      s  <- chroms[[tasks$seqid[i]]]
      ws <- max(1, tasks$start[i] - flank)
      we <- min(nchar(s), tasks$end[i] + flank)
      nupop_predict_window(substr(s, ws, we), tasks$start[i] - ws + 1, tasks$end[i] - ws + 1,
                           species = species, model = model, work_root = work_root)
    }, mc.cores = cores, mc.preschedule = FALSE)
  }

  chroms <- inputs$chroms
  coords <- inputs$coords
  occ    <- lapply(chroms, function(s) rep(NA_real_, nchar(s)))

  ## Node-local scratch space for every window's files, removed at the end
  work_root <- file.path(dirname(tempdir()), sprintf("nupop_work_%d", Sys.getpid()))
  dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_root, recursive = TRUE), add = TRUE)

  tasks <- do.call(rbind, lapply(names(chroms), function(sq) {
    n <- nchar(chroms[[sq]])
    s <- seq(1, n, by = window_bp)
    data.frame(seqid = sq, start = s, end = pmin(s + window_bp - 1, n), stringsAsFactors = FALSE)
  }))
  window_ok <- function(r, t) is.numeric(r) && length(r) == t$end - t$start + 1

  reduced <- NULL
  failed  <- NULL
  round   <- 0
  while (nrow(tasks) > 0) {
    round <- round + 1
    t0    <- Sys.time()
    res   <- nupop_run_windows(chroms, tasks, flank, species, model, cores, work_root)
    ok    <- vapply(seq_len(nrow(tasks)), function(i) window_ok(res[[i]], tasks[i, ]), logical(1))
    for (i in which(ok)) occ[[tasks$seqid[i]]][tasks$start[i]:tasks$end[i]] <- res[[i]]

    bad   <- tasks[!ok, , drop = FALSE]
    small <- bad[bad$end - bad$start + 1 <= min_core_bp, , drop = FALSE]
    big   <- bad[bad$end - bad$start + 1 >  min_core_bp, , drop = FALSE]

    if (nrow(small) > 0) {
      res2 <- nupop_run_windows(chroms, small, fallback_flank, species, model, cores, work_root)
      ok2  <- vapply(seq_len(nrow(small)), function(i) window_ok(res2[[i]], small[i, ]), logical(1))
      for (i in which(ok2)) occ[[small$seqid[i]]][small$start[i]:small$end[i]] <- res2[[i]]
      if (any(ok2))  reduced <- rbind(reduced, cbind(small[ok2, , drop = FALSE], flank = fallback_flank))
      if (any(!ok2)) failed  <- rbind(failed, small[!ok2, , drop = FALSE])
    }

    mid   <- (big$start + big$end) %/% 2
    tasks <- rbind(data.frame(seqid = big$seqid, start = big$start, end = mid,     stringsAsFactors = FALSE),
                   data.frame(seqid = big$seqid, start = mid + 1,   end = big$end, stringsAsFactors = FALSE))
    cat(sprintf("[%s] round %d: %d windows scored, %d split for the next round, %.1f min\n",
                format(Sys.time(), "%H:%M:%S"), round, sum(ok), nrow(big),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
  }

  out <- setNames(vector("list", nrow(coords)), coords$gene)
  for (i in seq_len(nrow(coords))) {
    s <- occ[[coords$seqid[i]]][coords$start[i]:coords$end[i]]
    out[[coords$gene[i]]] <- if (coords$strand[i] == "-") rev(s) else s
  }
  attr(out, "coords")                <- coords
  attr(out, "reduced_flank_regions") <- reduced
  attr(out, "failed_regions")        <- failed
  out
}

load("nupop_inputs.rda")     # NUPOP.INPUTS.SC, NUPOP.INPUTS.SE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

## Summarizes the regions that needed fallback flanks or stayed unscored,
## so the job log records them next to the progress lines
report_regions <- function(occ, label) {
  red <- attr(occ, "reduced_flank_regions")
  bad <- attr(occ, "failed_regions")
  cat(sprintf("%s: %d promoters, %d region(s) scored with fallback flanks, %d region(s) unscored (%d bp)\n",
              label, length(occ), if (is.null(red)) 0L else nrow(red),
              if (is.null(bad)) 0L else nrow(bad),
              if (is.null(bad)) 0L else sum(bad$end - bad$start + 1)))
  if (!is.null(bad)) print(bad, row.names = FALSE)
}

t0 <- Sys.time()
cat(sprintf("NuPoP start: %d cores\n", NUM.CORES))

cat(sprintf("Sc: %d chromosomes, %.1f Mb\n", length(NUPOP.INPUTS.SC$chroms), sum(nchar(NUPOP.INPUTS.SC$chroms)) / 1e6))
NUPOP.OCC.SC <- nupop_occupancy_cluster(NUPOP.INPUTS.SC, cores = NUM.CORES)
report_regions(NUPOP.OCC.SC, "Sc")

cat(sprintf("Se: %d chromosomes, %.1f Mb\n", length(NUPOP.INPUTS.SE$chroms), sum(nchar(NUPOP.INPUTS.SE$chroms)) / 1e6))
NUPOP.OCC.SE <- nupop_occupancy_cluster(NUPOP.INPUTS.SE, cores = NUM.CORES)
report_regions(NUPOP.OCC.SE, "Se")

save(NUPOP.OCC.SC, NUPOP.OCC.SE, file = "nupop_output.rda")

cat(sprintf("NuPoP done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
