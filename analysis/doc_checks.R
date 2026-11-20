###############################################################
### doc_checks.R
### Every check whose result the Methods document states or relies on, run from the checkpoints in one place. DOC.CHECKS maps each statement of the
### document to the checks that support it, where they run (a check script of analysis/checks/, a block of this script, or a table the pipeline writes
### into a checkpoint), and what they need. The runner sources each check script in DOC.PARTS whose checkpoints exist, runs the two simulation checks
### defined here (DOM1, M4), and writes results/checks/DOC_CHECKS.csv: one row per statement with the status and headline of each of its checks, read
### from results/checks/summary.csv or from the checkpoint that holds the pipeline's own table. A part whose inputs are missing is recorded with the
### reason instead of stopping the run, so the script can be rerun as checkpoints and cluster outputs arrive.
###   DOM1  bias of the dominance contrasts (mid, par) and of the full-cell total, cis and trans contrasts when the hybrid truly departs from its parents,
###         on the observed design (Section 2 checkpoint) or stand-in offsets, for gamma, lognormal and bimodal expression factors (dominance_bias_cell());
###   M4    whether the state-conditioned layer (analysis.R Section 2.5) flags the false calls that a rare state makes (check M3), with the state shared by
###         every group or present in the parents only, the two-culture case of the data (rare_cond_cell()); needs only the Section 1 checkpoint.
### The checks keep their own seeds and streams, so a part reruns identically on its own. Cluster-dependent parts read outputs already copied back:
### the sign-test array (sign_boot_test.R MERGE: P1) and the artefact array (artefact_tests.R T0, T8, T12).
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT   <- "analysis/doc_checks.R"
CORES          <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.DOC       <- 9100
## Settings; a value set before sourcing this script is kept, so a run can select parts or sizes.
if (!exists("DOC.PARTS"))      DOC.PARTS      <- c("sim", "data", "leak", "artefact", "mse", "state_groups", "state_conditioned", "sign_boot", "DOM1", "M4", "text")   # text last: it writes REPORT.md
if (!exists("DOC.CHECK.SIZE")) DOC.CHECK.SIZE <- 1     # replicate scale of checks_sim.R (1 = full size)
if (!exists("DOM1.G"))         DOM1.G         <- 300   # genes per DOM1 cell
if (!exists("M4.R"))           M4.R           <- 100   # datasets per M4 cell
if (!exists("M4.B"))           M4.B           <- 99    # bootstrap draws per dataset
set.seed(SEED.SECTION + 27)

## ---- What the document states and which checks support it ----
## section: the Methods section; statement: what the check supports; checks: check IDs (results/checks/<ID>.rds) or pipeline tables (CAPITALS, held in the
## checkpoint named in needs); part: the DOC.PARTS entry that runs it ("pipeline" for analysis.R itself); needs: checkpoints and cluster outputs.
DOC.CHECKS <- data.frame(rbind(
  c("Sequencing depth offset", "A gene's own counts in its library bias its dispersion little", "A7", "sim, data", "Section 1 (data arm)"),
  c("Estimating mu", "The ratio mean costs at most 1 + CV2(O) in variance; realized cost per gene", "G3, A5", "text, data", "Section 2"),
  c("Estimating phi", "Slope rule, optimizer, and fits at the ends of the search range", "A2, A3", "sim, data", "Section 2 (A3)"),
  c("Estimating phi", "Cox-Reid bias against the unadjusted and moment estimators; truncation at zero", "A1, A9", "sim", "none"),
  c("Estimating phi; Notes: non-gamma factors", "Non-gamma expression factors bias Cox-Reid; the data's third moment", "B4, A11, L4", "sim, data, leak", "Section 1"),
  c("Notes: Cox-Reid adjustment and estimator choice", "Bias, variance and MSE of both estimators per group and per contrast at the observed depths", "M1, M2", "mse", "Section 1"),
  c("Burst kinetics", "The dispersion is estimated far less precisely than the mean", "PILOT.PRECISION", "pipeline", "Section 2"),
  c("Burst kinetics", "Burst size and kinetic balance couple to the mean and dispersion only as their definitions force", "STRUCT.COUPLING", "pipeline", "Section 3"),
  c("Burst kinetics", "UMI capture thins counts binomially and leaves phi unchanged", "A8", "data", "Section 2"),
  c("Negative binomial model: CV2 and Fano", "CV2 and the Fano factor are reported for every gene and group, without contrasts or tests", "GENE.TABLE", "pipeline", "Section 2"),
  c("Negative binomial model: dropout", "UMI counts need no dropout term beyond the NB", "DROPOUT.AIC", "pipeline", "Section 2"),
  c("Transcriptome share corrections", "The observed share is binomial; the share-variance estimator", "B1, B2", "data, sim", "Section 2 (data arms)"),
  c("Transcriptome share corrections", "The share corrections close the loop; fits truncated by the correction", "B3, B6", "sim, data", "Section 2 (B6)"),
  c("Per-genome hybrid estimates", "Route A against the summed counts; the per-genome mean", "B5, F5, ROUTE.CHECK", "data, pipeline", "Section 2"),
  c("Per-genome hybrid estimates; dominance", "The dominance contrasts are unbiased when the hybrid departs from its parents", "DOM1", "DOM1", "none (observed design with Section 2)"),
  c("Contrasts", "Gene sets and exclusions per contrast and scale (observed-fit rules)", "E4", "data", "Section 2"),
  c("Contrasts: bootstrap intervals", "Coverage of the percentile intervals with zero refits at -Inf", "E3", "sim", "none"),
  c("Statistical power", "Partition fractions and sizes; cis variance falls as 1/n; equal precision; depth balance; call stability", "C1, C2, C3, C6", "data", "Section 2"),
  c("Model assumptions: large samples", "Every partition keeps about 800 cells", "PARTITION.REPORT", "pipeline", "Section 2"),
  c("Separate partitions for each quantity", "Mean and dispersion errors covary on shared cells by mu (E[(a - 1)^3] - 2 phi^2) / N", "A4, A10, C7", "sim, data", "Section 2 (data arms)"),
  c("Tests of contrasts", "The bootstrap sign test: type I error for Poisson and small equal dispersions, power, cost of ties on the null side", "P1, P1_local", "sign_boot", "sign_boot array outputs (P1); Section 1 (P1_local)"),
  c("Tests of contrasts", "The cis and trans sign tests keep their size when f is estimated from the same cells", "C4", "sim", "none"),
  c("Multiple testing", "Genes at the p-value floor 2 / (B + 1) and the reach of BH", "D9", "data", "Section 2"),
  c("Appendix A", "Bias of the moment estimator; convergence and agreement on the data", "F1", "sim, data", "Section 2 (data arm)"),
  c("Appendix B", "Asymptotics, the observed N and unequal intrinsic variances", "F2, F3, F4", "sim", "none"),
  c("Model assumptions (departures from the NB)", "Why the data depart from the NB: low-quality cells, ambient RNA, depth scaling, cell states, outlier cells", "T0, T1, T2, T3, T4, T6, T7, T8, T9, T10, T11, T12", "artefact", "Section 1; artefact array outputs (T0, T8, T12)"),
  c("Model assumptions (departures from the NB)", "Estimator comparison and the leak of mean differences into dispersion contrasts", "L1, L2, L3, L5, L6, L7, L8", "leak", "Section 1"),
  c("Appendix C: state-conditioned estimates", "A rare state inflates total and trans dispersion tests; conditioning flags those calls", "M3, M4", "mse, M4", "Section 1"),
  c("Appendix C: state-conditioned estimates", "Cell states, the program's frequency per culture, and state-conditioned dispersions", "ST1, SC1, SC2, SC3", "state_groups, state_conditioned", "Section 1"),
  c("Appendix C: state-conditioned estimates", "Calls that the program carries", "STATE.FLAG.SUMMARY", "pipeline", "Section 2 (gene_cond outputs)"),
  c("Every section", "Every number the text quotes, regenerated", "G1", "data", "Section 2"),
  c("References", "Every DOI resolves to its citation", "G2", "text", "network access to api.crossref.org")),
  stringsAsFactors = FALSE)
names(DOC.CHECKS) <- c("section", "statement", "checks", "part", "needs")

## Where each pipeline table lives
DOC.PIPELINE <- c(PILOT.PRECISION = 2, STRUCT.COUPLING = 3, GENE.TABLE = 2, DROPOUT.AIC = 2, ROUTE.CHECK = 2, PARTITION.REPORT = 2, STATE.FLAG.SUMMARY = 2)

## ---- The check scripts and what each needs ----
HAVE <- vapply(1:8, function(n) file.exists(ckpt_path(n)), logical(1))
DOC.SCRIPTS <- list(
  sim               = list(file = "analysis/checks/checks_sim.R",        needs = integer(0), settings = list(CHECK.SIZE = DOC.CHECK.SIZE)),
  data              = list(file = "analysis/checks/checks_data.R",       needs = 1:2,        settings = list()),
  leak              = list(file = "analysis/checks/estimator_leak.R",    needs = 1,          settings = list()),
  artefact          = list(file = "analysis/checks/artefact_tests.R",    needs = 1,          settings = list(ART.LOCAL = FALSE)),
  mse               = list(file = "analysis/checks/estimator_mse.R",     needs = 1,          settings = list(MSE.RUN = c("M1", "M2", "M3"))),
  state_groups      = list(file = "analysis/checks/state_groups.R",      needs = 1,          settings = list()),
  state_conditioned = list(file = "analysis/checks/state_conditioned.R", needs = 1,          settings = list()),
  sign_boot         = list(file = "analysis/checks/sign_boot_test.R",    needs = integer(0), settings = list(P1.RUN = c("INPUTS", "MERGE"))),
  text              = list(file = "analysis/checks/checks_text.R",       needs = integer(0), settings = list()))

## Each part runs in its own environment with its settings, so one part's objects never reach another; an error is recorded and the run goes on.
DOC.RUNS <- list()
for (part in setdiff(DOC.PARTS, c("DOM1", "M4", "text"))) {
  sc <- DOC.SCRIPTS[[part]]
  missing_ck <- sc$needs[!HAVE[sc$needs]]
  if (length(missing_ck)) { DOC.RUNS[[part]] <- data.frame(part = part, status = sprintf("not run: checkpoint %s missing", paste(missing_ck, collapse = ", ")), minutes = NA_real_); next }
  env <- new.env(parent = globalenv()); for (nm in names(sc$settings)) assign(nm, sc$settings[[nm]], envir = env)
  t0 <- Sys.time()
  msg <- tryCatch({ source(here(sc$file), local = env); "run" }, error = function(e) paste("error:", conditionMessage(e)))
  DOC.RUNS[[part]] <- data.frame(part = part, status = msg, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  cat(sprintf("[doc_checks] %s: %s (%.1f min)\n", part, msg, DOC.RUNS[[part]]$minutes))
}

##############################################################################
## DOM1. Bias of the dominance contrasts when the hybrid departs from its parents (Methods: Per-genome hybrid estimates; Notes: full-pipeline simulation)
if ("DOM1" %in% DOC.PARTS) {
  DOM1.ENV <- new.env()
  if (HAVE[2]) load(ckpt_path(2), envir = DOM1.ENV)
  DOM1.POOLS <- check_design_pools(env = DOM1.ENV)
  rm(DOM1.ENV)
  DOM1.GRID <- expand.grid(mean = c(0.5, 2, 8), phi = c(0.1, 0.5), family = c("gamma", "lognormal", "bimodal"), omega = 0.5, stringsAsFactors = FALSE)
  DOM1.STREAMS <- check_streams(nrow(DOM1.GRID), SEED.DOC + 1)
  t0 <- Sys.time()
  DOM1 <- do.call(rbind, mclapply(seq_len(nrow(DOM1.GRID)), dominance_bias_cell, grid = DOM1.GRID, pools = DOM1.POOLS, G = DOM1.G, streams = DOM1.STREAMS,
                                  mc.cores = CORES, mc.preschedule = FALSE))
  DOM1.RUN <- data.frame(part = "DOM1", status = "run", minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  DOM1.Z <- DOM1[abs(DOM1$bias) > 2 * DOM1$bias_se & is.finite(DOM1$bias_se), ]
  print(DOM1[DOM1$contrast %in% c("mid", "par_sc", "par_se"), c("mean", "phi", "family", "quantity", "contrast", "genes", "bias", "bias_se", "bias_over_sd")], digits = 3, row.names = FALSE)
  check_save("DOM1", list(cells = DOM1, grid = DOM1.GRID, design = DOM1.POOLS$source), "run",
             headline = c(sprintf("design %s, %d genes per cell", DOM1.POOLS$source, DOM1.G),
                          vapply(c("mid", "par_sc", "par_se", "cis", "trans"), function(md) { d <- DOM1[DOM1$contrast == md & DOM1$quantity %in% c("mean", "disp"), ]
                            sprintf("%s: largest |bias| / SD %.3f (%s, %s, mean %g, phi %g); |bias| > 2 MC SE in %d of %d cells", md, max(abs(d$bias_over_sd), na.rm = TRUE),
                                    d$quantity[which.max(abs(d$bias_over_sd))], d$family[which.max(abs(d$bias_over_sd))], d$mean[which.max(abs(d$bias_over_sd))],
                                    d$phi[which.max(abs(d$bias_over_sd))], sum(abs(d$bias) > 2 * d$bias_se, na.rm = TRUE), nrow(d)) }, "")),
             mc_error = "MC SE = SD over genes / sqrt(genes) per cell", script = CHECK.SCRIPT)
  DOC.RUNS$DOM1 <- DOM1.RUN
}

##############################################################################
## M4. Does the state-conditioned layer flag the false calls of a rare state? (Methods Appendix C)
if ("M4" %in% DOC.PARTS) {
  if (!HAVE[1]) {
    check_save("M4", NULL, "not run", headline = "the Section 1 checkpoint is missing", script = CHECK.SCRIPT)
    DOC.RUNS$M4 <- data.frame(part = "M4", status = "not run: checkpoint 1 missing", minutes = NA_real_)
  } else {
    M4.ENV <- new.env(); load(ckpt_path(1), envir = M4.ENV)
    M4.EXPOS <- with(M4.ENV, { ref <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
      list(Pc = unname(colSums(MIX.SC)) / ref, Pe = unname(colSums(MIX.SE)) / ref, H = unname(colSums(HYB.SC) + colSums(HYB.SE)) / ref / 2) })
    rm(M4.ENV)
    M4.HYB.SPLIT <- local({ ord <- order(M4.EXPOS$H, decreasing = TRUE); sp <- split_indices_by_depth(length(ord), 0.5); list(c = ord[sp$c], t = ord[sp$t]) })
    M4.GRID <- rbind(data.frame(mu = c(2, 8.9), phi = 0.29, p = 0, J = 1, scenario = "shared"),
                     expand.grid(mu = c(2, 8.9), phi = 0.29, p = c(0.003, 0.01), J = 50, scenario = c("shared", "parents"), stringsAsFactors = FALSE))
    M4.STREAMS <- check_streams(nrow(M4.GRID), SEED.DOC + 2)
    t0 <- Sys.time()
    M4 <- do.call(rbind, mclapply(seq_len(nrow(M4.GRID)), rare_cond_cell, grid = M4.GRID, expos = M4.EXPOS, hyb_split = M4.HYB.SPLIT, R = M4.R, B = M4.B,
                                  streams = M4.STREAMS, mc.cores = CORES, mc.preschedule = FALSE))
    DOC.RUNS$M4 <- data.frame(part = "M4", status = "run", minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    print(M4, digits = 3, row.names = FALSE)
    check_save("M4", list(cells = M4, grid = M4.GRID, settings = list(R = M4.R, B = M4.B)), "run",
               headline = sprintf("%s, %s, p %g, mu %g: type I %.3f (SE %.3f), flagged %.2f of false calls, type I after flags %.3f, flag rate %.3f",
                                  M4$scenario, M4$contrast, M4$p, M4$mu, M4$type1_05, M4$type1_05_se, M4$flagged_of_false_calls, M4$type1_05_unflagged, M4$flag_rate),
               mc_error = sprintf("binomial SE over the tested datasets of %d per cell; %d bootstrap draws per dataset", M4.R, M4.B), script = CHECK.SCRIPT)
  }
}

## The text audit last: it reads summary.csv and writes REPORT.md.
if ("text" %in% DOC.PARTS) {
  t0 <- Sys.time(); env <- new.env(parent = globalenv())
  msg <- tryCatch({ source(here(DOC.SCRIPTS$text$file), local = env); "run" }, error = function(e) paste("error:", conditionMessage(e)))
  DOC.RUNS$text <- data.frame(part = "text", status = msg, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

##############################################################################
## ---- One row per statement of the document, with the status of each of its checks ----
DOC.RUNS <- do.call(rbind, DOC.RUNS); rownames(DOC.RUNS) <- NULL
print(DOC.RUNS, row.names = FALSE)
SUMMARY <- if (file.exists(file.path(CHECK.DIR, "summary.csv"))) read.csv(file.path(CHECK.DIR, "summary.csv"), stringsAsFactors = FALSE) else data.frame(ID = character(0), status = character(0), headline = character(0))
CKPT.OBJECTS <- lapply(setNames(nm = unique(DOC.PIPELINE)), function(k) if (HAVE[k]) load(ckpt_path(k), envir = new.env()) else character(0))   # object names held by each checkpoint
DOC.STATUS <- DOC.CHECKS
DOC.STATUS$status <- DOC.STATUS$headline <- ""
for (i in seq_len(nrow(DOC.STATUS))) {
  ids <- trimws(strsplit(DOC.STATUS$checks[i], ",")[[1]])
  st <- character(0); hl <- character(0)
  for (id in ids) {
    if (id %in% names(DOC.PIPELINE)) {
      k <- DOC.PIPELINE[[id]]
      ok <- id %in% CKPT.OBJECTS[[as.character(k)]]
      st <- c(st, sprintf("%s: %s", id, if (ok) sprintf("in checkpoint %d", k) else if (!HAVE[k]) sprintf("checkpoint %d missing", k) else sprintf("not in checkpoint %d yet (rerun analysis.R)", k)))
      next
    }
    row <- SUMMARY[SUMMARY$ID == id | startsWith(SUMMARY$ID, paste0(id, "_")) | startsWith(SUMMARY$ID, paste0(id, " ")), , drop = FALSE]
    if (!nrow(row)) { st <- c(st, sprintf("%s: no result yet", id)); next }
    st <- c(st, sprintf("%s: %s (%s)", row$ID, row$status, if ("date" %in% names(row)) row$date else ""))
    hl <- c(hl, sprintf("%s: %s", row$ID, substr(row$headline, 1, 300)))
  }
  DOC.STATUS$status[i] <- paste(st, collapse = "; "); DOC.STATUS$headline[i] <- paste(hl, collapse = " | ")
}
dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(DOC.STATUS, file.path(CHECK.DIR, "DOC_CHECKS.csv"), row.names = FALSE)
cat(sprintf("\nDocument checks: %d statements; %d with every check run or in a checkpoint\n", nrow(DOC.STATUS),
            sum(!grepl("no result yet|not run|not in checkpoint|missing", DOC.STATUS$status))))
print(DOC.STATUS[, c("section", "checks", "status")], row.names = FALSE)
