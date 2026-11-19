###############################################################
### sign_boot_test.R
### P1: the proposed bootstrap sign test for every dispersion-family contrast (phi, theta, zeta, kappa, on the log2 and raw scales), tested before it
### replaces the gate, the projection, the shifted null and the likelihood-bound comparison. The procedure under test (sign_boot_cell(), sign_boot_p()):
###   1. Observed fits: a failed fit gives no contrast; a zero phi gives no log2 estimate, while the raw estimate and both tests stand.
###   2. Bootstrap: resample cells (each parent; the undivided hybrid cells, a copy keeping its cis or trans membership), refit every group, compute the
###      contrast per replicate with +-Inf allowed on log2, a failed refit counting as zero.
###   3. Interval: the 2.5% and 97.5% percentiles of the replicates.
###   4. p = min(1, 2 (1 + n0) / (1 + B)), n0 the replicates on the null side of zero (opposite the observed sign), ties and undefined values (0/0,
###      +Inf - Inf) counted on the null side. The observed sign is that of the log2 contrast (+-Inf allowed), or of the raw contrast when the log2 value
###      is undefined.
### Measured, at the observed depths (Section 1 checkpoint; stand-in offsets otherwise) and 1000 or 2000 cells per group:
###   type I error at 0.05 and 0.01 when both groups are Poisson (phi = 0) and when they share a small dispersion (0.01, 0.03, 0.1), at mean counts of
###   0.2, 1 and 5 per parental cell, for the total, cis and trans contrasts;
###   power when one group's phi is 1.5 or 2 times the others' (effect on Pc, which moves total and trans, or on the Sc allele, which moves cis and
###   trans), under the proposal and under a variant that splits ties and drops undefined replicates (the cost of counting them on the null side);
###   the coverage of the percentile interval, the share of datasets with a contrast and with a log2 estimate, and the shares of tied and undefined
###   replicates.
### Blocks (P1.RUN): INPUTS writes sign_boot_inputs.rda for the cluster array (cluster/scripts/sign_boot_grid.R, cluster/submit/sign_boot.sub); LOCAL runs
### a reduced grid here (P1_local, status "reduced"); MERGE reads the array outputs and writes P1, results/checks/P1_sign_boot.csv and
### figures/extra/check_P1_sign_boot.pdf. The hybrid has no transcriptome-share variation here (c = 0.5 for every cell), so the share corrections are not
### part of this test.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT <- "analysis/checks/sign_boot_test.R"
CORES        <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.P1      <- 9200
if (!exists("P1.RUN")) P1.RUN <- c("INPUTS", "LOCAL", "MERGE")
if (!exists("P1.B"))       P1.B       <- 499   # bootstrap replicates per dataset (smallest p 2 / 500 = 0.004)
if (!exists("P1.R.NULL"))  P1.R.NULL  <- 500   # datasets per null cell (binomial SE 0.010 at 0.05)
if (!exists("P1.R.ALT"))   P1.R.ALT   <- 200   # datasets per power cell
if (!exists("P1.LOCAL.R")) P1.LOCAL.R <- 40; if (!exists("P1.LOCAL.B")) P1.LOCAL.B <- 199   # reduced local run; a value set before sourcing is kept
set.seed(SEED.SECTION + 28)

## ---- Observed depths ----
P1.POOLS <- if (file.exists(ckpt_path(1))) local({
  e <- new.env(); load(ckpt_path(1), envir = e)
  ref <- median(c(colSums(e$MIX.SC), colSums(e$MIX.SE)))
  list(O_Pc = unname(colSums(e$MIX.SC)) / ref, O_Pe = unname(colSums(e$MIX.SE)) / ref, O_H = unname(colSums(e$HYB.SC) + colSums(e$HYB.SE)) / ref / 2, source = "observed")
}) else local({ p <- check_design_pools(env = new.env()); list(O_Pc = p$O_Pc, O_Pe = p$O_Pe, O_H = p$S / 2, source = "stand-in") })
cat(sprintf("P1 depths: %s (Pc %d, Pe %d, hybrid %d cells in the pools)\n", P1.POOLS$source, length(P1.POOLS$O_Pc), length(P1.POOLS$O_Pe), length(P1.POOLS$O_H)))

## ---- The grid: null cells (every true contrast zero) and power cells ----
P1.GRID <- rbind(
  transform(expand.grid(N = c(1000, 2000), mean = c(0.2, 1, 5), phi = c(0, 0.01, 0.03, 0.1), ratio = 1, target = "none", stringsAsFactors = FALSE), R = P1.R.NULL),
  transform(expand.grid(N = c(1000, 2000), mean = c(0.2, 1, 5), phi = c(0.03, 0.1), ratio = c(1.5, 2), target = c("total", "cis"), stringsAsFactors = FALSE), R = P1.R.ALT))
P1.GRID$task <- seq_len(nrow(P1.GRID))
P1.STREAMS <- check_streams(nrow(P1.GRID), SEED.P1)

if ("INPUTS" %in% P1.RUN) {
  save(P1.GRID, P1.B, P1.POOLS, P1.STREAMS, file = file.path(INPUT.DIR, "sign_boot_inputs.rda"))
  cat(sprintf("Saved sign_boot_inputs.rda: %d cells, %d datasets, B = %d\n", nrow(P1.GRID), sum(P1.GRID$R), P1.B))
}

## ---- Reduced local run: Poisson and small-phi nulls and one power cell at 1000 cells and a mean count of 1 ----
if ("LOCAL" %in% P1.RUN) {
  LOC <- P1.GRID[P1.GRID$N == 1000 & P1.GRID$mean == 1 & ((P1.GRID$target == "none" & P1.GRID$phi %in% c(0, 0.03)) | (P1.GRID$target == "total" & P1.GRID$phi == 0.1 & P1.GRID$ratio == 2)), ]
  LOC$R <- P1.LOCAL.R
  t0 <- Sys.time()
  P1.LOCAL <- do.call(rbind, mclapply(seq_len(nrow(LOC)), sign_boot_cell, grid = LOC, pools = P1.POOLS, B = P1.LOCAL.B, streams = P1.STREAMS[LOC$task],
                                      mc.cores = CORES, mc.preschedule = FALSE))
  cat(sprintf("P1 local: %d cells x %d datasets x %d replicates in %.1f min\n", nrow(LOC), P1.LOCAL.R, P1.LOCAL.B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(P1.LOCAL[P1.LOCAL$quantity %in% c("phi", "kappa"), c("phi", "target", "ratio", "contrast", "quantity", "scale", "tested", "log2_estimate", "reject_05", "reject_05_split", "coverage", "tie_share", "undefined_share")],
        digits = 3, row.names = FALSE)
  check_save("P1_local", P1.LOCAL, "reduced",
             headline = sprintf("%s %s.%s.%s, phi %g, ratio %g: reject at 0.05 %.3f (SE %.3f), ties split %.3f", P1.LOCAL$target, P1.LOCAL$contrast, P1.LOCAL$quantity, P1.LOCAL$scale,
                                P1.LOCAL$phi, P1.LOCAL$ratio, P1.LOCAL$reject_05, P1.LOCAL$reject_05_se, P1.LOCAL$reject_05_split),
             mc_error = sprintf("binomial SE over %d datasets; %d replicates per dataset", P1.LOCAL.R, P1.LOCAL.B), script = CHECK.SCRIPT)
}

## ---- The array outputs ----
if ("MERGE" %in% P1.RUN) {
  P1.FILES <- list.files(OUTPUT.DIR, "^sign_boot_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
  if (!length(P1.FILES)) {
    check_save("P1", NULL, "not run", headline = "no sign_boot_output files in OUTPUT.DIR (cluster job pending)", script = CHECK.SCRIPT)
  } else {
    P1 <- list(); for (f in P1.FILES) { e <- new.env(); load(f, envir = e); P1[[f]] <- e$P1.PART }
    P1 <- do.call(rbind, P1); rownames(P1) <- NULL
    stopifnot("sign_boot outputs do not cover every cell of the inputs" = setequal(unique(P1$task), P1.GRID$task))
    write.csv(P1, file.path(CHECK.DIR, "P1_sign_boot.csv"), row.names = FALSE)
    NUL <- P1[P1$target == "none", ]
    ## Power: the contrasts the effect moves, the effect on Pc for total and trans and on the Sc allele for cis
    POW <- P1[(P1$target == "total" & P1$contrast %in% c("total", "trans")) | (P1$target == "cis" & P1$contrast == "cis"), ]
    hl_null <- vapply(split(NUL, paste(NUL$contrast, NUL$quantity, NUL$scale)), function(d)
      sprintf("%s %s %s null: type I at 0.05 from %.3f to %.3f (largest at phi %g, mean %g, N %d), above 0.05 + 2 SE in %d of %d cells", d$contrast[1], d$quantity[1], d$scale[1],
              min(d$reject_05), max(d$reject_05), d$phi[which.max(d$reject_05)], d$mean[which.max(d$reject_05)], d$N[which.max(d$reject_05)],
              sum(d$reject_05 > 0.05 + 2 * d$reject_05_se, na.rm = TRUE), nrow(d)), "")
    hl_pow <- vapply(split(POW, paste(POW$contrast, POW$quantity, POW$scale)), function(d)
      sprintf("%s %s %s power at 0.05: median %.3f (ties on the null side) against %.3f (ties split); largest loss %.3f", d$contrast[1], d$quantity[1], d$scale[1],
              median(d$reject_05), median(d$reject_05_split), max(d$reject_05_split - d$reject_05)), "")
    check_save("P1", list(cells = P1, grid = P1.GRID, settings = list(B = P1.B, pools = P1.POOLS$source)), "run", headline = c(hl_null, hl_pow),
               mc_error = sprintf("binomial SE over %d datasets per null cell and %d per power cell; %d replicates per dataset", P1.R.NULL, P1.R.ALT, P1.B), script = CHECK.SCRIPT)

    ## Figure (figures/extra/check_P1_sign_boot.pdf). Page 1: type I error at 0.05 against the common phi (0 = Poisson), one panel per contrast (rows) and
    ## quantity-scale (columns), a line per mean count (single-hue ramp), solid at 2000 cells and dashed at 1000, the nominal 0.05 dashed grey. Page 2:
    ## power at 0.05 against the effect ratio at phi 0.1 and 2000 cells, solid under the proposal and dashed with ties split.
    QS <- c("phi.log2", "phi.raw", "zeta.log2", "zeta.raw", "kappa.log2", "kappa.raw")
    MC <- setNames(colorRampPalette(c("#9ECAE1", "#08306B"))(3), c(0.2, 1, 5))
    fig_pdf("extra/check_P1_sign_boot.pdf", 2.6 * length(QS), 2.6 * 3 + 0.6)
    par(mfrow = c(3, length(QS)), mar = c(3.6, 3.8, 2.2, 0.6), mgp = c(2.1, 0.6, 0), cex.axis = 0.8)
    PH <- sort(unique(NUL$phi))
    for (cn in c("total", "cis", "trans")) for (qs in QS) {
      d <- NUL[NUL$contrast == cn & paste(NUL$quantity, NUL$scale, sep = ".") == qs, ]
      plot(NA, xlim = c(1, length(PH)), ylim = c(0, max(0.1, d$reject_05 + d$reject_05_se, na.rm = TRUE)), xaxt = "n", xlab = "common phi", ylab = "type I error at 0.05",
           main = sprintf("%s, %s", cn, qs), cex.main = 0.85, las = 1)
      axis(1, at = seq_along(PH), labels = PH); abline(h = 0.05, lty = 2, col = COLOR.GREY[["mid"]])
      for (m in names(MC)) for (n in c(1000, 2000)) { s <- d[d$mean == as.numeric(m) & d$N == n, ]; s <- s[order(s$phi), ]
        lines(match(s$phi, PH), s$reject_05, col = MC[[m]], lty = if (n == 2000) 1 else 2, lwd = 1.5) }
    }
    par(mfrow = c(3, length(QS)))
    for (cn in c("total", "cis", "trans")) for (qs in QS) {
      d <- POW[POW$contrast == cn & paste(POW$quantity, POW$scale, sep = ".") == qs & POW$phi == 0.1 & POW$N == 2000, ]
      plot(NA, xlim = c(1.5, 2), ylim = c(0, 1), xlab = "effect ratio", ylab = "power at 0.05", main = sprintf("%s, %s", cn, qs), cex.main = 0.85, las = 1)
      for (m in names(MC)) { s <- d[d$mean == as.numeric(m), ]; s <- s[order(s$ratio), ]
        lines(s$ratio, s$reject_05, col = MC[[m]], lwd = 1.5); lines(s$ratio, s$reject_05_split, col = MC[[m]], lwd = 1.5, lty = 2) }
    }
    dev.off()
  }
}
