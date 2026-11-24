###############################################################
### dir_class_check.R
### Checks of the classification by direction (spec "Classification by direction", claude.ai artifact Y9Vn84qUaxXJVf8ZebEJEr, "Checks for Claude Code"):
###   definitions, identities, disjoint cell sets, split redraw, test formulas, gate logic, theta and wrap, p   deterministic checks of the functions
###   allocation, independence, class rule, no class bias, over and under                                     NB simulation through dir_contrasts_one()
### The simulation uses the cell numbers and depths of the Section 1 checkpoint when it exists (CHK.USE.CKPT1), stand-in values otherwise. Genes are simulated
### on the mean, which needs no dispersion fit, so the bootstrap runs locally; one gene runs the dispersion family to check its path. CHK.B replicates and
### CHK.GENES genes per simulated class set the run time (a few minutes at the defaults). The theta intervals are 95% intervals here (the pipeline uses
### 1 - 2 alpha_FCR), and the margins are the spec's delta_phi on the dispersion and CHK.DELTA.MEAN on the mean.
### Needs   : R/functions.R. Prints a table of checks; stops when a deterministic check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
if (!exists("CHK.B"))          CHK.B          <- 300
if (!exists("CHK.GENES"))      CHK.GENES      <- 12
if (!exists("CHK.DELTA.MEAN")) CHK.DELTA.MEAN <- 0.1
if (!exists("CHK.USE.CKPT1"))  CHK.USE.CKPT1  <- TRUE
set.seed(31)
RES <- list()
note <- function(id, pass, detail) { RES[[length(RES) + 1]] <<- data.frame(check = id, pass = pass, detail = detail); if (!isTRUE(pass)) message("FAILED: ", id, " ", detail) }

## ---- Contrast definitions and identities (no split: exact) ----
x <- list(Pc = 1.3, Pe = -0.4, Hc = 0.7, He = 0.2); H <- 0.45
ok <- all(abs(dir_parts("cis_trans", x)[c("A_all", "B_all", "S", "D")] - c(0.5, (1.3 + 0.4) - 0.5, 1.7, 0.5 - 1.2)) < 1e-12) &&
      all(abs(dir_parts("inheritance", x)[c("A_all", "B_all", "S", "D")] - c(H + 0.4, 1.3 - H, 1.7, 2 * H - 0.9)) < 1e-12) &&
      all(abs(dir_parts("trans_inh", x)[c("A_all", "B_all", "S", "D")] - c(0.2 + 0.4, 1.3 - 0.7, 1.2, 0.9 - 0.9)) < 1e-12)
note("Contrast definitions", ok, "A, B, S, D of the three schemes; second inheritance part x_Pc - x_H")
idn <- replicate(200, { x <- as.list(setNames(rnorm(4), c("Pc", "Pe", "Hc", "He")))
  max(vapply(.DIR.SCHEMES, function(sc) { r <- dir_parts(sc, x); max(abs(r[["S"]] - r[["A_all"]] - r[["B_all"]]), abs(r[["D"]] - r[["A_all"]] + r[["B_all"]])) }, numeric(1))) })
note("Identities (no split)", max(idn) < 1e-12, sprintf("max |S - (A + B)|, |D - (A - B)| = %.1e over 200 random designs", max(idn)))

## ---- Test formulas, including +-Inf and undefined replicates ----
xr <- c(0.5, 0.3, -0.1, Inf, -Inf, NA, 0.05)
t <- dir_margin_test(xr, 0.2)
nR <- min(sum(is.na(xr) | xr <= 0.2), sum(is.na(xr) | xr >= -0.2)); nE <- max(sum(is.na(xr) | xr >= 0.2), sum(is.na(xr) | xr <= -0.2))
note("Test formulas", t[["n0R"]] == nR && t[["n0E"]] == nE && t[["pR"]] == min(1, 2 * (1 + nR) / 8) && t[["pE"]] == (1 + nE) / 8 && nR == 4 && nE == 4,
     "p_R and p_E of Section 4; +Inf beyond +delta, -Inf beyond -delta, undefined on the null side of each test")

## ---- Gate logic over every combination of outcomes ----
comb <- as.matrix(expand.grid(A = 0:2, B = 0:2, S = 0:2, D = 0:2))       # 0 neither, 1 relevant, 2 equivalent
qR <- ifelse(comb == 1, 0.001, 0.5); qE <- ifelse(comb == 2, 0.001, 0.5); colnames(qR) <- colnames(qE) <- c("A", "B", "S", "D")
gt <- dir_gate(qR, qE)
exp_div <- rowSums(comb == 1) > 0; exp_cons <- comb[, "S"] == 2 & comb[, "D"] == 2
exp_gate <- ifelse(exp_div & exp_cons, "Inconclusive", ifelse(exp_div, "Divergent", ifelse(exp_cons, "Conserved", "Inconclusive")))
note("Gate logic", identical(gt$gate, exp_gate) && identical(gt$contradictory, exp_div & exp_cons),
     sprintf("81 combinations; %d contradictory cases flagged", sum(gt$contradictory)))

## ---- theta and wrap ----
ok <- dir_theta(1, 0) == 0 && dir_theta(0, 1) == 90 && dir_theta(1, 1) == 45 && dir_theta(-1, 1) == 135 && dir_theta(-1, -1) == 45 && dir_theta(-1, 0) == 0 &&
      abs(dir_recentre(179, 1) - (-1)) < 1e-12 && abs(dir_recentre(2, 178) - 182) < 1e-12
th0 <- 1; tb <- dir_recentre((th0 + rnorm(1000, 0, 3)) %% 180, th0)
ok <- ok && abs(median(tb) - th0) < 0.5 && diff(range(tb)) < 40
note("theta and wrap", ok, "atan2(B, A) mod 180; replicate angles re-expressed around a theta near 0 = 180 keep a narrow interval")

## ---- p ----
A <- runif(500, 0.1, 2); B <- runif(500, 0.1, 2)
ok <- max(abs(dir_share(A, B) - (1 - dir_theta(A, B) / 90))) < 1e-12
note("p", ok && is.na(dir_share(0.01, 0.01, 0.01, 0.01)), "p = 1 - theta / 90 for same-sign genes; undefined when neither part exceeds its noise")
## The correction for weak effects (SE 0.1): the bias of the mean p with and without it, over 20,000 normal draws of A and B per setting; it is checked to
## lower the bias where |A| and |B| sit within a few SE of zero and one part dominates. Where the raw p is already near unbiased it can add a little bias.
PB <- do.call(rbind, lapply(list(c(0.15, 0.05), c(0.3, 0.05), c(0.05, 0.15), c(0.05, 0.3), c(0.15, 0), c(0.3, 0.15)), function(ab) {
  A <- ab[1] + rnorm(20000, 0, 0.1); B <- ab[2] + rnorm(20000, 0, 0.1); tr <- dir_share(ab[1], ab[2])
  data.frame(A = ab[1], B = ab[2], true_p = tr, bias_raw = mean(dir_share(A, B), na.rm = TRUE) - tr, bias_corrected = mean(dir_share(A, B, 0.01, 0.01), na.rm = TRUE) - tr,
             undefined = mean(is.na(dir_share(A, B, 0.01, 0.01)))) }))
cat("Bias of p for weak effects (SE 0.1):\n"); print(PB, digits = 3, row.names = FALSE)
note("p bias correction", all(abs(PB$bias_corrected[1:5]) < abs(PB$bias_raw[1:5])),
     sprintf("|bias| raw %s; corrected %s", paste(sprintf("%.3f", abs(PB$bias_raw)), collapse = " "), paste(sprintf("%.3f", abs(PB$bias_corrected)), collapse = " ")))

## ---- Simulated design ----
if (CHK.USE.CKPT1 && file.exists(file.path("results", "checkpoints", "section1_checkpoint.rda"))) {
  e <- new.env(); load(file.path("results", "checkpoints", "section1_checkpoint.rda"), envir = e)
  D0 <- median(c(colSums(e$MIX.SC), colSums(e$MIX.SE))); OP <- list(Pc = colSums(e$MIX.SC) / D0, Pe = colSums(e$MIX.SE) / D0)
  LH <- colSums(e$HYB.SC) + colSums(e$HYB.SE); SH <- colSums(e$HYB.SC) / LH; DESIGN <- "Section 1 checkpoint"; rm(e)
} else {
  OP <- list(Pc = rlnorm(900, 0, 0.35), Pe = rlnorm(1000, 0, 0.35)); LH <- round(rlnorm(1500, log(4000), 0.35)); SH <- rbeta(1500, 600, 600); DESIGN <- "stand-in"
}
nH <- length(LH); LC <- rbinom(nH, LH, SH); D1 <- median(c(OP$Pc, OP$Pe)) * mean(LH) / 2
OH <- list(Hc = LC / D1 * 1, He = (LH - LC) / D1); ETA <- split_eta(LC / LH, LH)
## Classes on the inheritance plane of the mean: hybrid at the midparent (additive), at Pc or Pe (dominant), halfway toward Pc (partial), above both (over),
## below both (under); the parents differ by log2 1.2 in every class, so all classes carry the same total divergence.
sim_genes <- function(cls, n, mu0 = 4, d = 1.2, phi = 0.15) {
  h <- switch(cls, additive = 0, pc_dom = d / 2, pe_dom = -d / 2, partial_pc = d / 4, over = d, under = -d)
  data.frame(gene = paste0(cls, seq_len(n)), cls = cls, lpc = log2(mu0) + d / 2, lpe = log2(mu0) - d / 2, lh = log2(mu0) + h, phi = phi)
}
GS <- do.call(rbind, lapply(c("additive", "pc_dom", "pe_dom", "partial_pc", "over", "under"), sim_genes, n = CHK.GENES))
nbm <- function(lmu, o, phi) t(vapply(seq_along(lmu), function(i) rnbinom(length(o), mu = 2^lmu[i] * o, size = 1 / phi[i]), numeric(length(o))))
MATS <- list(MIX.SC = nbm(GS$lpc, OP$Pc, GS$phi), MIX.SE = nbm(GS$lpe, OP$Pe, GS$phi), HYB.SC = nbm(GS$lh, OH$Hc, GS$phi), HYB.SE = nbm(GS$lh, OH$He, GS$phi))
MATS <- lapply(MATS, function(m) { rownames(m) <- GS$gene; m })
EXPOS <- list(MIX.SC = OP$Pc, MIX.SE = OP$Pe, HYB.SC = OH$Hc, HYB.SE = OH$He); ETAS <- list(MIX.SC = 0, MIX.SE = 0, HYB.SC = ETA$Hc, HYB.SE = ETA$He)
FITS <- lapply(names(MATS), function(k) { f <- fit_counts_offset(MATS[[k]], EXPOS[[k]], eta = ETAS[[k]]); data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) }); names(FITS) <- names(MATS)

## Allocation from the per-cell variances of the ratio mean: Var(log2 mu) of a group is about 1 / (sum(y) ln2^2) times (1 + phi M), pooled over the genes
CHI <- c(mean = 0.4, disp = 0.35, kbal = 0.3, bsize = 0.33)
pv <- function(m, o) median(apply(m, 1, function(y) { M <- sum(y) / sum(o) * o; sum(M * (1 + 0.15 * M)) / sum(y)^2 }))
vpc <- pv(MATS$MIX.SC, OP$Pc); vpe <- pv(MATS$MIX.SE, OP$Pe); uhc <- nH * pv(MATS$HYB.SC, OH$Hc); uhe <- nH * pv(MATS$HYB.SE, OH$He)
AB <- lapply(CHI, function(ch) rbind(inh = dir_ab_fraction(vpe, vpc, 0.25 * (uhc + uhe), 0.25 * (uhc + uhe), (1 - ch) * nH), tinh = dir_ab_fraction(vpe, vpc, uhe, uhc, (1 - ch) * nH)))
PLAN <- make_dir_plan(c(Pc = length(OP$Pc), Pe = length(OP$Pe), H = nH), CHK.B, 77, CHI, AB)

## ---- Disjoint cell sets and the split redraw ----
d1 <- dir_draw(PLAN, 1); d2 <- dir_draw(PLAN, 2); d1b <- dir_draw(PLAN, 1)
okd <- all(vapply(names(CHI), function(p) { s <- d1[[p]]; pos <- function(v) match(v, d1$H)
  length(s$cis) == round(CHI[[p]] * nH) && length(s$inh_A) == round(AB[[p]]["inh", "f"] * length(s$trans)) && length(s$cis) + length(s$inh_A) + length(s$inh_B) == nH }, logical(1)))
## disjointness of the copies: the cis, A and B sets partition the resampled copies (as multisets of the drawn indices)
okm <- all(vapply(names(CHI), function(p) { s <- d1[[p]]; identical(sort(c(s$cis, s$inh_A, s$inh_B)), sort(d1$H)) && identical(sort(c(s$cis, s$tinh_A, s$tinh_B)), sort(d1$H)) }, logical(1)))
note("Disjoint cell sets", okd && okm, "cis, A and B copies partition each replicate's hybrid resample; sizes match chi and the stored A/B shares")
note("Split redraw", !identical(d1$disp$cis, d2$disp$cis) && identical(d1, d1b), "a new random partition in every replicate; the same draw from the stored seed")

## ---- The direction bootstrap on the simulated genes (mean family) ----
dl <- c(mean = CHK.DELTA.MEAN, disp = DIR.DELTA.PHI, alpha = DIR.DELTA.PHI, bsize = DIR.DELTA.PHI, kbal = 2 * DIR.DELTA.PHI)
t0 <- Sys.time()
TAB <- do.call(rbind, lapply(GS$gene, dir_contrasts_one, expos = EXPOS, fits = FITS, mats = MATS, etas = ETAS, plan = PLAN, delta = dl, family = "mean"))
TD <- dir_contrasts_one(GS$gene[1], EXPOS, FITS, MATS, ETAS, make_dir_plan(PLAN$ncells, 20, 78, CHI, AB), dl, "disp")
cat(sprintf("Direction bootstrap: %d genes x %d replicates (mean family) in %.1f min; design: %s\n", nrow(GS), CHK.B, as.numeric(difftime(Sys.time(), t0, units = "mins")), DESIGN))
note("Dispersion family path", nrow(TD) == 12 && all(abs(TD$S - TD$A - TD$B) < 1e-12 | !is.finite(TD$S)), "one gene, 4 quantities x 3 schemes, identities on the estimates")
inh <- TAB[TAB$scheme == "inheritance", ]; inh$cls <- GS$cls[match(inh$gene, GS$gene)]
## identities with the split: the replicate means of S - (A + B) and D - (A - B) near zero (A and B split, S and D all-cell), checked through the estimates of
## the split parts: the mean of the split replicates of A sits near the all-cell estimate
note("Identities (split)", median(abs(inh$A - (inh$A_lo + inh$A_hi) / 2) / ((inh$A_hi - inh$A_lo) / 3.92)) < 0.5,
     sprintf("median |all-cell A - centre of split interval| / SE = %.2f", median(abs(inh$A - (inh$A_lo + inh$A_hi) / 2) / ((inh$A_hi - inh$A_lo) / 3.92))))
note("Allocation", abs(log(median(inh$var_A / inh$var_B))) < log(1.5), sprintf("median Var A / Var B = %.2f (inheritance, mean)", median(inh$var_A / inh$var_B)))
note("Independence", abs(median(inh$cor_AB)) < 0.15 && median(inh$cor_AB_all) < -0.1,
     sprintf("median replicate cor(A, B): split %.2f, all cells %.2f", median(inh$cor_AB), median(inh$cor_AB_all)))

## ---- Class rule, no class bias, over and under ----
CL <- dir_classify_table(TAB); ci <- CL[CL$scheme == "inheritance", ]; ci$cls <- GS$cls[match(ci$gene, GS$gene)]
target <- c(additive = "Additive", pc_dom = "Pc dominant", pe_dom = "Pe dominant", partial_pc = "Partial toward Pc", over = "Overdominant", under = "Underdominant")
hit <- ci$class == target[ci$cls] | (ci$class == "Leaning" & (ci$lean_1 == target[ci$cls] | ci$lean_2 == target[ci$cls]))
rate <- tapply(hit, ci$cls, mean); exact <- tapply(ci$class == target[ci$cls], ci$cls, mean)
cat("Recovery on the inheritance plane (exact class / class or leaning toward it):\n"); print(round(rbind(exact = exact, with_leaning = rate), 2))
note("Class rule", all(rate >= 0.8), paste(sprintf("%s %.2f", names(rate), rate), collapse = ", "))
note("No class bias", abs(rate[["additive"]] - mean(rate[c("pc_dom", "pe_dom")])) < 0.2,
     sprintf("additive %.2f against dominant %.2f, equal total divergence", rate[["additive"]], mean(rate[c("pc_dom", "pe_dom")])))
ou <- ci[ci$cls %in% c("over", "under") & ci$gate %in% "Divergent", ]
note("Over and under", all(ifelse(ou$cls == "over", ou$D > 0, ou$D < 0)) && all(!(ou$class %in% c("Overdominant", "Underdominant")) | (ou$class == "Overdominant") == (ou$cls == "over")),
     "in the opposite-sign region the class follows the sign of D")
note("Outputs", all(c("gene", "quantity", "scheme", "A", "B", "S", "D", "pR_A_m2", "pE_D_m2", "theta", "theta_lo", "theta_hi", "class", "p", "p_lo", "p_hi", "k", "zero_fit", "failed_fit", "out_A") %in% names(CL)),
     "one row per gene, quantity and scheme with parts, outcomes, theta, class, p, k and flags")

CHECKS <- do.call(rbind, RES); print(CHECKS, row.names = FALSE)
stopifnot("a deterministic check failed" = all(CHECKS$pass[CHECKS$check %in% c("Contrast definitions", "Identities (no split)", "Test formulas", "Gate logic", "theta and wrap", "p", "Disjoint cell sets", "Split redraw")]))
