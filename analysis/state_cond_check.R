###############################################################
### state_cond_check.R
### Local checks of the state-conditioned layer (analysis.R Section 2.5, Methods appendix: State-conditioned estimates): the conditioned NB fit
### (fit_dispersion(z = )) reduces to the unconditioned fit for a constant covariate, its adjusted slope at zero matches the closed form, it recovers phi and the
### covariate's effect when expression follows a cell-state score, and it reaches the Poisson boundary at the expected rate; cond_contrasts_one() reproduces the
### unconditioned estimates of boot_contrasts_one() on the same draws, keeps its gate, leaves every contrast near unchanged when no state differs between groups,
### and, when a program occurs in one parent only, moves the trans and total dispersion contrasts towards zero with a change interval that excludes zero.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(21)

## ---- the conditioned fit ----
n <- 600; e <- rlnorm(n, 0, 0.5); z <- rnorm(n); y <- rnbinom(n, mu = 3 * e * exp(0.7 * z), size = 1 / 0.3)
## a constant covariate carries no information and gives the unconditioned fit
a <- fit_dispersion(y, e, 0.02, "NB", z = rep(1.5, n)); b <- fit_dispersion(y, e, 0.02, "NB")
stopifnot(isTRUE(all.equal(a[names(b)], b)), a[["beta1"]] == 0)
## the adjusted slope at phi = 0+ equals one half of sum((y - m)^2) - sum(m) + tr[(X' M X)^-1 X' M^2 X]
pl <- .pois_loglinear(y, e, z); ll <- .cr_loglik(y, pl$m, pl$X); h <- 1e-5; X <- pl$X; m <- pl$m
slope_num <- (4 * ll(h) - ll(2 * h) - 3 * ll(0)) / (2 * h)
slope_cf  <- 0.5 * (sum((y - m)^2) - sum(m) + sum(diag(solve(crossprod(X, m * X), crossprod(X, m^2 * X)))))
stopifnot(abs(slope_num / slope_cf - 1) < 1e-3, pl$converged == 1, abs(sum(pl$m) - sum(y)) < 1e-6 * sum(y))
## an intercept-only design gives the adjustment of the ratio estimate
stopifnot(isTRUE(all.equal(.cr_loglik(y, 2 * e, matrix(1, n))(0.3), .cr_loglik(y, 2 * e)(0.3))))
## recovery: the conditioned phi and beta match the truth on average, the unconditioned phi carries the state's variance
est <- t(replicate(200, { z <- rnorm(n); y <- rnbinom(n, mu = 3 * e * exp(0.7 * z), size = 1 / 0.3)
                          c(fit_dispersion(y, e, 0, "NB", z = z)[c("disp", "beta1", "mu")], fit_dispersion(y, e)[["disp"]]) }))
stopifnot(abs(mean(est[, 1]) / 0.3 - 1) < 0.05, abs(mean(est[, 2]) - 0.7) < 0.02, abs(mean(est[, 3]) / 3 - 1) < 0.03, mean(est[, 4]) > 2 * 0.3)
## Poisson counts with a state effect sit at the boundary about half the time
bd <- mean(replicate(400, { z <- rnorm(n); fit_dispersion(rpois(n, 3 * e * exp(0.5 * z)), e, 0, "NB", z = z)[["disp"]] == 0 }))
stopifnot(bd > 0.4, bd < 0.6)
## the weighted estimator takes no covariate
stopifnot(inherits(try(fit_dispersion(y, e, 0, "weighted", z = z), silent = TRUE), "try-error"))

## ---- a design with every dataset: parents in one tube, hybrid cells with the four partitions, a program score per cell ----
## prog_pc is the share of Pc cells in the program (expression times J); Pe and the hybrid cells never enter it.
G <- 4; genes <- paste0("Y", seq_len(G)); nP <- 300; nH <- 400; J <- 20; B <- 60
libP1 <- rlnorm(nP, 8, 0.3); libP2 <- rlnorm(nP, 8, 0.3); libH <- rlnorm(nH, 8.3, 0.3); chat <- pmin(pmax(rnorm(nH, 0.55, 0.03), 0.3), 0.8)
ref <- median(c(libP1, libP2)); mu <- runif(G, 1, 4); phi <- runif(G, 0.15, 0.4)
splits <- list(HYC = sort(sample(nH, nH / 2)), HYC.N = sort(sample(nH, nH / 2)), HYC.K = sort(sample(nH, nH / 2)), HYC.Z = sort(sample(nH, nH / 2)))
draws <- make_draws(c(MIX.SC = nP, MIX.SE = nP, HYB.COMB = nH), B, 3, splits)
st <- share_terms(chat, libH); sr <- share_r_groups(st$r_c, st$r_e)
designs <- list()
for (scen in c("none", "pc_only")) {
  inprog <- if (scen == "pc_only") runif(nP) < 0.06 else rep(FALSE, nP)
  fold <- ifelse(inprog, J, 1)
  P1 <- matrix(rnbinom(G * nP, mu = outer(mu, libP1 / ref * fold), size = 1 / phi), G, nP, dimnames = list(genes, paste0("pc", seq_len(nP))))
  P2 <- matrix(rnbinom(G * nP, mu = outer(mu * 1.3, libP2 / ref), size = 1 / phi), G, nP, dimnames = list(genes, paste0("pe", seq_len(nP))))
  HC <- matrix(rnbinom(G * nH, mu = outer(mu, libH / ref * mean(chat)), size = 1 / phi), G, nH, dimnames = list(genes, paste0("h", seq_len(nH))))
  HE <- matrix(rnbinom(G * nH, mu = outer(mu * 1.3, libH / ref * (1 - mean(chat))), size = 1 / phi), G, nH, dimnames = list(genes, paste0("h", seq_len(nH))))
  mats <- list(MIX.SC = P1, MIX.SE = P2, HYB.SC = HC, HYB.SE = HE)
  expos <- list(MIX.SC = libP1 / ref, MIX.SE = libP2 / ref, HYB.SC = libH / ref / 2, HYB.SE = libH / ref / 2)
  for (k in names(splits)) for (al in c("SC", "SE")) {
    suf <- sub("^HYC", "", k); cc <- splits[[k]]; tt <- setdiff(seq_len(nH), cc); src <- mats[[paste0("HYB.", al)]]
    mats[[paste0("HYC.", al, suf)]] <- src[, cc]; mats[[paste0("HYT.", al, suf)]] <- src[, tt]
    expos[[paste0("HYC.", al, suf)]] <- libH[cc] / ref / 2; expos[[paste0("HYT.", al, suf)]] <- libH[tt] / ref / 2
  }
  fits <- list()
  for (k in .CONTRAST.GROUPS) { f <- t(vapply(genes, function(g) fit_dispersion(mats[[k]][g, ], expos[[k]], sr[[k]], "NB")[c("mu", "disp")], numeric(2)))
                                fits[[k]] <- data.frame(MU = f[, "mu"], DISP = f[, "disp"], row.names = genes) }
  ## program module: 20 marker genes' summed count, raised J-fold in program cells
  state <- list(cells = list(Pc = list(total = setNames(rpois(nP, 60 * libP1 / ref * fold), colnames(P1)), lib = setNames(libP1 / ref, colnames(P1))),
                             Pe = list(total = setNames(rpois(nP, 60 * libP2 / ref), colnames(P2)), lib = setNames(libP2 / ref, colnames(P2))),
                             H  = list(total = setNames(rpois(nH, 120 * libH / ref), colnames(HC)), lib = setNames(libH / ref, colnames(HC)))),
                members = character(0), pseudo = 0.5)
  designs[[scen]] <- list(mats = mats, expos = expos, fits = fits, state = state)
}

## ---- cond_contrasts_one ----
d0 <- designs$none
cc0 <- do.call(rbind, lapply(genes, cond_contrasts_one, expos = d0$expos, mats = d0$mats, draws = draws, share_cells = list(chat = chat, lib = libH),
                             share_obs = list(cbar = st$cbar, r = sr), state = d0$state))
bc0 <- do.call(rbind, lapply(genes, boot_contrasts_one, expos = d0$expos, fits = d0$fits, mats = d0$mats, draws = draws, share_cells = list(chat = chat, lib = libH),
                             share_obs = list(cbar = st$cbar, r = sr)))
## the unconditioned estimates of the layer are those of the main bootstrap
for (q in c("mean", "disp", "rdisp", "bsize", "kbal", "rbsize", "rkbal")) for (md in .OUT_MODES)
  stopifnot(isTRUE(all.equal(cc0[[paste0("c", q, "_", md, "_est")]] + cc0[[paste0("d", q, "_", md, "_est")]], bc0[[paste0(q, "_", md, "_est")]], tolerance = 1e-8)))
## no state difference: the change intervals of the log2 dispersion contrasts cover zero
cover <- unlist(lapply(.OUT_MODES, function(md) cc0[[paste0("ddisp_", md, "_lo")]] <= 0 & cc0[[paste0("ddisp_", md, "_hi")]] >= 0))
stopifnot(mean(cover) >= 0.9)
## the gate sets the undefined contrasts of both outputs to NA
gates <- matrix(TRUE, G, 6 * length(.OUT_MODES), dimnames = list(genes, as.vector(outer(c("mean", "raw", "log", "bsize", "rbsize", "kbal"), .OUT_MODES, paste, sep = "_"))))
gates[1, "log_trans"] <- FALSE
cg <- cond_contrasts_one(genes[1], d0$expos, d0$mats, draws, list(chat = chat, lib = libH), list(cbar = st$cbar, r = sr), d0$state, gates)
stopifnot(is.na(cg$cdisp_trans_est), is.na(cg$ddisp_trans_lo), is.na(cg$cbsize_trans_est) == FALSE, !is.na(cg$cdisp_cis_est))
## a gene in the program module never enters its own score
st_in <- d0$state; st_in$members <- genes[1]
s_out <- state_scores(genes[1], d0$mats, d0$state); s_in <- state_scores(genes[1], d0$mats, st_in)
stopifnot(!isTRUE(all.equal(s_out$MIX.SC, s_in$MIX.SC)), all(abs(vapply(s_in, median, numeric(1))[c("MIX.SC", "MIX.SE", "HYB.SC")]) < 1e-12),
          identical(s_in$HYB.SC, s_in$HYB.SE), identical(s_in$HYC.SC, s_in$HYB.SC[splits$HYC]))

## program in Pc only: the trans and total dispersion contrasts are inflated, and conditioning brings them back with a change interval that excludes zero
d1 <- designs$pc_only
cc1 <- do.call(rbind, lapply(genes, cond_contrasts_one, expos = d1$expos, mats = d1$mats, draws = draws, share_cells = list(chat = chat, lib = libH),
                             share_obs = list(cbar = st$cbar, r = sr), state = d1$state))
for (md in c("total", "trans")) {
  u <- cc1[[paste0("cdisp_", md, "_est")]] + cc1[[paste0("ddisp_", md, "_est")]]
  stopifnot(all(u > 1), median(abs(cc1[[paste0("cdisp_", md, "_est")]])) < median(u) / 3, mean(cc1[[paste0("ddisp_", md, "_lo")]] > 0) >= 0.75,
            mean(cc1[[paste0("cdisp_", md, "_lo")]] <= 0 & cc1[[paste0("cdisp_", md, "_hi")]] >= 0) >= 0.75)
}
## cis shares the hybrid tube and stays unchanged
stopifnot(all(cc1$ddisp_cis_lo <= 0 & cc1$ddisp_cis_hi >= 0))
## the program's effect on the gene is recovered in Pc (log2 mean per log2 score: the gene and the module rise together)
stopifnot(all(abs(cc1$cbeta_MIX.SC - 1) < 0.25))

## ---- flags for the figures ----
BCx <- data.frame(gene = 1:3, stateflag_disp_trans = c(TRUE, NA, FALSE), stateflag_disp_cis = c(FALSE, FALSE, TRUE))
stopifnot(identical(state_flag_vec(BCx, "disp", c("cis", "trans")), c(TRUE, FALSE, TRUE)), !any(state_flag_vec(BCx, "mean", "cis")))

cat("state-conditioned checks passed\n")
