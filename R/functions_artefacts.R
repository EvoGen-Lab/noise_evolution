###############################################################
### functions_artefacts.R
### Functions of the tests of why the data depart from the NB (handoff "tests of why the data depart from the NB", tests T0 to T12),
### sourced after R/functions.R and R/functions_checks.R by analysis/checks/artefact_tests.R and cluster/scripts/artefact_null_grid.R.
### Groups are Pc, Pe, Hc and He; every fit is the pipeline's Cox-Reid NB fit (fit_dispersion()) with the group's share correction unless stated.
###############################################################

## ART.GROUPS: the four groups of every test, in the order of the Section 1 matrices.
ART.GROUPS <- c("Pc", "Pe", "Hc", "He")

## art_trim_idx: the cells of one fit left out by the trimming rule of the console run (estimator_pair_gene(trim = 0.005)): the top ceiling(frac x n) cells
## by scaled count y / O, at least one, ties broken by cell order.
art_trim_idx <- function(y, e, frac) order(y / e, decreasing = TRUE)[seq_len(max(1, ceiling(frac * length(y))))]

## art_gene_stats: the trimming, depth-stratum and thinning statistics of one gene in one group, from its counts y and offsets e (observed or simulated).
##   mu, mbar, phi        the Cox-Reid fit on all cells (phi = 0 is the Poisson limit);
##   log2_trim            log2(phi without the trimmed cells / phi), -Inf when the trimmed phi is 0 (kept and ranked), NA when phi is 0 (T0a);
##   excess_frac          the share of the excess sum((y - M)^2 - M) over all cells, M from the all-cells fit, carried by the trimmed cells; NA when the
##                        total excess is not positive (T0b);
##   strat1..3, mstrat1..3  log2(phi in depth tercile s / phi) and the expected count per cell in that tercile (terc gives each cell's tercile; T3);
##   thin_log2phi, thin_mu_ratio  log2(phi of the thinned counts / phi) and the ratio of the means, when y_thin and e_thin are given (T7a).
## The trimmed cell indices are returned as attribute "trim".
art_gene_stats <- function(y, e, eta, frac, terc, y_thin = NULL, e_thin = NULL) {
  f <- fit_dispersion(y, e, eta); phi <- f[["disp"]]; mu <- f[["mu"]]
  tr <- art_trim_idx(y, e, frac); ft <- fit_dispersion(y[-tr], e[-tr], eta)
  l2 <- function(a, b) if (!is.finite(a) || !is.finite(b) || b <= 0) NA_real_ else if (a <= 0) -Inf else log2(a / b)
  M <- mu * e; ex <- (y - M)^2 - M
  st <- vapply(1:3, function(s) l2(fit_dispersion(y[terc == s], e[terc == s], eta)[["disp"]], phi), numeric(1))
  out <- c(mu = mu, mbar = mu * mean(e), phi = phi, log2_trim = l2(ft[["disp"]], phi),
           excess_frac = if (sum(ex) > 0) sum(ex[tr]) / sum(ex) else NA_real_,
           strat1 = st[1], strat2 = st[2], strat3 = st[3], mstrat1 = mu * mean(e[terc == 1]), mstrat2 = mu * mean(e[terc == 2]), mstrat3 = mu * mean(e[terc == 3]),
           thin_log2phi = NA_real_, thin_mu_ratio = NA_real_)
  if (!is.null(y_thin)) { fh <- fit_dispersion(y_thin, e_thin, eta); out[["thin_log2phi"]] <- l2(fh[["disp"]], phi); out[["thin_mu_ratio"]] <- fh[["mu"]] / mu }
  attr(out, "trim") <- tr
  out
}

## art_nb_models: T4, three NB models of one gene in one group fitted by maximum likelihood with phi free: (a) E[y] = O mu, (b) E[y] = mu O^beta,
## (c) E[y] = mu (O + rho), an additive component lambda = mu rho outside the offset (rho = lambda / mu at unit offset). log phi is bounded to the
## pipeline's search range and log rho to [-15, 5]. Returns beta and rho with SEs from the Hessian, the three phi and log-likelihoods.
art_nb_models <- function(y, e) {
  na <- c(beta = NA, beta_se = NA, rho = NA, rho_se = NA, phi_a = NA, phi_b = NA, phi_c = NA, ll_a = NA, ll_b = NA, ll_c = NA)
  if (sum(y) == 0) return(na)
  lo <- log(e); m0 <- log(sum(y) / sum(e)); lp0 <- log(max(fit_dispersion(y, e)[["disp"]], 1e-3)); lpr <- -rev(LOG.THETA.RANGE)
  nll <- function(m, lphi) -sum(dnbinom(y, size = exp(-lphi), mu = m, log = TRUE))
  fa <- optim(c(m0, lp0), function(p) nll(exp(p[1]) * e, p[2]), method = "L-BFGS-B", lower = c(-Inf, lpr[1]), upper = c(Inf, lpr[2]))
  fb <- optim(c(m0, 1, lp0), function(p) nll(exp(p[1] + p[2] * lo), p[3]), method = "L-BFGS-B", lower = c(-Inf, -3, lpr[1]), upper = c(Inf, 5, lpr[2]))
  fc <- optim(c(m0, -3, lp0), function(p) nll(exp(p[1]) * (e + exp(p[2])), p[3]), method = "L-BFGS-B", lower = c(-Inf, -15, lpr[1]), upper = c(Inf, 5, lpr[2]))
  se <- function(f, fn, i) { h <- tryCatch(optimHess(f$par, fn), error = function(err) NULL); v <- if (is.null(h)) NULL else tryCatch(solve(h), error = function(err) NULL)
    if (is.null(v) || v[i, i] <= 0) NA_real_ else sqrt(v[i, i]) }
  rse <- se(fc, function(p) nll(exp(p[1]) * (e + exp(p[2])), p[3]), 2)
  c(beta = fb$par[2], beta_se = se(fb, function(p) nll(exp(p[1] + p[2] * lo), p[3]), 2), rho = exp(fc$par[2]), rho_se = exp(fc$par[2]) * rse,
    phi_a = exp(fa$par[2]), phi_b = exp(fb$par[3]), phi_c = exp(fc$par[3]), ll_a = -fa$value, ll_b = -fb$value, ll_c = -fc$value)
}

## art_influence: T11, the leave-one-out influence of the k cells with the largest scaled counts on phi of one gene in one group, for the Cox-Reid and the
## moment estimators. The influence of cell i is phi - phi without cell i; reported relative to phi: the largest single influence and the sum of the five
## largest. NA when phi is 0.
art_influence <- function(y, e, eta, k = 20) {
  top <- order(y / e, decreasing = TRUE)[seq_len(min(k, length(y)))]
  out <- c()
  for (m in c("NB", "weighted")) {
    phi <- fit_dispersion(y, e, eta, m)[["disp"]]
    inf <- vapply(top, function(i) phi - fit_dispersion(y[-i], e[-i], eta, m)[["disp"]], numeric(1))
    tag <- if (m == "NB") "cr" else "mom"
    out[paste0(c("max_share_", "top5_share_"), tag)] <- if (phi > 0) c(max(inf) / phi, sum(sort(inf, decreasing = TRUE)[1:min(5, length(inf))]) / phi) else c(NA, NA)
  }
  out
}

## art_null_chunk: T0 and the null parts of T1, T3, T4, T7a and T11 for a chunk of genes. For every gene and group, R counts are simulated as gamma NB
## with the gene's all-cells mu and Cox-Reid phi on the group's own offsets (Poisson when phi = 0), each gene from its own random stream, and every
## statistic of art_gene_stats() is computed on each replicate (thinned replicates keep each UMI with probability thin_p on the same offsets). The first
## R4 replicates also get art_nb_models() and the first R11 art_influence(). Returns
##   stats   one row per gene and group: the observed statistics, the null mean and 5/50/95% quantiles of log2_trim and excess_frac, the observed values'
##           upper-tail percentile in the null (share of null values >= observed), and the null medians of the stratum and thinning statistics;
##   burden  per group an R x cells matrix, the number of genes of the chunk for which each cell is trimmed in each replicate;
##   obs_burden  per group the observed count of genes for which each cell is trimmed;
##   t4, t11  the null rows of art_nb_models() and art_influence().
art_null_chunk <- function(genes, mats, expos, eta, frac, terc, R, R4, R11, thin_p, streams) {
  stats <- list(); t4 <- list(); t11 <- list()
  burden <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) matrix(0L, R, length(expos[[k]])))
  obs_burden <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) integer(length(expos[[k]])))
  for (g in genes) with_stream(streams[[g]], for (k in ART.GROUPS) {
    y <- as.numeric(mats[[k]][g, ]); e <- expos[[k]]; n <- length(e)
    ob <- art_gene_stats(y, e, eta[[k]], frac, terc[[k]], rbinom(n, y, thin_p), e)
    obs_burden[[k]][attr(ob, "trim")] <- obs_burden[[k]][attr(ob, "trim")] + 1L
    mu <- ob[["mu"]]; phi <- ob[["phi"]]
    if (!is.finite(mu) || mu <= 0) next
    sims <- matrix(NA_real_, R, length(ob)); colnames(sims) <- names(ob)
    for (r in seq_len(R)) {
      ys <- rpois(n, e * mu * rgamma1(n, if (is.finite(phi)) phi else 0))
      st <- art_gene_stats(ys, e, eta[[k]], frac, terc[[k]], rbinom(n, ys, thin_p), e)
      sims[r, ] <- st; burden[[k]][r, attr(st, "trim")] <- burden[[k]][r, attr(st, "trim")] + 1L
      if (r <= R4) t4[[length(t4) + 1]] <- c(group = match(k, ART.GROUPS), mbar = ob[["mbar"]], art_nb_models(ys, e))
      if (r <= R11) t11[[length(t11) + 1]] <- c(group = match(k, ART.GROUPS), mbar = ob[["mbar"]], art_influence(ys, e, eta[[k]]))
    }
    q <- function(v) { v <- v[!is.na(v)]; if (length(v) < 5) rep(NA_real_, 4) else c(mean(v[is.finite(v)]), quantile(v, c(0.05, 0.5, 0.95), type = 1, names = FALSE)) }
    pc <- function(o, v) { v <- v[!is.na(v)]; if (!is.finite(o) && !is.infinite(o) || !length(v)) NA_real_ else mean(v >= o) }
    stats[[length(stats) + 1]] <- data.frame(gene = g, group = k, t(ob),
      null_trim_mean = q(sims[, "log2_trim"])[1], null_trim_q05 = q(sims[, "log2_trim"])[2], null_trim_q50 = q(sims[, "log2_trim"])[3], null_trim_q95 = q(sims[, "log2_trim"])[4],
      null_exf_mean = q(sims[, "excess_frac"])[1], null_exf_q05 = q(sims[, "excess_frac"])[2], null_exf_q50 = q(sims[, "excess_frac"])[3], null_exf_q95 = q(sims[, "excess_frac"])[4],
      pct_trim = pc(ob[["log2_trim"]], sims[, "log2_trim"]), pct_exf = pc(ob[["excess_frac"]], sims[, "excess_frac"]),
      null_strat1 = median(sims[, "strat1"], na.rm = TRUE), null_strat2 = median(sims[, "strat2"], na.rm = TRUE), null_strat3 = median(sims[, "strat3"], na.rm = TRUE),
      null_thin_log2phi = median(sims[, "thin_log2phi"], na.rm = TRUE), null_thin_mu_ratio = median(sims[, "thin_mu_ratio"], na.rm = TRUE), row.names = NULL)
  })
  list(stats = do.call(rbind, stats), burden = burden, obs_burden = obs_burden,
       t4 = if (length(t4)) as.data.frame(do.call(rbind, t4)) else NULL, t11 = if (length(t11)) as.data.frame(do.call(rbind, t11)) else NULL)
}

## art_sim_dataset: T8 and T12, one simulated dataset of every gene and group from the observed fits (fits: genes x columns <group>.mu, <group>.phi of the
## all-cells Cox-Reid fits; hybrid allele means include their abundance factor), the observed offsets and cell counts, as gamma NB. scenario:
##   A           every group of a gene at one dispersion, the geometric mean of its positive fitted values (true Delta phi = 0, true Delta mu observed);
##   B / S0      every group at its fitted dispersion;
##   S1 rho      B plus an additive Poisson component of mean rho mu per cell (outside the offset);
##   S2 beta     B with mean mu O^beta;
##   S3 q, J     B with a fraction q of cells (the same cells for every gene of the group; the two hybrid alleles share them) whose expression factor is
##               multiplied by J, the others scaled so the mean is unchanged;
##   S4 g        B with every UMI replaced by 1 + Geometric(mean g) copies;
##   S5 f        B plus Poisson counts of mean f O mu in the 2% of cells with the lowest offsets (the hybrid alleles share them).
## Draws from the given random stream. Returns the four count matrices (genes x cells).
art_sim_dataset <- function(fits, expos, genes, scenario, par1 = NA, par2 = NA, stream) with_stream(stream, {
  nH <- length(expos$Hc)
  jack <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) if (scenario == "S3") sample.int(length(expos[[k]]), max(1, round(par1 * length(expos[[k]])))) else integer(0))
  low <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) if (scenario == "S5") order(expos[[k]])[seq_len(max(1, round(0.02 * length(expos[[k]]))))] else integer(0))
  jack$He <- jack$Hc; low$He <- low$Hc
  phiA <- apply(fits[genes, paste0(ART.GROUPS, ".phi"), drop = FALSE], 1, function(p) { p <- p[is.finite(p) & p > 0]; if (length(p)) exp(mean(log(p))) else 0 })
  out <- list()
  for (k in ART.GROUPS) {
    e <- expos[[k]]; n <- length(e)
    m <- matrix(0L, length(genes), n, dimnames = list(genes, NULL))
    for (j in seq_along(genes)) {
      mu <- fits[genes[j], paste0(k, ".mu")]; ph <- if (scenario == "A") phiA[[j]] else fits[genes[j], paste0(k, ".phi")]
      if (!is.finite(mu) || mu <= 0) next
      a <- rgamma1(n, if (is.finite(ph)) ph else 0)
      if (scenario == "S3") { a <- a / (1 + par1 * (par2 - 1)); a[jack[[k]]] <- a[jack[[k]]] * par2 }
      lam <- if (scenario == "S2") mu * e^par1 * a else mu * e * a
      y <- rpois(n, lam)
      if (scenario == "S1") y <- y + rpois(n, par1 * mu)
      if (scenario == "S4") { i <- y > 0; y[i] <- y[i] + rnbinom(sum(i), size = y[i], prob = 1 / (1 + par1)) }     # y UMIs each with Geometric(mean g) extra copies
      if (scenario == "S5") y[low[[k]]] <- y[low[[k]]] + rpois(length(low[[k]]), par1 * e[low[[k]]] * mu)
      m[j, ] <- y
    }
    out[[k]] <- m
  }
  out
})

## art_gene_core: the per-gene fits that the signature of a dataset is built from, for one gene: in every group the all-cells Cox-Reid fit (mu, mbar, phi)
## with bound_pos (1 when the lower likelihood bound on phi is above zero: phi > 0 and the boundary LRT above BOUND.CRIT, dispersion_bound_stats(); this
## needs no bootstrap and stands in for the pipeline's resolution, which these Section 1 checks cannot run), the trimming statistic and excess share (art_gene_stats() without strata or thinning), the linear term a of L6 with its SE,
## and the disjoint arm (mean on half A, phi on half B). Returns a named vector group.quantity.
art_gene_core <- function(g, mats, expos, eta, halves, frac) {
  out <- c()
  for (k in ART.GROUPS) {
    y <- as.numeric(mats[[k]][g, ]); e <- expos[[k]]; A <- halves[[k]]; B <- setdiff(seq_along(e), A)
    st <- art_gene_stats(y, e, eta[[k]], frac, rep(1:3, length.out = length(e)))
    cl <- as.numeric(dispersion_bound_stats(y, e, eta[[k]])[["disp_lower"]] > 0)
    le <- linear_excess_gene(1, matrix(y, 1), e)
    out[paste(k, c("mu", "mbar", "phi", "bound_pos", "log2_trim", "excess_frac", "a", "a_se", "mu_A", "phi_B"), sep = ".")] <-
      c(st[c("mu", "mbar", "phi")], cl, st[c("log2_trim", "excess_frac")], le[c("a", "a_se")], sum(y[A]) / sum(e[A]), fit_dispersion(y[B], e[B], eta[[k]])[["disp"]])
  }
  out
}

## art_signature: the signature of one dataset (observed or simulated) from art_gene_core() rows (core, genes x columns): by group and mean-count bin the
## median trimming statistic and excess share and the median a; per group the Spearman correlation of a with the mean and the trend of log2 phi on log2 mu
## across genes (phi > 0); and the own slopes of Delta log2 phi on Delta log2 mu for cis and trans in the same-cell and disjoint arms, on genes whose groups
## all have phi > 0 and a lower likelihood bound above zero (bound_pos, the bootstrap-free stand-in for the pipeline's resolution), with a gene-bootstrap interval (B draws, indices drawn first).
art_signature <- function(core, cbar, mbar_breaks, B = 200, seed = 1) {
  bins <- list(); grp <- list()
  for (k in ART.GROUPS) {
    v <- function(q) core[, paste0(k, ".", q)]
    mb <- cut(v("mbar"), mbar_breaks, right = FALSE)
    for (l in levels(mb)) { i <- mb == l & is.finite(v("mbar"))
      if (sum(i) < 10) next
      bins[[paste(k, l)]] <- data.frame(group = k, mean_count = l, genes = sum(i), median_log2_trim = median(v("log2_trim")[i], na.rm = TRUE),
                                        median_excess_frac = median(v("excess_frac")[i], na.rm = TRUE), median_a = median(v("a")[i], na.rm = TRUE), row.names = NULL) }
    ok <- is.finite(v("phi")) & v("phi") > 0 & v("mu") > 0
    grp[[k]] <- data.frame(group = k, spearman_a_mean = suppressWarnings(cor(v("a"), v("mbar"), method = "spearman", use = "complete.obs")),
                           trend_log2phi_on_log2mu = unname(coef(lm(log2(v("phi")[ok]) ~ log2(v("mu")[ok])))[2]), row.names = NULL)
  }
  lmu <- function(k, q) { m <- core[, paste0(k, ".", q)]; if (k == "Hc") m / (2 * cbar) else if (k == "He") m / (2 * (1 - cbar)) else m }
  slopes <- list()
  for (arm in c("same", "disjoint")) {
    qm <- if (arm == "same") "mu" else "mu_A"; qp <- if (arm == "same") "phi" else "phi_B"
    L <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) log2(lmu(k, qm)))
    P <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) ifelse(core[, paste0(k, ".", qp)] > 0 & core[, paste0(k, ".bound_pos")] %in% 1, log2(core[, paste0(k, ".", qp)]), NA_real_))
    for (md in c("cis", "trans")) {
      x <- if (md == "cis") L$Hc - L$He else (L$Pc - L$Pe) - (L$Hc - L$He); y <- if (md == "cis") P$Hc - P$He else (P$Pc - P$Pe) - (P$Hc - P$He)
      ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
      if (length(x) < 10) { slopes[[paste(arm, md)]] <- data.frame(arm = arm, contrast = md, genes = length(x), slope = NA, lo = NA, hi = NA); next }
      idx <- with_local_seed(seed, lapply(seq_len(B), function(b) sample.int(length(x), length(x), TRUE)))
      bs <- vapply(idx, function(i) unname(coef(lm(y[i] ~ x[i]))[2]), numeric(1))
      slopes[[paste(arm, md)]] <- data.frame(arm = arm, contrast = md, genes = length(x), slope = unname(coef(lm(y ~ x))[2]),
                                             lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)), row.names = NULL)
    }
  }
  list(bins = do.call(rbind, bins), groups = do.call(rbind, grp), slopes = do.call(rbind, slopes))
}

## art_scenario_task: one task of the scenario array (T8, T12): simulate a dataset (art_sim_dataset()), fit every gene (art_gene_core()) and return its
## signature (art_signature()) with the scenario labels.
art_scenario_task <- function(i, grid, fits, expos, eta, halves, genes, frac, cbar, mbar_breaks, streams, cores = 1) {
  t <- grid[i, ]
  mats <- art_sim_dataset(fits, expos, genes, t$scenario, t$par1, t$par2, streams[[i]])
  core <- do.call(rbind, parallel::mclapply(genes, art_gene_core, mats = mats, expos = expos, eta = eta, halves = halves, frac = frac, mc.cores = cores))
  rownames(core) <- genes
  sig <- art_signature(core, cbar, mbar_breaks, seed = i)
  lapply(sig, function(d) cbind(task = i, scenario = t$scenario, par1 = t$par1, par2 = t$par2, rep = t$rep, d))
}

## art_spearman_ci: T1, the Spearman correlation of a cell covariate with the excess burden, with a 95% interval from B resamples of cells (indices
## drawn first). NA when fewer than 10 cells have both values.
art_spearman_ci <- function(x, z, B = 500, seed = 1) {
  ok <- is.finite(x) & is.finite(z); x <- x[ok]; z <- z[ok]
  if (length(x) < 10) return(c(cells = length(x), rho = NA, lo = NA, hi = NA))
  idx <- with_local_seed(seed, lapply(seq_len(B), function(b) sample.int(length(x), length(x), TRUE)))
  bs <- vapply(idx, function(i) suppressWarnings(cor(x[i], z[i], method = "spearman")), numeric(1))
  c(cells = length(x), rho = cor(x, z, method = "spearman"), lo = unname(quantile(bs, 0.025, na.rm = TRUE)), hi = unname(quantile(bs, 0.975, na.rm = TRUE)))
}

## art_jaccard_mean: T10, the mean pairwise Jaccard index of a list of cell sets (the trimmed cells of each gene of a gene set). Genes that share the
## cells carrying their extreme counts give a high value; unrelated genes give about the overlap expected from cell depth alone.
art_jaccard_mean <- function(sets) {
  n <- length(sets); if (n < 2) return(NA_real_)
  pr <- utils::combn(n, 2)
  mean(vapply(seq_len(ncol(pr)), function(j) { a <- sets[[pr[1, j]]]; b <- sets[[pr[2, j]]]; length(intersect(a, b)) / length(union(a, b)) }, numeric(1)))
}
