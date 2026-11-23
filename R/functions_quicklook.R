###############################################################
### functions_quicklook.R
### Functions of the quick-look angle classification (analysis/quicklook_angles.R). Point estimates with model-based standard errors and delta-method
### intervals stand in for the bootstrap, so the whole classification runs locally from the Section 1 checkpoint in minutes. Sourced after R/functions.R,
### whose fit_dispersion() and .cr_loglik() give the pipeline's estimates.
###############################################################

## QL.STEP: the step in log(phi) of the central second difference that gives the observed information of the Cox-Reid adjusted likelihood.
QL.STEP <- 0.01

## ql_fit_gene: one gene's fit in one group, with the variances of its log2 mean and log2 dispersion from the fit alone (handoff Section 3).
##   mean        mu = sum(y) / sum(O) (the pipeline's ratio estimate). Its variance is that of the summed counts under the NB, sum M_i (1 + phi M_i), over
##               (sum y)^2, since mu is proportional to sum(y); divided by (ln 2)^2 for the log2 scale.
##   dispersion  phi* maximizes the Cox-Reid adjusted likelihood; its observed information comes from the second difference of the adjusted likelihood in
##               log(phi) at phi*, I_log = -d2 l / d(log phi)^2 = phi*^2 (-d2 l / d phi^2) at a maximum. The split correction delta is held fixed, so the
##               corrected phi = phi* - delta has the standard error of phi*, and V(log2 phi) = 1 / ((ln 2)^2 phi^2 (-d2 l / d phi^2)) = phi*^2 / ((ln 2)^2 phi^2 I_log).
## status: "fitted" (phi > 0), "zero" (phi = 0 from the boundary rule or the truncation after the split correction), "failed" (no finite phi) or
## "no counts" (every count zero, so mu = 0). Returns c(mu, phi, phi_app, v_mu, v_phi, status code 1 to 4).
QL.STATUS <- c("fitted", "zero", "failed", "no counts")
ql_fit_gene <- function(g, mat, expo, eta) {
  stopifnot(is.character(g), length(g) == 1, is.numeric(expo), all(expo > 0))
  y <- as.numeric(mat[g, ]); f <- fit_dispersion(y, expo, eta, "NB")
  out <- c(mu = f[["mu"]], phi = f[["disp"]], phi_app = f[["disp_app"]], v_mu = NA_real_, v_phi = NA_real_, status = NA_real_)
  if (!is.finite(f[["mu"]]) || sum(y) == 0) { out[["status"]] <- 4; return(out) }
  M <- f[["mu"]] * expo
  ph <- if (is.finite(f[["disp"]])) f[["disp"]] else 0
  out[["v_mu"]] <- sum(M * (1 + ph * M)) / sum(y)^2 / log(2)^2
  if (!is.finite(f[["disp"]])) { out[["status"]] <- 3; return(out) }
  if (f[["disp"]] == 0) { out[["status"]] <- 2; return(out) }
  ll <- .cr_loglik(y, M); lp <- log(f[["disp_app"]]); h <- QL.STEP
  i_log <- -(ll(exp(lp + h)) - 2 * ll(exp(lp)) + ll(exp(lp - h))) / h^2
  out[["v_phi"]] <- if (is.finite(i_log) && i_log > 0) f[["disp_app"]]^2 / (log(2)^2 * f[["disp"]]^2 * i_log) else NA_real_
  out[["status"]] <- 1
  out
}

## ql_contrasts: the four log2 contrasts of one quantity from the per-group log2 values x and their variances v (genes x c(Pc, Pe, Hc, He) matrices), with
## the variances and covariances that the classification planes need (handoff Section 4). The four groups are independent but for Hc and He, which share
## the hybrid cells; their covariance is set to zero here. Delta_trans = Delta_total - Delta_cis on the same hybrid cells, so Cov(cis, trans) = -Var(cis).
## Each contrast gets a two-sided Wald test, z = Delta / SE, and BH q-values over the genes with a defined contrast. Returns a data frame by gene.
ql_contrasts <- function(x, v) {
  stopifnot(identical(colnames(x), c("Pc", "Pe", "Hc", "He")), identical(dim(x), dim(v)))
  tot <- x[, "Pc"] - x[, "Pe"]; cis <- x[, "Hc"] - x[, "He"]; trans <- tot - cis; dom <- 0.5 * (x[, "Hc"] + x[, "He"]) - 0.5 * (x[, "Pc"] + x[, "Pe"])
  vt <- v[, "Pc"] + v[, "Pe"]; vc <- v[, "Hc"] + v[, "He"]
  out <- data.frame(total = tot, cis = cis, trans = trans, dom = dom,
                    var_total = vt, var_cis = vc, var_trans = vt + vc, var_dom = (vt + vc) / 4,
                    cov_cis_trans = -vc, cov_total_dom = 0.5 * (v[, "Pe"] - v[, "Pc"]),
                    cov_trans_dom = 0.5 * (v[, "Pe"] - v[, "Pc"]) + 0.5 * (v[, "He"] - v[, "Hc"]), row.names = rownames(x))
  for (k in c("total", "cis", "trans", "dom")) {
    z <- out[[k]] / sqrt(out[[paste0("var_", k)]]); p <- 2 * pnorm(-abs(z)); p[!is.finite(z)] <- NA
    out[[paste0("se_", k)]] <- sqrt(out[[paste0("var_", k)]]); out[[paste0("q_", k)]] <- p.adjust(p, "BH")
  }
  out
}

## QL.SCHEMES: the three classification planes (handoff Section 5). Each orients a gene by s, the sign of a parental difference, so the higher parent (or the
## higher trans environment) sits at +45 degrees for every gene; ref and other name the two contrasts of the plane, which also form the Conserved gate.
QL.SCHEMES <- list(
  cis_trans   = list(gate = c("cis", "trans"), landmarks = c(-45, 45, 90), circular = TRUE),
  inheritance = list(gate = c("total", "dom"), landmarks = c(-45, 0, 45), circular = FALSE),
  trans_inh   = list(gate = c("trans", "dom"), landmarks = c(-45, 0, 45), circular = FALSE))

## ql_scheme_table: one quantity's genes placed in one plane: h >= 0 and v, the angle theta = atan2(v, h) in degrees (-90 to 90), the magnitude, SE(theta) by
## the delta method, the 95% interval theta +- 1.96 SE, the landmarks it covers and the class (handoff Section 6).
##   cis_trans    s = sign(total), h = s (cis + trans) / 2, v = s (cis - trans) / 2. +-90 degrees are one direction (compensatory, total zero), so the angle
##                lives on a circle of 180 degrees: interval ends wrap and coverage uses the circular distance.
##   inheritance  s = sign(total), h = |total| / 2, v = Delta_dom; the interval is clamped to [-90, 90], whose ends are different directions.
##   trans_inh    s = sign(trans), h = |trans| / 2, v = Delta_dom.
## The gate first: a gene whose two contrasts both have q >= sig is Conserved ("No trans divergence" for trans_inh). Every other gene takes the class of the
## one landmark its interval covers, or, covering none or several, the class of the region holding theta. bad gives a class that replaces the angle for a
## gene whose fits leave the quantity undefined ("Zero fit", "Failed fit", "No counts"; NA when defined). Returns one row per gene.
ql_scheme_table <- function(scheme, ct, bad, sig = 0.05) {
  sc <- QL.SCHEMES[[scheme]]
  stopifnot(!is.null(sc), is.data.frame(ct), length(bad) == nrow(ct))
  ref <- if (scheme == "trans_inh") "trans" else "total"
  s <- sign(ct[[ref]]); s[s == 0] <- 1
  if (scheme == "cis_trans") {
    h <- 0.5 * s * (ct$cis + ct$trans); v <- 0.5 * s * (ct$cis - ct$trans)
    vh <- 0.25 * (ct$var_cis + ct$var_trans + 2 * ct$cov_cis_trans); vv <- 0.25 * (ct$var_cis + ct$var_trans - 2 * ct$cov_cis_trans)
    chv <- 0.25 * (ct$var_cis - ct$var_trans)
  } else {
    h <- 0.5 * abs(ct[[ref]]); v <- ct$dom
    vh <- 0.25 * ct[[paste0("var_", ref)]]; vv <- ct$var_dom; chv <- 0.5 * s * ct[[paste0("cov_", ref, "_dom")]]
  }
  theta <- atan2(v, h) * 180 / pi; mag <- sqrt(h^2 + v^2)
  se <- (180 / pi) * sqrt(pmax(v^2 * vh + h^2 * vv - 2 * h * v * chv, 0)) / (h^2 + v^2)
  hw <- 1.96 * se
  wrap <- function(a) ((a + 90) %% 180) - 90
  lo <- if (sc$circular) wrap(theta - hw) else pmax(theta - hw, -90); hi <- if (sc$circular) wrap(theta + hw) else pmin(theta + hw, 90)
  cover <- sapply(sc$landmarks, function(L) if (sc$circular) abs(wrap(theta - L)) <= hw | hw >= 90 else (theta - hw) <= L & L <= (theta + hw))
  cover <- matrix(cover, ncol = length(sc$landmarks))
  covers <- apply(cover, 1, function(r) paste(sc$landmarks[r %in% TRUE], collapse = ";"))
  hi_p <- ifelse(s > 0, "Pc", "Pe"); lo_p <- ifelse(s > 0, "Pe", "Pc")
  if (scheme == "cis_trans") {
    lm_class <- c("Trans only", "Cis only", "Compensatory")
    region <- ifelse(abs(theta) == 90, "Compensatory", ifelse(theta > 45, "Cis, opposed by trans", ifelse(theta == 45, "Cis only",
              ifelse(theta >= 0, "Cis + trans, cis larger", ifelse(theta > -45, "Cis + trans, trans larger", ifelse(theta == -45, "Trans only", "Trans, opposed by cis"))))))
  } else {
    above <- if (scheme == "inheritance") "Overdominant" else "Beyond both parents, above"
    below <- if (scheme == "inheritance") "Underdominant" else "Beyond both parents, below"
    lm_mat <- cbind(paste(lo_p, "dominant"), "Additive", paste(hi_p, "dominant"))
    region <- ifelse(theta > 45, above, ifelse(theta == 45, paste(hi_p, "dominant"), ifelse(theta > 0, paste("Partial toward", hi_p),
              ifelse(theta == 0, "Additive", ifelse(theta > -45, paste("Partial toward", lo_p), ifelse(theta == -45, paste(lo_p, "dominant"), below))))))
  }
  n_cov <- rowSums(cover, na.rm = TRUE); one <- which(n_cov == 1 & is.finite(theta))
  cls <- region
  if (length(one)) {
    j <- apply(cover[one, , drop = FALSE], 1, which)
    cls[one] <- if (scheme == "cis_trans") lm_class[j] else lm_mat[cbind(one, j)]
  }
  gq <- cbind(ct[[paste0("q_", sc$gate[1])]], ct[[paste0("q_", sc$gate[2])]])
  conserved <- rowSums(gq < sig, na.rm = TRUE) == 0 & rowSums(is.finite(gq)) == 2
  cls[conserved] <- if (scheme == "trans_inh") "No trans divergence" else "Conserved"
  cls[!is.finite(theta) & !conserved] <- NA
  cls[!is.na(bad)] <- bad[!is.na(bad)]
  merged <- sub("^Cis \\+ trans, (cis|trans) larger$", "Cis + trans", cls)
  keep <- is.na(bad)
  data.frame(gene = rownames(ct), scheme = scheme, contrast_1 = sc$gate[1], est_1 = ct[[sc$gate[1]]], se_1 = ct[[paste0("se_", sc$gate[1])]],
             q_1 = ct[[paste0("q_", sc$gate[1])]], contrast_2 = sc$gate[2], est_2 = ct[[sc$gate[2]]], se_2 = ct[[paste0("se_", sc$gate[2])]],
             q_2 = ct[[paste0("q_", sc$gate[2])]], s = s, h = h, v = v, magnitude = mag, theta = theta, se_theta = se, theta_lo = lo, theta_hi = hi,
             landmarks_covered = covers, class = cls, class_merged = merged, high_parent = ifelse(keep, hi_p, NA), fit_flag = ifelse(keep, "", bad),
             row.names = NULL, stringsAsFactors = FALSE)
}
