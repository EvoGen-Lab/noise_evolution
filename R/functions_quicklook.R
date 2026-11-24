###############################################################
### functions_quicklook.R
### Functions of the quick-look classification by direction (analysis/quicklook_angles.R). Point estimates with model-based standard errors and normal or
### delta-method intervals stand in for the bootstrap, so the classification runs locally from the Section 1 checkpoint in minutes. Sourced after R/functions.R,
### whose fit_dispersion() and .cr_loglik() give the pipeline's estimates and whose dir_classify_table() gives the gate and the classes.
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

## ql_dir_rows: one quantity's rows of the direction table for one scheme from the quick look's per-group log2 values x and variances v (genes x c(Pc, Pe, Hc, He))
## on all cells, in the columns dir_classify_table() reads, so the quick look gates and classes genes with the pipeline's rules. The parts and diagonals are those
## of dir_parts(); their variances and the covariance of A and B follow from the group variances with the Hc and He estimates taken as uncorrelated (the quick
## look has no split, so A and B share hybrid cells and covary). The margin tests are the normal versions of dir_margin_test(): #{x* <= delta} / B becomes
## Phi((delta - x) / se), and so on. theta's SE is the delta method's, sqrt(B^2 Var A + A^2 Var B - 2 A B Cov) / (A^2 + B^2) in radians; its interval is theta
## +- z se with z = qnorm(1 - tail) (tail per gene, or one value). p is dir_share() with the model variances. xi (by gene, or NULL) centres the inheritance
## schemes on the shared shift (scheme names with "_xi"). zero (by gene) flags an observed zero rate or dispersion, failed a failed fit.
ql_dir_rows <- function(q, scheme, x, v, delta, zero, failed, xi = NULL, tail = 0.025) {
  stopifnot(identical(colnames(x), c("Pc", "Pe", "Hc", "He")), identical(dim(x), dim(v)), scheme %in% c(.DIR.SCHEMES, "inheritance_xi", "trans_inh_xi"))
  sc <- sub("_xi$", "", scheme)
  if (grepl("_xi$", scheme)) { x[, c("Hc", "He")] <- x[, c("Hc", "He")] - .DIR.Q[[q]]$sh * xi }
  a <- lapply(setNames(colnames(x), colnames(x)), function(k) x[, k])
  pt <- sapply(seq_len(nrow(x)), function(i) dir_parts(sc, lapply(a, `[`, i)))
  A <- pt["A_all", ]; B <- pt["B_all", ]; S <- pt["S", ]; D <- pt["D", ]
  Vpc <- v[, "Pc"]; Vpe <- v[, "Pe"]; Vhc <- v[, "Hc"]; Vhe <- v[, "He"]
  vv <- switch(sc,
    cis_trans   = list(A = Vhc + Vhe, B = Vpc + Vpe + Vhc + Vhe, C = -(Vhc + Vhe), S = Vpc + Vpe),
    inheritance = list(A = 0.25 * (Vhc + Vhe) + Vpe, B = Vpc + 0.25 * (Vhc + Vhe), C = -0.25 * (Vhc + Vhe), S = Vpc + Vpe),
    trans_inh   = list(A = Vhe + Vpe, B = Vpc + Vhc, C = 0 * Vpc, S = Vpc + Vpe + Vhc + Vhe))
  vD <- if (sc == "trans_inh") vv$A + vv$B else vv$A + vv$B - 2 * vv$C
  se <- list(A = sqrt(vv$A), B = sqrt(vv$B), S = sqrt(vv$S), D = sqrt(vD)); est <- list(A = A, B = B, S = S, D = D)
  out <- data.frame(gene = rownames(x), quantity = q, scheme = scheme, failed_fit = failed, zero_fit = zero, A = A, B = B, S = S, D = D, stringsAsFactors = FALSE)
  for (k in c("A", "B", "S", "D")) {
    out[[paste0(k, "_se")]] <- se[[k]]; out[[paste0(k, "_lo")]] <- est[[k]] - 1.96 * se[[k]]; out[[paste0(k, "_hi")]] <- est[[k]] + 1.96 * se[[k]]
    for (m in seq_along(DIR.DELTA.MULT)) {
      dl <- DIR.DELTA.MULT[m] * delta; z1 <- pnorm((dl - est[[k]]) / se[[k]]); z2 <- pnorm((est[[k]] + dl) / se[[k]])
      out[[sprintf("pR_%s_m%d", k, m)]] <- pmin(1, 2 * pmin(z1, z2))
      out[[sprintf("pE_%s_m%d", k, m)]] <- pmax(pnorm((est[[k]] - dl) / se[[k]]), pnorm((-dl - est[[k]]) / se[[k]]))
    }
  }
  bad <- failed | zero | !is.finite(A) | !is.finite(B)
  th <- ifelse(bad, NA_real_, dir_theta(A, B))
  se_th <- (180 / pi) * sqrt(pmax(B^2 * vv$A + A^2 * vv$B - 2 * A * B * vv$C, 0)) / (A^2 + B^2)
  z <- qnorm(1 - tail)
  out$theta <- th; out$theta_se <- se_th; out$theta_tail <- tail; out$theta_lo <- th - z * se_th; out$theta_hi <- th + z * se_th
  out$p_raw <- ifelse(bad, NA_real_, dir_share(A, B)); out$p <- ifelse(bad, NA_real_, dir_share(A, B, vv$A, vv$B))
  out$cov_AB <- vv$C; out$xi <- if (grepl("_xi$", scheme)) xi else NA_real_
  out[zero, grep("^p[RE]_", names(out))] <- NA     # infinite estimates leave the normal tests undefined: the row is labelled Zero fit
  out[failed, setdiff(names(out), c("gene", "quantity", "scheme", "failed_fit", "zero_fit"))] <- NA
  out
}

## ql_dir_combo: the rows of one quantity and scheme (ql_dir_rows()) from the list Q of per-group log2 values and variances by quantity, with the quantity's
## margin, zero and failed flags (the mean family or the dispersion family) and, with alpha (the false coverage-statement rate by "<quantity>|<scheme>"), theta
## intervals at level 1 - 2 alpha; without it, or when no gene is divergent, 95% intervals.
ql_dir_combo <- function(q, sc, Q, delta, zero, failed, xi, alpha = NULL) {
  f <- if (q == "mean") "mean" else "disp"
  a <- if (is.null(alpha)) NA_real_ else alpha[[paste(q, sc, sep = "|")]]
  ql_dir_rows(q, sc, Q[[q]]$x, Q[[q]]$v, delta[[q]], zero[[f]], failed[[f]], xi, if (is.finite(a) && a > 0) a else 0.025)
}
