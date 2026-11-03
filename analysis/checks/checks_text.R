###############################################################
### checks_text.R
### Text audit of the methods document (docs/methods.md): G2, every DOI of the reference list resolved through Crossref and compared
### with its citation; G3, the displayed equations for the variance cap of the ratio mean, the balanced partition quadratic and its
### two limits checked against their closed forms; and REPORT.md, the table of every check ID with its headline numbers, written from
### results/checks/summary.csv and the mismatches. G1 (the quoted numbers) is in checks_data.R. G2 needs network access to
### api.crossref.org; without it the check is recorded as "not run".
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))

CHECK.SCRIPT <- "analysis/checks/checks_text.R"
SEED.CHECKS  <- 4026
set.seed(SEED.SECTION + 22)

##############################################################################
## G2. References: each DOI resolved and compared with its citation
# Each reference line reads "Authors (year) Title. Journal volume:pages. https://doi.org/<doi>"; the parse takes the authors before the year, the
# title up to the first period after it, the volume and first page from "volume:page", and the DOI from the link.
REF.LINES <- local({
  x <- readLines(here("docs", "methods.md"), encoding = "UTF-8")
  x <- x[(grep("^## References", x) + 1):(grep("^## Notes for revision", x) - 1)]
  x[grepl("https://doi.org/", x)]
})
REFS <- data.frame(authors = sub(" \\((19|20)[0-9]{2}\\).*", "", REF.LINES), year = sub("^.*?\\(((19|20)[0-9]{2})\\).*$", "\\1", REF.LINES),
                   title = sub("^.*?\\((19|20)[0-9]{2}\\) ([^.]+)\\..*$", "\\2", REF.LINES),
                   volume = ifelse(grepl(" [0-9]+:[0-9e]+", REF.LINES), sub("^.* ([0-9]+):([0-9e]+)[^0-9].*$", "\\1", REF.LINES), NA),
                   page = ifelse(grepl(" [0-9]+:[0-9e]+", REF.LINES), sub("^.* [0-9]+:(e?[0-9]+)[^0-9].*$", "\\1", REF.LINES), NA),
                   doi = sub("^.*https://doi.org/([^ ]+?)\\.?$", "\\1", REF.LINES), stringsAsFactors = FALSE)
REC <- do.call(rbind, lapply(REFS$doi, crossref_record))
if (!any(REC$retrieved)) {
  check_save("G2", list(citations = REFS, records = REC), "not run", headline = sprintf("no Crossref record could be read (%s); %d DOIs parsed", REC$error[1], nrow(REFS)), script = CHECK.SCRIPT)
} else {
  G2 <- do.call(rbind, lapply(seq_len(nrow(REFS)), function(i) reference_match(REFS[i, ], REC[i, ])))
  bad <- G2[!(G2$retrieved & G2$year_ok & G2$volume_ok & G2$page_ok & G2$authors_ok & G2$title_word_share >= 0.8), ]
  check_save("G2", G2, "run", headline = c(sprintf("%d of %d DOIs resolved; %d with a field that differs from the citation", sum(G2$retrieved), nrow(G2), nrow(bad)),
                                           sprintf("%s: year %s, volume %s, page %s, authors %s, title word share %.2f", bad$doi, bad$year_ok, bad$volume_ok, bad$page_ok, bad$authors_ok, bad$title_word_share)),
             mc_error = "exact", script = CHECK.SCRIPT)
}

##############################################################################
## G3. Displayed equations against their closed forms
# (a) The cap: N sum O^2 / (sum O)^2 = 1 + CV^2(O) holds exactly with the population CV (divisor N); with the sample CV (divisor N - 1) it does not.
# The ratio sum M (1 + phi M) sum M / (1 + phi M) / (sum M)^2 lies between 1 and the cap for every M and phi (Chebyshev's sum inequality).
# (b) The balanced partition: the minus root f = [(2 + gamma) - sqrt(gamma^2 + 4)] / (2 gamma) solves gamma f^2 - (2 + gamma) f + 1 = 0 and
# 1 / f = gamma + 1 / (1 - f), equals fstar_from_r(), and the plus root exceeds one. (c) The limits: f - (1/2 - gamma/8) is of order gamma^3
# (the next term of the expansion is gamma^3 / 128; read from gamma = 0.01 up, where rounding stays below that term), and gamma f - 1 is of order
# 1 / gamma (its leading term is -1 / gamma).
G3 <- with_local_seed(SEED.CHECKS, local({
  n_cap <- 10000; viol <- 0; cap_pop <- 0; cap_samp <- 0
  for (i in seq_len(n_cap)) {
    N <- sample(c(5, 50, 800), 1); O <- exp(rnorm(N, 0, runif(1, 0, 1.5))); phi <- exp(rnorm(1, 0, 2)); M <- O * exp(rnorm(1))
    r <- sum(M * (1 + phi * M)) * sum(M / (1 + phi * M)) / sum(M)^2; cap <- N * sum(O^2) / sum(O)^2
    viol <- viol + (r < 1 - 1e-12 || r > cap * (1 + 1e-12))
    cap_pop <- max(cap_pop, abs(cap - (1 + mean((O - mean(O))^2) / mean(O)^2))); cap_samp <- max(cap_samp, abs(cap - (1 + var(O) / mean(O)^2)))
  }
  g <- 10^seq(-4, 4, by = 0.25); f <- ((2 + g) - sqrt(g^2 + 4)) / (2 * g); fp <- ((2 + g) + sqrt(g^2 + 4)) / (2 * g)
  small <- g[g <= 0.1 & g >= 0.01]; fs <- ((2 + small) - sqrt(small^2 + 4)) / (2 * small); large <- g[g >= 10]; fl <- ((2 + large) - sqrt(large^2 + 4)) / (2 * large)
  list(cap_violations = viol, cap_cases = n_cap, cap_identity_err_population_cv = cap_pop, cap_identity_err_sample_cv = cap_samp,
       quad_residual = max(abs(g * f^2 - (2 + g) * f + 1)), balance_residual = max(abs(1 / f - g - 1 / (1 - f))), code_vs_doc = max(abs(fstar_from_r(g) - f)),
       plus_root_min = min(fp), small_limit = data.frame(gamma = small, err = fs - (0.5 - small / 8), err_over_gamma3 = (fs - (0.5 - small / 8)) / small^3),
       large_limit = data.frame(gamma = large, gamma_f_minus_1 = large * fl - 1, times_gamma = (large * fl - 1) * large))
}))
check_save("G3", G3, "code",
           headline = c(sprintf("(a) ratio outside [1, cap] in %d of %d random cases; cap identity error with the population CV %.2g, with the sample CV %.2g",
                                G3$cap_violations, G3$cap_cases, G3$cap_identity_err_population_cv, G3$cap_identity_err_sample_cv),
                        sprintf("(b) quadratic residual %.2g, balance residual %.2g, fstar_from_r() against the formula %.2g, smallest plus root %.4f", G3$quad_residual, G3$balance_residual, G3$code_vs_doc, G3$plus_root_min),
                        sprintf("(c) small gamma: (f - (1/2 - gamma/8)) / gamma^3 from %.5f to %.5f (1/128 = %.5f); large gamma: (gamma f - 1) gamma from %.4f to %.4f",
                                min(G3$small_limit$err_over_gamma3), max(G3$small_limit$err_over_gamma3), 1 / 128, min(G3$large_limit$times_gamma), max(G3$large_limit$times_gamma))),
           mc_error = "exact arithmetic to rounding", script = CHECK.SCRIPT)

##############################################################################
## REPORT.md: one row per check with its status and headline numbers, the checks whose numbers fall outside their criteria, the code and document
## mismatches, the checks that need the cluster or the data, and what cannot be tested with these data. Written last, from summary.csv.
SUMMARY <- read.csv(file.path(CHECK.DIR, "summary.csv"), stringsAsFactors = FALSE)
REPORT <- c("# Validation checks of the methods document", "",
            sprintf("Generated %s from results/checks/summary.csv by analysis/checks/checks_text.R. Headline numbers only; each check's full result is in results/checks/<ID>.rds.",
                    format(Sys.time(), "%Y-%m-%d %H:%M")), "",
            "| ID | Status | Headline numbers | Monte Carlo error | Script |", "|---|---|---|---|---|",
            sprintf("| %s | %s | %s | %s | %s |", SUMMARY$ID, SUMMARY$status, gsub("\\|", "/", gsub("; ", "<br>", SUMMARY$headline)), SUMMARY$mc_error, SUMMARY$script), "",
            "Not testable with these data (handoff section 7): species-specific ploidy effects (need diploid parental strains), capture differences between parental and",
            "hybrid cells, and the bursty-limit reading of theta and zeta.", "")
if (file.exists(file.path(CHECK.DIR, "mismatches.md"))) REPORT <- c(REPORT, readLines(file.path(CHECK.DIR, "mismatches.md")))
writeLines(REPORT, file.path(CHECK.DIR, "REPORT.md"))
cat("Wrote", file.path(CHECK.DIR, "REPORT.md"), "\n")
