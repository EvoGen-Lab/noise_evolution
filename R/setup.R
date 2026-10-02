###############################################################
### setup.R
### Session setup shared by analysis.R, power_analysis.R and
### external_validation.R: package check (renv), project paths and
### output folders, class names and color palettes, and the shared
### function library. Defines no functions. Each script loads it with
###   library(here); source(here("R", "setup.R"))
### and then loads its own packages and any analysis-specific
### function file.
###############################################################

# renv pins every package to the version recorded in renv.lock. Run
# renv::restore() once on a new machine
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::status()

# here() anchors every path to the project root so that functions and data load correctly
library(here)
stopifnot(file.exists(here("R", "functions.R")))
source(here("R", "functions.R"))

# data/ holds everything read from disk. cluster_inputs and cluster_outputs
# carry *_inputs.rda and *_output.rda to and from the SLURM cluster.
SC.DIR       <- here("data", "single_cell")
GENOME.DIR   <- here("data", "genomes")
EXTERNAL.DIR <- here("data", "external")
INPUT.DIR    <- here("data", "cluster_inputs")
OUTPUT.DIR   <- here("data", "cluster_outputs")

# results/ holds everything this script writes. checkpoints hold the
# section{N}_checkpoint.rda files, console holds a copy of what each section
# prints to the R console (section{N}_console.txt), and tables holds .csv files.
CHECKPOINT.DIR <- here("results", "checkpoints")
CONSOLE.DIR    <- here("results", "console")
TABLE.DIR      <- here("results", "tables")

# figures/ holds main text figures, extended data/supplement figures, and everything
# else (diagnostic and exploratory plots)
FIGURE.DIR <- here("figures")
for (d in c(INPUT.DIR, OUTPUT.DIR, CHECKPOINT.DIR, CONSOLE.DIR, TABLE.DIR,
            file.path(FIGURE.DIR, c("main", "extended", "extra")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Classes
REG.CLASS <- c("Conserved", "Cis", "Trans", "Cis + Trans", "Compensatory")
DOM.CLASS <- c("Conserved", "Sc.Dominant", "Se.Dominant", "Overdominant", "Underdominant", "Additive")

## Color palettes
# Species anchors
SPECIES.COLOR <- c(Sc = "#B5533C", Se = "#3E6B7A")
# Regulatory classes
COLOR.LIST.1 <- c("#D4D4CF","#2B2F42","#640B14","#E0A526","#5F8F5A"); names(COLOR.LIST.1) <- REG.CLASS
# Dominance classes
COLOR.LIST.2 <- c("#D4D4CF", SPECIES.COLOR[["Sc"]], SPECIES.COLOR[["Se"]], "#7A2E4E","#D9A7BF","#8C8C86"); names(COLOR.LIST.2) <- DOM.CLASS
# Diverging heatmap ramp (warm = enriched, cool = depleted
COLOR.LIST.3 <- colorRampPalette(c("#70453F","#A2623A","#C8893A","#DDB264","#F7F7F7","#D6ABD6","#BB82C3","#9D57AA","#6F3C79"))(11)
# Warm ramp with monotone lightness
SEQ.ANCHORS <- c("#F6F1E4","#E8C77A","#C98A3A","#8A4A3A","#33384A")
COLOR.SEQ <- colorRampPalette(SEQ.ANCHORS)(100)
# Cell-cycle phases
COLOR.PHASE <- c(G1 = "#A9CBF2", S = "#3F63BF", G2M = "#132057")
# Power-curve lines follow the same anchors as COLOR.SEQ, minus the cream
# end that vanishes on white
POWER.COLOR <- SEQ.ANCHORS[-1]
# Plum ramp for cluster labels
COLOR.CLUSTER <- c("#D7ACD7", "#A159AF", "#512B59")
# Neutral colors for thresholds and reference lines
COLOR.ACCENT <- "#3A3A3A"
COLOR.GREY <- c(light = "#DADAD5", mid = "#B0B0AB", dark = "#6E6E6A")

# Seed for the steps that carry no seed of their own: each section starts with
# set.seed(SEED.SECTION + n) for section n, so a section gives the same draws whether it runs
# from the top or is resumed from a checkpoint. Steps with their own seed use with_local_seed().
SEED.SECTION <- 100

# Record the package versions of this session beside the results and compare them with the records
# the cluster scripts write (pkg_versions_<script>.csv, copied back with each job's output). A
# difference is a warning, not a stop: check_pkg_versions(..., strict = TRUE) stops on a major or minor one.
write_pkg_versions("local", dir = dirname(CHECKPOINT.DIR))
check_pkg_versions(file.path(dirname(CHECKPOINT.DIR), "pkg_versions_local.csv"), OUTPUT.DIR)
