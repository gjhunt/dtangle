set.seed(3151990)
library("dtangle.data")
source("analyze_scripts/save_sim_analysis.R", chdir = TRUE)

load("./newman_pbmc.rda")
load("./newman_fl.rda")
dsets <- list(`Newman PBMC` = newman_pbmc, `Newman FL` = newman_fl)

# Analysis for main paper
deconv_params <- list(Quantile = 0, Marker.Method = "p.value", all_markers = TRUE)
paper_out <- save_sim_analysis("newman", deconv_params = deconv_params, Dataset = dsets)
