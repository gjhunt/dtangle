set.seed(3151990)
library("dtangle.data")
source("analyze_scripts/save_sim_analysis.R", chdir = TRUE)

# Analysis for main paper
deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 0)
paper_out <- save_sim_analysis("paper", deconv_params = deconv_params)
