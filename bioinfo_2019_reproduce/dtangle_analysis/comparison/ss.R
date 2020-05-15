set.seed(3151993)
library("dtangle.data")
source("analyze_scripts/save_sim_analysis.R", chdir = TRUE)

# Suppl Slope Analysis
deconv_params <- list(Quantile = c(0.85, 0.9, 0.95, 0.99), Marker.Method = c("p.value", 
    "Ratio"), gamma = seq(0.25, 2, 0.05))
slope_sens_out <- save_sim_analysis("supplement_s", deconv_params = deconv_params, 
    Deconv.Method = "dtangle")
