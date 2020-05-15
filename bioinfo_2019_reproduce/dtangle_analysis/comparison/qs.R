set.seed(3151992)
library("dtangle.data")
source("analyze_scripts/save_sim_analysis.R", chdir = TRUE)

# Quantile Sensitivity
deconv_params <- list(Quantile = seq(0.85, 1, 0.001), Marker.Method = c("p.value", 
    "Ratio"))
q_sens_out <- save_sim_analysis("q_sens", deconv_params = deconv_params, Deconv.Method = c("dtangle", 
    "Q Prog", "LS Fit", "CIBERSORT", "EPIC"))
