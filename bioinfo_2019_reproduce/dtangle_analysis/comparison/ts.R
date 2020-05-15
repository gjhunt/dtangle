set.seed(3151995)
library("dtangle.data")
source("analyze_scripts/save_sim_analysis.R", chdir = TRUE)

# Time Sensitivity
deconv_params <- list(Quantile = c(0.75, 0.85, 0.9, 0.95, 0.99), Marker.Method = c("p.value"))
t_sens_out <- save_sim_analysis("t_sens", deconv_params = deconv_params, )
