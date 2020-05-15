set.seed(3151991)

source("simulators.R")

for (err_struct in c("gaussian", "poisson")) {
    cat(err_struct, "\n")
    # Low Error
    cat("Low Error\n")
    deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, Scale = c("linear", 
        "log"))
    sim <- gen_sim(eps_sigma = 0.025)
    lowe_sim <- save_sim_analysis(paste0("sim_", err_struct, "_lowe"), deconv_params = deconv_params, 
        dset_table = NULL, Dataset = list(sim[[paste0("dset_", err_struct)]]), Deconv.Method = c("dtangle", 
            "linearRegression", "CIBERSORT", "EPIC", "LS Fit", "Q Prog"), outf = "./Analysis")
    
    # High Error
    cat("High Error\n")
    deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, Scale = c("linear", 
        "log"))
    sim <- gen_sim(eps_sigma = 0.75)
    highe_sim <- save_sim_analysis(paste0("sim_", err_struct, "_highe"), deconv_params = deconv_params, 
        dset_table = NULL, Dataset = list(sim[[paste0("dset_", err_struct)]]), Deconv.Method = c("dtangle", 
            "linearRegression", "CIBERSORT", "EPIC", "Q Prog", "LS Fit"), outf = "./Analysis")
    
    # Outliers
    cat("Outliers\n")
    deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, Scale = c("linear", 
        "log"))
    sim <- gen_sim(eps_sigma = 0.025, add_outlier = TRUE)
    outlier_sim <- save_sim_analysis(paste0("sim_", err_struct, "_outlier"), deconv_params = deconv_params, 
        dset_table = NULL, Dataset = list(sim[[paste0("dset_", err_struct)]]), Deconv.Method = c("dtangle", 
            "linearRegression", "CIBERSORT", "EPIC", "Q Prog", "LS Fit"), outf = "./Analysis")
    
    # Markers Near Zero
    cat("Markers Near Zero\n")
    deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, Scale = c("linear"))
    
    dset_params <- list(marker_quantile_low = seq(0.01, 1, 0.01), zero_markers = FALSE)
    dset_table <- expand.grid(dset_params)
    
    Dataset <- lapply(1:nrow(dset_table), function(i) {
        p <- dset_table[i, ]
        gen_sim(marker_quantile_low = p$marker_quantile_low, zero_markers = p$zero_markers)[[paste0("dset_", 
            err_struct)]]
    })
    names(Dataset) <- paste0("Sim", 1:length(Dataset))
    
    save_sim_analysis(paste0("sim_", err_struct, "_marker_value"), deconv_params = deconv_params, 
        dset_table = dset_table, Dataset = Dataset, Deconv.Method = "dtangle", outf = "./Analysis")
    
    # Number of Markers
    cat("Number of Markers\n")
    deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, Scale = c("linear"))
    
    dset_params <- list(pct_markers = seq(0.01, 0.2, 0.01))
    dset_table <- expand.grid(dset_params)
    
    Dataset <- lapply(1:nrow(dset_table), function(i) {
        gen_sim(pct_markers = dset_table[i, ])[[paste0("dset_", err_struct)]]
    })
    names(Dataset) <- paste0("Sim", 1:length(Dataset))
    
    save_sim_analysis(paste0("sim_", err_struct, "_num_markers"), deconv_params = deconv_params, 
        dset_table = dset_table, Dataset = Dataset, Deconv.Method = c("dtangle"), outf = "./Analysis")
}
