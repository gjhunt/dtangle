set.seed(3301991)

source("mix_simulators.R")
for (err_struct in c("gaussian", "poisson")) {
    cat(err_struct, "\n")
    for (dset in c("linsley", "parsons")) {
        cat(dset, "\n")
        # Low Error
        
        cat("Low Error\n")
        deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, 
            Scale = c("linear", "log"))
        sim <- gen_sim(eps_sigma = 0.025, dset = get(dset))
        lowe_sim <- save_sim_analysis(paste0("sim_mix_", err_struct, "_lowe_", dset), 
            deconv_params = deconv_params, dset_table = NULL, Dataset = list(sim[[paste0("dset_", 
                err_struct)]]), Deconv.Method = c("dtangle", "linearRegression", 
                "CIBERSORT", "EPIC", "LS Fit", "Q Prog"), outf = "./Analysis")
        
        # High Error
        cat("High Error\n")
        deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, 
            Scale = c("linear", "log"))
        sim <- gen_sim(eps_sigma = 0.75, dset = get(dset))
        highe_sim <- save_sim_analysis(paste0("sim_mix_", err_struct, "_highe_", 
            dset), deconv_params = deconv_params, dset_table = NULL, Dataset = list(sim[[paste0("dset_", 
            err_struct)]]), Deconv.Method = c("dtangle", "linearRegression", "CIBERSORT", 
            "EPIC", "Q Prog", "LS Fit"), outf = "./Analysis")
        
        # Outliers
        cat("Outliers\n")
        deconv_params <- list(Quantile = 0.9, Marker.Method = "p.value", gamma = 1, 
            Scale = c("linear", "log"))
        sim <- gen_sim(eps_sigma = 0.025, add_outlier = TRUE, outlier_multiplier = 1.25, 
            num_outlier = 5, , dset = get(dset))
        outlier_sim <- save_sim_analysis(paste0("sim_mix_", err_struct, "_outlier_", 
            dset), deconv_params = deconv_params, dset_table = NULL, Dataset = list(sim[[paste0("dset_", 
            err_struct)]]), Deconv.Method = c("dtangle", "linearRegression", "CIBERSORT", 
            "EPIC", "Q Prog", "LS Fit"), outf = "./Analysis")
        
    }
}
