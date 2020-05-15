library(doParallel)
cl <- makeCluster(max(detectCores() - 2, 1))
registerDoParallel(cl)

init_params <- function(param_list) {
    minimal_params <- list(Quantile = 0.9, Marker.Method = "Ratio", gamma = 0, Scale = "linear", 
        Normalize = TRUE, all_markers = FALSE)
    for (pname in names(minimal_params)) {
        if (is.null(param_list[[pname]])) {
            param_list[[pname]] <- minimal_params[[pname]]
        }
    }
    return(param_list)
}

save_sim_analysis <- function(fname, deconv_params, dset_table = NULL, Dataset = NULL, 
    Deconv.Method = NULL, outf = NULL) {
    
    deconv_params <- init_params(deconv_params)
    params <- expand.grid(deconv_params, stringsAsFactors = FALSE)
    
    if (is.null(outf)) 
        outf <- paste0("./Analysis/", fname, "/", fname, ".rds") else outf <- paste0(outf, "/", fname, "/", fname, ".rds")
    
    cat(toString(Sys.time()), "\n")
    out <- foreach(i = 1:nrow(params)) %dopar% {
        logdir <- paste0("../Analysis/", fname, "/")
        dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
        capture.output({
            p <- params[i, ]
            
            scl <- function(x) 2^x
            if (p$Scale == "log") 
                scl <- base::identity
            
            source("analyze.R")
            a <- analyze(p$Marker.Method, p$Quantile, p$gamma, dmeths = Deconv.Method, 
                normalize = p$Normalize, datasets = Dataset, scl = scl, all_markers = p$all_markers)
            saveRDS(list(a = a, p = p), file = paste0(logdir, paste0(params[i, ], 
                collapse = "_"), "_", i, ".rds"))
            a
        }, file = paste0(logdir, "log", i, ".txt"))
        return(a)
    }
    cat(toString(Sys.time()), "\n")
    
    print(outf)
    saveRDS(list(out = out, params = params, dsets = dset_table), file = outf)
    return(out)
}
