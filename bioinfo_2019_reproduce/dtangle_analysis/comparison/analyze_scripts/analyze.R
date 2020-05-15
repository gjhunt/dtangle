library("dtangle.data")
source("Utilities/util.R")
source("analyze_dset.R")

analyze <- function(method, q, gamma, dmeths = NULL, verb = TRUE, normalize = TRUE, 
                    datasets = NULL, scl = function(x) 2^x, all_markers = FALSE) {

    if (gamma == 0) 
        gamma <- NULL
    
    sig <- paste(method, q, gamma, "scl(1) = ", scl(1))
    
    if (is.null(dmeths)) {
        dmeths <- c("dtangle", "deconf", "ssFrobenius", "ssKL", "DSA", "Q Prog", 
            "LS Fit", "CIBERSORT", "EPIC")
    }
    
    if (is.null(datasets)) {
        datasets <- get_dtangle_data()
        dset_names <- c("Abbas", "Becht", "Gong", "Kuhn", "Linsley", "Liu", "Newman FL", 
            "Newman PBMC", "Shi", "Shen-Orr", "Parsons")
        datasets <- datasets[dset_names]
    }
    
    updt(paste(sig, "Starting."), init = TRUE)
    output <- lapply(datasets, function(dset) analyze_dset(dset, method, q, gamma, 
        dmeths, verb, normalize, scl, all_markers = all_markers))
    
    p_hat <- lapply(output, "[[", "p_hat")
    p_truth <- lapply(output, "[[", "p_truth")
    n_choose <- lapply(output, "[[", "n_choose")
    markers <- lapply(output, "[[", "markers")
    timing <- lapply(output, "[[", "timing")
    updt(paste(sig, "Complete."))
    
    return(list(p_hat = p_hat, p_truth = p_truth, n = n_choose, time = timing,markers=markers))
}
