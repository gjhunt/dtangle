source("Utilities/CIBERSORT_mod.R")

run_deconv_method <- function(method_name, data, to_deconv, pure_samples, markers, 
    sig_matrix, data_type, gamma, verb, sig) {
    NA_array <- array(NA, c(ncol(to_deconv), ncol(sig_matrix)))
    updt(paste(sig, "Running", method_name))
    out <- tryCatch({
        deconv_method_switch(method_name, data, to_deconv, pure_samples, markers, 
            sig_matrix, data_type, gamma, verb)
    }, error = function(cond) {
        message(red(paste0(">>>Failed ", method_name, ".")))
        cat(red(paste("Caught", cond)))
        return(list(out = NA_array, time = NA))
    })
    
    return(list(estimate = out$phat, time = out$tme))
}


deconv_method_switch <- function(method_name, data, to_deconv, pure_samples, markers, 
    sig_matrix, data_type, gamma, verb) {
    
    pure <- unlist(pure_samples)
    K <- length(pure_samples)
        
    if (!(method_name %in% c("dtangle", "deconf", "ssFrobenius", "ssKL", "DSA"))) {
        to_deconv <- to_deconv[unlist(markers), ]
        sig_matrix <- sig_matrix[unlist(markers), ]

        print(length(rownames(to_deconv)))
        print(length(rownames(sig_matrix)))
    }
    
    
    markers <- switch(method_name, dtangle = markers, MarkerList(markers))

    print(length(unlist(markers)))
    
    
    tme <- system.time(output <- switch(method_name, dtangle = dtangle(Y = data, 
        pure_samples = pure_samples, n_markers = NULL, data_type = data_type, markers = markers, 
        gamma = gamma)$estimates[-pure, ], deconf = t(coef(ged(object = to_deconv, 
        x = markers, method = "deconf", verbose = verb))), ssFrobenius = t(coef(ged(object = to_deconv, 
        x = markers, method = "ssFrobenius", verbose = verb))), ssKL = t(coef(ged(object = to_deconv, 
        x = markers, method = "ssKL", verbose = verb))), DSA = t(coef(ged(object = to_deconv, 
        x = markers, method = "DSA", verbose = verb))), `Q Prog` = t(coef(ged(object = to_deconv, 
        x = sig_matrix, method = "qprog", verbose = verb))), `LS Fit` = t(coef(ged(object = to_deconv, 
        x = sig_matrix, method = "lsfit", verbose = verb))), CIBERSORT = CIBERSORT(sig_matrix, 
        to_deconv)[, 1:K, drop = FALSE], logRegression = t(lm(log(to_deconv) ~ 0 + 
        log(sig_matrix))$coefficients), linearRegression = t(lm(to_deconv ~ 0 + sig_matrix)$coefficients), 
        EPIC = EPIC(bulk = to_deconv, reference = list(refProfiles = sig_matrix, 
            sigGenes = rownames(sig_matrix)))$cellFractions[, 1:K], stop("Method not found.")))
    
    return(list(phat = output, tme = tme))
}
