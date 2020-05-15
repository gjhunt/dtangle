library("CellMix")
library("nnls")
library("dtangle")
library("EPIC")
library("limma")
library("reshape2")
library("dtangle.data")

source("Utilities/util.R")
source("get_markers.R")
source("run_deconv_methods.R")

analyze_dset <- function(dset, method, q, gamma, dmeths, verb, normalize, scl, all_markers = FALSE) {
    sig <- paste(method, q, gamma, "scl(1) = ", scl(1), dset$name)
    
    p_hat <- list()
    timing <- list()
    mark_list <- list()
    
    data_type <- dset$annotation$data_type

    
    
    if (data_type == "rna-seq") {
        data <- log2(1 + as.matrix(dset$data$tpm))
    } else if (data_type == "microarray-gene") {
        data <- as.matrix(dset$data$log)
        data = data[, colnames(data) != "sample"]
        print(tail(colnames(data)))
        print(normalize)
        if(normalize)
            data <- t(normalizeBetweenArrays(t(data)))
    }
    
    

    colnames(data) <- make.unique(colnames(data))
    
    # Pure samples
    annotation <- dset$annotation
    pure_samples <- annotation$pure_samples
    pure <- unlist(pure_samples)
    K <- length(pure_samples)

    ## REMOVED SO THAT WE CAN FIT USING LOG SCALE
    ## Get Signature Matrix for Cell Types
    #sig_matrix <- sapply(pure_samples, function(x) colMeans(scl(data)[x, , drop = FALSE]))
    #to_deconv <- t(scl(data))[, -pure]

    sig_matrix <- sapply(pure_samples, function(x) colMeans(2^(data)[x, , drop = FALSE]))
    to_deconv <- t(2^(data))[, -pure]
    
    if (!all_markers) {
        minvar <- quantile(apply(sig_matrix, 1, var), 0.75)
        keep_vars <- which(apply(to_deconv, 1, var) > 0 & apply(sig_matrix, 1, var) > 
            minvar)
        to_deconv <- to_deconv[keep_vars, ]
        sig_matrix <- sig_matrix[keep_vars, ]
        data <- data[, keep_vars]
    }
    
    ## Get Marker Genes
    mrks <- get_comparison_markers(method, data, q, pure_samples, data_type, sig, 
                                   all_markers)

    cat("MARKERS:\n")
    print(mrks)

    ## Change scale if needed
    sig_matrix = scl(log2(sig_matrix))
    to_deconv = scl(log2(to_deconv))

    cat("MARKERS NAMES:\n")
    print(rownames(sig_matrix)[unlist(mrks)])
    print(rownames(to_deconv)[unlist(mrks)])
        
    ## Run Methods
    p_truth <- annotation$mixture[-pure, ]
    n_choose <- lengths(mrks)

    cat("MARKERS NUMBER:\n")
    cat(n_choose)
    
    for (mth in dmeths) {
        dcnv_out <- run_deconv_method(mth, data, to_deconv, pure_samples, mrks, sig_matrix, 
            data_type, gamma, verb, sig)
        p_hat[[mth]] <- dcnv_out$estimate
        timing[[mth]] <- dcnv_out$time
    }
    
    return(list(n_choose = n_choose, p_hat = p_hat, p_truth = p_truth, timing = timing, markers=mrks))
}

