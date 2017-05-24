#'Estimate the gene type proportions.
#' @param baseline_ests List of vectors (same structure as markers). One list entry for each cell type. Each list element is a vector of estimated offset for each marker of the respective type (output from \code{baseline_exprs}).
#' @inheritParams dtangle
#' @return Estimated matrix of mixing proportions.
#' @export
est_phats <- function(Y, markers, baseline_ests, gamma) {
    
    K <- length(markers)
    
    contrib_est <- function(i) {
        Y_i <- Y[, markers[[i]], drop = FALSE]
        baseline_adj <- sweep(Y_i, 2, baseline_ests[[i]])
        amt <- 2^rowMeans(baseline_adj/gamma)
        return(amt)
    }
    
    contribs <- sapply(1:K, contrib_est)
    phats <- t(apply(contribs, 1, function(x) {
        x/sum(x)
    }))
    
    return(phats)
}
