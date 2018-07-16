#'Estimate the gene type proportions.
#' @param baseline_ests List of vectors (same structure as markers). One list entry for each cell type. Each list element is a vector of estimated offset for each marker of the respective type (output from \code{baseline_exprs}).
#' @inheritParams dtangle
#' @param inv_scale Inverse scale transformation. Default to exponential as dtangle assumes data has been logarithmically transformed.
#' @return Estimated matrix of mixing proportions.
#' @examples 
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' markers = find_markers(Y=Y,pure_samples = pure_samples,
#' data_type='microarray-gene',marker_method='ratio')$L
#' K = length(pure_samples)
#' n_markers = rep(20,K)
#' mrkrs <- lapply(1:K, function(i) {
#'        markers[[i]][1:n_markers[i]]
#' })
#' baseline = dtangle:::baseline_exprs(Y, pure_samples, mrkrs)
#' phats <- dtangle:::est_phats(Y, mrkrs, baseline, gamma=.8)
est_phats <- function(Y, markers, baseline_ests, gamma, summary_fn = mean, inv_scale = function(x) 2^x) {
    
    K <- length(markers)
    markers <- get_marker_list(markers)
    
    contrib_est <- function(i) {
        Y_i <- Y[, markers[[i]], drop = FALSE]
        baseline_adj <- sweep(Y_i, 2, baseline_ests[[i]])
        amt <- inv_scale(apply(baseline_adj/gamma, 1, summary_fn))
        return(amt)
    }
    
    contribs <- sapply(1:K, contrib_est)
    phats <- t(apply(contribs, 1, function(x) {
        x/sum(x)
    }))
    
    colnames(phats) <- names(markers)
    
    return(phats)
}
