#'Estimate the gene type proportions.
#' @param baseline_ests List of vectors (same structure as markers). One list entry for each cell type. Each list element is a vector of estimated offset for each marker of the respective type (output from \code{baseline_exprs}).
#' @inheritParams dtangle
#' @param gamma Expression sensitivity parameter. A single positive number.
#' @param markers Marker gene indices. List of vectors. List should be same length as \code{pure_samples}, i.e. one element for each cell type. Each element of the top-level list is a vector of indicies (columns of Y) that will be considered markers of that particular type.
#' @param inv_scale Inverse scale transformation. Default to exponential as dtangle assumes data has been logarithmically transformed. 
#' @return Estimated matrix of mixing proportions.
#' @examples
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' markers = find_markers(Y,pure_samples,data_type='microarray-gene',marker_method='ratio')$L
#' K = length(pure_samples)
#' n_choose = rep(20,K)
#' mrkrs <- lapply(1:K, function(i) {
#'        markers[[i]][1:n_choose[i]]
#' })
#' baseline = baseline_exprs(Y, pure_samples, mrkrs)
#' phats <- est_phats(Y, mrkrs, baseline, gamma=.8)
#' @export
est_phats <- function(Y, markers, baseline_ests, gamma, inv_scale = function(x) 2^x) {
    
    K <- length(markers)
    markers <- get_marker_list(markers)
    
    contrib_est <- function(i) {
        Y_i <- Y[, markers[[i]], drop = FALSE]
        baseline_adj <- sweep(Y_i, 2, baseline_ests[[i]])
        amt <- inv_scale(rowMeans(baseline_adj/gamma))
        return(amt)
    }
    
    contribs <- sapply(1:K, contrib_est)
    phats <- t(apply(contribs, 1, function(x) {
        x/sum(x)
    }))
    
    return(phats)
}
