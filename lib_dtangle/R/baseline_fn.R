#' Estimate the offset terms.
#' @return List of vectors. Each vector is estimated estimated basline in pure samples of markers for each group, resp.
#' @inheritParams dtangle
#' @examples
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' markers = find_markers(Y=Y,
#' pure_samples = pure_samples,data_type='microarray-gene',marker_method='ratio')$L
#' K = length(pure_samples)
#' n_markers = rep(20,K)
#' mrkrs <- lapply(1:K, function(i) {
#'        markers[[i]][1:n_markers[i]]
#' })
#' dtangle:::baseline_exprs(Y, pure_samples, mrkrs)
baseline_exprs <- function(Y, pure_samples, markers, summary_fn = mean) {
    K <- length(pure_samples)
    markers <- get_marker_list(markers)
    
    bl <- list()
    for (i in 1:K) {
        compute_baseline <- function(x) {
            apply(Y[pure_samples[[i]], x, drop = FALSE], 2, summary_fn)
        }
        bl[[i]] <- base::sapply(markers[[i]], compute_baseline)
    }
    return(bl)
}
