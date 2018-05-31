#' Estimate the offset terms.
#' @return List of vectors. Each vector is estimated estimated basline in pure samples of markers for each group, resp.
#' @inheritParams dtangle
#' @param markers Marker gene indices. List of vectors. List should be same length as \code{pure_samples}, i.e. one element for each cell type. Each element of the top-level list is a vector of indicies (columns of Y) that will be considered markers of that particular type.
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
#' baseline_exprs(Y, pure_samples, mrkrs)
#' @export
baseline_exprs <- function(Y, pure_samples, markers) {
    K <- length(pure_samples)
    markers <- get_marker_list(markers)
    
    bl <- list()
    for (i in 1:K) {
        compute_baseline <- function(x) {
            colMeans(Y[pure_samples[[i]], x, drop = FALSE])
        }
        bl[[i]] <- base::sapply(markers[[i]], compute_baseline)
    }
    return(bl)
}
