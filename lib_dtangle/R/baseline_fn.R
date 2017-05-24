#' Estimate the offset terms.
#' @return List of vectors. Each vector is estimated estimated basline in pure samples of markers for each group, resp.
#' @inheritParams dtangle
#' @export
baseline_exprs <- function(Y, pure_samples, markers) {
    K <- length(pure_samples)
    
    bl <- list()
    for (i in 1:K) {
        compute_baseline <- function(x) {
            colMeans(Y[pure_samples[[i]], x, drop = FALSE])
        }
        bl[[i]] <- base::sapply(markers[[i]], compute_baseline)
    }
    return(bl)
}
