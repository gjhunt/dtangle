#' Find marker genes for each cell type.
#' @return List with two elements. ``L'' is respective ranked markers for each cell type and ``V'' is the corresponding values of the ranking method (higher are better) used to determine markers and sort them.
#' @inheritParams dtangle
#' @param marker_method The method used to determine which genes are markers. If not supplied defaults to ``ratio''. Options are
#' \itemize{
#' \item{'ratio'}{ selects and ranks markers by the ratio of the mean expression of each gene in each cell type to the mean of that gene in all other cell types.}
#' \item{'regression '}{ selects and ranks markers by estimated regression coefficients in a series of regressions with single covariate that is indicator of each type.}
#' \item{'diff'}{ selects and ranks markers based upon the difference, for each cell type, between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
#' \item{'p.value'}{ selects and ranks markers based upon the p-value of a t-test between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
#' }
#' @examples
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' find_markers(Y,pure_samples,data_type='microarray-gene',marker_method='ratio')
#' @export
find_markers <- function(Y, pure_samples, data_type = NULL, gamma = NULL, marker_method = "ratio") {
    if (any(lengths(pure_samples) == 1) & marker_method == "p.value") {
        message("Can't use p.value method. Using simple differences.")
        marker_method <- "diff"
    }
    if (is.null(gamma)) 
        gamma <- get_gamma(data_type)
    K <- length(pure_samples)
    N <- dim(Y)[2]
    pure <- unlist(pure_samples)
    C <- array(0, c(K, N))
    colnames(C) <- colnames(Y)
    if (marker_method == "ratio") {
        avg_exp_fn <- function(x) {
            colMeans(2^(Y[x, , drop = FALSE]))/gamma
        }
        eta_hats <- t(sapply(pure_samples, avg_exp_fn))
        C <- t(sapply(1:K, function(i) {
            eta_hats[i, ]/apply(eta_hats[-i, , drop = FALSE], 2, sum)
        }))
    } else if (marker_method == "regression") {
        for (i in 1:K) {
            X <- as.numeric(pure %in% pure_samples[[i]])
            Yp <- as.matrix(Y[pure, ])
            m <- stats::lm(Yp ~ 1 + X)
            cfdf <- data.frame(t(stats::coef(m)))
            C[i, ] <- cfdf$X
        }
    } else if (marker_method == "diff") {
        for (i in 1:K) {
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, stats::median)
        }
        less_second <- function(x) {
            x - sort(x, decreasing = TRUE)[2]
        }
        C <- apply(C, 2, less_second)
    } else if (marker_method == "p.value") {
        for (i in 1:K) {
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, stats::median)
        }
        calc_pvals <- function(i) {
            x <- C[, i]
            second <- order(x, decreasing = TRUE)[2]
            pvs <- rep(NA, length(x))
            for (j in 1:length(x)) {
                pvs[j] <- tryCatch({
                  tmp <- 1 - stats::t.test(Y[pure_samples[[j]], i], Y[pure_samples[[second]], 
                    i], alternative = "two.sided")$p.value
                  if (!is.finite(tmp)) 
                    stop("Non-finite")
                  tmp
                }, error = function(e) {
                  0
                })
            }
            return(pvs)
        }
        C <- sapply(1:ncol(C), calc_pvals)
    } else {
        stop("Marker method not found.")
    }
    pick_top <- function(x) {
        m <- which(x == max(x))
        if (length(m) > 1) 
            return(c(NA, NaN))
        return(c(m, x[m]))
    }
    M <- apply(C, 2, pick_top)
    M <- data.frame(t(M))
    colnames(M) <- c("top", "value")
    M$rn <- 1:N
    M <- M[stats::complete.cases(M), ]
    sM <- M[order(M$top, -M$value), ]
    L <- lapply(1:K, function(i) {
        sM[sM$top == i, "rn"]
    })
    V <- lapply(1:K, function(i) {
        sM[sM$top == i, "value"]
    })
    return(list(L = L, V = V))
}
