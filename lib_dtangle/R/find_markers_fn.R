#' Find marker genes for each cell type.
#' @return List with two elements. ``L'' is respective ranked markers for each cell type and ``V'' is the corresponding values of the ranking method (higher are better) used to determine markers and sort them.
#' @inheritParams dtangle
#' @export
find_markers <- function(Y, pure_samples, data_type = NULL, gamma = NULL, marker_method = "ratio") {
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
    }
    
    if (marker_method == "regression") {
        for (i in 1:K) {
            X <- as.numeric(pure %in% pure_samples[[i]])
            Yp <- as.matrix(Y[pure, ])
            m <- stats::lm(Yp ~ 1 + X)
            cfdf <- data.frame(t(stats::coef(m)))
            C[i, ] <- cfdf$X
        }
    }
    
    if (marker_method == "diff") {
        for (i in 1:K) {
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, stats::median)
        }
        less_second <- function(x) {
            x - sort(x, decreasing = TRUE)[2]
        }
        C <- apply(C, 2, less_second)
    }
    
    pick_top <- function(x) {
        m <- which.max(x)
        return(c(m, x[m]))
    }
    
    M <- apply(C, 2, pick_top)
    M <- data.frame(t(M))
    colnames(M) <- c("top", "value")
    M$rn <- 1:N
    
    sM <- M[order(M$top, -M$value), ]
    L <- lapply(1:K, function(i) {
        sM[sM$top == i, "rn"]
    })
    
    V <- lapply(1:K, function(i) {
        sM[sM$top == i, "value"]
    })
    
    return(list(L = L, V = V))
}
