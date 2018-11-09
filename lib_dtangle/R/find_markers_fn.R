#' Find marker genes for each cell type.
#' @return List with four elements. ``L'' is respective ranked markers for each cell type and ``V'' is the corresponding values of the ranking method (higher are better) used to determine markers and sort them, ``M'' is the matrix used to create the other two arguments after sorting and subsetting, and ``sM'' is a sorted version of M.
#' @inheritParams dtangle
#' @examples
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' find_markers(Y=Y,pure_samples=pure_samples,
#' data_type='microarray-gene',marker_method='ratio')
#' @export
find_markers <- function(Y, references = NULL, pure_samples = NULL, data_type = NULL, 
    gamma = NULL, marker_method = "ratio") {
    
    cmbd <- combine_Y_refs(Y, references, pure_samples)
    Y <- cmbd$Y
    pure_samples <- cmbd$pure_samples
    
    if (any(lengths(pure_samples) == 1) & marker_method == "p.value") {
        message("Can't use p.value method.")
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
            C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, mean)
        }
        calc_pvals <- function(i) {
            x <- C[, i]
            top <- which(x == max(x))[1]
            second <- order(x, decreasing = TRUE)[2]
            pvs <- rep(NA, length(x))
            for (j in 1:length(x)) {
                pvs[j] <- tryCatch({
                  x1 <- Y[pure_samples[[j]], i]
                  x2 <- Y[pure_samples[[second]], i]
                  n1 <- length(x1)
                  n2 <- length(x2)
                  sd1 <- stats::sd(x1)
                  sd2 <- stats::sd(x2)
                  sp <- sqrt(((n1 - 1) * sd1 + (n2 - 1) * sd2)/(n1 + n2 - 2))
                  t.value <- (mean(x1) - mean(x2))/(sp * sqrt((1/n1) + (1/n2)))
                  tmp <- stats::pt(abs(t.value), df = n1 + n2 - 2)
                  tmp
                })
            }
            pvs[-top] <- 0
            return(pvs)
        }
        C <- sapply(1:ncol(C), calc_pvals)
    } else {
        stop("Marker method not found.")
    }
    pick_top <- function(x) {
        m <- which(x == max(x, na.rm = TRUE))
        if (length(m) != 1) 
            return(c(NA, NaN))
        return(c(m, x[m]))
    }
    M <- apply(C, 2, pick_top)
    M <- data.frame(t(M))
    colnames(M) <- c("top", "value")
    M$rn <- 1:N
    rownames(M) <- colnames(Y)
    M$Cell.Type <- names(pure_samples)[M$top]
    if (marker_method == "p.value") {
        diffmm <- find_markers(Y = Y, pure_samples = pure_samples, data_type = data_type, 
            gamma = gamma, marker_method = "diff")$M
        M$diff <- diffmm$value
        iM <- M[stats::complete.cases(M), ]
        sM <- iM[order(iM$top, -iM$value, -iM$diff), ]
    } else {
        iM <- M[stats::complete.cases(M), ]
        sM <- iM[order(iM$top, -iM$value), ]
    }
    L <- lapply(1:K, function(i) {
        vals <- sM[sM$top == i, "rn"]
        names(vals) <- rownames(sM[sM$top == i, ])
        return(vals)
    })
    V <- lapply(1:K, function(i) {
        vals <- sM[sM$top == i, "value"]
        names(vals) <- rownames(sM[sM$top == i, ])
        return(vals)
    })
    names(L) <- names(pure_samples)
    names(V) <- names(pure_samples)
    return(list(L = L, V = V, M = M, sM = sM))
}

