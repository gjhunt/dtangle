find_markers2 <- function(Y, pure_samples, data_type = NULL, gamma = NULL, marker_method = "ratio") {
    if (any(lengths(pure_samples) == 1) & marker_method == "p.value") {
        message("Can't use p.value method.")
        marker_method <- "diff"
    }
    if (is.null(gamma)) 
        gamma <- dtangle:::get_gamma(data_type)
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
                    tmp <- NA
                    if(j!=second){
                        x1 <- Y[pure_samples[[j]], i]
                        x2 <- Y[pure_samples[[second]], i]
                        n1 <- length(x1)
                        n2 <- length(x2)
                        sd1 <- sd(x1)
                        sd2 <- sd(x2)
                        sp <- sqrt(((n1 - 1) * sd1 + (n2 - 1) * sd2)/(n1 + n2 - 2))
                        t.value <- (mean(x1) - mean(x2))/(sp * sqrt((1/n1) + (1/n2)))
                        tmp <- pt(abs(t.value), df = n1 + n2 - 2)
                        if(!is.finite(tmp)){
                        }
                    } else {
                        tmp <- 0
                    }
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
    if (marker_method == "p.value") {
        diffmm <- find_markers2(Y, pure_samples, data_type, gamma, marker_method = "diff")$M
        M$diff <- diffmm$value
        iM <- M[stats::complete.cases(M), ]
        sM <- iM[order(iM$top, -iM$value, -iM$diff), ]
    } else {
        iM <- M[stats::complete.cases(M), ]
        sM <- iM[order(iM$top, -iM$value), ]
    }
    L <- lapply(1:K, function(i) {
        sM[sM$top == i, "rn"]
    })
    V <- lapply(1:K, function(i) {
        sM[sM$top == i, "value"]
    })
    return(list(L = L, V = V, M = M, sM = sM))
}

get_comparison_markers <- function(method, data, q, pure_samples, data_type, sig, 
    all_markers = FALSE) {
    
    updt(paste(sig, "Finding Markers"))
    pure <- unlist(pure_samples)
    K <- length(pure_samples)
    n_choose <- floor((1 - q) * ncol(data)/K)
    
    if (method == "Ratio" || method == "Diff" || method == "p.value") {
        marker_list_ratio <- find_markers2(data, pure_samples, data_type = data_type, 
                                           marker_method = tolower(method))
        ## Choose number of markers
        q_ratio <- sapply(1:K, function(i) {
            quantile(marker_list_ratio$V[[i]], q)
        })
        med_q <- median(q_ratio)
        chs <- function(q_ratio, j) {
            qtn <- q_ratio[j]
            n <- max(which(marker_list_ratio$V[[j]] >= qtn), 1)
            return(n)
        }
        n_choose_ratio <- sapply(1:K, function(a) {
            chs(q_ratio, a)
        })
        marks_ratio <- marker_list_ratio$L
        mrks <- lapply(1:K, function(i) {
            marks_ratio[[i]][1:n_choose_ratio[i]]
        })
        if (all_markers) {
            mrks <- marks_ratio
        }
    } else if (method == "p-value") {
        factors <- unlist(lapply(1:K, function(x) rep(names(pure_samples)[x], length(pure_samples[[x]]))))
        ufactors <- names(pure_samples)
        markers_ged <- extractMarkers(t(data)[, pure], factors, method = "Abbas", 
            log = FALSE)
        for (nm in ufactors) {
            markers_ged[[nm]] <- markers_ged[[nm]][is.finite(markers_ged[[nm]])]
        }
        
        if (all(lengths(markers_ged) == 0)) {
            marker_mtx <- data.frame(markerScoreAbbas(t(data)[, pure], data = factor(factors), 
                log = FALSE))
            markers_ged <- split(marker_mtx$dm2, marker_mtx$top)
            L <- split(rownames(marker_mtx), marker_mtx$top)
            for (i in 1:length(markers_ged)) {
                names(markers_ged[[i]]) <- L[[i]]
            }
            markers_ged <- lapply(markers_ged, sort, decreasing = TRUE)
            names(markers_ged) <- sapply(factors, toString)
        }
        qs <- lapply(ufactors, function(x) {
            quantile(markers_ged[[x]], 1 - q)
        })
        mrks <- lapply(1:K, function(i) {
            which(colnames(data) %in% names(which(markers_ged[[ufactors[i]]] <= qs[i])))
        })
    }
    mrks <- lapply(mrks, function(x) x[1:min(n_choose, length(x))])
    names(mrks) <- names(pure_samples)
    return(mrks)
}
