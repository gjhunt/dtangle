get_errors <- function(p_hat, p_truth) {
    
    err <- function(dset_iter) {
        p_hat_i <- p_hat[[dset_iter]]
        if (length(p_hat_i) > 0) {
            p_truth_i <- p_truth[[dset_iter]]
            L <- lapply(p_hat_i, function(x) {
                as.matrix(x - p_truth_i)
            })
            Lm <- melt(L)
            colnames(Lm) <- c("Sample", "Cell.Type", "Error", "Deconv.Method")
            Lm[["Cell.Type"]] <- colnames(p_truth_i)[Lm[["Cell.Type"]]]
            Lm$Estimate <- melt(p_hat_i, measure.vars = 1:ncol(p_hat_i))$value
            Lm$Truth <- melt(p_truth_i, measure.vars = 1:ncol(p_truth_i))$value
            if (is.null(names(p_truth))) 
                names(p_truth) <- "Dset"
            Lm$Dataset <- names(p_truth)[dset_iter]
            return(Lm)
        } else {
            stop("Something went wrong.")
            return(NULL)
        }
    }
    
    E <- lapply(seq_along(p_hat), err)
    
    errors <- do.call(rbind, E)
    errors[["Deconv.Method"]] <- factor(errors[["Deconv.Method"]], levels = c("dtangle", 
        "CIBERSORT", "EPIC", "PERT", "LS Fit", "Q Prog", "linearRegression", "logRegression", 
        "DSA", "ssKL", "ssFrobenius", "deconf"), ordered = TRUE)
    types <- sapply(errors$Cell.Type, toString)
    types <- gsub("_", " ", types)
    types <- sapply(types, capitalize)
    errors$Cell.Type <- factor(types)
    errors$Dataset <- factor(errors$Dataset)
    
    return(errors)
}

all_errors <- function(results, separate_dtangle = FALSE) {
    outs <- results$out
    fn <- function(out) {
        get_errors(out$p_hat, out$p_truth)
    }
    E <- lapply(outs, fn)
    
    if (!is.null(results$dsets)) {
        eDsets <- results$dsets[rep(1:length(unique(levels(E[[1]]$Dataset))), each = table(E[[1]]$Dataset)[1]), 
            , drop = FALSE]
        E <- lapply(E, function(x) cbind(x, eDsets))
    }
    
    E <- lapply(1:length(E), function(i) cbind(E[[i]], results$params[i, ], row.names = NULL))
    allE <- do.call(rbind, E)
    allE$Deconv.Method <- factor(allE$Deconv.Method)
    
    if (separate_dtangle) {
        allE <- allE[-which(allE$Deconv.Method == "dtangle" & allE$Scale == "linear"), 
            ]
        allE[which(allE$Deconv.Method == "dtangle" & allE$Scale == "log"), "Scale"] <- ""
    }
    
    return(data.table(allE))
}

times <- function(res) {
    out <- res$out
    gett <- function(o) {
        ro <- melt(rapply(o$time[lengths(o$time) > 0], function(x) x["elapsed"], 
            how = "list"))
    }
    Ts <- lapply(out, gett)
    neach <- sapply(Ts, nrow)
    param_reps <- unlist(lapply(1:nrow(res$params), function(x) rep(x, times = neach[x])))
    T <- do.call(rbind, Ts)
    T <- cbind(T, res$params[param_reps, ])
    colnames(T)[1:3] <- c("Seconds", "Deconv.Method", "Dataset")
    T$Quantile <- 1 - T$Quantile
    T$Deconv.Method <- factor(T$Deconv.Method)
    T$Deconv.Method <- factor(T$Deconv.Method, levels = c("dtangle", "CIBERSORT", 
        "EPIC", "PERT", "LS Fit", "Q Prog", "linearRegression", "logRegression", 
        "DSA", "ssKL", "ssFrobenius", "deconf"), ordered = TRUE)
    T$Dataset <- factor(T$Dataset)
    T$LogMinutes <- log10(T$Seconds/60)
    return(data.table(T))
}

reduce_newman_res <- function(res, pbmc = "Newman PBMC", fl = "Newman FL") {
    
    if (!is.null(pbmc)) {
        
        newman_pbmc_hl_cts <- c("B", "Monocytes", "NK", "CD4", "CD8", "gamma")
        newman_pbmc_names <- colnames(res$out[[1]]$p_truth[[pbmc]])
        combine_newman_pbmc <- function(mtx) {
            mtx <- sapply(newman_pbmc_hl_cts, function(x) rowSums(mtx[, grepl(x, 
                newman_pbmc_names), drop = FALSE]))
            mtx/rowSums(mtx)
        }
        
        for (i in seq_along(res$out)) {
            res$out[[i]]$p_truth[[pbmc]] <- combine_newman_pbmc(res$out[[i]]$p_truth[[pbmc]])
            res$out[[i]]$p_hat[[pbmc]] <- lapply(res$out[[i]]$p_hat[[pbmc]], combine_newman_pbmc)
        }
    }
    
    if (!is.null(fl)) {
        
        
        newman_fl_hl_cts <- c("B", "CD4", "CD8")
        newman_fl_names <- colnames(res$out[[1]]$p_truth[[fl]])
        combine_newman_fl <- function(mtx) {
            mtx <- sapply(newman_fl_hl_cts, function(x) rowSums(mtx[, grepl(x, newman_fl_names), 
                drop = FALSE]))
            mtx/rowSums(mtx)
        }
        
        for (i in seq_along(res$out)) {
            res$out[[i]]$p_truth[[fl]] <- combine_newman_fl(res$out[[i]]$p_truth[[fl]])
            res$out[[i]]$p_hat[[fl]] <- lapply(res$out[[i]]$p_hat[[fl]], combine_newman_fl)
        }
    }
    return(res)
}

