library("dtangle.data")
library("dtangle")
library("gtools")
library("hitandrun")

source("./analyze_scripts/save_sim_analysis.R", chdir = TRUE)

tpm_fn <- function(counts, lengths) {
    rate <- counts/lengths
    rate/sum(rate) * 1e+06
}

gen_sim <- function(K = 3, pct_markers = 0.15, marker_quantile_top = 0.75, marker_quantile_low = 0.05, 
    zero_markers = TRUE, shared_markers = FALSE, n_samples = 50, eps_sigma = 0.025, 
    add_outlier = FALSE, seed = NULL) {
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    profile <- apply(parsons$data$count, 2, median)
    names(profile) <- make.unique(names(profile))
    N <- length(profile)
    N_markers <- floor(N * pct_markers/K)
    candidates <- which(profile > quantile(profile, marker_quantile_top))
    if (length(candidates) < N_markers) 
        stop("Not enough candidates.")
    all_markers <- sample(candidates, N_markers * K)
    
    roll_fn <- function(K) {
        if (shared_markers) {
            return(rbinom(K, 1, sel_prob) == 1)
        } else {
            rl <- runif(K)
            return(rl == min(rl))
        }
    }
    
    marker_list <- list()
    for (mark in all_markers) {
        sel_prob <- 1/K
        roll <- roll_fn(K)
        for (k in 1:K) {
            if (roll[k]) {
                marker_list <- tryCatch({
                  marker_list[[k]] <- c(marker_list[[k]], mark)
                  marker_list
                }, error = function(err) {
                  marker_list[[k]] <- mark
                  marker_list
                })
            }
        }
    }
    markers <- marker_list
    all_markers <- unique(unlist(markers))
    
    profiles <- rep(list(profile), K)
    profiles <- do.call(rbind, profiles)
    
    marker_value <- quantile(profile, marker_quantile_low)
    if (zero_markers) 
        marker_value <- 0
    
    for (mark in all_markers) {
        marks <- which(sapply(markers, function(x) mark %in% x) == TRUE)
        profiles[-marks, mark] <- marker_value
    }
    
    Z <- profiles
    
    p_rand <- simplex.sample(K, n_samples)$samples
    P <- t(cbind(diag(rep(1, K)), diag(rep(1, K)), t(p_rand)))
    
    Y <- P %*% Z
    
    Y_poisson <- apply(Y, c(1, 2), function(x) rpois(1, x))
    epsilon_poisson <- Y_poisson - Y
    
    eps_sigma <- sd(log2(1 + Z)) * eps_sigma
    epsilon_gaussian <- array(rnorm(prod(dim(Y)), 0, eps_sigma), dim(Y))
    Y_gaussian <- 2^(log2(Y) + epsilon_gaussian)
    
    pure_samples <- lapply(1:K, function(i) which(P[, i] == 1))
    names(pure_samples) <- LETTERS[1:K]
    colnames(P) <- names(pure_samples)
    rownames(Y_poisson) <- paste0("Sample", 1:nrow(Y))
    rownames(Y_gaussian) <- paste0("Sample", 1:nrow(Y))
    
    y_poisson <- log2(1 + Y_poisson)
    y_gaussian <- log2(1 + Y_gaussian)
    
    if (add_outlier) {
        multiplier <- 2
        outlier_value <- max(y_poisson) * multiplier
        num_outlier <- 10
        y_poisson[pure_samples[[1]], markers[[1]][1:num_outlier]] <- outlier_value
        y_poisson[pure_samples[[2]], markers[[2]][1:num_outlier]] <- outlier_value
        y_poisson[pure_samples[[3]], markers[[3]][1:num_outlier]] <- outlier_value
        outlier_value <- max(y_gaussian) * multiplier
        y_gaussian[pure_samples[[1]], markers[[1]][1:num_outlier]] <- outlier_value
        y_gaussian[pure_samples[[2]], markers[[2]][1:num_outlier]] <- outlier_value
        y_gaussian[pure_samples[[3]], markers[[3]][1:num_outlier]] <- outlier_value
    }
    
    Y_gaussian <- 2^y_gaussian - 1
    Y_poisson <- 2^y_poisson - 1
    
    tpm_gaussian <- t(apply(t(Y_gaussian), 2, function(x) tpm_fn(x, parsons$annotation$gene.length)))
    tpm_poisson <- t(apply(t(Y_poisson), 2, function(x) tpm_fn(x, parsons$annotation$gene.length)))
    
    dset_gaussian <- list(data = list(log = y_gaussian, count = Y_gaussian, tpm = tpm_gaussian), 
        annotation = list(mixture = P, pure_samples = pure_samples, data_type = "rna-seq"), 
        name = "Parsons")
    dset_poisson <- list(data = list(log = y_poisson, count = Y_poisson, tpm = tpm_poisson), 
        annotation = list(mixture = P, pure_samples = pure_samples, data_type = "rna-seq"), 
        name = "Parsons")
    
    
    return(list(Z = Z, Y = Y, P = P, markers = markers, pure_samples = pure_samples, 
        epsilon_gaussian = epsilon_gaussian, epsilon_poisson = epsilon_poisson, eps_sigma = eps_sigma, 
        dset_gaussian = dset_gaussian, dset_poisson = dset_poisson))
}
