library("dtangle.data")
library("dtangle")
library("gtools")
library("hitandrun")

source("./analyze_scripts/save_sim_analysis.R", chdir = TRUE)

tpm_fn <- function(counts, lengths) {
    rate <- counts/lengths
    rate/sum(rate) * 1e+06
}

gen_sim <- function(n_samples = 50, eps_sigma = 0.025, add_outlier = FALSE, seed = NULL, 
    outlier_multiplier = 2, num_outlier = 10, dset = parsons) {
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    pure_samples <- dset$annotation$pure_samples
    K <- length(pure_samples)
    Z <- t(sapply(dset$annotation$pure_samples, function(x) colMeans(dset$data$count[x, 
        , drop = FALSE])))
    
    p_rand <- simplex.sample(K, n_samples)$samples
    P <- t(cbind(diag(rep(1, K)), diag(rep(1, K)), t(p_rand)))
    P <- rbind(dset$annotation$mixture[unlist(pure_samples), ], P)
    
    pure_samples <- lapply(1:K, function(x) which(P[1:length(unlist(pure_samples)), 
        x] == 1))
    names(pure_samples) <- names(dset$annotation$pure_samples)
    
    Y <- P %*% Z
    
    Y_poisson <- apply(Y, c(1, 2), function(x) rpois(1, x))
    epsilon_poisson <- Y_poisson - Y
    
    
    eps_sigma <- sd(log2(1 + Z)) * eps_sigma
    epsilon_gaussian <- array(rnorm(prod(dim(Y)), 0, eps_sigma), dim(Y))
    Y_gaussian <- 2^(log2(Y) + epsilon_gaussian)
    
    rownames(Y_poisson) <- paste0("Sample", 1:nrow(Y))
    rownames(Y_gaussian) <- paste0("Sample", 1:nrow(Y))
    
    y_poisson <- log2(1 + Y_poisson)
    y_gaussian <- log2(1 + Y_gaussian)
    
    if (add_outlier) {
        outlier_value <- exp(log(max(y_poisson)) * outlier_multiplier)
        y_poisson[pure_samples[[1]], sample(ncol(y_poisson), num_outlier)] <- outlier_value
        y_poisson[pure_samples[[2]], sample(ncol(y_poisson), num_outlier)] <- outlier_value
        y_poisson[pure_samples[[3]], sample(ncol(y_poisson), num_outlier)] <- outlier_value
        outlier_value <- exp(log(max(y_gaussian)) * outlier_multiplier)
        y_gaussian[pure_samples[[1]], sample(ncol(y_gaussian), num_outlier)] <- outlier_value
        y_gaussian[pure_samples[[2]], sample(ncol(y_gaussian), num_outlier)] <- outlier_value
        y_gaussian[pure_samples[[3]], sample(ncol(y_gaussian), num_outlier)] <- outlier_value
    }
    
    Y_gaussian <- 2^y_gaussian - 1
    Y_poisson <- 2^y_poisson - 1
    
    tpm_gaussian <- t(apply(t(Y_gaussian), 2, function(x) tpm_fn(x, dset$annotation$gene.length)))
    tpm_poisson <- t(apply(t(Y_poisson), 2, function(x) tpm_fn(x, dset$annotation$gene.length)))
    
    dset_gaussian <- list(data = list(log = y_gaussian, count = Y_gaussian, tpm = tpm_gaussian), 
        annotation = list(mixture = P, pure_samples = pure_samples, data_type = "rna-seq"), 
        name = dset$name)
    dset_poisson <- list(data = list(log = y_poisson, count = Y_poisson, tpm = tpm_poisson), 
        annotation = list(mixture = P, pure_samples = pure_samples, data_type = "rna-seq"), 
        name = dset$name)
    
    
    return(list(Z = Z, Y = Y, P = P, pure_samples = pure_samples, epsilon_gaussian = epsilon_gaussian, 
        epsilon_poisson = epsilon_poisson, eps_sigma = eps_sigma, dset_gaussian = dset_gaussian, 
        dset_poisson = dset_poisson))
}
