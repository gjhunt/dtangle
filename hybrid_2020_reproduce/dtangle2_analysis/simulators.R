library("dtangle.data")
library("dtangle")
library("gtools")
library("hitandrun")

# source('./analyze_scripts/save_sim_analysis.R', chdir = TRUE)

tpm_fn <- function(counts, lengths) {
    rate <- counts/lengths
    rate/sum(rate) * 1e+06
}

gen_sim <- function(n_samples = 50, mu = 1, sigma = 1, seed = NULL, dset = parsons, 
    type = NULL, K = NULL, p_rand = NULL) {
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    pure_samples <- dset$annotation$pure_samples
    
    if (is.null(K)) {
        K <- length(pure_samples)
    }
    Z <- t(sapply(dset$annotation$pure_samples, function(x) colMeans(dset$data$count[x, 
        , drop = FALSE])))
    Z <- Z[1:K, ]
    colnames(Z) = make.unique(colnames(Z))
    
    if (is.null(p_rand)) 
        p_rand <- simplex.sample(K, n_samples)$samples
    P <- t(cbind(diag(rep(1, K)), diag(rep(1, K)), t(p_rand)))  # add at edges of simplex
    # P <- rbind(dset$annotation$mixture[unlist(pure_samples), ], P) # add pure
    # samples
    P <- rbind(diag(rep(1, K)), P)
    
    # pure_samples <- lapply(1:K, function(x) which(P[1:length(unlist(pure_samples)),
    # x] == 1))
    pure_samples = lapply(1:K, identity)
    names(pure_samples) <- paste0("ref", 1:K)  #names(dset$annotation$pure_samples)
    
    Y <- P %*% Z
    
    eps_sigma = median(sqrt(sapply(1:3,function(i)median(apply(parsons$data$log[parsons$annotation$pure_samples[[i]],],2,var)))))

    # Poisson
    if (type == "poisson") {
        Y_sim <- apply(Y, c(1, 2), function(x) rpois(1, mu * x))
    } else if (type == "nbinom") {
        # Neg Binom
        Y_sim <- apply(Y, c(1, 2), function(x) {
            x <- x * mu
            size = 1/sigma  #x/(sigma-1+1E-10)+1E-10
            rnbinom(1, mu = x, size = size)
        })
    } else if (type == "gaussian") {
        # Gaussian
        eps_sigma <- sqrt(sigma) * eps_sigma
        print(eps_sigma)
        epsilon_gaussian <- array(rnorm(prod(dim(Y)), 0, eps_sigma), dim(Y))
        Y_sim <- 2^(log2(1 + mu * Y) + epsilon_gaussian) - 1
    }
    Y_sim[1:K, ] <- mu * Z
    epsilon <- log(1 + Y_sim) - log(1 + Y)
    
    rownames(Y_sim) <- paste0("Sample", 1:nrow(Y))
    y <- log2(1 + Y_sim)
    
    SNR = apply(Z, 2, sd)/apply(epsilon, 2, sd)
    
    dset <- list(data = list(log = y, count = Y_sim), annotation = list(mixture = P, 
        pure_samples = pure_samples, data_type = "rna-seq"), name = dset$name)
    
    return(list(Z = Z, Y = Y, P = P, pure_samples = pure_samples, epsilon = epsilon, 
        eps_sigma = eps_sigma, mu = mu, sigma = sigma, dset = dset, snr = SNR))
}
