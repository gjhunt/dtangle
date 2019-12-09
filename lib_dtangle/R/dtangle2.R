#' Deconvolve cell type mixing proportions from gene expression data. 
#' @param Y Expression matrix.
#'
#' (Required) Two-dimensional numeric. Must implement \code{as.matrix}.
#'
#' Each row contains expression measurements for a particular sample. Each columm contains the measurements of the same gene over all individuals. Can either contain just the mixture samples to be deconvolved or both the mixture samples and the reference samples. See \code{pure_samples} and \code{references} for more details.
#' @param references Cell-type reference expression matrix.
#'
#' (Optional) Two-dimensional numeric. Must implement \code{as.matrix}. Must have same number of columns as \code{Y}. Columns must correspond to columns of \code{Y}.
#'
#' Each row contains expression measurements for a reference profile of a particular cell type. Columns contain measurements of reference profiles of a gene. Optionally may merge this matrix with \code{Y} and use \code{pure_samples} to indicate which rows of \code{Y} are pure samples. If \code{pure_samples} is not specified \code{references} must be specified. In this case each row of \code{references} is assumed to be a distinct cell-type. If both \code{pure_samples} and \code{references} are specified then \code{pure_samples} specifies to which cell-type each row of \code{references} corresponds. 
#' @param pure_samples The pure sample indicies.
#'
#' (Optional) List of one-dimensional integer. Must implement \code{as.list}.
#'
#' The i-th element of the top-level list is a vector of indicies (rows of \code{Y} or \code{references}) that are pure samples of type i. If \code{references} is not specified then this argument identifies which rows of \code{Y} correspond to pure reference samples of which cell-types. If \code{references} is specified then this makes same idenficiation but for the \code{references} matrix instead. 
#' @param n_markers Number of marker genes.
#'
#' (Optional) One-dimensional numeric.
#'
#' How many markers genes to use for deconvolution. Can either be a single integer, vector of integers (one for each cell type), or single or vector of percentages (numeric in 0 to 1). If a single integer then all cell types use that number of markers. If a vector then the i-th element determines how many marker genes are used for the i-th cell type. If single percentage (in 0 to 1) then that percentage of markers are used for all types. If vector of percentages then that percentage used for each type, respectively. If not specified then top 10\% of genes are used.
#' @param markers Marker gene indices.
#'
#' (Optional) List of one-dimensional integer.
#'
#' Top-level list should be same length as \code{pure_samples}, i.e. one element for each cell type. Each element of the top-level list is a vector of indicies (columns of \code{Y}) that will be considered markers of that particular type. If not supplied then \code{dtangle} finds markers internally using \code{find_markers}. Alternatively, one can supply the output of \code{find_markers} to the markers argument. 
#' @param marker_method Method used to rank marker genes.
#'
#' (Optional) One-dimensional string.
#' 
#' The method used to rank genes as markers. If not supplied defaults to ``ratio''. Only used if markers are not provided to argument ``markers''. Options are
#' \itemize{
#' \item{'ratio'}{ selects and ranks markers by the ratio of the mean expression of each gene in each cell type to the mean of that gene in all other cell types.}
#' \item{'regression '}{ selects and ranks markers by estimated regression coefficients in a series of regressions with single covariate that is indicator of each type.}
#' \item{'diff'}{ selects and ranks markers based upon the difference, for each cell type, between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
#' \item{'p.value'}{ selects and ranks markers based upon the p-value of a t-test between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
#' }
#' @param weights Weights for the genes.
#'
#' (Optional) String or one-dimensional numeric vector. 
#'
#' Weights for the genes in the optimization. If NULL (default) then does not weight genes differently. If 'variance' then inversely weights with the variance of the references. This only works if there is more than one reference per cell type so that the variance can be estimated. If a numeric then this uses whatever is specified as weights. They must be non-negative. 
#' @param sto Sum-to-one constraint.
#'
#' (Optional) Boolean.
#'
#' Re-normalize the estimates so that the cell-type proportions sum to one.
#' @param inv_scale Inverse scale transformation.
#'
#' (Optional) Function. 
#'
#' Defaults to 2^x. This is equivalent to assuming that the data has been log2-transformed. If another transformation has been applied to the data then this function should be used to specify the inverse of that transformation needed to put gene expressions on the linear scale.
#' @param fit_scale Transformation to used as part of optimization.
#'
#' (Optional) Function.
#'
#' Function to apply to gene expressions as part of optimization. Defaults to log.
#' @param loss_smry Loss summary function minimized to find estimated proportions.
#'
#' (Optional) String.
#'
#' Either 'var' (default) to minimze the (weighted) variance of the residuals or 'L2' to minimize the (weighted) sums of squares of the residuals. 
#' @param dtangle_init Optimization initialization.
#'
#' (Optional) Boolean.
#'
#' Boolean controlling if dtangle2 optimization should be initialized using dtangle1 estimates.
#' @param seed
#'
#' (Optional) Integer.
#'
#' Value at which to seed the random seed before estimating. Optimization initialization might change if this value is not specified.
#' @param verbose
#'
#' (Optional) Boolean.
#'
#' Controls if optimization output is printed or not.
#' @param optim_opts
#'
#' (Optional) List.
#'
#' Optimization options passed to DEoptimR controlling optimization. Options that may be set are
#' \itemize{
#' \item{'constr'} constraint to enforce. Either 'box' for 0-1 box constraints that proportions are between zero and one, 'ineq' for constraints that proportions sum to less than one, 'eq' for equality constraints that proportions sum to one, or 'eq_solve' to solve for one of the parameters in terms of the other and enfoce equality constraints using inequality on remaining parameters. Default and recommended is 'box'.
#' \item{'ninit'} number of randomly initalized points as part of the DEoptimR initial population.
#' \item{'tritter'} how often to print results if 'verbose=TRUE'.
#' \item{'maxiter'} maximum number of optimization iterations to use before exiting.
#' \item{'convtol'} tolerance for convergence tolerance stopping criterion. 
#' \item{'constrtol'} tolerance for constraint enforcement.
#' }
#' @return List.
#' \itemize{
#' \item{'estimates'}{ a matrix estimated mixing proportions. One row for each sample, one column for each cell type.}
#' \item{'markers'}{ list of vectors of marker used for each cell type. Each element of list is vector of columns of \code{Y} used as a marker for the i-th cell type.}
#' \item{'n_markers'}{ vector of number of markers used for each cell type.}
#' \item{'weights'}{ the weights used as part of the optimization.}
#' \item{'diag'}{ diagnostic values for the estimated proportions.
#' \code{resids_hat},\code{loss_hat}, and \code{p_hat} are the residuals, loss, and estimates for the proportions returned by dtangle2. Similarly, \code{resids_opt},\code{loss_opt} and \code{p_opt} are these values for the optimized value not re-scaled to enforce the STO constraint. 
#' }
#' }
#' @examples
#' \donttest{
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' n_markers = 20
#'
#' dtangle2(Y, pure_samples = pure_samples,
#' n_markers=n_markers)
#' }
#' 
#' @seealso \code{\link{find_markers}}
#' @export
dtangle2 <- function(Y, references = NULL, pure_samples = NULL, n_markers = NULL, 
    markers = NULL, marker_method = "ratio", weights = NULL, sto = TRUE, inv_scale = function(x) 2^x, 
    fit_scale = log, loss_smry = "var", dtangle_init = TRUE, seed = NULL, verbose = FALSE, 
    optim_opts = NULL) {
    
    stopifnot(all(n_markers > 0))
    stopifnot(!is.null(c(references, pure_samples)))
    
    cmbd <- combine_Y_refs(Y, references, pure_samples)
    Y <- cmbd$Y
    pure_samples <- cmbd$pure_samples
    
    # Removed arguments
    data_type <- NULL
    gamma <- 1
    
    prc <- process_markers(Y, pure_samples, n_markers, data_type, gamma, markers, 
        marker_method)
    n_markers <- prc$n_markers
    mrkrs <- prc$mrkrs
    gamma <- prc$gamma
    
    est <- dtangle2_est_phats(Y = Y, pure_samples = pure_samples, markers = mrkrs, 
        weights = weights, optim_opts = optim_opts, sto = sto, inv_scale = inv_scale, 
        dtangle_init = dtangle_init, fit_scale = fit_scale, loss_smry = loss_smry, 
        seed = seed, verbose = verbose)
    phats <- est$phats
    diag <- list(resids_hat = est$resids_hat, resids_opt = est$resids_opt, loss_hat = est$loss_hat, 
        loss_opt = est$loss_opt, p_opt = est$popt, p_hat = est$phats, opt = est$opt)
    
    return(list(estimates = phats, markers = mrkrs, n_markers = n_markers, weights = est$W, 
        diag = diag))
}

dtangle2_est_phats <- function(Y, pure_samples = NULL, markers = NULL, weights = NULL, 
    optim_opts = NULL, sto = TRUE, inv_scale = function(x) 2^x, dtangle_init = FALSE, 
    fit_scale = NULL, loss_smry = "var", seed = NULL, verbose = FALSE) {
    
    set.seed(seed)
    
    optim_opts <- set_optim_opts(optim_opts)
    
    pure <- unlist(pure_samples)
    mrks <- unlist(markers)
    
    X <- inv_scale(Y[-pure, mrks, drop = FALSE])
    
    Z <- t(sapply(pure_samples, function(x) apply(Y[x, mrks, drop = FALSE], 2, mean)))
    Z <- inv_scale(Z)
    
    if (is.null(weights)) {
        W <- rep(1, ncol(X))
    } else if (all(weights == "variance")) {
        Yss <- fit_scale(inv_scale(Y))
        SS <- t(sapply(pure_samples, function(x) apply(Yss[x, mrks, drop = FALSE], 
            2, function(x) sum((x - mean(x))^2))))
        SSS <- apply(SS, 2, sum)
        denom <- sum(lengths(pure_samples)) - length(pure_samples)
        VARS <- SSS/denom
        VARS[VARS < stats::median(VARS)] <- stats::median(VARS)
        W <- 1/VARS
    } else {
        W <- weights[mrks]
    }
    W <- W/sum(W, na.rm = TRUE)
    
    Kall <- nrow(Z)
    P0 <- array(1/Kall, c(nrow(X), Kall))
    if (dtangle_init) 
        P0 <- dtangle(Y = Y, pure_samples = pure_samples, markers = markers, gamma = 1)$estimates
    
    i <- 1
    dt_opt <- function(x) {
        if (verbose) 
            cat("====> ", i, "\n")
        opt <- optim_solve(x, Z, W = W, optim_opts = optim_opts, p0 = P0[i, ], fit_scale = fit_scale, 
            loss_smry = loss_smry, verbose = verbose)
        popt <- opt$phat
        phat <- popt
        lossf <- function(p) L(p, x, Z, W, sto = FALSE, fit_scale = fit_scale, loss_smry = loss_smry)
        residf <- function(p) resid(p, x, Z, W, sto = FALSE, fit_scale = fit_scale)
        if (optim_opts$constr == "eq_solve") {
            phat <- itob(phat)
            lossf <- function(p) L_ineq(p, x, Z, W, sto = FALSE, fit_scale = fit_scale, 
                loss_smry = loss_smry)
            residf <- function(p) resid_ineq(p, x, Z, W, sto = FALSE, fit_scale = fit_scale)
        }
        phat <- phat/sum(phat)
        resid_opt <- residf(popt)
        loss_opt <- lossf(popt)
        resid_hat <- resid(phat, x, Z = Z, W = W, fit_scale = fit_scale)
        loss_hat <- L(phat, x, Z, W, sto = FALSE, fit_scale = fit_scale, loss_smry = loss_smry)
        i <<- i + 1
        return(list(opt = opt, popt = popt, phat = phat, resid_opt = resid_opt, resid_hat = resid_hat, 
            loss_opt = loss_opt, loss_hat = loss_hat))
    }
    
    opts <- apply(X, 1, dt_opt)
    
    if (sto) 
        phats_all <- t(sapply(opts, "[[", "phat")) else phats_all <- t(sapply(opts, "[[", "popt"))
    
    colnames(phats_all) <- rownames(Z)
    
    resids_hat <- sapply(opts, "[[", "resid_hat")
    resids_opt <- sapply(opts, "[[", "resid_opt")
    loss_hat <- sapply(opts, "[[", "loss_hat")
    loss_opt <- sapply(opts, "[[", "loss_opt")
    popt <- t(sapply(opts, "[[", "popt"))
    opt <- t(sapply(opts, "[[", "opt"))
    
    return(list(phats = phats_all, W = W, resids_hat = resids_hat, resids_opt = resids_opt, 
        loss_hat = loss_hat, loss_opt = loss_opt, popt = popt, opt = opt))
}

optim_solve <- function(x, Z, W = NULL, optim_opts = NULL, maxiter = 5000, tol = 1e-10, 
    p0 = NULL, fit_scale = NULL, loss_smry = "var", verbose = TRUE) {
    
    Kall <- nrow(Z)
    if (is.null(p0)) 
        p0 <- rep(1/Kall, Kall)
    if (is.null(W)) 
        W <- rep(1, length(x))
    
    loss <- function(p) L(p, x, Z, W, sto = FALSE, fit_scale = fit_scale, loss_smry = loss_smry)
    if (optim_opts$constr == "eq_solve") {
        loss <- function(p) L_ineq(p, x, Z, W, sto = FALSE, fit_scale = fit_scale, 
            loss_smry = loss_smry)
        p0 <- p0[-length(p0)]
    }
    loss_deriv <- NULL
    
    optimize <- switch(optim_opts$pkg, DEoptimR = optim_solve_deoptimr(x = x, Z = Z, 
        loss = loss, loss_deriv = loss_deriv, W = W, optim_opts = optim_opts, maxiter = maxiter, 
        tol = tol, p0 = p0, verbose = verbose))
    
    return(list(phat = optimize$solution, opt = optimize))
}

optim_solve_deoptimr <- function(x, Z, loss, loss_deriv, W = NULL, optim_opts = NULL, 
    maxiter = 5000, tol = 1e-12, p0 = NULL, verbose = TRUE) {
    
    
    MEQ <- 0
    constrfn <- NULL  ## constr = 'box'
    if (optim_opts$constr == "eq") 
        MEQ <- 1
    if (optim_opts$constr %in% c("eq", "eq_solve", "ineq")) 
        constrfn <- constr
    
    out <- tryCatch({
        Kall <- length(p0)
        fs <- loss(p0)
        ADD <- cbind(diag(1, Kall), rep(1/Kall, Kall))
        RAND <- simplex_sample(N = optim_opts$ninit, Kall)
        ADD <- cbind(ADD, t(RAND))
        list(opt_out = DEoptimR::JDEoptim(lower = rep(0, Kall), upper = rep(1, Kall), 
            fn = loss, trace = verbose, triter = optim_opts$triter, maxiter = optim_opts$maxiter, 
            add_to_init_pop = ADD, tol = optim_opts$convtol, fnscale = fs, constr = constrfn, 
            meq = MEQ, eps = optim_opts$constrtol, details = TRUE), err = NULL)
    }, error = function(e) {
        list(opt_out = NULL, err = e)
    })
    out$solution <- out$opt_out$par
    return(out)
}

set_optim_opts <- function(optim_opts) {
    default <- list(pkg = "DEoptimR", maxiter = 5000, ninit = 20, convtol = 1e-12, 
        constr = "box", constrtol = 1e-06)
    for (nm in names(optim_opts)) {
        default[[nm]] <- optim_opts[[nm]]
    }
    if (is.null(optim_opts$triter)) 
        default$triter <- default$maxiter/10
    return(default)
}
