constr <- function(p) {
    sum(p) - 1
}

check_constr <- function(p) {
    if (abs(sum(p) - 1) > 1e-10) 
        return(FALSE)
    if (any(p < 0)) 
        return(FALSE)
    return(TRUE)
}

constr_deriv <- function(p) {
    rep(1, length(p))
}

p_to_mtx <- function(p, Kall, sto = FALSE) {
    mat <- matrix(p, ncol = Kall, byrow = TRUE)
    if (sto) 
        return(mat/rowSums(mat))
    return(mat)
}

p_to_vec <- function(p, sto = FALSE) {
    if (sto) 
        p <- p/rowSums(p)
    as.vector(t(p))
}

weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w)/sum(w)
    (sum.w/(sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na.rm)
}

L <- function(p, y, Z, W, sto = FALSE, fit_scale = NULL, loss_smry = NULL) {
    if (is.null(fit_scale)) 
        fit_scale <- log
    Kall <- nrow(Z)
    p <- p_to_mtx(p, Kall = Kall, sto = sto)
    if (any(y == 0)) 
        stop("Some values of Y are zero.")
    ly <- fit_scale(y)
    lpz <- fit_scale(p %*% Z)
    diff <- as.vector(ly - lpz)
    if (is.null(loss_smry) | loss_smry == "L2") {
        nrm <- sum((W * diff)^2)
    }
    if (loss_smry == "var") 
        nrm <- weighted.var(diff, W)
    if (!is.finite(nrm)) 
        return(Inf)
    return(nrm)
}

resid <- function(p, y, Z, W, sto = FALSE, fit_scale = NULL) {
    if (is.null(fit_scale)) 
        fit_scale <- log
    Kall <- nrow(Z)
    p <- p_to_mtx(p, Kall = Kall, sto = sto)
    ly <- fit_scale(y)
    lpz <- fit_scale(p %*% Z)
    diff <- as.vector(ly - lpz)
    return(diff)
}

## Ineq constraints
itob <- function(p) {
    return(c(p, max(1 - sum(p), 0)))
}

L_ineq <- function(pm1, ...) {
    L(itob(pm1), ...)
}

resid_ineq <- function(pm1, ...) {
    resid(itob(pm1), ...)
}

