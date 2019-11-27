simplex_sample <- function(N = 1, d = 2) {
    e <- array(stats::rexp(d * N, 1), c(N, d))
    s <- e/rowSums(e)
    return(s)
}
