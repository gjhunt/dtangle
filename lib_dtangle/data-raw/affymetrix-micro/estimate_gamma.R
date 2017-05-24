source("latin_data_read.R")
eps <- 1e-07

x <- log2(all_titration[, 1] + eps)
y <- sd[, 1]

zc <- 1:3

gfit <- sapply(1:dim(sdg)[2], function(i) {
    m <- lm(sdg[-zc, i] ~ 1 + x[-zc])
    return(coef(m)[2])
})

pfit <- sapply(1:dim(sd)[2], function(i) {
    m <- lm(sd[-zc, i] ~ 1 + x[-zc])
    return(coef(m)[2])
})

g_gamma <- median(gfit)
p_gamma <- median(pfit)
