library("dtangle")

Y <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), byrow = TRUE, nrow = 4)
pure_samples <- list(1, 3)
ml <- find_markers(Y, pure_samples, gamma = 1)$L
bes <- baseline_exprs(Y, pure_samples, ml)

test_that("phats are computed as expected", {
    ph <- est_phats(Y, ml, bes, gamma = 1)
    expect_equal(ph, matrix(c(2/3, 1/3, 2/3, 1/3, 1/3, 2/3, 1/3, 2/3), byrow = TRUE, 
        nrow = 4))
    ph <- est_phats(Y, ml, bes, gamma = 2)
    alpha <- 1/(1 + 2^(-0.5))
    expect_equal(ph, matrix(c(alpha, 1 - alpha, alpha, 1 - alpha, 1 - alpha, alpha, 
        1 - alpha, alpha), byrow = TRUE, nrow = 4))
})
