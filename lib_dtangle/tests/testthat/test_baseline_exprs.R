library("dtangle")

Y <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), byrow = TRUE, nrow = 4)
pure_samples <- list(1, 3)
ml <- find_markers(Y, pure_samples, gamma = 1)$L

test_that("baseline exprs are computed as expected", {
    bes <- baseline_exprs(Y, pure_samples, ml)
    expect_equal(bes[[1]], 1)
    expect_equal(bes[[2]], 1)
    
    bes <- baseline_exprs(5 * Y, pure_samples, ml)
    expect_equal(bes[[1]], 5)
    expect_equal(bes[[2]], 5)
})
