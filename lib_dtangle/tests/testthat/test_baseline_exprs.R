library("dtangle")

.runThisTest <- Sys.getenv("RunAllRRTests") == "yes"

if (.runThisTest) {
    
    Y <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), byrow = TRUE, nrow = 4)
    pure_samples <- list(1, 3)
    ml <- find_markers(Y = Y, pure_samples = pure_samples, gamma = 1)$L
    
    test_that("baseline exprs are computed as expected", {
        bes <- baseline_exprs(Y = Y, pure_samples = pure_samples, markers = ml)
        expect_equal(bes[[1]][[1]], 1)
        expect_equal(bes[[2]][[1]], 1)
        
        bes <- baseline_exprs(Y = 5 * Y, pure_samples = pure_samples, markers = ml)
        expect_equal(bes[[1]][[1]], 5)
        expect_equal(bes[[2]][[1]], 5)
    })
    
}
