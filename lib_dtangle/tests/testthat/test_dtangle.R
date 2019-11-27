library("dtangle")

.runThisTest <- Sys.getenv("RunAllRRTests") == "yes"

if (.runThisTest) {
    
    
    test_that("preset gamma values as estimated", {
        expect_equal(dtangle:::gma$ma_probe, 0.4522564, tolerance = 1e-07)
        expect_equal(dtangle:::gma$ma_gene, 0.6999978, tolerance = 1e-07)
        expect_equal(dtangle:::gma$rna_seq, 0.9433902, tolerance = 1e-07)
    })
    
}
