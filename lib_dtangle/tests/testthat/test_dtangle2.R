.runThisTest <- Sys.getenv("RunAllRRTests") == "yes"

if (.runThisTest) {
    
    
    ## Basic
    truth <- shen_orr_ex$annotation$mixture
    pure_samples <- lapply(1:3, function(i) {
        which(truth[, i] == 1)
    })
    Y <- shen_orr_ex$data$log
    n_markers <- 20
    
    # Basic: Y and References and pure_samples
    test_that("basic dtangle2 Y and refs and ps", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        dt_out <- dtangle2(Y, references = references, pure_samples = pure_samples, 
            n_markers = n_markers, marker_method = "ratio", seed = 4261992, dtangle_init = FALSE)
        expect_equal_to_reference(dt_out, file = "basic_dtangle2.rds")
    })
    
    # Basic: Y and References
    test_that("basic dtangle2 Y and refs", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic_dtangle2.rds")
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        combined_refs <- t(sapply(pure_samples, function(x) colMeans(references[x, 
            ])))
        dt_out <- dtangle2(Y, references = combined_refs, n_markers = n_markers, 
            markers = comp$markers, seed = 4261992, dtangle_init = FALSE)
        expect_equal(dt_out[-5], comp[-5], tolerance = 1e-05)
    })
    
    # Basic: Y and pure_samples
    test_that("basic dtangle2 Y and pure_samples", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic_dtangle2.rds")
        dt_out <- dtangle2(Y, pure_samples = pure_samples, n_markers = n_markers, 
            markers = comp$markers, seed = 4261992, dtangle_init = FALSE)
        expect_equal(dt_out[-5], comp[-5], tolerance = 1e-05)
        # expect_equal(dt_out$estimates[-unlist(pure_samples),],comp$estimates,tolerance=1E-5)
    })
    
    # Basic: n_markers
    test_that("basic dtangle2 n_markers", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic_dtangle2.rds")
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        combined_refs <- t(sapply(pure_samples, function(x) colMeans(references[x, 
            ])))
        dt_out <- dtangle2(Y, references = combined_refs, n_markers = NULL, markers = comp$markers, 
            seed = 4261992, dtangle_init = FALSE)
        expect_equal(dt_out, comp)
        dt_out <- dtangle2(Y, references = combined_refs, n_markers = c(10, 11, 12), 
            markers = comp$markers, seed = 4261992)
        expect_equal_to_reference(dt_out, "basic_dtangle2_markers.rds")
    })
    
    # Basic: markers
    test_that("basic dtangle2 markers", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        dt_out <- dtangle2(Y, pure_samples = pure_samples, n_markers = c(10, 11, 
            12), marker_method = "regression", seed = 4261992)
        expect_equal_to_reference(dt_out, "basic_dtangle2_marker_reg.rds")
    })
    
    # Basic: dtangle_init
    test_that("basic dtangle2 dtangle_init", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        dt_out <- dtangle2(Y, pure_samples = pure_samples, n_markers = n_markers, 
            marker_method = "ratio", seed = 4261992, dtangle_init = FALSE)
        expect_equal_to_reference(dt_out, "basic_dtangle2.rds")
        dt_out <- dtangle2(Y, pure_samples = pure_samples, n_markers = n_markers, 
            marker_method = "ratio", seed = 4261992, dtangle_init = TRUE)
        expect_equal_to_reference(dt_out, "basic_dtangle2_dt_init.rds")
    })
    
    # Basic: variance weights
    test_that("basic dtangle2 variance weights", {
        library("dtangle.data")
        Y <- kuhn$data$log
        pure_samples <- kuhn$annotation$pure_samples
        n_markers <- 20
        dt_out <- dtangle2(Y, pure_samples = pure_samples, n_markers = n_markers, 
            seed = 4261992, dtangle_init = FALSE, weights = "variance")
        expect_equal_to_reference(dt_out[-5], "basic_dtangle2_var_weights.rds")
    })
    
    
    # Basic: arbitrary weights
    test_that("basic dtangle2 arbitrary weights", {
        library("dtangle.data")
        Y <- kuhn$data$log
        pure_samples <- kuhn$annotation$pure_samples
        n_markers <- 20
        wts <- 1:ncol(Y)
        dt_out <- dtangle2(Y, pure_samples = pure_samples, n_markers = n_markers, 
            seed = 4261992, dtangle_init = FALSE, weights = wts)
        expect_equal_to_reference(dt_out[-5], "basic_dtangle2_arb_weights.rds")
    })
    
}
