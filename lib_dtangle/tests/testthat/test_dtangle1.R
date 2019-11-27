.runThisTest <- Sys.getenv("RunAllRRTests") == "yes"

if (.runThisTest) {
    
    ## Basic
    truth <- shen_orr_ex$annotation$mixture
    pure_samples <- lapply(1:3, function(i) {
        which(truth[, i] == 1)
    })
    Y <- shen_orr_ex$data$log
    n_markers <- 20
    
    # Basic: Y and pure_samples
    test_that("basic_deconv", {
        dt_out <- dtangle(Y, pure_samples = pure_samples, n_markers = n_markers, 
            data_type = "microarray-gene", marker_method = "ratio")
        unname(dt_out$markers)
        expect_equal_to_reference(dt_out, "basic.rds")
    })
    
    # Basic: Y and References and pure_samples
    test_that("basic deconv Y, refs, pure_samples", {
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        dt_out <- dtangle(Y, references = references, pure_samples = pure_samples, 
            n_markers = n_markers, data_type = "microarray-gene", marker_method = "ratio")
        comp <- readRDS("basic.rds")
        comp$estimates <- comp$estimates[-unlist(pure_samples), ]
        expect_equal(dt_out, comp)
    })
    
    # Basic: Y and References
    test_that("basic deconv Y and refs", {
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        combined_refs <- t(sapply(pure_samples, function(x) colMeans(references[x, 
            ])))
        comp <- readRDS("basic.rds")
        comp$estimates <- comp$estimates[-unlist(pure_samples), ]
        dt_out <- dtangle(Y, references = combined_refs, n_markers = n_markers, data_type = "microarray-gene", 
            markers = comp$markers)
        expect_equal(dt_out, comp)
    })
    
    # Basic: n_markers
    test_that("basic deconv n_markers", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        comp <- readRDS("basic.rds")
        comp$estimates <- comp$estimates[-unlist(pure_samples), ]
        references <- Y[unlist(pure_samples), ]
        Y <- Y[-unlist(pure_samples), ]
        combined_refs <- t(sapply(pure_samples, function(x) colMeans(references[x, 
            ])))
        dt_out <- dtangle(Y, references = combined_refs, n_markers = NULL, data_type = "microarray-gene", 
            markers = comp$markers)
        
        expect_equal(dt_out, comp)
        
        dt_out <- dtangle(Y, references = combined_refs, n_markers = c(10, 11, 12), 
            data_type = "microarray-gene", markers = comp$markers)
        expect_equal_to_reference(dt_out, "basic_markers.rds")
    })
    
    # Basic: gamma
    test_that("basic deconv gamma", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        dt_out <- dtangle(Y, pure_samples = pure_samples, n_markers = c(10, 11, 12), 
            marker_method = "ratio", gamma = 1)
        expect_equal_to_reference(dt_out, "basic_gamma_1.rds")
        
        dt_out <- dtangle(Y, pure_samples = pure_samples, n_markers = c(10, 11, 12), 
            marker_method = "ratio")
        expect_equal_to_reference(dt_out, "basic_gamma_1.rds")
        
        dt_out <- dtangle(Y, pure_samples = pure_samples, n_markers = c(10, 11, 12), 
            marker_method = "ratio", data_type = "rna-seq")
        expect_equal_to_reference(dt_out, "basic_gamma_seq.rds")
    })
    
    # Basic: markers
    test_that("basic deconv markers", {
        Y <- shen_orr_ex$data$log
        n_markers <- 20
        dt_out <- dtangle(Y, pure_samples = pure_samples, n_markers = c(10, 11, 12), 
            marker_method = "regression")
        expect_equal_to_reference(dt_out, "basic_marker_reg.rds")
    })
    
}
