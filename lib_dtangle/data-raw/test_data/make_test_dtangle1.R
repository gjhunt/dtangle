library("dtangle", lib.loc = "~/dev/dtangle/0.1.0")

# Basic
truth <- shen_orr_ex$annotation$mixture
pure_samples <- lapply(1:3, function(i) {
    which(truth[, i] == 1)
})
Y <- shen_orr_ex$data$log
n_choose <- 20

dt_out <- dtangle(Y, pure_samples = pure_samples, n_choose = n_choose, data_type = "microarray-gene", 
    marker_method = "ratio")

saveRDS(dt_out, "basic.rds")

# basic_markers
truth <- shen_orr_ex$annotation$mixture
pure_samples <- lapply(1:3, function(i) {
    which(truth[, i] == 1)
})
Y <- shen_orr_ex$data$log
n_choose <- c(10, 11, 12)

dt_out <- dtangle(Y, pure_samples = pure_samples, n_choose = n_choose, data_type = "microarray-gene", 
    marker_method = "ratio")

saveRDS(dt_out, "basic_markers.rds")

# gamma
truth <- shen_orr_ex$annotation$mixture
pure_samples <- lapply(1:3, function(i) {
    which(truth[, i] == 1)
})
Y <- shen_orr_ex$data$log
n_choose <- c(10, 11, 12)

dt_out <- dtangle(Y, pure_samples = pure_samples, n_choose = n_choose, data_type = "microarray-gene", 
    marker_method = "ratio", gamma = 1)

saveRDS(dt_out, "basic_gamma_1.rds")

dt_out <- dtangle(Y, pure_samples = pure_samples, n_choose = n_choose, data_type = "rna-seq", 
    marker_method = "ratio")
saveRDS(dt_out, "basic_gamma_seq.rds")

# markers
dt_out <- dtangle(Y, pure_samples = pure_samples, n_choose = n_choose, marker_method = "regression", 
    gamma = 1)
saveRDS(dt_out, "basic_marker_reg.rds")

