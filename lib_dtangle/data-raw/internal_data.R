set.seed(3251991)
source("./affymetrix-micro/estimate_gamma.R", local = TRUE, chdir = TRUE)
set.seed(3151990)
source("./illumina-seq/estimate_gamma.R", local = TRUE, chdir = TRUE)

gma <- list(ma_probe = p_gamma, ma_gene = g_gamma, rna_seq = s_gamma)
devtools::use_data(gma, internal = TRUE, overwrite = TRUE)
