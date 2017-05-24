ref <- read.csv("ERCC_Controls_Analysis.txt", sep = "\t", stringsAsFactors = FALSE)

data_dir <- dir("data/")

dtype <- "ILM"
fnames <- paste0("data/", data_dir[grepl(paste0("_", dtype, "_"), data_dir)])

read_data <- function(fname) {
    read.table(fname, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

d <- do.call(cbind, lapply(fnames, read_data))

which_series <- which(grepl("_E_", colnames(d)) | grepl("_F_", colnames(d)))
es <- d[, which_series]
rownames(es) <- gsub("_", "-", d$TranscriptID)

conc1 <- ref[, c(2, 4)]
colnames(conc1) <- c("id", "conc")

map <- sapply(1:dim(conc1)[1], function(i) {
    which(rownames(es) == conc1$id[i])
})

sp_es <- es[map, ]

log_conc <- log2(conc1$conc)
log_expr <- log2(1 + sp_es)

I <- array(1, c(1, dim(log_expr)[2]))
C <- kronecker(I, as.matrix(log_conc))

coefs <- rep(0, dim(log_expr)[2])
for (i in 1:dim(log_expr)[2]) {
    y <- as.vector(as.matrix(log_expr[, i]))
    x <- as.vector(C[, i])
    m <- lm(y ~ 1 + x)
    coefs[i] <- coef(m)[2]
}

s_gamma <- median(coefs)
