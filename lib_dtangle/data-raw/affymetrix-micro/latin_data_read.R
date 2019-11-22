## ---- latinread
library("affy")

affy_batch <- ReadAffy(celfile.path = "latin_data")

## probe-level data
ips <- indexProbes(affy_batch, "both")
d <- log2(t(intensity(affy_batch)))

# order according to experiment number
ord <- order(sapply(rownames(d), function(x) {
    as.integer(substr(strsplit(x, "_")[[1]][8], 5, 10))
}))
d <- d[ord, ]

# gene-level data
dg <- data.frame(rma(affy_batch, verbose = FALSE, normalize = TRUE, background = TRUE))
dg <- dg[ord, ]

# determine gene spike-in titrations
anno <- read.csv("anno.csv", stringsAsFactors = FALSE)
sp_genes <- sapply(anno[1, 2:15], function(x) {
    strsplit(x, "\n", fixed = TRUE)
})
titration <- as.matrix(anno[2:15, 2:15])
titration <- apply(titration, c(1, 2), as.numeric)

# rows are spiked-in genes, columns are experiment
all_titration <- array(0, c(dim(titration)[1] * 3, dim(titration)[2]))
for (i in 1:dim(all_titration)[1]) {
    all_titration[i, ] <- titration[ceiling(i/3), ]
}
all_titration_order <- apply(all_titration, 2, order)

group <- 1:14
gene <- 1:3
sd <- d
for (g in group) {
    probe_groups <- unlist(ips[unlist(lapply(sp_genes[g], "[", gene))])
    sd[, probe_groups] <- sd[all_titration_order[, g], probe_groups]
}
all_probes <- unlist(ips[unlist(lapply(sp_genes[group], "[", gene))])
sd <- sd[, all_probes]

sdg <- dg
for (g in group) {
    gns <- unlist(make.names(sp_genes[[g]]))
    sdg[, gns] <- sdg[all_titration_order[, g], gns]
}
all_genes <- make.names(unlist(sp_genes))
sdg <- sdg[, all_genes]





