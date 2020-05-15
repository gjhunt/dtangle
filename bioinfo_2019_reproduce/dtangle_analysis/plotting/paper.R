source("util/util.R", chdir = TRUE)

## Main Paper Plots Read in Error Data
paper_res <- readRDS("../comparison/Analysis/paper/paper.rds")
paper_res <- reduce_newman_res(paper_res)
paper_E <- all_errors(paper_res)

pdir <- "paper_plots/"

# Plots for each data set
dp_list <- list()
for (dset in levels(paper_E$Dataset)) {
    dp_list[[dset]] <- plots_of_dset(paper_E, dset, height = list(5, 5, 5, 5, 10), 
        width = c(5, 5, 5, 5, 10))
}

# Meta boxplot over datasets
ma_dsets <- c("Abbas", "Becht", "Gong", "Kuhn", "Newman FL", "Newman PBMC", "Shi", 
    "Shen-Orr")
seq_dsets <- c("Linsley", "Liu", "Parsons")
all_dsets <- c(ma_dsets, seq_dsets)

for (dset in c("all", "ma", "seq")) {
    E <- paper_E[paper_E$Dataset %in% switch(dset, all = all_dsets, ma = ma_dsets, 
        seq = seq_dsets), ]
    
    box <- (dset != "seq")
    
    pl <- list(make_meta_boxplot(E, useCorr = TRUE, box = box), make_meta_boxplot(E, 
        useCorr = TRUE, R2 = TRUE, box = box), make_meta_boxplot(E, box = box))
    save_plots(pl[[1]], dir = pdir, name = paste0("meta_boxplots_cor_", dset), subdir = "meta_plots/", 
        width = 10, height = 4)
    save_plots(pl[[2]], dir = pdir, name = paste0("meta_boxplots_r2_", dset), subdir = "meta_plots/", 
        width = 10, height = 4)
    save_plots(pl[[3]], dir = pdir, name = paste0("meta_boxplots_err_", dset), subdir = "meta_plots/", 
        width = 10, height = 4)
    pl[[1]] <- pl[[1]] + theme(legend.position = "none", legend.title = element_blank(), 
        axis.title.x = element_blank())
    pl[[3]] <- pl[[3]] + theme(legend.position = "none", legend.title = element_blank(), 
        axis.title.x = element_blank())
    leg <- g_legend(pl[[2]])
    
    pl[[2]] <- pl[[2]] + theme(legend.position = "none", legend.title = element_blank(), 
        axis.title.x = element_blank())
    
    gplts <- plot_grid(plotlist = pl, ncol = 3, labels = "AUTO")
    gplts <- plot_grid(gplts, leg, nrow = 2, rel_heights = c(1, 0.1))
    save_plots(gplts, dir = pdir, name = paste0("meta_boxplots_", dset), width = 12, 
        height = 4.5)
}

# Time Boxplots
T <- times(paper_res)
save_plots(make_time_boxplot(T), dir = pdir, name = "time_meta_boxplots", width = 10, 
    height = 4)

## Case-specific microarray mixture experiments
cin <- function(x, y) {
    tvals <- sapply(y, function(z) abs(z - x) < 1e-15)
    apply(tvals, 1, any)
}

## Quantile Assessment
sq_res <- readRDS("../comparison/Analysis/q_sens/q_sens.rds")
sq_res <- reduce_newman_res(sq_res)
sq_E_all <- all_errors(sq_res)
sq_E_all$Quantile <- 1 - sq_E_all$Quantile
sq_E_ma <- sq_E_all[sq_E_all$Dataset %in% ma_dsets, ]
sq_E_seq <- sq_E_all[sq_E_all$Dataset %in% seq_dsets, ]

byV <- "Quantile"
splitV <- "Deconv.Method"
for (dset in c("all", "ma", "seq")) {
    E <- get(paste0("sq_E_", dset))
    q_pl <- list()
    for (useCorr in c(FALSE, TRUE)) {
        for (R2 in c(TRUE, FALSE)) {
            if (!useCorr && R2) 
                next
            
            for (smry in c("mean", "median")) {
                
                q_pl[[length(q_pl) + 1]] <- by_var(E, byV = byV, splitV = splitV, 
                  ribbon = TRUE, useCorr = useCorr, R2 = R2, summ = switch(smry, 
                    mean = mn, median = mdn), title = "Error by Marker Selection", 
                  ylab = paste0("Grand ", capitalize(smry), " of ", ifelse(useCorr, 
                    ifelse(R2, "R2", "Cor."), "Errors")), xlab = "Quantile Cutoff for Markers", 
                  subtitle = switch(dset, all = "All Datasets", ma = "Microarray Datasets", 
                    seq = "Sequencing Datasets"), expandY = as.numeric(useCorr))
                save_plots(q_pl[[length(q_pl)]], dir = pdir, name = paste0(dset, 
                  "_", smry, "_", ifelse(useCorr, ifelse(R2, "r2", "cor"), "err"), 
                  "_by_quantile"), width = 12, height = 3.5, subdir = "quantile/")
            }
        }
    }
    q_pl <- lapply(q_pl, function(x) x + theme(plot.title = element_blank(), plot.subtitle = element_blank()))
    q_pl[1:(length(q_pl) - 1)] <- lapply(q_pl[1:(length(q_pl) - 1)], function(x) x + 
        theme(axis.title.x = element_blank(), legend.position = "none"))
    save_plots(plot_grid(plotlist = q_pl, labels = "AUTO", ncol = 1), dir = pdir, 
        name = dset, width = 14, height = 15, subdir = "quantile/")
}
