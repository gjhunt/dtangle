source("util/util.R", chdir = TRUE)

## Main Paper Plots Read in Error Data
paper_res <- readRDS("../comparison/Analysis/log_paper/log_paper.rds")
paper_res <- reduce_newman_res(paper_res)
paper_E <- all_errors(paper_res)

pdir <- "log_paper_plots/"

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
        height = 5)
}
