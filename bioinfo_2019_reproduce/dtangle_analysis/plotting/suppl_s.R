source("util/util.R", chdir = TRUE)

ma_dsets <- c("Abbas", "Becht", "Gong", "Kuhn", "Newman FL", "Newman PBMC", "Shi", 
    "Shen-Orr")
seq_dsets <- c("Linsley", "Liu", "Parsons")

ss_res <- readRDS("../comparison/Analysis/supplement_s/supplement_s.rds")
ss_res <- reduce_newman_res(ss_res)

ss_E_all <- all_errors(ss_res)
ss_E_ma <- ss_E_all[ss_E_all$Dataset %in% ma_dsets, ]
ss_E_seq <- ss_E_all[ss_E_all$Dataset %in% seq_dsets, ]

byV <- "gamma"
splitV <- "Quantile"

slope_pl <- list()
for (dset in c("all", "ma", "seq")) {
    E <- get(paste0("ss_E_", dset))
    s_pl <- list()
    for (useCorr in c(FALSE, TRUE)) {
        for (R2 in c(TRUE, FALSE)) {
            if (!useCorr && R2) 
                next
            
            for (smry in c("mean", "median")) {
                
                s_pl[[length(s_pl) + 1]] <- by_var(E, byV = byV, splitV = splitV, 
                  useCorr = useCorr, R2 = R2, ribbon = TRUE, summ = switch(smry, 
                    mean = mn, median = mdn), title = "Meta Error by Marker Slope", 
                  ylab = paste0(capitalize(smry), " of ", ifelse(useCorr, ifelse(R2, 
                    "R2", "Cor."), "Errors")), xlab = "Slope for Our Method", subtitle = switch(dset, 
                    all = "All Datasets", ma = "Microarray Datasets", seq = "Sequencing Datasets"), 
                  expandY = as.numeric(useCorr))
                save_plots(s_pl[[length(s_pl)]], dir = "ss/", name = paste0(dset, 
                  "_", smry, "_", ifelse(useCorr, ifelse(R2, "r2", "cor"), "err"), 
                  "_by_slope"), width = 10, height = 4)
            }
        }
    }
    s_pl <- lapply(s_pl, function(x) x + theme(plot.title = element_blank(), plot.subtitle = element_blank()))
    s_pl[1:(length(s_pl) - 1)] <- lapply(s_pl[1:(length(s_pl) - 1)], function(x) x + 
        theme(axis.title.x = element_blank(), legend.position = "none"))
    save_plots(plot_grid(plotlist = s_pl, labels = "AUTO", ncol = 1), dir = "ss/", 
        name = dset, width = 14, height = 15, subdir = NULL)
}
