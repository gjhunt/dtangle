source("util/util.R", chdir = TRUE)

ma_dsets <- c("Abbas", "Becht", "Gong", "Kuhn", "Newman FL", "Newman PBMC", "Shi", 
    "Shen-Orr")
seq_dsets <- c("Linsley", "Liu", "Parsons")
all_dsets <- c(ma_dsets, seq_dsets)

for (err_struct in c("gaussian", "poisson")) {
    for (err_level in c("low", "high")) {
        
        ## Quantile Assessment
        sq_res <- readRDS(paste0("../comparison/Analysis/markers_sim_", err_struct, 
            "_", err_level, "e/markers_sim_", err_struct, "_", err_level, "e.rds"))
        sq_E_all <- all_errors(sq_res)
        sq_E_all$Quantile <- 1 - sq_E_all$Quantile
        sq_E_ma <- sq_E_all[sq_E_all$Dataset %in% ma_dsets, ]
        sq_E_seq <- sq_E_all[sq_E_all$Dataset %in% seq_dsets, ]
        
        pdir <- paste0("./marker_sim/", err_struct, "/", err_level, "/")
        
        byV <- "Quantile"
        splitV <- "Deconv.Method"
        for (dset in c("all")) {
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
            q_pl <- lapply(q_pl, function(x) x + theme(plot.title = element_blank(), 
                plot.subtitle = element_blank()))
            q_pl[1:(length(q_pl) - 1)] <- lapply(q_pl[1:(length(q_pl) - 1)], function(x) x + 
                theme(axis.title.x = element_blank(), legend.position = "none"))
            save_plots(plot_grid(plotlist = q_pl, labels = "AUTO", ncol = 1), dir = pdir, 
                name = dset, width = 14, height = 15, subdir = "quantile/")
        }
    }
}
