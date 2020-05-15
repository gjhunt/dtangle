source("util/util.R", chdir = TRUE)
pdir <- "./sim"
dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
separate_dtangle <- TRUE

## Simple Simulation Plots
for (err_dist in c("gaussian", "poisson")) {
    for (nm in c("lowe", "highe", "outlier")) {
        long_nm <- switch(nm, lowe = "Low Error", highe = "High Error", outlier = "Outliers")
        f_nm <- gsub(" ", "_", tolower(long_nm))
        this_pdir <- paste0(paste(pdir, err_dist, f_nm, sep = "/"), "/")
        res <- readRDS(paste0("../comparison/Analysis/sim_", err_dist, "_", nm, "/sim_", 
            err_dist, "_", nm, ".rds"))
        E <- all_errors(res, separate_dtangle = separate_dtangle)
        ttl <- paste0("Linear Scale Mixing with ", long_nm)
        
        pl <- list()
        pl$err <- scale_boxplot(E, ttl = ttl)
        pl$cor <- corr_boxplot(E, ttl = ttl)
        pl$r2 <- corr_boxplot(E, ttl = ttl, R2 = TRUE)
        pl$scatter <- plot_scatter(E, title = ttl, methods = c("dtangle", "CIBERSORT", 
            "EPIC"), facet_formula = ".~Scale")
        pl$all <- plot_grid(plotlist = pl, nrow = 2, labels = "AUTO")
        lapply(seq_along(pl), function(i) save_plots(pl[[i]], dir = this_pdir, name = paste0(f_nm, 
            "_", err_dist, "_", names(pl)[[i]])))
    }
    
    this_pdir <- paste0(pdir, "/", err_dist, "/")
    
    # Markers Near Zero
    mnz <- list()
    mv_res <- readRDS(paste0("../comparison/Analysis/sim_", err_dist, "_marker_value/sim_", 
        err_dist, "_marker_value.rds"))
    mv_E <- all_errors(mv_res, separate_dtangle = FALSE)
    mnz[[1]] <- save_plots(by_var(mv_E, byV = "marker_quantile_low", splitV = "Deconv.Method", 
        ribbonV = "Cell.Type", pfrac = 0, title = "Accuacy of dtangle v. Expression of Marker Gene by Other Cell Types", 
        confidence = 1) + xlab("Expression of Marker Genes by Other Cell Types"), 
        this_pdir, name = paste0("zero_marker_assumption_", err_dist))
    mnz[[2]] <- save_plots(by_var(mv_E, byV = "marker_quantile_low", splitV = "Deconv.Method", 
        ribbonV = "Cell.Type", pfrac = 0, title = "Accuacy of dtangle v. Expression of Marker Gene by Other Cell Types", 
        confidence = 1, useCorr = TRUE, ylab = "Corr.") + xlab("Expression of Marker Genes by Other Cell Types"), 
        this_pdir, name = paste0("zero_marker_assumption_corr_", err_dist))
    mnz[[3]] <- save_plots(by_var(mv_E, byV = "marker_quantile_low", splitV = "Deconv.Method", 
        ribbonV = "Cell.Type", pfrac = 0, title = "Accuacy of dtangle v. Expression of Marker Gene by Other Cell Types", 
        confidence = 1, useCorr = TRUE, R2 = TRUE, ylab = "R-sqared") + xlab("Expression of Marker Genes by Other Cell Types"), 
        this_pdir, name = paste0("zero_marker_assumption_r2_", err_dist))
    save_plots(plot_grid(plotlist = mnz, labels = "AUTO", nrow = 3), dir = this_pdir, 
        name = "nmz")
    
    # Num of Markers
    nmk <- list()
    num_res <- readRDS(paste0("../comparison/Analysis/sim_", err_dist, "_num_markers/sim_", 
        err_dist, "_num_markers.rds"))
    num_E <- all_errors(num_res, separate_dtangle = FALSE)
    nmk[[1]] <- save_plots(by_var(num_E, byV = "pct_markers", splitV = "Deconv.Method", 
        ribbonV = "Cell.Type", pfrac = 0, confidence = 1, title = "Accuacy of dtangle v. Pct. of Genes that are Markers") + 
        xlab("Pct. of Genes that are Markers"), this_pdir, name = paste0("num_markers_", 
        err_dist))
    num_res <- readRDS(paste0("../comparison/Analysis/sim_", err_dist, "_num_markers/sim_", 
        err_dist, "_num_markers.rds"))
    num_E <- all_errors(num_res, separate_dtangle = FALSE)
    nmk[[2]] <- save_plots(by_var(num_E, byV = "pct_markers", splitV = "Deconv.Method", 
        ribbonV = "Cell.Type", confidence = 1, pfrac = 0, useCorr = TRUE, ylab = "Corr", 
        title = "Accuacy of dtangle v. Pct. of Genes that are Markers") + xlab("Pct. of Genes that are Markers"), 
        this_pdir, name = paste0("num_markers_corr_", err_dist))
    nmk[[3]] <- save_plots(by_var(num_E, byV = "pct_markers", splitV = "Deconv.Method", 
        ribbonV = "Cell.Type", confidence = 1, pfrac = 0, useCorr = TRUE, R2 = TRUE, 
        ylab = "R-squared", title = "Accuacy of dtangle v. Pct. of Genes that are Markers") + 
        xlab("Pct. of Genes that are Markers"), this_pdir, name = paste0("num_markers_r2_", 
        err_dist))
    save_plots(plot_grid(plotlist = nmk, labels = "AUTO", nrow = 3), dir = this_pdir, 
        name = "nmk")
    
}


# Mixture Simulation Plots
for (err_dist in c("mix_gaussian", "mix_poisson")) {
    for (dset in c("linsley", "parsons")) {
        for (nm in c("lowe", "highe", "outlier")) {
            tryCatch({
                long_nm <- switch(nm, lowe = "Low Error", highe = "High Error", outlier = "Outliers")
                f_nm <- gsub(" ", "_", tolower(long_nm))
                this_pdir <- paste0(paste(pdir, err_dist, dset, f_nm, sep = "/"), 
                  "/")
                res <- readRDS(paste0("../comparison/Analysis/sim_", err_dist, "_", 
                  nm, "_", dset, "/sim_", err_dist, "_", nm, "_", dset, ".rds"))
                E <- all_errors(res, separate_dtangle = separate_dtangle)
                ttl <- paste0("Linear Scale Mixing with ", long_nm)
                pl <- list()
                
                pl$err <- scale_boxplot(E, ttl = ttl)
                pl$cor <- corr_boxplot(E, ttl = ttl)
                pl$r2 <- corr_boxplot(E, ttl = ttl, R2 = TRUE)
                pl$scatter <- plot_scatter(E, title = ttl, methods = c("dtangle", 
                  "CIBERSORT", "EPIC"), facet_formula = ".~Scale")
                pl$all <- plot_grid(plotlist = pl, nrow = 2, labels = "AUTO")
                lapply(seq_along(pl), function(i) save_plots(pl[[i]], dir = this_pdir, 
                  name = paste0(f_nm, "_", err_dist, "_", names(pl)[[i]])))
                
            }, error = function(e) {
                cat(paste(e), "\n")
            })
        }
    }
}

