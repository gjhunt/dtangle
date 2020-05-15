plots_of_dset <- function(E, dset = NULL, save = TRUE, subdir = "datasets/", adjust = NULL, 
    height = list(5, 5, 5, 5, 10), width = c(5, 5, 5, 5, 10)) {
    if (length(unique(E$Dataset)) == 1) 
        dset <- unique(E$Dataset)
    E <- E[E$Dataset == dset, ]
    pl <- list()
    pl$err <- try(dataset_by_quantile_boxplot(E, dset, full_page = FALSE)[[1]])
    pl$cor <- try(corr_boxplot(E, dset))
    pl$r2 <- try(corr_boxplot(E, dset, R2 = TRUE))
    E <- E[E$Deconv.Method %in% c("dtangle", "CIBERSORT", "EPIC"), ]
    pl$scatter <- try(plot_scatter(E, dset))
    if (!is.null(adjust)) {
        if (length(adjust) == 1) 
            adjust <- rep(adjust, length(pl))
        nms <- names(pl)
        pl <- lapply(seq_along(pl), function(i) pl[[i]] + adjust[[i]])
        names(pl) <- nms
    }
    pl$all <- try(plot_grid(plotlist = pl, nrow = 2, labels = "AUTO"))
    lapply(seq_along(pl), function(i) save_plots(pl[[i]], dir = pdir, subdir = subdir, 
        name = paste0(gsub(" ", "_", dset), "_", names(pl)[[i]]), width = width[[i]], 
        height = height[[i]]))
    return(pl)
}
