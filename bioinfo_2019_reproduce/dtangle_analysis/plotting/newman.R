source("util/util.R", chdir = TRUE)

## Main Paper Plots Read in Error Data
newman_res <- readRDS("../comparison/Analysis/newman/newman.rds")
newman_res <- reduce_newman_res(newman_res)
newman_E <- all_errors(newman_res)

pdir <- "newman_markers_plots/"

# Plots for each data set
dp_list <- list()
for (dset in levels(newman_E$Dataset)) {
    dp_list[[dset]] <- plots_of_dset(newman_E, dset, height = list(5, 5, 5, 5, 10), 
        width = c(5, 5, 5, 5, 10))
}
