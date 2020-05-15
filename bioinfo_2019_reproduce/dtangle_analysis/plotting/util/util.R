library("cowplot")
library("data.table")
library("dtangle.data")
library("ggplot2")
library("grid")
library("Hmisc")
library("reshape2")
library("xtable")
source("./all_errors.R")
source("./by_var_fn.R")
source("./dset_plots.R")
source("./global_theme.R")
source("./plot_boxplots_fn.R")
source("./plot_scatter_fn.R")
source("./save_plots.R")

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
