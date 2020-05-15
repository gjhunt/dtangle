save_plots <- function(pl, dir, subdir = NULL, name, uid = NULL, width = NULL, height = NULL) {
    if (is.null(width)) 
        width <- 15
    if (is.null(height)) 
        height <- 10
    graphics.off()
    dir <- paste0(dir, subdir)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    
    nm <- paste0(dir, name, uid)
    
    dvce_print(pl, pdf, nm, "pdf", width, height)
    dvce_print(pl, cairo_ps, nm, "eps", width, height, fallback_resolution = 600)
    
    return(pl)
}

dvce_print <- function(pl, dvce, nm, extension, width, height, ...) {
    fname <- paste0(nm, ".", extension)
    cat(paste0("Saving ", fname, "\n"))
    dvce(fname, width = width, height = height, ...)
    if ("ggplot" %in% class(pl)) {
        print(pl)
    } else {
        for (p in pl) {
            tryCatch({
                grid.draw(p)
            }, error = function(e) {
                print(e)
            })
        }
    }
    dev.off()
}
