library("crayon")

sze <- function(unit = "Mb", envir = globalenv()) {
    os <- function(x) {
        object.size(get(x))
    }
    os_pretty <- function(x) {
        format(object.size(get(x)), units = unit)
    }
    ord <- order(sapply(ls(envir), os), decreasing = TRUE)
    print(sapply(ls(envir), os_pretty)[ord])
    cat("Total: ", sum(sapply(ls(envir), os)) * 1e-09, "Gb \n")
}

VERBOSE <- TRUE
TIME <- Sys.time()
updt <- function(msg, init = FALSE) {
    if (VERBOSE) {
        t2 <- Sys.time()
        if (init) {
            TIME <<- t2
            FIRST <<- t2
        }
        MSG <- paste(msg, " (", format(t2 - FIRST, digits = 1), " ellapsed)", sep = "")
        TIME <<- t2
        cat(blue(paste0("[", t2, "]")), MSG, "\n")
        return(MSG)
    }
    return("")
}
