#' Row-binds \code{Y} with \code{references} and generates \code{pure_samples}.
#' @inheritParams dtangle
combine_Y_refs <- function(Y, references, pure_samples) {
    
    if (is.null(pure_samples)) {
        pure_samples <- lapply(1:nrow(references), identity)
    }
    
    Y <- as.matrix(rbind(references, Y))
    
    return(list(Y = Y, pure_samples = pure_samples))
}
