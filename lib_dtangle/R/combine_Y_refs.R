#' Row-binds \code{Y} with \code{references} and generates \code{pure_samples}.
#' @inheritParams dtangle
combine_Y_refs <- function(Y, references, pure_samples) {
    
    if (is.null(pure_samples)) {
        pure_samples <- lapply(1:nrow(references), identity)
        names(pure_samples) <- rownames(references)
    }
    
    if (is.null(colnames(Y)) & !is.null(colnames(references))) {
        colnames(Y) <- colnames(references)
    }
    
    if (!is.null(references) & is.null(colnames(references))) 
        colnames(references) <- colnames(Y)
    
    if (!is.null(references)) 
        Y <- as.matrix(rbind(as.matrix(references), as.matrix(Y)))
    
    if (is.null(colnames(Y))) 
        colnames(Y) <- 1:ncol(Y)
    
    return(list(Y = Y, pure_samples = pure_samples))
}
