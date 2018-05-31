get_marker_list <- function(value) {
    typ <- typeof(value[[1]])
    if (typ == "list") 
        return(value$L)
    return(value)
}
