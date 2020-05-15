library("RColorBrewer")

plot_scatter <- function(E, title = "", subtitle = NA, ylab = NA, xlab = NA, methods = NULL, 
    facet_formula = NULL) {
    
    if (!is.null(methods)) {
        E <- E[E$Deconv.Method %in% methods, ]
    }
    
    net <- E$Truth + runif(length(E$Truth), 0, 0.02)
    E$Truth <- net
    
    p_size <- 2.5
    p_stroke <- 1
    p <- ggplot(E, aes(x = Truth, y = Estimate, shape = Deconv.Method, color = Deconv.Method))
    p <- p + ylim(0, 1) + xlim(0, 1)
    p <- p + geom_abline(intercept = 0, slope = 1, lwd = 0.5, color = "Black")
    p <- p + geom_point(size = p_size, stroke = p_stroke, alpha = 0.75)
    if (length(unique(E$Marker.Method)) > 1) {
        p <- p + facet_grid(. ~ Marker.Method)
    }
    
    p <- p + scale_shape_manual(values = 0:25, name = "Deconv\nMethod")
    
    if (!is.null(facet_formula)) 
        p <- p + facet_grid(as.formula(facet_formula), space = "free_x")
    
    group_list <- lapply(unique(E$Deconv.Method), function(x) E$Deconv.Method == 
        x)
    names(group_list) <- unique(E$Deconv.Method)
    
    if (!is.na(subtitle)) {
        p <- p + labs(title = title, subtitle = subtitle)
    } else {
        p <- p + labs(title = title)
    }
    
    if (!is.na(ylab)) 
        p <- p + ylab(ylab)
    if (!is.na(ylab)) 
        p <- p + xlab(xlab)
    
    lj <- length(unique(E$Cell.Type))
    
    p <- p + global_theme  # + theme(legend.title = element_blank(), #legend.position = c(.45,-.2),
    # plot.margin = unit(c(0,0,1,0), 'lines'))
    
    p <- p + scale_color_brewer("Deconv\nMethod", palette = "Dark2")
    
    return(p)
}
