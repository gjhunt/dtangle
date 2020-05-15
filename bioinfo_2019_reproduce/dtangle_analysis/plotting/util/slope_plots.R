mn <- function(x) {
    mean(abs(x[abs(x) < 1]), na.rm = TRUE)
}
tmn <- function(x) {
    mean(abs(x[abs(x) < 1]), trim = 0.25, na.rm = TRUE)
}
mdn <- function(x) {
    median(abs(x[abs(x) < 1]), na.rm = TRUE)
}
sdev <- function(x) (sd(abs(x[abs(x) < 1])))
iqr <- function(x) 0.5 * (quantile(x, 0.75) - quantile(x, 0.25))

by_var <- function(E, summ = mn, spread = sdev, title = "", subtitle = "", ylab = NA, 
    xlab = NA, ymax = NA, ribbon = TRUE, pfrac = 0.05, byV, splitV) {
    
    E[[splitV]] <- as.factor(E[[splitV]])
    
    a_vars <- colnames(E)[!(colnames(E) %in% c("Cell.Type", "Sample", "Error", "Estimate", 
        "Truth"))]
    aE <- aggregate(as.formula(paste0("Error ~ ", paste(a_vars, collapse = "+"))), 
        data = E, mn)
    t_vars <- a_vars[!(a_vars %in% c("Dataset"))]
    tE <- aggregate(as.formula(paste0("Error ~ ", paste(t_vars, collapse = "+"))), 
        data = aE, summ)
    
    uq <- unique(tE[[byV]])
    ss <- uq[seq(1, length(uq), max(1, round(length(uq) * pfrac)))]
    ss <- which(tE[[byV]] %in% ss)
    
    width <- aggregate(as.formula(paste0("Error ~ ", paste(t_vars, collapse = "+"))), 
        data = aE, spread)
    se_mean <- width$Error * abs(qnorm(0.1))/2 * 1/sqrt(length(unique(aE$Dataset)))
    width$lwr <- tE$Error - se_mean
    width$upr <- tE$Error + se_mean
    
    p <- ggplot(tE, aes_string(x = byV, y = "Error", color = splitV))
    
    if (ribbon) {
        p <- p + geom_ribbon(data = width, aes_string(fill = splitV, ymin = "lwr", 
            ymax = "upr"), linetype = 2, alpha = 0.1)
    }
    
    p <- p + geom_line(size = 1.5) + geom_point(data = tE[ss, ], aes_string(x = byV, 
        y = "Error", shape = splitV), size = 3, stroke = 1.5)
    p <- p + theme(plot.title = element_text(hjust = 0.5, size = 35), axis.text = element_text(size = 25), 
        axis.title = element_text(size = 30, face = "bold"), plot.subtitle = element_text(hjust = 0.5, 
            size = 20), legend.text = element_text(size = 18), legend.key.width = unit(2, 
            "line"), legend.key.height = unit(2, "line"), legend.title = element_text(size = 25), 
        legend.position = "bottom", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
        panel.spacing = unit(2, "lines"), strip.text = element_text(size = 25))
    p <- p + scale_shape_manual(values = c(15:18, 0:14))
    
    if (!is.na(ymax)) {
        p <- p + ylim(0, ymax)
    } else {
        p <- p + ylim(0, 0.25)
    }
    p <- p + labs(title = title, subtitle = subtitle)
    
    if (!is.na(ylab)) 
        p <- p + ylab(ylab)
    if (!is.na(ylab)) 
        p <- p + xlab(xlab)
    if (length(unique(E$Marker.Method)) > 1) {
        p <- p + facet_grid(. ~ Marker.Method)
    }
    return(p)
}
