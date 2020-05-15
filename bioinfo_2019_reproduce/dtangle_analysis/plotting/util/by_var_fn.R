
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
    xlab = NA, ribbon = TRUE, pfrac = 0.05, byV, splitV, ribbonV = "Dataset", facetV = "Marker.Method", 
    useCorr = FALSE, R2 = FALSE, expandY = 0, confidence = 0.95) {
    
    E[[splitV]] <- as.factor(E[[splitV]])
    E <- data.table(E)
    split_E <- E[, .(Error = mn(Error)), by = .(get(ribbonV), get(splitV), get(byV), 
        get(facetV))]
    if (useCorr) {
        sds <- E[, .(sd = sd(Truth)), by = .(get(ribbonV), get(splitV), get(byV), 
            get(facetV), Cell.Type)]
        if (R2) {
            split_E <- E[, .(Error = cor(Estimate, Truth)^2), by = .(get(ribbonV), 
                get(splitV), get(byV), get(facetV), Cell.Type)]
        } else {
            split_E <- E[, .(Error = cor(Estimate, Truth)), by = .(get(ribbonV), 
                get(splitV), get(byV), get(facetV), Cell.Type)]
        }
        
        split_E <- split_E[sds$sd != 0, ]
        split_E[is.na(split_E$Error), "Error"] <- 0
    }
    colnames(split_E)[1:4] <- c(ribbonV, splitV, byV, facetV)
    
    N <- E[, .(N = length(Error)), by = .(get(ribbonV), get(splitV), get(byV), get(facetV))]$N[1]
    multiplier <- abs(qnorm(1 - confidence))/2 * 1/sqrt(N)
    
    if (confidence < 1) {
        trend_E <- split_E[, .(Error = summ(Error), lwr = summ(Error) - multiplier * 
            spread(Error), upr = summ(Error) + multiplier * spread(Error)), by = .(get(splitV), 
            get(byV), get(facetV))]
    } else {
        trend_E <- split_E[, .(Error = summ(Error), lwr = min(Error), upr = max(Error)), 
            by = .(get(splitV), get(byV), get(facetV))]
        
    }
    colnames(trend_E) <- c(splitV, byV, facetV, "Error", "lwr", "upr")
    
    uq <- unique(trend_E[[byV]])
    ss <- uq[c(seq(1, length(uq), max(1, round(length(uq) * pfrac))), length(uq))]
    ss <- which(trend_E[[byV]] %in% ss)
    
    p <- ggplot(trend_E, aes_string(x = byV, y = "Error", color = splitV))
    if (ribbon) {
        p <- p + geom_ribbon(aes_string(fill = splitV, x = byV, ymin = "lwr", ymax = "upr"), 
            linetype = 2, alpha = 0.25, inherit.aes = FALSE)
    }
    p <- p + geom_line(size = 1, aes_string(linetype = splitV)) + geom_point(data = trend_E[ss, 
        ], aes_string(x = byV, y = "Error", shape = splitV), size = 3, stroke = 1)
    p <- p + facet_grid(as.formula(paste0(".~", facetV)))
    p <- p + scale_shape_manual(values = c(15:18, 0:14))
    p <- p + labs(title = title, subtitle = subtitle) + expand_limits(y = expandY)
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_fill_brewer(palette = "Set1")
    
    if (!is.na(ylab)) 
        p <- p + ylab(ylab)
    if (!is.na(ylab)) 
        p <- p + xlab(xlab)
    
    return(p + global_theme)
}
