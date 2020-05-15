

summ <- function(x) {
    mean(abs(x[abs(x) < 1]), na.rm = TRUE)
}

gen_boxplot <- function(aE, title = "", subtitle = "", ylab = NA, xlab = NA, jitter = "Dataset", 
    yvar = "Error", ylim = NULL, full_page = FALSE, facet_formula = NULL, box = TRUE, 
    jitter_width = 1) {
    bpf <- function(x) {
        r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
        names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
        r
    }
    
    bpf_mean <- function(x) {
        r <- rep(mean(x), 5)
        names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
        r
    }
    
    p <- ggplot(aE, aes_string("Deconv.Method", yvar))
    if (box) {
        p <- p + stat_summary(fun.data = bpf, geom = "boxplot", size = 1)
    }
    if (length(unique(aE[[jitter]])) > 1) {
        p <- p + stat_summary(fun.data = bpf_mean, geom = "boxplot", size = 1, color = "grey", 
            lty = 1)
    }
    
    print(jitter_width)
    print(jitter)
    
    p_size <- 3
    p_stroke <- 0.75
    if (jitter == "Dataset") {
        # nudges <- c(0, 0.15, -0.15, 0.05, -0.05, 0.2, -0.2, -0.1, 0.1, 0,
        # 0)[1:length(unique(aE$Dataset))] * jitter_width nudges <-
        # c(0,-.25,-.25,-.15,-.15,rep(.25,5))[1:length(unique(aE$Dataset))] print(nudges)
        # p <- p + geom_jitter(aes_string(color = jitter, shape = jitter), size = p_size,
        # stroke = p_stroke, position = position_nudge(x = nudges))
        p <- p + geom_jitter(aes_string(color = jitter, shape = jitter), stroke = p_stroke, 
            size = p_size, width = jitter_width * 0.2)
        # p = p + geom_segment(aes_string(color =
        # jitter,x='Deconv.Method',xend='Deconv.Method',y='MinError',yend='MaxError'),size=.5,position=position_nudge(x=nudges))
    } else {
        p <- p + geom_jitter(aes_string(color = jitter, shape = jitter), size = p_size, 
            stroke = p_stroke, wdith = jitter_width * 0.2)
    }
    p <- p + scale_shape_manual(values = 0:25)
    if (!any(is.null(ylim))) 
        p <- p + coord_cartesian(ylim = ylim)
    if (!is.na(subtitle)) {
        p <- p + labs(title = title, subtitle = subtitle)
    } else {
        p <- p + labs(title = title)
    }
    
    if (!is.na(ylab)) 
        p <- p + ylab(ylab)
    if (!is.na(ylab)) 
        p <- p + xlab(xlab)
    
    if (!is.null(facet_formula)) {
        p <- p + facet_grid(as.formula(facet_formula), scales = "free_x", space = "free_x")
    }
    
    p <- p + global_theme + theme(legend.text = element_text(size = global_size))
    if (jitter != "Dataset" & !full_page) {
        lj <- length(unique(aE[[jitter]]))
        if (lj <= 10) 
            p <- p + guides(shape = guide_legend(nrow = 1, byrow = TRUE), color = guide_legend(nrow = 1, 
                byrow = TRUE), legend.position = c(0.45, -0.33))
        if (lj > 10) 
            p <- p + theme(legend.position = "none")
    }
    
    return(p)
}

dataset_by_quantile_boxplot <- function(err, dsname, ttl = NA, full_page = FALSE, 
    yl = NULL) {
    err <- err[err$Dataset == dsname, ]
    uq <- unique(err$Quantile)
    um <- unique(err$Marker.Method)
    L <- list()
    for (q in uq) {
        E <- err[abs(err$Quantile - q) < 1e-10, ]
        E <- E[complete.cases(E), ]
        E$Error <- abs(E$Error)
        st <- paste0("Marker Method=", paste(um, collapse = ","), " Slope=Auto", 
            " Quantile=", q)
        if (is.na(ttl)) 
            ttl <- dsname
        yl <- c(0, 1)
        # if(max(E$Error) < .5) yl = c(0,.5)
        L[[toString(q)]] <- gen_boxplot(E, jitter = "Cell.Type", title = ttl, subtitle = "", 
            ylim = yl, full_page = full_page) + theme(legend.title = element_blank())
    }
    return(L)
}

scale_boxplot <- function(E, ttl = NA, full_page = FALSE, yl = NULL, separate_dtangle = FALSE) {
    E <- E[complete.cases(E), ]
    E$Error <- abs(E$Error)
    if (separate_dtangle) {
        E <- E[-which(E$Deconv.Method == "dtangle" & E$Scale == "linear"), ]
        E[which(E$Deconv.Method == "dtangle" & E$Scale == "log"), "Scale"] <- ""
    }
    yl <- c(0, 1)
    gen_boxplot(E, jitter = "Cell.Type", title = ttl, subtitle = "", ylim = yl, full_page = full_page, 
        facet_formula = ".~Scale")
}

make_meta_boxplot <- function(E, ylim = NULL, useCorr = FALSE, R2 = FALSE, jitter = "Dataset", 
    facet_formula = NULL, box = TRUE) {
    if (any(is.null(ylim))) 
        ylim <- c(0, 0.6)
    
    if (useCorr) {
        E_adj <- data.table(E)
        
        sds <- E_adj[, .(sd = sd(Truth)), by = .(Dataset, Deconv.Method, Marker.Method, 
            Quantile, Cell.Type)]
        if (R2) 
            aE <- E_adj[, .(Error = cor(Estimate, Truth)^2), by = .(Dataset, Deconv.Method, 
                Marker.Method, Quantile, Cell.Type)] else aE <- E_adj[, .(Error = cor(Estimate, Truth)), by = .(Dataset, Deconv.Method, 
            Marker.Method, Quantile, Cell.Type)]
        
        aE <- aE[sds$sd != 0, ]
        aE[is.na(aE$Error), "Error"] <- 0
        aE <- aE[, .(Error = median(Error, na.rm = TRUE), MinError = min(Error, na.rm = TRUE), 
            MaxError = max(Error, na.rm = TRUE)), by = .(Dataset, Deconv.Method, 
            Marker.Method, Quantile)]
        aE <- as.data.frame(aE)
        ylim <- c(-0.5, 1)
        if (R2) {
            ylab <- "R-squared"
            title <- "R-squared by Deconvolution Method"
            ylim <- c(0, 1)
        } else {
            ylab <- "Correlation"
            title <- "Correlation by Deconvolution Method"
        }
        
    } else {
        ylab <- "Mean Error"
        title <- "Mean Error by Deconvolution Method"
        E_adj <- data.table(E)
        aE <- E_adj[, .(Error = mean(abs(Error))), by = .(Dataset, Deconv.Method, 
            Marker.Method, Quantile, Cell.Type)]
        aE <- aE[, .(Error = mean(abs(Error), na.rm = TRUE), MinError = min(abs(Error), 
            na.rm = TRUE), MaxError = max(abs(Error), na.rm = TRUE)), by = .(Dataset, 
            Deconv.Method, Marker.Method, Quantile)]
        aE <- as.data.frame(aE)
    }
    
    gen_boxplot(aE, title = title, subtitle = "", ylab = ylab, xlab = "Deconvolution Method", 
        jitter = jitter, ylim = ylim, facet_formula = facet_formula, box = box, jitter_width = as.numeric(box))
}

corr_boxplot <- function(E, ttl = "", R2 = FALSE) {
    E_adj <- data.table(E)
    sds <- E_adj[, .(sd = sd(Truth)), by = .(Dataset, Deconv.Method, Marker.Method, 
        Quantile, Cell.Type, Scale)]
    if (R2) 
        aE <- E_adj[, .(Error = cor(Estimate, Truth)^2), by = .(Dataset, Deconv.Method, 
            Marker.Method, Quantile, Cell.Type, Scale)] else aE <- E_adj[, .(Error = cor(Estimate, Truth)), by = .(Dataset, Deconv.Method, 
        Marker.Method, Quantile, Cell.Type, Scale)]
    
    aE <- aE[sds$sd != 0, ]
    aE[is.na(aE$Error), "Error"] <- 0
    aE <- as.data.frame(aE)
    if (any(aE$Error < 0)) 
        ylim <- c(-1, 1) else ylim <- c(0, 1)
    if (length(unique(aE$Scale)) > 1) 
        ff <- ".~Scale" else ff <- NULL
    ylab <- "Correlation"
    if (R2) 
        ylab <- "R-squared"
    gen_boxplot(aE, title = ttl, subtitle = "", ylab = ylab, xlab = "Deconvolution Method", 
        ylim = ylim, facet_formula = ff, box = FALSE, jitter = "Cell.Type", jitter_width = 0.25)
}

make_time_boxplot <- function(T, ylim = NULL) {
    
    st <- ""
    aE <- T[, .(Error = mean((LogMinutes))), by = .(Dataset, Deconv.Method, Marker.Method, 
        Quantile)]
    aE <- aE[, .(Error = mean((Error), na.rm = TRUE), MinError = min((Error), na.rm = TRUE), 
        MaxError = max((Error), na.rm = TRUE)), by = .(Dataset, Deconv.Method, Marker.Method, 
        Quantile)]
    
    gen_boxplot(aE, title = "Computation Time by Deconvolution Method", subtitle = "", 
        ylab = "Computation Time Log10(Minutes)", xlab = "Deconvolution Method", 
        jitter = "Dataset", ylim = ylim) + scale_color_brewer(palette = "Set1") + 
        scale_fill_brewer(palette = "Set1")
}





