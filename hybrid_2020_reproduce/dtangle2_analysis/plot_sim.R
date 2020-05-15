library("data.table")
library("ggplot2")
library("cowplot")

L = fl$L
L = L[sapply(L, lengths) > 1]

L = lapply(L, function(x) {
    ret = x$err
    ret[method == "dtangle2", "method"] <- "Hybrid"
    ret[method == "LogRegression", "method"] <- "Log Reg."
    ret[method == "ds", "method"] <- "deconvSeq"
    ret[method == "icedt", "method"] <- "ICeDT"
    ret = ret[method %in% kp_mths, ]
    ret$method <- factor(ret$method, levels = unique(ret$method))
    return(list(err = ret))
})
type = fl$type
grd = fl$grd

smry = lapply(L, "[[", "err")
for (i in 1:length(smry)) {
    smry[[i]]$i <- i
}
smry = do.call(rbind, smry)
grd$i = 1:nrow(grd)
d = merge(smry, grd)
d = data.table(d)
d$n_markers = paste0("N. Marker = ", d$n_markers)
d$n_markers = ordered(d$n_markers, levels = unique(d$n_marker))
d$mu = paste0("Seq. Depth = ", d$mu)
d$mu = ordered(d$mu, levels = unique(d$mu))

# Meta boxplots for all cases
plt = ggplot(data = d, mapping = aes(x = as.factor(sigma), y = err, group = interaction(sigma, 
    method), color = method, shape = method)) + geom_boxplot(size = 1, position = "dodge")
if (length(levels(d$mu)) > 1) {
    plt = plt + facet_grid(n_markers ~ mu, scales = "fixed") + scale_y_sqrt()
} else {
    plt = plt + facet_grid(n_markers ~ ., scales = "fixed") + scale_y_sqrt()
}
plt = plt + labs(x = "Tau", y = "Error") + theme_bw()
ggsave(paste0(pdir, type, "_boxplots.pdf"), plt, width = 14, height = 10)

dzo = d[(truth < 0.05) | (truth > 0.95), ]
plt = ggplot(data = dzo, mapping = aes(x = as.factor(sigma), y = err, group = interaction(sigma, 
    method), color = method, shape = method)) + geom_boxplot(size = 1, position = "dodge")
if (length(levels(d$mu)) > 1) {
    plt = plt + facet_grid(n_markers ~ mu, scales = "fixed") + scale_y_sqrt()
} else {
    plt = plt + facet_grid(n_markers ~ ., scales = "fixed") + scale_y_sqrt()
}
plt = plt + labs(x = "Tau", y = "Error") + theme_bw()
ggsave(paste0(pdir, type, "_nz_boxplots.pdf"), plt, width = 14, height = 10)

NA_lm = function(fmla, data) {
    tryCatch({
        lm(formula = fmla, data = data)$coef
    }, error = function(e) {
        as.numeric(c(NA, NA))
    })
}
lm_summary = d[, .(b0 = NA_lm(truth ~ 1 + est, data = .SD)[1], b1 = NA_lm(truth ~ 
    1 + est, data = .SD)[2]), by = .(method, n_markers, mu, sigma)]
lm_summary$b1 = lm_summary$b1 - 1
plt = ggplot(data = lm_summary, mapping = aes(x = as.factor(sigma), y = abs(b1), 
    group = interaction(sigma, method), color = method, shape = method)) + geom_jitter(size = 5, 
    width = 0.1) + scale_y_sqrt()
if (length(levels(d$mu)) > 1) {
    plt = plt + facet_grid(n_markers ~ mu, scales = "fixed") + scale_y_sqrt()
} else {
    plt = plt + facet_grid(n_markers ~ ., scales = "fixed") + scale_y_sqrt()
}
ggsave(paste0(pdir, type, "_lm.pdf"), plt, width = 14, height = 10)

# Scatter plots
pdf(paste0(pdir, type, "_scatter.pdf"), width = 15, height = 7)
for (i in 1:length(L)) {
    E = L[[i]]$err
    plt = ggplot(data = E, mapping = aes(x = truth, y = est, color = method, shape = method)) + 
        geom_abline(intercept = 0, slope = 1, lty = 2) + geom_point() + facet_wrap(~method) + 
        xlim(0, 1) + ylim(0, 1)
    plt = plt + ggtitle(paste(names(grd[i, ]), grd[i, ], sep = "=", collapse = ", "))
    plt = plt + geom_smooth(method = "lm", color = "black") + scale_shape_manual(values = 0:14)
    print(plt)
}
dev.off()

sign_plt = function(S, N, M = 1) {
    smry = lapply(L, "[[", "err")
    for (i in 1:length(smry)) {
        smry[[i]]$i <- i
    }
    smry = do.call(rbind, smry)
    grd$i = 1:nrow(grd)
    d = merge(smry, grd)
    d = data.table(d)
    d = d[mu == M & sigma == S & n_markers %in% N, ]
    d$n_markers = paste0("N. Marker = ", d$n_markers)
    d$n_markers = ordered(d$n_markers, levels = unique(d$n_marker))
    d$mu = paste0("Seq. Depth = ", d$mu)
    d$mu = ordered(d$mu, levels = unique(d$mu))
    d$sigma = paste0("Var. Mult. = ", d$sigma)
    d$sigma = ordered(d$sigma, levels = unique(d$sigma))
    # levels(d$method)[2] <- 'Our Approach'
    plt = ggplot(data = d, mapping = aes(x = truth, y = est, group = interaction(sigma, 
        method, n_markers), color = method, shape = method)) + geom_point() + facet_grid(n_markers + 
        sigma ~ method) + scale_shape_manual(values = 0:15)
    plt = plt + theme_bw() + geom_abline(intercept = 0, slope = 1, lty = 2) + xlim(0, 
        1) + ylim(0, 1)
    plt = plt + labs(x = "Truth", y = "Estimate") + theme(axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 5), )
    return(plt)
}

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}


library("gridExtra")
pl = list(sign_plt(.5, 10),sign_plt(.5, 100),sign_plt(.5, 1000))
leg = g_legend(pl[[1]])
pl = lapply(pl, function(p) p + guides(color = FALSE, shape = FALSE)+coord_fixed())
pg = plot_grid(plotlist = pl, nrow = length(pl), labels = "AUTO")
 plt = grid.arrange(pg,leg,ncol=1,heights=c(5,1))
ggsave(paste0(pdir, type, "_4scatter_paper.pdf"), pg, width = 8, height = 6)

library("gridExtra")
pl = list(sign_plt(.2, 1),sign_plt(1, 1),sign_plt(5, 1),
          sign_plt(.2, 10),sign_plt(1, 10),sign_plt(5, 10),
          sign_plt(.2, 100),sign_plt(1, 100),sign_plt(5, 100),
          sign_plt(.2, 1000),sign_plt(1, 1000),sign_plt(5, 1000)
          )
# pl = list(sign_plt(1,10),sign_plt(1,100))
leg = g_legend(pl[[1]])
pl = lapply(pl, function(p) p + guides(color = FALSE, shape = FALSE)+coord_fixed())
pg = plot_grid(plotlist = pl, nrow = length(pl), labels = "AUTO")
# plt = grid.arrange(pg,leg,ncol=1,heights=c(5,1))
ggsave(paste0(pdir, type, "_4scatter_paper_wide.pdf"), pg, width = 9, height = 22)
