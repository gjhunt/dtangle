source("util/util.R", chdir = TRUE)

q_res <- readRDS("../comparison/Analysis/t_sens/t_sens.rds")
E <- all_errors(q_res)
T <- times(q_res)

T$Time <- T$Seconds
T$Time <- T$Time/60

Tm <- aggregate(Time ~ Deconv.Method + Marker.Method + Quantile, T, mean)
Tm$TimeMax <- aggregate(Time ~ Deconv.Method + Marker.Method + Quantile, T, max)$Time
Tm$TimeMin <- aggregate(Time ~ Deconv.Method + Marker.Method + Quantile, T, min)$Time

p <- ggplot(Tm, aes(y = Time, x = Quantile, color = Deconv.Method, shape = Deconv.Method)) + 
    geom_point(size = 1.5, stroke = 1) + geom_line(size = 1) + global_theme
p <- p + geom_ribbon(aes(ymin = TimeMin, ymax = TimeMax, fill = Deconv.Method), alpha = 0.05, 
    lty = 2, lwd = 0.35)
p <- p + scale_shape_manual(values = 0:14)
p <- p + labs(title = "Time Efficiency Across Methods", x = "Marker Quantile", y = "Mean of log10(Minutes)")
p <- p + scale_y_log10()
p <- p + scale_color_brewer(palette = "Set1")
p <- p + scale_fill_brewer(palette = "Set1")


save_plots(p, dir = "./time/", name = "time_plot", width = 5, height = 4)
