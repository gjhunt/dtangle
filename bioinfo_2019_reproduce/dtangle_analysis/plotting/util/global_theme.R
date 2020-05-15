library("ggthemes")

global_size <- 11

global_theme <- theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = global_size), 
    axis.text = element_text(size = global_size), axis.title = element_text(size = global_size, 
        face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 0.9 * global_size), 
    legend.text = element_text(size = 0.75 * global_size), legend.key.width = unit(0.3, 
        "line"), legend.key.height = unit(0.75, "line"), legend.title = element_text(size = global_size), 
    legend.position = "bottom", plot.margin = unit(c(0, 0, 1, 0), "lines"), panel.spacing = unit(2, 
        "lines"), strip.text = element_text(size = global_size), axis.line = element_line(size = 0.7, 
        color = "black"), axis.text.x = element_text(angle = 30, hjust = 1), legend.box = "vertical")
