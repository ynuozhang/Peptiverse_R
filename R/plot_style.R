# R/plot_style.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(showtext)
  library(sysfonts)
})

font_add_google("Ubuntu", "ubuntu")
showtext_auto()
theme_set(
  theme_bw(base_size = 15, base_family = "ubuntu") +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 20, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      strip.background = element_rect(fill = "#EDEDFD", color = NA),
      strip.text = element_text(face = "bold", color = "#013FB1", size = 20),
      legend.position = "bottom",
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),   # larger font
      legend.background = element_rect(
        fill = "#EDEDFD",
        color = "#013FB1",
        linewidth = 0.5
      ),
      legend.margin = margin(6, 8, 6, 8),
      legend.box = "vertical",
      legend.key.width = unit(1.6, "lines")
    )
)

PAL_MODEL <- c(
  "XGB" = "#01e3d2",
  "SVM" = "#cc5bf4",
  "ENET" = "#0098f1",
  "MLP" = "#f086ff",
  "CNN" = "#8662be",
  "Transformer" = "#f25e96"
)
