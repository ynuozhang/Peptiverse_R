# R/plot_style.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#EDEDFD", color = NA),
      strip.text = element_text(face = "bold", color = "#013FB1"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
)

PAL_MODEL <- c(
  xgb = "#01e3d2",
  svm_gpu = "#cc5bf4",
  enet_gpu = "#0098f1",
  mlp = "#f086ff",
  cnn = "#8662be",
  transformer = "#f25e96"
)
