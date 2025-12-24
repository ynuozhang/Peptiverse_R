#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(patchwork)
})
source("R/plot_style.R")
# -----------------------------
# Root + output
# -----------------------------
ROOT <- "~/Documents/Peptiverse_R/data_distribution"
OUT_DIR <- file.path(ROOT, "figures", "training_val_distributions")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAL_SPLIT <- c("train" = "#de9cec", "val" = "#48d1ec")

normalize_split <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("\\s+", "", x)
  dplyr::case_when(
    x %in% c("train", "training", "tr") ~ "train",
    x %in% c("val", "valid", "validation", "va", "dev") ~ "val",
    TRUE ~ x
  )
}

read_meta <- function(csv_path, split_col) {
  df <- read_csv(csv_path, show_col_types = FALSE)
  stopifnot(split_col %in% names(df))
  df %>%
    mutate(split = normalize_split(.data[[split_col]])) %>%
    filter(split %in% c("train", "val"))
}

plot_binary_bars <- function(csv_path, prop_name, split_col, label_col) {
  df <- read_meta(csv_path, split_col)
  stopifnot(label_col %in% names(df))
  
  df <- df %>%
    mutate(y = suppressWarnings(as.integer(.data[[label_col]]))) %>%
    filter(y %in% c(0, 1)) %>%
    mutate(
      y = factor(y, levels = c(0, 1), labels = c("0", "1")),
      split = factor(split, levels = c("train", "val"))
    )
  
  n_train <- sum(df$split == "train")
  n_val   <- sum(df$split == "val")
  
  counts <- df %>% count(split, y, name = "n")
  
  ggplot(counts, aes(x = y, y = n, fill = split)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    scale_fill_manual(values = PAL_SPLIT, breaks = c("train", "val")) +
    labs(x = NULL, y = NULL) +
    ggtitle(prop_name) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("train N=", n_train, "\nval N=", n_val),
      hjust = 1.08, vjust = 1.15,
      size = 5.0, family = "ubuntu"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "bottom"   # <-- Grid A legend lives here
    )
}

plot_continuous_overlay <- function(csv_path, prop_name, split_col, value_col,
                                    thresholds = numeric(0), x_lab = NULL,
                                    left_pad = 1L, right_pad = 0L) {
  df <- read_meta(csv_path, split_col)
  stopifnot(value_col %in% names(df))
  
  df <- df %>%
    mutate(v = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    filter(!is.na(v)) %>%
    mutate(split = factor(split, levels = c("train", "val")))
  
  n_train <- sum(df$split == "train")
  n_val   <- sum(df$split == "val")
  
  # raw data range
  xmin_raw <- min(df$v, na.rm = TRUE)
  xmax_raw <- max(df$v, na.rm = TRUE)
  
  # integer-rounded endpoints
  xmin0 <- floor(xmin_raw)
  xmax0 <- ceiling(xmax_raw)
  
  # 5 integer ticks => 4 equal integer segments
  span <- xmax0 - xmin0
  step_int <- max(1L, ceiling(span / 4L))
  
  xmin_tick <- xmin0
  xmax_tick <- xmin_tick + 4L * step_int
  xbreaks <- xmin_tick + step_int * (0:4)
  
  # extend view a bit left (and optionally right) WITHOUT dropping data
  x_view_min <- xmin_tick - as.integer(left_pad)
  x_view_max <- xmax_tick + as.integer(right_pad)
  
  p <- ggplot(df, aes(x = v, fill = split)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.55) +
    scale_fill_manual(values = PAL_SPLIT, breaks = c("train", "val"), guide = "none") +
    scale_x_continuous(
      breaks = xbreaks,
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(xlim = c(x_view_min, x_view_max)) +   # <-- zoom only, keep all data
    labs(x = x_lab, y = NULL) +
    theme(axis.title.x = element_text(size = 12)) +
    ggtitle(prop_name) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("train N=", n_train, "\nval N=", n_val),
      hjust = 1.08, vjust = 1.15,
      size = 5.0, family = "ubuntu"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "none"
    )
  
  if (length(thresholds) > 0) {
    for (thr in thresholds) {
      p <- p + geom_vline(xintercept = thr, linetype = "dashed", linewidth = 0.8)
    }
  }
  p
}



# ----- explicit specs -----
continuous_jobs <- tribble(
  ~prop,                       ~csv, ~split_col, ~value_col, ~thresholds,     ~x_lab,                     
  "Permeability (PAMPA)",      file.path(ROOT,"permeability_pampa","pampa_meta_with_split.csv"),
                               "split","PAMPA", list(c(-6)), "Log Pexp",
                               

  "Permeability (CACO2)",      file.path(ROOT,"permeability_caco2","caco2_meta_with_split.csv"),
                               "split","Caco2", list(c(-6)), "Log Pexp",
                               

  "Binding affinity (WT)",     file.path(ROOT,"binding_affinity","binding_affinity_wt_meta_with_split.csv"),
                               "split","affinity", list(c(7,9)), "Affinity Score",
                               

  "Binding affinity (SMILES)", file.path(ROOT,"binding_affinity","binding_affinity_smiles_meta_with_split.csv"),
                               "split","affinity", list(c(7,9)), "Affinity Score",
                              
)

binary_jobs <- tribble(
  ~prop,         ~csv,                                              ~split_col, ~label_col,
  "Solubility",  file.path(ROOT, "solubility", "sol_meta_with_split.csv"), "split", "label",
  "Hemolysis",   file.path(ROOT, "hemolysis",  "hemo_meta_with_split.csv"),"split", "label",
  "Nonfouling", file.path(ROOT, "nf",         "nf_meta_with_split.csv"),  "split", "label",
  "Toxicity", file.path(ROOT, "toxicity",         "tox_meta_with_split.csv"),  "split", "Label",
  "Permeability (Penetrance)", file.path(ROOT, "permeability_penetrance",         "data_full.csv"),  "split", "label"
)

plots_bin <- purrr::pmap(
  binary_jobs,
  \(prop, csv, split_col, label_col) plot_binary_bars(csv, prop, split_col, label_col)
)

plots_cont <- purrr::pmap(
  continuous_jobs,
  \(prop, csv, split_col, value_col, thresholds, x_lab)
  plot_continuous_overlay(
    csv, prop, split_col, value_col,
    thresholds = thresholds,
    x_lab = x_lab
  )
)

grid_A <- wrap_plots(plots_bin,  ncol = 3)   # legend comes from here
grid_B <- wrap_plots(plots_cont, ncol = 2)   # legends suppressed

final_plot <- (grid_A / grid_B) +
  plot_layout(guides = "collect") &          # collects only from plots that have legends
  theme(legend.position = "bottom")

print(final_plot)

