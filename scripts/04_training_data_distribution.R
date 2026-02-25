suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(patchwork)
})
source("R/plot_style.R")

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

read_meta <- function(csv_path, split_col = NA_character_) {
  df <- read_csv(csv_path, show_col_types = FALSE)
  
  if (is.null(split_col) || is.na(split_col) || !(split_col %in% names(df))) {
    return(df %>% mutate(split = factor("all", levels = "all")))
  }
  
  df %>%
    mutate(split = normalize_split(.data[[split_col]])) %>%
    filter(split %in% c("train", "val")) %>%
    mutate(split = factor(split, levels = c("train", "val")))
}

plot_binary_bars <- function(csv_path, prop_name, split_col, label_col) {
  df <- read_meta(csv_path, split_col)
  stopifnot(label_col %in% names(df))
  
  df <- df %>%
    mutate(y = suppressWarnings(as.integer(.data[[label_col]]))) %>%
    filter(y %in% c(0, 1)) %>%
    mutate(y = factor(y, levels = c(0, 1), labels = c("0", "1")))
  
  if (levels(df$split)[1] == "all") {
    n_all <- nrow(df)
    counts <- df %>% count(y, name = "n")
    
    return(
      ggplot(counts, aes(x = y, y = n)) +
        geom_col(width = 0.65) +
        labs(x = NULL, y = NULL) +
        ggtitle(prop_name) +
        annotate(
          "text", x = Inf, y = Inf,
          label = paste0("N=", n_all),
          hjust = 1.08, vjust = 1.15,
          size = 5.0, family = "ubuntu"
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 18),
          legend.position = "none"
        )
    )
  }
  
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
      legend.position = "bottom"
    )
}

plot_continuous_overlay <- function(csv_path, prop_name, split_col, value_col,
                                    thresholds = numeric(0), x_lab = NULL,
                                    left_pad = 1L, right_pad = 0L) {
  df <- read_meta(csv_path, split_col)
  stopifnot(value_col %in% names(df))
  
  df <- df %>%
    mutate(v = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    filter(!is.na(v))
  
  # counts
  if (levels(df$split)[1] == "all") {
    n_label <- paste0("N=", nrow(df))
  } else {
    n_train <- sum(df$split == "train")
    n_val   <- sum(df$split == "val")
    n_label <- paste0("train N=", n_train, "\nval N=", n_val)
  }
  
  xmin_raw <- min(df$v, na.rm = TRUE)
  xmax_raw <- max(df$v, na.rm = TRUE)
  xmin0 <- floor(xmin_raw)
  xmax0 <- ceiling(xmax_raw)
  span <- xmax0 - xmin0
  step_int <- max(1L, ceiling(span / 4L))
  xmin_tick <- xmin0
  xmax_tick <- xmin_tick + 4L * step_int
  xbreaks <- xmin_tick + step_int * (0:4)
  x_view_min <- xmin_tick - as.integer(left_pad)
  x_view_max <- xmax_tick + as.integer(right_pad)
  
  # plot
  p <- ggplot(df, aes(x = v)) +
    geom_histogram(bins = 50, alpha = 0.55, fill = PAL_SPLIT[["train"]]) +
    scale_x_continuous(
      breaks = xbreaks,
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(xlim = c(x_view_min, x_view_max)) +
    labs(x = x_lab, y = NULL) +
    theme(axis.title.x = element_text(size = 12)) +
    ggtitle(prop_name) +
    annotate(
      "text", x = Inf, y = Inf,
      label = n_label,
      hjust = 1.08, vjust = 1.15,
      size = 5.0, family = "ubuntu"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "none"
    )
  
  if (!(levels(df$split)[1] == "all")) {
    p <- ggplot(df, aes(x = v, fill = split)) +
      geom_histogram(bins = 50, position = "identity", alpha = 0.55) +
      scale_fill_manual(values = PAL_SPLIT, breaks = c("train", "val"), guide = "none") +
      scale_x_continuous(breaks = xbreaks, expand = expansion(mult = c(0, 0.02))) +
      coord_cartesian(xlim = c(x_view_min, x_view_max)) +
      labs(x = x_lab, y = NULL) +
      theme(axis.title.x = element_text(size = 12)) +
      ggtitle(prop_name) +
      annotate(
        "text", x = Inf, y = Inf,
        label = n_label,
        hjust = 1.08, vjust = 1.15,
        size = 5.0, family = "ubuntu"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 18),
        legend.position = "none"
      )
  }
  
  if (length(thresholds) > 0) {
    for (thr in thresholds) {
      p <- p + geom_vline(xintercept = thr, linetype = "dashed", linewidth = 0.8)
    }
  }
  p
}


continuous_jobs <- tribble(
  ~prop,                       ~csv, ~split_col, ~value_col,        ~thresholds,  ~x_lab,
  
  "Permeability (PAMPA)",      file.path(ROOT,"permeability_pampa","pampa_meta_with_split.csv"),
  "split","PAMPA",  list(c(-6)),       "Log Pexp",
  
  "Permeability (CACO2)",      file.path(ROOT,"permeability_caco2","caco2_meta_with_split.csv"),
  "split","Caco2",  list(c(-6)),       "Log Pexp",
  
  "Binding affinity (Amino Acids)",     file.path(ROOT,"binding_affinity","binding_affinity_wt_meta_with_split.csv"),
  "split","affinity", list(c(7,9)),    "Affinity Score",
  
  "Binding affinity (SMILES)", file.path(ROOT,"binding_affinity","binding_affinity_smiles_meta_with_split.csv"),
  "split","affinity", list(c(7,9)),    "Affinity Score",
  "Half-life (Amino Acids)",                 file.path(ROOT,"half_life","wt_halflife_merged_dedup.csv"),
  NA,    "half_life_hours", list(numeric(0)), "Half-life (hours)",
  "Half-life (SMILES)",                 file.path(ROOT,"half_life","halflife_merged_dedup.csv"),
  NA,    "half_life_hours", list(numeric(0)), "Half-life (hours)",
  "Stability (Amino Acids)", file.path(ROOT,"stability","stability_meta_with_split.csv"),
  "split","label", list(numeric(0)),    "Score",

)


binary_jobs <- tribble(
  ~prop,         ~csv,                                              ~split_col, ~label_col,
  "Solubility",  file.path(ROOT, "solubility", "sol_meta_with_split.csv"), "split", "label",
  "Hemolysis",   file.path(ROOT, "hemolysis",  "hemo_meta_with_split.csv"),"split", "label",
  "Nonfouling", file.path(ROOT, "nf",         "nf_meta_with_split.csv"),  "split", "label",
  "Toxicity", file.path(ROOT, "toxicity",         "tox_meta_with_split.csv"),  "split", "Label",
  "Permeability (Penetrance)", file.path(ROOT, "permeability_penetrance",         "data_full.csv"),  "split", "label",
)

plots_bin <- purrr::pmap(
  binary_jobs,
  \(prop, csv, split_col, label_col) plot_binary_bars(csv, prop, split_col, label_col)
)

plots_cont <- purrr::pmap(
  continuous_jobs,
  \(prop, csv, split_col, value_col, thresholds, x_lab) {
    plot_continuous_overlay(
      csv, prop, split_col, value_col,
      thresholds = thresholds,
      x_lab = x_lab
    )
  }
)


grid_A <- wrap_plots(plots_bin,  ncol = 3) 
grid_B <- wrap_plots(plots_cont, ncol = 2)

final_plot <- (grid_A / grid_B) +
  plot_layout(guides = "collect") &  
  theme(legend.position = "bottom")

print(final_plot)
out_png <- file.path(OUT_DIR, "distribution.png")
ggsave(out_png, final_plot, width = 10, height = 8, dpi = 500) 
message("Saved:\n  ", out_png)

# -----------------------------
# PRINT distributions 
# -----------------------------

read_meta_keep_all <- function(csv_path, split_col = NA_character_) {
  df <- read_csv(csv_path, show_col_types = FALSE)
  
  if (is.null(split_col) || is.na(split_col) || !(split_col %in% names(df))) {
    return(df %>% mutate(split = factor("all", levels = "all")))
  }
  
  df %>%
    mutate(split = normalize_split(.data[[split_col]])) %>%
    mutate(split = ifelse(split %in% c("train", "val"), split, NA_character_)) %>%
    filter(!is.na(split)) %>%
    mutate(split = factor(split, levels = c("train", "val")))
}

print_binary_distribution <- function(prop, csv, split_col, label_col) {
  cat("\n==============================\n")
  cat("BINARY:", prop, "\n")
  cat("File:", csv, "\n")
  
  df <- read_meta_keep_all(csv, split_col)
  if (!(label_col %in% names(df))) {
    cat("  [SKIP] label_col not found:", label_col, "\n")
    return(invisible(NULL))
  }
  
  df <- df %>%
    mutate(y = suppressWarnings(as.integer(.data[[label_col]]))) %>%
    filter(y %in% c(0, 1)) %>%
    mutate(y = factor(y, levels = c(0, 1), labels = c("0", "1")))
  
  # overall
  overall <- df %>% count(y, name = "n") %>% arrange(y)
  cat("Overall N =", nrow(df), "\n")
  cat("Counts: ", paste0(overall$y, "=", overall$n, collapse = ", "), "\n")
  
  # split (if present)
  if (!(levels(df$split)[1] == "all")) {
    by_split <- df %>% count(split, y, name = "n") %>% arrange(split, y)
    cat("By split:\n")
    for (sp in levels(df$split)) {
      d <- by_split %>% filter(split == sp)
      if (nrow(d) == 0) next
      cat("  ", sp, "N=", sum(d$n), " | ",
          paste0(d$y, "=", d$n, collapse = ", "),
          "\n", sep = "")
    }
  }
  invisible(NULL)
}

print_continuous_distribution <- function(prop, csv, split_col, value_col) {
  cat("\n==============================\n")
  cat("CONTINUOUS:", prop, "\n")
  cat("File:", csv, "\n")
  
  df <- read_meta_keep_all(csv, split_col)
  if (!(value_col %in% names(df))) {
    cat("  [SKIP] value_col not found:", value_col, "\n")
    return(invisible(NULL))
  }
  
  df <- df %>%
    mutate(v = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    filter(!is.na(v))
  
  summarise_one <- function(d) {
    tibble(
      N = nrow(d),
      min = min(d$v),
      median = median(d$v),
      mean = mean(d$v),
      max = max(d$v)
    )
  }
  
  overall <- summarise_one(df)
  cat("Overall:\n")
  cat("  N=", overall$N,
      " | min=", signif(overall$min, 4),
      " median=", signif(overall$median, 4),
      " mean=", signif(overall$mean, 4),
      " max=", signif(overall$max, 4),
      "\n", sep = "")
  
  if (!(levels(df$split)[1] == "all")) {
    cat("By split:\n")
    for (sp in levels(df$split)) {
      d <- df %>% filter(split == sp)
      if (nrow(d) == 0) next
      s <- summarise_one(d)
      cat("  ", sp, ": N=", s$N,
          " | min=", signif(s$min, 4),
          " median=", signif(s$median, 4),
          " mean=", signif(s$mean, 4),
          " max=", signif(s$max, 4),
          "\n", sep = "")
    }
  }
  
  invisible(NULL)
}

# ---- run prints ----
purrr::pwalk(
  binary_jobs,
  \(prop, csv, split_col, label_col) print_binary_distribution(prop, csv, split_col, label_col)
)

purrr::pwalk(
  continuous_jobs,
  \(prop, csv, split_col, value_col, thresholds, x_lab) print_continuous_distribution(prop, csv, split_col, value_col)
)

