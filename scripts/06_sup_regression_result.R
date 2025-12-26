ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(scales)
  library(dplyr)
})

RAW_DIR <- "raw"
OUT_DIR <- "figures/supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


if (!exists("parse_meta_from_path")) {
  stop("parse_meta_from_path() not found. Please source the script where you defined it, then rerun.")
}

read_csv_fast <- function(path) data.table::fread(path, data.table = FALSE, showProgress = FALSE)


norm_task <- function(x) str_replace_all(tolower(x), "[^a-z0-9]+", "")

task_to_prop <- function(task_norm) {
  dplyr::case_when(
    str_detect(task_norm, "halflife") ~ "Half-life",
    str_detect(task_norm, "bindingaffinity") ~ "Binding_affinity",
    str_detect(task_norm, "permeabilitypampa") ~ "Permeability_PAMPA",
    str_detect(task_norm, "permeabilitycaco2") ~ "Permeability_CACO2",
    TRUE ~ NA_character_
  )
}

# -------------------------
# Regression metrics on a val_predictions.csv (y_true/y_pred)
# We compute: pearson r, spearman rho, r2(lm), rmse, mae
# -------------------------
compute_reg_metrics_from_val <- function(df, quiet = TRUE) {
  nm <- names(df)
  
  pick_first <- function(cands) {
    hit <- intersect(cands, nm)
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
  pred_col  <- pick_first(c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
  
  if (is.na(truth_col) || is.na(pred_col)) {
    stop("Regression val file missing y_true/y_pred cols. Found: ", paste(nm, collapse = ", "))
  }
  
  yt <- suppressWarnings(as.numeric(df[[truth_col]]))
  yp <- suppressWarnings(as.numeric(df[[pred_col]]))
  ok <- is.finite(yt) & is.finite(yp)
  yt <- yt[ok]; yp <- yp[ok]
  
  if (!quiet) message(sprintf("[reg-metrics] truth=%s pred=%s n=%d", truth_col, pred_col, length(yt)))
  
  pear <- suppressWarnings(cor(yt, yp, method = "pearson"))
  rho  <- suppressWarnings(cor(yt, yp, method = "spearman"))
  
  lmfit <- tryCatch(lm(yt ~ yp), error = function(e) NULL)
  r2 <- if (!is.null(lmfit)) summary(lmfit)$r.squared else NA_real_
  
  rmse <- sqrt(mean((yt - yp)^2))
  mae  <- mean(abs(yt - yp))
  
  tibble(
    n = length(yt),
    pearson_r = as.numeric(pear),
    spearman_rho = as.numeric(rho),
    r2_lm = as.numeric(r2),
    rmse = as.numeric(rmse),
    mae  = as.numeric(mae)
  )
}

# -------------------------
# Plot style: y_pred vs y_true with lm + stats text
# -------------------------
make_scatter <- function(df, stats_df, facet_var = "panel") {
  ggplot(df, aes(x = y_true, y = y_pred)) +
    geom_point(
      shape = 16,
      color = "#dba5d6",
      alpha = 0.70,
      size = 2.0
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "#5dd8df",
      fill  = "#5dd8df",
      alpha = 0.18,
      linewidth = 0.9
    ) +
    facet_wrap(as.formula(paste("~", facet_var)), scales = "free") +
    labs(x = "y_true", y = "y_pred") +
    theme(legend.position = "none") +
    geom_text(
      data = stats_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.05,
      size = 4.8
    )
}

# =========================================================
# 1) BEST MODEL PER TASK×REPR from raw/**/val_predictions.csv
#    (excluding binding affinity; we special-case it)
# =========================================================
val_files <- list.files(RAW_DIR, pattern = "val_predictions\\.csv$", recursive = TRUE, full.names = TRUE)

val_meta <- purrr::map_dfr(val_files, function(f) {
  meta <- parse_meta_from_path(f, RAW_DIR)  # must return task, repr, model, run
  meta %>% mutate(path = f)
}) %>%
  mutate(
    repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr)),
    task_norm = norm_task(task),
    prop = task_to_prop(task_norm)
  ) %>%
  filter(!is.na(prop)) %>%
  # drop binding affinity here; we add it from your explicit file below
  filter(prop != "Binding_affinity")

val_metrics <- val_meta %>%
  mutate(
    m = purrr::map(path, ~ compute_reg_metrics_from_val(read_csv_fast(.x), quiet = TRUE))
  ) %>%
  unnest(m)

# Pick best run per prop×repr by max spearman rho (ties: higher r2, then lower rmse)
best_runs <- val_metrics %>%
  group_by(prop, repr) %>%
  arrange(desc(spearman_rho), desc(r2_lm), rmse) %>%
  slice(1) %>%
  ungroup()

message("[Best regressors picked]")
print(best_runs %>% select(prop, repr, task, model, run, spearman_rho, r2_lm, rmse, mae, path))

best_scatter_df <- best_runs %>%
  mutate(
    df = purrr::map(path, function(p) {
      d <- read_csv_fast(p)
      nm <- names(d)
      
      pick_first <- function(cands) { hit <- intersect(cands, nm); if (length(hit)==0) NA_character_ else hit[1] }
      truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
      pred_col  <- pick_first(c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
      
      yt <- suppressWarnings(as.numeric(d[[truth_col]]))
      yp <- suppressWarnings(as.numeric(d[[pred_col]]))
      ok <- is.finite(yt) & is.finite(yp)
      
      tibble(
        y_true = yt[ok],
        y_pred = yp[ok]
      )
    }),
    repr_lab = if_else(repr == "wt", "WT", if_else(repr == "smiles", "SMILES", toupper(repr))),
    panel = sprintf("%s (%s)\nBest: %s", prop, repr_lab, toupper(model))
  ) %>%
  select(prop, repr, model, run, panel, df) %>%
  unnest(df)

best_stats <- best_runs %>%
  mutate(
    repr_lab = if_else(repr == "wt", "WT", if_else(repr == "smiles", "SMILES", toupper(repr))),
    panel = sprintf("%s (%s)\nBest: %s", prop, repr_lab, toupper(model)),
    label = sprintf(
      "N=%d\nPearson r=%.3f\nSpearman ρ=%.3f\nLinear regression R²=%.3f",
      n, pearson_r, spearman_rho, r2_lm
    )
  ) %>%
  select(panel, label)

# =========================================================
# 2) Binding affinity: force WT_SMILES_UNPOOLED val file
#    path: binding_affinity/wt_smiles_unpooled/val_smiles_unpooled.csv
# =========================================================
BA_VAL_SMI <- file.path("raw/binding_affinity", "wt_smiles_unpooled", "val_smiles_unpooled.csv")
BA_VAL_WT  <- file.path("raw/binding_affinity", "wt_wt_unpooled",     "val_wt_unpooled.csv")

if (file.exists(BA_VAL)) {
  # =========================================================
  # 2) Binding affinity: force WT_UNPOOLED + SMILES_UNPOOLED val files
  # =========================================================
  read_ba_val <- function(path, panel_label) {
    dba <- readr::read_csv(path, show_col_types = FALSE)
    
    nm <- names(dba)
    pick_first <- function(cands) { hit <- intersect(cands, nm); if (length(hit)==0) NA_character_ else hit[1] }
    truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
    pred_col  <- pick_first(c("y_pred_reg","pred","y_hat","prediction","pred_mean","estimate"))
    
    if (is.na(truth_col) || is.na(pred_col)) {
      stop("Binding affinity val file missing y_true/y_pred cols: ", path)
    }
    
    yt <- suppressWarnings(as.numeric(dba[[truth_col]]))
    yp <- suppressWarnings(as.numeric(dba[[pred_col]]))
    ok <- is.finite(yt) & is.finite(yp)
    
    df <- tibble(
      y_true = yt[ok],
      y_pred = yp[ok],
      panel  = panel_label
    )
    
    # stats
    pear <- suppressWarnings(cor(df$y_true, df$y_pred, method = "pearson"))
    rho  <- suppressWarnings(cor(df$y_true, df$y_pred, method = "spearman"))
    r2   <- tryCatch(summary(lm(y_pred ~ y_true, data = df))$r.squared, error = function(e) NA_real_)
    
    st <- tibble(
      panel = panel_label,
      label = sprintf(
        "N=%d\nPearson r=%.3f\nSpearman ρ=%.3f\nLinear regression R²=%.3f",
        nrow(df), pear, rho, r2
      )
    )
    
    list(df = df, st = st)
  }
  
  # SMILES unpooled
  if (file.exists(BA_VAL_SMI)) {
    res <- read_ba_val(BA_VAL_SMI, "Binding_affinity (SMILES)\nBest: TRANSFORMER")
    best_scatter_df <- bind_rows(best_scatter_df, res$df)
    best_stats      <- bind_rows(best_stats,      res$st)
    message("Included binding affinity SMILES panel from: ", BA_VAL_SMI)
  } else {
    message("Binding affinity SMILES val file not found (skipping): ", BA_VAL_SMI)
  }
  
  # WT unpooled
  if (file.exists(BA_VAL_WT)) {
    res <- read_ba_val(BA_VAL_WT, "Binding_affinity (WT)\nBest: TRANSFORMER")
    best_scatter_df <- bind_rows(best_scatter_df, res$df)
    best_stats      <- bind_rows(best_stats,      res$st)
    message("Included binding affinity WT panel from: ", BA_VAL_WT)
  } else {
    message("Binding affinity WT val file not found (skipping): ", BA_VAL_WT)
  }
} else {
  message("Binding affinity val file not found (skipping): ", BA_VAL_WT)
}

# =========================================================
# 3) Plot + save
# =========================================================
PANEL_ORDER <- c(
  "Permeability_PAMPA (SMILES)\nBest: XGB",
  "Permeability_CACO2 (SMILES)\nBest: XGB",
  "Half-life (SMILES)\nBest: XGB",
  "Binding_affinity (SMILES)\nBest: TRANSFORMER"
)

best_scatter_df <- best_scatter_df %>%
  mutate(panel = factor(panel, levels = unique(c(PANEL_ORDER, sort(unique(as.character(panel)))))))

p <- make_scatter(best_scatter_df, best_stats, facet_var = "panel") +
  theme(
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )

print(p)

out_png <- file.path(OUT_DIR, "supp_best_regressors_val_scatter.png")
out_pdf <- file.path(OUT_DIR, "supp_best_regressors_val_scatter.pdf")

ggsave(out_png, p, width = 12, height = 3.8 * max(1, ceiling(nlevels(best_scatter_df$panel) / 2)), dpi = 400)
ggsave(out_pdf, p, width = 12, height = 3.8 * max(1, ceiling(nlevels(best_scatter_df$panel) / 2)))

message("Saved:\n  ", out_png, "\n  ", out_pdf)
