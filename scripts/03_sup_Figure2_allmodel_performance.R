ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(yardstick)
  library(stringr)
  library(scales)
  library(dplyr)
})

RAW_DIR <- "raw"
OUT_DIR <- "figures/supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv_fast <- function(path) data.table::fread(path, data.table = FALSE, showProgress = FALSE)

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  dplyr::case_when(
    str_detect(x, "^xgb$|xgboost") ~ "XGB",
    str_detect(x, "svm")          ~ "SVM",
    str_detect(x, "enet|elastic|logreg") ~ "ENET",
    str_detect(x, "mlp")          ~ "MLP",
    str_detect(x, "cnn")          ~ "CNN",
    str_detect(x, "transformer|bert|esm") ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}
prop_spec <- tribble(
  ~prop,                      ~task_key,                 ~kind,          ~include_wt, ~include_smiles, ~placeholder_ok,
  "Hemolysis",                "hemolysis",               "classifier",    TRUE,        TRUE,            TRUE,
  "Nonfouling",               "nonfouling",              "classifier",    TRUE,        TRUE,            TRUE,
  "Halflife",                 "halflife",                "regression",    TRUE,        TRUE,            TRUE,
  "Toxicity",                 "toxicity",                "classifier",    TRUE,        TRUE,            TRUE,
  
  # DEFINITELY missing in SMILES -> do not include in SMILES at all
  "Solubility",               "solubility",              "classifier",    TRUE,        FALSE,           FALSE,
  
  # WT-only umbrella permeability (if exists)
  "Permeability",             "permeability",            "regression",    TRUE,        FALSE,           TRUE,
  
  # SMILES-only permeability assays
  "Permeability (PAMPA)",     "permeability_pampa",      "regression",    FALSE,       TRUE,            TRUE,
  "Permeability (CACO2)",     "permeability_caco2",      "regression",    FALSE,       TRUE,            TRUE,
  
  # Binding affinity regression in BOTH panels (placeholder until you add rho results)
  "Binding affinity",         "binding",                 "regression",    TRUE,        TRUE,            TRUE
)
task_matches <- function(task_col, key) {
  task_col == key | str_detect(task_col, fixed(key, ignore_case = TRUE))
}

metrics_plot_named <- prop_spec %>%
  select(prop, task_key, include_wt, include_smiles, placeholder_ok) %>%
  mutate(prop = factor(prop, levels = prop_spec$prop)) %>%
  # join by "task matches task_key"
  tidyr::crossing(task = unique(metrics_plot$task)) %>%
  filter(task_matches(task, task_key)) %>%
  distinct(task, prop, include_wt, include_smiles, placeholder_ok) %>%
  right_join(metrics_plot, by = "task") %>%
  mutate(
    # ensure prop exists even if task didn't match any key
    prop = if_else(is.na(prop), str_to_title(task), as.character(prop)),
    prop = factor(prop, levels = c(prop_spec$prop, setdiff(unique(as.character(prop)), prop_spec$prop)))
  )

# Drop disallowed panels (e.g., Solubility SMILES)
metrics_plot_named <- metrics_plot_named %>%
  filter(
    (repr == "WT"     & (is.na(include_wt)     | include_wt)) |
      (repr == "SMILES" & (is.na(include_smiles) | include_smiles))
  )
# ---- use your parse_meta_from_path() here (same as you already have) ----
# parse_meta_from_path <- function(path, raw_dir) { ... }

# --- metrics from a single val_predictions.csv ---
compute_metrics_from_val <- function(df) {
  d <- df %>%
    transmute(
      truth    = factor(as.integer(y_true), levels = c(0,1)),
      estimate = factor(as.integer(y_pred), levels = c(0,1)),
      .pred_1  = as.numeric(y_prob)
    )
  
  has_both <- dplyr::n_distinct(d$truth) == 2
  
  f1 <- tryCatch(f_meas(d, truth, estimate, event_level = "second")$.estimate, error = \(e) NA_real_)
  acc <- tryCatch(accuracy(d, truth, estimate)$.estimate, error = \(e) NA_real_)
  mccv <- tryCatch(mcc(d, truth, estimate, event_level = "second")$.estimate, error = \(e) NA_real_)
  
  aucv <- if (has_both) {
    tryCatch(roc_auc(d, truth, .pred_1, event_level = "second")$.estimate, error = \(e) NA_real_)
  } else NA_real_
  
  tibble(
    metric = c("F1", "AUC", "MCC", "Accuracy"),
    value  = c(f1, aucv, mccv, acc)
  )
}

# -------------------------
# Discover val files
# -------------------------
pred_files <- list.files(RAW_DIR, pattern = "val_predictions\\.csv$", recursive = TRUE, full.names = TRUE)

preds_meta <- purrr::map_dfr(pred_files, function(f) {
  meta <- parse_meta_from_path(f, RAW_DIR)
  meta %>% mutate(path = f, split = "val")
})

# -------------------------
# Compute metrics per run (file)
# -------------------------
metrics_long <- preds_meta %>%
  mutate(
    repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr)),
    model_lbl = clean_model_label(model),
    tbl = purrr::map(path, ~ compute_metrics_from_val(read_csv_fast(.x)))
  ) %>%
  unnest(tbl) %>%
  select(task, repr, model_lbl, run, metric, value)

# -------------------------
# Aggregate across runs -> mean ± SE
# -------------------------
metrics_sum <- metrics_long %>%
  group_by(task, repr, model_lbl, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = if_else(n > 0, sd / sqrt(n), NA_real_),
    .groups = "drop"
  )

write_csv(metrics_sum, file.path(OUT_DIR, "supp_per_task_metrics_summary.csv"))

# -------------------------
# Placeholders so each metric shows all 6 models
# -------------------------
ALL_METRICS <- c("F1","AUC","MCC","Accuracy")
ALL_REPR    <- c("wt","smiles")
MODEL_ORDER <- c("XGB", "ENET", "SVM", "MLP", "CNN", "Transformer")

PAL_MODEL  <- PAL_MODEL[MODEL_ORDER]
ALL_MODELS <- MODEL_ORDER

metrics_plot <- metrics_sum %>%
  filter(repr %in% ALL_REPR) %>%
  mutate(
    model_lbl = clean_model_label(model_lbl),   # only if needed; otherwise remove
    model_lbl = factor(model_lbl, levels = ALL_MODELS),
    metric    = factor(metric, levels = ALL_METRICS),
    repr      = factor(repr, levels = ALL_REPR, labels = c("WT", "SMILES"))
  ) %>%
  tidyr::complete(
    task,
    repr,
    metric    = factor(ALL_METRICS, levels = ALL_METRICS),
    model_lbl = factor(ALL_MODELS,  levels = ALL_MODELS),
    fill = list(mean = NA_real_, se = NA_real_, n = 0)
  ) %>%
  mutate(
    available = !is.na(mean),
    mean_plot = if_else(available, mean, 0),
    se_plot   = if_else(available, se, 0),
    metric_id = as.numeric(metric),
    x_metric  = metric_id * 1.6          # <-- GAP between metrics (tune 1.4–1.9)
  )


# -------------------------
# Plot: one figure per task, panels = WT vs SMILES
# Bars packed tightly (no padding) and metrics adjacent
# -------------------------
# order tasks nicely (optional)
task_levels <- c("hemolysis","nonfouling","toxicity","solubility","halflife",
                 "permeability","permeability_pampa","permeability_caco2","binding")

metrics_plot2 <- metrics_plot %>%
  mutate(
    task = factor(task, levels = intersect(task_levels, unique(task)))
  )

p_big <- ggplot(metrics_plot_named, aes(x = x_metric, y = mean_plot, fill = model_lbl)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(
    data = metrics_plot_named %>% filter(available & !is.na(se_plot)),
    aes(ymin = pmax(mean_plot - se_plot, 0), ymax = pmin(mean_plot + se_plot, 1)),
    position = position_dodge(width = 0.9),
    width = 0.12, linewidth = 0.30
  ) +
  geom_text(
    data = metrics_plot_named %>% filter(!available),
    aes(label = "N/A", y = 0.03),
    position = position_dodge(width = 0.9),
    size = 3
  ) +
  facet_grid(prop ~ repr, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PAL_MODEL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_x_continuous(
    breaks = (1:length(ALL_METRICS)) * 1.6,
    labels = ALL_METRICS,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.01)) +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  )


print(p_big)

# Save ONE file
ggsave(file.path(OUT_DIR, "supp_all_tasks_per_metric_bars.png"),
       p_big, width = 12, height = 2.2 * nlevels(metrics_plot2$task), dpi = 500)
ggsave(file.path(OUT_DIR, "supp_all_tasks_per_metric_bars.pdf"),
       p_big, width = 12, height = 2.2 * nlevels(metrics_plot2$task))


