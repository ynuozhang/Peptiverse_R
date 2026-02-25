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

if (!exists("parse_meta_from_path")) {
  stop("parse_meta_from_path() not found. Please source the script where you defined it, then rerun.")
}

read_csv_fast <- function(path) data.table::fread(path, data.table = FALSE, showProgress = FALSE)

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  dplyr::case_when(
    str_detect(x, "^xgb$|xgboost")              ~ "XGB",
    str_detect(x, "svm")                       ~ "SVM",
    str_detect(x, "enet|elastic|logreg")        ~ "ENET",
    str_detect(x, "\\bmlp\\b")                 ~ "MLP",
    str_detect(x, "\\bcnn\\b")                 ~ "CNN",
    str_detect(x, "transformer|bert|esm")      ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}

prop_spec <- tribble(
  ~prop,        ~task_key,        ~include_wt, ~include_smiles, ~placeholder_ok,
  "Hemolysis",  "hemolysis",      TRUE,        TRUE,            TRUE,
  "Non-Fouling", "nonfouling",     TRUE,        TRUE,            TRUE,
  "Toxicity",   "toxicity",       FALSE,       TRUE,            TRUE,
  "Solubility", "solubility",     TRUE,        FALSE,           FALSE,
  "Permeability","permeabilitypenetrance",  TRUE,        FALSE,           TRUE
)

PROP_ORDER <- prop_spec$prop

compute_metrics_from_val <- function(df, threshold = 0.5, quiet = TRUE) {
  nm <- names(df)
  
  pick_first <- function(cands) {
    hit <- intersect(cands, nm)
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
  pred_col  <- pick_first(c("y_pred","pred","pred_label","estimate","yhat","prediction"))
  prob_col  <- pick_first(c("y_prob","prob","p","pred_prob","pred_probability","y_score","score",".pred_1",".pred1","prob_1","p1"))
  
  if (is.na(truth_col)) {
    stop("No truth column found. Expected one of: y_true, truth, label, y, target, gt. Found: ",
         paste(nm, collapse = ", "))
  }
  
  prob <- NULL
  if (!is.na(prob_col)) prob <- suppressWarnings(as.numeric(df[[prob_col]]))
  
  if (!is.na(pred_col)) {
    pred <- suppressWarnings(as.integer(df[[pred_col]]))
  } else if (!is.null(prob)) {
    pred <- as.integer(prob >= threshold)
  } else {
    stop("No predicted label or probability found. Need one of y_pred/pred/... or y_prob/prob/.... Found: ",
         paste(nm, collapse = ", "))
  }
  
  truth <- suppressWarnings(as.integer(df[[truth_col]]))
  
  if (!quiet) {
    msg <- sprintf("[metrics] truth=%s pred=%s prob=%s", truth_col,
                   ifelse(is.na(pred_col), "<thresholded>", pred_col),
                   ifelse(is.na(prob_col), "<none>", prob_col))
    message(msg)
  }
  
  d <- tibble(
    truth    = factor(truth, levels = c(0,1)),
    estimate = factor(pred,  levels = c(0,1)),
    .pred_1  = if (is.null(prob)) NA_real_ else prob
  )
  
  has_both <- dplyr::n_distinct(d$truth) == 2
  
  f1  <- tryCatch(yardstick::f_meas(d, truth, estimate, event_level = "second")$.estimate, error = \(e) NA_real_)
  acc <- tryCatch(yardstick::accuracy(d, truth, estimate)$.estimate, error = \(e) NA_real_)
  mccv<- tryCatch(yardstick::mcc(d, truth, estimate, event_level = "second")$.estimate, error = \(e) NA_real_)
  
  aucv <- if (has_both && !all(is.na(d$.pred_1))) {
    tryCatch(yardstick::roc_auc(d, truth, .pred_1, event_level = "second")$.estimate, error = \(e) NA_real_)
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
  meta <- parse_meta_from_path(f, RAW_DIR)   # MUST return at least: task, repr, model, run
  meta %>% mutate(path = f, split = "val")
})

# -------------------------
# Compute metrics per run (file)
# -------------------------
metrics_long <- preds_meta %>%
  mutate(
    repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr)),
    model_lbl = clean_model_label(model),
    tbl = purrr::map(path, ~ compute_metrics_from_val(read_csv_fast(.x), quiet = FALSE))
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
# Keep only classifier tasks we want (from prop_spec)
# -------------------------
KEEP_TASKS <- unique(prop_spec$task_key)
metrics_sum <- metrics_sum %>% filter(task %in% KEEP_TASKS, repr %in% c("wt", "smiles"))

# -------------------------
# Placeholders so every prop×repr×metric shows all models
# -------------------------
ALL_METRICS <- c("F1","AUC","MCC","Accuracy")
ALL_REPR    <- c("wt","smiles")
MODEL_ORDER <- c("XGB", "ENET", "SVM", "MLP", "CNN", "Transformer")

if (!exists("PAL_MODEL")) {
  PAL_MODEL <- setNames(RColorBrewer::brewer.pal(6, "Set2"), MODEL_ORDER)
} else {
  PAL_MODEL <- PAL_MODEL[MODEL_ORDER]
}

metrics_plot <- metrics_sum %>%
  mutate(
    metric    = factor(metric, levels = ALL_METRICS),
    model_lbl = factor(model_lbl, levels = MODEL_ORDER),
    repr      = factor(repr, levels = ALL_REPR, labels = c("Amino Acids","SMILES"))
  ) %>%
  tidyr::complete(
    task,
    repr,
    metric    = factor(ALL_METRICS, levels = ALL_METRICS),
    model_lbl = factor(MODEL_ORDER,  levels = MODEL_ORDER),
    fill = list(mean = NA_real_, se = NA_real_, n = 0)
  ) %>%
  mutate(
    available = !is.na(mean),
    mean_plot = if_else(available, mean, 0),
    se_plot   = if_else(available, se, 0)
  )

metrics_plot_named <- metrics_plot %>%
  left_join(prop_spec %>% select(prop, task_key, include_wt, include_smiles, placeholder_ok),
            by = c("task" = "task_key")) %>%
  mutate(
    prop = factor(prop, levels = PROP_ORDER)
  ) %>%
  filter(
    (repr == "Amino Acids"     & include_wt) |
      (repr == "SMILES" & include_smiles)
  )

# -------------------------
# Plot
# x-axis = metrics (F1/AUC/MCC/Accuracy)
# bars = models
# facets = prop × repr
# -------------------------
p_big <- ggplot(metrics_plot_named, aes(x = metric, y = mean_plot, fill = model_lbl)) +
  geom_col(position = position_dodge(width = 0.90), width = 0.85) +
  geom_errorbar(
    data = metrics_plot_named %>% filter(available & !is.na(se_plot)),
    aes(ymin = pmax(mean_plot - se_plot, 0), ymax = pmin(mean_plot + se_plot, 1)),
    position = position_dodge(width = 0.90),
    width = 0.18, linewidth = 0.30
  ) +
  geom_text(
    data = metrics_plot_named %>% filter(!available),
    aes(label = "N/A", y = 0.03),
    position = position_dodge(width = 0.90),
    size = 3
  ) +
  facet_grid(prop ~ repr, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PAL_MODEL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.01)) +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  )

print(p_big)

n_props <- dplyr::n_distinct(metrics_plot_named$prop)
out_png <- file.path(OUT_DIR, "supp_all_tasks_per_metric_bars.png")
out_pdf <- file.path(OUT_DIR, "supp_all_tasks_per_metric_bars.pdf")

ggsave(out_png, p_big, width = 12, height = 2.2 * n_props, dpi = 500)
ggsave(out_pdf, p_big, width = 12, height = 2.2 * n_props)

message("Saved:\n  ", out_png, "\n  ", out_pdf)
