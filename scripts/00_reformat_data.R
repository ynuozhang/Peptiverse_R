suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(yardstick)
  library(stringr)
})
# -------------------------
# Config
# -------------------------
RAW_DIR  <- "raw"
OUT_DIR  <- "data_processed"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
read_csv_fast <- function(path) {
  data.table::fread(path, data.table = FALSE, showProgress = FALSE)
}

# Parse task/model from the file path.
# Assumes structure like: .../<task>/<model_or_run>/<file>
parse_meta_from_path <- function(path, raw_dir) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  
  rel <- sub(paste0("^", r, "/?"), "", p)
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  
  # expected common layout: <task>/<run>/<file>
  task_raw <- if (length(parts) >= 1) parts[1] else NA_character_
  run_raw  <- if (length(parts) >= 2) parts[2] else "default"
  file     <- tail(parts, 1)
  
  task_l <- tolower(task_raw)
  # normalize task names
  task <- dplyr::case_when(
    task_l %in% c("nf", "nonfouling", "non_fouling", "non-fouling", "nonfoul") ~ "nonfouling",
    task_l %in% c("halflife", "half_life", "half-life", "t12") ~ "halflife",
    TRUE ~ str_replace_all(task_l, "[^a-z0-9]+", "")
  )
  
  run_l <- tolower(run_raw)
  
  # repr can be embedded in run folder name like enet_gpu_smiles / svm_gpu_wt
  repr <- dplyr::case_when(
    str_detect(run_l, "(^|[_-])smiles($|[_-])") ~ "smiles",
    str_detect(run_l, "(^|[_-])wt($|[_-])")     ~ "wt",
    TRUE ~ NA_character_
  )
  
  # model = run with repr stripped off
  model <- run_raw
  if (!is.na(repr)) {
    model <- str_replace_all(model, "(^|[_-])smiles($|[_-])", "\\1")
    model <- str_replace_all(model, "(^|[_-])wt($|[_-])", "\\1")
    model <- str_replace_all(model, "[-_]+$", "")
  }
  
  # normalize model names (expand as needed)
  model_norm <- dplyr::case_when(
    str_detect(tolower(model), "svm_gpu|cuml.*svc|cusvc") ~ "svm_gpu",
    str_detect(tolower(model), "enet_gpu|elastic|cuml.*logreg|culogreg|logreg") ~ "enet_gpu",
    str_detect(tolower(model), "\\bxgb\\b|xgboost") ~ "xgb",
    str_detect(tolower(model), "\\bmlp\\b") ~ "mlp",
    str_detect(tolower(model), "\\bcnn\\b") ~ "cnn",
    str_detect(tolower(model), "transformer|esm|bert") ~ "transformer",
    TRUE ~ model
  )
  
  tibble::tibble(
    path = path, filename = file,
    task = task, run = run_raw, model = model_norm, repr = repr
  )
}


# Safe metric computations from y_true / y_prob / y_pred
compute_metrics_from_preds <- function(df) {
  # normalize header names (handles BOM + whitespace)
  names(df) <- names(df) %>%
    stringr::str_replace("^\ufeff", "") %>%  # remove BOM if present
    stringr::str_trim()
  nms <- names(df)
  
  # -------- find y_true ----------
  y_true_col <- intersect(nms, c("y_true", "y", "label", "target"))[1]
  if (is.na(y_true_col)) {
    stop("Missing y_true column. Expected one of: y_true, y, label, target")
  }
  y_true <- suppressWarnings(as.numeric(df[[y_true_col]]))
  
  # -------- find prediction columns ----------
  prob_col <- intersect(nms, c("y_prob", "prob", "p", "p_hat", "pred_prob", ".pred_1"))[1]
  pred_col <- intersect(nms, c("y_pred", "pred", "y_hat", "prediction", "pred_mean"))[1]
  
  y_prob <- if (!is.na(prob_col)) suppressWarnings(as.numeric(df[[prob_col]])) else NA_real_
  y_pred <- if (!is.na(pred_col)) suppressWarnings(as.numeric(df[[pred_col]])) else NA_real_
  
  # -------- decide classification vs regression ----------
  # classification if y_true is 0/1 AND we have probabilities
  is_binary <- all(na.omit(y_true) %in% c(0, 1)) && !is.na(prob_col)
  
  if (is_binary) {
    if (is.na(pred_col)) {
      y_hat <- as.integer(y_prob >= 0.5)
    } else {
      y_hat <- as.integer(y_pred)
    }
    
    d <- tibble::tibble(
      truth    = factor(as.integer(y_true), levels = c(0, 1)),
      estimate = factor(y_hat, levels = c(0, 1)),
      .pred_1  = y_prob
    )
    
    has_both <- dplyr::n_distinct(d$truth) == 2
    
    out <- list(
      f1       = tryCatch(yardstick::f_meas(d, truth, estimate, beta = 1)$.estimate, error = function(e) NA_real_),
      accuracy = tryCatch(yardstick::accuracy(d, truth, estimate)$.estimate,         error = function(e) NA_real_),
      mcc      = tryCatch(yardstick::mcc(d, truth, estimate)$.estimate,              error = function(e) NA_real_)
    )
    
    if (has_both) {
      out$roc_auc <- tryCatch(yardstick::roc_auc(d, truth, .pred_1)$.estimate, error = function(e) NA_real_)
      out$auprc   <- tryCatch(yardstick::pr_auc(d, truth, .pred_1)$.estimate,  error = function(e) NA_real_)
    } else {
      out$roc_auc <- NA_real_
      out$auprc   <- NA_real_
    }
    
    return(tibble::tibble(metric = names(out), value = as.numeric(out)))
  }
  
  # -------- regression branch ----------
  if (is.na(pred_col)) {
    stop("Regression preds: missing prediction column. Expected one of: y_pred, pred, y_hat, prediction, pred_mean")
  }
  
  # basic regression metrics + Spearman rho
  ok <- is.finite(y_true) & is.finite(y_pred)
  yt <- y_true[ok]; yp <- y_pred[ok]
  
  rho  <- suppressWarnings(cor(yt, yp, method = "spearman"))
  rmse <- sqrt(mean((yt - yp)^2))
  mae  <- mean(abs(yt - yp))
  
  tibble::tibble(
    metric = c("spearman_rho", "rmse", "mae"),
    value  = c(as.numeric(rho), as.numeric(rmse), as.numeric(mae))
  )
}


# Best-so-far curve from study_trials.csv
make_best_so_far <- function(trials_df) {
  # trial index column
  idx_col <- dplyr::case_when(
    "number" %in% names(trials_df) ~ "number",
    "trial_number" %in% names(trials_df) ~ "trial_number",
    TRUE ~ names(trials_df)[1]
  )
  
  # choose what should populate "f1"
  # priority: explicit CV spearman (half-life XGB), else Optuna "value"
  score_col <- dplyr::case_when(
    "user_attrs_cv_spearman_rho" %in% names(trials_df) ~ "user_attrs_cv_spearman_rho",
    "value" %in% names(trials_df) ~ "value",
    TRUE ~ NA_character_
  )
  
  if (is.na(score_col)) {
    trials_df <- trials_df %>% mutate(value = NA_real_)
    score_col <- "value"
  }
  
  trials_df %>%
    arrange(.data[[idx_col]]) %>%
    mutate(
      trial_index = as.integer(.data[[idx_col]]),
      # keep the column name f1 for compatibility (even if it's rho)
      f1 = as.numeric(.data[[score_col]]),
      best_so_far = cummax(replace_na(f1, -Inf))
    )
}



# -------------------------
# 1) Discover files
# -------------------------
study_files <- list.files(RAW_DIR, pattern = "study_trials\\.csv$", recursive = TRUE, full.names = TRUE)
pred_files <- list.files(
  RAW_DIR,
  pattern = "((train|val)_predictions|oof_predictions)\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
# Exclude half-life SMILES prediction files (schema differs, and we don't use them here)
pred_files <- pred_files[!(
  str_detect(tolower(pred_files), "(^|/)half_life(/|$)") &
    str_detect(tolower(pred_files), "(^|/)([^/]*smiles[^/]*)/")  # run folder contains "smiles"
)]

sum_files   <- list.files(RAW_DIR, pattern = "optimization_summary\\.txt$", recursive = TRUE, full.names = TRUE)

message("Found: ",
        length(study_files), " study_trials.csv | ",
        length(pred_files),  " *_predictions.csv | ",
        length(sum_files),   " optimization_summary.txt")

# -------------------------
# 2) Read + standardize trials
# -------------------------
trials_long <- purrr::map_dfr(study_files, function(f) {
  meta <- parse_meta_from_path(f, RAW_DIR)
  df <- read_csv_fast(f)
  
  df2 <- make_best_so_far(df)
  
  # carry auc/ap if optuna stored as user attrs
  auc_col <- names(df2)[str_detect(names(df2), "user_attrs.*auc")]
  ap_col  <- names(df2)[str_detect(names(df2), "user_attrs.*ap")]
  thr_col <- names(df2)[str_detect(names(df2), "user_attrs.*threshold")]
  
  kind_hint <- if (!is.na(meta$task) && meta$task == "halflife") {
    "regression"
  } else if (any(str_detect(names(df), "^user_attrs_cv_(spearman_rho|rmse|mae|r2)$"))) {
    "regression"
  } else {
    NA_character_
  }
  
  df2 %>%
    mutate(
      task  = meta$task,
      repr  = meta$repr,
      model = meta$model,
      run   = meta$run,
      kind_hint = kind_hint,
      auc = if (length(auc_col) > 0) as.numeric(df2[[auc_col[1]]]) else NA_real_,
      ap  = if (length(ap_col)  > 0) as.numeric(df2[[ap_col[1]]])  else NA_real_,
      threshold = if (length(thr_col) > 0) as.numeric(df2[[thr_col[1]]]) else NA_real_
    ) %>%
    select(task, repr, model, run, trial_index, f1, best_so_far, auc, ap, threshold, kind_hint)
  
})

write_csv(trials_long, file.path(OUT_DIR, "trials_long.csv"))

# -------------------------
# 3) Read + standardize predictions and compute metrics
# -------------------------
preds_meta <- purrr::map_dfr(pred_files, function(f) {
  meta <- parse_meta_from_path(f, RAW_DIR)
  split <- dplyr::case_when(
    str_detect(basename(f), "^train_") ~ "train",
    str_detect(basename(f), "^val_")   ~ "val",
    str_detect(basename(f), "^oof_") | basename(f) == "oof_predictions.csv" ~ "oof",
    TRUE ~ "unknown"
  )
  meta %>% mutate(split = split)
})

# For each (task, repr, model, run), compute metrics per split
results_long <- preds_meta %>%
  distinct(path, task, repr, model, run, split) %>% 
  mutate(
    metrics_tbl = purrr::map(path, function(f) {
      message("METRICS FILE: ", f)
      df <- read_csv_fast(f)
      compute_metrics_from_preds(df)
    })
  ) %>%
  unnest(metrics_tbl) %>%
  arrange(task, repr, model, run, split, metric)


write_csv(results_long, file.path(OUT_DIR, "results_long.csv"))

# -------------------------
# 4) Split outputs into classification vs regression CSVs
# -------------------------

results_long2 <- results_long %>% mutate(repr2 = dplyr::coalesce(repr, "unknown"))
trials_long2  <- trials_long  %>% mutate(repr2 = dplyr::coalesce(repr, "unknown"))

kind_from_results <- results_long2 %>%
  group_by(task, repr2, model, run) %>%
  summarise(
    kind = dplyr::case_when(
      any(metric %in% c("spearman_rho", "rmse", "mae")) ~ "regression",
      any(metric %in% c("f1", "accuracy", "mcc", "roc_auc", "auprc")) ~ "classification",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

kind_from_trials <- trials_long2 %>%
  group_by(task, repr2, model, run) %>%
  summarise(
    kind = dplyr::first(stats::na.omit(kind_hint)),
    .groups = "drop"
  )

run_kind <- dplyr::full_join(
  kind_from_results, kind_from_trials,
  by = c("task", "repr2", "model", "run"),
  suffix = c("_res", "_trial")
) %>%
  transmute(
    task, repr2, model, run,
    kind = dplyr::coalesce(kind_res, kind_trial, "unknown")
  )

results_classification <- results_long2 %>%
  left_join(run_kind, by = c("task", "repr2", "model", "run")) %>%
  filter(kind == "classification") %>%
  select(-kind, -repr2)

results_regression <- results_long2 %>%
  left_join(run_kind, by = c("task", "repr2", "model", "run")) %>%
  filter(kind == "regression") %>%
  select(-kind, -repr2)

trials_classification <- trials_long2 %>%
  left_join(run_kind, by = c("task", "repr2", "model", "run")) %>%
  filter(kind == "classification") %>%
  select(-kind, -repr2, -kind_hint)

trials_regression <- trials_long2 %>%
  left_join(run_kind, by = c("task", "repr2", "model", "run")) %>%
  filter(kind == "regression") %>%
  select(-kind, -repr2, -kind_hint)


write_csv(results_classification, file.path(OUT_DIR, "results_classification.csv"))
write_csv(results_regression,     file.path(OUT_DIR, "results_regression.csv"))
write_csv(trials_classification,  file.path(OUT_DIR, "trials_classification.csv"))
write_csv(trials_regression,      file.path(OUT_DIR, "trials_regression.csv"))

