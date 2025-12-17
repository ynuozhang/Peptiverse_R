suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(yardstick)
  library(stringr)
})
# -------------------------
# Config
# -------------------------
RAW_DIR  <- "raw"         # where you untar the server results
OUT_DIR  <- "data_processed"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
read_csv_fast <- function(path) {
  data.table::fread(path, data.table = FALSE, showProgress = FALSE)
}

# Parse task/model from the file path.
# Assumes structure like: .../<task>/<model_or_run>/<file>
# If your folder layout differs, this still gives useful defaults.
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
  # expects y_true, y_prob, y_pred columns
  df <- df %>%
    mutate(
      y_true = as.integer(y_true),
      y_pred = as.integer(y_pred),
      y_prob = as.numeric(y_prob)
    )
  
  # yardstick wants factors for class metrics
  d <- df %>%
    transmute(
      truth = factor(y_true, levels = c(0,1)),
      estimate = factor(y_pred, levels = c(0,1)),
      .pred_1 = y_prob
    )
  
  # Handle edge cases: AUC/AUPRC require both classes present
  has_both <- n_distinct(d$truth) == 2
  
  out <- list(
    f1       = tryCatch(f_meas(d, truth, estimate, beta = 1)$.estimate, error = function(e) NA_real_),
    accuracy = tryCatch(accuracy(d, truth, estimate)$.estimate,          error = function(e) NA_real_),
    mcc      = tryCatch(mcc(d, truth, estimate)$.estimate,               error = function(e) NA_real_)
  )
  
  if (has_both) {
    out$roc_auc <- tryCatch(roc_auc(d, truth, .pred_1)$.estimate, error = function(e) NA_real_)
    out$auprc   <- tryCatch(pr_auc(d, truth, .pred_1)$.estimate,  error = function(e) NA_real_)
  } else {
    out$roc_auc <- NA_real_
    out$auprc   <- NA_real_
  }
  
  tibble(metric = names(out), value = as.numeric(out))
}
# Extract "best threshold" if present in optimization_summary.txt
read_threshold_from_summary <- function(path) {
  txt <- readLines(path, warn = FALSE)
  # looks like: "Best threshold (picked on val): 0.1234"
  line <- txt[str_detect(txt, "Best threshold")]
  if (length(line) == 0) return(NA_real_)
  x <- str_extract(line[1], "[0-9]*\\.?[0-9]+")
  as.numeric(x)
}

# Best-so-far curve from study_trials.csv
make_best_so_far <- function(trials_df) {
  # optuna trials_dataframe() often has columns:
  # number, value, params_..., user_attrs_threshold, user_attrs_auc, user_attrs_ap
  # Sometimes trial index is "number" or "trial_number"
  idx_col <- dplyr::case_when(
    "number" %in% names(trials_df) ~ "number",
    "trial_number" %in% names(trials_df) ~ "trial_number",
    TRUE ~ names(trials_df)[1]
  )
  
  val_col <- dplyr::case_when(
    "value" %in% names(trials_df) ~ "value",
    TRUE ~ NA_character_
  )
  
  if (is.na(val_col)) {
    trials_df <- trials_df %>% mutate(value = NA_real_)
    val_col <- "value"
  }
  
  trials_df %>%
    arrange(.data[[idx_col]]) %>%
    mutate(
      trial_index = as.integer(.data[[idx_col]]),
      f1 = as.numeric(.data[[val_col]]),
      best_so_far = cummax(replace_na(f1, -Inf))
    )
}

# -------------------------
# 1) Discover files
# -------------------------
study_files <- list.files(RAW_DIR, pattern = "study_trials\\.csv$", recursive = TRUE, full.names = TRUE)
pred_files  <- list.files(RAW_DIR, pattern = "(train|val)_predictions\\.csv$", recursive = TRUE, full.names = TRUE)
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
  
  # optional: carry auc/ap if optuna stored as user attrs
  # these columns usually are named like "user_attrs_auc", "user_attrs_ap", "user_attrs_threshold"
  auc_col <- names(df2)[str_detect(names(df2), "user_attrs.*auc")]
  ap_col  <- names(df2)[str_detect(names(df2), "user_attrs.*ap")]
  thr_col <- names(df2)[str_detect(names(df2), "user_attrs.*threshold")]
  
  df2 %>%
    mutate(
      task  = meta$task,
      repr  = meta$repr,   # <-- add this
      model = meta$model,
      run   = meta$run,
      auc = if (length(auc_col) > 0) as.numeric(df2[[auc_col[1]]]) else NA_real_,
      ap  = if (length(ap_col)  > 0) as.numeric(df2[[ap_col[1]]])  else NA_real_,
      threshold = if (length(thr_col) > 0) as.numeric(df2[[thr_col[1]]]) else NA_real_
    ) %>%
    select(task, repr, model, run, trial_index, f1, best_so_far, auc, ap, threshold)
})

write_csv(trials_long, file.path(OUT_DIR, "trials_long.csv"))

# -------------------------
# 3) Read + standardize predictions and compute metrics
# -------------------------
preds_meta <- purrr::map_dfr(pred_files, function(f) {
  meta <- parse_meta_from_path(f, RAW_DIR)
  split <- ifelse(str_detect(basename(f), "^train_"), "train", "val")
  meta %>% mutate(split = split)
})

# For each (task, model, run), we’ll compute metrics per split from predictions csv
results_long <- preds_meta %>%
  distinct(path, task, model, run, split) %>%
  mutate(
    metrics_tbl = purrr::map(path, function(f) {
      df <- read_csv_fast(f)
      compute_metrics_from_preds(df)
    })
  ) %>%
  unnest(metrics_tbl) %>%
  arrange(task, model, run, split, metric)

write_csv(results_long, file.path(OUT_DIR, "results_long.csv"))

