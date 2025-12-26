suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
})

RAW_DIR  <- "raw"         
OUT_DIR  <- "data_processed"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv_fast <- function(path) {
  data.table::fread(path, data.table = FALSE, showProgress = FALSE)
}

make_best_so_far <- function(trials_df) {
  idx_col <- case_when(
    "number" %in% names(trials_df) ~ "number",
    "trial_number" %in% names(trials_df) ~ "trial_number",
    "trial" %in% names(trials_df) ~ "trial",
    TRUE ~ names(trials_df)[1]
  )
  
  val_col <- if ("value" %in% names(trials_df)) "value" else NA_character_
  if (is.na(val_col)) {
    trials_df <- trials_df %>% mutate(value = NA_real_)
    val_col <- "value"
  }
  
  trials_df %>%
    arrange(.data[[idx_col]]) %>%
    mutate(
      trial_index = as.integer(.data[[idx_col]]),
      rho = suppressWarnings(as.numeric(.data[[val_col]])),
      best_so_far = cummax(replace_na(rho, -Inf))
    )
}

# Binding affinity-specific path parser
# Expected layouts commonly look like:
#   raw/binding_affinity/<setting>/<mode>/optuna_trials.csv
# or
#   raw/binding_affinity/<setting>_<mode>/optuna_trials.csv
parse_binding_affinity_meta <- function(path, raw_dir) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  
  rel <- sub(paste0("^", r, "/?"), "", p)
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  
  # Find any token that looks like pooled/unpooled
  parts_l <- tolower(parts)
  mode <- NA_character_
  if (any(parts_l %in% c("pooled", "unpooled"))) {
    mode <- parts[which(parts_l %in% c("pooled", "unpooled"))[1]]
    mode <- tolower(mode)
  }
  
  # Identify a "setting" token:
  # - If a folder includes "setting" or "set" or looks like "s1/s2/s3/s4"
  # - Else fallback to the immediate parent folder of the CSV
  parent <- if (length(parts) >= 2) parts[length(parts) - 1] else "unknown"
  parent_l <- tolower(parent)
  
  setting <- case_when(
    str_detect(parent_l, "setting[ _-]*[0-9]+") ~ str_extract(parent_l, "setting[ _-]*[0-9]+"),
    str_detect(parent_l, "set[ _-]*[0-9]+")     ~ str_replace(str_extract(parent_l, "set[ _-]*[0-9]+"), "^set", "setting"),
    parent_l %in% c("s1","s2","s3","s4")        ~ str_replace(parent_l, "^s", "setting"),
    TRUE ~ parent_l
  )
  setting <- str_replace_all(setting, "[ _-]+", "_")
  
  # If mode wasn't found as its own directory, try extracting from parent folder name
  if (is.na(mode)) {
    if (str_detect(parent_l, "unpooled")) mode <- "unpooled"
    if (str_detect(parent_l, "pooled"))   mode <- "pooled"
  }
  
  # repr: keep as wt/smiles if embedded, else NA
  repr <- case_when(
    str_detect(tolower(rel), "(^|/|[_-])wt($|/|[_-])") ~ "wt",
    str_detect(tolower(rel), "(^|/|[_-])smiles($|/|[_-])") ~ "smiles",
    TRUE ~ NA_character_
  )
  
  # binding affinity task id
  task <- "binding_affinity"
  
  model <- "transformer"
  
  # run: use parent folder as run tag (keeps your 4-setting structure intact)
  run <- parent
  
  tibble(
    path = path,
    task = task,
    setting = setting,
    mode = mode,
    repr = repr,
    model = model,
    run = run
  )
}

# -------------------------
# 1) Discover Optuna trials files (binding affinity only)
# -------------------------
# primary file name in cross-attn script: optuna_trials.csv
trial_files1 <- list.files(RAW_DIR, pattern = "optuna_trials\\.csv$", recursive = TRUE, full.names = TRUE)
# fallback older name
trial_files2 <- list.files(RAW_DIR, pattern = "study_trials\\.csv$", recursive = TRUE, full.names = TRUE)

trial_files <- unique(c(trial_files1, trial_files2))

trial_files <- trial_files[str_detect(tolower(trial_files), "binding")]

message("Found binding-affinity trial files: ", length(trial_files))

if (length(trial_files) == 0) {
  stop("No binding-affinity optuna trial files found under RAW_DIR. ",
       "Expected something like raw/**/binding/**/optuna_trials.csv")
}

# -------------------------
# 2) Read + standardize
# -------------------------
ba_trials_long <- purrr::map_dfr(trial_files, function(f) {
  meta <- parse_binding_affinity_meta(f, RAW_DIR)
  df <- read_csv_fast(f)
  
  # Standardize objective value -> rho
  df2 <- make_best_so_far(df)
  
  df2 %>%
    mutate(
      task = meta$task,
      setting = meta$setting,
      mode = meta$mode,
      repr = meta$repr,
      model = meta$model,
      run = meta$run
    ) %>%
    select(task, setting, mode, repr, model, run, trial_index, rho, best_so_far)
})

# make repr explicit if missing (for your downstream code)
ba_trials_long <- ba_trials_long %>%
  mutate(
    repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr)),
    mode = case_when(is.na(mode) | mode == "" ~ "unknown", TRUE ~ tolower(mode)),
    setting = case_when(is.na(setting) | setting == "" ~ "unknown", TRUE ~ tolower(setting))
  )

out_csv <- file.path(OUT_DIR, "binding_affinity_trials_long.csv")
readr::write_csv(ba_trials_long, out_csv)
message("Saved: ", out_csv)

print(ba_trials_long %>%
        count(setting, mode, repr, sort = TRUE))
print(ba_trials_long %>%
        group_by(setting, mode, repr) %>%
        summarise(best_rho = max(rho, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(best_rho)))
