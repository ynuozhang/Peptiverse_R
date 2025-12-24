suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
})

RAW_DIR <- "raw"                 # run this from the directory that contains nf/ hemolysis/ solubility/
OUT_PNG <- "supp_selected_hparams_cards.png"
OUT_PDF <- "supp_selected_hparams_cards.pdf"

read_csv_fast <- function(path) {
  data.table::fread(path, data.table = FALSE, showProgress = FALSE)
}

# Parse: <task>/<run>/<file>  e.g. nf/cnn_wt/study_trials.csv
parse_meta_simple <- function(path, raw_dir = RAW_DIR) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  
  rel <- sub(paste0("^", r, "/?"), "", p)
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  
  task_raw <- if (length(parts) >= 1) parts[1] else NA_character_
  run_raw  <- if (length(parts) >= 2) parts[2] else NA_character_
  
  # normalize task label similar to your other script
  task_l <- tolower(task_raw)
  task <- dplyr::case_when(
    task_l %in% c("nf", "nonfouling", "non_fouling", "non-fouling", "nonfoul") ~ "nonfouling",
    TRUE ~ str_replace_all(task_l, "[^a-z0-9]+", "")
  )
  
  run_l <- tolower(run_raw)
  repr <- dplyr::case_when(
    str_detect(run_l, "(^|[_-])smiles($|[_-])") ~ "smiles",
    str_detect(run_l, "(^|[_-])wt($|[_-])")     ~ "wt",
    TRUE ~ "unspecified"
  )
  
  tibble(path = path, task = task, run = run_raw, repr = repr)
}

# Extract best trial's params_* as long dataframe
read_best_params <- function(study_csv) {
  df <- read_csv_fast(study_csv)
  
  # objective column
  obj_col <- dplyr::case_when(
    "value" %in% names(df) ~ "value",
    "values_0" %in% names(df) ~ "values_0",
    TRUE ~ NA_character_
  )
  if (is.na(obj_col)) return(tibble())
  
  # prefer COMPLETE if state exists
  if ("state" %in% names(df)) df <- df[df$state == "COMPLETE", , drop = FALSE]
  df <- df[!is.na(df[[obj_col]]), , drop = FALSE]
  if (nrow(df) == 0) return(tibble())
  
  best_row <- df[which.max(df[[obj_col]]), , drop = FALSE]
  pcols <- names(best_row)[str_detect(names(best_row), "^params_")]
  if (length(pcols) == 0) return(tibble())
  
  tibble(best_row) %>%
    select(all_of(pcols)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "param",
      values_to = "chosen",
      values_transform = list(chosen = as.character)
    ) %>%
    mutate(
      param = sub("^params_", "", param),
      chosen = str_trim(chosen),
      # common xgb naming normalization
      param = recode(param,
                     "reg_lambda" = "lambda",
                     "reg_alpha"  = "alpha")
    )
}

# Nicely order params (optional)
param_order_hint <- c(
  # XGB-ish
  "learning_rate","max_depth","min_child_weight","subsample","colsample_bytree",
  "gamma","lambda","alpha","num_boost_round","early_stopping_rounds",
  # linear / svm / general
  "C","l1_ratio","tol","max_iter","kernel","class_weight",
  # NN-ish
  "lr","weight_decay","dropout","batch_size","layers","hidden","channels",
  "d_model","ff","nhead"
)

# -------------------------
# Discover study_trials.csv
# -------------------------
study_files <- list.files(RAW_DIR, pattern = "study_trials\\.csv$", recursive = TRUE, full.names = TRUE)
message("Found ", length(study_files), " study_trials.csv")

meta <- purrr::map_dfr(study_files, parse_meta_simple)

# -------------------------
# Build "cards" table
# -------------------------
best_params <- meta %>%
  mutate(params = purrr::map(path, read_best_params)) %>%
  unnest(params)

if (nrow(best_params) == 0) {
  stop("No params extracted. Check that study_trials.csv contains params_* columns.")
}

# Create a readable model label from run folder name
best_params <- best_params %>%
  mutate(
    run_l = tolower(run),
    model_lbl = case_when(
      str_detect(run_l, "xgb") ~ "XGB",
      str_detect(run_l, "enet|elastic") ~ "ENET",
      str_detect(run_l, "svm") ~ "SVM",
      str_detect(run_l, "mlp") ~ "MLP",
      str_detect(run_l, "cnn") ~ "CNN",
      str_detect(run_l, "transformer|bert|esm") ~ "Transformer",
      TRUE ~ run
    ),
    repr_lbl = case_when(
      repr == "wt" ~ "WT",
      repr == "smiles" ~ "SMILES",
      TRUE ~ repr
    ),
    panel = paste0(model_lbl, " (", repr_lbl, ")")
  )

# order params for prettier display
best_params <- best_params %>%
  mutate(
    param = factor(param, levels = unique(c(param_order_hint, sort(unique(as.character(param))))))
  ) %>%
  arrange(task, panel, run, param)

cards <- best_params %>%
  group_by(task, panel, run) %>%
  summarise(
    label = paste0(as.character(param), " = ", chosen, collapse = "\n"),
    .groups = "drop"
  ) %>%
  mutate(
    task = factor(task, levels = c("solubility","hemolysis","nonfouling")),
    panel = factor(panel, levels = c("XGB (WT)","ENET (WT)","SVM (WT)","MLP (WT)","CNN (WT)","Transformer (WT)",
                                     "XGB (SMILES)","ENET (SMILES)","SVM (SMILES)","MLP (SMILES)","CNN (SMILES)","Transformer (SMILES)"))
  )

# -------------------------
# Plot: facet "cards"
# -------------------------
p <- ggplot(cards, aes(x = 1, y = 1, label = label)) +
  geom_label(
    hjust = 0, vjust = 1,
    size = 2.6, lineheight = 0.95,
    label.size = 0.25
  ) +
  facet_grid(task ~ panel, scales = "free", space = "free") +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

print(p)

ggsave(OUT_PNG, p, width = 16, height = 9, dpi = 300)
ggsave(OUT_PDF, p, width = 16, height = 9)
message("Saved: ", OUT_PNG, " and ", OUT_PDF)
