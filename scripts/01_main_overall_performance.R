ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(scales)
  library(dplyr)
  library(forcats)
  library(jsonlite)
})

# ============================================================
# Configuration
# ============================================================
OUT_DIR   <- "figures/main"
RAW_DIR   <- "raw"
BA_CSV    <- "data_processed/binding_affinity_trials_long.csv"
IN_TRIALS_CLS <- "data_processed/trials_classification.csv"
IN_TRIALS_REG <- "data_processed/trials_regression.csv"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAPER_REG_TASKS <- c("halflife", "permeabilitypampa", "permeabilitycaco2")

# ============================================================
# Helper Functions
# ============================================================

read_csv_fast <- function(path) {
  data.table::fread(path, data.table = FALSE, showProgress = FALSE)
}

norm_task <- function(x) {
  stringr::str_replace_all(tolower(x %||% ""), "[^a-z0-9]+", "")
}

canon_task <- function(task) {
  tn <- norm_task(task)
  dplyr::case_when(
    stringr::str_detect(tn, "halflife|halflifesecond|halflifehour|half_life") ~ "halflife",
    stringr::str_detect(tn, "permeabilitypampa") ~ "permeabilitypampa",
    stringr::str_detect(tn, "permeabilitycaco2|permeabilitycaco") ~ "permeabilitycaco2",
    stringr::str_detect(tn, "bindingaffinity") ~ "bindingaffinity",
    TRUE ~ tn
  )
}

keep_paper_reg_task <- function(task) {
  canon_task(task) %in% PAPER_REG_TASKS
}

task_matches <- function(task_col, key) {
  norm <- function(x) tolower(x) %>% str_replace_all("[^a-z0-9]+", "")
  str_detect(norm(task_col), fixed(norm(key)))
}

task_matches_eq <- function(task_col, key) {
  norm <- function(x) tolower(x) %>% stringr::str_replace_all("[^a-z0-9]+", "")
  norm(task_col) == norm(key)
}

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  x <- str_remove(x, "_log$")
  x <- str_remove(x, "_raw$")
  x <- str_remove(x, "_unpooled$")
  x <- str_remove(x, "_pooled$")
  x <- str_remove(x, "_wt$")  
  x <- str_remove(x, "_smiles$") 
  
  dplyr::case_when(
    str_detect(x, "svr") ~ "SVR",
    str_detect(x, "svm") ~ "SVM",
    str_detect(x, "xgb|xgboost") ~ "XGB",
    str_detect(x, "enet|elastic|logreg") ~ "ENET",
    str_detect(x, "\\bmlp\\b") ~ "MLP",
    str_detect(x, "\\bcnn\\b") ~ "CNN",
    str_detect(x, "transformer|bert|esm") ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}

parse_meta_from_path <- function(path, raw_dir) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  
  rel <- sub(paste0("^", r, "/?"), "", p)
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  
  task_raw <- if (length(parts) >= 1) parts[1] else NA_character_
  run_raw  <- if (length(parts) >= 2) parts[2] else "default"
  file     <- tail(parts, 1)
  
  task_l <- tolower(task_raw)
  task <- dplyr::case_when(
    task_l %in% c("nf", "nonfouling", "non_fouling", "non-fouling", "nonfoul") ~ "nonfouling",
    task_l %in% c("halflife", "half_life", "half-life", "t12") ~ "halflife",
    TRUE ~ str_replace_all(task_l, "[^a-z0-9]+", "")
  )
  
  run_l <- tolower(run_raw)
  
  repr <- dplyr::case_when(
    str_detect(run_l, "_smiles") ~ "smiles",
    str_detect(run_l, "_wt") ~ "wt",
    TRUE ~ NA_character_
  )
  
  model <- run_raw
  if (!is.na(repr)) {
    model <- str_replace(model, "_wt.*$", "")
    model <- str_replace(model, "_smiles.*$", "")
  }
  
  model_norm <- dplyr::case_when(
    str_detect(tolower(model), "svm") ~ "svm",
    str_detect(tolower(model), "xgb") ~ "xgb",
    str_detect(tolower(model), "transformer") ~ "transformer",
    str_detect(tolower(model), "cnn") ~ "cnn",
    str_detect(tolower(model), "mlp") ~ "mlp",
    str_detect(tolower(model), "enet|elastic") ~ "enet",
    TRUE ~ model
  )
  
  tibble::tibble(
    path = path, filename = file,
    task = task, run = run_raw, model = model_norm, repr = repr
  )
}

pick_first <- function(nm, cands) {
  hit <- intersect(cands, nm)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ============================================================
# Binding Affinity Override
# ============================================================

ba_best <- readr::read_csv(BA_CSV, show_col_types = FALSE) %>%
  mutate(
    setting_l = tolower(setting),
    repr2 = case_when(
      str_detect(setting_l, "wt_smiles") ~ "smiles",
      str_detect(setting_l, "wt_wt")     ~ "wt",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(repr2)) %>%
  group_by(repr2) %>%
  summarise(
    mean_best = max(rho, na.rm = TRUE),
    model_lbl = "Transformer",
    .groups = "drop"
  ) %>%
  transmute(
    repr = repr2,
    mean_best = mean_best,
    model_lbl = model_lbl,
    task = "binding_affinity",
    kind = "regression",
    available = TRUE
  )

print(ba_best)

# ============================================================
# Read Trials Data
# ============================================================

read_trials_one <- function(path, kind_tag) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  
  rho_candidates <- c("spearman", "spearman_rho", "rho", "spearman_corr", "corr_spearman")
  rho_hit <- intersect(rho_candidates, names(df))
  
  if (kind_tag == "regression" && length(rho_hit) == 0 && "f1" %in% names(df)) {
    message("[regression] No rho column found; treating 'f1' as Spearman rho for: ", path)
    df <- df %>% mutate(rho = suppressWarnings(as.numeric(f1)))
  } else if (length(rho_hit) > 0) {
    df <- df %>% mutate(rho = suppressWarnings(as.numeric(.data[[rho_hit[1]]])))
  } else {
    df <- df %>% mutate(rho = NA_real_)
  }
  
  obj_col <- if (kind_tag == "regression") "rho" else if ("f1" %in% names(df)) "f1" else "best_so_far"
  
  message("[", kind_tag, "] objective = ", obj_col, "  (", path, ")")
  
  df %>%
    mutate(
      kind = kind_tag,
      task  = as.character(task),
      model = as.character(model),
      run   = as.character(run),
      repr  = as.character(repr),
      trial_index = as.integer(trial_index),
      value = suppressWarnings(as.numeric(.data[[obj_col]]))
    ) %>%
    mutate(
      repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr))
    )
}

trials_cls <- read_trials_one(IN_TRIALS_CLS, "classifier")
trials_reg <- read_trials_one(IN_TRIALS_REG, "regression")

# ============================================================
# Data Collection Functions
# ============================================================

spearman_from_predfile <- function(path) {
  d <- read_csv_fast(path)
  nm <- names(d)
  
  truth_col <- pick_first(nm, c("y_true","truth","label","y","target","gt"))
  pred_col  <- pick_first(nm, c("y_pred","pred","y_hat","prediction","pred_mean","estimate","y_pred_reg"))
  
  if (is.na(truth_col) || is.na(pred_col)) return(NA_real_)
  
  yt <- suppressWarnings(as.numeric(d[[truth_col]]))
  yp <- suppressWarnings(as.numeric(d[[pred_col]]))
  ok <- is.finite(yt) & is.finite(yp)
  if (sum(ok) < 3) return(NA_real_)
  
  as.numeric(suppressWarnings(cor(yt[ok], yp[ok], method = "spearman")))
}

read_cv_oof_spearman <- function(path) {
  x <- jsonlite::read_json(path, simplifyVector = TRUE)
  
  cand <- c(
    suppressWarnings(as.numeric(x[["oof_metrics"]][["spearman_rho"]])),
    suppressWarnings(as.numeric(x[["spearman_rho"]])),
    suppressWarnings(as.numeric(x[["spearman"]])),
    suppressWarnings(as.numeric(x[["best_objective_cv_spearman"]]))
  )
  cand <- cand[is.finite(cand)]
  if (length(cand) == 0) return(NA_real_)
  cand[1]
}

read_opt_summary_spearman <- function(path) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  all_txt <- paste(txt, collapse = "\n")
  
  model <- NA_character_
  i_mod <- grep("^\\s*MODEL\\s*:\\s*", txt, ignore.case = TRUE)
  if (length(i_mod) > 0) {
    model <- trimws(sub("^\\s*MODEL\\s*:\\s*", "", txt[i_mod[1]], ignore.case = TRUE))
  } else {
    folder <- basename(dirname(path))
    model <- str_replace(folder, "_wt.*$", "")
    model <- str_replace(model, "_smiles.*$", "")
  }
  
  sp <- NA_real_
  
  # Format 1: Inside OOF metrics JSON block
  m <- regexpr("\"spearman_rho\"\\s*:\\s*([-+0-9\\.eE]+)", all_txt, perl = TRUE)
  if (m[1] != -1) {
    hit <- regmatches(all_txt, m)
    val <- suppressWarnings(as.numeric(sub(".*:\\s*", "", hit)))
    if (is.finite(val)) sp <- val
  }
  
  # Format 2: For XGB style
  if (!is.finite(sp)) {
    m <- regexpr("\"best_objective_cv_spearman\"\\s*:\\s*([-+0-9\\.eE]+)", all_txt, perl = TRUE)
    if (m[1] != -1) {
      hit <- regmatches(all_txt, m)
      val <- suppressWarnings(as.numeric(sub(".*:\\s*", "", hit)))
      if (is.finite(val)) sp <- val
    }
  }
  
  # Format 3: ML format
  if (!is.finite(sp)) {
    i_line <- grep("Val\\s+Spearman\\s+rho\\s*\\(objective\\)\\s*:", txt, ignore.case = TRUE)
    if (length(i_line) > 0) {
      s <- txt[i_line[1]]
      val <- suppressWarnings(as.numeric(sub(".*:\\s*([-+0-9\\.eE]+)\\s*$", "\\1", s)))
      if (is.finite(val)) sp <- val
    }
  }
  
  tibble(model = model, spearman = sp)
}

# A) Collect from val_predictions.csv
collect_reg_from_val_preds <- function(RAW_DIR) {
  val_files <- list.files(RAW_DIR, pattern = "val_predictions\\.csv$", recursive = TRUE, full.names = TRUE)
  
  purrr::map_dfr(val_files, function(f) {
    meta <- parse_meta_from_path(f, RAW_DIR)
    
    if (is.na(meta$task) || !keep_paper_reg_task(meta$task)) return(tibble())
    
    repr2 <- dplyr::case_when(
      !is.na(meta$repr) ~ tolower(meta$repr),
      stringr::str_detect(tolower(f), "(^|/|_)smiles($|/|_)") ~ "smiles",
      stringr::str_detect(tolower(f), "(^|/|_)wt($|/|_)") ~ "wt",
      TRUE ~ "unspecified"
    )
    
    sp <- spearman_from_predfile(f)
    if (!is.finite(sp)) return(tibble())
    
    tibble(
      source    = "val_predictions",
      kind      = "regression",
      task      = canon_task(meta$task),
      repr      = repr2,
      model     = as.character(meta$model),
      model_lbl = clean_model_label(meta$model),
      run       = as.character(meta$run),
      score     = sp
    )
  })
}

# B) Collect from cv_oof_summary.json
collect_reg_from_oof_summary <- function(RAW_DIR) {
  js <- list.files(RAW_DIR, pattern = "cv_oof_summary\\.json$", recursive = TRUE, full.names = TRUE)
  
  purrr::map_dfr(js, function(f) {
    meta <- parse_meta_from_path(f, RAW_DIR)
    
    if (is.na(meta$task) || !keep_paper_reg_task(meta$task)) return(tibble())
    
    repr2 <- dplyr::case_when(
      stringr::str_detect(tolower(f), "(^|/|_)smiles($|/|_)") ~ "smiles",
      stringr::str_detect(tolower(f), "(^|/|_)wt($|/|_)") ~ "wt",
      TRUE ~ NA_character_
    )
    if (is.na(repr2)) return(tibble())
    
    sp <- read_cv_oof_spearman(f)
    if (!is.finite(sp)) return(tibble())
    
    tibble(
      source    = "cv_oof_summary",
      kind      = "regression",
      task      = canon_task(meta$task),
      repr      = repr2,
      model     = as.character(meta$model),
      model_lbl = clean_model_label(meta$model),
      run       = as.character(meta$run),
      score     = sp
    )
  })
}

# C) Collect from opt_summary*.txt
collect_reg_from_opt_summaries <- function(RAW_DIR) {
  sum_files <- list.files(RAW_DIR, pattern = "opt_summary\\.txt$|summary\\.txt$|optuna_summary\\.txt$", 
                          recursive = TRUE, full.names = TRUE)
  
  purrr::map_dfr(sum_files, function(f) {
    meta <- parse_meta_from_path(f, RAW_DIR)
    if (is.na(meta$task) || !keep_paper_reg_task(meta$task)) return(tibble())
    
    repr2 <- dplyr::case_when(
      !is.na(meta$repr) ~ tolower(meta$repr),
      stringr::str_detect(tolower(f), "_smiles") ~ "smiles",
      stringr::str_detect(tolower(f), "_wt") ~ "wt",
      stringr::str_detect(tolower(f), "/smiles/") ~ "smiles", 
      stringr::str_detect(tolower(f), "/wt/") ~ "wt",
      TRUE ~ NA_character_
    )
    if (is.na(repr2)) return(tibble())
    
    sp_tbl <- read_opt_summary_spearman(f)
    sp <- sp_tbl$spearman[1]
    if (!is.finite(sp)) return(tibble())
    
    tibble(
      source    = "opt_summary",
      kind      = "regression",
      task      = canon_task(meta$task),
      repr      = repr2,
      model     = sp_tbl$model[1],
      model_lbl = clean_model_label(sp_tbl$model[1]),
      run       = as.character(meta$run),
      score     = sp
    )
  })
}

# D) Collect from trials_regression.csv
collect_reg_from_trials <- function(IN_TRIALS_REG) {
  df <- readr::read_csv(IN_TRIALS_REG, show_col_types = FALSE)
  
  df <- df %>% mutate(task = as.character(task))
  df <- df %>% mutate(task_c = canon_task(task)) %>% filter(task_c %in% PAPER_REG_TASKS)
  
  rho_candidates <- c("spearman", "spearman_rho", "rho", "spearman_corr", "corr_spearman")
  rho_hit <- intersect(rho_candidates, names(df))
  
  if (length(rho_hit) > 0) {
    df <- df %>% mutate(rho = suppressWarnings(as.numeric(.data[[rho_hit[1]]])))
  } else if ("f1" %in% names(df)) {
    df <- df %>% mutate(rho = suppressWarnings(as.numeric(f1)))
  } else {
    df <- df %>% mutate(rho = NA_real_)
  }
  
  df %>%
    mutate(
      source    = "trials_regression",
      kind      = "regression",
      task      = task_c,
      repr      = dplyr::case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(as.character(repr))),
      model     = as.character(model),
      model_lbl = clean_model_label(model),
      run       = as.character(run),
      score     = suppressWarnings(as.numeric(rho))
    ) %>%
    filter(is.finite(score)) %>%
    select(source, kind, task, repr, model, model_lbl, run, score)
}

# ============================================================
# Pick Best Regressors
# ============================================================

pick_best_regressors <- function(RAW_DIR, IN_TRIALS_REG) {
  cand_val <- collect_reg_from_val_preds(RAW_DIR)
  cand_oof <- collect_reg_from_oof_summary(RAW_DIR)
  cand_sum <- collect_reg_from_opt_summaries(RAW_DIR)
  cand_tr  <- collect_reg_from_trials(IN_TRIALS_REG)
  
  source_rank <- c(
    opt_summary      = 1L,
    cv_oof_summary   = 1L,
    val_predictions  = 2L,
    trials_regression= 3L
  )
  
  cands <- bind_rows(cand_val, cand_oof, cand_sum, cand_tr) %>%
    mutate(src_rank = source_rank[source]) %>%
    filter(!is.na(src_rank))
  
  model_tbl <- cands %>%
    group_by(task, repr, model_lbl, source, src_rank) %>%
    summarise(
      mean_best = max(score, na.rm = TRUE),
      med_best  = median(score, na.rm = TRUE),
      n_runs    = n_distinct(run),
      .groups   = "drop"
    )
  
  model_tbl %>%
    group_by(task, repr) %>%
    arrange(src_rank, desc(mean_best), desc(n_runs), desc(med_best)) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      kind      = "regression",
      task,
      repr,
      model_lbl,
      mean_best,
      med_best,
      n_runs,
      available = TRUE
    )
}

best_reg <- pick_best_regressors(RAW_DIR, IN_TRIALS_REG)

best_reg %>%
  arrange(task, repr) %>%
  select(task, repr, model_lbl, mean_best, n_runs) %>%
  print(n = 200)

# ============================================================
# Best Model Processing
# ============================================================

best_by_run_metric <- function(df) {
  df %>%
    filter(!is.na(value)) %>%
    group_by(kind, task, repr, model, run) %>%
    summarise(best_metric = max(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(model_lbl = clean_model_label(model))
}

best_model_per_task <- function(best_by_run_df) {
  best_by_run_df %>%
    group_by(kind, task, repr, model_lbl) %>%
    summarise(
      mean_best = mean(best_metric),
      med_best  = median(best_metric),
      n_runs    = n_distinct(run),
      sd_best   = sd(best_metric),
      se_best   = sd_best / sqrt(n_runs),
      .groups = "drop"
    ) %>%
    arrange(kind, task, repr, desc(mean_best), desc(med_best), desc(n_runs)) %>%
    group_by(kind, task, repr) %>%
    slice(1) %>%
    ungroup()
}

best_model_per_task_all <- best_by_run_metric(trials_cls) %>%
  best_model_per_task() %>%
  bind_rows(best_reg)

# ============================================================
# Property Specification
# ============================================================

prop_spec <- tribble(
  ~prop,                      ~task_key,                 ~kind,          ~include_wt, ~include_smiles, ~placeholder_ok,
  "Hemolysis",                "hemolysis",               "classifier",    TRUE,        TRUE,            TRUE,
  "Non-Fouling",              "nonfouling",              "classifier",    TRUE,        TRUE,            TRUE,
  "Toxicity",                 "toxicity",                "classifier",    FALSE,       TRUE,            TRUE,
  "Solubility",               "solubility",              "classifier",    TRUE,        FALSE,           FALSE,
  "Permeability",             "permeabilitypenetrance",  "classifier",    TRUE,        FALSE,           TRUE,
  "Half-Life",                "halflife",                "regression",    TRUE,        TRUE,            TRUE,
  "Permeability (PAMPA)",     "permeabilitypampa",       "regression",    FALSE,       TRUE,            TRUE,
  "Permeability (CACO2)",     "permeabilitycaco2",       "regression",    FALSE,       TRUE,            TRUE,
  "Binding Affinity",         "bindingaffinity",         "regression",    TRUE,        TRUE,            TRUE
)

# ============================================================
# Build Panel Data
# ============================================================

build_panel_df <- function(repr_key, panel_label) {
  spec2 <- prop_spec %>%
    mutate(include = if (repr_key == "wt") include_wt else include_smiles) %>%
    filter(include) %>%
    mutate(prop = as.character(prop)) %>%
    arrange(prop)
  
  out <- purrr::pmap_dfr(
    spec2 %>% select(prop, task_key, kind, placeholder_ok),
    function(prop, task_key, kind_spec, placeholder_ok) {
      
      if (task_matches_eq(task_key, "bindingaffinity")) {
        dba <- ba_best %>% filter(repr == repr_key)
        
        if (nrow(dba) == 0) {
          if (!placeholder_ok) return(NULL)
          return(tibble(
            repr = repr_key,
            panel = panel_label,
            prop = as.character(prop),
            kind = "regression",
            available = FALSE,
            value = NA_real_,
            model_lbl = "Missing"
          ))
        }
        
        return(tibble(
          repr = repr_key,
          panel = panel_label,
          prop = as.character(prop),
          kind = "regression",
          available = TRUE,
          value = dba$mean_best[1],
          model_lbl = dba$model_lbl[1]
        ))
      }
      
      d <- best_model_per_task_all %>%
        filter(repr == repr_key) %>%
        filter(kind == kind_spec) %>%
        filter(task_matches_eq(task, task_key))
      
      if (nrow(d) == 0) {
        if (!placeholder_ok) return(NULL)
        return(tibble(
          repr = repr_key,
          panel = panel_label,
          prop = as.character(prop),
          kind = kind_spec,
          available = FALSE,
          value = NA_real_,
          model_lbl = "Missing"
        ))
      }
      
      d <- d %>% arrange(desc(mean_best)) %>% slice(1)
      
      tibble(
        repr = repr_key,
        panel = panel_label,
        prop = as.character(prop),
        kind = kind_spec,
        available = TRUE,
        value = d$mean_best,
        model_lbl = d$model_lbl
      )
    }
  )
  
  out
}

df_wt     <- build_panel_df("wt",     "Amino Acids Predictors")
df_smiles <- build_panel_df("smiles", "SMILES Predictors")

plot_df <- bind_rows(df_wt, df_smiles) %>%
  mutate(panel = factor(panel, levels = c("Amino Acids Predictors", "SMILES Predictors")))

# ============================================================
# Panel-specific ordering
# ============================================================

order_wt <- c(
  "Hemolysis",
  "Non-Fouling",
  "Solubility",
  "Permeability",
  "Binding Affinity",
  "Half-Life"
)

order_smiles <- c(
  "Hemolysis",
  "Non-Fouling",
  "Toxicity",
  "Binding Affinity",
  "Half-Life",
  "Permeability (PAMPA)",
  "Permeability (CACO2)"
)

# ============================================================
# Prepare Plot Data
# ============================================================

# Dual-axis mapping
rho_min <- 0
rho_max <- 1
rho_to_primary <- function(rho) (rho - rho_min) / (rho_max - rho_min)
primary_to_rho <- function(y) y * (rho_max - rho_min) + rho_min

plot_df <- plot_df %>%
  mutate(
    prop = as.character(prop),
    prop = case_when(
      panel == "Amino Acids Predictors" ~ factor(prop, levels = order_wt),
      panel == "SMILES Predictors" ~ factor(prop, levels = order_smiles),
      TRUE ~ factor(prop)
    ),
    y_plot = case_when(
      kind == "classifier" ~ value,
      kind == "regression" ~ rho_to_primary(value),
      TRUE ~ NA_real_
    )
  )

# Create x-axis key for faceting
x_levels <- c(
  paste("Amino Acids Predictors", order_wt, sep = "||"),
  paste("SMILES Predictors", order_smiles, sep = "||")
)

plot_df2 <- plot_df %>%
  mutate(
    x_key = paste(panel, as.character(prop), sep = "||"),
    x_key = factor(x_key, levels = x_levels),
    x_lab = as.character(prop)
  )

# ============================================================
# Create Plot
# ============================================================

PAL_KIND <- c(
  "classifier" = "#5dd8df",
  "regression" = "#dba5d6"
)

p_main <- ggplot(plot_df2, aes(x = x_key, y = y_plot)) +
  geom_col(aes(fill = kind, alpha = available), width = 0.75) +
  geom_text(
    data = subset(plot_df2, available),
    aes(label = model_lbl, size = ifelse(model_lbl == "Transformer", 4, 5)),
    vjust = -0.55
  ) +
  scale_size_identity() +
  geom_text(
    data = subset(plot_df2, !available),
    aes(y = 0.55, label = "N/A"),
    size = 3.2,
    alpha = 0.9
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  scale_x_discrete(labels = setNames(plot_df2$x_lab, plot_df2$x_key), drop = TRUE) +
  scale_fill_manual(values = PAL_KIND,
                    breaks = c("classifier", "regression"),
                    labels = c("Classifier (F1)", "Regression (Spearman ρ)")) +
  scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.25), guide = "none") +
  scale_y_continuous(
    limits = c(0, 1),
    name = "Best validation F1",
    labels = scales::number_format(accuracy = 0.01),
    sec.axis = sec_axis(~ primary_to_rho(.),
                        name = "Best validation Spearman ρ",
                        labels = scales::number_format(accuracy = 0.01))
  ) +
  labs(x = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

print(p_main)

# ============================================================
# Save Output
# ============================================================

out_png <- file.path(OUT_DIR, "main_best_model_dualaxis_f1_rho.png")
ggsave(out_png, p_main, width = 15, height = 6, dpi = 500)

out_pdf <- file.path(OUT_DIR, "main_best_model_dualaxis_f1_rho.pdf")
ggsave(out_pdf, p_main, width = 15, height = 6)

message("Saved:\n  ", out_png, "\n  ", out_pdf)