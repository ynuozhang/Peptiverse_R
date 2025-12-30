ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(scales)
  library(dplyr)
  library(grid)
})
# -----------------------------
# Regression traces: Spearman rho vs trial
# -----------------------------
IN_TRIALS <- "data_processed/trials_regression.csv"
OUT_DIR   <- "figures/supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

trials <- readr::read_csv(IN_TRIALS, show_col_types = FALSE) %>%
  mutate(
    task  = as.character(task),
    model = as.character(model),
    run   = as.character(run),
    repr  = as.character(repr),
    trial_index = as.integer(trial_index),
    rho = as.numeric(f1)   # <-- use f1 column as rho
  ) %>%
  mutate(
    repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr))
  )

# -----------------------------
# expects: data_processed/binding_affinity_trials_long.csv
# -----------------------------
BA_TRIALS <- "data_processed/binding_affinity_trials_long.csv"

if (file.exists(BA_TRIALS)) {
  ba <- readr::read_csv(BA_TRIALS, show_col_types = FALSE) %>%
    mutate(setting_l = tolower(setting)) %>%
    transmute(
      task = "binding_affinity",
      repr = case_when(
        str_detect(setting_l, "wt_smiles") ~ "smiles",
        str_detect(setting_l, "wt_wt")     ~ "wt",
        TRUE ~ NA_character_
      ),
      model = "transformer",
      run = coalesce(run, "binding_affinity"),
      trial_index = as.integer(trial_index),
      f1 = as.numeric(rho),
      best_so_far = as.numeric(best_so_far),
      auc = NA_real_, ap = NA_real_, threshold = NA_real_
    ) %>%
    filter(!is.na(repr))   # keep only wt/smiles
  
  
  trials <- bind_rows(trials, ba) %>%
    mutate(
      rho = as.numeric(f1),
      repr = case_when(is.na(repr) | repr == "" ~ "unspecified", TRUE ~ tolower(repr))
    )
  
  message("Appended binding affinity trials: ", nrow(ba))
} else {
  message("Binding affinity trials not found: ", BA_TRIALS)
}

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  x <- str_remove(x, "_log$")
  x <- str_remove(x, "_raw$")
  
  dplyr::case_when(
    str_detect(x, "svr")          ~ "SVR",
    str_detect(x, "^svm$|svm")    ~ "SVM",
    str_detect(x, "^xgb$|xgboost") ~ "XGB",
    str_detect(x, "enet|elastic") ~ "ENET",
    str_detect(x, "mlp")          ~ "MLP",
    str_detect(x, "cnn")          ~ "CNN",
    str_detect(x, "transformer|bert|esm") ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}

task_matches <- function(task_col, key) {
  task_col == key | stringr::str_detect(task_col, stringr::fixed(key, ignore_case = TRUE))
}

MODEL_ORDER <- c("XGB","ENET","SVR","MLP","CNN","Transformer")
COLOR_LEVELS <- c("XGB","ENET","SVM","MLP","CNN","Transformer")  # note SVM here only as a COLOR KEY

y_offset <- 0.20 
REG_TASK_KEYS <- c("permeabilitypampa", "permeabilitycaco2", "halflife", "bindingaffinity")
trials_reg <- trials %>%
  mutate(task_norm = str_replace_all(tolower(task), "[^a-z0-9]+", "")) %>%  # normalize
  filter(str_detect(task_norm, paste(REG_TASK_KEYS, collapse = "|"))) %>%
  mutate(
    task_label = case_when(
      str_detect(task_norm, "permeabilitypampa") ~ "Permeability_PAMPA",
      str_detect(task_norm, "permeabilitycaco2") ~ "Permeability_CACO2",
      str_detect(task_norm, "halflife")          ~ "Half-life",
      str_detect(task_norm, "bindingaffinity")   ~ "Binding_affinity",
      TRUE ~ task
    ),
    repr_label = case_when(
      repr == "wt" ~ "WT",
      repr == "smiles" ~ "SMILES",
      TRUE ~ repr
    ),
    panel = paste0(task_label, " (", repr_label, ")"),
    model_lbl = factor(clean_model_label(model), levels = c("XGB","ENET","SVM","SVR","MLP","CNN","Transformer"))
  )

# quick sanity check
message("Rows in trials_reg: ", nrow(trials_reg))
print(unique(trials_reg$task_norm))


best_pts <- trials_reg %>%
  filter(!is.na(rho)) %>%
  group_by(task, repr, model, run, panel, model_lbl) %>%
  arrange(desc(rho), trial_index) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    # draw a longer arrow so it doesn't look like a tiny triangle
    y_top = pmin(rho + 0.12, 0.995),   # start above
    y_bot = rho + 0.02                # end near the point
  )

targets <- tibble::tribble(
  ~panel_id, ~task_key,            ~repr_key, ~panel_label,
  "caco2_s", "permeabilitycaco2",  "smiles",  "Permeability_CACO2 (SMILES)",
  "pampa_s", "permeabilitypampa",  "smiles",  "Permeability_PAMPA (SMILES)",
  "hl_wt",   "halflife",           "wt",      "Half-Life (Amino Acids)",
  "ba_all",  "bindingaffinity",    "any",     "Binding Affinity (Amino Acids+SMILES)"
)


panel_dfs <- purrr::pmap_dfr(targets, function(panel_id, task_key, repr_key, panel_label) {
  
  d <- trials_reg %>%
    filter(str_detect(task_norm, fixed(task_key))) %>%   # <-- robust
    { if (repr_key == "any") . else filter(., tolower(repr) == repr_key) } %>%
    mutate(panel = panel_label)
  
  if (nrow(d) == 0) {
    message(sprintf("[Regression traces] Missing data for: %s -> placeholder facet", panel_label))
    tibble(
      task = task_key, repr = repr_key,
      model = "MISSING", run = "MISSING",
      trial_index = 0L, f1 = NA_real_, rho = NA_real_,
      best_so_far = NA_real_, auc = NA_real_, ap = NA_real_, threshold = NA_real_,
      panel = panel_label
    )
  } else {
    d
  }
})

panel_dfs <- panel_dfs %>%
  mutate(panel = factor(panel, levels = targets$panel_label))

trials_plot <- panel_dfs %>%
  mutate(
    model_lbl = factor(clean_model_label(model), levels = MODEL_ORDER),
    # map SVR to SVM color key so it shares SVM’s hex
    model_color = case_when(
      model_lbl == "SVR" ~ "SVM",
      TRUE ~ as.character(model_lbl)
    ),
    model_color = factor(model_color, levels = COLOR_LEVELS)
  )
trials_plot <- trials_plot %>%
  filter(

    # Binding affinity: keep Transformer only
    !(str_detect(str_replace_all(tolower(task), "[^a-z0-9]+", ""), "bindingaffinity") &
        model_lbl != "Transformer")
  )

best_pts2 <- trials_plot %>%
  filter(!is.na(rho), model != "MISSING") %>%
  group_by(task, repr, model, run, panel, model_lbl, model_color) %>%
  arrange(desc(rho), trial_index) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    x_mark = trial_index,
    y_mark = pmin(rho + y_offset, 0.98)  # long enough to see the shaft
  )

p_rho <- ggplot(trials_plot, aes(
  trial_index, rho,
  group = interaction(task, repr, model, run),
  color = model_color
)) +
  geom_line(linewidth = 0.6, alpha = 0.35, na.rm = TRUE) +
  
  geom_segment(
    data = best_pts2,
    aes(x = x_mark, xend = trial_index, y = y_mark, yend = rho, color = model_color),
    inherit.aes = FALSE,
    linewidth = 0.45,
    arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed"),
    show.legend = FALSE
  ) +
  
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  
  scale_color_manual(
    values = PAL_MODEL,                                  # uses SVM color for SVR via model_color
    breaks = c("XGB","ENET","SVM","MLP","CNN","Transformer"),
    labels = c("XGB","ENET","SVR","MLP","CNN","Transformer"), 
    drop = FALSE
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 1.8, alpha = 1, linetype = 1),
      nrow = 1
    )
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    x = "Optuna trial index",
    y = "Spearman ρ (val)"
  )


print(p_rho)


out_png <- file.path(OUT_DIR, "supp_S?_rho_vs_trial.png")
ggsave(out_png, p_rho, width = 10, height = 8, dpi = 300)

out_pdf <- file.path(OUT_DIR, "supp_S?_rho_vs_trial.pdf")
ggsave(out_pdf, p_rho, width = 20, height = 16)
