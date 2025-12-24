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

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  
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

# choose which regression tasks to include (edit this vector)
REG_TASK_KEYS <- c("permeability_pampa", "permeability_caco2", "binding_affinity")

task_matches <- function(task_col, key) {
  task_col == key | stringr::str_detect(task_col, stringr::fixed(key, ignore_case = TRUE))
}

MODEL_ORDER <- c("XGB","ENET","SVR","MLP","CNN","Transformer")
COLOR_LEVELS <- c("XGB","ENET","SVM","MLP","CNN","Transformer")  # note SVM here only as a COLOR KEY

y_offset <- 0.20 
# regression tasks to include (normalized form)
REG_TASK_KEYS <- c("permeabilitypampa", "permeabilitycaco2")  # add bindingaffinity if needed

trials_reg <- trials %>%
  mutate(task_norm = str_replace_all(tolower(task), "[^a-z0-9]+", "")) %>%  # normalize
  filter(str_detect(task_norm, paste(REG_TASK_KEYS, collapse = "|"))) %>%
  mutate(
    task_label = case_when(
      str_detect(task_norm, "permeabilitypampa") ~ "Permeability_PAMPA",
      str_detect(task_norm, "permeabilitycaco2") ~ "Permeability_CACO2",
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
  "pampa_s", "permeabilitypampa",  "smiles",  "Permeability_PAMPA (SMILES)"
)

panel_dfs <- purrr::pmap_dfr(targets, function(panel_id, task_key, repr_key, panel_label) {
  
  d <- trials %>%
    mutate(task_norm = str_replace_all(tolower(task), "[^a-z0-9]+", "")) %>%
    filter(task_norm == task_key) %>%
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
    labels = c("XGB","ENET","SVR","MLP","CNN","Transformer"),  # <-- rename legend entry
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
ggsave(out_png, p_rho, width = 10, height = 8, dpi = 300)   # recommend 300

out_pdf <- file.path(OUT_DIR, "supp_S?_rho_vs_trial.pdf")
ggsave(out_pdf, p_rho, width = 20, height = 16)
