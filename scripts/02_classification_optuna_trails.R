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
IN_TRIALS <- "data_processed/trials_classification.csv"
OUT_DIR   <- "figures/supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  
  dplyr::case_when(
    str_detect(x, "^xgb$|xgboost") ~ "XGB",
    str_detect(x, "svm")          ~ "SVM",
    str_detect(x, "enet|elastic") ~ "ENET",
    str_detect(x, "mlp")          ~ "MLP",
    str_detect(x, "cnn")          ~ "CNN",
    str_detect(x, "transformer|bert|esm") ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}
# -----------------------------
# Load trials
# -----------------------------
trials <- readr::read_csv(IN_TRIALS, show_col_types = FALSE) %>%
  mutate(
    task  = as.character(task),
    model = as.character(model),
    run   = as.character(run),
    repr  = as.character(repr),
    trial_index = as.integer(trial_index),
    f1 = as.numeric(f1),
    best_so_far = as.numeric(best_so_far)
  )

# Normalize repr a bit (in case NA or "")
trials <- trials %>%
  mutate(
    repr = case_when(
      is.na(repr) | repr == "" ~ "unspecified",
      TRUE ~ tolower(repr)
    )
  )

y_offset <- 0.2 
trials_plot <- trials %>%
  mutate(
    task_label = task %>%
      str_replace_all("permeabilitypenetrance", "permeability_penetrance") %>% 
      str_replace_all("nonfouling", "non-fouling") %>% 
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_replace_all(" ", "_"),
    
    repr_label = case_when(
      repr == "wt" ~ "Amino Acids",
      repr == "smiles" ~ "SMILES",
      TRUE ~ repr
    ),
    
    panel = paste0(task_label, " (", repr_label, ")"),
    model_lbl = clean_model_label(model)
  )


best_pts <- trials_plot %>%
  group_by(task, repr, model, run) %>%
  arrange(desc(f1), trial_index) %>%
  slice(1) %>%
  ungroup()
best_pts2 <- best_pts %>%
  mutate(
    y_mark = pmin(f1 + y_offset, 0.995),   # float above, avoid clipping at 1.0
    x_mark = trial_index
  )
# best points df must also have model_lbl

targets <- tribble(
  ~panel_id, ~task_key,      ~repr_key,        ~panel_label,
  "sol",     "solubility",   "any",            "Solubility (all repr)",
  "hemo_wt", "hemolysis",    "wt",             "Hemolysis (Amino Acids)",
  "hemo_sm", "hemolysis",    "smiles",         "Hemolysis (SMILES)",
  "nf_wt",   "nonfouling",   "wt",             "Nonfouling (Amino Acids)",
  "nf_sm",   "nonfouling",   "smiles",         "Nonfouling (SMILES)",
  "tox_wt",   "toxicity",   "smiles",         "Toxicity (SMILES)",
  "perm_wt",   "permeabilitypenetrance",   "wt",         "Permeability_penetrance (WT)"
)
MODEL_ORDER <- c("XGB", "ENET", "SVM", "MLP", "CNN", "Transformer")
best_pts2 <- best_pts2 %>%
  mutate(model_lbl = clean_model_label(model))
trials_plot <- trials_plot %>%
  mutate(
    model_lbl = factor(
      clean_model_label(model),
      levels = MODEL_ORDER
    )
  )
best_pts2 <- best_pts2 %>%
  mutate(
    model_lbl = factor(
      clean_model_label(model),
      levels = MODEL_ORDER
    )
  )
task_matches <- function(task_col, key) {
  task_col == key | str_detect(task_col, fixed(key, ignore_case = TRUE))
}

panel_dfs <- purrr::pmap_dfr(targets, function(panel_id, task_key, repr_key, panel_label) {
  d <- trials %>%
    filter(task_matches(task, task_key)) %>%
    { if (repr_key == "any") . else filter(., repr == repr_key) } %>%
    mutate(panel = panel_label)
  
  if (nrow(d) == 0) {
    message(sprintf("[Supp S1] Missing data for: %s (task=%s, repr=%s) -> placeholder facet",
                    panel_label, task_key, repr_key))
    tibble(
      task = task_key,
      repr = repr_key,
      model = "MISSING",
      run = "MISSING",
      trial_index = 0L,
      f1 = NA_real_,
      best_so_far = NA_real_,
      auc = NA_real_,
      ap = NA_real_,
      threshold = NA_real_,
      panel = panel_label
    )
  } else {
    d
  }
})

panel_dfs <- panel_dfs %>%
  mutate(panel = factor(panel, levels = targets$panel_label))

# -----------------------------
# Plot: raw F1 (light) + best-so-far (bold)
# -----------------------------

p <- ggplot(trials_plot, aes(
  trial_index, f1,
  group = interaction(task, repr, model, run),
  color = model_lbl
)) +
  geom_line(linewidth = 0.6, alpha = 0.35) +
  
  geom_segment(
    data = best_pts2,
    aes(x = x_mark, xend = trial_index, y = y_mark, yend = f1, color = model_lbl),
    inherit.aes = FALSE,
    linewidth = 0.4,
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed"),
    show.legend = FALSE
  ) +
  
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  
  scale_color_manual(values = PAL_MODEL, drop = FALSE) +
  
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 1.8, alpha = 1, linetype = 1),
      nrow = 1
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    x = "Optuna trial index",
    y = "F1 (val)",
  )


print(p)

out_png <- file.path(OUT_DIR, "supp_S1_f1_vs_trial.png")
ggsave(out_png, p, width = 10, height =8, dpi = 300, device = ragg::agg_png)
message("Saved: ", out_png)

out_pdf <- file.path(OUT_DIR, "supp_S1_f1_vs_trial.pdf")
ggsave(out_pdf, p, width = 20, height = 16)
message("Saved: ", out_pdf)


