ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(scales)
  library(dplyr)
})

IN_TRIALS <- "data_processed/trials_long.csv"
OUT_DIR   <- "figures/main"
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

task_matches <- function(task_col, key) {
  task_col == key | str_detect(task_col, fixed(key, ignore_case = TRUE))
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
    f1 = suppressWarnings(as.numeric(f1))
  ) %>%
  mutate(
    repr = case_when(
      is.na(repr) | repr == "" ~ "unspecified",
      TRUE ~ tolower(repr)
    )
  )

# Detect Spearman rho column if present
rho_candidates <- c("spearman", "spearman_rho", "rho", "spearman_corr", "corr_spearman")
rho_col <- intersect(rho_candidates, names(trials))
if (length(rho_col) > 0) {
  trials <- trials %>% mutate(rho = suppressWarnings(as.numeric(.data[[rho_col[1]]])))
} else {
  trials <- trials %>% mutate(rho = NA_real_)
}

# -----------------------------
# Property spec
# kind = classifier => use F1
# kind = regression => use rho
# include_* controls whether prop exists at all in that panel
# placeholder_ok controls whether missing data show N/A
# -----------------------------
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

prop_order <- prop_spec$prop

# -----------------------------
# Helpers: best per run for a metric, then best model per task/repr
# -----------------------------
best_by_run_metric <- function(df, metric_col) {
  df %>%
    filter(!is.na(.data[[metric_col]])) %>%
    group_by(task, repr, model, run) %>%
    summarise(best_metric = max(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(model_lbl = clean_model_label(model))
}

best_model_per_task <- function(best_by_run_df) {
  best_by_run_df %>%
    group_by(task, repr, model_lbl) %>%
    summarise(
      mean_best = mean(best_metric),
      med_best  = median(best_metric),
      n_runs    = n_distinct(run),
      sd_best   = sd(best_metric),
      se_best   = sd_best / sqrt(n_runs),
      .groups = "drop"
    ) %>%
    arrange(task, repr, desc(mean_best), desc(med_best), desc(n_runs)) %>%
    group_by(task, repr) %>%
    slice(1) %>%
    ungroup()
}

# Classification summary (F1)
best_by_run_f1 <- best_by_run_metric(trials, "f1")
best_model_per_task_f1 <- best_model_per_task(best_by_run_f1)

# Regression summary (rho)
best_by_run_rho <- best_by_run_metric(trials, "rho")
best_model_per_task_rho <- best_model_per_task(best_by_run_rho)

# -----------------------------
# Build panel DF (drops disallowed props, keeps placeholders where allowed)
# -----------------------------
build_panel_df <- function(repr_key, panel_label) {
  spec2 <- prop_spec %>%
    mutate(include = if (repr_key == "wt") include_wt else include_smiles) %>%
    filter(include) %>%
    mutate(prop = factor(prop, levels = prop_order)) %>%
    arrange(prop)
  
  out <- purrr::pmap_dfr(
    spec2 %>% select(prop, task_key, kind, placeholder_ok),
    function(prop, task_key, kind, placeholder_ok) {
      
      d <- if (kind == "classifier") {
        best_model_per_task_f1 %>%
          filter(repr == repr_key) %>%
          filter(task_matches(task, task_key))
      } else {
        best_model_per_task_rho %>%
          filter(repr == repr_key) %>%
          filter(task_matches(task, task_key))
      }
      
      if (nrow(d) == 0) {
        if (!placeholder_ok) return(NULL)
        return(tibble(
          repr = repr_key,
          panel = panel_label,
          prop = factor(prop, levels = prop_order),
          kind = kind,
          available = FALSE,
          value = NA_real_,
          se = NA_real_,
          model_lbl = "Missing"
        ))
      }
      
      d <- d %>% arrange(desc(mean_best)) %>% slice(1)
      
      tibble(
        repr = repr_key,
        panel = panel_label,
        prop = factor(prop, levels = prop_order),
        kind = kind,
        available = TRUE,
        value = d$mean_best,
        se = d$se_best,
        model_lbl = d$model_lbl
      )
    }
  )
  
  out
}

df_wt     <- build_panel_df("wt",     "WT Predictors")
df_smiles <- build_panel_df("smiles", "SMILES Predictors")

plot_df <- bind_rows(df_wt, df_smiles) %>%
  mutate(panel = factor(panel, levels = c("WT Predictors", "SMILES Predictors"))) %>%
  group_by(panel) %>%
  mutate(prop = droplevels(prop)) %>%   # per-panel x categories
  ungroup()

# -----------------------------
# Dual-axis mapping
# Primary axis: F1 in [0,1]
# Secondary axis: rho (mapped linearly into [0,1] for plotting)
# -----------------------------
rho_vals <- plot_df %>%
  filter(kind == "regression", available, !is.na(value)) %>%
  pull(value)

rho_min <- if (length(rho_vals) > 0) min(rho_vals) else 0
rho_max <- if (length(rho_vals) > 0) max(rho_vals) else 1

if (isTRUE(all.equal(rho_min, rho_max))) {
  rho_min <- rho_min - 0.05
  rho_max <- rho_max + 0.05
}

rho_to_primary <- function(rho) (rho - rho_min) / (rho_max - rho_min)
primary_to_rho <- function(y) y * (rho_max - rho_min) + rho_min

plot_df <- plot_df %>%
  mutate(
    y_plot = case_when(
      kind == "classifier" ~ value,
      kind == "regression" ~ rho_to_primary(value),
      TRUE ~ NA_real_
    ),
    se_plot = case_when(
      kind == "classifier" ~ se,
      kind == "regression" ~ se / (rho_max - rho_min),
      TRUE ~ NA_real_
    )
  )

# -----------------------------
# Two colors only
# -----------------------------
PAL_KIND <- c(
  "classifier" = "#5dd8df",
  "regression" = "#dba5d6"
)

p_main <- ggplot(plot_df, aes(x = prop, y = y_plot)) +
  geom_col(aes(fill = kind, alpha = available), width = 0.75) +
  
  geom_errorbar(
    data = subset(plot_df, available & !is.na(se_plot) & !is.na(y_plot)),
    aes(ymin = pmax(y_plot - se_plot, 0), ymax = pmin(y_plot + se_plot, 1)),
    width = 0.18,
    linewidth = 0.4
  ) +
  
  geom_text(
    data = subset(plot_df, available),
    aes(label = model_lbl),
    vjust = -0.55,
    size = 5
  ) +
  
  geom_text(
    data = subset(plot_df, !available),
    aes(y = 0.55, label = "N/A"),
    size = 3.2,
    alpha = 0.9
  ) +
  
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  scale_x_discrete(drop = TRUE) +
  
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

out_png <- file.path(OUT_DIR, "main_best_model_dualaxis_f1_rho.png")
ggsave(out_png, p_main, width = 15, height = 6, dpi = 500)

out_pdf <- file.path(OUT_DIR, "main_best_model_dualaxis_f1_rho.pdf")
ggsave(out_pdf, p_main, width = 15, height = 6)

message("Saved:\n  ", out_png, "\n  ", out_pdf)
