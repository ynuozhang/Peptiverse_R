ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(scales)
  library(dplyr)
  library(forcats)
})

# -----------------------------
# Binding affinity override:
# pick best across pooled/unpooled separately for wt and smiles
# -----------------------------
BA_CSV <- "data_processed/binding_affinity_trials_long.csv"

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
    mean_best = max(rho, na.rm = TRUE),  # best across pooled/unpooled
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


OUT_DIR   <- "figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
IN_TRIALS_CLS <- "data_processed/trials_classification.csv"
IN_TRIALS_REG <- "data_processed/trials_regression.csv"
#IN_TRIALS_REG <- "data_processed/trials_regression_plus_binding.csv"

pick_obj_col <- function(df) {
  # prefer f1; if absent, fall back to best_so_far; else NA
  if ("f1" %in% names(df)) return("f1")
  if ("best_so_far" %in% names(df)) return("best_so_far")
  return(NA_character_)
}

rho_candidates <- c("spearman", "spearman_rho", "rho", "spearman_corr", "corr_spearman")

read_trials_one <- function(path, kind_tag) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  
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
trials <- bind_rows(trials_cls, trials_reg)
# sanity checks (should NOT error)
print(names(trials))
print(dplyr::count(trials, kind, model, sort = TRUE))
stopifnot("value" %in% names(trials), "kind" %in% names(trials))

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  
  dplyr::case_when(
    str_detect(x, "svr")          ~ "SVR",
    str_detect(x, "svm")          ~ "SVM",
    str_detect(x, "^xgb$|xgboost") ~ "XGB",
    str_detect(x, "enet|elastic") ~ "ENET",
    str_detect(x, "\\bmlp\\b") ~ "MLP",
    str_detect(x, "\\bcnn\\b") ~ "CNN",
    str_detect(x, "transformer|bert|esm") ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}



task_matches <- function(task_col, key) {
  norm <- function(x) {
    x %>% tolower() %>% str_replace_all("[^a-z0-9]+", "")
  }
  str_detect(norm(task_col), fixed(norm(key)))
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
  "Toxicity",                 "toxicity",                "classifier",    FALSE,        TRUE,            TRUE,
  "Solubility",               "solubility",              "classifier",    TRUE,        FALSE,           FALSE,
  "Permeability",             "permeabilitypenetrance",            "classifier",    TRUE,        FALSE,           TRUE,
  "Halflife",                 "halflife",                "regression",    TRUE,        TRUE,            TRUE,
  "Permeability (PAMPA)",     "permeability_pampa",      "regression",    FALSE,       TRUE,            TRUE,
  "Permeability (CACO2)",     "permeability_caco2",      "regression",    FALSE,       TRUE,            TRUE,
  "Binding affinity",         "binding",                 "regression",    TRUE,        TRUE,            TRUE
)

prop_order <- prop_spec$prop

# -----------------------------
# Helpers: best per run for a metric, then best model per task/repr
# -----------------------------
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

best_by_run_all <- best_by_run_metric(trials)
best_model_per_task_all %>%
  filter(repr == "wt", kind == "classifier") %>%
  filter(task_matches(task, "permeabilitypenetrance")) %>%
  select(task, repr, kind, model_lbl, mean_best, n_runs) %>%
  print(n = 50)

best_model_per_task_all <- best_model_per_task(best_by_run_all)

# -----------------------------
# Build panel DF
# -----------------------------
build_panel_df <- function(repr_key, panel_label) {
  spec2 <- prop_spec %>%
    mutate(include = if (repr_key == "wt") include_wt else include_smiles) %>%
    filter(include) %>%
    mutate(prop = as.character(prop)) %>%
    arrange(prop)
  
  out <- purrr::pmap_dfr(
    spec2 %>% select(prop, task_key, kind, placeholder_ok),
    function(prop, task_key, kind_spec, placeholder_ok) {
      
      # ---- SPECIAL CASE: Binding affinity comes from ba_best, not best_model_per_task_all ----
      if (task_matches("binding", task_key)) {
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
      
      
      # ---- DEFAULT: everything else uses best_model_per_task_all ----
      d <- best_model_per_task_all %>%
        filter(repr == repr_key) %>%
        filter(kind == kind_spec) %>%
        filter(task_matches(task, task_key))
      
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

df_wt     <- build_panel_df("wt",     "WT Predictors")
df_smiles <- build_panel_df("smiles", "SMILES Predictors")

plot_df <- bind_rows(df_wt, df_smiles) %>%
  mutate(panel = factor(panel, levels = c("WT Predictors", "SMILES Predictors")))

# -----------------------------
# Panel-specific x ordering
# -----------------------------
order_wt <- c(
  "Hemolysis",
  "Nonfouling",
  "Solubility",
  "Permeability",
  "Binding affinity",
  "Halflife"
)

order_smiles <- c(
  "Hemolysis",
  "Nonfouling",
  "Toxicity",
  "Binding affinity",
  "Halflife",
  "Permeability (PAMPA)",
  "Permeability (CACO2)"
)

plot_df <- plot_df %>%
  mutate(
    prop_chr = as.character(prop),
    prop = case_when(
      panel == "WT Predictors" ~ factor(prop_chr, levels = order_wt),
      panel == "SMILES Predictors" ~ factor(prop_chr, levels = order_smiles),
      TRUE ~ factor(prop_chr)
    )
  ) %>%
  select(-prop_chr)
# -----------------------------
# Dual-axis mapping
# Primary axis: F1 in [0,1]
# Secondary axis: rho (mapped linearly into [0,1] for plotting)
# -----------------------------
rho_vals <- plot_df %>%
  filter(kind == "regression", available, !is.na(value)) %>%
  pull(value)

rho_min <- 0
rho_max <- 1

rho_to_primary <- function(rho) (rho - rho_min) / (rho_max - rho_min)  # == rho
primary_to_rho <- function(y) y * (rho_max - rho_min) + rho_min        # == y


plot_df <- plot_df %>%
  mutate(
    y_plot = case_when(
      kind == "classifier" ~ value,
      kind == "regression" ~ rho_to_primary(value),
      TRUE ~ NA_real_
    )
  )
# -----------------------------
# Panel-specific x ordering
# -----------------------------
order_wt <- c(
  "Hemolysis",
  "Nonfouling",
  "Solubility",
  "Permeability",
  "Binding affinity",
  "Halflife"
)

order_smiles <- c(
  "Hemolysis",
  "Nonfouling",
  "Toxicity",
  "Binding affinity",
  "Halflife",
  "Permeability (PAMPA)",
  "Permeability (CACO2)"
)

plot_df <- plot_df %>%
  mutate(
    prop = as.character(prop),
    prop = case_when(
      panel == "WT Predictors"     ~ factor(prop, levels = order_wt),
      panel == "SMILES Predictors" ~ factor(prop, levels = order_smiles),
      TRUE ~ factor(prop)
    )
  )

# -----------------------------
# Two colors only
# -----------------------------
PAL_KIND <- c(
  "classifier" = "#5dd8df",
  "regression" = "#dba5d6"
)

x_levels <- c(
  paste("WT Predictors", order_wt, sep = "||"),
  paste("SMILES Predictors", order_smiles, sep = "||")
)

plot_df2 <- plot_df %>%
  mutate(
    x_key = paste(panel, as.character(prop), sep = "||"),
    x_key = factor(x_key, levels = x_levels),
    x_lab = as.character(prop)
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

out_png <- file.path(OUT_DIR, "main_best_model_dualaxis_f1_rho.png")
ggsave(out_png, p_main, width = 15, height = 6, dpi = 500)

out_pdf <- file.path(OUT_DIR, "main_best_model_dualaxis_f1_rho.pdf")
ggsave(out_pdf, p_main, width = 15, height = 6)

message("Saved:\n  ", out_png, "\n  ", out_pdf)
