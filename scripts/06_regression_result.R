ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(scales)
  library(dplyr)
})

RAW_DIR <- "raw"
OUT_DIR <- "figures/supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("parse_meta_from_path")) {
  stop("parse_meta_from_path() not found. Please source the script where you defined it, then rerun.")
}

read_csv_fast <- function(path) data.table::fread(path, data.table = FALSE, showProgress = FALSE)

norm_task <- function(x) str_replace_all(tolower(x), "[^a-z0-9]+", "")

task_to_prop <- function(task_norm) {
  dplyr::case_when(
    str_detect(task_norm, "halflife") ~ "Half-life",
    str_detect(task_norm, "bindingaffinity") ~ "Binding_affinity",
    str_detect(task_norm, "permeabilitypampa") ~ "Permeability_PAMPA",
    str_detect(task_norm, "permeabilitycaco2") ~ "Permeability_CACO2",
    TRUE ~ NA_character_
  )
}

clean_model_label <- function(x) {
  x <- tolower(x)
  x <- str_remove(x, "_gpu$")
  x <- str_remove(x, "_log$")
  x <- str_remove(x, "_raw$")
  x <- str_remove(x, "_unpooled$")
  x <- str_remove(x, "_pooled$")
  
  case_when(
    str_detect(x, "svr") ~ "SVR",
    str_detect(x, "svm") ~ "SVM",
    str_detect(x, "^xgb$|xgboost") ~ "XGB",
    str_detect(x, "enet|elastic") ~ "ENET",
    str_detect(x, "\\bmlp\\b") ~ "MLP",
    str_detect(x, "\\bcnn\\b") ~ "CNN",
    str_detect(x, "transformer|bert|esm") ~ "Transformer",
    TRUE ~ str_to_title(x)
  )
}

pick_first <- function(nm, cands) {
  hit <- intersect(cands, nm)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_cv_oof_spearman <- function(path) {
  x <- jsonlite::read_json(path, simplifyVector = TRUE)
  sp <- suppressWarnings(as.numeric(x[["spearman"]]))
  if (length(sp) == 0 || !is.finite(sp[1])) return(NA_real_)
  sp[1]
}

# Function to format p-values
format_pval <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

# -------------------------
# Regression metrics on a val_predictions.csv (y_true/y_pred)
# We compute: pearson r, spearman rho, r2(lm), rmse, mae WITH p-values
# -------------------------
compute_reg_metrics_from_val <- function(df, quiet = TRUE) {
  nm <- names(df)
  
  pick_first <- function(cands) {
    hit <- intersect(cands, nm)
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
  pred_col  <- pick_first(c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
  
  if (is.na(truth_col) || is.na(pred_col)) {
    stop("Regression val file missing y_true/y_pred cols. Found: ", paste(nm, collapse = ", "))
  }
  
  yt <- suppressWarnings(as.numeric(df[[truth_col]]))
  yp <- suppressWarnings(as.numeric(df[[pred_col]]))
  ok <- is.finite(yt) & is.finite(yp)
  yt <- yt[ok]; yp <- yp[ok]
  
  if (!quiet) message(sprintf("[reg-metrics] truth=%s pred=%s n=%d", truth_col, pred_col, length(yt)))
  
  # Pearson with p-value
  pearson_test <- suppressWarnings(cor.test(yt, yp, method = "pearson"))
  pear <- pearson_test$estimate
  pear_p <- pearson_test$p.value
  
  # Spearman with p-value
  spearman_test <- suppressWarnings(cor.test(yt, yp, method = "spearman"))
  rho <- spearman_test$estimate
  rho_p <- spearman_test$p.value
  
  # Linear regression with p-value
  lmfit <- tryCatch(lm(yt ~ yp), error = function(e) NULL)
  if (!is.null(lmfit)) {
    r2 <- summary(lmfit)$r.squared
    fstat <- summary(lmfit)$fstatistic
    lm_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else {
    r2 <- NA_real_
    lm_p <- NA_real_
  }
  
  rmse <- sqrt(mean((yt - yp)^2))
  mae  <- mean(abs(yt - yp))
  
  tibble(
    n = length(yt),
    pearson_r = as.numeric(pear),
    pearson_p = as.numeric(pear_p),
    spearman_rho = as.numeric(rho),
    spearman_p = as.numeric(rho_p),
    r2_lm = as.numeric(r2),
    lm_p = as.numeric(lm_p),
    rmse = as.numeric(rmse),
    mae  = as.numeric(mae)
  )
}

make_scatter <- function(df, stats_df, facet_var = "panel") {
  ggplot(df, aes(x = y_true, y = y_pred)) +
    geom_point(shape = 16, color = "#dba5d6", alpha = 0.70, size = 2.0) +
    geom_smooth(method = "lm", se = TRUE,
                color = "#5dd8df", fill = "#5dd8df",
                alpha = 0.18, linewidth = 0.9) +
    facet_wrap(as.formula(paste("~", facet_var)), scales = "free", ncol = 3) +
    labs(x = "y_true", y = "y_pred") +
    theme(legend.position = "none") +
    geom_text(
      data = stats_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.05,
      size = 4.8
    )
}

# =========================================================
# 1) BEST MODEL PER TASK×REPR from raw/**/val_predictions.csv
# =========================================================
val_files <- list.files(RAW_DIR, pattern = "val_predictions\\.csv$", recursive = TRUE, full.names = TRUE)

val_meta <- purrr::map_dfr(val_files, function(f) {
  meta <- parse_meta_from_path(f, RAW_DIR)
  meta %>% mutate(path = f)
}) %>%
  mutate(
    repr = case_when(
      !is.na(repr) ~ tolower(repr),
      str_detect(tolower(path), "(^|/|_)smiles($|/|_)") ~ "smiles",
      str_detect(tolower(path), "(^|/|_)wt($|/|_)") ~ "wt",
      TRUE ~ "unspecified"
    ),
    task_norm = norm_task(task),
    prop = task_to_prop(task_norm),
    repr = if_else(str_detect(task_norm, "permeabilitypampa|permeabilitycaco2"), "smiles", repr)
  ) %>%
  filter(!is.na(prop)) %>%
  filter(prop != "Binding_affinity")

val_metrics <- val_meta %>%
  mutate(
    m = purrr::map(path, ~ compute_reg_metrics_from_val(read_csv_fast(.x), quiet = TRUE))
  ) %>%
  unnest(m)

# Pick best run per prop×repr by max spearman rho
best_runs <- val_metrics %>%
  group_by(prop, repr) %>%
  arrange(desc(spearman_rho), desc(r2_lm), rmse) %>%
  slice(1) %>%
  ungroup()

message("[Best regressors picked]")
print(best_runs %>% select(prop, repr, task, model, run, spearman_rho, spearman_p, r2_lm, lm_p, rmse, mae, path))

best_scatter_df <- best_runs %>%
  mutate(
    df = purrr::map(path, function(p) {
      d <- read_csv_fast(p)
      nm <- names(d)
      
      pick_first <- function(cands) { hit <- intersect(cands, nm); if (length(hit)==0) NA_character_ else hit[1] }
      truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
      pred_col  <- pick_first(c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
      
      yt <- suppressWarnings(as.numeric(d[[truth_col]]))
      yp <- suppressWarnings(as.numeric(d[[pred_col]]))
      ok <- is.finite(yt) & is.finite(yp)
      
      tibble(
        y_true = yt[ok],
        y_pred = yp[ok]
      )
    }),
    repr_lab = case_when(
      repr == "wt" ~ "Amino acids",
      repr == "smiles" ~ "SMILES",
      TRUE ~ toupper(repr)
    ),
    panel = sprintf("%s (%s)\nBest: %s", prop, repr_lab, toupper(clean_model_label(model)))
  ) %>%
  select(prop, repr, model, run, panel, df) %>%
  unnest(df)

best_stats <- best_runs %>%
  mutate(
    repr_lab = case_when(
      repr == "wt" ~ "Amino acids",
      repr == "smiles" ~ "SMILES",
      TRUE ~ toupper(repr)
    ),
    panel = sprintf("%s (%s)\nBest: %s", prop, repr_lab, toupper(clean_model_label(model))),
    label = sprintf(
      "N=%d\nPearson r=%.3f (p=%s)\nSpearman ρ=%.3f (p=%s)\nLinear R²=%.3f (p=%s)",
      n, pearson_r, format_pval(pearson_p),
      spearman_rho, format_pval(spearman_p),
      r2_lm, format_pval(lm_p)
    )
  ) %>%
  select(panel, label)

# =========================================================
# 2) Binding affinity: separate WT and SMILES panels
# =========================================================
BA_VAL_SMI <- file.path("raw/binding_affinity", "wt_smiles_unpooled", "val_smiles_unpooled.csv")
BA_VAL_WT  <- file.path("raw/binding_affinity", "wt_wt_unpooled",     "val_wt_unpooled.csv")

read_ba_val <- function(path, panel_label) {
  dba <- readr::read_csv(path, show_col_types = FALSE)
  nm <- names(dba)
  
  pick_first <- function(cands) { hit <- intersect(cands, nm); if (length(hit)==0) NA_character_ else hit[1] }
  truth_col <- pick_first(c("y_true","truth","label","y","target","gt"))
  pred_col  <- pick_first(c("y_pred","y_pred_reg","pred","y_hat","prediction","pred_mean","estimate"))
  
  if (is.na(truth_col) || is.na(pred_col)) {
    stop("Binding affinity val file missing truth/pred cols: ", path,
         "\nFound: ", paste(nm, collapse = ", "))
  }
  
  yt <- suppressWarnings(as.numeric(dba[[truth_col]]))
  yp <- suppressWarnings(as.numeric(dba[[pred_col]]))
  ok <- is.finite(yt) & is.finite(yp)
  
  df <- tibble(y_true = yt[ok], y_pred = yp[ok], panel = panel_label)
  
  # Calculate correlations with p-values
  pearson_test <- suppressWarnings(cor.test(df$y_true, df$y_pred, method = "pearson"))
  pear <- pearson_test$estimate
  pear_p <- pearson_test$p.value
  
  spearman_test <- suppressWarnings(cor.test(df$y_true, df$y_pred, method = "spearman"))
  rho <- spearman_test$estimate
  rho_p <- spearman_test$p.value
  
  lmfit <- tryCatch(lm(y_pred ~ y_true, data = df), error = function(e) NULL)
  if (!is.null(lmfit)) {
    r2 <- summary(lmfit)$r.squared
    fstat <- summary(lmfit)$fstatistic
    lm_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else {
    r2 <- NA_real_
    lm_p <- NA_real_
  }
  
  st <- tibble(
    panel = panel_label,
    label = sprintf(
      "N=%d\nPearson r=%.3f (p=%s)\nSpearman ρ=%.3f (p=%s)\nLinear R²=%.3f (p=%s)",
      nrow(df), pear, format_pval(pear_p),
      rho, format_pval(rho_p),
      r2, format_pval(lm_p)
    )
  )
  
  list(df = df, st = st)
}

BA_PANEL_WT  <- "Binding Affinity (Amino acids)\nBest: TRANSFORMER on unpooled"
BA_PANEL_SMI <- "Binding Affinity (SMILES)\nBest: TRANSFORMER on unpooled"

if (file.exists(BA_VAL_WT)) {
  res <- read_ba_val(BA_VAL_WT,  BA_PANEL_WT)
  best_scatter_df <- bind_rows(best_scatter_df, res$df)
  best_stats      <- bind_rows(best_stats,      res$st)
}

if (file.exists(BA_VAL_SMI)) {
  res <- read_ba_val(BA_VAL_SMI, BA_PANEL_SMI)
  best_scatter_df <- bind_rows(best_scatter_df, res$df)
  best_stats      <- bind_rows(best_stats,      res$st)
}

# =========================================================
# 3) Half-life
# =========================================================
hl_wt_oof <- list.files(
  file.path(RAW_DIR, "half_life"),
  pattern = "^oof_predictions\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
hl_wt_oof <- hl_wt_oof[
  str_detect(tolower(dirname(hl_wt_oof)), "wt") &
    !str_detect(tolower(dirname(hl_wt_oof)), "smiles")
]

if (length(hl_wt_oof) > 0) {
  hl_wt_tbl <- purrr::map_dfr(hl_wt_oof, function(f) {
    meta <- parse_meta_from_path(f, RAW_DIR)
    
    d <- read_csv_fast(f)
    nm <- names(d)
    
    truth_col <- pick_first(nm, c("y_true","truth","label","y","target","gt"))
    pred_col  <- pick_first(nm, c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
    
    if (is.na(truth_col) || is.na(pred_col)) return(tibble())
    
    yt <- suppressWarnings(as.numeric(d[[truth_col]]))
    yp <- suppressWarnings(as.numeric(d[[pred_col]]))
    ok <- is.finite(yt) & is.finite(yp)
    if (sum(ok) < 3) return(tibble())
    
    sp <- suppressWarnings(cor(yt[ok], yp[ok], method = "spearman"))
    
    tibble(
      path = f,
      model = meta$model,
      run = meta$run,
      spearman = as.numeric(sp)
    )
  })
  
  if (nrow(hl_wt_tbl) > 0) {
    hl_wt_best <- hl_wt_tbl %>% arrange(desc(spearman)) %>% slice(1)
    best_m <- clean_model_label(hl_wt_best$model[1])
    
    dhl <- read_csv_fast(hl_wt_best$path[1])
    nm <- names(dhl)
    truth_col <- pick_first(nm, c("y_true","truth","label","y","target","gt"))
    pred_col  <- pick_first(nm, c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
    
    yt <- suppressWarnings(as.numeric(dhl[[truth_col]]))
    yp <- suppressWarnings(as.numeric(dhl[[pred_col]]))
    ok <- is.finite(yt) & is.finite(yp)
    
    df <- tibble(
      y_true = yt[ok],
      y_pred = yp[ok],
      panel  = sprintf("Half-life (Amino acids)\nBest: %s", toupper(best_m))
    )
    best_scatter_df <- bind_rows(best_scatter_df, df)
    
    # Calculate with p-values
    pearson_test <- suppressWarnings(cor.test(df$y_true, df$y_pred, method = "pearson"))
    pear <- pearson_test$estimate
    pear_p <- pearson_test$p.value
    
    spearman_test <- suppressWarnings(cor.test(df$y_true, df$y_pred, method = "spearman"))
    rho <- spearman_test$estimate
    rho_p <- spearman_test$p.value
    
    lmfit <- tryCatch(lm(y_pred ~ y_true, data = df), error = function(e) NULL)
    if (!is.null(lmfit)) {
      r2 <- summary(lmfit)$r.squared
      fstat <- summary(lmfit)$fstatistic
      lm_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    } else {
      r2 <- NA_real_
      lm_p <- NA_real_
    }
    
    best_stats <- bind_rows(
      best_stats,
      tibble(
        panel = unique(df$panel),
        label = sprintf(
          "N=%d\nPearson r=%.3f (p=%s)\nSpearman ρ=%.3f (p=%s)\nLinear R²=%.3f (p=%s)",
          nrow(df), pear, format_pval(pear_p),
          rho, format_pval(rho_p),
          r2, format_pval(lm_p)
        )
      )
    )
    
    message("Included half-life WT panel from: ", hl_wt_best$path[1])
  } else {
    message("No usable half-life WT oof_predictions.csv found.")
  }
} else {
  message("No half-life WT oof_predictions.csv found.")
}

# Half-life SMILES processing with p-values
hl_sum_json <- list.files(
  file.path(RAW_DIR, "half_life"),
  pattern = "cv_oof_summary\\.json$",
  recursive = TRUE,
  full.names = TRUE
)

hl_sum_json <- hl_sum_json[str_detect(tolower(hl_sum_json), "(^|/|_)smiles($|/|_)")]

if (length(hl_sum_json) > 0) {
  
  hl_tbl <- purrr::map_dfr(hl_sum_json, function(f) {
    meta <- parse_meta_from_path(f, RAW_DIR)
    if (is.na(meta$task) || meta$task != "halflife") return(tibble())
    
    is_smiles_run <- str_detect(tolower(f), "(^|/|_)smiles($|/|_)")
    if (!is_smiles_run) return(tibble())
    
    sp <- read_cv_oof_spearman(f)
    if (!is.finite(sp)) return(tibble())
    
    pred_path <- file.path(dirname(f), "cv_oof_predictions.csv")
    if (!file.exists(pred_path)) return(tibble())
    
    tibble(
      sum_path  = f,
      pred_path = pred_path,
      model     = meta$model,
      spearman  = sp
    )
  })
  
  if (nrow(hl_tbl) > 0) {
    hl_best <- hl_tbl %>% arrange(desc(spearman)) %>% slice(1)
    
    best_m <- clean_model_label(hl_best$model[1])
    panel_lbl <- sprintf("Half-life (SMILES)\nBest: %s", toupper(best_m))
    
    dhl <- read_csv_fast(hl_best$pred_path[1])
    nm <- names(dhl)
    
    truth_col <- pick_first(nm, c("y_true","truth","label","y","target","gt"))
    pred_col  <- pick_first(nm, c("y_pred","pred","y_hat","prediction","pred_mean","estimate"))
    
    if (!is.na(truth_col) && !is.na(pred_col)) {
      yt <- suppressWarnings(as.numeric(dhl[[truth_col]]))
      yp <- suppressWarnings(as.numeric(dhl[[pred_col]]))
      ok <- is.finite(yt) & is.finite(yp)
      
      df <- tibble(
        y_true = yt[ok],
        y_pred = yp[ok],
        panel  = panel_lbl
      )
      best_scatter_df <- bind_rows(best_scatter_df, df)
      
      # Calculate with p-values
      pearson_test <- suppressWarnings(cor.test(df$y_true, df$y_pred, method = "pearson"))
      pear <- pearson_test$estimate
      pear_p <- pearson_test$p.value
      
      spearman_test <- suppressWarnings(cor.test(df$y_true, df$y_pred, method = "spearman"))
      rho <- spearman_test$estimate
      rho_p <- spearman_test$p.value
      
      lmfit <- tryCatch(lm(y_pred ~ y_true, data = df), error = function(e) NULL)
      if (!is.null(lmfit)) {
        r2 <- summary(lmfit)$r.squared
        fstat <- summary(lmfit)$fstatistic
        lm_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      } else {
        r2 <- NA_real_
        lm_p <- NA_real_
      }
      
      best_stats <- bind_rows(
        best_stats,
        tibble(
          panel = panel_lbl,
          label = sprintf(
            "N=%d\nPearson r=%.3f (p=%s)\nSpearman ρ=%.3f (p=%s)\nLinear R²=%.3f (p=%s)",
            nrow(df), pear, format_pval(pear_p),
            rho, format_pval(rho_p),
            r2, format_pval(lm_p)
          )
        )
      )
      
      message("Included half-life SMILES panel from: ", hl_best$pred_path[1])
      message("Half-life SMILES best model: ", best_m, " (spearman=", round(hl_best$spearman[1], 4), ")")
    } else {
      message("Half-life SMILES cv_oof_predictions.csv missing truth/pred cols: ", hl_best$pred_path[1],
              "\nFound cols: ", paste(nm, collapse = ", "))
    }
  } else {
    message("No usable half-life SMILES found (skipping).")
  }
} else {
  message("No half-life SMILES cv_oof_summary.json found (skipping).")
}

# =========================================================
# 4) Plot + save
# =========================================================
panel_base <- function(p) stringr::str_replace(p, "\\nBest:.*$", "")

panel_key <- function(p) {
  p %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "")
}

PANEL_BASE <- c(
  "Permeability_PAMPA (SMILES)",
  "Half-Life (Amino acids)",
  "Binding Affinity (Amino acids)\nBest: TRANSFORMER on unpooled",
  "Permeability_CACO2 (SMILES)",
  "Half-Life (SMILES)",
  "Binding Affinity (SMILES)\nBest: TRANSFORMER on unpooled"
)

PANEL_BASE_BASE <- panel_base(PANEL_BASE)
PANEL_BASE_KEYS <- panel_key(PANEL_BASE_BASE)

panel_levels <- best_scatter_df %>%
  distinct(panel) %>%
  mutate(
    base = panel_base(panel),
    key  = panel_key(base),
    ord  = match(key, PANEL_BASE_KEYS),
    ord2 = dplyr::if_else(is.na(ord), 999L, ord)
  ) %>%
  arrange(ord2, panel) %>%
  pull(panel)

best_scatter_df <- best_scatter_df %>% mutate(panel = factor(panel, levels = panel_levels))
best_stats      <- best_stats      %>% mutate(panel = factor(panel, levels = panel_levels))

p <- make_scatter(best_scatter_df, best_stats, facet_var = "panel") +
  theme(
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )

print(p)

out_png <- file.path(OUT_DIR, "supp_best_regressors_val_scatter.png")
out_pdf <- file.path(OUT_DIR, "supp_best_regressors_val_scatter.pdf")

ggsave(out_png, p, width = 12, height = 3.8 * max(1, ceiling(nlevels(best_scatter_df$panel) / 2)), dpi = 400)
ggsave(out_pdf, p, width = 12, height = 3.8 * max(1, ceiling(nlevels(best_scatter_df$panel) / 2)))

message("Saved:\n  ", out_png, "\n  ", out_pdf)