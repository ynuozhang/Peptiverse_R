ROOT <- rprojroot::find_rstudio_root_file()
source(file.path(ROOT, "R", "plot_style.R"))

suppressPackageStartupMessages({
  library(tidyverse)
})

WT_CSV     <- "iptm_scatter_csv_all/wt_iptm_affinity_all.csv"
SMILES_CSV <- "iptm_scatter_csv_all/smiles_iptm_affinity_all.csv"
OUT_DIR    <- "iptm_scatter_csv_all/plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_one <- function(path, branch) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(branch = branch) %>%
    mutate(
      affinity = as.numeric(affinity),
      iptm     = as.numeric(iptm),
      split    = if ("split" %in% names(.)) as.character(split) else "all"
    ) %>%
    filter(!is.na(affinity), !is.na(iptm))
}

df <- bind_rows(
  read_one(WT_CSV, "WT"),
  read_one(SMILES_CSV, "SMILES")
)

stats <- df %>%
  group_by(branch) %>%
  summarise(
    n = n(),
    pearson_r    = cor(iptm, affinity, method = "pearson"),
    spearman_rho = cor(iptm, affinity, method = "spearman"),
    r2_lm        = summary(lm(affinity ~ iptm, data = cur_data()))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "N=%d\nPearson r=%.3f\nSpearman ρ=%.3f\nLinear regression R²=%.3f",
      n, pearson_r, spearman_rho, r2_lm
    )
  )

# ---- plot (facet WT vs SMILES) ----
p <- ggplot(df, aes(x = iptm, y = affinity)) +
  geom_point(
    shape = 16,
    color = "#dba5d6",
    alpha = 0.70,
    size = 2.0
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#5dd8df",
    fill  = "#5dd8df",
    alpha = 0.18,        # subtle CI band
    linewidth = 0.9
  ) +
  facet_wrap(~ branch, scales = "free_y") +
  labs(
    x = "ipTM score",
    y = "Affinity"
  ) +
  theme(
    legend.position = "none"
  ) +
  geom_text(
    data = stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.05,
    size = 4.8
  )

print(p)
ggsave(file.path(OUT_DIR, "iptm_vs_affinity_facet.png"),
       p, width = 10, height = 4.2, dpi = 300)


for (b in unique(df$branch)) {
  d  <- df %>% filter(branch == b)
  st <- stats %>% filter(branch == b)
  
  pp <- ggplot(d, aes(x = iptm, y = affinity)) +
    geom_point(aes(shape = split), alpha = 0.55, size = 1.8, stroke = 0.2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "ipTM score",
      y = "Affinity"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = st$label,
             hjust = -0.05, vjust = 1.05, size = 3.4)
  
  ggsave(file.path(OUT_DIR, paste0("iptm_vs_affinity_", tolower(b), ".png")),
         pp, width = 5.2, height = 4.2, dpi = 300)
}

message("Wrote plots to: ", OUT_DIR)
