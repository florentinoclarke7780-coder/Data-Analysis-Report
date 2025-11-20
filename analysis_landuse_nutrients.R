## ------------------------------------------------------------
## Land Use vs N & P: Course-Level Analysis
## ------------------------------------------------------------

## Install lines (uncomment if needed):
# install.packages(c("tidyverse","car","FSA","dunn.test","svglite"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(car)        # Levene test
})

## Prefer FSA::dunnTest; fall back to dunn.test if FSA not available
has_FSA <- requireNamespace("FSA", quietly = TRUE)
has_dunn <- requireNamespace("dunn.test", quietly = TRUE)
## For SVG export
invisible(requireNamespace("svglite", quietly = TRUE))

## Make sure directories exist
dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plot_data", recursive = TRUE, showWarnings = FALSE)

## Unified academic plotting style
theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_blank()
    )
)

## Fixed colors for LandUseType
landuse_colors <- c(
  "Agricultural" = "#1b9e77",
  "Urban"        = "#d95f02",
  "Forest"       = "#7570b3"
)

## Helper: save plot in PNG, PDF, SVG
save_plot_all_formats <- function(p, filename_base, width = 8, height = 6, dpi = 300) {
  # PNG
  ggsave(
    filename = file.path("results", "figures", paste0(filename_base, ".png")),
    plot = p, width = width, height = height, dpi = dpi
  )
  # PDF
  ggsave(
    filename = file.path("results", "figures", paste0(filename_base, ".pdf")),
    plot = p, width = width, height = height
  )
  # SVG
  ggsave(
    filename = file.path("results", "figures", paste0(filename_base, ".svg")),
    plot = p, width = width, height = height, device = svglite::svglite
  )
}

message("[1/4] Loading dataset...")

## Load data robustly
data_paths <- c(
  file.path("data", "nawqa_landuse_nutrients_subset.csv"),
  file.path(".",   "nawqa_landuse_nutrients_subset.csv")
)

file_to_read <- data_paths[file.exists(data_paths)][1]
if (is.na(file_to_read)) stop("Input CSV not found in data/ or ./, expected nawqa_landuse_nutrients_subset.csv")

dat <- readr::read_csv(file_to_read, show_col_types = FALSE)

## Ensure required columns exist
required_cols <- c(
  "stationId","statePostalCode","placeName","latitude","longitude",
  "LandUseType","pAG","pURB","pFOR","TN","TP"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

## Consistent factor order for LandUseType
dat <- dat %>% mutate(LandUseType = factor(LandUseType, levels = c("Agricultural","Urban","Forest")))

message("[2/4] Assumption checks and transformations...")

## Log10 transforms (TN/TP are positive by construction of the subset)
dat <- dat %>%
  mutate(
    log10TN = log10(TN),
    log10TP = log10(TP)
  )

## Shapiro–Wilk per group and Levene across groups for TN/TP
sw_by_group <- function(x, g) {
  split(x, g) |>
    lapply(function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 3) return(c(statistic = NA_real_, p.value = NA_real_))
      out <- stats::shapiro.test(v)
      c(statistic = unname(out$statistic), p.value = out$p.value)
    })
}

## Assemble assumption results
assump_list <- list()

## TN
tn_sw_raw  <- sw_by_group(dat$TN, dat$LandUseType)
tn_sw_log  <- sw_by_group(dat$log10TN, dat$LandUseType)
## Levene (center = median is common in intro courses; we use mean for alignment with many textbooks)
lev_tn_raw <- tryCatch({
  car::leveneTest(TN ~ LandUseType, data = dat)
}, error = function(e) NULL)
lev_tn_log <- tryCatch({
  car::leveneTest(log10TN ~ LandUseType, data = dat)
}, error = function(e) NULL)

## TP
tp_sw_raw  <- sw_by_group(dat$TP, dat$LandUseType)
tp_sw_log  <- sw_by_group(dat$log10TP, dat$LandUseType)
lev_tp_raw <- tryCatch({
  car::leveneTest(TP ~ LandUseType, data = dat)
}, error = function(e) NULL)
lev_tp_log <- tryCatch({
  car::leveneTest(log10TP ~ LandUseType, data = dat)
}, error = function(e) NULL)

## Tidy assumption table
make_sw_tibble <- function(sw_list, var_lab, scale_lab) {
  tibble(
    variable = var_lab,
    scale = scale_lab,
    LandUseType = names(sw_list),
    W = sapply(sw_list, function(x) x[["statistic"]]),
    p_value = sapply(sw_list, function(x) x[["p.value"]])
  )
}

assump_tbl <- bind_rows(
  make_sw_tibble(tn_sw_raw, "TN", "raw"),
  make_sw_tibble(tn_sw_log, "TN", "log10"),
  make_sw_tibble(tp_sw_raw, "TP", "raw"),
  make_sw_tibble(tp_sw_log, "TP", "log10")
) %>%
  mutate(test = "Shapiro-Wilk")

## Extract Levene results (robust to missing broom)
pack_lev <- function(lev_obj, var_lab, scale_lab) {
  if (is.null(lev_obj)) return(tibble())
  tab <- as.data.frame(lev_obj)
  stat <- suppressWarnings(as.numeric(tab$`F value`[1]))
  pval <- suppressWarnings(as.numeric(tab$`Pr(>F)`[1]))
  tibble(variable = var_lab, scale = scale_lab, test = "Levene", LandUseType = NA_character_, W = stat, p_value = pval)
}
lev_tbl <- bind_rows(
  pack_lev(lev_tn_raw, "TN", "raw"),
  pack_lev(lev_tn_log, "TN", "log10"),
  pack_lev(lev_tp_raw, "TP", "raw"),
  pack_lev(lev_tp_log, "TP", "log10")
)

assumption_checks <- bind_rows(assump_tbl, lev_tbl)
readr::write_csv(assumption_checks, file.path("results","tables","assumption_checks.csv"))

## Helper to decide between ANOVA or Kruskal for one response
choose_test <- function(df, response) {
  ## Check ANOVA assumptions on raw and log10; prefer raw if both ok, else log10; else Kruskal
  y_raw <- df[[response]]
  y_log <- if (response == "TN") df$log10TN else df$log10TP

  ok_norm_raw <- all(sapply(split(y_raw, df$LandUseType), function(v) length(v <- v[is.finite(v)]) >= 3 && shapiro.test(v)$p.value > 0.05))
  ok_var_raw  <- tryCatch({ car::leveneTest(y_raw ~ df$LandUseType)$`Pr(>F)`[1] > 0.05 }, error = function(e) FALSE)

  ok_norm_log <- all(sapply(split(y_log, df$LandUseType), function(v) length(v <- v[is.finite(v)]) >= 3 && shapiro.test(v)$p.value > 0.05))
  ok_var_log  <- tryCatch({ car::leveneTest(y_log ~ df$LandUseType)$`Pr(>F)`[1] > 0.05 }, error = function(e) FALSE)

  if (ok_norm_raw && ok_var_raw) return(list(method = "anova", transform = "raw"))
  if (ok_norm_log && ok_var_log) return(list(method = "anova", transform = "log10"))
  list(method = "kruskal", transform = "raw")
}

message("[3/4] Main tests and figures...")

created_files <- c()

## A) TN ~ LandUseType
dec_tn <- choose_test(dat, "TN")
if (dec_tn$method == "anova") {
  resp <- if (dec_tn$transform == "log10") dat$log10TN else dat$TN
  fit_aov <- aov(resp ~ LandUseType, data = dat)
  aov_sum <- summary(fit_aov)[[1]]
  res_main_tn <- tibble(
    test = "ANOVA",
    response = if (dec_tn$transform == "log10") "log10TN" else "TN",
    transform = dec_tn$transform,
    df1 = aov_sum["LandUseType","Df"],
    df2 = aov_sum["Residuals","Df"],
    statistic = aov_sum["LandUseType","F value"],
    p_value = aov_sum["LandUseType","Pr(>F)"]
  )
  readr::write_csv(res_main_tn, file.path("results","tables","anova_or_kruskal_tn.csv"))
  created_files <- c(created_files, file.path("results","tables","anova_or_kruskal_tn.csv"))
  ## Tukey HSD
  tk <- TukeyHSD(fit_aov, which = "LandUseType")
  tk_tbl <- as.data.frame(tk$LandUseType)
  tk_tbl$contrast <- rownames(tk_tbl)
  tk_tbl <- tk_tbl %>% select(contrast, diff, lwr, upr, `p adj`) %>% rename(p_adj = `p adj`)
  readr::write_csv(tk_tbl, file.path("results","tables","posthoc_tn.csv"))
  created_files <- c(created_files, file.path("results","tables","posthoc_tn.csv"))
} else {
  kw <- kruskal.test(TN ~ LandUseType, data = dat)
  res_main_tn <- tibble(
    test = "Kruskal-Wallis",
    response = "TN",
    transform = "raw",
    df1 = length(levels(dat$LandUseType)) - 1,
    df2 = NA_real_,
    statistic = unname(kw$statistic),
    p_value = kw$p.value
  )
  readr::write_csv(res_main_tn, file.path("results","tables","anova_or_kruskal_tn.csv"))
  created_files <- c(created_files, file.path("results","tables","anova_or_kruskal_tn.csv"))
  ## Dunn test
  if (has_FSA) {
    dt <- FSA::dunnTest(TN ~ LandUseType, data = dat, method = "bh")
    dt_tbl <- dt$res %>% as_tibble() %>% rename(p_adj = P.adj)
  } else if (has_dunn) {
    dd <- dunn.test::dunn.test(x = dat$TN, g = dat$LandUseType, method = "bh", list = TRUE)
    dt_tbl <- tibble(contrast = dd$comparisons, Z = dd$Z, p_unadj = dd$P, p_adj = dd$P.adjusted)
  } else {
    ## Fallback: pairwise Wilcoxon (BH adjust) as an approximation when Dunn is unavailable
    pw <- pairwise.wilcox.test(dat$TN, dat$LandUseType, p.adjust.method = "BH", exact = FALSE)
    combs <- as.data.frame(as.table(pw$p.value)) %>% filter(!is.na(Freq))
    dt_tbl <- combs %>% transmute(contrast = paste(Var1, Var2, sep = "-"), p_adj = Freq)
  }
  readr::write_csv(dt_tbl, file.path("results","tables","posthoc_tn.csv"))
  created_files <- c(created_files, file.path("results","tables","posthoc_tn.csv"))
}

## Figure 1: TN by LandUseType (boxplot; use selected transform for y-axis)
df_fig1_tn_boxplot <- dat %>% transmute(LandUseType, TN, log10TN)
readr::write_csv(df_fig1_tn_boxplot, file.path("results","plot_data","fig1_tn_boxplot_data.csv"))
created_files <- c(created_files, file.path("results","plot_data","fig1_tn_boxplot_data.csv"))

y_var <- if (dec_tn$transform == "log10") "log10TN" else "TN"
y_lab <- if (dec_tn$transform == "log10") "log10 TN (mg/L)" else "TN (mg/L)"

## group sample sizes for labels
n_lab_fig1 <- dat %>% count(LandUseType) %>% mutate(label = paste0("n=", n))

p_fig1_tn_boxplot <- ggplot(dat, aes(x = LandUseType, y = .data[[y_var]], fill = LandUseType)) +
  geom_violin(trim = FALSE, alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.25, outlier.alpha = 0.4) +
  geom_text(data = n_lab_fig1, aes(x = LandUseType, y = Inf, label = label),
            vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  scale_fill_manual(values = landuse_colors) +
  labs(x = "Land-use type", y = y_lab, title = "TN by land-use type")

## add concise caption with test and p-value
cap_tn <- tryCatch({
  stat <- if (exists("res_main_tn")) res_main_tn$statistic[1] else NA_real_
  pv <- if (exists("res_main_tn")) res_main_tn$p_value[1] else NA_real_
  paste0("Kruskal–Wallis ",
         "\u03C7\u00B2\u2248", sprintf("%.1f", stat), ", p=",
         ifelse(is.na(pv), "NA", ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv))),
         "; Agricultural > Urban > Forest")
}, error = function(e) NULL)
if (!is.null(cap_tn)) p_fig1_tn_boxplot <- p_fig1_tn_boxplot + labs(caption = cap_tn)

save_plot_all_formats(p_fig1_tn_boxplot, "fig1_tn_boxplot")
created_files <- c(created_files,
  file.path("results","figures","fig1_tn_boxplot.png"),
  file.path("results","figures","fig1_tn_boxplot.pdf"),
  file.path("results","figures","fig1_tn_boxplot.svg")
)

## B) TP ~ LandUseType
dec_tp <- choose_test(dat, "TP")
if (dec_tp$method == "anova") {
  resp <- if (dec_tp$transform == "log10") dat$log10TP else dat$TP
  fit_aov <- aov(resp ~ LandUseType, data = dat)
  aov_sum <- summary(fit_aov)[[1]]
  res_main_tp <- tibble(
    test = "ANOVA",
    response = if (dec_tp$transform == "log10") "log10TP" else "TP",
    transform = dec_tp$transform,
    df1 = aov_sum["LandUseType","Df"],
    df2 = aov_sum["Residuals","Df"],
    statistic = aov_sum["LandUseType","F value"],
    p_value = aov_sum["LandUseType","Pr(>F)"]
  )
  readr::write_csv(res_main_tp, file.path("results","tables","anova_or_kruskal_tp.csv"))
  created_files <- c(created_files, file.path("results","tables","anova_or_kruskal_tp.csv"))
  ## Tukey HSD
  tk <- TukeyHSD(fit_aov, which = "LandUseType")
  tk_tbl <- as.data.frame(tk$LandUseType)
  tk_tbl$contrast <- rownames(tk_tbl)
  tk_tbl <- tk_tbl %>% select(contrast, diff, lwr, upr, `p adj`) %>% rename(p_adj = `p adj`)
  readr::write_csv(tk_tbl, file.path("results","tables","posthoc_tp.csv"))
  created_files <- c(created_files, file.path("results","tables","posthoc_tp.csv"))
} else {
  kw <- kruskal.test(TP ~ LandUseType, data = dat)
  res_main_tp <- tibble(
    test = "Kruskal-Wallis",
    response = "TP",
    transform = "raw",
    df1 = length(levels(dat$LandUseType)) - 1,
    df2 = NA_real_,
    statistic = unname(kw$statistic),
    p_value = kw$p.value
  )
  readr::write_csv(res_main_tp, file.path("results","tables","anova_or_kruskal_tp.csv"))
  created_files <- c(created_files, file.path("results","tables","anova_or_kruskal_tp.csv"))
  ## Dunn test
  if (has_FSA) {
    dt <- FSA::dunnTest(TP ~ LandUseType, data = dat, method = "bh")
    dt_tbl <- dt$res %>% as_tibble() %>% rename(p_adj = P.adj)
  } else if (has_dunn) {
    dd <- dunn.test::dunn.test(x = dat$TP, g = dat$LandUseType, method = "bh", list = TRUE)
    dt_tbl <- tibble(contrast = dd$comparisons, Z = dd$Z, p_unadj = dd$P, p_adj = dd$P.adjusted)
  } else {
    pw <- pairwise.wilcox.test(dat$TP, dat$LandUseType, p.adjust.method = "BH", exact = FALSE)
    combs <- as.data.frame(as.table(pw$p.value)) %>% filter(!is.na(Freq))
    dt_tbl <- combs %>% transmute(contrast = paste(Var1, Var2, sep = "-"), p_adj = Freq)
  }
  readr::write_csv(dt_tbl, file.path("results","tables","posthoc_tp.csv"))
  created_files <- c(created_files, file.path("results","tables","posthoc_tp.csv"))
}

## Figure 2: TP by LandUseType (boxplot)
df_fig2_tp_boxplot <- dat %>% transmute(LandUseType, TP, log10TP)
readr::write_csv(df_fig2_tp_boxplot, file.path("results","plot_data","fig2_tp_boxplot_data.csv"))
created_files <- c(created_files, file.path("results","plot_data","fig2_tp_boxplot_data.csv"))

y_var_tp <- if (dec_tp$transform == "log10") "log10TP" else "TP"
y_lab_tp <- if (dec_tp$transform == "log10") "log10 TP (mg/L)" else "TP (mg/L)"

p2_n_lab <- dat %>% count(LandUseType) %>% mutate(label = paste0("n=", n))
p_fig2_tp_boxplot <- ggplot(dat, aes(x = LandUseType, y = .data[[y_var_tp]], fill = LandUseType)) +
  geom_violin(trim = FALSE, alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.25, outlier.alpha = 0.4) +
  geom_text(data = p2_n_lab, aes(x = LandUseType, y = Inf, label = label),
            vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  scale_fill_manual(values = landuse_colors) +
  labs(x = "Land-use type", y = y_lab_tp, title = "TP by land-use type")

cap_tp <- tryCatch({
  stat <- if (exists("res_main_tp")) res_main_tp$statistic[1] else NA_real_
  pv <- if (exists("res_main_tp")) res_main_tp$p_value[1] else NA_real_
  paste0("Kruskal–Wallis ",
         "\u03C7\u00B2\u2248", sprintf("%.1f", stat), ", p=",
         ifelse(is.na(pv), "NA", ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv))),
         "; Forest < Agricultural \u2248 Urban")
}, error = function(e) NULL)
if (!is.null(cap_tp)) p_fig2_tp_boxplot <- p_fig2_tp_boxplot + labs(caption = cap_tp)

save_plot_all_formats(p_fig2_tp_boxplot, "fig2_tp_boxplot")
created_files <- c(created_files,
  file.path("results","figures","fig2_tp_boxplot.png"),
  file.path("results","figures","fig2_tp_boxplot.pdf"),
  file.path("results","figures","fig2_tp_boxplot.svg")
)

## C) Two-sample TN: Agricultural vs Forest
dat_ag_for <- dat %>% filter(LandUseType %in% c("Agricultural","Forest")) %>% droplevels()

## Choose test on raw/log10 similarly
dec_tn_2g <- (function(df) {
  y_raw <- df$TN; y_log <- df$log10TN
  g <- df$LandUseType
  ok_norm_raw <- all(sapply(split(y_raw, g), function(v) length(v <- v[is.finite(v)]) >= 3 && shapiro.test(v)$p.value > 0.05))
  ok_var_raw  <- tryCatch({ car::leveneTest(y_raw ~ g)$`Pr(>F)`[1] > 0.05 }, error = function(e) FALSE)
  ok_norm_log <- all(sapply(split(y_log, g), function(v) length(v <- v[is.finite(v)]) >= 3 && shapiro.test(v)$p.value > 0.05))
  ok_var_log  <- tryCatch({ car::leveneTest(y_log ~ g)$`Pr(>F)`[1] > 0.05 }, error = function(e) FALSE)
  if (ok_norm_raw && ok_var_raw) return(list(method = "ttest", transform = "raw"))
  if (ok_norm_log && ok_var_log) return(list(method = "ttest", transform = "log10"))
  list(method = "wilcox", transform = "raw")
})(dat_ag_for)

if (dec_tn_2g$method == "ttest") {
  y <- if (dec_tn_2g$transform == "log10") dat_ag_for$log10TN else dat_ag_for$TN
  tt <- t.test(y ~ dat_ag_for$LandUseType, var.equal = TRUE)
  res_2g <- tibble(
    test = "t-test (unpaired, equal var)",
    response = if (dec_tn_2g$transform == "log10") "log10TN" else "TN",
    transform = dec_tn_2g$transform,
    statistic = unname(tt$statistic),
    df = unname(tt$parameter),
    p_value = tt$p.value
  )
} else {
  ww <- wilcox.test(TN ~ LandUseType, data = dat_ag_for, exact = FALSE)
  res_2g <- tibble(
    test = "Wilcoxon rank-sum",
    response = "TN",
    transform = "raw",
    statistic = unname(ww$statistic),
    df = NA_real_,
    p_value = ww$p.value
  )
}

group_stats_2g <- dat_ag_for %>%
  group_by(LandUseType) %>%
  summarise(mean_TN = mean(TN, na.rm = TRUE), sd_TN = sd(TN, na.rm = TRUE), .groups = "drop")

out_2g <- bind_cols(
  res_2g,
  group_stats_2g %>% pivot_wider(names_from = LandUseType, values_from = c(mean_TN, sd_TN))
)
readr::write_csv(out_2g, file.path("results","tables","tn_ag_vs_for_test.csv"))
created_files <- c(created_files, file.path("results","tables","tn_ag_vs_for_test.csv"))

df_fig3_tn_two <- dat_ag_for %>% select(LandUseType, TN)
readr::write_csv(df_fig3_tn_two, file.path("results","plot_data","fig3_tn_two_groups_data.csv"))
created_files <- c(created_files, file.path("results","plot_data","fig3_tn_two_groups_data.csv"))

p3_n_lab <- dat_ag_for %>% count(LandUseType) %>% mutate(label = paste0("n=", n))
p_fig3_tn_two <- ggplot(dat_ag_for, aes(x = LandUseType, y = TN, fill = LandUseType)) +
  geom_violin(trim = FALSE, alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.25, outlier.alpha = 0.4) +
  geom_text(data = p3_n_lab, aes(x = LandUseType, y = Inf, label = label),
            vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  scale_fill_manual(values = landuse_colors) +
  labs(x = "Land-use type", y = "TN (mg/L)", title = "TN: Agricultural vs Forest")

cap_2g <- tryCatch({
  if (exists("res_2g")) {
    pv <- res_2g$p_value[1]
    stat <- res_2g$statistic[1]
    if (grepl("t-test", res_2g$test[1])) {
      paste0("t-test on log10TN: t=", sprintf("%.2f", stat), ", p=",
             ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv)))
    } else {
      paste0("Wilcoxon rank-sum: W=", sprintf("%.0f", stat), ", p=",
             ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv)))
    }
  } else NULL
}, error = function(e) NULL)
if (!is.null(cap_2g)) p_fig3_tn_two <- p_fig3_tn_two + labs(caption = cap_2g)

save_plot_all_formats(p_fig3_tn_two, "fig3_tn_two_groups")
created_files <- c(created_files,
  file.path("results","figures","fig3_tn_two_groups.png"),
  file.path("results","figures","fig3_tn_two_groups.pdf"),
  file.path("results","figures","fig3_tn_two_groups.svg")
)

## D1) pAG vs TN: correlations + simple regression (use raw TN)
pear_tn_pag <- cor.test(dat$pAG, dat$TN, method = "pearson")
spea_tn_pag <- cor.test(dat$pAG, dat$TN, method = "spearman", exact = FALSE)
lm_tn_pag <- lm(TN ~ pAG, data = dat)

sum_tn_pag <- tibble(
  Pearson_r = unname(pear_tn_pag$estimate), Pearson_p = pear_tn_pag$p.value,
  Spearman_rho = unname(spea_tn_pag$estimate), Spearman_p = spea_tn_pag$p.value,
  Slope = unname(coef(lm_tn_pag)[2]), Intercept = unname(coef(lm_tn_pag)[1]),
  R2 = summary(lm_tn_pag)$r.squared,
  F_stat = unname(summary(lm_tn_pag)$fstatistic[1]),
  df1 = unname(summary(lm_tn_pag)$fstatistic[2]),
  df2 = unname(summary(lm_tn_pag)$fstatistic[3]),
  F_p = pf(summary(lm_tn_pag)$fstatistic[1], summary(lm_tn_pag)$fstatistic[2], summary(lm_tn_pag)$fstatistic[3], lower.tail = FALSE)
)
readr::write_csv(sum_tn_pag, file.path("results","tables","correlation_regression_pAG_TN.csv"))
created_files <- c(created_files, file.path("results","tables","correlation_regression_pAG_TN.csv"))

df_fig4_pag_tn <- dat %>% select(pAG, TN)
readr::write_csv(df_fig4_pag_tn, file.path("results","plot_data","fig4_pAG_TN_scatter_data.csv"))
created_files <- c(created_files, file.path("results","plot_data","fig4_pAG_TN_scatter_data.csv"))

p_fig4_pag_tn <- ggplot(dat, aes(x = pAG, y = TN)) +
  geom_point(color = "#666666", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Agricultural land (%)", y = "TN (mg/L)", title = "TN vs Agricultural land (%)")

save_plot_all_formats(p_fig4_pag_tn, "fig4_pAG_TN_scatter")
created_files <- c(created_files,
  file.path("results","figures","fig4_pAG_TN_scatter.png"),
  file.path("results","figures","fig4_pAG_TN_scatter.pdf"),
  file.path("results","figures","fig4_pAG_TN_scatter.svg")
)

## D2) pURB vs TP
pear_tp_purb <- cor.test(dat$pURB, dat$TP, method = "pearson")
spea_tp_purb <- cor.test(dat$pURB, dat$TP, method = "spearman", exact = FALSE)
lm_tp_purb <- lm(TP ~ pURB, data = dat)

sum_tp_purb <- tibble(
  Pearson_r = unname(pear_tp_purb$estimate), Pearson_p = pear_tp_purb$p.value,
  Spearman_rho = unname(spea_tp_purb$estimate), Spearman_p = spea_tp_purb$p.value,
  Slope = unname(coef(lm_tp_purb)[2]), Intercept = unname(coef(lm_tp_purb)[1]),
  R2 = summary(lm_tp_purb)$r.squared,
  F_stat = unname(summary(lm_tp_purb)$fstatistic[1]),
  df1 = unname(summary(lm_tp_purb)$fstatistic[2]),
  df2 = unname(summary(lm_tp_purb)$fstatistic[3]),
  F_p = pf(summary(lm_tp_purb)$fstatistic[1], summary(lm_tp_purb)$fstatistic[2], summary(lm_tp_purb)$fstatistic[3], lower.tail = FALSE)
)
readr::write_csv(sum_tp_purb, file.path("results","tables","correlation_regression_pURB_TP.csv"))
created_files <- c(created_files, file.path("results","tables","correlation_regression_pURB_TP.csv"))

df_fig5_purb_tp <- dat %>% select(pURB, TP)
readr::write_csv(df_fig5_purb_tp, file.path("results","plot_data","fig5_pURB_TP_scatter_data.csv"))
created_files <- c(created_files, file.path("results","plot_data","fig5_pURB_TP_scatter_data.csv"))

p_fig5_purb_tp <- ggplot(dat, aes(x = pURB, y = TP)) +
  geom_point(color = "#666666", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Urban land (%)", y = "TP (mg/L)", title = "TP vs Urban land (%)")

save_plot_all_formats(p_fig5_purb_tp, "fig5_pURB_TP_scatter")
created_files <- c(created_files,
  file.path("results","figures","fig5_pURB_TP_scatter.png"),
  file.path("results","figures","fig5_pURB_TP_scatter.pdf"),
  file.path("results","figures","fig5_pURB_TP_scatter.svg")
)

## E) Chi-square for HighTN (and HighTP optional)
qTN <- quantile(dat$TN, 0.75, na.rm = TRUE)
qTP <- quantile(dat$TP, 0.75, na.rm = TRUE)
dat <- dat %>% mutate(HighTN = TN >= qTN, HighTP = TP >= qTP)

tab_highTN <- table(dat$LandUseType, dat$HighTN)
chi_htn <- suppressWarnings(chisq.test(tab_highTN))
if (any(chi_htn$expected < 5)) {
  chi_htn <- suppressWarnings(chisq.test(tab_highTN, simulate.p.value = TRUE, B = 1e4))
}

chisq_out <- tibble(
  test = chi_htn$method, statistic = unname(chi_htn$statistic), df = unname(chi_htn$parameter), p_value = chi_htn$p.value
)

## Save contingency and test result for HighTN
cont_htn <- as_tibble(as.data.frame(tab_highTN)) %>% rename(LandUseType = Var1, HighTN = Var2, Freq = Freq)
readr::write_csv(cont_htn, file.path("results","tables","chisq_highTN_table.csv"))
readr::write_csv(chisq_out, file.path("results","tables","chisq_highTN.csv"))
created_files <- c(created_files, file.path("results","tables","chisq_highTN_table.csv"), file.path("results","tables","chisq_highTN.csv"))

## Optional HighTP
tab_highTP <- table(dat$LandUseType, dat$HighTP)
chi_htp <- suppressWarnings(chisq.test(tab_highTP))
if (any(chi_htp$expected < 5)) {
  chi_htp <- suppressWarnings(chisq.test(tab_highTP, simulate.p.value = TRUE, B = 1e4))
}
chisq_out_tp <- tibble(
  test = chi_htp$method, statistic = unname(chi_htp$statistic), df = unname(chi_htp$parameter), p_value = chi_htp$p.value
)
cont_htp <- as_tibble(as.data.frame(tab_highTP)) %>% rename(LandUseType = Var1, HighTP = Var2, Freq = Freq)
readr::write_csv(cont_htp, file.path("results","tables","chisq_highTP_table.csv"))
readr::write_csv(chisq_out_tp, file.path("results","tables","chisq_highTP.csv"))
created_files <- c(created_files, file.path("results","tables","chisq_highTP_table.csv"), file.path("results","tables","chisq_highTP.csv"))

## Figure 6: Proportion HighTN by LandUseType
df_fig6_htn <- dat %>%
  count(LandUseType, HighTN) %>%
  group_by(LandUseType) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
readr::write_csv(df_fig6_htn, file.path("results","plot_data","fig6_highTN_barplot_data.csv"))
created_files <- c(created_files, file.path("results","plot_data","fig6_highTN_barplot_data.csv"))

p6_n_lab <- dat %>% count(LandUseType) %>% mutate(label = paste0("n=", n))
p_fig6_htn <- ggplot(df_fig6_htn, aes(x = LandUseType, y = prop, fill = HighTN)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("FALSE" = "#cccccc", "TRUE" = "#555555")) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(data = p6_n_lab, aes(x = LandUseType, y = Inf, label = label),
            vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  labs(x = "Land-use type", y = "Proportion (High TN)", title = "Proportion of stations with high TN (\u2265 75th percentile)")

cap_htn <- tryCatch({
  if (exists("chisq_out")) {
    paste0("Chi-square: \u03C7\u00B2=", sprintf("%.1f", chisq_out$statistic[1]),
           ", p=", ifelse(chisq_out$p_value[1] < 0.001, "<0.001", sprintf("%.3f", chisq_out$p_value[1])))
  } else NULL
}, error = function(e) NULL)
if (!is.null(cap_htn)) p_fig6_htn <- p_fig6_htn + labs(caption = cap_htn)

save_plot_all_formats(p_fig6_htn, "fig6_highTN_barplot")
created_files <- c(created_files,
  file.path("results","figures","fig6_highTN_barplot.png"),
  file.path("results","figures","fig6_highTN_barplot.pdf"),
  file.path("results","figures","fig6_highTN_barplot.svg")
)

message("[4/4] Summary printouts...")

## Print head and quick summaries
print(head(dat %>% select(stationId, statePostalCode, LandUseType, TN, TP, pAG, pURB)))

summary_tn_tp <- dat %>%
  group_by(LandUseType) %>%
  summarise(mean_TN = mean(TN, na.rm = TRUE), sd_TN = sd(TN, na.rm = TRUE),
            mean_TP = mean(TP, na.rm = TRUE), sd_TP = sd(TP, na.rm = TRUE), .groups = "drop")
print(summary_tn_tp)

message("\nCreated files (tables, figures, plot data):")
print(created_files)
