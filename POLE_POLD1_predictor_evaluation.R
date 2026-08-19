#!/usr/bin/env Rscript

# ================================================================================= #
# Evaluation of in-silico pathogenicity predictors in POLE and POLD1
#
# Analyses reported in the manuscript:
#   1. Violin plots of the predictor scores in the Benign/Pathogenic controls, in the
#      VUS and in the gnomAD variants, with the ClinGen PP3/BP4 Supporting cut-offs.
#      -> Figure 1A, Supplementary Figures 1-4
#   2. ROC curves per predictor, AUC with 95% CI (stratified BCa bootstrap) and
#      pairwise AUC comparison (paired DeLong test).
#      -> Figures 1B-C
#   3. Sensitivity analysis at the ClinGen thresholds for the AM, REVEL and Hybrid
#      (PP3 from AM, BP4 from REVEL) evidence rules.
#      -> Figure 1D
#   4. Violin plots of the AM and REVEL scores in the VUS and gnomAD variants against
#      the full published PP3/BP4 threshold ladder of each predictor.
#      -> Figures 1E-H
#   5. Distribution of the control, VUS and gnomAD variants across the AM and REVEL
#      evidence intervals defined by those thresholds.
#      -> Supplementary Tables 4-5
#
# ================================================================================= #

# General settings --------------------------------------------------------

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(pROC)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(boot)))
suppressWarnings(suppressPackageStartupMessages(library(Hmisc)))
suppressWarnings(suppressPackageStartupMessages(library(openxlsx)))

# Input / output. Paths are relative to the repository root.
input.excel <- "data/POLE_POLD1_predictors.xlsx"
results.dir <- "results"
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

# Minimum and maximum values defined for the predictors whose score is not in 0-1.
# Used ONLY to rescale the scores for the multi-predictor violin panels; every
# threshold-based analysis works on the raw score scale.
bayesdel_range <- c(-1.29333, 0.75731)
cadd_range     <- c(1, 48)

# Shared ggplot theme
gg_theme <- theme_bw() + theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10), axis.title = element_text(size = 12),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 12), legend.title = element_text(size = 12),
    strip.background = element_rect(color = "black", fill = "grey95"),
    strip.text = element_text(size = 12))

# Canonical predictor set, used for column selection, factor ordering and the
# threshold-free analyses (score distributions, ROC / AUC, DeLong).
# PrimateAI-3D (PAI3D) is included here but has NO ClinGen-calibrated PP3/BP4
# thresholds, so it carries no threshold guide-lines and takes no part in the
# threshold-based analyses.
predictors <- c("REVEL", "AM", "BD", "CADD", "PAI", "PAI3D")

# Reproducibility. `violin_seed` governs the point jitter of the violin plots;
# `boot_seed` is set ONCE before the per-predictor AUC loop and governs both the
# stratified bootstrap of the CIs and the label permutations, which share a single RNG
# stream (see the note above that loop).
violin_seed <- 1233
boot_seed   <- 20250707
n_boot      <- 5000L    # stratified bootstrap replicates for the AUC BCa CIs
n_perm      <- 20000L   # label permutations for the AUC > 0.5 test


# ClinGen thresholds -------------------------------------------------------

# Operating-point thresholds: ClinGen-calibrated PP3 / BP4 cut-offs per predictor, as a
# NAMED list of named numeric vectors holding EVERY evidence strength the predictor
# reached. Values are on the raw score scale (BD and CADD NOT rescaled).
#
# Each vector is an ordered LADDER, from the strongest benign to the strongest pathogenic
# cut-off. Every entry is the inclusive boundary of its tier, read in the direction of its
# criterion: `bp4_*` means "score <= value", `pp3_*` means "score >= value". A tier is
# reached only if the next-stronger one is not, e.g. for REVEL a score of 0.100 is
# <= bp4_moderate_b (0.183) but > bp4_moderate_a (0.052), so it carries BP4_Moderate-B.
#   - `bp4` / `pp3` (no suffix) are the SUPPORTING cut-offs.
#   - A strength that a predictor never reached ("-" in the source tables) is simply
#     ABSENT from its vector.
#   - Moderate-A is STRONGER than Moderate-B (Bergquist 2025 splits the single Pejaver
#     2022 Moderate tier in two). Predictors calibrated only by Pejaver keep a single
#     `*_moderate` tier.
#
# Sources:
#   - REVEL / AM / BD : Bergquist et al. 2025 (doi:10.1016/j.gim.2025.101402), which refines
#     Pejaver 2022 with the Moderate-A/B split and agrees with it on every shared cut-off.
#   - CADD / PAI      : Pejaver et al. 2022 (doi:10.1016/j.ajhg.2022.10.013); not
#     recalibrated by Bergquist 2025.
op_thresholds <- list(
    REVEL = c(
        bp4_very_strong =  0.003,   # Pejaver 2022 ONLY; the Bergquist 2025 ladder stops at Strong
        bp4_strong      =  0.016,
        bp4_moderate_a  =  0.052,
        bp4_moderate_b  =  0.183,
        bp4             =  0.290,   # BP4_Supporting
        pp3             =  0.644,   # PP3_Supporting
        pp3_moderate_b  =  0.773,
        pp3_moderate_a  =  0.879,
        pp3_strong      =  0.932),
    AM = c(
        bp4_moderate_a  =  0.070,
        bp4_moderate_b  =  0.099,
        bp4             =  0.169,
        pp3             =  0.792,
        pp3_moderate_b  =  0.906,
        pp3_moderate_a  =  0.972,
        pp3_strong      =  0.990),
    BD = c(                         # raw BayesDel scale
        bp4_moderate_a  = -0.520,
        bp4_moderate_b  = -0.360,
        bp4             = -0.180,
        pp3             =  0.130,
        pp3_moderate_b  =  0.270,
        pp3_moderate_a  =  0.410,
        pp3_strong      =  0.500),
    CADD = c(                       # raw CADD phred scale
        bp4_strong      =  0.150,
        bp4_moderate    = 17.300,
        bp4             = 22.700,
        pp3             = 25.300,
        pp3_moderate    = 28.100),
    PAI = c(
        bp4_moderate    =  0.362,
        bp4             =  0.483,
        pp3             =  0.790,
        pp3_moderate    =  0.867))


# Helper functions --------------------------------------------------------

#' Rescale BD / CADD scores to the 0-1 range used in the multi-predictor violin panels
#' (REVEL, AM, PAI and PAI3D are already 0-1). Applied to both the plotted scores and the
#' threshold guide-lines so the two share one definition.
rescale_score <- function(pred, x) {
    ifelse(pred == "BD",   (x - bayesdel_range[1]) / (bayesdel_range[2] - bayesdel_range[1]),
    ifelse(pred == "CADD", (x - cadd_range[1])     / (cadd_range[2] - cadd_range[1]),
           x))
}

#' Collapse the raw Classification terms into the plotting subgroups used in
#' Supplementary Figure 1. VUS variants map to NA and are excluded from that figure.
#'   Benign     : B, LB                     -> "B/LB"
#'   Pathogenic : LP, P                     -> "P/LP"
#'                hotVUS                    -> "hotVUS"
#'                LP - somatic, P - somatic -> "Somatic"
classification_group <- function(x) {
    x <- trimws(x)
    dplyr::case_when(
        x %in% c("B", "LB")                     ~ "B/LB",
        x %in% c("LP", "P")                     ~ "P/LP",
        x == "hotVUS"                           ~ "hotVUS",
        x %in% c("LP - somatic", "P - somatic") ~ "Somatic",
        TRUE                                    ~ NA_character_)
}

# Classification subgroups (Supplementary Figure 1).
classification_levels <- c("B/LB", "P/LP", "hotVUS", "Somatic")
classification_pal <- c(
    "B/LB"    = "#2166AC",   # blue   - benign
    "P/LP"    = "#B2182B",   # red    - pathogenic (germline)
    "hotVUS"  = "#9970AB",   # purple - hotspot VUS
    "Somatic" = "#E08214")   # orange - somatic-designated pathogenic

# Gene subgroups (Supplementary Figure 2).
gene_levels <- c("POLE", "POLD1")
gene_pal    <- c("POLE" = "#1B7837", "POLD1" = "#762A83")   # green / purple

# Long table of the PP3/BP4 Supporting thresholds for the multi-predictor violin panels
# (one row per predictor x criterion), rescaled to match the plotted BD / CADD scores.
# geom_hline reads this per facet to draw a threshold line for each predictor.
op_thresholds_df <- lapply(names(op_thresholds), function(p) {
    data.frame(
        Predictor = p,
        criterion = c("PP3_Supporting", "BP4_Supporting"),
        value     = c(op_thresholds[[p]][["pp3"]], op_thresholds[[p]][["bp4"]]),
        stringsAsFactors = FALSE)
}) %>%
    bind_rows() %>%
    mutate(
        value     = rescale_score(Predictor, value),
        Predictor = factor(Predictor, levels = predictors),
        criterion = factor(criterion, levels = c("PP3_Supporting", "BP4_Supporting")))

# Reusable ggplot layers: PP3/BP4 Supporting guide-lines, one per predictor facet.
threshold_layers <- list(
    geom_hline(data = op_thresholds_df,
               aes(yintercept = value, color = criterion),
               linetype = "dashed", linewidth = 0.3),
    scale_color_manual(
        name = NULL,
        values = c(PP3_Supporting = "red", BP4_Supporting = "#2780E3"),
        labels = c(PP3_Supporting = "PP3 Supporting", BP4_Supporting = "BP4 Supporting")))

# The same guide-lines with FIXED colours and no colour scale, for the plots that already
# map `color` to something else (the classification / gene subgroup figures). Keeping the
# colour aesthetic free avoids two competing colour scales in one plot.
threshold_layers_fixed <- list(
    geom_hline(data = op_thresholds_df %>% filter(criterion == "PP3_Supporting"),
               aes(yintercept = value), color = "red", linetype = "dashed",
               linewidth = 0.3),
    geom_hline(data = op_thresholds_df %>% filter(criterion == "BP4_Supporting"),
               aes(yintercept = value), color = "#2780E3", linetype = "dashed",
               linewidth = 0.3))

# Label, ACMG points and drawing style of each tier of a threshold ladder. Single source
# of truth for both the ladder figures (ladder_layers(), Figures 1E-1H) and the tier
# tables (tier_interval_table(), Supplementary Tables 4-5).
# Rows are in drawing order (strongest pathogenic -> strongest benign); the Bergquist
# Moderate-A / Moderate-B tiers are labelled "Moderate-1" / "Moderate-2" throughout.
# `points` are the ACMG/Tavtigian 2018 points claimed by the tier, as an unsigned
# magnitude (the tables sign them: pathogenic positive, benign negative).
# Only the AM and REVEL ladders are drawn / tabulated, so the single unsplit `*_moderate`
# tier of the Pejaver-only predictors (CADD, PAI) is not listed here.
ladder_style <- data.frame(
    tier = c("pp3_strong", "pp3_moderate_a", "pp3_moderate_b", "pp3",
             "bp4", "bp4_moderate_b", "bp4_moderate_a", "bp4_strong", "bp4_very_strong"),
    label = c("PP3 Strong", "PP3 Moderate-1", "PP3 Moderate-2", "PP3 Supporting",
              "BP4 Supporting", "BP4 Moderate-2", "BP4 Moderate-1", "BP4 Strong",
              "BP4 Very Strong"),
    points = c(4, 3, 2, 1,
               1, 2, 3, 4, 8),
    linetype = c("twodash", "dotdash", "dashed", "dotted",
                 "dotted", "dashed", "dotdash", "twodash", "longdash"),
    stringsAsFactors = FALSE)
ladder_style$color <- ifelse(startsWith(ladder_style$tier, "pp3"), "red", "#2780E3")

# Named lookups over the ladder, for code that indexes by tier name.
tier_label  <- setNames(ladder_style$label,  ladder_style$tier)
tier_points <- setNames(ladder_style$points, ladder_style$tier)

#' ggplot layers drawing the FULL published ladder of one predictor (Figures 1E-1H):
#' one line per tier, the tier name to the right of the violin and its cut-off value to
#' the left. The layers inherit the plot data, so a label is drawn once per plotted
#' variant, exactly as in the published figures.
#' @param pred predictor name; must have a ladder in `op_thresholds`.
#' @param label_y named vector overriding the y position of a tier's labels. AM's
#'   PP3_Strong cut-off (0.990) sits too close to the panel top, so its labels are
#'   pushed to y = 1.
ladder_layers <- function(pred, label_y = NULL) {
    thr <- op_thresholds[[pred]]
    sty <- ladder_style[ladder_style$tier %in% names(thr), ]
    unlist(lapply(seq_len(nrow(sty)), function(i) {
        tier  <- sty$tier[i]
        value <- unname(thr[[tier]])
        y_lab <- if (!is.null(label_y) && tier %in% names(label_y)) label_y[[tier]] else value - 0.01
        list(
            geom_hline(yintercept = value, color = sty$color[i],
                       linetype = sty$linetype[i], linewidth = 0.2),
            geom_text(label = sty$label[i], nudge_x =  0.50, y = y_lab, size = 1.5),
            geom_text(label = value,        nudge_x = -0.55, y = y_lab, size = 1.5))
    }), recursive = FALSE)
}


# Load data ---------------------------------------------------------------

# Pathogenic / Benign controls. `Categorization` is the classification made WITHOUT the
# PP3/BP4 and PM1_Supporting (exonuclease / DNA-binding domain) evidence codes, so the
# labels are independent of the predictor scores evaluated here: the control set carries
# no circularity by construction.
# `Classification` carries the specific term behind each Categorization (Benign: B / LB;
# Pathogenic: LP / P / hotVUS / LP - somatic / P - somatic; VUS is identical in both
# columns). It is used only to colour Supplementary Figure 1.
controls_vars <- read_xlsx(input.excel, na = c("NA", "", "n.a."), sheet = "controls") %>%
    select(all_of(predictors), Categorization, Classification, GENE)

# VUS: same columns.
vus_vars <- read_xlsx(input.excel, na = c("NA", "", "n.a."), sheet = "VUS") %>%
    select(all_of(predictors), Categorization, Classification, GENE)

# gnomAD missense variants of the exonuclease domains (variants also seen in TCGA are
# already excluded from this sheet).
gnomad_vars <- read_xlsx(input.excel, na = c("NA", "", "n.a."), sheet = "gnomAD") %>%
    select(all_of(predictors), GENE)

# Controls and VUS together, for the figures that show the three groups side by side.
input_vars <- bind_rows(controls_vars, vus_vars)

# Long format, with the factor ordering used throughout
df_long <- input_vars %>%
    pivot_longer(all_of(predictors), names_to = "Predictor", values_to = "Score") %>%
    mutate(
        Predictor      = factor(Predictor, levels = predictors),
        Categorization = factor(Categorization, levels = c("Benign", "VUS", "Pathogenic")),
        Classification = factor(classification_group(Classification),
                                levels = classification_levels),
        GENE           = factor(GENE, levels = gene_levels))

# Control table with a binary label y (Pathogenic = 1, Benign = 0), used by the ROC / AUC
# and the evidence-rule analyses. NA scores are dropped per predictor inside each helper,
# so n_P / n_B are reported per predictor.
dat_controls <- controls_vars %>%
    mutate(y = as.integer(Categorization == "Pathogenic"))

dat_vus <- vus_vars



# ========================================================================= #
# 1. PREDICTOR SCORE DISTRIBUTIONS
# ------------------------------------------------------------------------- #
# Violin + boxplot + jittered points of every predictor score, faceted by predictor,
# with the ClinGen PP3/BP4 Supporting cut-offs as dashed guide-lines. BD and CADD are
# rescaled to 0-1 (scores and guide-lines alike) so all predictors share one y axis.
# ========================================================================= #

benign_N     <- sum(controls_vars$Categorization == "Benign")
pathogenic_N <- sum(controls_vars$Categorization == "Pathogenic")
vus_N        <- nrow(vus_vars)

plot_vars <- df_long %>%
    mutate(Score = rescale_score(Predictor, Score),
           group = case_when(
               Categorization == "Benign"     ~ paste0("Benign (n=", benign_N, ")"),
               Categorization == "Pathogenic" ~ paste0("Pathogenic (n=", pathogenic_N, ")"),
               TRUE                           ~ paste0("VUS (n=", vus_N, ")"))) %>%
    arrange(Categorization)
plot_vars$group <- factor(plot_vars$group, levels = unique(plot_vars$group))


# --- Figure 1A: Benign / Pathogenic controls -----------------------------
set.seed(violin_seed)
pdf(file.path(results.dir, "Figure1A_predictor_scores_B-P.pdf"), width = 7, height = 3.5)
ggplot(plot_vars %>% filter(Categorization != "VUS"),
       aes(x = Categorization, y = Score, fill = group)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.2, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
    threshold_layers +
    gg_theme +
    theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Supplementary Figure 1: controls, points coloured by classification subgroup ----
# The violins / boxplots still summarise the Benign and Pathogenic groups; the jittered
# points are coloured by `Classification` to show how the subgroups distribute within
# each group. The guide-lines use the fixed-colour variant so the colour scale stays
# free for the classification terms.
set.seed(violin_seed)
pdf(file.path(results.dir, "SupplementaryFigure1_predictor_scores_B-P_classification.pdf"),
    width = 7.5, height = 4.2)
ggplot(plot_vars %>% filter(Categorization != "VUS"),
       aes(x = Categorization, y = Score)) +
    geom_violin(aes(fill = group), linewidth = 0.3, alpha = 0.55) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(aes(color = Classification), size = 0.6, width = 0.2, alpha = 0.9) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores", fill = NULL, color = NULL,
         caption = "Dashed guide-lines: red = PP3_Supporting, blue = BP4_Supporting") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
    scale_color_manual(values = classification_pal, drop = TRUE) +
    threshold_layers_fixed +
    gg_theme +
    guides(fill  = guide_legend(order = 1, nrow = 1),
           color = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.2))) +
    theme(
        legend.margin = margin(t = 0, b = 0), legend.spacing.y = unit(1, "pt"),
        plot.caption = element_text(size = 7, color = "grey30"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Supplementary Figure 2: controls, points coloured by gene -----------
# Same construction as Supplementary Figure 1, but the colour aesthetic carries GENE.
# Unlike Classification, GENE crosses BOTH categories, so each violin mixes POLE and
# POLD1 points.
set.seed(violin_seed)
pdf(file.path(results.dir, "SupplementaryFigure2_predictor_scores_B-P_gene.pdf"),
    width = 7.5, height = 4.2)
ggplot(plot_vars %>% filter(Categorization != "VUS"),
       aes(x = Categorization, y = Score)) +
    geom_violin(aes(fill = group), linewidth = 0.3, alpha = 0.55) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(aes(color = GENE), size = 0.6, width = 0.2, alpha = 0.9) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores", fill = NULL, color = NULL,
         caption = "Dashed guide-lines: red = PP3_Supporting, blue = BP4_Supporting") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
    scale_color_manual(values = gene_pal, drop = TRUE) +
    threshold_layers_fixed +
    gg_theme +
    guides(fill  = guide_legend(order = 1, nrow = 1),
           color = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.2))) +
    theme(
        legend.margin = margin(t = 0, b = 0), legend.spacing.y = unit(1, "pt"),
        plot.caption = element_text(size = 7, color = "grey30"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Supplementary Figure 3: VUS -----------------------------------------
set.seed(violin_seed)
pdf(file.path(results.dir, "SupplementaryFigure3_predictor_scores_VUS.pdf"),
    width = 6, height = 3)
ggplot(plot_vars %>% filter(Categorization == "VUS"),
       aes(x = Categorization, y = Score, fill = group)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.4, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = c("#FFDB6D")) +
    threshold_layers +
    gg_theme +
    theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Supplementary Figure 4: gnomAD variants -----------------------------
plot_gnomad <- gnomad_vars %>%
    pivot_longer(all_of(predictors), names_to = "Predictor", values_to = "Score") %>%
    mutate(Predictor = factor(Predictor, levels = predictors),
           X = "",
           Score = rescale_score(Predictor, Score))

set.seed(violin_seed)
pdf(file.path(results.dir, "SupplementaryFigure4_predictor_scores_gnomAD.pdf"),
    width = 7, height = 3.5)
ggplot(plot_gnomad, aes(x = X, y = Score)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.4, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    threshold_layers +
    gg_theme +
    theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()



# ========================================================================= #
# 2. AM AND REVEL SCORES AGAINST THEIR FULL THRESHOLD LADDER
# ------------------------------------------------------------------------- #
# The same violins restricted to AM and to REVEL, the two predictors carrying a complete
# ClinGen ladder, with EVERY published PP3/BP4 tier drawn and labelled. Scores are shown
# on their native 0-1 scale (no rescaling is needed for AM or REVEL).
# ========================================================================= #

# --- Figure 1E: VUS, AM ladder -------------------------------------------
set.seed(violin_seed)
pdf(file.path(results.dir, "Figure1E_VUS_AM_thresholds.pdf"), width = 4, height = 3)
ggplot(plot_vars %>% filter(Categorization == "VUS", Predictor == "AM") %>%
           mutate(Predictor = paste0("AM (", n(), " VUS variants)")),
       aes(x = Categorization, y = Score, fill = group)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.4, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ladder_layers("AM", label_y = c(pp3_strong = 1)) +
    scale_fill_manual(values = c("#FFDB6D")) +
    gg_theme +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Figure 1F: VUS, REVEL ladder ----------------------------------------
set.seed(violin_seed)
pdf(file.path(results.dir, "Figure1F_VUS_REVEL_thresholds.pdf"), width = 4, height = 3)
ggplot(plot_vars %>% filter(Categorization == "VUS", Predictor == "REVEL") %>%
           mutate(Predictor = paste0("REVEL (", n(), " VUS variants)")),
       aes(x = Categorization, y = Score, fill = group)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.4, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ladder_layers("REVEL") +
    scale_fill_manual(values = c("#FFDB6D")) +
    gg_theme +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Figure 1G: gnomAD, AM ladder ----------------------------------------
set.seed(violin_seed)
pdf(file.path(results.dir, "Figure1G_gnomAD_AM_thresholds.pdf"), width = 4, height = 3)
ggplot(plot_gnomad %>% filter(Predictor == "AM") %>%
           mutate(Predictor = paste0("AM (", n(), " gnomAD variants)")),
       aes(x = X, y = Score)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.4, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ladder_layers("AM", label_y = c(pp3_strong = 1)) +
    gg_theme +
    theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


# --- Figure 1H: gnomAD, REVEL ladder -------------------------------------
set.seed(violin_seed)
pdf(file.path(results.dir, "Figure1H_gnomAD_REVEL_thresholds.pdf"), width = 4, height = 3)
ggplot(plot_gnomad %>% filter(Predictor == "REVEL") %>%
           mutate(Predictor = paste0("REVEL (", n(), " gnomAD variants)")),
       aes(x = X, y = Score)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.1, outliers = FALSE, lwd = 0.2, fill = "grey") +
    geom_jitter(color = "black", size = 0.6, pch = 21, width = 0.4, stroke = 0.3) +
    facet_wrap(. ~ Predictor, nrow = 1) +
    labs(y = "Prediction scores") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ladder_layers("REVEL") +
    gg_theme +
    theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()



# ========================================================================= #
# 3. ROC CURVES, AUC WITH 95% CI AND PAIRWISE AUC COMPARISON
# ------------------------------------------------------------------------- #
# Discrimination between the Pathogenic and Benign controls, per predictor:
#   - ROC curve of every predictor on one panel (Figure 1B),
#   - AUC with a 95% stratified BCa bootstrap CI, and a label-permutation test of
#     AUC > 0.5 (Figure 1B),
#   - pairwise AUC comparison with the paired DeLong test (Figure 1C).
#
# All predictor scores are oriented "higher = more pathogenic"; the ROC direction of the
# CI / permutation / DeLong analyses is pinned to "<" (controls = Benign < cases =
# Pathogenic) so that the resamples cannot be flipped by pROC's `direction = "auto"`,
# which would force every replicate to AUC >= 0.5 and bias the null distributions and
# the intervals upward.
# ========================================================================= #

input_roc <- df_long %>%
    filter(Categorization != "VUS") %>%
    mutate(Categorization = factor(Categorization, levels = c("Benign", "Pathogenic")))

roc_metrics <- lapply(predictors, function(p) {
    aux_df    <- input_roc %>% filter(Predictor == p)
    aux_rocobj <- roc(aux_df, "Categorization", "Score", auc = TRUE, ci = TRUE)
    data.frame(
        PREDICTOR     = p,
        ROC_THRESHOLD = aux_rocobj$thresholds,
        TPR           = aux_rocobj$sensitivities,
        FPR           = 1 - aux_rocobj$specificities,
        AUC           = round(aux_rocobj$auc, 4),
        CI_LOW        = aux_rocobj$ci[1],
        CI_HIGH       = aux_rocobj$ci[3])
}) %>% bind_rows()


# --- Figure 1B: ROC curves of all predictors on one panel ----------------
plot_roc <- roc_metrics %>%
    mutate(PREDICTOR = factor(PREDICTOR, levels = predictors)) %>%
    arrange(PREDICTOR, FPR, TPR) %>%
    mutate(LABEL = paste0(PREDICTOR, " (AUC=", AUC, ")"))
plot_roc$LABEL <- factor(plot_roc$LABEL, levels = unique(plot_roc$LABEL))

pdf(file.path(results.dir, "Figure1B_ROC_curves.pdf"), width = 7, height = 4)
ggplot(plot_roc, aes(x = FPR, y = TPR, group = LABEL, color = LABEL)) +
    geom_line(linewidth = 0.4) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    labs(
        title = "ROC curves by predictor",
        x = "FPR (1 - Specificity)",
        y = "TPR (Sensitivity)") +
    gg_theme +
    guides(color = guide_legend(nrow = 2)) +
    theme(legend.title = element_blank())
dev.off()


# --- AUC helpers ---------------------------------------------------------

#' Point AUC with fixed orientation.
#' @param y integer 0/1 labels (1 = Pathogenic).
#' @param s numeric predictor scores (higher = more pathogenic).
#' @return AUC as a plain numeric, or NA if a class is empty / all scores are NA.
auc_point <- function(y, s) {
    ok <- !is.na(s) & !is.na(y)
    y <- y[ok]; s <- s[ok]
    if (length(unique(y)) < 2L) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(response = y, predictor = s,
                                   levels = c(0, 1), direction = "<",
                                   quiet = TRUE)))
}

#' Stratified BCa bootstrap CI for the AUC.
#' Resamples Pathogenic and Benign SEPARATELY (boot strata = y) so n_P and n_B are
#' preserved in every replicate. The tryCatch is a safety net that returns NA bounds
#' instead of aborting the run if boot.ci fails.
#' @return list(auc, ci_low, ci_high).
boot_auc_bca <- function(y, s, R = n_boot) {
    ok <- !is.na(s) & !is.na(y)
    y <- y[ok]; s <- s[ok]
    obs <- auc_point(y, s)
    boot_stat <- function(data, idx) auc_point(data$y[idx], data$s[idx])
    bt <- boot::boot(data = data.frame(y = y, s = s), statistic = boot_stat,
                     R = R, strata = factor(y))
    ci <- tryCatch(boot::boot.ci(bt, type = "bca", conf = 0.95),
                   error = function(e) NULL)
    if (is.null(ci) || is.null(ci$bca)) {
        warning("BCa CI not estimable; reporting NA bounds", call. = FALSE)
        return(list(auc = obs, ci_low = NA_real_, ci_high = NA_real_))
    }
    list(auc = obs, ci_low = ci$bca[4], ci_high = ci$bca[5])
}

#' Permutation test for AUC > 0.5.
#' Permutes the labels `n_perm` times, recomputes the AUC, and returns the empirical
#' one-sided p = (1 + #{perm AUC >= observed}) / (n_perm + 1). At the floor the p is
#' reported as "< 1/(N+1)".
#' @return list(obs_auc, p_value, p_label, n_perm).
perm_auc_p <- function(y, s, reps = n_perm) {
    ok <- !is.na(s) & !is.na(y)
    y <- y[ok]; s <- s[ok]
    obs <- auc_point(y, s)
    if (is.na(obs)) return(list(obs_auc = NA_real_, p_value = NA_real_,
                                p_label = NA_character_, n_perm = reps))
    perm <- vapply(seq_len(reps),
                   function(i) auc_point(sample(y), s),
                   numeric(1))
    p <- (1 + sum(perm >= obs, na.rm = TRUE)) / (reps + 1)
    floor_p <- 1 / (reps + 1)
    p_label <- if (p <= floor_p) paste0("< ", signif(floor_p, 2)) else signif(p, 3)
    list(obs_auc = obs, p_value = p, p_label = as.character(p_label), n_perm = reps)
}

#' Pairwise paired DeLong test across predictors (shared control set).
#' For each predictor pair the test runs on the rows where BOTH scores are non-NA
#' (the paired = TRUE requirement). Returns an upper-triangular matrix of p-values.
pairwise_delong <- function(wide, preds = predictors) {
    m <- matrix(NA_real_, nrow = length(preds), ncol = length(preds),
                dimnames = list(preds, preds))
    for (i in seq_along(preds)) {
        for (j in seq_along(preds)) {
            if (j <= i) next
            pi <- preds[i]; pj <- preds[j]
            ok <- !is.na(wide[[pi]]) & !is.na(wide[[pj]])
            y  <- wide$y[ok]
            if (length(unique(y)) < 2L) next
            r1 <- pROC::roc(response = y, predictor = wide[[pi]][ok],
                            levels = c(0, 1), direction = "<", quiet = TRUE)
            r2 <- pROC::roc(response = y, predictor = wide[[pj]][ok],
                            levels = c(0, 1), direction = "<", quiet = TRUE)
            m[i, j] <- tryCatch(
                pROC::roc.test(r1, r2, paired = TRUE, method = "delong")$p.value,
                error = function(e) NA_real_)
        }
    }
    m
}


# --- Figure 1B: AUC with 95% BCa CI and permutation p, per predictor -----
# The bootstrap and the permutation test draw from ONE RNG stream seeded once, before
# the loop. Their order within an iteration (bootstrap first, then permutation) is
# therefore part of the reproducibility contract: reordering them, or adding/removing a
# random draw anywhere in the loop, shifts every predictor after the first one.
set.seed(boot_seed)
auc_results <- lapply(predictors, function(p) {
    s  <- dat_controls[[p]]; y <- dat_controls$y
    ok <- !is.na(s) & !is.na(y)
    bca  <- boot_auc_bca(y, s)
    perm <- perm_auc_p(y, s)
    data.frame(
        predictor = p,
        n_P       = sum(y[ok] == 1),
        n_B       = sum(y[ok] == 0),
        AUC       = round(bca$auc, 4),
        CI_low    = bca$ci_low,
        CI_high   = bca$ci_high,
        CI_95     = sprintf("%.3f-%.3f", bca$ci_low, bca$ci_high),
        perm_p    = perm$p_label,
        stringsAsFactors = FALSE)
}) %>% bind_rows() %>%
    mutate(predictor = factor(predictor, levels = predictors)) %>%
    arrange(predictor)

# --- Figure 1C: pairwise paired DeLong p-value matrix --------------------
delong_mat <- pairwise_delong(dat_controls)
delong_df  <- cbind(predictor = rownames(delong_mat),
                    as.data.frame(delong_mat), row.names = NULL)

write.xlsx(
    list(AUC_results      = auc_results,   # one row per predictor
         pairwiseDeLong_p = delong_df),    # paired DeLong p-value matrix
    file = file.path(results.dir, "POLE_POLD1_AUC_results.xlsx"),
    overwrite = TRUE)

# --- Console report ------------------------------------------------------
cat("\n", strrep("=", 73), "\n", sep = "")
cat("AUC ON THE P/B CONTROLS\n")
cat("  seed =", boot_seed, "| bootstrap reps =", n_boot, "| permutations =", n_perm, "\n")
cat(strrep("=", 73), "\n", sep = "")
cat("\nAUC (95% BCa CI) and permutation p:\n")
print(auc_results[, c("predictor", "n_P", "n_B", "AUC", "CI_95", "perm_p")],
      row.names = FALSE)

cat("\n-- Is discrimination better than chance (CI lower bound > 0.5)? --\n")
for (p in predictors) {
    i <- auc_results[auc_results$predictor == p, ]
    verdict <- if (!is.na(i$CI_low) && i$CI_low > 0.5) {
        "YES - CI lower bound > 0.5"
    } else {
        "NOT confirmed - CI includes 0.5 (or not estimable)"
    }
    cat(sprintf("  %-6s AUC=%.3f (%s)  => %s\n", p, i$AUC, i$CI_95, verdict))
}

cat("\n-- Pairwise paired DeLong test (p-values) --\n")
print(delong_df, row.names = FALSE)
cat(strrep("=", 73), "\n", sep = "")



# ========================================================================= #
# 4. SENSITIVITY ANALYSIS AT THE ClinGen THRESHOLDS: AM, REVEL AND HYBRID
# ------------------------------------------------------------------------- #
# Operating characteristics of three evidence rules, all of them using the PUBLISHED
# ClinGen PP3/BP4 Supporting cut-offs (nothing is re-fitted):
#   - AM     : both codes from AlphaMissense,
#   - REVEL  : both codes from REVEL,
#   - Hybrid : PP3 from AM, BP4 from REVEL, conflicts left uncoded.
#
# The hybrid is motivated by the two OPPOSITE mis-calibrations of the ClinVar-derived
# thresholds in these genes: REVEL's pathogenic scores run low (PP3 rarely fires) but it
# never mis-fires BP4 on a pathogenic control, whereas AM's benign scores run high (BP4
# rarely fires) but PP3 fires on every pathogenic control. The pathogenic evidence is
# therefore taken from AM and the benign evidence from REVEL. When the two conflict
# (AM says PP3 and REVEL says BP4) the variant gets NO in-silico code, the conservative
# choice.
#
# Metrics on the labelled controls (exact Clopper-Pearson 95% CIs):
#   sensitivity = TP / n_P  (PP3 yield among true Pathogenic)
#   specificity = TN / n_B  (BP4 yield among true Benign)
#   PPV = TP / (TP + FP)    (among PP3-called variants)
#   NPV = TN / (TN + FN)    (among BP4-called variants)
# Indeterminate and Conflict variants are a genuine "no call": they are excluded from the
# PPV / NPV denominators and reported as separate counts, NOT folded into the benign class.
#
# CAVEAT for the manuscript: the hybrid rule is SELECTED and EVALUATED on the same
# controls, so its advantage here is mildly optimistic and needs independent validation.
# ========================================================================= #

call_levels <- c("BP4", "Indeterminate", "Conflict", "PP3")

#' PP3 / BP4 / Indeterminate call of one predictor at its Supporting cut-offs.
#' Reaching ANY tier of the pathogenic ladder is equivalent to clearing its weakest
#' (Supporting) cut-off, and likewise on the benign side; since pp3 > bp4 in every
#' ladder the two rules are mutually exclusive.
#' @return character vector, NA where the score is missing.
evidence_call <- function(s, thr) {
    ifelse(is.na(s), NA_character_,
    ifelse(s >= thr[["pp3"]], "PP3",
    ifelse(s <= thr[["bp4"]], "BP4", "Indeterminate")))
}

#' Evidence call of each variant under one of the three candidate rules.
#' A missing score simply cannot produce its side's code, so a variant without an AM
#' score can still receive BP4 from REVEL, and vice versa.
#' @return factor over `call_levels`, one element per row of `d`.
rule_calls <- function(d, rule = c("AM", "REVEL", "Hybrid")) {
    rule  <- match.arg(rule)
    am    <- evidence_call(d$AM,    op_thresholds$AM)
    revel <- evidence_call(d$REVEL, op_thresholds$REVEL)
    out <- switch(
        rule,
        AM     = ifelse(is.na(am),    "Indeterminate", am),
        REVEL  = ifelse(is.na(revel), "Indeterminate", revel),
        Hybrid = {
            is_pp3 <- !is.na(am)    & am    == "PP3"      # pathogenic evidence: AM
            is_bp4 <- !is.na(revel) & revel == "BP4"      # benign evidence: REVEL
            ifelse(is_pp3 & is_bp4, "Conflict",
            ifelse(is_pp3, "PP3",
            ifelse(is_bp4, "BP4", "Indeterminate")))
        })
    factor(out, levels = call_levels)
}

#' Operating characteristics of a rule on the labelled controls, with exact 95% CIs.
rule_metrics <- function(calls, truth) {
    fmt <- function(x, n) {
        if (n == 0) return(NA_character_)
        cc <- Hmisc::binconf(x, n, method = "exact")
        sprintf("%.3f [%.3f-%.3f]", cc[1], cc[2], cc[3])
    }
    n_P <- sum(truth == "Pathogenic"); n_B <- sum(truth == "Benign")
    TP <- sum(calls == "PP3" & truth == "Pathogenic")
    FP <- sum(calls == "PP3" & truth == "Benign")
    TN <- sum(calls == "BP4" & truth == "Benign")
    FN <- sum(calls == "BP4" & truth == "Pathogenic")
    data.frame(
        n_PP3 = TP + FP, n_BP4 = TN + FN,
        n_indeterminate = sum(calls == "Indeterminate"),
        n_conflict      = sum(calls == "Conflict"),
        called_rate     = mean(calls %in% c("PP3", "BP4")),
        sensitivity = fmt(TP, n_P),        # PP3 yield among true Pathogenic
        specificity = fmt(TN, n_B),        # BP4 yield among true Benign
        PPV = fmt(TP, TP + FP), NPV = fmt(TN, TN + FN),
        stringsAsFactors = FALSE)
}

rules <- c("AM", "REVEL", "Hybrid")

# Controls: accuracy + yield, one row per rule
rule_controls <- lapply(rules, function(rl) {
    cbind(rule = rl,
          rule_metrics(rule_calls(dat_controls, rl), dat_controls$Categorization))
}) %>% bind_rows()

write.xlsx(
    list(controls = rule_controls),   # accuracy and yield on the P/B controls
    file = file.path(results.dir, "POLE_POLD1_evidence_rules.xlsx"),
    overwrite = TRUE)

# --- Console report (Figure 1D) ------------------------------------------
cat("\n", strrep("=", 73), "\n", sep = "")
cat("EVIDENCE RULES AT THE ClinGen SUPPORTING THRESHOLDS\n")
cat("  AM | REVEL | Hybrid (PP3 from AM, BP4 from REVEL, conflicts uncoded)\n")
cat(strrep("=", 73), "\n", sep = "")
cat("\n-- Controls: accuracy and yield by rule --\n")
print(rule_controls, row.names = FALSE)
cat(strrep("=", 73), "\n", sep = "")



# ========================================================================= #
# 5. TIER-INTERVAL DISTRIBUTION OF THE CONTROL, VUS AND gnomAD VARIANTS
#    Supplementary Table 4 (AM) and Supplementary Table 5 (REVEL)
# ------------------------------------------------------------------------- #
# For AM and REVEL, the WHOLE score range is partitioned into the intervals cut by the
# published ClinGen ladder (strongest PP3 tier down to the strongest BP4 tier, with the
# Indeterminate gap in the middle), and each variant set is counted in each interval:
# the P/B controls, the VUS and the gnomAD variants. The CONTROLS column is where the
# tiers can be read against a known answer, and is the reference for the two unlabelled
# sets beside it.
#
# The intervals are a strict partition - every scored variant falls in exactly one, so
# the counts of a column sum to the number of scored variants in that set. Each
# pathogenic tier is bounded BELOW by its own cut-off (inclusive) and ABOVE by the
# next-stronger cut-off (exclusive), mirrored on the benign side; the outermost tiers are
# closed at the limits of the score scale. THRESHOLDS is printed LOW -> HIGH in standard
# interval notation, the bracket marking each bound's inclusivity, e.g. AM PP3 Strong =
# "[0.990, 1.000]" and PP3 Moderate-1 = "[0.972, 0.990)". Scores are compared RAW, as
# everywhere else in the script, so these counts agree with the evidence rules above.
#
# POINTS are SIGNED ACMG/Tavtigian 2018 points: pathogenic evidence positive, benign
# evidence negative, Indeterminate 0.
# ========================================================================= #

#' Native score range of a predictor, used to close the outermost tiers.
predictor_range <- function(p) {
    switch(p,
           BD   = bayesdel_range,
           CADD = cadd_range,
           c(0, 1))          # REVEL / AM / PAI / PAI3D are 0-1 scores
}

#' Full tier partition of a predictor's score range + per-set counts.
#' @param pred predictor name (must have a ladder in `op_thresholds`).
#' @param sets named list of data.frames, each holding a `pred` column.
#' @param drop_tiers tier names to leave out of the ladder for this table. The partition
#'   stays complete: the neighbouring tier simply absorbs the dropped one's interval.
#' @return data.frame(EVIDENCE, POINTS, THRESHOLDS, <one column per set>).
tier_interval_table <- function(pred, sets, drop_tiers = character(0)) {
    thr <- op_thresholds[[pred]]
    thr <- thr[setdiff(names(thr), drop_tiers)]
    rng <- predictor_range(pred)
    pp3_t <- names(thr)[startsWith(names(thr), "pp3")]
    bp4_t <- names(thr)[startsWith(names(thr), "bp4")]
    pp3_t <- pp3_t[order(thr[pp3_t], decreasing = TRUE)]   # strongest -> weakest
    bp4_t <- bp4_t[order(thr[bp4_t], decreasing = TRUE)]   # weakest  -> strongest

    rows <- list()
    # Pathogenic tiers, strongest first: [own cut-off, next-stronger cut-off)
    for (i in seq_along(pp3_t)) {
        tier <- pp3_t[i]
        rows[[length(rows) + 1]] <- data.frame(
            EVIDENCE = tier_label[[tier]],
            POINTS   = tier_points[[tier]],
            lo = thr[[tier]], lo_inc = TRUE,
            hi = if (i == 1) rng[2] else thr[[pp3_t[i - 1]]], hi_inc = (i == 1),
            stringsAsFactors = FALSE)
    }
    # The gap between the two Supporting cut-offs: no evidence either way.
    rows[[length(rows) + 1]] <- data.frame(
        EVIDENCE = "Indeterminate", POINTS = 0,
        lo = thr[["bp4"]], lo_inc = FALSE,
        hi = thr[["pp3"]], hi_inc = FALSE, stringsAsFactors = FALSE)
    # Benign tiers, weakest first: (next-stronger cut-off, own cut-off]
    for (i in seq_along(bp4_t)) {
        tier <- bp4_t[i]
        rows[[length(rows) + 1]] <- data.frame(
            EVIDENCE = tier_label[[tier]],
            POINTS   = -tier_points[[tier]],
            lo = if (i == length(bp4_t)) rng[1] else thr[[bp4_t[i + 1]]],
            lo_inc = (i == length(bp4_t)),
            hi = thr[[tier]], hi_inc = TRUE, stringsAsFactors = FALSE)
    }
    tab <- bind_rows(rows)

    # Count each set's scored variants per interval, as "n (pct%)".
    count_in <- function(x, lo, hi, lo_inc, hi_inc) {
        x  <- x[!is.na(x)]
        ge <- if (lo_inc) x >= lo else x > lo
        le <- if (hi_inc) x <= hi else x < hi
        sum(ge & le)
    }
    for (set_name in names(sets)) {
        x <- sets[[set_name]][[pred]]
        n_scored <- sum(!is.na(x))
        n <- vapply(seq_len(nrow(tab)),
                    function(i) count_in(x, tab$lo[i], tab$hi[i], tab$lo_inc[i], tab$hi_inc[i]),
                    numeric(1))
        tab[[set_name]] <- sprintf("%d (%.1f%%)", n, 100 * n / n_scored)
    }
    # Standard interval notation, written LOW -> HIGH, with the bracket type marking each
    # bound's true inclusivity: "[" / "]" closed, "(" / ")" open. So "score >= 0.972 and
    # < 0.990" prints as "[0.972, 0.990)" - exact against the raw scores counted above.
    tab$THRESHOLDS <- sprintf("%s%.3f, %.3f%s",
                              ifelse(tab$lo_inc, "[", "("), tab$lo,
                              tab$hi, ifelse(tab$hi_inc, "]", ")"))
    tab[, c("EVIDENCE", "POINTS", "THRESHOLDS", names(sets))]
}

# One column per set, in this order. The counts are of variants with a non-missing score,
# so each set's own denominator drives its percentages.
tier_dist_sets <- list(CONTROLS = dat_controls, VUS = dat_vus, gnomAD = gnomad_vars)

# Supplementary Table 4: AM
tier_dist_AM <- tier_interval_table("AM", tier_dist_sets)

# Supplementary Table 5: REVEL.
# REVEL's bp4_very_strong (<= 0.003) is a Pejaver 2022 tier only - the Bergquist 2025
# ladder reported in this table stops at Strong, so it is dropped here and BP4 Strong
# extends down to the bottom of the scale. It stays in `op_thresholds`, and is still
# drawn in Figures 1F / 1H.
tier_dist_REVEL <- tier_interval_table("REVEL", tier_dist_sets, drop_tiers = "bp4_very_strong")

write.xlsx(
    list(SupplTable4_AM    = tier_dist_AM,     # controls / VUS / gnomAD counts per AM tier
         SupplTable5_REVEL = tier_dist_REVEL), # controls / VUS / gnomAD counts per REVEL tier
    file = file.path(results.dir, "POLE_POLD1_tier_distribution.xlsx"),
    overwrite = TRUE)

# --- Console report ------------------------------------------------------
cat("\n", strrep("=", 73), "\n", sep = "")
cat("TIER-INTERVAL DISTRIBUTION OF THE CONTROL, VUS AND gnomAD VARIANTS\n")
cat(strrep("=", 73), "\n", sep = "")
for (m in c("AM", "REVEL")) {
    cat(sprintf("\nDenominators (%s): %s (variants with a score)\n", m,
                paste(sprintf("%s n=%d", names(tier_dist_sets),
                              vapply(tier_dist_sets, function(d) sum(!is.na(d[[m]])), integer(1))),
                      collapse = ", ")))
}
cat("\n--- Supplementary Table 4: AM ---\n");    print(tier_dist_AM,    row.names = FALSE)
cat("\n--- Supplementary Table 5: REVEL ---\n"); print(tier_dist_REVEL, row.names = FALSE)

cat("\nOutputs written to:", normalizePath(results.dir), "\n")
cat(strrep("=", 73), "\n", sep = "")


# Session information for reproducibility
cat("\n"); print(sessionInfo())
