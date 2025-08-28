# ==========================================================
# Limonene SOA: Paper vs SIRIUS - Tables, Colored Plots, Structures
# Robust matching, white BG, accuracy captions, structure rendering
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(rcdk)     # requires rcdklibs
  library(grid)     # for grid.raster()
})

# ----------------
# User parameters
# ----------------
sirius_csv      <- "C:/Users/91962/Predicting Chemicals/Real Sample/Sirius_predictionLimoneneNegative.csv"
sirius_project  <- "C:/Users/91962/Predicting Chemicals/LimoneneNegative"   # SIRIUS task folder (with Summaries)
out_dir         <- "C:/Users/91962/Predicting Chemicals/Real Sample/Structures"

# outputs
out_matches_csv <- file.path(out_dir, "limonene_paper_vs_sirius_matches.csv")
out_rollup_csv  <- file.path(out_dir, "limonene_per_formula_rollup.csv")
fig_bar_png     <- file.path(out_dir, "Fig_bar_median_prob_by_formula.png")
fig_hist_png    <- file.path(out_dir, "Fig_hist_prob_medians.png")
fig_rt_png      <- file.path(out_dir, "Fig_scatter_RT_vs_prob.png")
struct_dir      <- file.path(out_dir, "structures")  # structure PNGs

# ROC-optimised threshold from your model
THRESH <- 0.4805

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(struct_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------
# Model metrics to show on the plots
# ----------------------------------
acc_internal    <- 0.9025  # 90.25%
sens_internal   <- 0.9370  # 93.70%
spec_internal   <- 0.8951  # 89.51%
balacc_internal <- 0.9161  # 91.61%
acc_external    <- 0.7236  # 72.36%
acc_caption <- sprintf(
  "Model (thr = 0.4805): Internal Acc=%.2f, Sens=%.2f, Spec=%.2f, BalAcc=%.2f; External spectral Acc=%.2f",
  acc_internal, sens_internal, spec_internal, balacc_internal, acc_external
)
cat(acc_caption, "\n")

# ---------------------------
# Helpers
# ---------------------------
canonicalize_formula <- function(x) {
  x <- toupper(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("\\[.*?\\]", "", x)     # strip [M+H]+ etc.
  x <- gsub("[^A-Z0-9]", "", x)     # keep letters/digits only
  x
}

# ---------------------------
# 1) Reference table (paper)
# ---------------------------
ref_df <- tibble::tribble(
  ~Reference_compound, ~Elemental_formula, ~Nominal_mass_Da, ~Retention_time_min,
  "Aromatic compound",                                  "C11H10O2", 174,  4.2,
  "Aromatic compound",                                  "C11H10O2", 174,  4.6,
  "Limononic acid",                                     "C10H16O3", 184, 27.7,
  "Limonic acid",                                        "C9H14O4", 186, 22.5,
  "Ketolimononic acid",                                  "C9H14O4", 186, 12.2,
  "Ketolimonic acid",                                    "C8H12O5", 188,  6.8,
  "Ketolimonic acid (isomer)",                           "C8H12O5", 188,  7.4,
  "Dicarbonyl derivative of limononic acid",            "C10H14O4", 198, 17.9,
  "Dicarbonyl derivative of limononic acid",            "C10H14O4", 198, 23.1,
  "2-Hydroxy limononic acid (tentative)",               "C10H16O4", 200, 13.0,
  "Non-acidic isomer of 7-hydroxy limononic acid",      "C10H16O4", 200, 15.9,
  "Non-acidic isomer of 7-hydroxy limononic acid",      "C10H16O4", 200, 18.7,
  "7-hydroxy limononic acid",                           "C10H16O4", 200, 19.5,
  "2-hydroperoxy limononic acid",                       "C10H16O5", 216, 13.2,
  "Derivative of limononic acid",                       "C10H18O6", 234, 25.2,
  "Oligomer: 7-hydroxy limonoaldehyde + norlimonoaldehyde","C19H30O5",338,36.3,
  "Oligomer: keto-limononic + norlimonoaldehyde",       "C18H28O6", 340, 31.3,
  "Oligomer: keto-limononic + norlimonoaldehyde",       "C18H28O6", 340, 32.0,
  "Oligomer: limonic + 7-hydroxy limononic acid (iso 1)","C19H28O7",368, 33.4,
  "Oligomer: limonic + 7-hydroxy limononic acid (iso 2)","C19H28O7",368, 33.9,
  "Oligomer: limonic + 7-hydroxy limononic acid (iso 3)","C19H28O7",368, 34.9
) %>%
  mutate(
    Elemental_formula = str_trim(Elemental_formula),
    formula_clean     = canonicalize_formula(Elemental_formula)
  )

# ---------------------------------
# 2) Read your SIRIUS predictions
# ---------------------------------
sirius <- read_csv(sirius_csv, show_col_types = FALSE) %>%
  rename_with(~str_trim(.x)) %>%
  mutate(
    predform   = str_trim(as.character(predform)),
    predform_clean = canonicalize_formula(predform),
    Carcinogenicity_predicted = suppressWarnings(as.numeric(Carcinogenicity_predicted)),
    exactMass  = suppressWarnings(as.numeric(exactMass))
  )

if (!"predform" %in% names(sirius))  message("Warning: 'predform' column not found.")
if (!"Carcinogenicity_predicted" %in% names(sirius)) message("Warning: 'Carcinogenicity_predicted' column not found.")
if (!"predion" %in% names(sirius))   message("Note: 'predion' column not found (example ions will be blank).")
if (!"exactMass" %in% names(sirius)) message("Note: 'exactMass' column not found (mass range will be blank).")

# --------------------------------------------
# 3) Match by formula (canonical) and summarise per entry
# --------------------------------------------
summarise_one <- function(formula_clean, rt_min, name) {
  sub <- sirius %>% filter(predform_clean == formula_clean)
  n <- nrow(sub)
  
  if (n == 0) {
    tibble(
      Reference_compound = name,
      formula_clean = formula_clean,
      Retention_time_min = rt_min,
      SIRIUS_matches = 0L,
      Pred_prob_min = NA_real_,
      Pred_prob_median = NA_real_,
      Pred_prob_max = NA_real_,
      Percent_predicted_Active_at_0.4805 = NA_real_,
      Example_predion_values = "",
      ExactMass_range_in_SIRIUS = NA_character_
    )
  } else {
    probs <- sub$Carcinogenicity_predicted
    prob_min <- suppressWarnings(min(probs, na.rm = TRUE))
    prob_med <- suppressWarnings(median(probs, na.rm = TRUE))
    prob_max <- suppressWarnings(max(probs, na.rm = TRUE))
    if (!is.finite(prob_min)) prob_min <- NA_real_
    if (!is.finite(prob_med)) prob_med <- NA_real_
    if (!is.finite(prob_max)) prob_max <- NA_real_
    
    actives <- mean(probs >= THRESH, na.rm = TRUE) * 100
    if (!is.finite(actives)) actives <- NA_real_
    
    ions <- if ("predion" %in% names(sub)) paste(head(unique(na.omit(sub$predion)), 5), collapse = "; ") else ""
    
    mass_range <- if ("exactMass" %in% names(sub) && any(is.finite(sub$exactMass))) {
      rng <- range(sub$exactMass, na.rm = TRUE)
      sprintf("%.4f-%.4f", rng[1], rng[2])
    } else NA_character_
    
    tibble(
      Reference_compound = name,
      formula_clean = formula_clean,
      Retention_time_min = rt_min,
      SIRIUS_matches = n,
      Pred_prob_min = prob_min,
      Pred_prob_median = prob_med,
      Pred_prob_max = prob_max,
      Percent_predicted_Active_at_0.4805 = actives,
      Example_predion_values = ions,
      ExactMass_range_in_SIRIUS = mass_range
    )
  }
}

matches_df <- purrr::pmap_dfr(
  list(ref_df$formula_clean, ref_df$Retention_time_min, ref_df$Reference_compound),
  summarise_one
) %>%
  # Join back the pretty formula + nominal mass
  left_join(
    ref_df %>% select(formula_clean, Elemental_formula, Nominal_mass_Da, Reference_compound, Retention_time_min),
    by = c("formula_clean", "Reference_compound", "Retention_time_min")
  ) %>%
  relocate(Elemental_formula, .after = formula_clean) %>%
  relocate(Nominal_mass_Da, .after = Elemental_formula)

# -------------------------------------
# 4) Per-formula rollup (compact view)
# -------------------------------------
rollup <- matches_df %>%
  group_by(Elemental_formula) %>%
  summarise(
    n_rows = n(),
    total_SIRIUS_matches = sum(SIRIUS_matches, na.rm = TRUE),
    median_prob = suppressWarnings(median(Pred_prob_median, na.rm = TRUE))
  ) %>%
  arrange(desc(median_prob)) %>%
  ungroup()

# --------------------
# 5) Write to files
# --------------------
write_csv(matches_df, out_matches_csv)
write_csv(rollup, out_rollup_csv)
message("Done. CSVs written:\n - ", out_matches_csv, "\n - ", out_rollup_csv)

# --------------------
# 6) Visualisations (colored, high-DPI, white bg)
# --------------------

base_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# Color settings
grad_cols  <- c("#D0F0FD", "#74C0E3", "#1F78B4", "#08306B")           # bar fill gradient
class_cols <- c("< threshold" = "#9E9E9E", "??? threshold" = "#D81B60")  # scatter point colors
thresh_col <- "#E63946"

# Guard: if all median_prob are NA, avoid blank plot by setting to 0
if (all(!is.finite(rollup$median_prob))) {
  rollup$median_prob <- 0
}

# (a) Bar plot: median probability by formula
rollup_plot_df <- rollup %>%
  filter(is.finite(median_prob)) %>%
  mutate(Elemental_formula = factor(Elemental_formula, levels = Elemental_formula[order(median_prob, decreasing = TRUE)]))

ymax_bar <- max(0.05, max(rollup_plot_df$median_prob, na.rm = TRUE) * 1.10)

p_bar <- ggplot(rollup_plot_df, aes(x = Elemental_formula, y = median_prob, fill = median_prob)) +
  geom_col(width = 0.7, color = "white") +
  scale_fill_gradientn(colors = grad_cols, name = "Median\nprob.") +
  geom_text(aes(label = sprintf("%.2f", median_prob)), hjust = -0.1, size = 4, color = "#222222") +
  coord_flip(ylim = c(0, ymax_bar)) +
  labs(
    title = "Median predicted p53-activity probability by formula",
    x = "Elemental formula",
    y = "Median probability",
    caption = acc_caption
  ) + base_theme

ggsave(fig_bar_png, p_bar, width = 10, height = 8, dpi = 600, bg = "white")

# (b) Histogram: distribution of per-entry median probabilities
p_hist <- ggplot(matches_df, aes(x = Pred_prob_median)) +
  geom_histogram(binwidth = 0.05, boundary = 0, closed = "left",
                 fill = "#5E60CE", color = "white", na.rm = TRUE) +
  geom_vline(xintercept = THRESH, linetype = "dashed", linewidth = 1, color = thresh_col) +
  labs(
    title = "Distribution of predicted probabilities (per paper entry)",
    x = "Median probability for paper entry",
    y = "Count",
    caption = acc_caption
  ) + base_theme

ggsave(fig_hist_png, p_hist, width = 10, height = 6, dpi = 600, bg = "white")

# (c) Scatter: Retention time vs probability (colored by class vs threshold)
matches_col <- matches_df %>%
  mutate(Class_048 = ifelse(Pred_prob_median >= THRESH, "??? threshold", "< threshold"))

p_rt <- ggplot(matches_col, aes(x = Retention_time_min, y = Pred_prob_median, color = Class_048, label = Elemental_formula)) +
  geom_point(size = 3.5, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(values = class_cols, name = "Class @ 0.4805") +
  geom_hline(yintercept = THRESH, linetype = "dashed", linewidth = 1, color = thresh_col) +
  ggrepel::geom_text_repel(size = 4, min.segment.length = 0, max.overlaps = 20, color = "#1A1A1A", na.rm = TRUE) +
  labs(
    title = "Retention time vs predicted p53-activity probability",
    x = "Retention time (min)",
    y = "Median probability for paper entry",
    caption = acc_caption
  ) + base_theme

ggsave(fig_rt_png, p_rt, width = 10, height = 6, dpi = 600, bg = "white")

message("Figures saved:\n - ", fig_bar_png, "\n - ", fig_hist_png, "\n - ", fig_rt_png)

# ---------------------------------------------------
# 7) Read + standardize SIRIUS candidates; draw structures (white bg, robust draw)
# ---------------------------------------------------

read_candidates_safely <- function(path) {
  df <- tryCatch(
    readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()), na = c("", "NA")),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  names(df) <- tolower(names(df))
  if (!"formula" %in% names(df) && "molecularformula" %in% names(df)) df <- dplyr::rename(df, formula = molecularformula)
  if (!"smiles"  %in% names(df) && "smiles2d"        %in% names(df)) df <- dplyr::rename(df, smiles  = smiles2d)
  if (!"name"    %in% names(df) && "candidate"       %in% names(df)) df <- dplyr::rename(df, name    = candidate)
  
  req <- c("formula","adduct","rank","score","name","smiles","inchi","inchikey")
  for (r in req) if (!r %in% names(df)) df[[r]] <- NA_character_
  
  df %>%
    mutate(
      rank  = suppressWarnings(as.numeric(rank)),
      score = suppressWarnings(as.numeric(score))
    ) %>%
    select(any_of(c("formula","adduct","rank","score","name","smiles","inchi","inchikey"))) %>%
    mutate(source_file = path)
}

cand_files <- list.files(sirius_project, pattern = "candidates.tsv$", recursive = TRUE, full.names = TRUE)
cands <- purrr::map_dfr(cand_files, read_candidates_safely)

if (nrow(cands) == 0) {
  message("No SIRIUS candidate files found or readable. Skipping structure drawing.")
} else {
  cands <- cands %>% mutate(formula_clean = canonicalize_formula(formula))
  
  # top per (formula, source file) - no slice(), dplyr-version-safe
  top_per_file <- cands %>%
    filter(!is.na(smiles), !is.na(formula_clean)) %>%
    group_by(formula_clean, source_file) %>%
    arrange(rank, desc(score), .by_group = TRUE) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()
  
  # choose best per formula among files
  formulas_of_interest <- unique(ref_df$formula_clean)
  chosen <- top_per_file %>%
    filter(formula_clean %in% formulas_of_interest) %>%
    group_by(formula_clean) %>%
    arrange(rank, desc(score), .by_group = TRUE) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()
  
  # robust depiction that actually renders into PNG
  depict_and_save <- function(smiles, outfile, width = 900, height = 600, res = 150) {
    mol <- try(rcdk::parse.smiles(smiles)[[1]], silent = TRUE)
    if (inherits(mol, "try-error") || is.null(mol)) return(FALSE)
    
    # prepare molecule: coords + aromaticity
    try(rcdk::convert.implicit.to.explicit(mol), silent = TRUE)
    try(rcdk::do.aromaticity(mol),                silent = TRUE)
    try(rcdk::generate.2d.coordinates(mol),      silent = TRUE)
    
    # Try to get a raster
    img <- try(rcdk::view.image.2d(mol, width = width, height = height), silent = TRUE)
    
    png(outfile, width = width, height = height, res = res, bg = "white")
    grid::grid.newpage()
    if (!inherits(img, "try-error") && !is.null(img)) {
      grid::grid.raster(img, interpolate = TRUE)
    } else {
      # Fallback: draw directly onto the device
      par(bg = "white", mar = c(1,1,1,1))
      rcdk::view.image.2d(mol)
    }
    dev.off()
    TRUE
  }
  
  for (i in seq_len(nrow(chosen))) {
    form_clean <- chosen$formula_clean[i]
    smi        <- chosen$smiles[i]
    nm         <- ifelse(!is.na(chosen$name[i]) && nchar(chosen$name[i]) > 0, chosen$name[i], "TopCandidate")
    # Use pretty formula from ref_df if available
    pretty_form <- ref_df$Elemental_formula[match(form_clean, ref_df$formula_clean)]
    if (is.na(pretty_form)) pretty_form <- form_clean
    
    file_safe <- gsub("[^A-Za-z0-9_\\-]+", "_", paste0(pretty_form, "_", nm))
    outfile   <- file.path(struct_dir, paste0("struct_", file_safe, ".png"))
    ok <- depict_and_save(smi, outfile)
    message(if (ok) paste("Saved structure:", outfile) else paste("Failed to depict:", pretty_form))
  }
  message("Structure images saved in: ", struct_dir)
}

# ----------------------------
# 8) Quick console preview
# ----------------------------
print(head(matches_df, 10))
print(rollup)
