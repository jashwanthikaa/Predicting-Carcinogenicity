# ===============================
# Limonene SOA - SIRIUS vs Paper
# With output directory
# ===============================

# ---- 0) Settings ----
# Change these two paths as needed:
in_csv  <- "C:/Users/91962/Predicting Chemicals/Real Sample/Sirius_predictionLimoneneNegative.csv"
out_dir <- "C:/Users/91962/Predicting Chemicals/Real Sample/"   # e.g., "C:/Users/91962/Predicting Chemicals/outputs_limonene"

# Create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 1) Packages ----
pkgs <- c("tidyverse")
to_install <- pkgs[!suppressWarnings(sapply(pkgs, requireNamespace, quietly = TRUE))]
if (length(to_install)) install.packages(to_install)
library(tidyverse)

# ---- 2) Load SIRIUS CSV ----
df <- readr::read_csv(in_csv, show_col_types = FALSE)

# Standardise column names used below
df <- df %>%
  rename(
    Formula         = dplyr::any_of("predform"),
    Ion             = dplyr::any_of("predion"),
    Carcinogenicity = dplyr::any_of("Carcinogenicity_predicted")
  )

# Basic checks
req_cols <- c("Formula","Ion","Carcinogenicity","foldernumber")
missing_cols <- setdiff(req_cols, names(df))
if (length(missing_cols)) {
  stop("Missing required columns in input CSV: ", paste(missing_cols, collapse = ", "))
}

# ---- 3) Literature (paper) anchor formulas ----
literature_map <- tibble::tribble(
  ~Formula,    ~literature_name,
  "C9H14O4",   "Limonic/limononic acid",
  "C10H16O4",  "7-hydroxy-limononic acid",
  "C9H14O5",   "Carbonyl related to hydroxy-limononic",
  "C10H16O3",  "Limononaldehyde (candidate)",
  "C8H12O4",   "Keto-limonic acid",
  "C8H12O5",   "7-hydroxy-keto-limononic acid",
  "C18H28O6",  "Oligomer (ester/aldol/hemiacetal)",
  "C19H28O7",  "Oligomer (ester/aldol/hemiacetal)",
  "C19H30O5",  "Oligomer (ester/aldol/hemiacetal)",
  "C19H30O7",  "Oligomer (ester/aldol/hemiacetal)",
  "C20H34O9",  "Oligomer (ester/aldol/hemiacetal)"
)

# ---- 4) Build summary table ----
join_df <- df %>%
  filter(Formula %in% literature_map$Formula) %>%
  left_join(literature_map, by = "Formula")

join_df$Ion <- as.character(join_df$Ion)

summary_tbl <- join_df %>%
  group_by(Formula, literature_name) %>%
  summarise(
    n_detections  = n(),
    n_active      = sum(Carcinogenicity == "Active",   na.rm = TRUE),
    n_inactive    = sum(Carcinogenicity == "Inactive", na.rm = TRUE),
    example_foldernumbers = paste(unique(foldernumber)[seq_len(min(5, dplyr::n()))],
                                  collapse = ", "),
    ions          = paste(unique(na.omit(Ion)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(`Found in file?` = if_else(n_detections > 0, "Yes", "No"))

# Save table to output directory
tbl_path <- file.path(out_dir, "Literature_vs_SIRIUS_summary.csv")
readr::write_csv(summary_tbl, tbl_path)
print(summary_tbl)
message("Saved table: ", tbl_path)

# ---- 5) Plot: detections per literature formula ----
plot1_df <- summary_tbl %>% filter(n_detections > 0) %>%
  arrange(desc(n_detections)) %>%
  mutate(Formula = factor(Formula, levels = unique(Formula)))

p1 <- ggplot(plot1_df, aes(x = Formula, y = n_detections)) +
  geom_col(fill = "#3182bd") +
  labs(
    title = "Detections per literature formula (found in SIRIUS)",
    x = "Formula", y = "Count of detections"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1_path <- file.path(out_dir, "plot_counts_by_formula.png")
ggsave(p1_path, p1, width = 9, height = 5, dpi = 300)
print(p1)
message("Saved figure: ", p1_path)

# ---- 6) Plot: stacked Active vs Inactive per formula ----
stacked <- summary_tbl %>%
  filter(n_detections > 0) %>%
  transmute(
    Formula,
    Active   = n_active,
    Inactive = n_inactive
  ) %>%
  pivot_longer(-Formula, names_to = "Prediction", values_to = "Count") %>%
  mutate(
    Formula    = factor(Formula, levels = unique(Formula)),
    Prediction = factor(Prediction, levels = c("Inactive","Active"))
  )

p2 <- ggplot(stacked, aes(x = Formula, y = Count, fill = Prediction)) +
  geom_col() +
  scale_fill_manual(values = c("Inactive" = "#3182bd", "Active" = "#e6550d")) +
  labs(
    title = "Carcinogenicity predictions by formula (SIRIUS features)",
    x = "Formula", y = "Count of detections"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2_path <- file.path(out_dir, "plot_active_inactive_by_formula.png")
ggsave(p2_path, p2, width = 9, height = 5, dpi = 300)
print(p2)
message("Saved figure: ", p2_path)

# ---- 7) OPTIONAL: Same formula, different RT (runs only if RT column exists) ----
rt_col <- names(df)[grepl("^rt$|retention", tolower(names(df)))]
if (length(rt_col)) {
  rt_col <- rt_col[1]
  key_formulas <- c("C10H16O4","C9H14O4","C19H28O7")
  df_rt <- df %>%
    filter(Formula %in% key_formulas) %>%
    mutate(rt_min = suppressWarnings(as.numeric(.data[[rt_col]]))) %>%
    filter(!is.na(rt_min))
  
  if (nrow(df_rt) > 0) {
    p3 <- ggplot(df_rt, aes(x = rt_min, y = factor(Formula, levels = key_formulas))) +
      geom_point(aes(shape = Carcinogenicity), size = 2, alpha = 0.9) +
      labs(
        title = "Same formula, different retention times (isomeric features)",
        x = "Retention time (min)", y = "Formula"
      ) +
      theme_classic(base_size = 14)
    
    p3_path <- file.path(out_dir, "plot_same_formula_different_rt.png")
    ggsave(p3_path, p3, width = 9, height = 4, dpi = 300)
    print(p3)
    message("Saved figure: ", p3_path)
  } else {
    message("RT column detected (", rt_col, ") but no numeric RT values available.")
  }
} else {
  message("No RT/retention time column found - skipping the RT plot.")
}

message("\nAll outputs are in: ", normalizePath(out_dir, winslash = "/"))
