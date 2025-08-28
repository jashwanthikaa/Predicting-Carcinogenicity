# ============================================================
# FULL PIPELINE: Align features -> predict -> SHAP -> Explanation figure
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(caret)       # ok even if model is a raw Booster
  library(xgboost)
  library(shapviz)
  library(patchwork)
  library(cowplot)
  library(httr)
})

# -------------------- 1) Paths & settings --------------------
in_xlsx  <- "C:/Users/91962/Predicting Chemicals/FinalData (1).xlsx"
in_model <- "C:/Users/91962/Predicting Chemicals/xgb_final_model_trained_on_all_dataws.rds"
out_dir  <- "C:/Users/91962/Predicting Chemicals/fig_explain_final"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Which row to explain: "maxprob" or numeric index (e.g., 715)
row_pick <- "maxprob"
# How many features to show explicitly in the waterfall & cards
waterfall_k <- 9

# Optional: feature???SMILES/SMARTS mapping for bottom cards
use_feature_mapping <- FALSE
mapping_sheet       <- "Mapping"  # requires columns: Feature, Label (SMILES/SMARTS)

# -------------------- 2) Load data --------------------
df <- read_excel(in_xlsx)

# Keep a copy of SMILES aligned to rows we actually use later
if (!"SMILES" %in% names(df)) stop("SMILES column not found in workbook.")
smiles_all <- df$SMILES

# Start with ALL numeric feature columns. We'll subset to model features later.
# (This avoids guessing patterns like ^Un or exactMass)
is_num <- sapply(df, function(z) is.numeric(z) || is.integer(z))
X_all  <- df[, is_num, drop = FALSE]

# Remove rows with any NA in these numeric features (to mirror usual training cleaning)
non_na_rows <- complete.cases(X_all)
X_all   <- X_all[non_na_rows, , drop = FALSE]
smiles_all <- smiles_all[non_na_rows]
df_clean   <- df[non_na_rows, , drop = FALSE]

# -------------------- 3) Load model & figure out expected features --------------------
model_obj <- readRDS(in_model)

is_caret <- inherits(model_obj, "train")
if (is_caret) {
  booster <- model_obj$finalModel
  # Most caret xgb models retain feature names here:
  model_feats <- booster$feature_names
  if (is.null(model_feats)) {
    # Fallback: try to infer from trainingData colnames minus outcome column
    trn <- model_obj$trainingData
    model_feats <- setdiff(colnames(trn), ".outcome")
  }
} else if (inherits(model_obj, "xgb.Booster")) {
  booster <- model_obj
  model_feats <- booster$feature_names
} else {
  stop("Unknown model type in RDS. Expected caret::train or xgb.Booster.")
}

if (is.null(model_feats) || length(model_feats) == 0) {
  stop("Could not extract feature names from the model. The model must be trained with column names.")
}

# Intersect current data columns with model features & order exactly as in the model
missing_in_data <- setdiff(model_feats, colnames(X_all))
if (length(missing_in_data) > 0) {
  stop(sprintf(
    "Your Excel lacks %d model feature(s), e.g.: %s",
    length(missing_in_data), paste(head(missing_in_data, 10), collapse = ", ")
  ))
}

X <- X_all[, model_feats, drop = FALSE]

# -------------------- 4) Predict probabilities --------------------
# caret::train uses predict(..., type="prob"); raw Booster uses xgboost::predict
if (is_caret) {
  pred_prob <- predict(model_obj, newdata = as.matrix(X), type = "prob")[, "Active"]
} else {
  # For Booster, we must ensure matrix with colnames==model_feats
  pred_prob <- predict(booster, as.matrix(X))
  # If it's binary: xgb returns P(class=1). Good.
}

# -------------------- 5) SHAP with shapviz --------------------
sv <- shapviz(
  booster,
  X      = as.matrix(X),
  X_pred = as.matrix(X),
  feature_names = model_feats
)
S_mat <- sv$S
stopifnot(nrow(S_mat) > 0, ncol(S_mat) > 0)

# -------------------- 6) Pick row to explain --------------------
row_id <- if (identical(row_pick, "maxprob")) which.max(pred_prob) else as.integer(row_pick)
if (is.na(row_id) || row_id < 1 || row_id > nrow(X)) stop("row_id out of range for current data.")

smiles_row <- smiles_all[row_id]
if (is.na(smiles_row) || !nzchar(smiles_row)) smiles_row <- "C"  # harmless fallback (methane)

# -------------------- 7) Hardened SMILES -> PNG depiction --------------------
depict_smiles <- function(smi, out_png, w = 900, h = 900, label = NULL, timeout_sec = 12) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  enc <- utils::URLencode(smi, reserved = TRUE)
  
  # Try NCI Cactus first
  url1 <- sprintf("https://cactus.nci.nih.gov/chemical/structure/%s/image?format=png&w=%d&h=%d",
                  enc, w, h)
  # SMARTS.plus fallback
  url2 <- sprintf("https://smarts.plus/render/mol/png?smi=%s&w=%d&h=%d", enc, w, h)
  
  for (u in c(url1, url2)) {
    resp <- try(httr::GET(u, httr::timeout(timeout_sec)), silent = TRUE)
    if (!(inherits(resp, "try-error") || httr::http_error(resp))) {
      ct <- httr::headers(resp)[["content-type"]]
      if (!is.null(ct) && grepl("^image/png", ct, ignore.case = TRUE)) {
        bin <- httr::content(resp, "raw")
        # Validate PNG signature
        if (length(bin) >= 8 && identical(unname(bin[1:8]),
                                          as.raw(c(0x89,0x50,0x4E,0x47,0x0D,0x0A,0x1A,0x0A)))) {
          writeBin(bin, out_png)
          return(out_png)
        }
      }
    }
  }
  
  # Safe placeholder (always valid PNG)
  png(out_png, width = w, height = h, res = 150, bg = "white")
  par(mar = c(0,0,0,0)); plot.new()
  text(0.5, 0.6, "Couldn't render SMILES", cex = 1.2)
  lab <- if (is.null(label)) smi else label
  text(0.5, 0.4, substr(lab, 1, 70), cex = 0.9)
  dev.off()
  out_png
}

# -------------------- 8) Left panel: compound structure --------------------
left_png <- file.path(out_dir, "panel_assets", sprintf("compound_row_%s.png", row_id))
depict_smiles(smiles_row, left_png, w = 1000, h = 1000)

left_panel <- ggdraw() +
  draw_image(left_png, x = 0.5, y = 0.56, scale = 0.92) +
  draw_label(
    sprintf("ROW #%s\nPred P(active)=%.3f\nSMILES:\n%s",
            row_id, pred_prob[row_id], substr(smiles_row, 1, 95)),
    x = 0.5, y = 0.06, size = 10, fontface = "bold", lineheight = 0.9
  )

# -------------------- 9) Waterfall (top-right) --------------------
p_water <- sv_waterfall(sv, row_id = row_id, max_display = waterfall_k) +
  ggtitle("Top feature contributions") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# -------------------- 10) Bottom feature cards --------------------
# Build mapping Feature -> SMILES/SMARTS if you have it
feat_to_smiles <- NULL
if (use_feature_mapping && mapping_sheet %in% readxl::excel_sheets(in_xlsx)) {
  mdf <- readxl::read_excel(in_xlsx, sheet = mapping_sheet)
  names(mdf) <- gsub("\\s+", "", names(mdf))
  if (all(c("Feature","Label") %in% names(mdf))) {
    feat_to_smiles <- setNames(trimws(mdf$Label), trimws(mdf$Feature))
  }
}

shaps_row <- as.numeric(S_mat[row_id, ])
names(shaps_row) <- colnames(S_mat)
top_feats <- names(sort(abs(shaps_row), decreasing = TRUE))[1:waterfall_k]

assets_dir <- file.path(out_dir, "panel_assets")
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

feat_cards <- vector("list", length(top_feats))
for (i in seq_along(top_feats)) {
  f <- top_feats[i]
  label_text <- paste0(f, "\nSHAP=", sprintf("%.3f", shaps_row[f]))
  
  have_map <- !is.null(feat_to_smiles) && !is.na(feat_to_smiles[f]) && nzchar(feat_to_smiles[f])
  
  if (have_map) {
    fpng <- file.path(assets_dir, sprintf("feat_%02d.png", i))
    depict_smiles(feat_to_smiles[f], fpng, w = 720, h = 520, label = f)
    g <- ggdraw() +
      draw_image(fpng, x = 0.5, y = 0.58, scale = 1) +
      draw_label(label_text, x = 0.5, y = 0.06, size = 9, fontface = "bold")
  } else {
    # Clean text card when no depictable mapping (e.g., Un#### / exactMass)
    g <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.62, label = "Fingerprint feature", size = 4.2) +
      annotate("text", x = 0.5, y = 0.42, label = label_text, size = 3.6, fontface = 2)
  }
  feat_cards[[i]] <- g
}
cards_panel <- wrap_plots(feat_cards, ncol = 3)

# -------------------- 11) Assemble & save --------------------
top_row <- left_panel + p_water + plot_layout(widths = c(1.1, 1.9))
final_panel <- (top_row) / cards_panel + plot_layout(heights = c(1.05, 1.0))

ggsave(file.path(out_dir, sprintf("Model_explanation_row_%s.png", row_id)),
       final_panel, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(out_dir, sprintf("Model_explanation_row_%s.pdf", row_id)),
       final_panel, width = 14, height = 12, device = cairo_pdf)

cat("\nSaved figure(s) in:", out_dir,
    "\nFile prefix: ", sprintf("Model_explanation_row_%s.*", row_id), "\n")
