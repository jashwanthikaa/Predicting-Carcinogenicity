# =============================================
# FINAL XGBOOST TRAINING ON TRAIN + TEST SET
# =============================================

# Load libraries
library(readxl)
library(caret)
library(dplyr)
library(xgboost)
library(Matrix)
library(pROC)

# ---------------------------------------------
# 1. Load and Prepare the Original Dataset
# ---------------------------------------------

# Load dataset
df <- read_excel("C:/Users/91962/Predicting Chemicals/FinalData (1).xlsx")

# Select fingerprint and exactMass columns
fp_cols <- grep("^Un", names(df), value = TRUE)
exactMass_col <- "exactMass"
X <- df %>% select(all_of(c(fp_cols, exactMass_col)))
y <- ifelse(tolower(df$`HIT.CALL`) == "active", 1, 0)

# Remove rows with missing values
non_na_rows <- complete.cases(X)
X <- X[non_na_rows, ]
y <- y[non_na_rows]
df_clean <- df[non_na_rows, ]

# Remove near-zero variance fingerprint features
nzv <- nearZeroVar(X)
if (length(nzv) > 0) {
  X <- X[, -nzv]
}

# ---------------------------------------------
# 2. Recreate Train/Test Split (for reference)
# ---------------------------------------------
set.seed(123)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
y_train <- y[train_index]
X_test <- X[-train_index, ]
y_test <- y[-train_index]

# ---------------------------------------------
# 3. Combine Train + Test Sets
# ---------------------------------------------
X_full <- rbind(X_train, X_test)
y_full <- c(y_train, y_test)
y_full_factor <- factor(ifelse(y_full == 1, "Active", "Inactive"))
X_full_matrix <- data.matrix(X_full)

# ---------------------------------------------
# 4. Class Imbalance Handling
# ---------------------------------------------
scale_weight <- sum(y_full == 0) / sum(y_full == 1)

# ---------------------------------------------
# 5. Training Control and Tuning Grid
# ---------------------------------------------
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

xgb_grid <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.01,
  gamma = 5,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  subsample = 0.6
)

# ---------------------------------------------
# 6. Train Final XGBoost Model
# ---------------------------------------------
set.seed(123)
xgb_final_model <- train(
  x = X_full_matrix,
  y = y_full_factor,
  trControl = train_control,
  tuneGrid = xgb_grid,
  method = "xgbTree",
  metric = "ROC",
  verbose = TRUE,
  scale_pos_weight = scale_weight
)

cat("\n Final Model Trained on All Data!\n")
print(xgb_final_model$bestTune)

# Save the model
saveRDS(xgb_final_model, "C:/Users/91962/Predicting Chemicals/xgb_final_model_trained_on_all_dataws.rds")

# ---------------------------------------------
# 7. Predict and Evaluate on Full Data
# ---------------------------------------------
xgb_pred_prob <- predict(xgb_final_model, newdata = X_full_matrix, type = "prob")[, "Active"]

# Find optimal threshold from ROC
roc_obj <- roc(y_full, xgb_pred_prob)
best_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))
optimal_threshold <- as.numeric(best_coords["threshold"])
cat("\n Optimal Threshold from ROC:", optimal_threshold, "\n")

# Final predictions using optimal threshold
xgb_pred <- ifelse(xgb_pred_prob > optimal_threshold, 1, 0)

# Confusion Matrix
y_full_factor_bin <- factor(y_full, levels = c(0, 1))
xgb_pred_factor_bin <- factor(xgb_pred, levels = c(0, 1))
xgb_cm <- confusionMatrix(xgb_pred_factor_bin, y_full_factor_bin, positive = "1")

cat("\n Confusion Matrix on Full Data:\n")
print(xgb_cm)

# Plot ROC Curve
png("ROC_curve_final_model.png", width = 800, height = 600)
plot(roc_obj, col = "blue", main = "ROC Curve - Final XGBoost Model (Train + Test Combined)")
dev.off()
cat("\n??? ROC Curve saved as 'ROC_curve_final_model6.png'\n")

# ---------------------------------------------
# 8. Save Predictions with Metadata
# ---------------------------------------------
metadata_cols <- c("Name", "CASRN", "SMILES", "exactMass", "DTXSID")
metadata_full <- df_clean[, metadata_cols]

results <- data.frame(
  Name = metadata_full$Name,
  DTXSID = metadata_full$DTXSID,
  CASRN = metadata_full$CASRN,
  SMILES = metadata_full$SMILES,
  exactMass = metadata_full$exactMass,
  Actual = y_full,
  Predicted = xgb_pred,
  Probability = xgb_pred_prob
)

write.csv(results, "C:/Users/91962/Predicting Chemicals/xgb_predictions_trained_on_allws.csv", row.names = FALSE)
cat("\n Predictions saved to CSV!\n")
-----------------------------------------------------------------------------------------
#Thershold Optimizations 
thresh_metrics <- coords(roc_obj, x = "all", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
df_thresh <- as.data.frame(thresh_metrics)
ggplot(df_thresh, aes(x = threshold)) +
  geom_line(aes(y = sensitivity, color = "Sensitivity")) +
  geom_line(aes(y = specificity, color = "Specificity")) +
  geom_line(aes(y = accuracy, color = "Accuracy")) +
  labs(title = "Threshold Optimization", y = "Metric Value") +
  scale_color_manual(values = c("red", "blue", "green"))
-------------------------------------------------------------------------------
# ================================
# Libraries
# ================================
library(readxl)
library(dplyr)
library(caret)
library(xgboost)
library(pROC)
library(ggplot2)
library(reshape2)
library(shapviz)       # SHAP for xgboost/caret
# For structure drawings (install once if needed):
# install.packages(c("rcdk","rcdklibs"), repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(rcdk)
  library(rcdklibs)
})
# ================================
# 1) Paths & output
# ================================
in_xlsx  <- "C:/Users/91962/Predicting Chemicals/FinalData (1).xlsx"
in_model <- "C:/Users/91962/Predicting Chemicals/xgb_final_model_trained_on_all_dataws.rds"
out_dir  <- "C:/Users/91962/Predicting Chemicals/figs_final"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ================================
# 2) Load data (mirror your training prep)
# ================================
df <- read_excel(in_xlsx)

fp_cols <- grep("^Un", names(df), value = TRUE) # fingerprint columns
exactMass_col <- "exactMass"

X <- df %>% dplyr::select(all_of(c(fp_cols, exactMass_col)))
y <- ifelse(tolower(df$`HIT.CALL`) == "active", 1, 0)

# Clean NAs exactly the same way
non_na_rows <- complete.cases(X)
X <- X[non_na_rows, , drop = FALSE]
y <- y[non_na_rows]
df_clean <- df[non_na_rows, ]

# Remove near-zero-variance features
nzv <- caret::nearZeroVar(X)
if (length(nzv) > 0) X <- X[, -nzv, drop = FALSE]

feature_names <- colnames(X)
X_mat <- data.matrix(X)

# ================================
# 3) Load trained model & predictions
# ================================
xgb_final_model <- readRDS(in_model)          # caret train object
xgb_native <- xgb_final_model$finalModel      # xgboost booster

pred_prob <- predict(xgb_final_model, newdata = X_mat, type = "prob")[, "Active"]

# ROC, threshold, binary preds
roc_obj <- pROC::roc(y, pred_prob)
best_coords <- pROC::coords(roc_obj, "best",
                            ret = c("threshold", "specificity", "sensitivity"))
opt_thr <- as.numeric(best_coords["threshold"])
pred_bin <- ifelse(pred_prob > opt_thr, 1, 0)

# Confusion matrix
cm <- caret::confusionMatrix(factor(pred_bin, levels = c(0,1)),
                             factor(y, levels = c(0,1)), positive = "1")
print(cm)

# ================================
# 4) Save ROC (white background)
# ================================
png(file.path(out_dir, "ROC_curve_final_model.png"),
    width = 1800, height = 1400, res = 300, bg = "white")
plot(roc_obj, col = "blue",
     main = "ROC Curve - Final XGBoost Model (Train+Test Combined)")
abline(a = 0, b = 1, lty = 2)
legend("bottomright",
       legend = paste0("AUC = ", round(pROC::auc(roc_obj), 3)),
       bty = "n")
dev.off()

# ================================
# 5) Confusion-matrix heatmap
# ================================
cm_df <- as.data.frame(cm$table)
colnames(cm_df) <- c("Pred", "True", "Freq")
p_cm <- ggplot(cm_df, aes(Pred, True, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 6) +
  scale_fill_gradient(low = "#c7e9c0", high = "#005a32") +
  labs(title = "Confusion Matrix (Final Model)", x = "Predicted", y = "Actual") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "Confusion_Matrix_heatmap.png"),
       p_cm, width = 7, height = 6, dpi = 300, bg = "white")

# ================================
# 6) SHAP values via shapviz
# ================================
# IMPORTANT: pass both X and X_pred
sv <- shapviz(
  xgb_native,
  X = X_mat,           # background data
  X_pred = X_mat,      # compute SHAP for these rows
  feature_names = feature_names
)

# Quick sanity check
stopifnot(nrow(sv$S) > 0, ncol(sv$S) > 0)

# --- (a) Bar plot of mean(|SHAP|) ---
p_bar <- sv_importance(sv, kind = "bar")
ggsave(file.path(out_dir, "SHAP_bar_mean_abs.png"),
       p_bar, width = 7, height = 6, dpi = 300, bg = "white")

# --- (b) Beeswarm (summary) plot ---
p_bee <- sv_importance(sv, kind = "bee", max_display = 20)
ggsave(file.path(out_dir, "SHAP_beeswarm_top20.png"),
       p_bee, width = 8, height = 7, dpi = 300, bg = "white")

# --- (c) Dependence plots for top features ---
top_feats_idx <- head(order(colMeans(abs(sv$S)), decreasing = TRUE), 8)
top_feats <- feature_names[top_feats_idx]
dep_dir <- file.path(out_dir, "SHAP_dependence")
if (!dir.exists(dep_dir)) dir.create(dep_dir)
for (f in top_feats) {
  p_dep <- sv_dependence(sv, v = f, color_var = f)
  ggsave(file.path(dep_dir, paste0("DEP_", gsub("[^A-Za-z0-9_]+","_", f), ".png")),
         p_dep, width = 6.5, height = 5.5, dpi = 300, bg = "white")
}

# --- (c2) SHAP dependence panel: top 15 in one figure -----------------
# Requires patchwork for easy layout
# install.packages("patchwork")  # run once if needed
library(patchwork)

top15 <- names(sort(colMeans(abs(S_mat)), decreasing = TRUE))[1:15]

dep_plots <- lapply(top15, function(f) {
  sv_dependence(sv, v = f, color_var = f) +
    ggtitle(f) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8)
    )
})

panel <- wrap_plots(dep_plots, ncol = 3)

# Save high-res PNG (white background)
ggsave(
  filename = file.path(out_dir, "SV_dependence_top15_panel.png"),
  plot = panel, width = 14, height = 16, dpi = 300, bg = "white"
)

# Also save a PDF (vector)
ggsave(
  filename = file.path(out_dir, "SV_dependence_top15_panel.pdf"),
  plot = panel, width = 14, height = 16, dpi = 300, device = cairo_pdf
)

# --- (d) Per-class mean |SHAP| bars (fixed implementation) ---
S_mat <- sv$S
cls <- factor(ifelse(y == 1, "Active", "Inactive"))
mean_abs_by_class_list <- lapply(split(as.data.frame(S_mat), cls), function(d) {
  colMeans(abs(as.matrix(d)), na.rm = TRUE)
})
mean_abs_by_class_mat <- do.call(rbind, mean_abs_by_class_list)   # rows = classes
mean_abs_by_class_df <- cbind(
  Class = rownames(mean_abs_by_class_mat),
  as.data.frame(mean_abs_by_class_mat, check.names = FALSE)
)
mean_abs_by_class_long <- reshape2::melt(mean_abs_by_class_df, id.vars = "Class",
                                         variable.name = "Feature",
                                         value.name = "MeanAbsSHAP")
global_top <- names(sort(colMeans(abs(S_mat)), decreasing = TRUE))[1:20]
mean_abs_by_class_long <- subset(mean_abs_by_class_long, Feature %in% global_top)
p_cls <- ggplot(mean_abs_by_class_long,
                aes(x = reorder(Feature, MeanAbsSHAP, FUN = mean),
                    y = MeanAbsSHAP, fill = Class)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = NULL, y = "mean(|SHAP|)",
       title = "Per-class mean(|SHAP|) for Top 20 Features") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "SHAP_per_class_bar_top20.png"),
       p_cls, width = 8.5, height = 7, dpi = 300, bg = "white")

# --- (e) Waterfall for a single observation ---
row_id <- which.max(pred_prob)  # choose any row you want to showcase
p_wf <- sv_waterfall(sv, row_id = row_id, max_display = 10)
ggsave(file.path(out_dir, sprintf("Waterfall_row_%s.png", row_id)),
       p_wf, width = 7.5, height = 6, dpi = 300, bg = "white")




# ================================
# 7) Export top-20 SHAP table
# ================================
imp_tbl <- sort(colMeans(abs(S_mat)), decreasing = TRUE)[1:20]
imp_df <- data.frame(
  Feature = names(imp_tbl),
  MeanAbsSHAP = as.numeric(imp_tbl),
  stringsAsFactors = FALSE
)
write.csv(imp_df, file.path(out_dir, "Top20_SHAP_features.csv"), row.names = FALSE)

cat("\nAll figures and tables saved to:\n", out_dir, "\n")

---------------------------------------------------------------------------
TOP_N <- 25  # try 20-30 for papers

imp_all <- sort(colMeans(abs(S_mat)), decreasing = TRUE)
top_vars <- names(imp_all)[1:TOP_N]

# optional: pretty/short names for labels
shorten <- function(x, width = 28) {
  x <- gsub("^Un", "Bit ", x)                 # 'Un357' -> 'Bit 357'
  x <- gsub("exactMass", "Exact mass", x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width-1), "."), x)
}

-------------------------------------------------------------------------------
TOP_N <- 25
imp_all <- sort(colMeans(abs(S_mat)), decreasing = TRUE)
top_vars <- names(imp_all)[1:TOP_N]

# shorten function
shorten <- function(x, width = 28) {
  x <- gsub("^Un", "Bit ", x)
  x <- gsub("exactMass", "Exact mass", x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width-1), "."), x)
}

# Subset shapviz object to top features only
sv_top <- sv
sv_top$S <- sv$S[, top_vars, drop = FALSE]
sv_top$X <- sv$X[, top_vars, drop = FALSE]

# Beeswarm plot
p_bee <- sv_importance(sv_top, kind = "bee", max_display = TOP_N) +
  ggtitle(sprintf("SHAP Summary (Top %d Features)", TOP_N)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(10, 14, 10, 14)
  )

# Make sure y-axis uses shortened names
p_bee <- p_bee + scale_y_discrete(
  limits = rev(top_vars),  # preserve order
  labels = shorten(rev(top_vars))
)

ggsave(file.path(out_dir, sprintf("SHAP_beeswarm_top%d.png", TOP_N)),
       p_bee, width = 12, height = 9, dpi = 350, bg = "white")


---------------------------------------------------------------------------------
library(readxl)
library(stringr)

# --- 0) Where to write the auto map
auto_map_csv <- file.path(out_dir, "feature_labels_auto.csv")

# --- 1) Helper to normalise names
normcol <- function(x) tolower(gsub("[^a-z0-9]+", "_", x))

# --- 2) Scan all sheets for a mapping-like table
build_feature_label_map <- function(xlsx_path) {
  sheets <- excel_sheets(xlsx_path)
  candidates <- list()
  
  for (sh in sheets) {
    df_sh <- suppressMessages(read_excel(xlsx_path, sheet = sh))
    if (nrow(df_sh) == 0) next
    colnames(df_sh) <- normcol(colnames(df_sh))
    
    # Likely feature-id column: values like Un1234, Un_123 etc., or exactMass
    feat_cols <- names(df_sh)[
      sapply(df_sh, function(col) {
        v <- as.character(col)
        any(grepl("^\\s*un\\s*[_-]?\\d+\\s*$", v, ignore.case = TRUE)) ||
          any(grepl("^\\s*exact\\s*mass\\s*$", v, ignore.case = TRUE)) ||
          "un" %in% unique(substr(v, 1, 2))
      })
    ]
    
    # Add obvious name variants too
    feat_cols <- union(feat_cols, intersect(c("feature","id","bit","un_id","un"), names(df_sh)))
    
    # Label-ish columns we could show on plots (priority order)
    label_priority <- c("smarts","smiles","pattern","formula","molecular_formula",
                        "label","name","description")
    label_cols <- names(df_sh)[colnames(df_sh) %in% label_priority |
                                 grepl(paste(label_priority, collapse="|"), colnames(df_sh))]
    
    if (length(feat_cols) == 0 || length(label_cols) == 0) next
    
    # Try all combos and pick the one with most usable rows
    best <- NULL
    best_n <- 0
    for (fc in feat_cols) {
      for (lc in label_cols) {
        tmp <- df_sh[, c(fc, lc)]
        colnames(tmp) <- c("Feature", "Label")
        tmp$Feature <- trimws(as.character(tmp$Feature))
        tmp$Label   <- trimws(as.character(tmp$Label))
        tmp <- tmp[nchar(tmp$Feature) > 0 & !is.na(tmp$Feature), ]
        # Keep rows that look like Un#### or exactMass or already SMARTS/SMILES text
        keep <- grepl("^un[_-]?\\d+$", tolower(tmp$Feature)) |
          tolower(tmp$Feature) %in% c("exactmass","monoisotopic.mass","monoisotopic_mass")
        tmp <- tmp[keep, , drop = FALSE]
        if (nrow(tmp) > best_n) { best <- tmp; best_n <- nrow(tmp) }
      }
    }
    if (!is.null(best) && nrow(best) > 0) candidates[[sh]] <- best
  }
  
  if (length(candidates) == 0) return(NULL)
  
  # Merge candidates; prefer first non-empty Label per Feature
  all_map <- do.call(rbind, candidates)
  all_map <- all_map[!duplicated(paste0(all_map$Feature, "_", all_map$Label)), ]
  # Collapse duplicates to first non-empty
  all_map <- aggregate(Label ~ Feature, data = all_map, function(v) v[which.max(nchar(v))])
  
  # Clean common variants
  all_map$Feature <- gsub("\\s+", "", all_map$Feature)
  all_map$Feature <- gsub("-", "_", all_map$Feature)
  all_map$Feature <- gsub("^UN", "Un", all_map$Feature, ignore.case = TRUE)
  
  # Nice default for exactMass if present
  all_map$Label[ tolower(all_map$Feature) %in% c("exactmass","monoisotopic.mass","monoisotopic_mass") ] <- "Exact mass"
  
  all_map
}

feature_label_map <- build_feature_label_map(in_xlsx)

if (!is.null(feature_label_map) && nrow(feature_label_map) > 0) {
  write.csv(feature_label_map, auto_map_csv, row.names = FALSE)
  message("??? Auto-built feature???label map: ", auto_map_csv)
} else {
  message("??? No mapping-like table found in the workbook. Using raw feature names.")
  feature_label_map <- data.frame(Feature = feature_names,
                                  Label   = feature_names, check.names = FALSE)
}

# --- 3) A lookup that always returns a label
get_feature_label <- function(feat) {
  hit <- which(feature_label_map$Feature == feat)
  if (length(hit)) feature_label_map$Label[hit[1]] else feat
}
