# =============================================
# FINAL XGBOOST TRAINING ON TRAIN + TEST SET
# =============================================

# ---- 0) Libraries ----
suppressPackageStartupMessages({
  library(readxl)
  library(caret)
  library(dplyr)
  library(xgboost)
  library(Matrix)
  library(pROC)
})

# ---- 1) Paths ----
data_path   <- "C:/Users/91962/Predicting Chemicals/FinalData (1).xlsx"
out_dir     <- "C:/Users/91962/Predicting Chemicals"
model_path  <- file.path(out_dir, "xgb_final_model_trained_on_all_data11.rds")
roc_png     <- file.path(out_dir, "ROC_curve_final_model.png")
cm_txt      <- file.path(out_dir, "confusion_matrix_final_model1.txt")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 2) Load & Prepare Data ----
df <- read_excel(data_path)

# Fingerprints start with "Un" + add exactMass
fp_cols <- grep("^Un", names(df), value = TRUE)
exactMass_col <- "exactMass"

if (!exactMass_col %in% names(df)) {
  stop(sprintf("Column '%s' not found in the data.", exactMass_col))
}

X <- df %>% dplyr::select(all_of(c(fp_cols, exactMass_col)))
y <- ifelse(tolower(df$`HIT.CALL`) == "active", 1, 0)

# Remove rows with missing values
non_na_rows <- complete.cases(X)
X <- X[non_na_rows, , drop = FALSE]
y <- y[non_na_rows]
df_clean <- df[non_na_rows, , drop = FALSE]

# Sanity checks
if (length(unique(y)) < 2) {
  stop("After cleaning, y has only one class. Cannot train or compute ROC.")
}

# Remove near-zero variance features
nzv_idx <- nearZeroVar(X)
if (length(nzv_idx) > 0) {
  X <- X[, -nzv_idx, drop = FALSE]
}

# ---- 3) Create a reproducible Train/Test split (for reference only) ----
set.seed(123)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, , drop = FALSE]
y_train <- y[train_index]
X_test  <- X[-train_index, , drop = FALSE]
y_test  <- y[-train_index]

# ---- 4) Combine Train + Test for final training ----
X_full <- rbind(X_train, X_test)
y_full <- c(y_train, y_test)

# Explicit factor: make "Active" the positive class (first level for twoClassSummary)
y_full_factor <- factor(ifelse(y_full == 1, "Active", "Inactive"),
                        levels = c("Active", "Inactive"))

# Matrix for xgboost / caret
X_full_matrix <- data.matrix(X_full)

# ---- 5) Handle class imbalance ----
scale_weight <- sum(y_full == 0) / sum(y_full == 1)

# ---- 6) Training control & grid ----
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,   # computes ROC/Sens/Spec for first level = "Active"
  allowParallel = TRUE
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

# ---- 7) Train final model on all data ----
set.seed(123)
xgb_final_model <- train(
  x = X_full_matrix,
  y = y_full_factor,
  trControl = train_control,
  tuneGrid = xgb_grid,
  method = "xgbTree",
  metric = "ROC",
  verbose = TRUE,
  scale_pos_weight = scale_weight    # passes through to xgboost
)

cat("\nFinal model trained on all data.\n")
print(xgb_final_model$bestTune)

# Save model
saveRDS(xgb_final_model, model_path)
cat(sprintf("Model saved to: %s\n", model_path))

# ---- 8) Predict probabilities on full data ----
pred_prob_df <- predict(xgb_final_model, newdata = X_full_matrix, type = "prob")
if (!"Active" %in% colnames(pred_prob_df)) {
  stop("Predicted probability column 'Active' not found. Check factor levels for y.")
}
xgb_pred_prob <- pred_prob_df[, "Active"]

# ---- 9) Safe ROC / Optimal threshold ----
# Basic checks before ROC
if (anyNA(y_full) || anyNA(xgb_pred_prob)) {
  stop("ROC cannot be computed: NA values present in y_full or predictions.")
}
if (length(unique(y_full)) < 2) {
  stop("ROC cannot be computed: only one class present in y_full.")
}

roc_obj <- pROC::roc(
  response  = y_full,           # numeric 0/1
  predictor = xgb_pred_prob,
  levels    = c(0, 1),          # 0 = control, 1 = case
  direction = "<",              # higher score => more likely case=1
  quiet     = TRUE
)

auc_val <- as.numeric(pROC::auc(roc_obj))
if (is.na(auc_val)) {
  stop("AUC is NA; check predictions and labels.")
}
cat(sprintf("AUC: %.4f\n", auc_val))

best_coords <- pROC::coords(
  roc_obj, "best",
  ret = c("threshold", "specificity", "sensitivity"),
  best.method = "youden"        # Youden's J
)
optimal_threshold <- as.numeric(best_coords["threshold"])
cat(sprintf("Optimal threshold (Youden): %.4f\n", optimal_threshold))
cat(sprintf("At threshold: Sens=%.3f, Spec=%.3f\n",
            as.numeric(best_coords["sensitivity"]),
            as.numeric(best_coords["specificity"])))

# ---- 10) Final class predictions & Confusion Matrix ----
xgb_pred <- ifelse(xgb_pred_prob > optimal_threshold, 1, 0)

y_full_factor_bin   <- factor(y_full, levels = c(0, 1))
xgb_pred_factor_bin <- factor(xgb_pred, levels = c(0, 1))

xgb_cm <- caret::confusionMatrix(
  data = xgb_pred_factor_bin,
  reference = y_full_factor_bin,
  positive = "1"
)

cat("\nConfusion Matrix on Full Data:\n")
print(xgb_cm)

# Save CM to text
sink(cm_txt)
cat("Confusion Matrix (Final Model on Full Data)\n\n")
print(xgb_cm)
sink()
cat(sprintf("Confusion matrix saved to: %s\n", cm_txt))

# ---- 11) Plot ROC Curve ----
png(filename = roc_png, width = 900, height = 700)
plot(roc_obj,
     main = "ROC Curve - Final XGBoost Model (Train + Test Combined)")
dev.off()
cat(sprintf("ROC curve saved to: %s\n", roc_png))

# ---- 12) Quick summary line ----
cat(sprintf("\nSUMMARY:\nAUC=%.4f | Threshold=%.4f | Sens=%.3f | Spec=%.3f\n",
            auc_val,
            optimal_threshold,
            as.numeric(best_coords["sensitivity"]),
            as.numeric(best_coords["specificity"])))
------------------------------------------------------------------
  # ============================
# Figures tied to your model (robust)
# ============================
out_dir <- "C:/Users/91962/Predicting Chemicals/Plots(Modelsb)"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

need <- c("ggplot2","pROC","xgboost","reshape2","gridExtra")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(ggplot2)
library(pROC)
library(xgboost)
library(reshape2)
library(gridExtra)

# ---- This script expects these in memory:
# X_train_matrix, X_test_matrix, y_test, xgb_model (caret::train object)

# --- 0) Sanity checks on inputs ---
stopifnot(exists("xgb_model"))
stopifnot(exists("X_test_matrix"))
stopifnot(exists("y_test"))

if (!is.matrix(X_test_matrix)) {
  stop("X_test_matrix must be a numeric matrix.")
}

# --- 1) Recompute predicted probabilities on the *same* rows used for test ---
pred_df <- predict(xgb_model, newdata = X_test_matrix, type = "prob")
# Make sure we grab the right column name for event class
prob_col <- if ("Active" %in% colnames(pred_df)) "Active" else colnames(pred_df)[1]
xgb_pred_prob <- as.numeric(pred_df[, prob_col])

# --- 2) Coerce y_test to a clean 0/1 numeric vector aligned to X_test_matrix ---
y_vec <- y_test
# if it's a data frame / matrix, take the first column
if (is.data.frame(y_vec) || is.matrix(y_vec)) y_vec <- y_vec[, 1]

if (is.factor(y_vec)) {
  # common cases: levels = c("Inactive","Active") or c("0","1") or c("1","0")
  lv <- levels(y_vec)
  if (all(c("Inactive","Active") %in% lv)) {
    y_vec <- as.numeric(y_vec == "Active")
  } else {
    # try numeric conversion through character
    y_vec <- suppressWarnings(as.numeric(as.character(y_vec)))
  }
} else if (is.character(y_vec)) {
  y_low <- tolower(y_vec)
  if (all(y_low %in% c("active","inactive"))) {
    y_vec <- as.numeric(y_low == "active")
  } else {
    y_vec <- suppressWarnings(as.numeric(y_vec))
  }
}

# --- 3) Remove any NAs in either vector and ensure equal length ---
ok_idx <- which(!is.na(y_vec) & !is.na(xgb_pred_prob))
y_vec <- y_vec[ok_idx]
xgb_pred_prob <- xgb_pred_prob[ok_idx]

if (length(y_vec) != length(xgb_pred_prob)) {
  stop(sprintf("After alignment, lengths still differ: y=%d, prob=%d",
               length(y_vec), length(xgb_pred_prob)))
}
if (length(unique(y_vec)) < 2) {
  stop("Test labels contain only one class after cleaning; ROC cannot be computed.")
}

# --- 4) Feature names (for importance heatmap) ---
feature_names <- colnames(X_test_matrix)
if (is.null(feature_names)) {
  feature_names <- paste0("f", seq_len(ncol(X_test_matrix)))
  colnames(X_test_matrix) <- feature_names
}

# ============================
# Panel A - probability distribution (test set)
# ============================
roc_best <- pROC::roc(response = y_vec, predictor = xgb_pred_prob,
                      levels = c(0,1), direction = "<", quiet = TRUE)
thr_best <- as.numeric(pROC::coords(roc_best, "best", ret = "threshold"))

df_prob <- data.frame(
  Prob = xgb_pred_prob,
  Class = factor(ifelse(y_vec == 1, "Active", "Inactive"),
                 levels = c("Inactive","Active"))
)

pA <- ggplot(df_prob, aes(x = Prob, fill = Class)) +
  geom_histogram(position = "identity", alpha = 0.55, bins = 40, color = "white") +
  geom_vline(xintercept = thr_best, linetype = "dashed") +
  labs(title = "A", subtitle = "Predicted probability (test set)",
       x = "Predicted probability of Active", y = "Count",
       caption = paste0("ROC-optimal threshold = ", sprintf("%.3f", thr_best),
                        " | AUC = ", sprintf("%.3f", pROC::auc(roc_best)))) +
  scale_fill_manual(values = c("Inactive" = "#6baed6", "Active" = "#fb6a4a")) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "PanelA_ProbDist_Test.png"), pA, width = 7, height = 5, dpi = 300)

# ============================
# Panel B - correlation of TOP-K most important features
# ============================
TOP_K <- 12  # tweak as desired

# Pull importance from the final booster inside caret model
imp <- tryCatch(
  xgboost::xgb.importance(model = xgb_model$finalModel, feature_names = feature_names),
  error = function(e) {
    stop("Could not extract xgb.importance(): ", e$message)
  }
)

if (nrow(imp) == 0) stop("xgb.importance returned zero rows.")
imp_top <- head(imp, TOP_K)
top_feats <- imp_top$Feature

# subset test matrix to top features (only those that exist)
top_feats <- intersect(top_feats, colnames(X_test_matrix))
if (length(top_feats) < 2) stop("Fewer than 2 top features available in X_test_matrix.")

X_top <- as.matrix(X_test_matrix[, top_feats, drop = FALSE])

# Correlation and ordering
C <- suppressWarnings(cor(X_top, use = "pairwise.complete.obs"))
C[is.na(C)] <- 0; diag(C) <- 1
D <- as.dist(pmax(0, 1 - C)); D[!is.finite(D)] <- 0
ord <- tryCatch(hclust(D, method = "average")$order, error = function(e) seq_len(ncol(C)))
Ck <- C[ord, ord]

# Upper triangle only
Ck[lower.tri(Ck)] <- NA
hm <- reshape2::melt(Ck, varnames = c("Var1","Var2"), value.name = "corr", na.rm = TRUE)

pB <- ggplot(hm, aes(Var2, Var1, fill = corr)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0,1), low = "white", high = "#08519C", name = "Correlation") +
  coord_fixed() +
  labs(title = "B", subtitle = paste0("Correlation among top-", length(top_feats), " important features"),
       x = NULL, y = NULL,
       caption = "Importance ranked by XGBoost gain on the training split") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank())

ggsave(file.path(out_dir, "PanelB_TopFeatures_Correlation.png"), pB, width = 7.8, height = 6, dpi = 300)

# ============================
# Combine A + B into one figure
# ============================
combo <- gridExtra::arrangeGrob(pA, pB, ncol = 2, widths = c(1, 1.15))
ggsave(file.path(out_dir, "Figure_ProbDist_and_TopFeatureCorrelation.png"),
       combo, width = 14.5, height = 6, dpi = 300)

cat("Saved:\n",
    file.path(out_dir, "PanelA_ProbDist_Test.png"), "\n",
    file.path(out_dir, "PanelB_TopFeatures_Correlation.png"), "\n",
    file.path(out_dir, "Figure_ProbDist_and_TopFeatureCorrelation.png"), "\n", sep = "")
