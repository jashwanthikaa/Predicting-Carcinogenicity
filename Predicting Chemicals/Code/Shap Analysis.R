# =============================================
# FINAL XGBOOST (caret) + SHAP EXPLANATIONS
# =============================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(readxl)
  library(caret)
  library(dplyr)
  library(xgboost)
  library(Matrix)
  library(pROC)
  library(ggplot2)
  library(SHAPforxgboost)
  library(gridExtra)
  library(data.table)
})

# ---------------------------------------------
#Paths
# ---------------------------------------------
data_path <- "C:/Users/91962/Predicting Chemicals/FinalData (1).xlsx"
out_dir   <- "C:/Users/91962/Predicting Chemicals/SHAP1/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------
# 1. Load and Prepare the Dataset
# ---------------------------------------------
df <- read_excel(data_path)

# Select fingerprint and exactMass columns
fp_cols <- grep("^Un", names(df), value = TRUE)
exactMass_col <- "exactMass"
X <- df %>% dplyr::select(all_of(c(fp_cols, exactMass_col)))

# Binary label from HIT.CALL
y <- ifelse(tolower(df$`HIT.CALL`) == "active", 1, 0)

# Remove rows with missing values
non_na_rows <- complete.cases(X)
X <- X[non_na_rows, ]
y <- y[non_na_rows]
df_clean <- df[non_na_rows, ]

# Remove near-zero variance features (only among fingerprints/mass)
nzv <- caret::nearZeroVar(X)
if (length(nzv) > 0) {
  X <- X[, -nzv]
}

# ---------------------------------------------
# 2. Create a Reference Train/Test Split (same seed)
# ---------------------------------------------
set.seed(123)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
y_train <- y[train_index]
X_test  <- X[-train_index, ]
y_test  <- y[-train_index]

# ---------------------------------------------
# 3. Combine Train + Test for Final Training
# ---------------------------------------------
X_full <- rbind(X_train, X_test)
y_full <- c(y_train, y_test)

# Factor with positive class first for caret::twoClassSummary
y_full_factor <- factor(ifelse(y_full == 1, "Active", "Inactive"),
                        levels = c("Active", "Inactive"))

X_full_matrix <- data.matrix(X_full)

# ---------------------------------------------
# 4. Class Imbalance Handling (case weights)
# ---------------------------------------------
pos_wt <- sum(y_full == 0) / sum(y_full == 1)
case_wts <- ifelse(y_full_factor == "Active", pos_wt, 1)

# ---------------------------------------------
# 5. Training Control and Tuning Grid
# ---------------------------------------------
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
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
# 6. Train Final XGBoost Model (caret::xgbTree)
# ---------------------------------------------
set.seed(123)
xgb_final_model <- train(
  x = X_full_matrix,
  y = y_full_factor,
  method = "xgbTree",
  trControl = train_control,
  tuneGrid  = xgb_grid,
  metric    = "ROC",
  verbose   = TRUE,
  weights   = case_wts
)

cat("\n Final Model Trained on All Data!\n")
print(xgb_final_model$bestTune)

# Save the model
saveRDS(xgb_final_model, file.path(out_dir, "xgb_final_model_trained_on_all_dataws.rds"))

# ---------------------------------------------
# 7. Predict and Evaluate on Full Data
# ---------------------------------------------
xgb_pred_prob <- predict(xgb_final_model, newdata = X_full_matrix, type = "prob")[, "Active"]

# ROC + optimal threshold
roc_obj <- pROC::roc(response = as.numeric(y_full_factor == "Active"), predictor = xgb_pred_prob)
best_coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))
optimal_threshold <- as.numeric(best_coords["threshold"])
cat("\n Optimal Threshold from ROC:", optimal_threshold, "\n")

# Final label predictions
xgb_pred <- ifelse(xgb_pred_prob > optimal_threshold, 1, 0)
y_full_factor_bin <- factor(y_full, levels = c(0, 1))
xgb_pred_factor_bin <- factor(xgb_pred, levels = c(0, 1))

# Confusion Matrix
xgb_cm <- caret::confusionMatrix(xgb_pred_factor_bin, y_full_factor_bin, positive = "1")
cat("\n Confusion Matrix on Full Data:\n")
print(xgb_cm)

# ROC plot
png(file.path(out_dir, "ROC_curve_final_model.png"), width = 900, height = 700)
plot(roc_obj, main = "ROC Curve - Final XGBoost Model (Train + Test Combined)")
dev.off()
cat("\nROC curve saved to 'ROC_curve_final_model.png'\n")

# ---------------------------------------------
# 8. Save Predictions with Metadata
# ---------------------------------------------
meta_cols <- intersect(c("Name", "CASRN", "SMILES", "exactMass", "DTXSID"), names(df_clean))
metadata_full <- df_clean[, meta_cols, drop = FALSE]

results <- data.frame(
  metadata_full,
  Actual = y_full,
  Predicted = xgb_pred,
  Probability = xgb_pred_prob,
  stringsAsFactors = FALSE
)

write.csv(results, file.path(out_dir, "xgb_predictions_trained_on_all_dataws.csv"), row.names = FALSE)
cat("\nPredictions saved to 'xgb_predictions_trained_on_all_data.csv'\n")

# ---------------------------------------------
# 9. SHAP EXPLANATIONS
# ---------------------------------------------
booster <- xgb_final_model$finalModel     # xgb.Booster used by caret
X_used  <- X_full_matrix

# 9.1 Compute SHAP values
shap_vals <- SHAPforxgboost::shap.values(xgb_model = booster, X_train = X_used)

global_imp <- data.frame(
  feature = names(shap_vals$mean_shap_score),
  mean_abs_shap = as.numeric(shap_vals$mean_shap_score),
  row.names = NULL
)

write.csv(global_imp, file.path(out_dir, "shap_global_importance.csv"), row.names = FALSE)

# 9.2 Prepare long data for plotting
shap_long <- SHAPforxgboost::shap.prep(xgb_model = booster, X_train = X_used)

# Summary (sina/violin-style)
p_sum <- SHAPforxgboost::shap.plot.summary(data_long = shap_long) +
  ggtitle("SHAP Summary (All Features)")
ggsave(file.path(out_dir, "SHAP_summary.png"), p_sum, width = 9, height = 7, dpi = 300)
ggsave(file.path(out_dir, "SHAP_summary.svg"), p_sum, width = 9, height = 7)

# Bar plot (Top features by mean |SHAP|)
p_bar <- SHAPforxgboost::shap.plot.summary(
  data_long = shap_long,
  scientific = FALSE,
  kind = "bar"
) + ggtitle("Top Features by mean(|SHAP|)")
ggsave(file.path(out_dir, "SHAP_bar_top.png"), p_bar, width = 9, height = 7, dpi = 300)
ggsave(file.path(out_dir, "SHAP_bar_top.svg"), p_bar, width = 9, height = 7)

# 9.3 Dependence plots (Top 10)
top_k <- 10
top_feats <- head(global_imp$feature[order(-global_imp$mean_abs_shap)], top_k)

make_dep_plot <- function(x_feat, y_feat = x_feat, color_feat = x_feat, smooth = TRUE) {
  SHAPforxgboost::shap.plot.dependence(
    data_long = shap_long,
    x = x_feat,
    y = y_feat,
    color_feature = color_feat,
    smooth = smooth
  ) + ggtitle(paste0("SHAP Dependence: ", y_feat, " vs ", x_feat))
}

dep_plots <- lapply(top_feats, function(f) make_dep_plot(f, f, f, smooth = TRUE))

png(file.path(out_dir, "SHAP_dependence_top10.png"), width = 1800, height = 1600)
gridExtra::grid.arrange(grobs = dep_plots, ncol = 2)
dev.off()

# Save individual dependence plots
dep_dir <- file.path(out_dir, "SHAP_dependence_plots")
if (!dir.exists(dep_dir)) dir.create(dep_dir, recursive = TRUE)
for (f in top_feats) {
  p <- make_dep_plot(f, f, f, smooth = TRUE)
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", f)
  ggsave(file.path(dep_dir, paste0("SHAP_dep_", safe_name, ".png")),
         p, width = 8, height = 6, dpi = 300)
}

cat("\nSHAP outputs saved:\n",
    "- shap_global_importance.csv\n",
    "- SHAP_summary.(png|svg)\n",
    "- SHAP_bar_top.(png|svg)\n",
    "- SHAP_dependence_top10.png + individual plots in 'SHAP_dependence_plots/'\n")
----------------------------------------------------------------------------------------
# ============================
# READABLE SHAP PLOTS (Top-20)
# ============================
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(gridExtra)
library(SHAPforxgboost)

# 1) Choose how many top features to show
TOP_N <- 20

# Importance table (mean |SHAP|)
imp_tbl <- global_imp %>%
  arrange(desc(mean_abs_shap)) %>%
  mutate(rank = row_number())

top_feats <- imp_tbl$feature[1:min(TOP_N, nrow(imp_tbl))]

# 2) Nice, readable Top-N bar plot
imp_top <- imp_tbl %>% slice_head(n = TOP_N)

p_bar20 <- ggplot(imp_top, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Avg(|SHAP|)", title = "Top Features by mean(|SHAP|)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave(file.path(out_dir, "SHAP_bar_top20.png"), p_bar20, width = 8.5, height = 10, dpi = 350)
ggsave(file.path(out_dir, "SHAP_bar_top20.svg"), p_bar20, width = 8.5, height = 10)

# 3) Identify binary vs continuous features from the actual training matrix
X_df <- as.data.frame(X_used)
is_binary <- vapply(X_df, function(col) {
  u <- unique(col)
  length(u) <= 2 && all(u %in% c(0, 1))
}, logical(1))

# For binary features, compute prevalence (share of 1s)
bin_prevalence <- sapply(names(X_df), function(nm) {
  if (!is_binary[[nm]]) return(NA_real_)
  mean(X_df[[nm]] == 1)
})

# 4) Save a CSV with Top-N + prevalence
imp_top %>%
  mutate(
    type = ifelse(feature %in% names(is_binary)[is_binary], "binary", "continuous"),
    prevalence_1 = ifelse(type == "binary", round(bin_prevalence[feature], 4), NA_real_)
  ) %>%
  write.csv(file.path(out_dir, "shap_top20_with_prevalence.csv"), row.names = FALSE)

# 5) Dependence plots: split by feature type
top_bin  <- top_feats[top_feats %in% names(is_binary)[is_binary]]
top_cont <- setdiff(top_feats, top_bin)

# Helper to make a readable file name
safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# Create an output folder
dep_dir <- file.path(out_dir, "SHAP_dependence_top20(1)")
if (!dir.exists(dep_dir)) dir.create(dep_dir, recursive = TRUE)

# 5a) CONTINUOUS features: keep smoothing
for (f in top_cont) {
  p <- SHAPforxgboost::shap.plot.dependence(
    data_long = shap_long,
    x = f, y = f, color_feature = f,
    smooth = TRUE
  ) +
    ggtitle(paste0("SHAP Dependence (continuous): ", f))
  ggsave(file.path(dep_dir, paste0("DEP_CONT_", safe_name(f), ".png")),
         p, width = 10, height = 8, dpi = 350)
}

# 5b) BINARY features: no smoothing, jitter to separate points; show prevalence in title
for (f in top_bin) {
  prev <- bin_prevalence[[f]]
  p <- SHAPforxgboost::shap.plot.dependence(
    data_long = shap_long,
    x = f, y = f, color_feature = f,
    smooth = FALSE,              # <<< no smoothing for 0/1 bits
    jitter_width = 0.12,
    jitter_height = 0.0,
    alpha = 0.8,
    size0 = 0.9
  ) +
    ggtitle(paste0("SHAP Dependence (binary): ", f,
                   "   |   prevalence(1) = ", sprintf("%.1f%%", 100*prev)))
  ggsave(file.path(dep_dir, paste0("DEP_BIN_", safe_name(f), ".png")),
         p, width = 10, height = 8, dpi = 350)
}

# 6) Also create two contact sheets for quick viewing
# (one for continuous, one for binary) - only if there are multiple plots
if (length(top_cont) > 0) {
  cont_plots <- lapply(top_cont, function(f) {
    SHAPforxgboost::shap.plot.dependence(
      data_long = shap_long, x = f, y = f, color_feature = f, smooth = TRUE
    ) + ggtitle(paste0("Continuous: ", f))
  })
  png(file.path(out_dir, "SHAP_dependence_CONT_contact.png"), width = 2200, height = 1800)
  gridExtra::grid.arrange(grobs = cont_plots, ncol = 2)
  dev.off()
}

if (length(top_bin) > 0) {
  bin_plots <- lapply(top_bin, function(f) {
    prev <- bin_prevalence[[f]]
    SHAPforxgboost::shap.plot.dependence(
      data_long = shap_long, x = f, y = f, color_feature = f,
      smooth = FALSE, jitter_width = 0.12, alpha = 0.8, size0 = 0.9
    ) + ggtitle(paste0("Binary: ", f, "  |  p(1)=", sprintf("%.1f%%", 100*prev)))
  })
  png(file.path(out_dir, "SHAP_dependence_BIN_contact.png"), width = 2200, height = 1800)
  gridExtra::grid.arrange(grobs = bin_plots, ncol = 2)
  dev.off()
}

cat("\nClean SHAP visuals saved:\n",
    " - SHAP_bar_top20.(png|svg)\n",
    " - shap_top20_with_prevalence.csv\n",
    " - SHAP_dependence_top20/DEP_CONT_*.png & DEP_BIN_*.png\n",
    " - SHAP_dependence_CONT_contact.png (if any continuous)\n",
    " - SHAP_dependence_BIN_contact.png (if any binary)\n")
----------------------------------------------------------------------------------------------
#Panel of SHAP dependence plots (clean, high-res)
#library(gridExtra)
  
# Make dependence plots for top 9 features
TOP_N <- 9
imp_ord <- global_imp$feature[order(-global_imp$mean_abs_shap)]
top_feats <- head(imp_ord, TOP_N)

plots <- lapply(top_feats, function(f) {
  SHAPforxgboost::shap.plot.dependence(
    data_long = shap_long,
    x = f, y = f, color_feature = f,
    smooth = TRUE
  ) + labs(title = f)
})

# Save as large PNG for clarity
png(file.path(out_dir, "SHAP_dep_panel_top9.png"), width = 3600, height = 3000, res = 400)
grid.arrange(grobs = plots, ncol = 3)
dev.off()

# Also save each individually (super clear)
for (f in top_feats) {
  p <- SHAPforxgboost::shap.plot.dependence(
    data_long = shap_long,
    x = f, y = f, color_feature = f,
    smooth = TRUE
  ) + labs(title = f)
  ggsave(file.path(out_dir, paste0("SHAP_dep_", gsub("[^A-Za-z0-9]", "_", f), ".png")),
         p, width = 8, height = 6, dpi = 400)
}

----------------------------------------------------------------------------------
#"Textbook" SHAP visuals (beeswarm + waterfall for one chemical)

# Install once: install.packages("shapviz")
install.packages("shapviz")
library(shapviz)

# booster: xgb_final_model$finalModel
# X_used: numeric matrix you trained on (same features, same order)

sv <- shapviz(
  xgb_model = xgb_final_model$finalModel,
  X = X_used,
  X_pred = X_used,  # explain all rows
  X_names = colnames(X_used)
)

# (a) Global importance (bar) - mean(|SHAP|)
p_imp <- sv_importance(sv, kind = "bar")
ggsave(file.path(out_dir, "SV_importance_bar.png"), p_imp, width = 8, height = 6, dpi = 300)

# (b) Beeswarm (summary) - like the classic SHAP summary with colors = feature value
p_bee <- sv_importance(sv, kind = "beeswarm")
ggsave(file.path(out_dir, "SV_beeswarm.png"), p_bee, width = 9, height = 7, dpi = 300)

# (c) Dependence for a specific feature
p_dep <- sv_dependence(sv, v = "exactMass", color_var = "exactMass")
ggsave(file.path(out_dir, "SV_dependence_exactMass.png"), p_dep, width = 7, height = 6, dpi = 300)

# (d) Waterfall for a single observation (choose an index you care about)
i <- 1  # e.g., first row
p_wf <- sv_waterfall(sv, row_index = i, max_display = 15)
ggsave(file.path(out_dir, sprintf("SV_waterfall_obs_%03d.png", i)), p_wf, width = 8, height = 6, dpi = 300)

--------------------------------------------------------------------------------------------
# ===============================
# ULTRA-CRISP SHAP PLOTS (no magick)
# ===============================
suppressPackageStartupMessages({
  library(shapviz)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(svglite)
})

# --- Inputs you already have ---
booster <- xgb_final_model$finalModel
X_used  <- X_full_matrix

# --- Output folder ---
out_dir <- "C:/Users/91962/Predicting Chemicals/SHAP Analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: save PNG(600dpi) + SVG + PDF ---
save_all <- function(p, name, w=10, h=8, dpi=600) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width=w, height=h, dpi=dpi)
  ggsave(file.path(out_dir, paste0(name, ".svg")), p, width=w, height=h, device = svglite)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width=w, height=h, device = cairo_pdf)
}

# --- Build shapviz object ---
sv <- shapviz(
  object  = booster,
  X       = X_used,
  X_pred  = X_used,
  X_names = colnames(X_used)
)

# ---------- 1) Global importance ----------
p_bar <- sv_importance(sv, kind = "bar") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold"))
save_all(p_bar, "SV_importance_bar", w=11, h=8)

# ---------- 2) Beeswarm ----------
p_bee <- sv_importance(sv, kind = "beeswarm") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold"))
save_all(p_bee, "SV_beeswarm", w=12, h=9)

# ---------- 3) Dependence: Top-15 ----------
S <- as.matrix(sv$S)              
abs_means <- colMeans(abs(S))
ord <- names(sort(abs_means, decreasing = TRUE))
top_feats <- head(ord, 15)

dep_plots <- list()
for (f in top_feats) {
  p <- sv_dependence(sv, v = f, color_var = f) +
    labs(title = paste0("SHAP Dependence: ", f),
         x = f, y = "SHAP value (impact)") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(face = "bold"))
  fname <- paste0("SV_dependence_", gsub("[^A-Za-z0-9_]+","_", f))
  save_all(p, fname, w=11, h=8)   # individual high-res + vector
  dep_plots[[f]] <- p
}

# ---------- 4) Large multi-plot panel ----------
# Arrange into a big PNG/SVG/PDF using gridExtra
panel <- gridExtra::grid.arrange(grobs = dep_plots, ncol = 3)
ggsave(file.path(out_dir, "SV_dependence_top15_panel.png"), panel, width = 20, height = 16, dpi = 300)
ggsave(file.path(out_dir, "SV_dependence_top15_panel.svg"), panel, width = 20, height = 16)
ggsave(file.path(out_dir, "SV_dependence_top15_panel.pdf"), panel, width = 20, height = 16)

# ---------- 5) Example: exactMass ----------
if ("exactMass" %in% colnames(X_used)) {
  p_mass <- sv_dependence(sv, v = "exactMass", color_var = "exactMass") +
    labs(title = "SHAP Dependence: exactMass",
         x = "exactMass", y = "SHAP value (impact)") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(face = "bold"))
  save_all(p_mass, "SV_dependence_exactMass", w=11, h=8)
}

message("All plots saved in: ", normalizePath(out_dir))
