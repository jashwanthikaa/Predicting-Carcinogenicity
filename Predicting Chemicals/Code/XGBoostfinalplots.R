# ============================
#  Setup
# ============================
out_dir <- "C:/Users/91962/Predicting Chemicals/Plots(Models)"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

library(xgboost)
library(pROC)
library(dplyr)
library(caret)
library(ggplot2)
library(reshape2)
library(fmsb)
library(scales)

# ============================
# 1. Train XGBoost models for 3 nrounds
# ============================
params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.01,
  gamma = 5,
  max_depth = 6,
  subsample = 0.6,
  colsample_bytree = 0.6,
  scale_pos_weight = scale_weight
)

dtrain <- xgb.DMatrix(data = X_train_matrix, label = y_train)
dtest  <- xgb.DMatrix(data = X_test_matrix,  label = y_test)

rounds <- c(100, 200, 300)
models <- list()
preds  <- list()
rocs   <- list()
aucs   <- numeric(length(rounds))
names(aucs) <- paste0("n", rounds)

set.seed(123)
for (i in seq_along(rounds)) {
  nr <- rounds[i]
  models[[i]] <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = nr,
    watchlist = list(train = dtrain),
    verbose = 0
  )
  preds[[i]] <- predict(models[[i]], dtest)
  rocs[[i]] <- roc(response = y_test, predictor = preds[[i]], levels = c(0,1), direction = "<")
  aucs[i] <- as.numeric(auc(rocs[[i]]))
}

# ============================
# 2. Metrics Table
# ============================
summary_list <- list()
for (i in seq_along(rounds)) {
  best <- coords(rocs[[i]], "best",
                 ret = c("threshold","specificity","sensitivity","precision","recall"))
  thr <- as.numeric(best["threshold"])
  pred_bin <- ifelse(preds[[i]] > thr, 1, 0)
  cm <- confusionMatrix(
    factor(pred_bin, levels = c(0,1)),
    factor(y_test,  levels = c(0,1)),
    positive = "1"
  )
  summary_list[[i]] <- data.frame(
    nrounds = rounds[i],
    AUC = aucs[i],
    threshold = thr,
    Accuracy = cm$overall["Accuracy"] * 100,
    Kappa = cm$overall["Kappa"],
    Sensitivity = cm$byClass["Sensitivity"] * 100,
    Specificity = cm$byClass["Specificity"] * 100,
    Precision = cm$byClass["Pos Pred Value"] * 100,
    NPV = cm$byClass["Neg Pred Value"] * 100,
    Balanced_Accuracy = cm$byClass["Balanced Accuracy"] * 100
  )
}
summary_df <- bind_rows(summary_list)
write.csv(summary_df, file.path(out_dir, "ROC_comparison_summary_100_200_300.csv"), row.names = FALSE)

# ----------------------------
# 3. ROC Curves (combined) - FIXED
# ----------------------------
roc_file <- file.path(out_dir, "ROC_Models_100_200_300.png")
png(filename = roc_file, width = 1400, height = 900, res = 150)

# First ROC curve
plot.roc(rocs[[1]], col = "blue", legacy.axes = TRUE,
         main = "ROC Curves - XGBoost (nrounds = 100, 200, 300)",
         lwd = 2, print.auc = FALSE)

# Add the others
plot.roc(rocs[[2]], col = "red", add = TRUE, lwd = 2, print.auc = FALSE)
plot.roc(rocs[[3]], col = "darkgreen", add = TRUE, lwd = 2, print.auc = FALSE)

legend("bottomright",
       legend = c(
         paste0("n=100 (AUC=", sprintf("%.3f", aucs[1]), ")"),
         paste0("n=200 (AUC=", sprintf("%.3f", aucs[2]), ")"),
         paste0("n=300 (AUC=", sprintf("%.3f", aucs[3]), ")")
       ),
       col = c("blue","red","darkgreen"), lwd = 2, bty = "n")

dev.off()

# ============================
# 4. Grouped Bar Chart
# ============================
metrics_df <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity", "Precision", "NPV", "Balanced Accuracy", "ROC AUC"),
  n100 = c(79.02, 66.67, 81.51, 42.11, 92.38, 74.09, 85.0),
  n200 = c(64.34, 87.50, 59.66, 30.43, 95.95, 73.58, 86.0),
  n300 = c(66.43, 83.33, 63.03, 31.25, 94.94, 73.18, 87.0)
)

metrics_melt <- melt(metrics_df, id.vars = "Metric", variable.name = "Model", value.name = "Value")

p_bar <- ggplot(metrics_melt, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f", Value)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = "Performance Metrics by Model (nrounds)", y = "Value (%)", x = "") +
  scale_fill_manual(values = c("n100" = "blue", "n200" = "red", "n300" = "darkgreen")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "Grouped_Bar_Metrics.png"), p_bar, width = 10, height = 6)

# ============================
# 5. Radar Chart
# ============================
radar_df <- metrics_df
rownames(radar_df) <- radar_df$Metric
radar_df <- radar_df[, -1]
radar_df <- rbind(max = rep(100, ncol(radar_df)),
                  min = rep(0, ncol(radar_df)),
                  radar_df)
png(file.path(out_dir, "Radar_Metrics.png"), width = 1200, height = 800, res = 150)
radarchart(radar_df,
           axistype = 1,
           pcol = c("blue", "red", "darkgreen"),
           pfcol = alpha(c("blue", "red", "darkgreen"), 0.3),
           plwd = 2, plty = 1,
           cglcol = "grey", cglty = 1,
           axislabcol = "grey", caxislabels = seq(0,100,20),
           vlcex = 0.8,
           title = "Radar Chart of Model Performance")
legend("topright", legend = c("n100", "n200", "n300"),
       bty = "n", pch = 20, col = c("blue", "red", "darkgreen"), text.col = "black", cex = 1.2, pt.cex = 1.5)
dev.off()

# ============================
# 6. Metric Trend Lines
# ============================
metrics_long <- metrics_melt
metrics_long$Model <- as.numeric(gsub("n", "", metrics_long$Model))
p_trend <- ggplot(metrics_long, aes(x = Model, y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "Metric Trends Across nrounds", x = "nrounds", y = "Value (%)") +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(100, 200, 300))
ggsave(file.path(out_dir, "Metric_Trends.png"), p_trend, width = 10, height = 6)

# ============================
# 7. Confusion Matrix Heatmaps
# ============================
cm_data <- list(
  n100 = data.frame(Actual = rep(c("Positive","Negative"), each = 2),
                    Predicted = rep(c("Positive","Negative"), times = 2),
                    Count = c(16, 8, 22, 97)),
  n200 = data.frame(Actual = rep(c("Positive","Negative"), each = 2),
                    Predicted = rep(c("Positive","Negative"), times = 2),
                    Count = c(21, 3, 48, 71)),
  n300 = data.frame(Actual = rep(c("Positive","Negative"), each = 2),
                    Predicted = rep(c("Positive","Negative"), times = 2),
                    Count = c(20, 4, 44, 75))
)

plot_cm_heatmap <- function(df, model_label, fill_color, save_path) {
  p <- ggplot(df, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count), size = 5, color = "black") +
    scale_fill_gradient(low = "white", high = fill_color) +
    labs(title = paste("Confusion Matrix -", model_label),
         x = "Predicted", y = "Actual") +
    theme_minimal(base_size = 14)
  ggsave(save_path, p, width = 5, height = 4)
}

plot_cm_heatmap(cm_data$n100, "nrounds = 100", "blue", file.path(out_dir, "ConfusionMatrix_n100.png"))
plot_cm_heatmap(cm_data$n200, "nrounds = 200", "red", file.path(out_dir, "ConfusionMatrix_n200.png"))
plot_cm_heatmap(cm_data$n300, "nrounds = 300", "darkgreen", file.path(out_dir, "ConfusionMatrix_n300.png"))

cat("??? All plots & metrics saved in:", out_dir, "\n")



--------------------------------------------------------------------------------------------
#FINAL Model plots 
  # ============================================
# Best vs Final - Publication-style Figure
# ============================================

# ---- 0) Setup & packages ----
out_dir <- "C:/Users/91962/Predicting Chemicals/Plots(Models6)"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

req <- c("pROC","ggplot2","grid","gridExtra","png","scales","caret")
to_install <- setdiff(req, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, dependencies = TRUE)

# Try to use fmsb for radar; if not installed, we'll fall back to ggplot radar
has_fmsb <- "fmsb" %in% rownames(installed.packages())
if(!has_fmsb) install.packages("fmsb", dependencies = TRUE)
has_fmsb <- "fmsb" %in% rownames(installed.packages())

library(pROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(png)
library(scales)
library(caret)
if (has_fmsb) library(fmsb)

# ---- helpers ----
binify <- function(y) {
  if (is.factor(y)) return(as.numeric(as.character(y)))
  if (is.logical(y)) return(as.numeric(y))
  as.numeric(y)
}

plot_cm_gg <- function(cm, title_txt, hi_col = "steelblue") {
  df <- as.data.frame(cm$table)
  names(df) <- c("Predicted","Actual","Count")
  ggplot(df, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count), size = 6) +
    scale_fill_gradient(low = "white", high = hi_col) +
    labs(title = title_txt, x = "Predicted", y = "Actual") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
}

to_grob <- function(path_png) rasterGrob(png::readPNG(path_png), interpolate = TRUE)

# ============================================
# 1) ROC objects (namespaced pROC::roc) + thresholds
# ============================================
roc_best  <- pROC::roc(response = binify(y_test), predictor = xgb_pred_prob, quiet = TRUE)
roc_final <- pROC::roc(response = binify(y_full), predictor = xgb_pred_prob_full, quiet = TRUE)
auc_best  <- as.numeric(pROC::auc(roc_best))
auc_final <- as.numeric(pROC::auc(roc_final))

thr_best  <- as.numeric(pROC::coords(roc_best,  "best", ret = "threshold"))
thr_final <- as.numeric(pROC::coords(roc_final, "best", ret = "threshold"))

# ============================================
# 2) Confusion matrices (rebuild if missing)
# ============================================
if (!exists("xgb_cm")) {
  pred_best_bin <- ifelse(xgb_pred_prob > thr_best, 1, 0)
  xgb_cm <- caret::confusionMatrix(
    factor(pred_best_bin, levels = c(0,1)),
    factor(binify(y_test), levels = c(0,1)),
    positive = "1"
  )
}

if (!exists("xgb_cm_full")) {
  pred_final_bin <- ifelse(xgb_pred_prob_full > thr_final, 1, 0)
  xgb_cm_full <- caret::confusionMatrix(
    factor(pred_final_bin, levels = c(0,1)),
    factor(binify(y_full), levels = c(0,1)),
    positive = "1"
  )
}

# ============================================
# 3) Metrics (compute from CMs so they stay consistent)
# ============================================
met_from_cm <- function(cm) {
  data.frame(
    Accuracy        = as.numeric(cm$overall["Accuracy"]) * 100,
    Sensitivity     = as.numeric(cm$byClass["Sensitivity"]) * 100,
    Specificity     = as.numeric(cm$byClass["Specificity"]) * 100,
    Precision       = as.numeric(cm$byClass["Pos Pred Value"]) * 100,
    NPV             = as.numeric(cm$byClass["Neg Pred Value"]) * 100,
    BalancedAcc     = as.numeric(cm$byClass["Balanced Accuracy"]) * 100,
    check.names = FALSE
  )
}
m_best  <- met_from_cm(xgb_cm)
m_final <- met_from_cm(xgb_cm_full)

metrics_df <- data.frame(
  Metric = c("Accuracy","Sensitivity","Specificity","Precision","NPV","BalancedAcc"),
  Best   = unlist(m_best[1, ]),
  Final  = unlist(m_final[1, ]),
  row.names = NULL,
  check.names = FALSE
)

# ============================================
# 4) Draw individual panels to PNG
# ============================================

## (A) Combined ROC
file_roc <- file.path(out_dir, "TMP_panel_ROC.png")
png(file_roc, width = 1800, height = 1400, res = 200)
pROC::plot.roc(roc_best,  col = "blue",     lwd = 3, legacy.axes = TRUE,
               main = "(A) ROC - Best (test) vs Final (full)")
pROC::plot.roc(roc_final, col = "firebrick", lwd = 3, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "grey50", lwd = 2)
legend("bottomright",
       legend = c(
         paste0("Best  (AUC=", sprintf("%.3f", auc_best),  ", thr=", sprintf("%.3f", thr_best),  ")"),
         paste0("Final (AUC=", sprintf("%.3f", auc_final), ", thr=", sprintf("%.3f", thr_final), ")")
       ),
       col = c("blue","firebrick"), lwd = 3, bty = "n", cex = 1.1)
dev.off()

## (B) Radar chart - try fmsb, else ggplot fallback
file_radar <- file.path(out_dir, "TMP_panel_Radar.png")
if (has_fmsb) {
  rad <- metrics_df
  rownames(rad) <- rad$Metric
  rad <- rad[, -1]
  rad <- rbind(max = rep(100, ncol(rad)),
               min = rep(0,   ncol(rad)),
               rad)
  png(file_radar, width = 1800, height = 1400, res = 200)
  fmsb::radarchart(rad,
                   axistype = 1,
                   pcol  = c("blue","firebrick"),
                   pfcol = alpha(c("blue","firebrick"), 0.25),
                   plwd  = 3, plty = 1,
                   cglcol = "grey60", cglty = 1, cglwd = 1,
                   axislabcol = "grey30", caxislabels = seq(0,100,20),
                   vlcex = 1.1,
                   title = "(B) Metrics Radar - Best vs Final")
  legend("topright", legend = c("Best","Final"),
         bty = "n", pch = 20, col = c("blue","firebrick"),
         text.col = "black", cex = 1.2, pt.cex = 2)
  dev.off()
} else {
  # ggplot fallback radar (polar)
  rad_long <- reshape2::melt(metrics_df, id.vars = "Metric",
                             variable.name = "Model", value.name = "Value")
  rad_long$Metric <- factor(rad_long$Metric, levels = metrics_df$Metric)
  p_radar <- ggplot(rad_long, aes(x = Metric, y = Value, group = Model, color = Model)) +
    geom_polygon(aes(fill = Model), alpha = 0.2, linewidth = 1.2) +
    geom_point(size = 3) +
    coord_polar() +
    scale_color_manual(values = c("Best" = "blue", "Final" = "firebrick")) +
    scale_fill_manual(values  = c("Best" = alpha("blue", 0.25), "Final" = alpha("firebrick", 0.25))) +
    labs(title = "(B) Metrics Radar - Best vs Final", x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(file_radar, p_radar, width = 12, height = 9, dpi = 200)
}

## (C) Confusion matrix - Best
p_cm_best  <- plot_cm_gg(xgb_cm, "(C) Confusion Matrix - Best (Test)", "dodgerblue")
file_cm_best <- file.path(out_dir, "TMP_panel_CM_best.png")
ggsave(file_cm_best, p_cm_best, width = 8, height = 6, dpi = 200)

## (D) Confusion matrix - Final
p_cm_final <- plot_cm_gg(xgb_cm_full, "(D) Confusion Matrix - Final (Full)", "tomato")
file_cm_final <- file.path(out_dir, "TMP_panel_CM_final.png")
ggsave(file_cm_final, p_cm_final, width = 8, height = 6, dpi = 200)

# ============================================
# 5) Assemble 2×2 figure
# ============================================
gA <- to_grob(file_roc)
gB <- to_grob(file_radar)
gC <- to_grob(file_cm_best)
gD <- to_grob(file_cm_final)

file_final <- file.path(out_dir, "Figure_Summary_Best_vs_Final.png")
png(file_final, width = 3200, height = 2400, res = 250)
grid.arrange(gA, gB, gC, gD, ncol = 2)
dev.off()

cat("??? Saved combined multi-panel figure to:\n", file_final, "\n")
cat(sprintf("AUCs - Best: %.3f | Final: %.3f\n", auc_best, auc_final))
cat(sprintf("Thresholds - Best: %.3f | Final: %.3f\n", thr_best, thr_final))

------------------------------------------------------------------------------------------
  #Heatmap
  # ============================
# Signed, annotated correlation heatmap (Top-K features)
# ============================

# What this uses from your session:
#   xgb_model         # caret::train object (method="xgbTree")
#   X_train_matrix    # training matrix used to fit
#   X_test_matrix     # test matrix for evaluation
# If any are missing, load them first.

# ---- packages ----
need <- c("xgboost","ggplot2","reshape2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(xgboost)
library(ggplot2)
library(reshape2)

# ---- settings ----
TOP_K       <- 20             # how many top features to show (try 12/20/30)
REMOVE_MASS <- FALSE          # TRUE to drop MONOISOTOPIC.MASS from the plot
OUT_DIR     <- "C:/Users/91962/Predicting Chemicals/Plots(Models)"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---- importance ??? pick top-K feature names ----
feature_names <- colnames(X_train_matrix)
imp <- xgb.importance(model = xgb_model$finalModel, feature_names = feature_names)
imp_top <- head(imp, TOP_K)
top_feats <- imp_top$Feature
if (REMOVE_MASS) top_feats <- setdiff(top_feats, c("MONOISOTOPIC.MASS","MONOISOTOPIC.MASS."))

# ---- build correlation on TEST split (what you report elsewhere) ----
X_top <- as.matrix(X_test_matrix[, top_feats, drop = FALSE])

# Pearson (for binary fingerprints this equals phi). Keep SIGNED correlations.
C <- suppressWarnings(cor(X_top, use = "pairwise.complete.obs", method = "pearson"))
C[is.na(C)] <- 0
diag(C) <- 1

# Order variables for nicer blocks (cluster by absolute similarity)
D   <- as.dist(pmax(0, 1 - abs(C)))
ord <- tryCatch(hclust(D, method = "average")$order,
                error = function(e) seq_len(ncol(C)))
C   <- C[ord, ord]

# Long format + text color based on magnitude (white text on strong colors)
m <- reshape2::melt(C, varnames = c("Var1","Var2"), value.name = "r")
m$txt <- sprintf("%.2f", m$r)
m$txt_col <- ifelse(abs(m$r) > 0.5, "white", "black")

# ---- plot ----
p <- ggplot(m, aes(Var2, Var1, fill = r)) +
  geom_tile() +
  geom_text(aes(label = txt, color = txt_col), size = 3.2, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradient2(limits = c(-1, 1),
                       low = "#2b8cbe", mid = "white", high = "#cb181d",
                       midpoint = 0, name = "Correlation r") +
  coord_fixed() +
  labs(
    title    = paste0("Correlation heatmap - Top-", length(top_feats), " important features"),
    subtitle = "Signed Pearson correlations (computed on test set)",
    x = NULL, y = NULL,
    caption  = "Red = positive correlation; Blue = negative; values inside cells"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid  = element_blank()
  )

out_png <- file.path(OUT_DIR, sprintf("Heatmap_SignedCorr_Top%d_Features.png", length(top_feats)))
ggsave(out_png, p, width = 10, height = 8.5, dpi = 300, bg = "white")
cat("??? Saved heatmap to:", out_png, "\n")

# (optional) also save the importance table used in this figure
write.csv(imp_top, file.path(OUT_DIR, sprintf("Top%d_Features_Importance.csv", length(top_feats))), row.names = FALSE)
--------------------------------------------------------------------------------------------
  # ============================
# Figures tied to your model
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

# ---- 0) Inputs this uses (already in your session) ----
# X_train_matrix, X_test_matrix, y_test, xgb_model (caret::train object)
# If xgb_pred_prob not in memory, compute it:
if (!exists("xgb_pred_prob")) {
  xgb_pred_prob <- predict(xgb_model, newdata = X_test_matrix, type = "prob")[, "Active"]
}

# Feature names for importance
feature_names <- colnames(X_train_matrix)

# ============================
# Panel A - probability distribution (test set)
# ============================
roc_best <- pROC::roc(response = y_test, predictor = xgb_pred_prob, quiet = TRUE)
thr_best <- as.numeric(pROC::coords(roc_best, "best", ret = "threshold"))

df_prob <- data.frame(
  Prob = xgb_pred_prob,
  Class = factor(ifelse(y_test == 1, "Active", "Inactive"), levels = c("Inactive","Active"))
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
TOP_K <- 12  # change to 20/25 if you want a denser heatmap

# importance from the final xgboost booster inside caret model
imp <- xgboost::xgb.importance(model = xgb_model$finalModel, feature_names = feature_names)
imp_top <- head(imp, TOP_K)
top_feats <- imp_top$Feature

# correlate on the TEST matrix (what produced the probs)
X_top <- as.matrix(X_test_matrix[, top_feats, drop = FALSE])

# pearson on binaries ??? phi; clean and cluster
C <- suppressWarnings(cor(X_top, use = "pairwise.complete.obs"))
C[is.na(C)] <- 0; diag(C) <- 1
D <- as.dist(pmax(0, 1 - C)); D[!is.finite(D)] <- 0
ord <- tryCatch(hclust(D, method = "average")$order, error = function(e) seq_len(ncol(C)))
Ck <- C[ord, ord]

# upper triangle only for readability
Ck[lower.tri(Ck)] <- NA
hm <- reshape2::melt(Ck, varnames = c("Var1","Var2"), value.name = "corr", na.rm = TRUE)

pB <- ggplot(hm, aes(Var2, Var1, fill = corr)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0,1), low = "white", high = "#08519C", name = "Correlation") +
  coord_fixed() +
  labs(title = "B", subtitle = paste0("Correlation among top-", TOP_K, " important features"),
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

cat("??? Saved:\n",
    file.path(out_dir, "PanelA_ProbDist_Test.png"), "\n",
    file.path(out_dir, "PanelB_TopFeatures_Correlation.png"), "\n",
    file.path(out_dir, "Figure_ProbDist_and_TopFeatureCorrelation.png"), "\n", sep = "")


--------------------------------------------------------
# ================================================================
# XGBoost classification: Calibration (Original vs Platt vs Isotonic)
# ================================================================

out_dir <- "C:/Users/91962/Predicting Chemicals/Plots(Modelsb)"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(pROC); library(gridExtra)
})

# ---------- helpers ----------
as01 <- function(y){
  if (is.data.frame(y) || is.matrix(y)) y <- y[,1]
  if (is.factor(y)) {
    lv <- tolower(levels(y)); yc <- tolower(as.character(y))
    if (all(c("active","inactive") %in% lv)) return(as.numeric(yc=="active"))
    return(suppressWarnings(as.numeric(as.character(y))))
  }
  if (is.character(y)) { ylc <- tolower(y); 
  if (all(ylc %in% c("active","inactive"))) return(as.numeric(ylc=="active"))
  return(suppressWarnings(as.numeric(y)))
  }
  as.numeric(y)
}

# equal-count (quantile) bins -> stable points
bin_equal_count <- function(p, y, bins = 10){
  brks <- unique(quantile(p, probs = seq(0,1,length.out=bins+1), na.rm=TRUE))
  if (length(brks) < 3) brks <- seq(0,1,length.out=3)
  b <- cut(p, breaks = brks, include.lowest = TRUE, right = TRUE)
  tibble(prob=p, y=y, bin=b) |>
    group_by(bin) |>
    summarise(mean_pred = mean(prob),
              obs_rate  = mean(y),
              n = n(), .groups="drop")
}

# Expected Calibration Error (ECE, 10 bins)
ece10 <- function(y, p, k=10){
  d <- bin_equal_count(p, y, bins = k)
  sum(d$n/length(y) * abs(d$mean_pred - d$obs_rate))
}

# Isotonic calibration fitted on train
fit_isotonic <- function(p_train, y_train){
  ir <- stats::isoreg(p_train, y_train)  # monotone fit of y on p
  df <- data.frame(x = ir$x, y = ir$yf)
  df <- df[order(df$x), ]
  df <- df[!duplicated(df$x), ]
  f  <- stats::approxfun(df$x, pmin(pmax(df$y,0),1), method="linear", rule=2)
  function(p) pmin(pmax(f(p),0),1)
}

# Single panel with Original points + Platt + Isotonic lines
make_panel <- function(y, p_orig, p_platt, p_iso, title_text, band_col="#1f78b4", bins=10, base_size=18){
  d_orig  <- bin_equal_count(p_orig,  y, bins)
  d_platt <- bin_equal_count(p_platt, y, bins)
  d_iso   <- bin_equal_count(p_iso,   y, bins)
  
  grid  <- data.frame(x = seq(0,1,length.out=1000))
  p <- ggplot() +
    geom_ribbon(data=grid, aes(x, ymin=pmax(0,x-0.30), ymax=pmin(1,x+0.30)), fill=band_col, alpha=.18) +
    geom_ribbon(data=grid, aes(x, ymin=pmax(0,x-0.20), ymax=pmin(1,x+0.20)), fill=band_col, alpha=.26) +
    geom_ribbon(data=grid, aes(x, ymin=pmax(0,x-0.10), ymax=pmin(1,x+0.10)), fill=band_col, alpha=.36) +
    geom_abline(slope=1, intercept=0, linewidth=1.1) +
    # Original points/line
    geom_point(data=d_orig, aes(obs_rate, mean_pred), size=3.4, color="black") +
    geom_path (data=d_orig, aes(obs_rate, mean_pred), linewidth=.9, color="black") +
    # Calibrated lines
    geom_line (data=d_platt, aes(obs_rate, mean_pred, color="Platt"), linewidth=1.1) +
    geom_line (data=d_iso,   aes(obs_rate, mean_pred, color="Isotonic"), linewidth=1.1, linetype="22") +
    scale_color_manual("", values=c("Platt"="#0055FF","Isotonic"="#D81B60")) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    labs(title=title_text, x="Observed Active rate (binned)", y="Predicted probability (Active)") +
    annotate("text", x=.14, y=.03, label="Overpredicted", angle=45, alpha=.7, size=5) +
    annotate("text", x=.03, y=.14, label="Underpredicted", angle=45, alpha=.7, size=5) +
    theme_classic(base_size = base_size) +
    theme(panel.border = element_rect(color="black", fill=NA, linewidth=.7),
          plot.title = element_text(face="bold"), legend.position = "top")
  p
}

# ---------- predictions & calibration fits ----------
stopifnot(exists("xgb_model"), exists("X_train_matrix"), exists("X_test_matrix"),
          exists("y_train"), exists("y_test"))

p_train <- predict(xgb_model, newdata=X_train_matrix, type="prob")[,"Active"]
p_test  <- predict(xgb_model, newdata=X_test_matrix,  type="prob")[,"Active"]

y_tr <- as01(y_train); y_te <- as01(y_test)

# Platt (sigmoid) fit on train
platt_fit <- glm(y_tr ~ p_train, family = binomial())
p_train_platt <- predict(platt_fit, newdata = data.frame(p_train=p_train), type="response")
p_test_platt  <- predict(platt_fit, newdata = data.frame(p_train=p_test),  type="response")

# Isotonic fit on train
iso_fun <- fit_isotonic(p_train, y_tr)
p_train_iso <- iso_fun(p_train)
p_test_iso  <- iso_fun(p_test)

# ---------- metrics ----------
metrics <- function(y, p){
  auc <- as.numeric(pROC::auc(pROC::roc(y, p, levels=c(0,1), direction="<", quiet=TRUE)))
  brier <- mean((p - y)^2)
  ece <- ece10(y, p, k=10)
  c(AUC = auc, Brier = brier, ECE10 = ece)
}
cat("\nVALIDATION metrics\n",
    "\nOriginal :", round(metrics(y_te, p_test), 4),
    "\nPlatt    :", round(metrics(y_te, p_test_platt), 4),
    "\nIsotonic :", round(metrics(y_te, p_test_iso), 4), "\n")

# ---------- panels ----------
p_train_panel <- make_panel(y_tr, p_train, p_train_platt, p_train_iso,
                            title_text="Train", band_col="#1f78b4", bins=10)
p_test_panel  <- make_panel(y_te, p_test, p_test_platt, p_test_iso,
                            title_text="Validation", band_col="#33a02c", bins=10)

combo <- gridExtra::arrangeGrob(p_train_panel, p_test_panel, ncol=2)
ggsave(file.path(out_dir, "Calibration_Train_vs_Validation_XGB_Orig-Platt-Iso.png"),
       combo, width=14.5, height=6.5, dpi=320)

cat("Saved: ", file.path(out_dir, "Calibration_Train_vs_Validation_XGB_Orig-Platt-Iso.png"), "\n")
--------------------------------------------------------------------------------------------------
  
  # ================================================================
# Reliability (Calibration) Figure for XGBoost Classification
# - Predicted on x; Observed on y (standard reading)
# - Continuous Platt + Isotonic lines
# - Decile points with n and 95% Wilson CIs
# - AUC, Brier, ECE(10) annotated
# - Validation panel shows your chosen threshold (vertical line)
# ================================================================
out_dir <- "C:/Users/91962/Predicting Chemicals/Plots(Modelsb)"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(pROC); library(gridExtra)
})

# -------- helpers --------
as01 <- function(y){
  if (is.data.frame(y) || is.matrix(y)) y <- y[,1]
  if (is.factor(y)) {
    lv <- tolower(levels(y)); yc <- tolower(as.character(y))
    if (all(c("active","inactive") %in% lv)) return(as.numeric(yc=="active"))
    return(suppressWarnings(as.numeric(as.character(y))))
  }
  if (is.character(y)) {
    ylc <- tolower(y)
    if (all(ylc %in% c("active","inactive"))) return(as.numeric(ylc=="active"))
    return(suppressWarnings(as.numeric(y)))
  }
  as.numeric(y)
}

wilson_ci <- function(k, n, conf = 0.95){
  if (n == 0) return(c(NA,NA))
  z <- qnorm(1 - (1 - conf)/2)
  p <- k/n; denom <- 1 + z^2/n
  centre <- (p + z^2/(2*n)) / denom
  half   <- z * sqrt((p*(p-1)/(-n) + z^2/(4*n^2))) / denom
  c(max(0, centre - half), min(1, centre + half))
}

bin_equal_count <- function(p, y, bins = 10){
  brks <- unique(quantile(p, probs = seq(0,1,length.out=bins+1), na.rm=TRUE))
  if (length(brks) < 3) brks <- seq(0,1,length.out=3)
  cut(p, breaks = brks, include.lowest = TRUE, right = TRUE)
}

ece10 <- function(y, p, k=10){
  b <- bin_equal_count(p, y, bins=k)
  d <- tibble(p=p, y=y, b=b) |>
    group_by(b) |>
    summarise(pp = mean(p), oo = mean(y), n=n(), .groups="drop")
  sum(d$n/length(y) * abs(d$pp - d$oo))
}

fit_isotonic <- function(p_train, y_train){
  ir <- stats::isoreg(p_train, y_train)
  df <- data.frame(x = ir$x, y = ir$yf)
  df <- df[order(df$x), ]; df <- df[!duplicated(df$x), ]
  f  <- stats::approxfun(df$x, pmin(pmax(df$y,0),1), method="linear", rule=2)
  function(p) pmin(pmax(f(p),0),1)
}

panel_reliability <- function(y, p_raw, platt_fun, iso_fun,
                              title_text, band_col="#1f78b4",
                              bins=10, base_size=18, vline = NA,
                              metrics_label = NULL){
  
  # bin points (x = mean predicted, y = observed rate)
  b <- bin_equal_count(p_raw, y, bins)
  pts <- tibble(p=p_raw, y=y, b=b) |>
    group_by(b) |>
    summarise(x = mean(p), y = mean(y), k = sum(y), n = n(), .groups="drop") |>
    rowwise() |>
    mutate(ci = list(wilson_ci(k, n, 0.95)),
           y_lo = ci[1], y_hi = ci[2]) |>
    ungroup()
  
  # continuous curves across 0..1
  grid <- data.frame(x = seq(0,1,length.out=1001))
  grid$platt <- platt_fun(grid$x)
  grid$iso   <- iso_fun(grid$x)
  
  # agreement ribbons
  band <- data.frame(x = grid$x)
  p <- ggplot() +
    geom_ribbon(data=band, aes(x, ymin=pmax(0, x-0.30), ymax=pmin(1, x+0.30)),
                fill = band_col, alpha=.18) +
    geom_ribbon(data=band, aes(x, ymin=pmax(0, x-0.20), ymax=pmin(1, x+0.20)),
                fill = band_col, alpha=.26) +
    geom_ribbon(data=band, aes(x, ymin=pmax(0, x-0.10), ymax=pmin(1, x+0.10)),
                fill = band_col, alpha=.36) +
    geom_abline(slope=1, intercept=0, linewidth=1.1, color="black") +
    # CI bars + points (deciles)
    geom_errorbar(data=pts, aes(x=x, ymin=y_lo, ymax=y_hi), width=0.012, linewidth=0.7) +
    geom_point(data=pts, aes(x=x, y=y), size=3.4, color="black") +
    geom_text(data=pts, aes(x=x, y=y, label=n), vjust=-1.1, size=4) +
    # Platt / Isotonic continuous lines
    geom_line(data=grid, aes(x=x, y=platt, color="Platt"), linewidth=1.1) +
    geom_line(data=grid, aes(x=x, y=iso,   color="Isotonic"), linewidth=1.1, linetype="22") +
    scale_color_manual("", values=c("Platt"="#0055FF","Isotonic"="#D81B60")) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    labs(title = title_text,
         x = "Predicted probability (Active)",
         y = "Observed Active rate") +
    theme_classic(base_size = base_size) +
    theme(panel.border = element_rect(color="black", fill=NA, linewidth=.7),
          plot.title   = element_text(face="bold"),
          legend.position = "top")
  
  if (!is.na(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = "dashed", linewidth = 0.9) +
      annotate("text", x = vline, y = 0.04, label = sprintf("Threshold = %.2f", vline),
               angle = 90, vjust = -0.5, size = 4)
  }
  
  if (!is.null(metrics_label)) {
    p <- p + annotate("label", x = 0.02, y = 0.92, hjust = 0, vjust = 1,
                      label = metrics_label, size = 4.2, label.size = 0, alpha = 0.95)
  }
  p
}

# -------- predictions + calibration fits --------
stopifnot(exists("xgb_model"), exists("X_train_matrix"), exists("X_test_matrix"),
          exists("y_train"), exists("y_test"))

p_tr <- predict(xgb_model, newdata = X_train_matrix, type="prob")[,"Active"]
p_te <- predict(xgb_model, newdata = X_test_matrix,  type="prob")[,"Active"]
y_tr <- as01(y_train); y_te <- as01(y_test)

# Platt (sigmoid) fitted on train
platt_fit <- glm(y_tr ~ p_tr, family = binomial())
platt_fun <- function(p){ predict(platt_fit, newdata = data.frame(p_tr=p), type="response") }

# Isotonic fitted on train
iso_fun <- fit_isotonic(p_tr, y_tr)

# Metrics on VALIDATION
thr_val <- as.numeric(coords(roc(y_te, p_te, levels=c(0,1), direction="<", quiet=TRUE),
                             "best", ret="threshold"))
m_vec <- function(y, p){
  auc <- as.numeric(auc(roc(y, p, levels=c(0,1), direction="<", quiet=TRUE)))
  brier <- mean((p - y)^2)
  ece <- ece10(y, p, k=10)
  c(AUC = auc, Brier = brier, ECE10 = ece)
}
orig_v <- m_vec(y_te, p_te)
plat_v <- m_vec(y_te, platt_fun(p_te))
iso_v  <- m_vec(y_te, iso_fun(p_te))

label_v <- sprintf("Validation metrics:\nOriginal  AUC=%.3f  Brier=%.3f  ECE10=%.3f\nPlatt     AUC=%.3f  Brier=%.3f  ECE10=%.3f\nIsotonic  AUC=%.3f  Brier=%.3f  ECE10=%.3f",
                   orig_v[1], orig_v[2], orig_v[3],
                   plat_v[1], plat_v[2], plat_v[3],
                   iso_v[1],  iso_v[2],  iso_v[3])

# (optional) Train metrics label if you want both
orig_t <- m_vec(y_tr, p_tr); plat_t <- m_vec(y_tr, platt_fun(p_tr)); iso_t <- m_vec(y_tr, iso_fun(p_tr))
label_t <- sprintf("Train metrics:\nOriginal  AUC=%.3f  Brier=%.3f  ECE10=%.3f\nPlatt     AUC=%.3f  Brier=%.3f  ECE10=%.3f\nIsotonic  AUC=%.3f  Brier=%.3f  ECE10=%.3f",
                   orig_t[1], orig_t[2], orig_t[3],
                   plat_t[1], plat_t[2], plat_t[3],
                   iso_t[1],  iso_t[2],  iso_t[3])

# -------- build panels --------
p_train <- panel_reliability(y_tr, p_tr, platt_fun, iso_fun,
                             title_text = "Train", band_col = "#1f78b4",
                             bins = 10, base_size = 19, vline = NA,
                             metrics_label = label_t)

p_valid <- panel_reliability(y_te, p_te, platt_fun, iso_fun,
                             title_text = "Validation", band_col = "#33a02c",
                             bins = 10, base_size = 19, vline = thr_val,
                             metrics_label = label_v)

combo <- gridExtra::arrangeGrob(p_train, p_valid, ncol = 2)
ggsave(file.path(out_dir, "Reliability_Train_vs_Validation_CALIBRATED.png"),
       combo, width = 15, height = 6.6, dpi = 320)

cat("Saved: ", file.path(out_dir, "Reliability_Train_vs_Validation_CALIBRATED.png"), "\n")

# -------- SUGGESTED CAPTION (copy into your dissertation) --------
cat("\nFIGURE CAPTION SUGGESTION:\n",
    'Reliability of the XGBoost classifier on the training (left) and validation (right) sets.
Each dot is a decile of predicted probability (x-axis) with the corresponding observed Active rate (y-axis); numbers above the dots are bin sizes and horizontal bars show 95% Wilson confidence intervals.
The diagonal line indicates perfect calibration, with shaded bands for ±0.10/0.20/0.30 deviation.
Coloured curves are probability-calibration mappings learned on the training set: Platt (blue, logistic) and Isotonic (magenta, monotone).
The validation panel also marks the selected decision threshold (dashed vertical line; Youden from ROC).
Text insets report AUC (discrimination), Brier score (overall accuracy of probabilities), and ECE(10) (average calibration error over 10 equal-count bins).
Curves closer to the diagonal and lower Brier/ECE indicate better-calibrated probabilities without necessarily changing AUC.',
    "\n")
