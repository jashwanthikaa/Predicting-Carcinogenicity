# ============================
# Correlation heatmap of TOP-K important features
# ============================

# 1) Packages
need <- c("xgboost","ggplot2","reshape2")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(xgboost)
library(ggplot2)
library(reshape2)

# 2) Settings
TOP_K   <- 12   # change to 20/25 for a denser map
ABS_COR <- TRUE # TRUE ??? show |correlation| on 0..1 scale (looks like your example)
out_dir <- "C:/Users/91962/Predicting Chemicals/Plots(Models)"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_png <- file.path(out_dir, sprintf("Heatmap_Top%d_FeatureCorrelation.png", TOP_K))

# 3) Get importance from the trained model (caret::train with method='xgbTree')
feature_names <- colnames(X_train_matrix)
imp <- xgb.importance(model = xgb_model$finalModel, feature_names = feature_names)
imp_top <- head(imp, TOP_K)
top_feats <- imp_top$Feature

# 4) Build correlation on the TEST matrix (the data you used for evaluation)
X_top <- as.matrix(X_test_matrix[, top_feats, drop = FALSE])

# Pearson on binaries ??? phi; clean & cluster
C <- suppressWarnings(cor(X_top, use = "pairwise.complete.obs", method = "pearson"))
if (ABS_COR) C <- abs(C)
C[is.na(C)] <- 0
diag(C) <- 1

# Cluster for nicer blocks
D <- as.dist(pmax(0, 1 - C)); D[!is.finite(D)] <- 0
ord <- tryCatch(hclust(D, method = "average")$order, error = function(e) seq_len(ncol(C)))
Ck <- C[ord, ord]

# Upper triangle only (readable)
Ck[lower.tri(Ck)] <- NA
hm <- reshape2::melt(Ck, varnames = c("Var1","Var2"), value.name = "corr", na.rm = TRUE)

# 5) Plot & save
p <- ggplot(hm, aes(Var2, Var1, fill = corr)) +
  geom_tile() +
  (if (ABS_COR) scale_fill_gradient(limits = c(0,1), low = "white", high = "#08519C",
                                    name = if (ABS_COR) "Abs. corr" else "Correlation")
   else          scale_fill_gradient2(limits = c(-1,1), low = "#67000D", mid = "white",
                                      high = "#08519C", name = "Correlation")) +
  coord_fixed() +
  labs(title = sprintf("Correlation among top-%d important features", TOP_K),
       subtitle = "Ranked by XGBoost gain",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank())

ggsave(out_png, p, width = 8.5, height = 7, dpi = 300)
cat("??? Saved heatmap to:", out_png, "\n")

# (Optional) also save the importance table used
write.csv(imp_top, file.path(out_dir, sprintf("Top%d_Features_Importance.csv", TOP_K)), row.names = FALSE)
------------------------------------------------
  
library(xgboost); library(ggplot2); library(reshape2)

TOP_K   <- 12
ABS_COR <- TRUE                 # FALSE = show signed correlations (-1..1)
out_png <- "C:/Users/91962/Predicting Chemicals/Plots(Models)/Heatmap_Top12_FeatureCorrelation_labeled.png"

# Importance ??? top features
feature_names <- colnames(X_train_matrix)
imp      <- xgb.importance(model = xgb_model$finalModel, feature_names = feature_names)
imp_top  <- head(imp, TOP_K)
top_feats <- imp_top$Feature

# OPTIONAL: drop mass from the heatmap if you prefer
top_feats <- setdiff(top_feats, c("MONOISOTOPIC.MASS", "MONOISOTOPIC.MASS."))

# Correlations on the TEST data
X_top <- as.matrix(X_test_matrix[, top_feats, drop = FALSE])
C <- suppressWarnings(cor(X_top, use = "pairwise.complete.obs", method = "pearson"))
if (ABS_COR) C <- abs(C)
C[is.na(C)] <- 0; diag(C) <- 1

# Cluster for nicer blocks
D <- as.dist(pmax(0, 1 - C)); D[!is.finite(D)] <- 0
ord <- tryCatch(hclust(D, method = "average")$order, error = function(e) seq_len(ncol(C)))
Ck  <- C[ord, ord]
Ck[lower.tri(Ck)] <- NA   # upper triangle only

hm <- reshape2::melt(Ck, varnames = c("Var1","Var2"), value.name = "corr", na.rm = TRUE)

p <- ggplot(hm, aes(Var2, Var1, fill = corr)) +
  geom_tile() +
  # add numeric labels (remove this line if you want a cleaner look)
  geom_text(aes(label = sprintf(if (ABS_COR) "%.2f" else "%.2f", corr)), size = 3, na.rm = TRUE) +
  (if (ABS_COR)
    scale_fill_gradient(limits = c(0,1), low = "white", high = "#08519C", name = "Abs. corr")
   else
     scale_fill_gradient2(limits = c(-1,1), low = "#67000D", mid = "white", high = "#08519C", name = "Correlation")) +
  coord_fixed() +
  labs(title = sprintf("Correlation among top-%d important features", length(top_feats)),
       subtitle = "Ranked by XGBoost gain (computed on test set)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank())

ggsave(out_png, p, width = 9, height = 8, dpi = 300, bg = "white")
cat("Saved:", out_png, "\n")
---------------------------------------------------
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
