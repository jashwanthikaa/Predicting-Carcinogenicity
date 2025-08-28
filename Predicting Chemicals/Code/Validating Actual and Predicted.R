# ================================================
# Confusion Matrix Plot: Actual vs Predicted
# ================================================

# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)

# -----------------------------
# 1. Load the data
# -----------------------------
actual_df <- read_excel("C:/Users/91962/Predicting Chemicals/Validation1.xlsx")
predicted_df <- read.csv("C:/Users/91962/Predicting Chemicals/Sirius_predictions_names.csv", stringsAsFactors = FALSE)

# -----------------------------
# 2. Clean and standardise keys
# -----------------------------
actual_df$MOLECULAR_FORMULA <- toupper(trimws(actual_df$MOLECULAR_FORMULA))
predicted_df$predform <- toupper(trimws(predicted_df$predform))

# Clean outcome labels to avoid mismatches
actual_df$HIT.CALL <- trimws(actual_df$HIT.CALL)
predicted_df$Carcinogenicity_predicted <- trimws(predicted_df$Carcinogenicity_predicted)

# -----------------------------
# 3. Merge datasets
# -----------------------------
merged_df <- merge(actual_df, predicted_df,
                   by.x = "MOLECULAR_FORMULA",
                   by.y = "predform",
                   all = FALSE)

# -----------------------------
# 4. Create confusion matrix
# -----------------------------
lvl <- c("Active", "Inactive")
y_true <- factor(merged_df$HIT.CALL, levels = lvl)
y_pred <- factor(merged_df$Carcinogenicity_predicted, levels = lvl)

tab <- table(Predicted = y_pred, Actual = y_true)
print(tab)

# Extract counts
TP <- tab["Active",   "Active"]
FN <- tab["Inactive", "Active"]
FP <- tab["Active",   "Inactive"]
TN <- tab["Inactive", "Inactive"]

# Accuracy
accuracy <- (TP + TN) / sum(tab)
cat(sprintf("\nAccuracy: %.2f%%\n", 100 * accuracy))

# -----------------------------
# 5. Function to plot confusion matrix
# -----------------------------
plot_conf_mat <- function(TP, FN, FP, TN,
                          x_title = "Fingerprint Present / Absent",
                          y_title = "rcdk calculated Fingerprints",
                          plot_title = "SIRIUS Fingerprints") {
  
  df <- tibble::tibble(
    Predicted = factor(c("Fingerprint Present","Fingerprint Absent",
                         "Fingerprint Present","Fingerprint Absent"),
                       levels = c("Fingerprint Present","Fingerprint Absent")),
    Actual    = factor(c("Fingerprint Present","Fingerprint Present",
                         "Fingerprint Absent","Fingerprint Absent"),
                       levels = c("Fingerprint Present","Fingerprint Absent")),
    Count     = c(TP, FN, FP, TN),
    Kind      = c("True Positive","False Negative","False Positive","True Negative")
  ) %>%
    mutate(Fill = ifelse(Kind %in% c("True Positive","True Negative"), "correct", "error"))
  
  gg <- ggplot(df, aes(Predicted, Actual, fill = Fill)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = paste0(Kind, "\n", Count)),
              size = 7, fontface = "bold") +
    scale_fill_manual(values = c(correct = "#86e686", error = "#f27d7d"), guide = "none") +
    labs(x = x_title, y = y_title, title = plot_title) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  print(gg)
}

# -----------------------------
# 6. Plot confusion matrix
# -----------------------------
plot_conf_mat(TP, FN, FP, TN,
              x_title = "Fingerprint Present / Absent",
              y_title = "rcdk calculated Fingerprints",
              plot_title = "SIRIUS Fingerprints")
