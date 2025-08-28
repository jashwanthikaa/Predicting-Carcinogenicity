# Load libraries
library(readxl)
library(caret)
library(dplyr)
library(xgboost)
library(Matrix)
library(pROC)

# Load dataset
df <- read_excel("C:/Users/91962/Predicting Chemicals/FinalData.xlsx")

# Select fingerprint and mass columns
fp_cols <- grep("^Un", names(df), value = TRUE)
mass_col <- "MONOISOTOPIC.MASS"
X <- df %>% select(all_of(c(fp_cols, mass_col)))
y <- ifelse(tolower(df$`HIT.CALL`) == "active", 1, 0)

# Remove rows with missing values
non_na_rows <- complete.cases(X)
X <- X[non_na_rows, ]
y <- y[non_na_rows]
df_clean <- df[non_na_rows, ]

# Remove near-zero variance features
nzv <- nearZeroVar(X)
if (length(nzv) > 0) {
  X <- X[, -nzv]
}

# Train/test split
set.seed(123)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
y_train <- y[train_index]
X_test <- X[-train_index, ]
y_test <- y[-train_index]

# Convert to matrix
X_train_matrix <- data.matrix(X_train)
X_test_matrix <- data.matrix(X_test)
y_train_factor <- factor(ifelse(y_train == 1, "Active", "Inactive"))

# Metadata for test set
metadata_cols <- c("Name", "CASRN", "SMILES", "MONOISOTOPIC.MASS", "DTXSID")
metadata_test <- df_clean[-train_index, metadata_cols]

# Define train control with repeated CV
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Calculate class imbalance weight
scale_weight <- sum(y_train == 0) / sum(y_train == 1)

# Define tuning grid
xgb_grid <- expand.grid(
  nrounds = 300,
  max_depth = c(6, 9, 12),
  eta = c(0.01, 0.05, 0.1),
  gamma = c(0, 1, 5),
  colsample_bytree = c(0.6, 0.8, 1),
  min_child_weight = c(1, 3, 5),
  subsample = c(0.6, 0.8, 1)
)

# Train the model
set.seed(123)
xgb_model <- train(
  x = X_train_matrix,
  y = y_train_factor,
  trControl = train_control,
  tuneGrid = xgb_grid,
  method = "xgbTree",
  metric = "ROC",
  verbose = TRUE,
  scale_pos_weight = scale_weight
)

# ----------------------------
# Best hyperparameters
# ----------------------------
best_params <- xgb_model$bestTune
cat("\n Best Hyperparameters:\n")
print(best_params)

# ----------------------------
# Make predictions on test set
# ----------------------------
xgb_pred_prob <- predict(xgb_model, newdata = X_test_matrix, type = "prob")[, "Active"]

# Find optimal threshold from ROC
roc_obj <- roc(y_test, xgb_pred_prob)
best_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))
optimal_threshold <- as.numeric(best_coords["threshold"])
cat("\n Optimal threshold based on ROC:", optimal_threshold, "\n")

# Final binary predictions using optimized threshold
xgb_pred <- ifelse(xgb_pred_prob > optimal_threshold, 1, 0)

# ----------------------------
# Confusion Matrix
# ----------------------------
y_test_factor <- factor(y_test, levels = c(0, 1))
xgb_pred_factor <- factor(xgb_pred, levels = c(0, 1))

xgb_cm <- confusionMatrix(xgb_pred_factor, y_test_factor, positive = "1")

cat("\n Confusion Matrix:\n")
print(xgb_cm)

# ----------------------------
# Class distribution breakdown
# ----------------------------
cm_table <- xgb_cm$table

actual_active <- sum(y_test == 1)
actual_inactive <- sum(y_test == 0)
predicted_active <- sum(xgb_pred == 1)
predicted_inactive <- sum(xgb_pred == 0)

cat("\n Confusion Matrix Counts (Test Set):\n")
print(cm_table)

cat("\n Class Distribution:\n")
cat("  Actual Active     :", actual_active, "\n")
cat("  Actual Inactive   :", actual_inactive, "\n")
cat("  Predicted Active  :", predicted_active, "\n")
cat("  Predicted Inactive:", predicted_inactive, "\n")

# ----------------------------
# Save prediction results
# ----------------------------
results <- data.frame(
  Name = metadata_test$Name,
  DTXSID = metadata_test$DTXSID,
  CASRN = metadata_test$CASRN,
  SMILES = metadata_test$SMILES,
  MONOISOTOPIC_MASS = metadata_test$MONOISOTOPIC.MASS,
  Actual = y_test,
  Predicted = xgb_pred,
  Probability = xgb_pred_prob
)

write.csv(results, "xgb_model_predictions_optimized300(1).csv", row.names = FALSE)

# Plot ROC
plot(roc_obj, col = "blue", main = "ROC Curve for XGBoost Model")
