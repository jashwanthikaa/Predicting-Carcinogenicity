# Load libraries
library(readxl)
library(caret)
library(dplyr)
library(xgboost)
library(Matrix)
library(pROC)

# Load dataset
df <- read_excel("C:/Users/91962/Predicting Chemicals/FinalData.xlsx")

# Select fingerprint and MONOISOTOPIC.MASS columns
fp_cols <- grep("^Un", names(df), value = TRUE)
MONOISOTOPIC.MASS_col <- "MONOISOTOPIC.MASS"
X <- df %>% select(all_of(c(fp_cols, MONOISOTOPIC.MASS_col)))
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

# Metadata for test results
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

# Calculate scale_pos_weight for imbalance
scale_weight <- sum(y_train == 0) / sum(y_train == 1)

# Expanded tuning grid
xgb_grid <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.01,
  gamma = 5,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  subsample = 0.6
)

# Train XGBoost
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

# Best model details
cat("\nBest Hyperparameters:\n")
print(xgb_model$bestTune)

# Predict probabilities on test data
xgb_pred_prob <- predict(xgb_model, newdata = X_test_matrix, type = "prob")[, "Active"]

# ROC-based threshold optimization
roc_obj <- roc(y_test, xgb_pred_prob)
best_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))
optimal_threshold <- as.numeric(best_coords["threshold"])
cat("\nOptimal threshold based on ROC:", optimal_threshold, "\n")

# Final prediction using optimized threshold
xgb_pred <- ifelse(xgb_pred_prob > optimal_threshold, 1, 0)

# Confusion matrix (fixing levels)
y_test_factor <- factor(y_test, levels = c(0, 1))
xgb_pred_factor <- factor(xgb_pred, levels = c(0, 1))

xgb_cm <- confusionMatrix(xgb_pred_factor, y_test_factor, positive = "1")
cat("\nConfusion Matrix:\n")
print(xgb_cm)

# Save results
results <- data.frame(
  Name = metadata_test$Name,
  DTXSID = metadata_test$DTXSID,
  CASRN = metadata_test$CASRN,
  SMILES = metadata_test$SMILES,
  MONOISOTOPIC.MASS = metadata_test$MONOISOTOPIC.MASS,
  Actual = y_test,
  Predicted = xgb_pred,
  Probability = xgb_pred_prob
)

write.csv(results, "xgb_model_predictionsBT", row.names = FALSE)

# Plot ROC curve
plot(roc_obj, col = "blue", main = "ROC Curve for XGBoost Model")
