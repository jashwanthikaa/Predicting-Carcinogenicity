# Load libraries
library(readxl)
library(caret)
library(dplyr)
library(xgboost)
library(Matrix)
library(pROC)

# Load dataset
df <- read_excel("C:/Users/91962/LessData.xlsx")

# Check exact column names:
print(colnames(df))  # Check this to confirm actual column names

# Select fingerprint columns and mass
fp_cols <- grep("^Un", names(df), value = TRUE)
mass_col <- "MONOISOTOPIC.MASS"

# Prepare input and output
X <- df %>% select(all_of(c(fp_cols, mass_col)))
y <- df$HIT_CALL  # Assumed to be 0 or 1 numeric

# Metadata for saving predictions
metadata_test <- df %>% select(Name, DTXSID, CASRN, SMILES, MONOISOTOPIC.MASS)
# If 'Name' doesn't exist, replace with correct column name, e.g. 'PREFERRED.NAME' or similar

# Train/test split
set.seed(123)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
X_test <- X[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]
metadata_test <- metadata_test[-train_index, ]

# Convert to matrix (for xgboost and caret)
X_train_mat <- as.matrix(X_train)
X_test_mat <- as.matrix(X_test)

# Prepare factor for caret with Yes/No labels for training
y_train_factor <- factor(ifelse(y_train == 1, "Yes", "No"), levels = c("No", "Yes"))

# Train control with 5-fold CV and ROC metric
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Define tuning grid for xgboost
xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(4, 6, 8),
  eta = c(0.05, 0.1),
  gamma = 0,
  colsample_bytree = c(0.6, 0.8),
  min_child_weight = 1,
  subsample = 0.8
)

# Train the model with caret using xgbTree method
set.seed(123)
xgb_model <- train(
  x = X_train_mat,
  y = y_train_factor,
  method = "xgbTree",
  metric = "ROC",
  trControl = train_control,
  tuneGrid = xgb_grid
)

# Best hyperparameters
cat("\nBest hyperparameters:\n")
print(xgb_model$bestTune)

# Predict probabilities on test set
xgb_pred_prob <- predict(xgb_model, newdata = X_test_mat, type = "prob")[, "Yes"]

# Determine optimal threshold based on ROC
roc_obj <- roc(response = factor(ifelse(y_test == 1, "Yes", "No"), levels = c("No", "Yes")),
               predictor = xgb_pred_prob)

optimal_threshold <- coords(roc_obj, "best", ret = "threshold")
cat("\nOptimal threshold based on ROC:", optimal_threshold, "\n")

# Predict classes using optimal threshold and convert to Yes/No
xgb_pred <- ifelse(xgb_pred_prob > optimal_threshold, "Yes", "No")
actual_label <- ifelse(y_test == 1, "Yes", "No")

# Convert to factors with same levels for confusionMatrix
xgb_pred_factor <- factor(xgb_pred, levels = c("No", "Yes"))
actual_label_factor <- factor(actual_label, levels = c("No", "Yes"))

# Confusion matrix
cm <- confusionMatrix(xgb_pred_factor, actual_label_factor, positive = "Yes")
cat("\nConfusion Matrix:\n")
print(cm)

# Save predictions with metadata
results <- data.frame(
  Name = metadata_test$Name,
  DTXSID = metadata_test$DTXSID,
  CASRN = metadata_test$CASRN,
  SMILES = metadata_test$SMILES,
  MONOISOTOPIC_MASS = metadata_test$MONOISOTOPIC_MASS,
  Actual = actual_label_factor,
  Predicted = xgb_pred_factor,
  Probability = xgb_pred_prob
)

write.xlsx(results, "xgboost_predictions_yes_no.xlsx", rowNames = FALSE)
cat("\nPredictions saved to 'xgboost_predictions_yes_no.xlsx'\n")

# Plot ROC curve
plot(roc_obj, col = "blue", main = "ROC Curve for XGBoost Model")
