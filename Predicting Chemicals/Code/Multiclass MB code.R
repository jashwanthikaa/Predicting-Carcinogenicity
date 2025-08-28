# Load libraries
library(readxl)
library(caret)
library(dplyr)
library(xgboost)
library(Matrix)
library(ggplot2)

# Load dataset
df <- read_excel("C:/Users/91962/Predicting Chemicals/FinalData.xlsx")

# Select fingerprint and mass columns
fp_cols <- grep("^Un", names(df), value = TRUE)
mass_col <- "MONOISOTOPIC.MASS"
X <- df %>% select(all_of(c(fp_cols, mass_col)))

# Target: multiclass HIT.CALL
y <- factor(tolower(df$`HIT.CALL`))  # Make sure it has 3+ unique levels

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

# Metadata for test results
metadata_cols <- c("Name", "CASRN", "SMILES", "MONOISOTOPIC.MASS", "DTXSID")
metadata_test <- df_clean[-train_index, metadata_cols]

# Define train control for multiclass (no ROC or classProbs)
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  verboseIter = TRUE
)

# Expanded tuning grid
xgb_grid <- expand.grid(
  nrounds = c(100, 200, 300),
  max_depth = c(6, 9, 12),
  eta = c(0.01, 0.05, 0.1),
  gamma = c(0, 1, 5),
  colsample_bytree = c(0.6, 0.8, 1),
  min_child_weight = c(1, 3, 5),
  subsample = c(0.6, 0.8, 1)
)

# Train XGBoost model
set.seed(123)
xgb_model <- train(
  x = X_train_matrix,
  y = y_train,  # Multiclass target (factor)
  trControl = train_control,
  tuneGrid = xgb_grid,
  method = "xgbTree",
  metric = "Accuracy",
  verbose = TRUE
)

# Best hyperparameters
cat("\nBest Hyperparameters:\n")
print(xgb_model$bestTune)

# Predict on test set
xgb_pred <- predict(xgb_model, newdata = X_test_matrix)

# Confusion matrix
xgb_cm <- confusionMatrix(xgb_pred, y_test)
cat("\nConfusion Matrix:\n")
print(xgb_cm)

# Save results
results <- data.frame(
  Name = metadata_test$Name,
  DTXSID = metadata_test$DTXSID,
  CASRN = metadata_test$CASRN,
  SMILES = metadata_test$SMILES,
  MONOISOTOPIC_MASS = metadata_test$MONOISOTOPIC.MASS,
  Actual = y_test,
  Predicted = xgb_pred
)

write.csv(results, "C:/Users/91962/Predicting Chemicals/xgb_multiclass_predictions.csv", row.names = FALSE)

# Plot confusion matrix
cm_df <- as.data.frame(xgb_cm$table)
ggplot(cm_df, aes(Prediction, Reference, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 4) +
  labs(title = "Confusion Matrix", x = "Predicted", y = "Actual") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal()
