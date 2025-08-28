# MODEL DEVELOPMENT PIPELINE
# Baselines: Logistic Regression, Decision Tree
# Advanced: XGBoost with extensive hyperparameter tuning

library(readxl)
library(caret)
library(dplyr)
library(rpart)
library(xgboost)
library(Matrix)

# -------------------------------------------------------------
# 1. Load data
# -------------------------------------------------------------
df <- read_excel("C:/Users/91962/All_319/Filtered_compounds_with_Set1.xlsx")

# Select fingerprint columns (starting with "Un") and MONOISOTOPIC.MASS
fp_cols <- grep("^Un", names(df), value = TRUE)
mass_col <- "MONOISOTOPIC.MASS"

# Input features: fingerprints + mass
X <- df %>%
  select(all_of(c(fp_cols, mass_col))) %>%
  as.data.frame()

# Output labels: convert "active"/"inactive" to binary 1/0
y <- ifelse(tolower(df$`HIT.CALL`) == "active", 1, 0)

# Remove rows with missing data
non_na_rows <- complete.cases(X)
X <- X[non_na_rows, ]
y <- y[non_na_rows]

df_clean <- df[non_na_rows, ]

# -------------------------------------------------------------
# 2. Train-test split
# -------------------------------------------------------------
set.seed(123)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[train_index, ]
y_train <- y[train_index]

X_test <- X[-train_index, ]
y_test <- y[-train_index]

# -------------------------------------------------------------
# 3. Baseline Model: Logistic Regression
# -------------------------------------------------------------
logistic_model <- glm(y_train ~ ., data = X_train, family = binomial)
logistic_pred_prob <- predict(logistic_model, newdata = X_test, type = "response")
logistic_pred <- ifelse(logistic_pred_prob > 0.5, 1, 0)

logistic_cm <- confusionMatrix(
  factor(logistic_pred),
  factor(y_test),
  positive = "1"
)
cat("\nLogistic Regression Performance:\n")
print(logistic_cm)

# -------------------------------------------------------------
# 4. Baseline Model: Decision Tree
# -------------------------------------------------------------
tree_model <- rpart(factor(y_train) ~ ., data = X_train, method = "class")
tree_pred <- predict(tree_model, newdata = X_test, type = "class")

tree_cm <- confusionMatrix(
  factor(tree_pred),
  factor(y_test),
  positive = "1"
)
cat("\nDecision Tree Performance:\n")
print(tree_cm)

# -------------------------------------------------------------
# 5. Advanced Model: XGBoost with expanded hyperparameter tuning
# -------------------------------------------------------------
# Convert to matrix format for XGBoost
X_train_matrix <- data.matrix(X_train)
X_test_matrix <- data.matrix(X_test)

# Convert y_train to factor for caret
y_train_factor <- factor(ifelse(y_train == 1, "Active", "Inactive"))

# Define training control with 5-fold CV
train_control <- trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Define expanded parameter grid
xgb_grid <- expand.grid(
  nrounds = c(50, 100, 200, 300, 500),
  max_depth = c(3, 5, 7, 9),
  eta = c(0.01, 0.05, 0.1, 0.2, 0.3),
  gamma = c(0, 0.1, 0.2),
  colsample_bytree = c(0.6, 0.8, 1.0),
  min_child_weight = c(1, 3, 5),
  subsample = c(0.6, 0.8, 1.0)
)

# Train XGBoost model with grid search
set.seed(123)
xgb_caret_model <- train(
  x = X_train_matrix,
  y = y_train_factor,
  trControl = train_control,
  tuneGrid = xgb_grid,
  method = "xgbTree",
  metric = "ROC"
)

# Print best model parameters
cat("\nBest XGBoost Parameters:\n")
print(xgb_caret_model$bestTune)

# -------------------------------------------------------------
# 6. Evaluate XGBoost on test set
# -------------------------------------------------------------
xgb_pred_prob <- predict(xgb_caret_model, newdata = X_test_matrix, type = "prob")[, "Active"]
xgb_pred <- ifelse(xgb_pred_prob > 0.5, 1, 0)

xgb_cm <- confusionMatrix(
  factor(xgb_pred),
  factor(y_test),
  positive = "1"
)
cat("\nXGBoost Performance:\n")
print(xgb_cm)

# -------------------------------------------------------------
# 7. Export XGBoost predictions with metadata
# -------------------------------------------------------------
metadata_cols <- c("Name", "CASRN", "SMILES", "MONOISOTOPIC.MASS", "DTXSID")
metadata_test <- df_clean[-train_index, metadata_cols]

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

write.csv(results, "model_comparison_predictions.csv", row.names = FALSE)
cat("\nPredictions saved to 'model_comparison_predictions.csv'.\n")
