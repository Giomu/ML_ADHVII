################### Second Analysis: Unaware Infection Identification ####################
# This script compares different classification models (k-NN, Random Forest, SVM-RBF) 
# for predicting asymptomatic infections using a labelled dataset (data2). 
# The models are trained with 5-fold cross-validation and their performance is evaluated  
# based on multiple metrics. A majority based vote consensus approach is then used for final  
# predictions on the unlabelled dataset (data3).
#
# Steps:
# 1. Load and explore synthetic data
# 2. Define cross-validation settings
# 3. Train k-NN, Random Forest, and SVM models
# 4. Compare model performance using several metrics
# 5. Predict asymptomatic infection status using models consensus
#
##########################################################################################

# Load necessary libraries
library(caret)
library(dplyr)
library(foreach)
library(doParallel)
library(MLmetrics)
library(stringr)
library(purrr)

# Load the labeled dataset (data2) for Model Construction phase
load("data/data2.rda")
# Ensure the target variable is a factor
data2$Class <- as.factor(data2$Class)
# Use ID as row names
rownames(data2) <- data2$ID
data2$ID <- NULL
# Transform the data: log2 transformation and scaling
data2[, c(1:8, 13)] <- log2(data2[, c(1:8, 13)] + 1)
# scale the data
data2[, c(1:13)] <- scale(data2[, c(1:13)], center = T, scale = T)

# Display basic dataset structure and summary
str(data2)
summary(data2)

##############################
# Model Construction and Evaluation
##############################

# Set a global seed
set.seed(1939) 

# Create 5 cross-validation folds to use for all models
cv_folds <- caret::createFolds(data2$Class, k = 5, returnTrain = TRUE)

# Define training control settings
control <- caret::trainControl(
  method = "cv", number = 5, 
  classProbs = TRUE, 
  summaryFunction = multiClassSummary, 
  verboseIter = TRUE, 
  index = cv_folds
)

# Define a function for hyper parameter tuning and training models
train_model <- function(method, tL) {
  caret::train(Class ~ ., data = data2, method = method, 
               trControl = control, tuneLength = tL)
}

# Train models with hyperparameter tuning
knn_model <- train_model("knn", tL = 5)
rf_model  <- train_model("rf", tL = 5)
svm_model <- train_model("svmRadial", tL = 5)

# Collect and compare model results
results <- resamples(list(kNN = knn_model, RF = rf_model, SVM = svm_model))

# Select metrics of interest
selected_metrics <- results$values %>% 
  select(contains(c("Accuracy", "Precision", "Recall", "F1")))

# Display performance summary
summary(selected_metrics)

# Visualization of model performance
bwplot(results)
dotplot(results)

##############################
# Variable Importance analysis
##############################
# Define a function to compute permutation-based feature importance for SVM and k-NN
permute_importance <- function(model, data, target_col, metric = "Accuracy", 
                               n_permutations = 10, parallel = TRUE) {
  set.seed(0306)
  
  y <- data[[target_col]]
  if (!is.factor(y)) y <- as.factor(y)  # Converte il target in fattore se necessario
  X <- data[, colnames(data) != target_col, drop = FALSE]
  
  # Check input data
  stopifnot(is.data.frame(data))
  stopifnot(target_col %in% colnames(data))
  stopifnot(nrow(data) > 10) 
  
  # Compute original accuracy
  original_preds <- predict(model, newdata = X)
  if (is.numeric(original_preds)) {
    original_preds <- ifelse(original_preds > 0.5, levels(y)[2], levels(y)[1])}
  original_acc <- mean(original_preds == y)
  
  # Prepare data frame to store importances
  importances <- data.frame(Feature = colnames(X), Importance = 0)
  
  # Allow parallel computation
  if (parallel) {
    registerDoParallel(cores = detectCores() - 1)}
  
  # Loop on each feature and compute importances
  results <- foreach(feature = colnames(X), .combine = rbind, .packages = "caret") %dopar% {
    acc_drops <- numeric(n_permutations)
    
    for (i in 1:n_permutations) {
      X_permuted <- X
      X_permuted[[feature]] <- sample(na.omit(X_permuted[[feature]]), replace = TRUE)
      
      permuted_preds <- predict(model, newdata = X_permuted)
      permuted_acc <- mean(permuted_preds == y)
      
      acc_drops[i] <- original_acc - permuted_acc}
    data.frame(Feature = feature, Importance = median(acc_drops))}
  
  # End parallel computation
  if (parallel) {
    stopImplicitCluster()}
  
  # Organize results
  results <- results[order(-results$Importance), ]
  
  # Plot importances
  p <- ggplot(results, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Feature Importance via Permutation",
         x = "Feature",
         y = "Importance (Drop in Accuracy)") +
    theme_minimal()
  print(p)
  
  return(results)
}

# Compute permutation-based importances for SVM-Radial
importance_svm <- permute_importance(svm_model, data2, "Class", metric = "Accuracy", n_permutations = 10)
print(importance_svm)
# Compute permutation-based importances for k-NN
importance_knn <- permute_importance(knn_model, data2, "Class", metric = "Accuracy", n_permutations = 30)
print(importance_knn)
# Compute vip based importances for Random Forest
importance_rf <- caret::varImp(rf_model, scale = TRUE)
importance_rf <- importance_rf$importance
importance_rf$Feature <- rownames(importance_rf)
ggplot(importance_rf, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance via Permutation",
       x = "Feature",
       y = "Importance (Drop in Accuracy)") +
  theme_minimal()
print(importance_rf)

##############################
# Model Application 
##############################
# Load unlabeled data (data3) for Model Application phase
load("data/data3.rda")
rownames(data3) <- data3$ID
data3$ID <- NULL
# save ID and self-declared class for later
self_df <- data.frame(self_status = data3$Infection_0_1pre_2post,
                      row.names = rownames(data3))
data3$Infection_0_1pre_2post <- NULL

data3[, c(1:8, 13)] <- log2(data3[, c(1:8, 13)] + 1)
# scale the data
data3[, c(1:13)] <- scale(data3[, c(1:13)], center = T, scale = T)

# Ensure feature consistency between data2 and data3
common_features <- intersect(names(data2), names(data3))
data3 <- data3[, common_features, drop = FALSE]

# Make predictions using trained models
self_df$kNN <- predict(knn_model, data3)
self_df$RF  <- predict(rf_model, data3)
self_df$SVM <- predict(svm_model, data3)

# Compare model predictions
table(self_df$kNN)
table(self_df$RF)
table(self_df$SVM)

# Clean up the predictions by removing "mc" prefix and find consensus
self_df <- self_df %>%
  mutate(
    kNN = str_remove(kNN, "mc"),
    RF = str_remove(RF, "mc"),
    SVM = str_remove(SVM, "mc"),
    consensus = pmap_chr(list(kNN, RF, SVM), function(a, b, c) {
      votes <- c(a, b, c)
      vote_table <- table(votes)
      winner <- names(vote_table)[which.max(vote_table)]
      if (max(vote_table) >= 2) winner else NA_character_
    })
  )

# Compare new predictions with self-declared status
self_df <- self_df %>%
  # convert self_status labels
  mutate(self_status = ifelse(self_status %in% c(1, 2), "I", "NI")) %>%
  
  # Assign new class based on self-declared status and consensus results
  mutate(new_class = case_when(
    self_status == "NI" & consensus == "NI" ~ "NI",
    self_status == "NI" & consensus == "I"  ~ "UI",
    self_status == "I"  & consensus == "I"  ~ "I",
    self_status == "I"  & consensus == "NI" ~ "exclude",
    TRUE ~ NA_character_  # for NAs or other cases
  ))

# Display the final classification results
table(self_df$new_class)
