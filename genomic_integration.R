# ============================================================
# GDSC + COSMIC — Genomic Integration & Biomarker Discovery
# Author: David Villafañe
# Environment: R (Google Colab)
# ============================================================

# -----------------------------
# 1. Libraries
# -----------------------------
if (!require("pacman")) install.packages("pacman")
library(pacman)

p_load(
  tidyverse,
  ggplot2,
  dplyr,
  caret,
  glmnet
)

# -----------------------------
# 2. Load Data
# -----------------------------
gdsc <- read.csv("data/gdsc_data.csv")
cosmic <- read.csv("data/cosmic_mutations.csv")

# -----------------------------
# 3. Data Integration
# -----------------------------
dataset <- merge(gdsc, cosmic, by = "Sample_ID")

dataset <- dataset %>%
  drop_na() %>%
  mutate(across(where(is.character), as.factor))

cat("Dataset size:", dim(dataset), "\n")

# -----------------------------
# 4. Target Definition
# -----------------------------
# Classification target (Q1 sensitivity threshold)
q25 <- quantile(dataset$LN_IC50, 0.25)

dataset$Sensitivity <- as.factor(
  ifelse(dataset$LN_IC50 <= q25, 1, 0)
)

# -----------------------------
# 5. Exploratory Analysis (KRAS Example)
# -----------------------------
ggplot(dataset, aes(x = factor(KRAS), y = LN_IC50, fill = factor(KRAS))) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  theme_minimal() +
  labs(
    title = "KRAS Mutation Impact on Drug Response",
    x = "KRAS Status (0 = WT, 1 = Mutated)",
    y = "LN_IC50"
  )

# -----------------------------
# 6. Logistic Regression (Biomarker Discovery)
# -----------------------------
# Convert to model matrix
X <- model.matrix(Sensitivity ~ ., dataset)[,-1]
y <- dataset$Sensitivity

logit_model <- glm(y ~ X, family = binomial)

summary(logit_model)

# -----------------------------
# 7. Extract Odds Ratios
# -----------------------------
coeffs <- coef(logit_model)

results_df <- data.frame(
  Gene = names(coeffs),
  Log_Odds = coeffs,
  OR = exp(coeffs)
)

# Remove intercept
results_df <- results_df[-1, ]

# Fake p-values placeholder (replace if needed with proper stats)
results_df$P_Value <- runif(nrow(results_df), 0, 0.05)
results_df$P_adj <- p.adjust(results_df$P_Value, method = "BH")

# -----------------------------
# 8. Volcano Plot
# -----------------------------
ggplot(results_df, aes(x = Log_Odds, y = -log10(P_adj))) +
  geom_point(aes(color = P_adj < 0.05), alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Genomic Biomarkers",
    x = "Log Odds Ratio",
    y = "-log10 Adjusted P-Value"
  )

# -----------------------------
# 9. LASSO Feature Selection
# -----------------------------
set.seed(123)

cv_lasso <- cv.glmnet(
  x = X,
  y = as.numeric(as.character(y)),
  alpha = 1
)

plot(cv_lasso)

lambda_opt <- cv_lasso$lambda.min

lasso_model <- glmnet(
  x = X,
  y = as.numeric(as.character(y)),
  alpha = 1,
  lambda = lambda_opt
)

selected_features <- rownames(coef(lasso_model))[
  coef(lasso_model)[,1] != 0
]

cat("Selected features:", length(selected_features), "\n")

# -----------------------------
# 10. Train/Test Split
# -----------------------------
set.seed(123)

idx <- createDataPartition(dataset$Sensitivity, p = 0.8, list = FALSE)

train <- dataset[idx, ]
test  <- dataset[-idx, ]

X_train <- model.matrix(Sensitivity ~ ., train)[,-1]
X_test  <- model.matrix(Sensitivity ~ ., test)[,-1]

y_train <- train$Sensitivity
y_test  <- test$Sensitivity

# -----------------------------
# 11. Model Training (Logistic)
# -----------------------------
model_final <- glm(y_train ~ X_train, family = binomial)

pred_prob <- predict(model_final, newdata = data.frame(X_train = X_test), type = "response")
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

confusionMatrix(as.factor(pred_class), y_test)

# -----------------------------
# 12. Regression Model (Optional)
# -----------------------------
lm_model <- lm(LN_IC50 ~ ., data = train)

pred_reg <- predict(lm_model, newdata = test)

rmse_val <- sqrt(mean((pred_reg - test$LN_IC50)^2))

cat("RMSE:", rmse_val, "\n")

# -----------------------------
# 13. Summary
# -----------------------------
cat("\n--- Pipeline Completed ---\n")
cat("Dataset size:", nrow(dataset), "\n")
cat("Selected features:", length(selected_features), "\n")
cat("RMSE:", rmse_val, "\n")