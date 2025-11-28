rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)
library(viridis)
library(caret)
library(MASS)

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
df_plot             = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R')

# Linear Discriminant Analysis (LDA)
# Finds linear combinations that maximize class separation
df_lda = df_plot %>% dplyr::select(all_of(param_names), tregs_better_cohens)

lda_model = lda(tregs_better_cohens ~ ., data = df_lda)

# Project onto discriminant axes
lda_proj = predict(lda_model, df_lda)

lda_df = data.frame(
  LD1 = lda_proj$x[, 1],
  treg_outcome = df_lda$tregs_better_cohens
)

lda_df$treg_outcome = as.factor(lda_df$treg_outcome)
ggplot(lda_df, aes(x = LD1, y = 0, color = treg_outcome)) +
  geom_jitter(height = 0.2, alpha = 0.7, size = 3) +
  scale_color_manual(values = c("0" = "gray70", "1" = "blue")) +
  theme_minimal() +
  labs(title = "LDA Projection (1D)", 
       subtitle = "Optimized to separate classes")

# Check separation quality
cat("Between-class variance explained:", 
    round(lda_model$svd^2 / sum(lda_model$svd^2) * 100, 2), "%\n")

# Check the actual separation quality
lda_proj = predict(lda_model, df_lda)

# Look at the overlap between classes
lda_full_df <- data.frame(
  LD1 = lda_proj$x[, 1],
  treg_outcome = df_lda$tregs_better_cohens,
  predicted = lda_proj$class
)

# Density plot to see overlap
lda_full_df$treg_outcome = as.factor(lda_full_df$treg_outcome)
ggplot(lda_full_df, aes(x = LD1, fill = treg_outcome)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "gray70", "1" = "blue")) +
  theme_minimal() +
  labs(title = "LDA: Class Separation on Discriminant Axis",
       subtitle = "How much do the classes overlap?")

# Classification accuracy
confusion <- table(Actual = lda_full_df$treg_outcome, 
                   Predicted = lda_full_df$predicted)
print(confusion)
accuracy <- sum(diag(confusion)) / sum(confusion)
cat("LDA Classification Accuracy:", round(accuracy * 100, 1), "%\n")

# Specifically for class 1
sensitivity <- confusion["1", "1"] / sum(confusion["1", ])
cat("Sensitivity (finding 1s):", round(sensitivity * 100, 1), "%\n")
