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

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
df_plot             = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R')
dir_name = '/Users/burcutepekule/Dropbox/tregs_ubelix/pca_plots'

inj_type_select = 'pathogenic'
df_pca          = df_plot %>% filter(injury_type == "pathogenic")
params_data     = df_pca[, param_names]
pca_result      = prcomp(params_data, scale. = TRUE, center = TRUE)

pca_df = data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  tregs_better_cohens = df_pca$tregs_better_cohens,
  param_set_id = df_pca$param_set_id
)

var_explained = summary(pca_result)$importance[2, 1:2] * 100

# Create the plot
p = ggplot(pca_df, aes(x = PC1, y = PC2, color = tregs_better_cohens)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis(option = "plasma", name = "Tregs Better\n(Cohen's d)") +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
    title = "PCA of Parameter Space",
    subtitle = "Colored by Treg Effect Size (Cohen's d)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "right"
  )

print(p)

# UMAP (better for non-linear relationships)
library(umap)
umap_result = umap(params_data)
umap_df = data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  tregs_better_cohens = df_pca$tregs_better_cohens
)

umap_df = umap_df %>% dplyr::filter(tregs_better_cohens!=0)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = tregs_better_cohens)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis(option = "plasma") + theme_minimal()
