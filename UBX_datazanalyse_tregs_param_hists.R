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

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
df_plot             = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
smoothing_param     = 0.01

inj_type_select     = c('sterile','pathogenic')
df_comparisons_plot = df_plot %>% dplyr::filter(injury_type %in% inj_type_select)
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R')

dir_name = '/Users/burcutepekule/Dropbox/tregs_ubelix/tregs_better_worse_histograms'
if(length(unique(df_comparisons_plot$tregs_better_cohens))==3){ # still no significant tregs worse case!
  df_comparisons_plot = df_comparisons_plot %>%
    dplyr::filter(tregs_better_cohens %in% c(-1, 0, 1)) %>%
    dplyr::mutate(tregs_better_cohens = factor(tregs_better_cohens, labels = c("Tregs worse", "Tregs DM", "Tregs better")))
  df_comparisons_plot = df_comparisons_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
  
  dir.create(dir_name, showWarnings = FALSE)
  for (param in param_names) {
    p = ggplot(df_comparisons_plot, aes(x = .data[[param]], fill = tregs_better_cohens)) +
      geom_density(alpha = 0.5, bw = smoothing_param) +
      scale_fill_manual(values = c("Tregs better" = "blue", "Tregs worse" = "red", "Tregs DM"='gray')) +
      labs(x = param, title = paste("Density of", param)) +
      theme_minimal()
    ggsave(
      filename = paste0(dir_name,"/",paste(inj_type_select, collapse = "_"),"_treg_better_",param,"_density.png"),
      plot = p,
      width = 9,
      height = 6,
      dpi = 300
    )
  }
}else{
  df_comparisons_plot = df_comparisons_plot %>%
    dplyr::filter(tregs_better_cohens %in% c(0, 1)) %>%
    dplyr::mutate(tregs_better_cohens = factor(tregs_better_cohens, labels = c("Tregs DM", "Tregs better")))
  df_comparisons_plot = df_comparisons_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
  
  dir.create(dir_name, showWarnings = FALSE)
  for (param in param_names) {
    p = ggplot(df_comparisons_plot, aes(x = .data[[param]], fill = tregs_better_cohens)) +
      geom_density(alpha = 0.5, bw = smoothing_param) +
      scale_fill_manual(values = c("Tregs better" = "blue", "Tregs DM"='gray')) +
      labs(x = param, title = paste("Density of", param)) +
      theme_minimal()
    ggsave(
      filename = paste0(dir_name,"/",paste(inj_type_select, collapse = "_"),"_treg_better_",param,"_density.png"),
      plot = p,
      width = 9,
      height = 6,
      dpi = 300
    )
  }
}
