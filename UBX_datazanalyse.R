rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_results_keep = readRDS('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/data_cpp_read.rds')
length(unique(df_results_keep$param_set_id))
df_params  = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)

# ss_start_0: sterile==0 & tregs_on==0 & randomize_tregs==0
# ss_start_1: sterile==0 & tregs_on==1 & randomize_tregs==0
# ss_start_2: sterile==0 & tregs_on==1 & randomize_tregs==1
# ss_start_3: sterile==1 & tregs_on==0 & randomize_tregs==0
# ss_start_4: sterile==1 & tregs_on==1 & randomize_tregs==0
# ss_start_5: sterile==1 & tregs_on==1 & randomize_tregs==1

df_results = df_results_keep %>% dplyr::filter(comparison=='Treg_ON_OFF' & injury_type=='sterile')

df_plot    = df_results[c('param_set_id','replicate_id',
                   'ss_start_0','ss_start_1', # pick the relevant ones
                   'mean_treg_on','mean_treg_off')]

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
param_id_all_below = df_plot %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start_tregs_off <= ss_start_threshold & ss_start_tregs_on <= ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_plot = df_plot %>% dplyr::filter(param_set_id %in% param_id_all_below)


