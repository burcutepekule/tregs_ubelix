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
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R')

df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
df_plot             = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
df_plot_keep        = df_plot

# --- STERILE
inj_type_select = 'sterile'
df_plot_pick    = df_plot_keep %>% dplyr::filter(injury_type==inj_type_select)
df_tregs_better = df_plot_pick %>% dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
df_tregs_worse  = df_plot_pick %>% dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)

df_results_pick_better = df_results_keep %>% dplyr::filter(param_set_id %in% df_tregs_better & injury_type==inj_type_select)
df_results_pick_worse  = df_results_keep %>% dplyr::filter(param_set_id %in% df_tregs_worse & injury_type==inj_type_select)


df_results_pick_better = left_join(df_results_pick_better, df_plot_keep[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd')], 
                                    by=c('param_set_id','injury_type'))
df_results_pick_better_s = distinct(df_results_pick_better[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd',"d_43","d_54")]) #d_43 and d_54 for sterile

df_results_pick_worse = left_join(df_results_pick_worse, df_plot_keep[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd')], 
                                  by=c('param_set_id','injury_type'))
df_results_pick_worse_s = distinct(df_results_pick_worse[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd',"d_43","d_54")]) #d_43 and d_54 for sterile


# --- PATHOGENIC
inj_type_select = 'pathogenic'
df_plot_pick    = df_plot_keep %>% dplyr::filter(injury_type==inj_type_select)
df_tregs_better = df_plot_pick %>% dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
df_tregs_worse  = df_plot_pick %>% dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)

df_results_pick_better = df_results_keep %>% dplyr::filter(param_set_id %in% df_tregs_better & injury_type==inj_type_select)
df_results_pick_worse  = df_results_keep %>% dplyr::filter(param_set_id %in% df_tregs_worse & injury_type==inj_type_select)


df_results_pick_better = left_join(df_results_pick_better, df_plot_keep[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd')], 
                                   by=c('param_set_id','injury_type'))
df_results_pick_better_p = distinct(df_results_pick_better[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd',"d_10","d_21")]) #d_10 and d_21 for pathogenic

df_results_pick_worse = left_join(df_results_pick_worse, df_plot_keep[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd')], 
                                   by=c('param_set_id','injury_type'))
df_results_pick_worse_p = distinct(df_results_pick_worse[c('param_set_id','injury_type','diff_tregs_on_minus_off','diff_notrnd_minus_rnd',"d_10","d_21")]) #d_10 and d_21 for pathogenic


# --- MERGE
colnames(df_results_pick_better_s)[5:6] = c('cohens_d_on_off','cohens_d_nrnd_rnd')
colnames(df_results_pick_better_p)[5:6] = c('cohens_d_on_off','cohens_d_nrnd_rnd')
colnames(df_results_pick_worse_s)[5:6] = c('cohens_d_on_off','cohens_d_nrnd_rnd')
colnames(df_results_pick_worse_p)[5:6] = c('cohens_d_on_off','cohens_d_nrnd_rnd')

df_results_tregs_on_better = rbind(df_results_pick_better_s, df_results_pick_better_p)
df_results_tregs_on_worse  = rbind(df_results_pick_worse_s, df_results_pick_worse_p)















