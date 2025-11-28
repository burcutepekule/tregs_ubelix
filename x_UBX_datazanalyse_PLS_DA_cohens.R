rm(list=ls())
sd_tol_in   = Inf 
cohens_d_th = 0.5
tol_in      = NA
M1_M2_diff  = 0

source('~/Dropbox/tregs_ubelix/UBX_datazanalyse_cohens_SD_regions.R')
df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
if(M1_M2_diff==1){
  df_plot = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
}
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R') #M1_M2_diff adjusts params as well

inj_type= 'sterile'
inj_type= 'pathogenic'
# inj_type= 'pooled'

if(inj_type!='pooled'){
  df_plot = df_plot %>% dplyr::filter(injury_type==inj_type)
}
df_lda  = df_plot %>% dplyr::select(all_of(param_names), tregs_better_cohens)
classes = unique(df_lda$tregs_better_cohens)

# Choose confidence levels for visualization
level_plus1  = 0.5 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
level_minus1 = 0.5 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)

if(identical(classes, c(0,1))){
  source('~/Dropbox/tregs_ubelix/sub_UBX_datazanalyse_PLS_DA_2classes.R')
}else{
  source('~/Dropbox/tregs_ubelix/sub_UBX_datazanalyse_PLS_DA_3classes.R')
  source('~/Dropbox/tregs_ubelix/sub_UBX_datazanalyse_LDA_3classes.R')
}

## =====================================

df_params = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)

neg_ids = df_plot %>% dplyr::filter(tregs_better_cohens==-1) %>% pull(param_set_id)
pos_ids = df_plot %>% dplyr::filter(tregs_better_cohens==+1) %>% pull(param_set_id)
net_ids = df_plot %>% dplyr::filter(tregs_better_cohens==0) %>% pull(param_set_id)

df_params_pick = df_params %>% dplyr::filter(param_set_id %in% neg_ids)

summary(df_params_pick$epith_recovery_chance)
