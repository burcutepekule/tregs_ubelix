rm(list=ls())
sd_tol_in   = Inf 
jsd_th      = 0.3
tol_in      = 125*0.2
M1_M2_diff  = 0

source('~/Dropbox/tregs_ubelix/UBX_datazanalyse_jsd_regions.R')

df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
if(M1_M2_diff==1){
  df_plot = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
}
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R') #M1_M2_diff adjusts params as well

inj_type= 'sterile'
inj_type= 'pathogenic'
inj_type= 'pooled'

if(inj_type!='pooled'){
  df_plot = df_plot %>% dplyr::filter(injury_type==inj_type)
}
df_lda  = df_plot %>% dplyr::select(all_of(param_names), tregs_better_cohens)
classes = unique(df_lda$tregs_better_cohens)

# Choose confidence levels for visualization
level_plus1   = 0.95 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
level_minus1  = 0.95 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
sigma_coeff   = 1.5
n_new_samples = 10000
n_arrows      = 8 # 24 for all params
source('~/Dropbox/tregs_ubelix/sub_UBX_datazanalyse_PLS_DA_3classes_only_worse.R')

