rm(list=ls())
sd_tol_in   = Inf 
jsd_th      = 0.3
tol_in      = 125*0.2
# jsd_th      = 0.4
# tol_in      = 125*0.25
M1_M2_diff  = 1
exclude_param_id = c(204505) # exclude this for now, better for sterile, worse for pathogenic. Only one case!
filter_M1_better = 1

source('~/Dropbox/tregs_ubelix/UBX_datazanalyse_jsd_regions_merged_dataset.R')
source('~/Dropbox/tregs_ubelix/UBX_datazanalyse_jsd_regions_merged_dataset_nolabel.R')

table(df_comparisons_plot$tregs_better_cohens)

id_better = sort(df_comparisons_plot %>% dplyr::filter(tregs_better_cohens == 1) %>% pull(param_set_id))
id_better_both = data.frame(id = id_better) %>%
  count(id) %>%
  filter(n == 2) %>%
  pull(id)

df_params_merged_better_both = df_params_merged %>% dplyr::filter(param_set_id %in% id_better_both)

df_comparisons_plot_uniform = df_comparisons_plot %>% dplyr::filter(param_set_id <= max_param_id_uniform)
table(df_comparisons_plot_uniform$tregs_better_cohens)

df_comparisons_plot_targetted = df_comparisons_plot %>% dplyr::filter(param_set_id > max_param_id_uniform)
table(df_comparisons_plot_targetted$tregs_better_cohens)

#==== do I lose the benefit when spatial randomization is on?===================

df_tregs_better = df_comparisons_plot %>% dplyr::filter(tregs_better_cohens==1)
## when d_21 for sterile (or d_54 for pathogenic) approaches 1, that meant I lose the benefit of Tregs. 
## There is a few cases like that, so that means spatial heterogeneity is necessary (but not all the time)
## It would be interesting to see for which parameter sets spatial heterogeneity is NOT necessary

df_tregs_better_sterile = df_tregs_better %>% dplyr::filter(injury_type=='sterile')
plot(df_tregs_better_sterile$diff_tregs_on_minus_off, df_tregs_better_sterile$d_21)

df_tregs_pathogenic = df_tregs_better %>% dplyr::filter(injury_type=='pathogenic')
plot(df_tregs_pathogenic$diff_tregs_on_minus_off, df_tregs_pathogenic$d_54)

#==== do I lose the benefit when spatial randomization is on?===================

df_plot = merge(df_comparisons_plot, df_params_merged, by='param_set_id')

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
level_plus1  = 0.75 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
level_minus1 = 0.75 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)

n_arrows  = 12 #24 to select all 
source('~/Dropbox/tregs_ubelix/sub_UBX_datazanalyse_PLS_DA.R')

violin_on = 0
source('~/Dropbox/tregs_ubelix/sub_UBX_param_hists_category.R')


