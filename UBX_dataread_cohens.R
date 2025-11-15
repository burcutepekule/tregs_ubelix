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

# Load files
results_merged             = c()
sterile_comparison_keep    = c()
pathogenic_comparison_keep = c()


path      = "/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/"
files_0_0_0   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_tregs_0_trnd_0.rds$", full.names = TRUE)
files_0_1_0   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_tregs_1_trnd_0.rds$", full.names = TRUE)
files_0_1_1   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_tregs_1_trnd_1.rds$", full.names = TRUE)
files_1_0_0   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_1_tregs_0_trnd_0.rds$", full.names = TRUE)
files_1_1_0   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_1_tregs_1_trnd_0.rds$", full.names = TRUE)
files_1_1_1   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_1_tregs_1_trnd_1.rds$", full.names = TRUE)


indices_0_0_0 = str_extract(basename(files_0_0_0), "\\d+") |> as.numeric()
indices_0_1_0 = str_extract(basename(files_0_1_0), "\\d+") |> as.numeric()
indices_0_1_1 = str_extract(basename(files_0_1_1), "\\d+") |> as.numeric()
indices_1_0_0 = str_extract(basename(files_1_0_0), "\\d+") |> as.numeric()
indices_1_1_0 = str_extract(basename(files_1_1_0), "\\d+") |> as.numeric()
indices_1_1_1 = str_extract(basename(files_1_1_1), "\\d+") |> as.numeric()

indices = Reduce(intersect, list(
  indices_0_0_0,
  indices_0_1_0,
  indices_0_1_1,
  indices_1_0_0,
  indices_1_1_0,
  indices_1_1_1
))

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/ids_read_cohens.rds')){
  inds_read = readRDS('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/ids_read_cohens.rds')
}else{
  inds_read = c()
}

inds2read = sort(setdiff(indices,inds_read))
length(inds2read)

if(length(inds2read)>0){
  for (i in inds2read){
    message("Processing param_set_", i)
    results_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_tregs_0_trnd_0.rds'))
    results_0_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_tregs_1_trnd_0.rds'))
    results_0_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_tregs_1_trnd_1.rds'))
    
    results_1_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_1_tregs_0_trnd_0.rds'))
    results_1_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_1_tregs_1_trnd_0.rds'))
    results_1_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_1_tregs_1_trnd_1.rds'))
    
    results = rbind(
      results_0_0_0,
      results_0_1_0,
      results_0_1_1,
      results_1_0_0,
      results_1_1_0,
      results_1_1_1
    )
    
    # Merge
    results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                          5*epithelial_inj_1+
                                          4*epithelial_inj_2+
                                          3*epithelial_inj_3+
                                          2*epithelial_inj_4+
                                          1*epithelial_inj_5)
    
    full_data_comparison = results %>% dplyr::select(param_set_id, sterile, tregs_on, randomize_tregs, rep_id, t, epithelial_score)
    min_reps  = min(full_data_comparison$rep_id)
    max_reps  = max(full_data_comparison$rep_id)
    t_max_ind = max(full_data_comparison$t)
    
    scores_0_keep = c()
    scores_1_keep = c()
    scores_2_keep = c()
    scores_3_keep = c()
    scores_4_keep = c()
    scores_5_keep = c()

    all_comparison_results_reps = data.frame()
    
    for (rep in min_reps:max_reps){  
      
      #### PATHOGENIC INJURY
      # tregs OFF
      full_data_comparison_scores_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & tregs_on==0 & randomize_tregs==0)
      # tregs ON
      full_data_comparison_scores_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & tregs_on==1 & randomize_tregs==0)
      # tregs ON, BUT ARE random
      full_data_comparison_scores_2 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & tregs_on==1 & randomize_tregs==1)
      
      #### STERILE INJURY
      # tregs OFF
      full_data_comparison_scores_3 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on==0 & randomize_tregs==0)
      # tregs ON
      full_data_comparison_scores_4 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on==1 & randomize_tregs==0)
      # tregs ON, BUT ARE random
      full_data_comparison_scores_5 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on==1 & randomize_tregs==1)
      
      # --- Steady-state detection ---
      time_ss_0 = steady_state_idx(full_data_comparison_scores_0$epithelial_score)
      time_ss_1 = steady_state_idx(full_data_comparison_scores_1$epithelial_score)
      time_ss_2 = steady_state_idx(full_data_comparison_scores_2$epithelial_score)
      time_ss_3 = steady_state_idx(full_data_comparison_scores_3$epithelial_score)
      time_ss_4 = steady_state_idx(full_data_comparison_scores_4$epithelial_score)
      time_ss_5 = steady_state_idx(full_data_comparison_scores_5$epithelial_score)
      
      time_ss_vec = c(time_ss_0, time_ss_1, time_ss_2, time_ss_3, time_ss_4, time_ss_5)
      
      if(!any(is.na(time_ss_vec))){
        
        scores_0 = full_data_comparison_scores_0$epithelial_score[time_ss_0:t_max_ind]
        scores_1 = full_data_comparison_scores_1$epithelial_score[time_ss_1:t_max_ind]
        scores_2 = full_data_comparison_scores_2$epithelial_score[time_ss_2:t_max_ind]
        scores_3 = full_data_comparison_scores_0$epithelial_score[time_ss_3:t_max_ind]
        scores_4 = full_data_comparison_scores_1$epithelial_score[time_ss_4:t_max_ind]
        scores_5 = full_data_comparison_scores_2$epithelial_score[time_ss_5:t_max_ind]
        
        scores_0_keep = c(scores_0_keep, scores_0)
        scores_1_keep = c(scores_1_keep, scores_1)
        scores_2_keep = c(scores_2_keep, scores_2)
        scores_3_keep = c(scores_3_keep, scores_3)
        scores_4_keep = c(scores_4_keep, scores_4)
        scores_5_keep = c(scores_5_keep, scores_5)
        
        # --- Tabulate all comparisons ---
        comparison_results = data.frame(
          param_set_id = i,
          replicate_id = rep,
          injury_type  = c("pathogenic","pathogenic","pathogenic","sterile","sterile","sterile"),
          tregs_on     = c(0, 1, 1, 0, 1, 1),
          tregs_rnd    = c(0, 0, 1, 0, 0, 1),
          ss_start     = c(time_ss_0, time_ss_1, time_ss_3, time_ss_3, time_ss_4, time_ss_5),
          mean_score   = c(mean(scores_0), mean(scores_1), mean(scores_2), mean(scores_3), mean(scores_4), mean(scores_5))
        )
        
        # Append to global results
        all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
      }
    }
    
    all_comparison_results_reps$d_10 = cohens_d(scores_1_keep, scores_0_keep)
    all_comparison_results_reps$d_21 = cohens_d(scores_2_keep, scores_1_keep)
    all_comparison_results_reps$d_43 = cohens_d(scores_4_keep, scores_3_keep)
    all_comparison_results_reps$d_54 = cohens_d(scores_5_keep, scores_4_keep) 
    
    all_comparison_results = bind_rows(all_comparison_results, all_comparison_results_reps)
    
  }
  message("num of param_id successfully added: ", length(inds2read))
  
  #---- read previous file
  if(!file.exists('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/ids_read_cohens.rds')){
    saveRDS(inds2read, '/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/ids_read_cohens.rds')
  }else{
    inds2read_old = readRDS('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/ids_read_cohens.rds')
    inds2read     = c(inds2read_old, inds2read)
    saveRDS(inds2read, '/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/ids_read_cohens.rds')
  }
  
  #---- read previous file
  if(!file.exists('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/data_cpp_read_cohens.rds')){
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/data_cpp_read_cohens.rds')
  }else{
    all_comparison_results_old = readRDS('/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/data_cpp_read_cohens.rds')
    all_comparison_results     = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/data_cpp_read_cohens.rds')
  }
}else{
  message("No new pts added.")
}



