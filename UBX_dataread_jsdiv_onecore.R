rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(philentropy)

setwd('/Users/burcutepekule/Dropbox/tregs_ubelix')
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

if(file.exists('/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')){
  inds_read = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
}else{
  inds_read = c()
}

inds2read = sort(setdiff(indices,inds_read))
length(inds2read)

colnames_insert_fixed = c("t","sterile","tregs_on","randomize_tregs","param_set_id","rep_id","epithelial_score","time_ss",
                          'epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4',
                          'epithelial_inj_5','phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3',
                          'phagocyte_M1_L_4','phagocyte_M1_L_5','phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2',
                          'phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5','commensal','pathogen','treg_resting',
                          'treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

if(length(inds2read)>0){
  # Track indices that have been successfully processed
  processed_indices = c()
  
  # for (i_idx in seq_along(inds2read)){
  for (i_idx in seq_along(inds2read)){
    
    i = inds2read[i_idx]
    message("Processing param_set_", i)
    # Check file sizes for this parameter set
    files_to_check = c(
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_tregs_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_tregs_1_trnd_1.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_1_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_1_tregs_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_1_tregs_1_trnd_1.rds')
    )
    if(any(file.info(files_to_check)$size<150000)){
      processed_indices      = c(processed_indices, i) #add and skip
      message("Skipped one")
    }else{
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
      colnames(results) = colnames_insert_fixed # this is fixing the mistake in UBX_datagen_cpp.R -> 7:37 should have been 9:39
      full_data_comparison = results %>% dplyr::select(param_set_id, sterile, tregs_on, randomize_tregs, rep_id, t, time_ss, epithelial_score)
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
        time_ss_0 = unique(full_data_comparison_scores_0$time_ss)
        time_ss_1 = unique(full_data_comparison_scores_1$time_ss)
        time_ss_2 = unique(full_data_comparison_scores_2$time_ss)
        time_ss_3 = unique(full_data_comparison_scores_3$time_ss)
        time_ss_4 = unique(full_data_comparison_scores_4$time_ss)
        time_ss_5 = unique(full_data_comparison_scores_5$time_ss)
        
        time_ss_vec = c(time_ss_0, time_ss_1, time_ss_2, time_ss_3, time_ss_4, time_ss_5)
        
        if(!any(is.na(time_ss_vec))){
          
          scores_0 = full_data_comparison_scores_0$epithelial_score[time_ss_0:t_max_ind]
          scores_1 = full_data_comparison_scores_1$epithelial_score[time_ss_1:t_max_ind]
          scores_2 = full_data_comparison_scores_2$epithelial_score[time_ss_2:t_max_ind]
          scores_3 = full_data_comparison_scores_3$epithelial_score[time_ss_3:t_max_ind]
          scores_4 = full_data_comparison_scores_4$epithelial_score[time_ss_4:t_max_ind]
          scores_5 = full_data_comparison_scores_5$epithelial_score[time_ss_5:t_max_ind]
          
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
            mean_score_ones_ts = c(mean(scores_0), mean(scores_1), mean(scores_2), mean(scores_3), mean(scores_4), mean(scores_5))
          )
          
          # Append to global results
          all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
        }
      }
      
      all_comparison_results_reps$d_10 = calculate_js_divergence(scores_1_keep, scores_0_keep, method = "fd")[1]
      all_comparison_results_reps$d_21 = calculate_js_divergence(scores_2_keep, scores_1_keep, method = "fd")[1]
      all_comparison_results_reps$d_43 = calculate_js_divergence(scores_4_keep, scores_3_keep, method = "fd")[1]
      all_comparison_results_reps$d_54 = calculate_js_divergence(scores_5_keep, scores_4_keep, method = "fd")[1]
      
      all_comparison_results_reps$n_10 = calculate_js_divergence(scores_1_keep, scores_0_keep, method = "fd")[2]
      all_comparison_results_reps$n_21 = calculate_js_divergence(scores_2_keep, scores_1_keep, method = "fd")[2]
      all_comparison_results_reps$n_43 = calculate_js_divergence(scores_4_keep, scores_3_keep, method = "fd")[2]
      all_comparison_results_reps$n_54 = calculate_js_divergence(scores_5_keep, scores_4_keep, method = "fd")[2]
      
      all_comparison_results_reps$mean_score_all_ts_0 = mean(scores_0_keep)
      all_comparison_results_reps$mean_score_all_ts_1 = mean(scores_1_keep)
      all_comparison_results_reps$mean_score_all_ts_2 = mean(scores_2_keep)
      all_comparison_results_reps$mean_score_all_ts_3 = mean(scores_3_keep)
      all_comparison_results_reps$mean_score_all_ts_4 = mean(scores_4_keep)
      all_comparison_results_reps$mean_score_all_ts_5 = mean(scores_5_keep)
      
      all_comparison_results_reps$sd_score_all_ts_0 = sd(scores_0_keep)
      all_comparison_results_reps$sd_score_all_ts_1 = sd(scores_1_keep)
      all_comparison_results_reps$sd_score_all_ts_2 = sd(scores_2_keep)
      all_comparison_results_reps$sd_score_all_ts_3 = sd(scores_3_keep)
      all_comparison_results_reps$sd_score_all_ts_4 = sd(scores_4_keep)
      all_comparison_results_reps$sd_score_all_ts_5 = sd(scores_5_keep)
      
      all_comparison_results = bind_rows(all_comparison_results, all_comparison_results_reps)
      processed_indices      = c(processed_indices, i)
    }
    # Save after every 10 parameter sets (if total is > 10)
    if(length(inds2read) > 10 && i_idx %% 10 == 0){
      message("Saving intermediate results after ", i_idx, " parameter sets...")
      
      # Update the list of read indices
      if(!file.exists('/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')){
        saveRDS(processed_indices, '/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
      }else{
        inds_read_old = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
        inds_read_updated = c(inds_read_old, processed_indices)
        saveRDS(inds_read_updated, '/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
      }
      
      # Update the results data
      if(!file.exists('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')){
        saveRDS(all_comparison_results, '/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
      }else{
        all_comparison_results_old = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
        all_comparison_results_combined = rbind(all_comparison_results_old, all_comparison_results)
        saveRDS(all_comparison_results_combined, '/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
      }
      
      # Reset for next batch
      all_comparison_results = data.frame()
      processed_indices = c()
      
      message("Intermediate save complete with 10 more param_ids.")
    }
  }
  
  # Final save for any remaining results
  message("Saving final results. Total param_sets processed: ", length(inds2read))
  
  if(!file.exists('/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')){
    saveRDS(processed_indices, '/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
  }else{
    inds_read_old = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
    inds_read_final = c(inds_read_old, processed_indices)
    saveRDS(inds_read_final, '/Users/burcutepekule/Dropbox/tregs_ubelix/ids_read_jsd.rds')
  }
  
  if(!file.exists('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')){
    saveRDS(all_comparison_results, '/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
  }else{
    all_comparison_results_old = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
    all_comparison_results_final = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results_final, '/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd.rds')
  }
  
  message("Reprinting space.")
  
  sd_tol_in   = Inf 
  jsd_th      = 0.3
  tol_in      = 125*0.2
  M1_M2_diff  = 1
  source('~/Dropbox/tregs_ubelix/UBX_datazanalyse_jsd_regions.R')
  
}else{
  message("No new pts added.")
}




