library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)

setwd('/Users/burcutepekule/Dropbox/tregs_ubelix')
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/tregs_better_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_jsd_tregs_better.rds')
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
keep_param_id = reps_df %>% dplyr::filter(Freq==600) %>% dplyr::pull(Var1) # 600 = 100 reps per scenario, 6 scenarios
df_results    = df_results_keep %>% filter(param_set_id %in% keep_param_id)
unique(table(df_results$param_set_id)) #66, sanity check
length(unique(df_results$param_set_id))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
length(param_id_all_below)/length(unique(df_results$param_set_id)) # >99%!
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
max(df_results$ss_start)<ss_start_threshold # TRUE, sanity check
unique(table(df_results$param_set_id)) #600, sanity check
length(unique(df_results$param_set_id))

df_comparisons_mean = df_results %>%
  # Create a key combining the conditions
  mutate(condition = case_when(
    tregs_on == 0 & tregs_rnd == 0 ~ "tregs_off_mean",
    tregs_on == 1 & tregs_rnd == 0 ~ "tregs_on_no_rnd_mean",
    tregs_on == 1 & tregs_rnd == 1 ~ "tregs_on_rnd_mean"
  )) 

df_comparisons = distinct(df_comparisons_mean %>% dplyr::select(param_set_id, injury_type,
                                                                mean_score_all_ts_0, mean_score_all_ts_1, 
                                                                mean_score_all_ts_2, mean_score_all_ts_3,
                                                                mean_score_all_ts_4, mean_score_all_ts_5,
                                                                sd_score_all_ts_0, sd_score_all_ts_1, sd_score_all_ts_2,
                                                                sd_score_all_ts_3, sd_score_all_ts_4, sd_score_all_ts_5))

df_comparisons = df_comparisons %>% dplyr::mutate(diff_tregs_on_minus_off = ifelse(injury_type=='pathogenic', 
                                                                                   mean_score_all_ts_1-mean_score_all_ts_0,
                                                                                   mean_score_all_ts_4-mean_score_all_ts_3))

df_comparisons = df_comparisons %>% dplyr::mutate(diff_notrnd_minus_rnd = ifelse(injury_type=='pathogenic', 
                                                                                 mean_score_all_ts_1-mean_score_all_ts_2,
                                                                                 mean_score_all_ts_4-mean_score_all_ts_5))

df_comparisons = df_comparisons %>% dplyr::mutate(tregs_off_sd = ifelse(injury_type=='pathogenic', sd_score_all_ts_0, sd_score_all_ts_3))
df_comparisons = df_comparisons %>% dplyr::mutate(tregs_on_no_rnd_sd = ifelse(injury_type=='pathogenic', sd_score_all_ts_1, sd_score_all_ts_4))
df_comparisons = df_comparisons %>% dplyr::mutate(tregs_on_rnd_sd = ifelse(injury_type=='pathogenic', sd_score_all_ts_2, sd_score_all_ts_5))

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, diff_tregs_on_minus_off, diff_notrnd_minus_rnd, 
                                                  tregs_off_sd, tregs_on_no_rnd_sd, tregs_on_rnd_sd)
df_comparisons = merge(df_comparisons, distinct(df_results[c('param_set_id','d_10','d_21','d_43','d_54')]), by='param_set_id')

#----- filter based on cohen's d since you are looking for differences!
df_comparisons_below= df_comparisons %>% dplyr::filter((injury_type=='pathogenic' & abs(d_10)<jsd_th)|(injury_type=='sterile' & abs(d_43)<jsd_th))
# tol_in              = 125*0.15 # this is NOT over epithelial cells (so don't use 25*stg), it's over epithelial score! max diff is 6*25-1*25=125!
# or get tol_in based on the range that is not significant!
if(is.na(tol_in)){
  tol_in = max(max(df_comparisons_below$diff_tregs_on_minus_off), abs(min(df_comparisons_below$diff_tregs_on_minus_off)))
}

# # === CUTOFF FOR VARIANCE? STABLE PARAMETER SETS?
# log10_variance = log10(c(df_comparisons$tregs_off_sd, df_comparisons$tregs_on_no_rnd_sd, df_comparisons$tregs_on_rnd_sd))
# d = density(log10_variance)
# y = d$y
# x = d$x
# # Local maxima (peaks)
# peaks_idx = which(diff(sign(diff(y))) == -2) + 1
# peaks_x   = x[peaks_idx]
# peaks_y   = y[peaks_idx]
# # Local minima (valleys)
# valleys_idx = which(diff(sign(diff(y))) == 2) + 1
# valleys_x   = x[valleys_idx]
# valleys_y   = y[valleys_idx]
# 
# plot(x, y, type="l")
# points(peaks_x, peaks_y, col="red", pch=19, cex=1.5)
# points(valleys_x, valleys_y, col="blue", pch=19, cex=1.5)
# 
# sd_tol_in = 10^valleys_x[2] # take the second, less conservative?

# Your variable (epithelial score) spans 100 units: between 25 to 125 
# For a distribution on that scale:
# SD ≈ 5–10 → very low variation
# SD ≈ 10–20 → moderate variation
# SD ≈ 20–30+ → high variation

#----- cases where randomizing tregs have the opposite effect to turning them on, and significant?

df_comparisons_keep  = df_comparisons
df_comparisons_sense = df_comparisons %>% dplyr::mutate(makes_sense=ifelse(sign(diff_tregs_on_minus_off)!=sign(diff_notrnd_minus_rnd) 
                                                                           & abs(diff_notrnd_minus_rnd)>tol_in
                                                                           & (abs(d_21)>=jsd_th | abs(d_54)>=jsd_th) ,0,1)) #either for sterile or pathogenic
df_comparisons_sense = df_comparisons_sense %>% dplyr::filter(makes_sense==0) 
dim(df_comparisons_sense)[1] #0! -> ALL MAKES SENSE # IF (sign(diff_tregs_on_minus_off)!=sign(diff_notrnd_minus_rnd), it's insignificant, meaning d_21<0.5 (or d_54<0.5)!

#----- Find significant cases

df_comparisons_plot            = df_comparisons_keep
# --- PATHOGENIC ---
df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(tregs_better = ifelse(diff_tregs_on_minus_off > tol_in, 1,
                                                                                                        ifelse(diff_tregs_on_minus_off < -1*tol_in,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(tregs_better_cohens = ifelse(abs(d_10)>jsd_th,
                                                                                                               # & tregs_on_no_rnd_sd<sd_tol_in
                                                                                                               # & tregs_off_sd<sd_tol_in,
                                                                                                               tregs_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-d_43)
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)=='d_10')]='cohens_d'

# --- STERILE ---
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(tregs_better = ifelse(diff_tregs_on_minus_off > tol_in, 1,
                                                                                                  ifelse(diff_tregs_on_minus_off < -1*tol_in,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(tregs_better_cohens = ifelse(abs(d_43)>jsd_th,
                                                                                                         # & tregs_on_no_rnd_sd<sd_tol_in
                                                                                                         # & tregs_off_sd<sd_tol_in,
                                                                                                         tregs_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-d_10)
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)=='d_43')]='cohens_d'

# ============CONFLICT?=========================================================
tregs_better_sterile = df_comparisons_plot_sterile %>% dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_sterile  = df_comparisons_plot_sterile %>% dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)
tregs_better_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_pathogenic  = df_comparisons_plot_pathogenic %>% dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)

s_better_p_better = intersect(tregs_better_sterile, tregs_better_pathogenic)
s_better_p_worse  = intersect(tregs_better_sterile, tregs_worse_pathogenic)
s_worse_p_better  = intersect(tregs_worse_sterile, tregs_better_pathogenic)
s_worse_p_worse   = intersect(tregs_worse_sterile, tregs_worse_pathogenic)
# ==============================================================================

df_comparisons_plot = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)
saveRDS(df_comparisons_plot,'./df_comparisons_plot.rds')

cohens_th = jsd_th # example threshold for x
e_th      = tol_in # example threshold for y

dfp    = df_comparisons_plot %>%
  mutate(
    x = abs(cohens_d),
    y = diff_tregs_on_minus_off
  )

dfp = dfp %>%
  mutate(color_group = case_when(
    x >= cohens_th & y >= e_th & tregs_better_cohens!=0 ~ as.character(param_set_id),  # blue region
    x >= cohens_th & y <= -e_th & tregs_better_cohens!=0 ~ as.character(param_set_id), # pink region
    TRUE ~ "other"  # everything else
  ))

# Create separate color assignments for each region
dfp = dfp %>%
  mutate(region = case_when(
    x >= cohens_th & y >= e_th & tregs_better_cohens!=0 ~ "blue_region",  # blue region
    x >= cohens_th & y <= -e_th & tregs_better_cohens!=0 ~ "pink_region", # pink region
    TRUE ~ "other"
  ))

# Get unique param_set_ids in each region
blue_params = unique(dfp$param_set_id[dfp$region == "blue_region"])
pink_params = unique(dfp$param_set_id[dfp$region == "pink_region"])

# Create color palettes
library(colorspace)

blue_colors = setNames(
  colorRampPalette(c("#1A2A80", "#3B38A0", "#637AB9", "#4FB7B3", "#016B61","#2F5755","#313647"))(length(blue_params)), 
  as.character(blue_params)
)

pink_colors = setNames(
  colorRampPalette(c("#FF9B17","#FA812F","#E83F25","#F75270", "#DC143C","#BF092F","#7D0A0A"))(length(pink_params)), 
  as.character(pink_params)
)

# Combine all colors
all_colors = c("other" = "gray70", blue_colors, pink_colors)

# Update color_group to use param_set_id for coloring
dfp = dfp %>%
  mutate(color_group = ifelse(region == "other", "other", as.character(param_set_id)))

# Before creating the plot, use the order from all_colors
dfp$color_group = factor(dfp$color_group, levels = names(all_colors))

p = ggplot(dfp, aes(x = x, y = y, shape = injury_type, color = color_group)) +
  annotate("rect",
           xmin = cohens_th, xmax = Inf,
           ymin = -Inf, ymax = -1*e_th,
           fill = "pink", alpha = 0.3) +
  annotate("rect",
           xmin = cohens_th, xmax = Inf,
           ymin = e_th, ymax = Inf,
           fill = "lightblue", alpha = 0.3) +
  annotate("rect",
           xmin = -Inf, xmax = cohens_th,
           ymin = -Inf, ymax = Inf,
           fill = "gray80", alpha = 0.15) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -1*e_th, ymax = e_th,
           fill = "gray80", alpha = 0.15) +
  geom_point(size = 3) +
  geom_text_repel(data = dfp %>% filter(region != "other"),
                  aes(label = param_set_id),
                  size = 3,
                  max.overlaps = 20,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.2) +
  scale_shape_manual(name = "Injury Type", values = c("sterile" = 16, "pathogenic" = 2)) +
  scale_color_manual(name = "Parameter Set ID", values = all_colors) +
  geom_vline(xintercept = cohens_th, linetype = "dashed") +
  geom_hline(yintercept = e_th, linetype = "dashed") +
  geom_hline(yintercept = -1*e_th, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid", col='red') +
  theme_minimal() +
  labs(x = paste0("jensen-shannon distance (threshold = ", cohens_th, ")"),
       y = paste0("Epithelial score for tregs_on - tregs_off (threshold = ±", round(e_th), " for max score of 125.)"),
       title = paste0('Num of parameter sets: ',length(unique(dfp$param_set_id))),
       shape = "Injury Type",
       color = "Parameter Set ID") + 
  scale_y_continuous(
    breaks = sort(c(seq(round(min(dfp$diff_tregs_on_minus_off))-5,
                        round(max(dfp$diff_tregs_on_minus_off))+5, by=10), 0))
  ) +
  scale_x_continuous(
    breaks = sort(c(seq(0, 1, by=0.1), cohens_th))
  ) +
  guides( # 3 columns for legends  # <-- REMOVED the + here
    color = guide_legend(ncol = 3),
    shape = guide_legend(ncol = 3)
  )

ggsave(
  filename = paste0("/Users/burcutepekule/Dropbox/tregs_ubelix/treg_better_JSD_tregs_better.png"),
  plot = p,
  width = 9,
  height = 6,
  dpi = 300,
  bg='white'
)

# print(p)

print(length(unique(dfp$param_set_id)))
dim(df_comparisons_sense)[1] #0! -> ALL MAKES SENSE # IF (sign(diff_tregs_on_minus_off)!=sign(diff_notrnd_minus_rnd), it's insignificant, meaning d_21<0.5 (or d_54<0.5)!

dfp_p = dfp %>% dplyr::filter(injury_type=='pathogenic') %>% dplyr::select(param_set_id, tregs_better_cohens)
dfp_s = dfp %>% dplyr::filter(injury_type=='sterile') %>% dplyr::select(param_set_id, tregs_better_cohens)
colnames(dfp_p)[2]='tregs_better_cohens_p'
colnames(dfp_s)[2]='tregs_better_cohens_s'
dfp_sp = merge(dfp_p, dfp_s, by='param_set_id')
dfp_sp = dfp_sp %>% dplyr::mutate(treg_better_cohens_both = ifelse(tregs_better_cohens_p+tregs_better_cohens_s==2, 1, 0))

100*length(which(dfp_sp$tregs_better_cohens_p==1))/dim(dfp_sp)[1] # % better for pathogenic
100*length(which(dfp_sp$tregs_better_cohens_p==-1))/dim(dfp_sp)[1] # % worse for pathogenic

100*length(which(dfp_sp$tregs_better_cohens_s==1))/dim(dfp_sp)[1] # % better for sterile
100*length(which(dfp_sp$tregs_better_cohens_s==-1))/dim(dfp_sp)[1] # % worse for sterile

100*sum(dfp_sp$treg_better_cohens_both)/dim(dfp_sp)[1] # % better for both

## check sds?

df_comparisons_plot_nz = df_comparisons_plot %>% dplyr::filter(tregs_better_cohens!=0)

# check conflict?
print(c(s_better_p_worse, s_worse_p_better))
