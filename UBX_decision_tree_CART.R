rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(mgcv)
library(factoextra)
library(cluster)

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params           = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep     = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
df_plot             = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
df_plot             = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
source('~/Dropbox/tregs_ubelix/LOAD_PARAM_VECTOR.R')
# # ---------- Conditional Inference Trees (CTree): These use statistical tests for splitting and handle interactions better:
library(partykit)
inj_type      = 'pathogenic'
df_plot       = df_plot %>% dplyr::filter(injury_type==inj_type)     
df_clustering = distinct(df_plot[c('tregs_better_cohens',param_names)])
sort(unique(df_clustering$tregs_better_cohens))
df_clustering = df_clustering %>% mutate(treg_outcome = factor(tregs_better_cohens, levels = c(-1,0,1)))
table(df_clustering$tregs_better_cohens)

library(rpart)
library(rpart.plot)

# Calculate class weights (inverse of frequency)
class_weights = 1 / prop.table(table(df_clustering$treg_outcome))

cart_purity = rpart(
  treg_outcome ~ .,
  data = df_clustering %>% select(treg_outcome, all_of(param_names)),
  method = "class",
  parms = list(split = "gini"),  # Gini explicitly optimizes purity
  weights = ifelse(df_clustering$treg_outcome == "0", 
                   class_weights["0"],
                   class_weights["1"]),
  control = rpart.control(
    cp = 0.001,
    minsplit = 5,
    minbucket = 5,
    maxdepth = 5
  )
)

# Save a very high resolution plot to see all details
png(paste0("./tree_tregs_CART_",inj_type,".png"), width = 3000, height = 3000, res = 300)
rpart.plot(
  cart_purity,
  type = 4,
  extra = 104,
  under = TRUE,      # put labels underneath
  faclen = 0,
  cex = 0.6,
  tweak = 1.2,
  fallen.leaves = TRUE,
  nn = TRUE          # <— shows node numbers on the plot
)
dev.off()


# Get full frame info
frame_info = cart_purity$frame
# Add node IDs
frame_info$node_id = as.numeric(row.names(frame_info))

# note that rpart.plot shows weighted class proportions — NOT the raw ones! so the very right terminal
# leaf will have 50% of the data although it's all class 1, which is not the case 
terminal_nodes = frame_info %>%
  dplyr::filter(var == "<leaf>") %>%
  dplyr::mutate(
    Node_ID = as.numeric(row.names(.)),
    Predicted_Class = ifelse(yval == 1, "High Variance", "Low Variance"),
    N_Samples = n,
    P_HighVar = yval2[, 4],
    P_LowVar  = yval2[, 5],
    Node_Percent = yval2[, 6]
  ) %>%
  dplyr::select(Node_ID, Predicted_Class, N_Samples, P_HighVar, P_LowVar, Node_Percent)

library(partykit)

ctree_model = ctree(treg_outcome ~ .,
                    data = df_clustering %>% select(treg_outcome, all_of(param_names)),
                    control = ctree_control(
                      # === Splitting criteria ===
                      mincriterion = 0.90, # 1-p-value threshold (higher = stricter splits)
                      # === Node size controls ===
                      minsplit = 5,               # Min observations to attempt a split
                      minbucket = 5,              # Min observations in terminal node
                      maxdepth = 10,               # Max tree depth (0 = unlimited)
                      testtype = "Univariate"    # Multiple testing correction
                    ))

df_clustering = df_clustering %>% mutate(ctree_node = predict(ctree_model, type = "node"))

library(ggparty)

get_split_info = function(tree_model) {
  # Get the party object
  party_obj = as.party(tree_model)

  # Function to extract split info for each node
  extract_node_info = function(node_id) {
    node = party_obj[[node_id]]

    if(is.terminal(node)) {
      return(paste0("Node ", node_id))
    } else {
      # Get the split variable and threshold
      split_var = names(node$split$varid)[1]
      if(is.null(split_var)) {
        split_var = names(data_in_node(party_obj, node_id))[node$split$varid]
      }

      # For numeric splits
      if(!is.null(node$split$breaks)) {
        split_val = round(node$split$breaks, 2)
        return(paste0(split_var, "\n≤ ", split_val))
      } else {
        # For categorical splits
        return(split_var)
      }
    }
  }

  # Apply to all nodes
  node_labels = sapply(1:length(party_obj), extract_node_info)
  return(node_labels)
}

p = ggparty(ctree_model) +
  geom_edge() +
  geom_edge_label(
    aes(label = breaks_label),  # This shows the threshold values (≤ or >)
    size = 3
  ) +
  geom_node_label(
    aes(label = ifelse(
      is.na(splitvar),  # Terminal nodes
      paste0("Node ", id, "\nN = ", nodesize),
      paste0(splitvar, "\np = ", round(p.value, 4), "\nN = ", nodesize)  # Internal nodes with split variable
    )),
    size = 3.5,
    fontface = "bold",
    nudge_y = 0.02
  ) +
  theme_minimal(base_size = 14) +
  geom_node_plot(
    gglist = list(
      geom_bar(
        aes(
          x = treg_outcome,
          y = after_stat(count / tapply(count, PANEL, sum)[PANEL]),
          fill = treg_outcome
        ),
        show.legend = FALSE
      ),
      scale_fill_manual(values = c("1" = "blue","0" = "gray","-1"='red')),
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)),
      theme_minimal(base_size = 10) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    ),
    shared_axis_labels = FALSE
  )

ggsave(paste0("./tree_tregs_CTREE_",inj_type,".png"), p, width = 28, height = 16, dpi = 300)

node_summary_counts = df_clustering %>%
  group_by(ctree_node, treg_outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ctree_node) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  select(ctree_node, treg_outcome, frac) %>%
  pivot_wider(
    names_from = treg_outcome,
    values_from = frac,
    values_fill = 0
  )

