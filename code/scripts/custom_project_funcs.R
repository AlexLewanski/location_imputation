################################################
################################################
### FUNCTIONS FOR SPATIAL IMPUTATION PROJECT ###
################################################
################################################

################
### OVERVIEW ###
################
#This script contains various custom functions for the spatial imputation project
#outside of the core functions that implement the location inference method. Those
#functions are housed in location_est_code.R



#############################
### SIMULATION PROCESSING ###
#############################

pedigree_founder_info <- function(indiv_id, admix_node_df, offspring_df, parent_loc_df) {
  
  pedig_id <- admix_node_df[admix_node_df[,'indiv'] == indiv_id,]$ped_id[1]
  
  counter <- 1
  ancestor_list <- list()
  focal_vec <- pedig_id
  ancestor_list[[counter]] <- focal_vec
  while(length(focal_vec) > 0) {
    focal_vec <- unlist(offspring_df[offspring_df$ped_id %in% focal_vec, c('p1_ped_id', 'p2_ped_id')])
    
    if (length(focal_vec) > 0) ancestor_list[[counter + 1]] <- focal_vec
    counter <- counter + 1
  }
  
  moved_indivs <- ancestor_list[[length(ancestor_list) - 1]]
  
  return(parent_loc_df[parent_loc_df$pedigree_id %in% moved_indivs,])
}