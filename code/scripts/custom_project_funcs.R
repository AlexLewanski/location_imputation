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


#assign ancestry for each node in each tree based in genome interval info for a
#node's ancestry. This is used when to assign ancestry info in simplified trees
#when the ancestry info is extracted from the unsimplified trees (quick, rough
#explanation ... will update later)
assign_ancestor <- function(tree_interval_info, 
                            ancestry_interval_info, 
                            ancestor_locations,
                            ancestor_location_join_col) {
  
  #checking to see if dataframes have the required columns
  if (!all(c("tree_index", "left_interval", "right_interval", "span") %in% colnames(tree_interval_info)))
    stop('tree_interval_info must have been the following columns: "tree_index", "left_interval", "right_interval", "span"')
  
  if (!all(c("focal_node_simplified", "focal_node_original", "ancestor", "left_interval", "right_interval") %in% colnames(ancestry_interval_info)))
    stop('ancestry_interval_info must have been the following columns: "focal_node_simplified", "focal_node_original", "ancestor", "left_interval", "right_interval"')
  
  
  simplified_node_vec <- unique(ancestry_interval_info$focal_node_simplified)
  
  df_list <- list()
  
  for (focal_node in simplified_node_vec) {
    tree_interval_info_copy <- tree_interval_info
    anc_interv_subset <- ancestry_interval_info[ancestry_interval_info$focal_node_simplified == focal_node,]
    
    tree_interval_info_copy[,'focal_node'] <- focal_node
    tree_interval_info_copy[,'ancestor_id'] <- NA
    for (i in seq_len(nrow(tree_interval_info_copy))) {
      left_index <- max(which(anc_interv_subset$left_interval <= tree_interval_info_copy[i,]$left_interval))
      right_index <- min(which(anc_interv_subset$right_interval >= tree_interval_info_copy[i,]$right_interval))
      
      if (left_index == right_index) {
        tree_interval_info_copy[i,'ancestor_id'] <- anc_interv_subset[left_index,]$ancestor
      } else {
        tree_interval_info_copy[i,'ancestor_id'] <- 'multianc'
      }
    }
    
    df_list[[paste0('node', focal_node)]] <- left_join(tree_interval_info_copy, 
                                                       ancestor_locations, 
                                                       by = c('ancestor_id' = ancestor_location_join_col))
  }
  
  return(df_list)
  
}



########################################################
### VISUALIZATION (AND PROCESSING FOR VISUALIZATION) ###
########################################################

summarize_by_bins <- function(true_locs, inferred_locs, bin_size_vec) {
  #recover()
  
  ### FUNCTION OVERVIEW: ###
  #Using the mean and inferred location for each individual, the function creates
  #a set of bins (based on the user-defined bin size) and then summarizes the
  #true/estimated location pairs in two ways for each bin. First, it calculates
  #the mean Euclidean true/estimated distance. Second, it quantifies the bias in
  #the direction of true --> inferred locations by summing the vectors and dividing
  #by the total number of vectors.
  
  ### THINGS TO POTENTIALLY ADD ###
  #-check to see if the true_locs and inferred_locs have different individuals
  
  
  ### MINOR INPUT CHECKS ###
  if (!all(c('x', 'y', 'indiv') %in% colnames(true_locs)))
    stop('true_locs must have x, y, and indiv columns')
  
  if (!all(c('x', 'y', 'indiv') %in% colnames(inferred_locs)))
    stop('true_locs must have x, y, and indiv columns')
  
  if (!all(inferred_locs$indiv %in% true_locs$indiv))
    stop('Not all the indivs with inferred locations have true locations.')
  
  ### INITIAL PROCESSING OF INPUTS ###
  colnames(true_locs)[colnames(true_locs) == 'x'] <- 'x_true'
  colnames(true_locs)[colnames(true_locs) == 'y'] <- 'y_true'
  colnames(inferred_locs)[colnames(inferred_locs) == 'x'] <- 'x_est'
  colnames(inferred_locs)[colnames(inferred_locs) == 'y'] <- 'y_est'
  
  combined_loc_df <- left_join(true_locs, inferred_locs, by = 'indiv')
  
  
  ### SUMMARIZE BY BINS ###
  bin_summary_list <- list()
  mean_bin_list <- list()
  prop_bin_list <- list()
  
  for (BIN in bin_size_vec) {
    #current attempt at creating nice bins for summarizing based on the user-supplied bin sizes
    x_seq <- seq(floor(min(combined_loc_df$x_true) / BIN) * BIN,ceiling(max(combined_loc_df$x_true) / BIN) * BIN,by = BIN)
    y_seq <- seq(floor(min(combined_loc_df$y_true) / BIN) * BIN,ceiling(max(combined_loc_df$y_true) / BIN) * BIN,by = BIN)
    combined_loc_df$x_bin <- cut(combined_loc_df$x_true, x_seq) #cut(combined_loc_df$x_true, breaks = seq(0, 100, BIN))
    combined_loc_df$y_bin <- cut(combined_loc_df$y_true, y_seq) #cut(combined_loc_df$y_true, breaks = seq(0, 100, BIN))
    combined_loc_df$dist <- apply(combined_loc_df[,c('x_est', 'y_est', 'x_true', 'y_true')], 1, function(m) sqrt( (m[1] - m[3])^2 + (m[2] - m[4])^2 ))
    
    bin_summary_list[[paste0('bin_size', BIN)]] <- combined_loc_df %>% 
      group_by(x_bin, y_bin) %>% 
      mutate(x_disp = x_est - x_true,
             y_disp = y_est - y_true) %>% 
      summarize(mean_dist = mean(dist),
                mean_x = sum(x_disp)/n(),
                mean_y = sum(y_disp)/n(),
                .groups = 'drop') %>% 
      mutate(xmin = as.numeric(sub("^[\\[\\(]([-+]?[0-9]*\\.?[0-9]+),.*", "\\1", x_bin)),
             xmax = as.numeric(sub("^[\\[\\(][^,]+,\\s*([-+]?[0-9]*\\.?[0-9]+).*", "\\1", x_bin)),
             ymin = as.numeric(sub("^[\\[\\(]([-+]?[0-9]*\\.?[0-9]+),.*", "\\1", y_bin)),
             ymax = as.numeric(sub("^[\\[\\(][^,]+,\\s*([-+]?[0-9]*\\.?[0-9]+).*", "\\1", y_bin))) %>% 
      mutate(summed_x_start = (xmax + xmin)/2, #middle of bin in x direction
             summed_y_start = (ymax + ymin)/2) %>% #middle of bin in y direction
      mutate(summed_x_end = summed_x_start + mean_x,
             summed_y_end = summed_y_start + mean_y)
    
    bin_count <- (length(x_seq) - 1)*(length(y_seq) - 1)
    mean_bin_list[[paste0('bin_size', BIN)]] <- nrow(combined_loc_df)/bin_count
    prop_bin_list[[paste0('bin_size', BIN)]] <- nrow(bin_summary_list[[paste0('bin_size', BIN)]])/bin_count
  }
  
  return(
    list(inferred_true_join_df = combined_loc_df, #merged df with inferred and true locations
         bin_summary_list = bin_summary_list, #list of dataframes with bin summaries
         mean_bin_points_vec = unlist(mean_bin_list), #mean number of true locations in each bin
         prop_bins_vec = unlist(prop_bin_list)) #prop of bins with a location in it
  )
  
  ### EXAMPLE BASE VIZ FROM THIS FUNCTION ###
  # ggplot(data = test_summarize$bin_summary_list$bin_size15) +
  #   geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = mean_dist)) +
  #   geom_segment(aes(x = summed_x_start, y = summed_y_start,xend = summed_x_end, yend = summed_y_end),
  #                arrow = arrow())
  
}
