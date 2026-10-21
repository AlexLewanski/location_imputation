#############################################################
#############################################################
### EXAMPLES OF ESTIMATING LOCATIONS IN ADMIXED SCENARIOS ###
#############################################################
#############################################################


################################################
### 2 LOADING AND PROCESSING SIMULATED DATA  ###
################################################

########################
### LOADING PACKAGES ###
########################
library(here)
library(dplyr)
library(ggplot2)
library(ape)
library(cowplot)

#source(here('code', 'scripts', 'location_est_code.R'))

#devtools::install("/Users/alexlewanski/Documents/r_packages/point")
library(point)
source(here('code', 'scripts', 'custom_project_funcs.R'))



#################################
### PROCESSING INFO FOR K = 2 ###
#################################

base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME10000_DISP1.5_CC0.17_K2_Kreps10_sampsize150'
base_file_name_minussamp <- gsub("_sampsize1[0-9]*", "", base_file_name)

files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_files_names <- files_vec[grep(".*nwk.gz$", files_vec)]
names(tree_files_names) <- sub(".*(INDIV[0-9]+).*", "\\1", tree_files_names)

tree_list <- lapply(tree_files_names, function(x) read_trees_from_gz(x))


metadat_files_names <- files_vec[grep(".*txt$", files_vec)]
metadat_list <- lapply(metadat_files_names, function(x) read.delim(x, header = TRUE))
names(metadat_list) <- sapply(metadat_files_names, function(x) gsub("\\.txt", "",gsub(".*sampsize[0-9]+_", "", x)))
metadat_list$location_info$n_node <- paste0('n', metadat_list$location_info$node_id_postsimpl)
metadat_list$focal_inds_node_map$n_node <- paste0('n', metadat_list$focal_inds_node_map$node_id_postsimpl)
metadat_list[['offspring_info']] <- read.delim(here('simulation_output', 'postburnin_sims', paste0('offspring_info_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
metadat_list[['parent_locs']] <- read.delim(here('simulation_output', 'postburnin_sims', paste0('parents_locs_ped_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")


processed_location_mat_list <- list()
admix_indivs <- unique(metadat_list$focal_inds_node_map$indiv_id)
for (INDIV in admix_indivs) {
  focal_nodes <- metadat_list$focal_inds_node_map$node_id_postsimpl[metadat_list$focal_inds_node_map$indiv_id == INDIV]
  
  miss_info <- metadat_list$location_info[metadat_list$location_info$focal_indiv_id == INDIV & !(metadat_list$location_info$node_id_postsimpl %in% focal_nodes),]
  rownames(miss_info) <- miss_info$n_node
  
  processed_location_mat_list[[paste0('INDIV',INDIV)]] <- as.matrix(miss_info[,c('x', 'y')])
}




### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
root_rate_info_list <- list()
admix_indivs <- unique(metadat_list$focal_inds_node_map$indiv_id)

for (INDIV in admix_indivs) {
  focal_nodes <- metadat_list$focal_inds_node_map$n_node[metadat_list$focal_inds_node_map$indiv_id == INDIV]
  
  root_rate_info_list[[paste0('INDIV',INDIV)]] <- root_rate_multitree(tree = tree_list[[paste0('INDIV',INDIV)]], 
                                                                      processed_location_mat_list[[paste0('INDIV',INDIV)]], 
                                                                      focal_nodes)
  message('Finished root and Brownian rate estimation for ', INDIV)
}




k2_est_list <- list()
k2_plot_list <- list()

#color_vec <- c("#ccbadc", "#bba3d0", "#9975b9", "#7647a2", "#551a8b")
color_vec <- rep(c("#ccbadc", "#ccbadc", "#bba3d0", "#bba3d0", "#7647a2"), each = 2)

for (MISSING_INDIV in admix_indivs[1:5]) {
  focal_nodes <- metadat_list$focal_inds_node_map$n_node[metadat_list$focal_inds_node_map$indiv_id == MISSING_INDIV]
  
  k2_process <- process_trees(
    tree = tree_list[[paste0('INDIV',MISSING_INDIV)]], 
    trait = processed_location_mat_list[[paste0('INDIV',MISSING_INDIV)]], 
    brown_rate_mat = root_rate_info_list[[paste0('INDIV',MISSING_INDIV)]]$mean_brown_rate, 
    root_coords_list = root_rate_info_list[[paste0('INDIV',MISSING_INDIV)]]$root_list,
    missing_inds = focal_nodes
  )
  
  
  cond_mean_list <- lapply(k2_process, function(x) x$cond_distr$condmean)
  cond_var_list <- lapply(k2_process, function(x) x$cond_distr$condvar)
  
  k2_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 5,
                                                                       cond_mean = cond_mean_list, 
                                                                       cond_var = cond_var_list, 
                                                                       k = 2, 
                                                                       max_stp = 40, 
                                                                       conv_thresh = 1e-5, 
                                                                       sp_bounds = list(x = c(-5, 105), 
                                                                                        y = c(-5, 105)),
                                                                       starting_location_seed = NULL)
  
  cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
  
  plot_color <- color_vec[index]
  
  ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
                                         admix_node_df = metadat_list$focal_inds_parent_info, 
                                         offspring_df = metadat_list$offspring_info, 
                                         parent_loc_df = metadat_list$parent_locs)
  
  k2_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
    annotate("rect", 
             xmin = 0, xmax = 100,       # Start and end on X axis
             ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
             fill = "#e5e5e5", #"#fac9c3",       # Background color
             alpha = 0.4) +
    geom_point(data = cond_mean_df,
               aes(x = x, y = y), shape = 8,
               color = plot_color, size = 0.75, alpha = 0.4) +
    geom_point(data = as.data.frame(processed_location_mat_list[[paste0('INDIV',MISSING_INDIV)]]),
               aes(x = x, y = y),
               color = '#999999', size = 0.7, alpha = 0.3) +
    geom_point(data = k2_est_list[[paste0('IND', MISSING_INDIV)]]$max_lik_xmat_list[[1]],
               aes(x, y), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
    geom_point(data = ancestor_locs,
               aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
          plot.title = element_text(hjust = 0.5, size = 15))
  
  message('FINISHED ', MISSING_INDIV)
  
}


k2_plot_title <- title <- cowplot::ggdraw() + cowplot::draw_label("K = 2", fontface = 'bold', x = 0.5, hjust = 0.5)
k2_multipan <- cowplot::plot_grid(plotlist = k2_plot_list[c(2,3)], nrow = 1)

k2_multipan_withtitle <- cowplot::plot_grid(k2_plot_title, k2_multipan, ncol = 1, rel_heights = c(0.1, 1))





cowplot::plot_grid(plotlist = k2_plot_list, nrow = 2)






#################################
### PROCESSING INFO FOR K = 3 ###
#################################

base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME10000_DISP1.5_CC0.17_K3_Kreps10_sampsize150'
base_file_name_minussamp <- gsub("_sampsize1[0-9]*", "", base_file_name)

files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_files_names <- files_vec[grep(".*nwk.gz$", files_vec)]
names(tree_files_names) <- sub(".*(INDIV[0-9]+).*", "\\1", tree_files_names)

tree_list <- lapply(tree_files_names, function(x) read_trees_from_gz(x))


metadat_files_names <- files_vec[grep(".*txt$", files_vec)]
metadat_list <- lapply(metadat_files_names, function(x) read.delim(x, header = TRUE))
names(metadat_list) <- sapply(metadat_files_names, function(x) gsub("\\.txt", "",gsub(".*sampsize[0-9]+_", "", x)))
metadat_list$location_info$n_node <- paste0('n', metadat_list$location_info$node_id_postsimpl)
metadat_list$focal_inds_node_map$n_node <- paste0('n', metadat_list$focal_inds_node_map$node_id_postsimpl)
metadat_list[['offspring_info']] <- read.delim(here('simulation_output', 'postburnin_sims', paste0('offspring_info_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
metadat_list[['parent_locs']] <- read.delim(here('simulation_output', 'postburnin_sims', paste0('parents_locs_ped_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")


processed_location_mat_list <- list()
admix_indivs <- unique(metadat_list$focal_inds_node_map$indiv_id)
for (INDIV in admix_indivs) {
  focal_nodes <- metadat_list$focal_inds_node_map$node_id_postsimpl[metadat_list$focal_inds_node_map$indiv_id == INDIV]
  
  miss_info <- metadat_list$location_info[metadat_list$location_info$focal_indiv_id == INDIV & !(metadat_list$location_info$node_id_postsimpl %in% focal_nodes),]
  rownames(miss_info) <- miss_info$n_node
  
  processed_location_mat_list[[paste0('INDIV',INDIV)]] <- as.matrix(miss_info[,c('x', 'y')])
}




### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
root_rate_info_list <- list()
admix_indivs <- unique(metadat_list$focal_inds_node_map$indiv_id)

for (INDIV in admix_indivs) {
  focal_nodes <- metadat_list$focal_inds_node_map$n_node[metadat_list$focal_inds_node_map$indiv_id == INDIV]
  
  root_rate_info_list[[paste0('INDIV',INDIV)]] <- root_rate_multitree(tree = tree_list[[paste0('INDIV',INDIV)]], 
                                                                      processed_location_mat_list[[paste0('INDIV',INDIV)]], 
                                                                      focal_nodes)
  message('Finished root and Brownian rate estimation for ', INDIV)
}




k3_est_list <- list()
k3_plot_list <- list()

#color_vec <- c("#ccbadc", "#bba3d0", "#9975b9", "#7647a2", "#551a8b")
color_vec <- rep(c("#ccbadc", "#ccbadc", "#bba3d0", "#bba3d0", "#7647a2"), each = 2)

for (MISSING_INDIV in admix_indivs[1:6]) {
  focal_nodes <- metadat_list$focal_inds_node_map$n_node[metadat_list$focal_inds_node_map$indiv_id == MISSING_INDIV]
  
  k3_process <- process_trees(
    tree = tree_list[[paste0('INDIV',MISSING_INDIV)]], 
    trait = processed_location_mat_list[[paste0('INDIV',MISSING_INDIV)]], 
    brown_rate_mat = root_rate_info_list[[paste0('INDIV',MISSING_INDIV)]]$mean_brown_rate, 
    root_coords_list = root_rate_info_list[[paste0('INDIV',MISSING_INDIV)]]$root_list,
    missing_inds = focal_nodes
  )
  
  
  cond_mean_list <- lapply(k3_process, function(x) x$cond_distr$condmean)
  cond_var_list <- lapply(k3_process, function(x) x$cond_distr$condvar)
  
  k3_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 5,
                                                                       cond_mean = cond_mean_list, 
                                                                       cond_var = cond_var_list, 
                                                                       k = 3, 
                                                                       max_stp = 40, 
                                                                       conv_thresh = 1e-5, 
                                                                       sp_bounds = list(x = c(-5, 105), 
                                                                                        y = c(-5, 105)),
                                                                       starting_location_seed = NULL)
  
  cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
  
  plot_color <- color_vec[index]
  
  ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
                                         admix_node_df = metadat_list$focal_inds_parent_info, 
                                         offspring_df = metadat_list$offspring_info, 
                                         parent_loc_df = metadat_list$parent_locs)
  
  k3_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
    annotate("rect", 
             xmin = 0, xmax = 100,       # Start and end on X axis
             ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
             fill = "#e5e5e5", #"#fac9c3",       # Background color
             alpha = 0.4) +
    geom_point(data = cond_mean_df,
               aes(x = x, y = y), shape = 8,
               color = plot_color, size = 0.75, alpha = 0.4) +
    geom_point(data = as.data.frame(processed_location_mat_list[[paste0('INDIV',MISSING_INDIV)]]),
               aes(x = x, y = y),
               color = '#999999', size = 0.7, alpha = 0.3) +
    geom_point(data = k3_est_list[[paste0('IND', MISSING_INDIV)]]$max_lik_xmat_list[[1]],
               aes(x, y), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
    geom_point(data = ancestor_locs,
               aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
          plot.title = element_text(hjust = 0.5, size = 15))
  
  message('FINISHED ', MISSING_INDIV)
  
}


k3_plot_title <- title <- cowplot::ggdraw() + cowplot::draw_label("K = 3", fontface = 'bold', x = 0.5, hjust = 0.5)
k3_multipan <- cowplot::plot_grid(plotlist = k3_plot_list, nrow = 2)

k3_multipan_withtitle <- cowplot::plot_grid(k3_plot_title, k3_multipan, ncol = 2, rel_heights = c(0.1, 1))






#### IGNORE BELOW HERE ###




### LOADING TREES AND OTHER INFO ###
# base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME10000_DISP2.0_CC0.17_K2_Kreps10_sampsize150'
# base_file_name_minussamp <- gsub("_sampsize1[0-9]*", "", base_file_name)
# 
# files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# 
# tree_files_vec <- files_vec[grep(".*nwk$", files_vec)]
# tree_number_numeric <- as.numeric(gsub(paste0('.*', base_file_name, '|\\.nwk'),"", tree_files_vec))
# 
# #tree_list_init <- lapply(files_vec[grep(".*nwk$", files_vec)],
# #                         function(x) ape::read.tree(x))
# #tree_list <- tree_list_init[order(tree_number_numeric)]
# 
# tree_list <- read_trees_from_gz(file_name = here('simulation_output', 'processed_treeseq', base_file_name, paste0(base_file_name, '.nwk.gz')))
# 
# 
# #con <- gzfile(here('simulation_output', 'processed_treeseq', base_file_name, paste0(base_file_name, '.nwk.gz')), open = "rt")
# #tree_strings <- readLines(con)
# #close(con)
# 
# #tree_list_alt <- lapply(
# #  tree_strings,
# #  function(x) ape::read.tree(text = x)
# #)
# 
# 
# indiv_node_info <- read.delim(here('simulation_output',
#                                    'processed_treeseq', 
#                                    base_file_name,
#                                    paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
# indiv_node_info$n_node <- paste0('n', indiv_node_info$node)
# 
# 
# 
# offspring_info <- read.delim(here('simulation_output', 
#                                   'postburnin_sims',
#                                   paste0('offspring_info_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# 
# parent_locs <- read.delim(here('simulation_output', 
#                                'postburnin_sims',
#                                paste0('parents_locs_ped_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# 
# admix_node_info <- read.delim(here('simulation_output',
#                                    'processed_treeseq',
#                                    base_file_name,
#                                    paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)
# 
# 
# 
# ### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
# admix_indivs <- unique(admix_node_info$indiv)
# admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]
# 
# miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
# rownames(miss_info) <- miss_info$n_node
# miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# 
# #tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# tree_list_subset1 <- tree_list
# root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
#                                       miss_info_mat, 
#                                       admix_nodes)
# 
# 
# ### LOCATION ESTIMATION AND PLOTTING ###
# 
# k2_est_list <- list()
# k2_plot_list <- list()
# 
# #color_vec <- c("#ccbadc", "#bba3d0", "#9975b9", "#7647a2", "#551a8b")
# color_vec <- c("#ccbadc", "#ccbadc", "#bba3d0", "#bba3d0", "#7647a2")
# 
# for (MISSING_INDIV in admix_indivs) {
#   index <- which(admix_indivs == MISSING_INDIV)
# 
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
#   REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
# 
# 
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# 
#   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# 
#   k2_process <- process_trees(
#     tree = tree_list_subset,
#     trait = miss_info_mat,
#     brown_rate_mat = root_rate_info$mean_brown_rate,
#     root_coords_list = root_rate_info$root_list,
#     missing_inds = MISSING_NODES
#   )
# 
# 
#   cond_mean_list <- lapply(k2_process, function(x) x$cond_distr$condmean)
#   cond_var_list <- lapply(k2_process, function(x) x$cond_distr$condvar)
# 
#   k2_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 5,
#                                                                        cond_mean = cond_mean_list,
#                                                                        cond_var = cond_var_list,
#                                                                        k = 2,
#                                                                        max_stp = 40,
#                                                                        conv_thresh = 1e-5,
#                                                                        sp_bounds = list(x = c(-5, 105),
#                                                                                         y = c(-5, 105)),
#                                                                        starting_location_seed = NULL)
# 
#   cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
# 
#   plot_color <- color_vec[index]
# 
#   ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV,
#                         admix_node_df = admix_node_info,
#                         offspring_df = offspring_info,
#                         parent_loc_df = parent_locs)
# 
#   k2_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
#     annotate("rect",
#              xmin = 0, xmax = 100,       # Start and end on X axis
#              ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
#              fill = "#e5e5e5", #"#fac9c3",       # Background color
#              alpha = 0.4) +
#     geom_point(data = cond_mean_df,
#                aes(x = x, y = y), shape = 8,
#                color = plot_color, size = 0.75, alpha = 0.4) +
#     geom_point(data = as.data.frame(miss_info_mat),
#                aes(x = x, y = y),
#                color = '#999999', size = 0.7, alpha = 0.3) +
#     geom_point(data = k2_est_list[[paste0('IND', MISSING_INDIV)]]$max_lik_xmat_list[[1]],
#                aes(x, y), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
#     geom_point(data = ancestor_locs,
#                aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
#           panel.grid.major = element_blank(),
#           panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
#           plot.title = element_text(hjust = 0.5, size = 15))
# 
#   message('FINISHED ', MISSING_INDIV)
# 
# }
# 
# 
# k2_plot_title <- title <- cowplot::ggdraw() + cowplot::draw_label("K = 2", fontface = 'bold', x = 0.5, hjust = 0.5)
# k2_multipan <- cowplot::plot_grid(plotlist = k2_plot_list[c(2,3)], nrow = 1)
# 
# k2_multipan_withtitle <- cowplot::plot_grid(k2_plot_title, k2_multipan, ncol = 1, rel_heights = c(0.1, 1))
# 
# 
# 
# 
# 
# ### LOADING TREES AND OTHER INFO ###
# base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME10000_DISP2.0_CC0.17_K2_Kreps10_sampsize150'
# base_file_name_minussamp <- gsub("_sampsize1[0-9]*", "", base_file_name)
# 
# files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# 
# tree_files_vec <- files_vec[grep(".*nwk$", files_vec)]
# tree_number_numeric <- as.numeric(gsub(paste0('.*', base_file_name, '|\\.nwk'),"", tree_files_vec))
# 
# #tree_list_init <- lapply(files_vec[grep(".*nwk$", files_vec)],
# #                         function(x) ape::read.tree(x))
# #tree_list <- tree_list_init[order(tree_number_numeric)]
# 
# tree_list <- read_trees_from_gz(file_name = here('simulation_output', 'processed_treeseq', base_file_name, paste0(base_file_name, '.nwk.gz')))
# 
# 
# #con <- gzfile(here('simulation_output', 'processed_treeseq', base_file_name, paste0(base_file_name, '.nwk.gz')), open = "rt")
# #tree_strings <- readLines(con)
# #close(con)
# 
# #tree_list_alt <- lapply(
# #  tree_strings,
# #  function(x) ape::read.tree(text = x)
# #)
# 
# 
# indiv_node_info <- read.delim(here('simulation_output',
#                                    'processed_treeseq', 
#                                    base_file_name,
#                                    paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
# indiv_node_info$n_node <- paste0('n', indiv_node_info$node)
# 
# 
# 
# offspring_info <- read.delim(here('simulation_output', 
#                                   'postburnin_sims',
#                                   paste0('offspring_info_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# 
# parent_locs <- read.delim(here('simulation_output', 
#                                'postburnin_sims',
#                                paste0('parents_locs_ped_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# 
# admix_node_info <- read.delim(here('simulation_output',
#                                    'processed_treeseq',
#                                    base_file_name,
#                                    paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)
# 
# 
# 
# ### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
# admix_indivs <- unique(admix_node_info$indiv)
# admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]
# 
# miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
# rownames(miss_info) <- miss_info$n_node
# miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# 
# #tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# tree_list_subset1 <- tree_list
# root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
#                                       miss_info_mat, 
#                                       admix_nodes)
# 
# 
# ### LOCATION ESTIMATION AND PLOTTING ###
# 
# k2_est_list <- list()
# k2_plot_list <- list()
# 
# #color_vec <- c("#ccbadc", "#bba3d0", "#9975b9", "#7647a2", "#551a8b")
# color_vec <- c("#ccbadc", "#ccbadc", "#bba3d0", "#bba3d0", "#7647a2")
# 
# for (MISSING_INDIV in admix_indivs) {
#   index <- which(admix_indivs == MISSING_INDIV)
#   
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
#   REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
#   
#   
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
#   
#   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
#   
#   k2_process <- process_trees(
#     tree = tree_list_subset,
#     trait = miss_info_mat,
#     brown_rate_mat = root_rate_info$mean_brown_rate,
#     root_coords_list = root_rate_info$root_list,
#     missing_inds = MISSING_NODES
#   )
#   
#   
#   cond_mean_list <- lapply(k2_process, function(x) x$cond_distr$condmean)
#   cond_var_list <- lapply(k2_process, function(x) x$cond_distr$condvar)
#   
#   k2_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 5,
#                                                                        cond_mean = cond_mean_list,
#                                                                        cond_var = cond_var_list,
#                                                                        k = 2,
#                                                                        max_stp = 40,
#                                                                        conv_thresh = 1e-5,
#                                                                        sp_bounds = list(x = c(-5, 105),
#                                                                                         y = c(-5, 105)),
#                                                                        starting_location_seed = NULL)
#   
#   cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
#   
#   plot_color <- color_vec[index]
#   
#   ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV,
#                                          admix_node_df = admix_node_info,
#                                          offspring_df = offspring_info,
#                                          parent_loc_df = parent_locs)
#   
#   k2_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
#     annotate("rect",
#              xmin = 0, xmax = 100,       # Start and end on X axis
#              ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
#              fill = "#e5e5e5", #"#fac9c3",       # Background color
#              alpha = 0.4) +
#     geom_point(data = cond_mean_df,
#                aes(x = x, y = y), shape = 8,
#                color = plot_color, size = 0.75, alpha = 0.4) +
#     geom_point(data = as.data.frame(miss_info_mat),
#                aes(x = x, y = y),
#                color = '#999999', size = 0.7, alpha = 0.3) +
#     geom_point(data = k2_est_list[[paste0('IND', MISSING_INDIV)]]$max_lik_xmat_list[[1]],
#                aes(x, y), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
#     geom_point(data = ancestor_locs,
#                aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
#           panel.grid.major = element_blank(),
#           panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
#           plot.title = element_text(hjust = 0.5, size = 15))
#   
#   message('FINISHED ', MISSING_INDIV)
#   
# }
# 
# 
# k2_plot_title <- title <- cowplot::ggdraw() + cowplot::draw_label("K = 2", fontface = 'bold', x = 0.5, hjust = 0.5)
# k2_multipan <- cowplot::plot_grid(plotlist = k2_plot_list[c(2,3)], nrow = 1)
# 
# k2_multipan_withtitle <- cowplot::plot_grid(k2_plot_title, k2_multipan, ncol = 1, rel_heights = c(0.1, 1))
# 
# 
# 
# 
# 
# 
# 
# # #################################
# # ### PROCESSING INFO FOR K = 3 ###
# # #################################
# # 
# # ### LOADING TREES AND OTHER INFO ###
# # base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME10000_DISP2.0_CC0.17_K3_Kreps10_sampsize150'
# # base_file_name_minussamp <- gsub("_sampsize1[0-9]*", "", base_file_name)
# # 
# # #files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# # #tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
# # #                    function(x) ape::read.tree(x)) 
# # tree_list <- read_trees_from_gz(file_name = here('simulation_output', 'processed_treeseq', base_file_name, paste0(base_file_name, '.nwk.gz')))
# # 
# # 
# # indiv_node_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq', 
# #                                    base_file_name,
# #                                    paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
# # indiv_node_info$n_node <- paste0('n', indiv_node_info$node)
# # 
# # 
# # 
# # offspring_info <- read.delim(here('simulation_output', 
# #                                   'postburnin_sims',
# #                                   paste0('offspring_info_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# # 
# # parent_locs <- read.delim(here('simulation_output', 
# #                                'postburnin_sims',
# #                                paste0('parents_locs_ped_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# # 
# # admix_node_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq',
# #                                    base_file_name,
# #                                    paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)
# # 
# # 
# # # files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# # # tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
# # #                     function(x) ape::read.tree(x)) 
# # 
# # 
# # 
# # ### ESTIMATION ###
# # admix_indivs <- unique(admix_node_info$indiv)
# # admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]
# # 
# # miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
# # rownames(miss_info) <- miss_info$n_node
# # miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# # 
# # #tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# # tree_list_subset1 <- tree_list
# # root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
# #                                       miss_info_mat, 
# #                                       admix_nodes)
# # 
# # 
# # ### LOCATION ESTIMATION AND PLOTTING ###
# # 
# # k3_est_list <- list()
# # k3_plot_list <- list()
# # 
# # color_vec <- c("#c7e7c1", "#abdba2", "#8fcf83", "#74c365", "#68af5a")
# # 
# # for (MISSING_INDIV in admix_indivs[2:3]) {
# #   index <- which(admix_indivs == MISSING_INDIV)
# #   
# #   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
# #   REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
# #   
# #   
# #   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
# #   rownames(miss_info) <- miss_info$n_node
# #   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# #   
# #   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# #   
# #   k3_process <- process_trees(
# #     tree = tree_list_subset, 
# #     trait = miss_info_mat, 
# #     brown_rate_mat = root_rate_info$mean_brown_rate, 
# #     root_coords_list = root_rate_info$root_list,
# #     missing_inds = MISSING_NODES
# #   )
# #   
# #   
# #   cond_mean_list <- lapply(k3_process, function(x) x$cond_distr$condmean)
# #   cond_var_list <- lapply(k3_process, function(x) x$cond_distr$condvar)
# #   
# #   k3_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 10,
# #                                                                        cond_mean = cond_mean_list, 
# #                                                                        cond_var = cond_var_list, 
# #                                                                        k = 3, 
# #                                                                        max_stp = 25, 
# #                                                                        conv_thresh = 1e-5, 
# #                                                                        sp_bounds = list(x = c(-5, 105), 
# #                                                                                         y = c(-5, 105)),
# #                                                                        starting_location_seed = NULL)
# #   
# #   cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
# #   
# #   plot_color <- color_vec[index]
# #   
# #   ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
# #                                          admix_node_df = admix_node_info, 
# #                                          offspring_df = offspring_info, 
# #                                          parent_loc_df = parent_locs)
# #   
# #   
# #   k3_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
# #     annotate("rect", 
# #              xmin = 0, xmax = 100,       # Start and end on X axis
# #              ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #              fill = "#e5e5e5", #"#fac9c3",       # Background color
# #              alpha = 0.4) +
# #     geom_point(data = cond_mean_df,
# #                aes(x = x, y = y), shape = 8,
# #                color = plot_color, size = 0.75, alpha = 0.4) +
# #     geom_point(data = as.data.frame(miss_info_mat),
# #                aes(x = x, y = y),
# #                color = '#999999', size = 0.7, alpha = 0.3) +
# #     geom_point(data = k3_est_list[[paste0('IND', MISSING_INDIV)]]$max_lik_xmat_list[[1]],
# #                aes(x, y), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
# #     geom_point(data = ancestor_locs,
# #                aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
# #     theme_bw() +
# #     theme(panel.grid.minor = element_blank(),
# #           #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #           panel.grid.major = element_blank(),
# #           panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #           plot.title = element_text(hjust = 0.5, size = 15))
# #   message('FINISHED ', MISSING_INDIV)
# #   
# # }
# # 
# # k3_plot_title <- title <- cowplot::ggdraw() + cowplot::draw_label("K = 3", fontface = 'bold', x = 0.5, hjust = 0.5)
# # k3_multipan <- cowplot::plot_grid(plotlist = k3_plot_list, nrow = 1)
# # 
# # k3_multipan_withtitle <- cowplot::plot_grid(k3_plot_title, k3_multipan, ncol = 1, rel_heights = c(0.1, 1))
# # 
# # 
# # 
# # 
# # #################################
# # ### PROCESSING INFO FOR K = 4 ###
# # #################################
# # 
# # ### LOADING TREES AND OTHER INFO ###
# # base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME5000_DISP1.5_CC0.15_K4_Kreps4_sampsize150'
# # base_file_name_minussamp <- gsub("_sampsize1[0-9]*", "", base_file_name)
# # 
# # #files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# # #tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
# # #                    function(x) ape::read.tree(x)) 
# # 
# # tree_list <- read_trees_from_gz(file_name = here('simulation_output', 'processed_treeseq', base_file_name, paste0(base_file_name, '.nwk.gz')))
# # 
# # 
# # indiv_node_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq', 
# #                                    base_file_name,
# #                                    paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
# # indiv_node_info$n_node <- paste0('n', indiv_node_info$node)
# # 
# # 
# # 
# # offspring_info <- read.delim(here('simulation_output', 
# #                                   'postburnin_sims',
# #                                   paste0('offspring_info_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# # 
# # parent_locs <- read.delim(here('simulation_output', 
# #                                'postburnin_sims',
# #                                paste0('parents_locs_ped_', base_file_name_minussamp, '.txt')), header = TRUE, sep = " ")
# # 
# # admix_node_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq',
# #                                    base_file_name,
# #                                    paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)
# # 
# # # files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# # # tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
# # #                     function(x) ape::read.tree(x)) 
# # 
# # 
# # 
# # ### ESTIMATION ###
# # admix_indivs <- unique(admix_node_info$indiv)
# # admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]
# # 
# # miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
# # rownames(miss_info) <- miss_info$n_node
# # miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# # 
# # #tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# # tree_list_subset1 <- tree_list
# # root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
# #                                       miss_info_mat, 
# #                                       admix_nodes)
# # 
# # 
# # ### LOCATION ESTIMATION AND PLOTTING ###
# # 
# # k4_est_list <- list()
# # k4_plot_list <- list()
# # 
# # color_vec <- c("#ffcdb3", "#ffb48e", "#ff9b68", "#ff8243", "#e5753c")
# # 
# # for (MISSING_INDIV in admix_indivs) {
# #   index <- which(admix_indivs == MISSING_INDIV)
# #   
# #   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
# #   REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
# #   
# #   
# #   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
# #   rownames(miss_info) <- miss_info$n_node
# #   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# #   
# #   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# #   
# #   k4_process <- process_trees(
# #     tree = tree_list_subset, 
# #     trait = miss_info_mat, 
# #     brown_rate_mat = root_rate_info$mean_brown_rate, 
# #     root_coords_list = root_rate_info$root_list,
# #     missing_inds = MISSING_NODES
# #   )
# #   
# #   
# #   cond_mean_list <- lapply(k4_process, function(x) x$cond_distr$condmean)
# #   cond_var_list <- lapply(k4_process, function(x) x$cond_distr$condvar)
# #   
# #   k4_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 3,
# #                                                                        cond_mean = cond_mean_list, 
# #                                                                        cond_var = cond_var_list, 
# #                                                                        k = 4, 
# #                                                                        max_stp = 25, 
# #                                                                        conv_thresh = 1e-5, 
# #                                                                        sp_bounds = list(x = c(-5, 105), 
# #                                                                                         y = c(-5, 105)),
# #                                                                        starting_location_seed = NULL)
# #   
# #   cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
# #   
# #   plot_color <- color_vec[index]
# #   
# #   ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
# #                                          admix_node_df = admix_node_info, 
# #                                          offspring_df = offspring_info, 
# #                                          parent_loc_df = parent_locs)
# #   
# #   
# #   k4_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
# #     annotate("rect", 
# #              xmin = 0, xmax = 100,       # Start and end on X axis
# #              ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #              fill = "#e5e5e5", #"#fac9c3",       # Background color
# #              alpha = 0.4) +
# #     geom_point(data = cond_mean_df,
# #                aes(x = x, y = y), shape = 8,
# #                color = plot_color, size = 0.75, alpha = 0.4) +
# #     geom_point(data = as.data.frame(miss_info_mat),
# #                aes(x = x, y = y),
# #                color = '#999999', size = 0.7, alpha = 0.3) +
# #     geom_point(data = k4_est_list[[paste0('IND', MISSING_INDIV)]]$max_lik_xmat_list[[1]],
# #                aes(x, y), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
# #     geom_point(data = ancestor_locs,
# #                aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
# #     theme_bw() +
# #     theme(panel.grid.minor = element_blank(),
# #           #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #           panel.grid.major = element_blank(),
# #           panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #           plot.title = element_text(hjust = 0.5, size = 15))
# #   message('FINISHED ', MISSING_INDIV)
# #   
# # }
# # 
# # k4_plot_title <- title <- cowplot::ggdraw() + cowplot::draw_label("k4", fontface = 'bold', x = 0.5, hjust = 0.5)
# # k4_multipan <- cowplot::plot_grid(plotlist = k4_plot_list, nrow = 1)
# # 
# # k4_multipan_withtitle <- cowplot::plot_grid(k4_plot_title, k4_multipan, ncol = 1, rel_heights = c(0.1, 1))
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ancestry_df_full <- left_join(ancestry_prop, parent_locs, by = c("ancestor_id" = "pedigree_id")) %>% 
# #   group_by(indiv, location_id) %>% 
# #   summarize(subpop = first(subpop),
# #             x = first(x),
# #             y = first(y),
# #             count = sum(count), .groups = 'drop') %>% 
# #   group_by(indiv) %>% 
# #   mutate(prop = count/sum(count)) %>% 
# #   ungroup()
# # 
# # k3_plot_list$IND152 + 
# #   geom_point(data = ancestry_df_full %>% 
# #                filter(subpop == 4),
# #              aes(x = x, y = y, size = prop)) +
# #   scale_size_area(max_size = 10, limits = c(0, 1))
# # 
# # 
# # true_loc <- ancestry_df_full %>% 
# #   filter(subpop == 4) %>% 
# #   select('x', 'y') %>% 
# #   as.matrix()
# # 
# # estimation_mat <- matrix(k3_est_list$IND152$top_x, 
# #        ncol = 2, 
# #        byrow = TRUE)
# # 
# # 
# # ranked_match_up(ref_location = true_loc, 
# #                 compare_location = estimation_mat)
# #   
# #  
# # 
# # 
# # #######
# # #######
# # #######
# # 
# # k3cond_mean_list <- list()
# # k3cond_var_list <- list()
# # 
# # for (MISSING_INDIV in admix_indivs) {
# #   index <- which(admix_indivs == MISSING_INDIV)
# #   
# #   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
# #   REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
# #   
# #   
# #   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
# #   rownames(miss_info) <- miss_info$n_node
# #   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# #   
# #   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# #   
# #   k3_process <- process_trees(
# #     tree = tree_list_subset, 
# #     trait = miss_info_mat, 
# #     brown_rate_mat = root_rate_info$mean_brown_rate, 
# #     root_coords_list = root_rate_info$root_list,
# #     missing_inds = MISSING_NODES
# #   )
# #   
# #   
# #   k3cond_mean_list[[paste0('ind', MISSING_INDIV)]] <- lapply(k3_process, function(x) x$cond_distr$condmean)
# #   k3cond_var_list[[paste0('ind', MISSING_INDIV)]] <- lapply(k3_process, function(x) x$cond_distr$condvar)
# # }
# # 
# # k3_est_list <- list()
# # for (MISSING_INDIV in names(k3cond_mean_list)) {
# #   k3_est_list[[MISSING_INDIV]] <- em_hardassign_wrapper(times = 5,
# #                                                                        cond_mean = k3cond_mean_list[[MISSING_INDIV]], 
# #                                                                        cond_var = k3cond_var_list[[MISSING_INDIV]],
# #                                                                        k = 3, 
# #                                                                        max_stp = 25, 
# #                                                                        conv_thresh = 1e-5, 
# #                                                                        sp_bounds = list(x = c(-5, 105), 
# #                                                                                         y = c(-5, 105)),
# #                                                                        starting_location_seed = NULL)
# #   
# # }
# # 
# # 
# # k3_est_list$ind152$extra_info$max_lik_group_membership_list_nonreorder
# # 
# # k3_max_lik_xmat_list <- lapply(k3_est_list, function(x) x$max_lik_xmat_list)
# # k3_max_lik_group_membership_list <- lapply(k3_est_list, function(x) x$extra_info$max_lik_group_membership_list_nonreorder)
# # #k3_max_lik_group_membership_list <- lapply(k3_est_list, function(x) x$max_lik_group_membership_list)
# # 
# # k3_max_lik_xmat_list
# # 
# # ancestor_locs_summarize_list <- list()
# # for (MISSING_INDIV in admix_indivs) {
# #   ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
# #                                          admix_node_df = admix_node_info, 
# #                                          offspring_df = offspring_info, 
# #                                          parent_loc_df = parent_locs)
# #   
# #   ancestor_locs_summarize_list[[paste0('ind', MISSING_INDIV)]] <- ancestor_locs %>% 
# #     mutate(location_id_plus = location_id + 1) %>% 
# #     group_by(location_id_plus) %>%
# #     summarize(x = mean(x),
# #               y = mean(y),
# #               .groups = 'drop')
# # }
# # 
# # 
# # 
# # 
# # true_location_matchup_list <- list()
# # for (MISSING_INDIV in names(k3_max_lik_xmat_list)) {
# #   true_location_matchup_list[[MISSING_INDIV]] <- ranked_match_up_multi(list(as.matrix(ancestor_locs_summarize_list[[MISSING_INDIV]][,c('x', 'y')]),
# #                                                                             k3_max_lik_xmat_list[[MISSING_INDIV]][[1]][,c(1,2)])
# #                                                                        )
# #   
# # }
# # 
# # 
# # k3_max_lik_group_membership_list_relabel_list <- list()
# # k3_max_lik_group_membership_list_relabel_list_df <- list()
# # for (MISSING_INDIV in names(true_location_matchup_list)) {
# #   k3_max_lik_group_membership_list_relabel_list[[MISSING_INDIV]] <- convert_labels(old_label = k3_max_lik_group_membership_list[[MISSING_INDIV]][[1]], 
# #                                                              label_map_original = true_location_matchup_list[[MISSING_INDIV]][[2]]$original_group, 
# #                                                              label_map_new = true_location_matchup_list[[MISSING_INDIV]][[2]]$ref_group)
# #   
# #   tree_ind <- as.numeric(gsub('_n.*', '', names(k3_max_lik_group_membership_list_relabel_list[[MISSING_INDIV]])))
# #   node <- as.numeric(gsub('[0-9]*_n', '', names(k3_max_lik_group_membership_list_relabel_list[[MISSING_INDIV]])))
# #   k3_max_lik_group_membership_list_relabel_list_df[[MISSING_INDIV]] <- data.frame(indiv = as.numeric(gsub('ind', '', MISSING_INDIV)),
# #                                                                                   focal_node = node,
# #             tree_ind = tree_ind,
# #             group = k3_max_lik_group_membership_list_relabel_list[[MISSING_INDIV]]) 
# # }
# # 
# # 
# # k3_max_lik_group_membership_list_relabel_list_df
# # 
# # 
# # 
# # tree_interval_info <- read.delim(here('simulation_output',
# #                                       'processed_treeseq',
# #                                       base_file_name,
# #                                       'tree_interval_info.txt'), header = TRUE)
# # 
# # 
# # node_founder_ancestry_intervals <- read.delim(here('simulation_output',
# #                                                    'processed_treeseq',
# #                                                    base_file_name,
# #                                                    'node_founder_ancestry_intervals.txt'), header = TRUE)
# # 
# # 
# # test_assign <- assign_ancestor(tree_interval_info = tree_interval_info, 
# #                                ancestry_interval_info = node_founder_ancestry_intervals,
# #                                ancestor_locations = parent_locs,
# #                                ancestor_location_join_col = 'pedigree_id')
# # 
# # 
# # 
# # true_node_anc <- assign_ancestor(tree_interval_info = tree_interval_info, 
# #                                ancestry_interval_info = node_founder_ancestry_intervals,
# #                                ancestor_locations = parent_locs,
# #                                ancestor_location_join_col = 'pedigree_id')
# # true_node_anc_combine <- do.call(rbind, test_assign) %>% 
# #   left_join(., admix_node_info, by = c("focal_node" = "node"))
# # 
# # # ancestry_prop <- test_assign_combine %>% 
# # #   group_by(indiv, ancestor_id) %>% 
# # #   summarize(count = n()) %>% 
# # #   ungroup() %>% 
# # #   group_by(indiv) %>% 
# # #   mutate(prop = count/sum(count)) %>% 
# # #   ungroup()
# # 
# # 
# # 
# # true_node_anc_combine %>% 
# #   filter(indiv == 152) %>% 
# #   group_by(focal_node) %>% 
# #   arrange(tree_index) %>% 
# #   mutate(right = cumsum(span),
# #          left = lag(right, default = 0),
# #          .groups = 'drop') %>% 
# #   #filter(focal_node == 305) %>% 
# #   ggplot(aes(x = left, xend = right,
# #              y = focal_node, yend = focal_node,
# #              color = as.character(location_id))) +
# #   geom_segment(linewidth = 15, lineend = "butt")# +
# #   #scale_x_continuous(breaks = seq(0, 40, 5)) +
# #   #scale_y_continuous(NULL, breaks = NULL) +
# #   #theme_minimal() +
# #   #labs(x = "Position", color = "Value")
# # 
# # tree_interval_info$tree_ind <- seq_len(nrow(tree_interval_info))
# # 
# # k3_max_lik_group_membership_list_relabel_list_df$ind152 %>% 
# #   left_join(., tree_interval_info, by = 'tree_ind') %>% 
# #   group_by(focal_node) %>% 
# #   arrange(tree_ind) %>% 
# #   mutate(right = cumsum(span),
# #          left = lag(right, default = 0),
# #          .groups = 'drop') %>% 
# #   #filter(focal_node == 305) %>% 
# #   ggplot(aes(x = left, xend = right,
# #              y = focal_node, yend = focal_node,
# #              color = as.character(group))) +
# #   geom_segment(linewidth = 15, lineend = "butt")#+
# #   #scale_x_continuous(breaks = seq(0, 40, 5)) +
# #   #scale_y_continuous(NULL, breaks = NULL) +
# #   #theme_minimal() +
# #   #labs(x = "Position", color = "Value")
# #   
# # 
# # 
# # 
# # 
# # 
# # true_node_anc_combine
# # 
# # 
# # df <- data.frame(
# #   start = c(1, 2, 5, 5.5, 1, 3, 6, 8),
# #   end   = c(2, 5, 5.5, 10, 3, 6, 8, 10),
# #   group = c("a", "b", "a", "c", "d", "e", "d", "f"),
# #   broad_group = c(1, 1, 1, 1, 2, 2, 2, 2)
# # )
# # 
# # ggplot(df) +
# #   geom_segment(
# #     aes(
# #       x = start, xend = end,
# #       y = broad_group, yend = broad_group,
# #       color = group
# #     ),
# #     linewidth = 10,
# #     lineend = "butt"
# #   ) +
# #   scale_y_continuous(
# #     breaks = c(1, 2),
# #     labels = c("Group 1", "Group 2")
# #   ) +
# #   scale_x_continuous(breaks = 0:10) +
# #   labs(
# #     x = "Length",
# #     y = NULL,
# #     color = "Segment group"
# #   ) +
# #   theme_minimal()
# # 
# # 
# # 
# # 
# # k3_max_lik_group_membership_list_relabel_list_df
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # k3_max_lik_group_membership_list_relabel_list$ind151
# # 
# # 
# # 
# # true_location_matchup_list$ind151
# # 
# # k3_max_lik_group_membership_list$ind151[[1]][1:10]
# # k3_max_lik_group_membership_list_relabel_list$ind151[1:10]
# # 
# # 
# # 
# # k3_max_lik_group_membership_list_relabel
# # 
# # 
# # #######
# # #######
# # #######
# # 
# # MISSING_INDIV <- admix_indivs[3]
# # 
# # index <- which(admix_indivs == MISSING_INDIV)
# # 
# # MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
# # REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
# # 
# # 
# # miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
# # rownames(miss_info) <- miss_info$n_node
# # miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# # 
# # tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# # 
# # k3_process <- process_trees(
# #   tree = tree_list_subset, 
# #   trait = miss_info_mat, 
# #   brown_rate_mat = root_rate_info$mean_brown_rate, 
# #   root_coords_list = root_rate_info$root_list,
# #   missing_inds = MISSING_NODES
# # )
# # 
# # 
# # cond_mean_list <- lapply(k3_process, function(x) x$cond_distr$condmean)
# # cond_var_list <- lapply(k3_process, function(x) x$cond_distr$condvar)
# # 
# # test_hard_assign <- em_hardassign_wrapper(times = 5,
# #                                                                      cond_mean = cond_mean_list, 
# #                                                                      cond_var = cond_var_list, 
# #                                                                      k = 3, 
# #                                                                      max_stp = 25, 
# #                                                                      conv_thresh = 1e-5, 
# #                                                                      sp_bounds = list(x = c(-5, 105), 
# #                                                                                       y = c(-5, 105)),
# #                                                                      starting_location_seed = NULL)
# # 
# # test_hard_assign$max_lik_xmat_list
# # test_hard_assign$maxlik_indices
# # 
# # cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
# # 
# # plot_color <- color_vec[index]
# # 
# # ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
# #                                        admix_node_df = admix_node_info, 
# #                                        offspring_df = offspring_info, 
# #                                        parent_loc_df = parent_locs)
# # 
# # ancestor_locs_summarize <- ancestor_locs %>% 
# #   mutate(location_id_plus = location_id + 1) %>% 
# #   group_by(location_id_plus) %>%
# #   summarize(x = mean(x),
# #             y = mean(y),
# #             .groups = 'drop')
# # 
# # 
# # test <- ranked_match_up_multi(list(as.matrix(ancestor_locs_summarize[,c('x', 'y')]),
# #                            test_hard_assign$max_lik_xmat_list[[1]][,c(1,2)])
# #                       )
# # 
# # convert_labels <- function(old_label, label_map_original, label_map_new) {
# #   
# #   new_label <- setNames(rep(NA, length(old_label)), nm = names(old_label))
# #   
# #   for (i in seq_along(label_map_original)) {
# #     original_label <- label_map_original[i]
# #     updated_label <- label_map_new[i]
# #     new_label[old_label == original_label] <- updated_label
# #   }
# #   
# #   return(new_label)
# # }
# # 
# # convert_labels(old_label = test_hard_assign$max_lik_group_membership_list[[1]], 
# #                label_map_original = test[[2]]$original_group, 
# #                label_map_new = test[[2]]$ref_group)
# # 
# # 
# # test_hard_assign$max_lik_group_membership_list[[1]]
# # 
# # new_label <- rep(NA, length(test_hard_assign$max_lik_group_membership_list[[1]]))
# # names(new_label) <- names(test_hard_assign$max_lik_group_membership_list[[1]])
# # 
# # for (i in seq_len(nrow(test[[2]]))) {
# #   original_label <- test[[2]]$original_group[i]
# #   updated_label <- test[[2]]$ref_group[i]
# #   new_label[test_hard_assign$max_lik_group_membership_list[[1]] == original_label] <- updated_label
# # }
# # 
# # test
# # test_hard_assign$max_lik_group_membership_list[[1]][1:12]
# # new_label[1:12]
# # 
# # 
# # test_hard_assign$max_lik_group_membership_list[[1]]
# # 
# # new_label
# # 
# # 
# # ancestor_locs_summarize
# # 
# # test_hard_assign$max_lik_xmat_list
# # 
# # 
# # 
# # matchup_df <- ranked_match_up(ref_location = as.matrix(ancestor_locs_summarize[,c(2,3)]), 
# #                 compare_location = matrix(k3_est_list$IND152$top_x, ncol = 2, byrow = TRUE)
# #                 )
# # 
# # true_locs <- as.data.frame(ancestor_locs_summarize[,c(2,3)])
# # true_locs[matchup_df$ref_group,'group'] <- 1:3
# # true_locs[,'value'] <- 'ref'
# # 
# # 
# # est_locs <- as.data.frame(matrix(k3_est_list$IND152$top_x, ncol = 2, byrow = TRUE))
# # colnames(est_locs) <- c('x', 'y')
# # est_locs[matchup_df$compare_group,'group'] <- 1:3
# # est_locs[,'value'] <- 'compare'
# # 
# # 
# # test_emhardassign <- em_hardassign_wrapper(times = 10,
# #                       cond_mean = cond_mean_list[1:400], 
# #                       cond_var = cond_var_list[1:400], 
# #                       k = 4, 
# #                       max_stp = 25, 
# #                       conv_thresh = 1e-5, 
# #                       sp_bounds = list(x = c(-5, 105), 
# #                                        y = c(-5, 105)),
# #                       starting_location_seed = NULL)
# # 
# # 
# # 
# # test_emhardassign$max_lik_xmat_list
# # 
# # test_emsoftassign <- em_softassign_wrapper(times = 5,
# #                                            cond_mean = cond_mean_list[1:400], 
# #                                            cond_var = cond_var_list[1:400], 
# #                                            k = 3, 
# #                                            max_stp = 10, 
# #                                            conv_thresh = 1e-4, 
# #                                            sp_bounds = list(x = c(-5, 105), 
# #                                                             y = c(-5, 105)),
# #                                            starting_location_seed = NULL)
# # 
# # test_emsoftassign$max_lik_xmat_list
# # 
# # 
# # test_emsoftassign$max_lik_xmat_list
# # 
# # test_emhardassign$max_lik_xmat_list
# # 
# # test_hard_assign <- em_hardassign(cond_mean = cond_mean_list, 
# #                       cond_var = cond_var_list, 
# #                       k = 3, 
# #                       max_stp = 25, 
# #                       conv_thresh = 1e-5, 
# #                       sp_bounds = list(x = c(-5, 105), 
# #                                        y = c(-5, 105)),
# #                       starting_location_seed = NULL)
# # 
# # test_soft_assign <- em_softassign(cond_mean = cond_mean_list, 
# #               cond_var = cond_var_list, 
# #               k = 3, 
# #               max_stp = 25, 
# #               conv_thresh = 1e-5, 
# #               sp_bounds = list(x = c(-5, 105), 
# #                                y = c(-5, 105)),
# #               starting_location_seed = NULL)
# # 
# # 
# # 
# # em_hardassign(cond_mean = cond_mean_list, 
# #                       cond_var = cond_var_list, 
# #                       k = 3, 
# #                       max_stp = 25, 
# #                       conv_thresh = 1e-8, 
# #                       sp_bounds = list(x = c(-5, 105), 
# #                                        y = c(-5, 105)),
# #                       starting_location_seed = NULL)
# # 
# # 
# # 
# # 
# # k3_est_list$IND152$output_list[[1]]$final_x
# # 
# # 
# # k3_est_list$IND150$output_list
# # 
# # 
# # 
# # 
# # 
# # k3_est_list[[paste0('IND', MISSING_INDIV)]]
# # 
# # 
# # 
# # k3_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = plot_color, size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = k3_est_list[[paste0('IND', MISSING_INDIV)]]$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
# #   geom_point(data = ancestor_locs,
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = plot_color, stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   geom_point(data = rbind(matchup_df$compare_group_reorder %>% mutate(type = 'compare'), 
# #                           matchup_df$ref_group_reorder %>% mutate(type = 'ref')),
# #              aes(x,y, color = as.character(ref_group), shape = type)) #+
# #   #geom_line(data = rbind(matchup_df$compare_group_reorder, 
# #   #                       matchup_df$ref_group_reorder),
# #   #          aes(x = x, y = y, group = ref_order), color = "gray50") # Connects points sharing
# # 
# # message('FINISHED ', MISSING_INDIV)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ancestor_locs
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# #   
# #   
# # 
# # ### COMBINING PLOTS ###
# # infer_examples_multipan <- cowplot::plot_grid(k2_multipan_withtitle, k3_multipan_withtitle, nrow = 2)
# # 
# # ggsave(infer_examples_multipan,
# #        filename = here('plots', 'tmp', 'infer_examples_multipan.png'),
# #        width = 5*1.2, height = 5*1.2,
# #        bg = 'white')
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # tree_interval_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq',
# #                                    base_file_name,
# #                                    'tree_interval_info.txt'), header = TRUE)
# # 
# # 
# # node_founder_ancestry_intervals <- read.delim(here('simulation_output',
# #                                       'processed_treeseq',
# #                                       base_file_name,
# #                                       'node_founder_ancestry_intervals.txt'), header = TRUE)
# # 
# # 
# # parent_locs
# # focal_node <- 300
# # 
# # focal_node
# # 
# # tree_interval_info
# # 
# # 
# # 
# # 
# # 
# # location_summary <- parent_locs %>% 
# #   group_by(subpop, location_id) %>% 
# #   summarize(x_mean = mean(x),
# #             y_mean = mean(y),
# #             .groups = 'drop')
# # 
# # test_assign <- assign_ancestor(tree_interval_info = tree_interval_info, 
# #                                ancestry_interval_info = node_founder_ancestry_intervals,
# #                                ancestor_locations = parent_locs,
# #                                ancestor_location_join_col = 'pedigree_id')
# # test_assign_combine <- do.call(rbind, test_assign) %>% 
# #   left_join(., admix_node_info, by = c("focal_node" = "node"))
# # 
# # ancestry_prop <- test_assign_combine %>% 
# #   group_by(indiv, ancestor_id) %>% 
# #   summarize(count = n()) %>% 
# #   ungroup() %>% 
# #   group_by(indiv) %>% 
# #   mutate(prop = count/sum(count)) %>% 
# #   ungroup()
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # tree_interval_info_copy <- tree_interval_info
# # 
# # 
# # left_index <- max(which(anc_interv_subset$left_interval <= tree_interval_subset$left_interval))
# # right_index <- min(which(anc_interv_subset$right_interval >= tree_interval_subset$right_interval))
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # location_summary <- parent_locs %>% 
# #   group_by(subpop, location_id) %>% 
# #   summarize(x_mean = mean(x),
# #             y_mean = mean(y),
# #             .groups = 'drop')
# # 
# # 
# # location_pop2 <- location_summary[location_summary$subpop == 2,]
# # 
# # 
# # 
# # test_emhardassign <- em_hardassign_wrapper(times = 20,
# #                                            cond_mean = cond_mean_list, 
# #                                            cond_var = cond_var_list, 
# #                                            k = 4, 
# #                                            max_stp = 25, 
# #                                            conv_thresh = 1e-6, 
# #                                            sp_bounds = list(x = c(-5, 105), 
# #                                                             y = c(-5, 105)),
# #                                            starting_location_seed = NULL)
# # 
# # 
# # process_em_output(
# #   test_emhardassign$output_list
# # )
# # 
# # 
# # test_emhardassign <- em_softassign_wrapper(5,
# #                                            cond_mean = cond_mean_list, 
# #                                            cond_var = cond_var_list, 
# #                                            k = 3, 
# #                                            max_stp = 25, 
# #                                            conv_thresh = 1e-5, 
# #                                            sp_bounds = list(x = c(-5, 105), 
# #                                                             y = c(-5, 105)),
# #                                            starting_location_seed = NULL)
# # 
# # 
# # k3_est_list$IND152
# # 
# # 
# # k3_est_list[['IND152']]$output_list[[2]]$group_membership_list
# # 
# # test_assign$node304
# # 
# # node_membership_305 <- test_emhardassign$max_lik_group_membership_list[[1]][grepl(pattern = 'n304', names(test_emhardassign$max_lik_group_membership_list[[1]]))]
# # 
# # 
# # 
# # 
# # test_assign$node304 %>% 
# #   mutate(right = cumsum(span),
# #          left = lag(right, default = 0)) %>% 
# #   mutate(inferred_vals = node_membership_305) %>% 
# #   ggplot(aes(x = left, xend = right,
# #                  y = 1, yend = 1,
# #                  color = as.character(location_id))) +
# #   geom_segment(linewidth = 15, lineend = "butt") +
# #   #scale_x_continuous(breaks = seq(0, 40, 5)) +
# #   scale_y_continuous(NULL, breaks = NULL) +
# #   theme_minimal() +
# #   labs(x = "Position", color = "Value")
# # 
# # test_assign$node305 %>% 
# #   mutate(right = cumsum(span),
# #          left = lag(right, default = 0)) %>% 
# #   mutate(inferred_vals = node_membership_305) %>% 
# #   ggplot(aes(x = left, xend = right,
# #              y = 1, yend = 1,
# #              color = as.character(location_id))) +
# #   geom_segment(linewidth = 15, lineend = "butt") +
# #   #scale_x_continuous(breaks = seq(0, 40, 5)) +
# #   scale_y_continuous(NULL, breaks = NULL) +
# #   theme_minimal() +
# #   labs(x = "Position", color = "Value")
# # 
# # 
# # test_assign$node304 %>% 
# #   mutate(right = cumsum(span),
# #          left = lag(right, default = 0)) %>% 
# #   mutate(inferred_vals = node_membership_305) %>% 
# #   ggplot(aes(x = left, xend = right,
# #              y = 1, yend = 1,
# #              color = as.character(inferred_vals))) +
# #   geom_segment(linewidth = 15, lineend = "butt") +
# #   #scale_x_continuous(breaks = seq(0, 40, 5)) +
# #   scale_y_continuous(NULL, breaks = NULL) +
# #   theme_minimal() +
# #   labs(x = "Position", color = "Value")
# # 
# # 
# # 
# # 
# # test_assign$node304 %>% 
# #   mutate(right = cumsum(span),
# #          left = lag(right, default = 0)) %>% 
# #   ggplot(aes(x = left_interval, xend = right_interval,
# #              y = 1, yend = 1,
# #              color = as.character(location_id))) +
# #   geom_segment(data = node_founder_ancestry_intervals,
# #                aes(x = left_interval, xend = right_interval,
# #                    y = 1, yend = 1,
# #                    color = as.character(location_id))) +
# #   geom_segment(linewidth = 15, lineend = "butt") +
# #   #scale_x_continuous(breaks = seq(0, 40, 5)) +
# #   scale_y_continuous(NULL, breaks = NULL) +
# #   theme_minimal() +
# #   labs(x = "Position", color = "Value")
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # test_assign <- assign_ancestor(tree_interval_info = tree_interval_info, 
# #                 ancestry_interval_info = node_founder_ancestry_intervals,
# #                 ancestor_locations = parent_locs,
# #                 ancestor_location_join_col = 'pedigree_id')
# # test_assign_combine <- do.call(rbind, test_assign) %>% 
# #   left_join(., admix_node_info, by = c("focal_node" = "node"))
# # 
# # ancestry_prop <- test_assign_combine %>% 
# #   group_by(indiv, ancestor_id) %>% 
# #   summarize(count = n()) %>% 
# #   ungroup() %>% 
# #   group_by(indiv) %>% 
# #   mutate(prop = count/sum(count)) %>% 
# #   ungroup()
# # 
# # 
# # 
# # 
# # 
# # 
# # left_join(test_assign[[1]], parent_locs, by = c('ancestor_id' = 'pedigree_id'))
# # 
# # all(tree_interval_info_copy$`300` %in% parent_locs$pedigree_id)
# # 
# # 
# # tree_interval_info_copy[,paste0(focal_node)] <- NA
# # for (i in seq_len(nrow(tree_interval_info_copy))) {
# #   left_index <- max(which(anc_interv_subset$left_interval <= tree_interval_info_copy[i,]$left_interval))
# #   right_index <- min(which(anc_interv_subset$right_interval >= tree_interval_info_copy[i,]$right_interval))
# #   
# #   if (left_index == right_index) {
# #     tree_interval_info_copy[i,paste0(focal_node)] <- anc_interv_subset[left_index,]$ancestor
# #   } else {
# #     tree_interval_info_copy[i,paste0(focal_node)] <- 'multianc'
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # true_locs <- data.frame(
# #   group = 1:4,
# #   x = c(10, 20, 40, 5),
# #   y = c(11, 22, 38, 5)
# # )
# # 
# # 
# # estimate_locs <- data.frame(
# #   x = c(18.4, 23.1, 50.0, 9.1),
# #   y = c(12.2, 18.2, 25.9, 0.4)
# # )
# # 
# # estimate_locs <- estimate_locs[sample(nrow(estimate_locs)),]
# # estimate_locs$group <- 1:4
# # 
# # 
# # 
# # dist_mat <- matrix(data = NA, nrow = nrow(true_locs), ncol = nrow(estimate_locs))
# # 
# # for (REAL in 1:nrow(true_locs)) {
# #   for (EST in 1:nrow(estimate_locs)) {
# #     dist_mat[REAL,EST] <- dist(rbind(true_locs[REAL,], estimate_locs[EST,]))
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # matchup_df <- data.frame(ref_group = rep(NA, ncol(dist_mat)),
# #                          compare_group = rep(NA, ncol(dist_mat))
# #                          )
# # 
# # total_length <- nrow(dist_mat)
# # row_ind <- 1:total_length
# # col_ind <- 1:total_length
# # 
# # for (i in 1:total_length) {
# #   #find the smallest distance
# #   #i <- 2
# #   min_coord <- which(dist_mat == min(dist_mat), arr.ind = TRUE)
# #   
# #   ref_ind <- row_ind[min_coord[1]]
# #   compare_ind <- col_ind[min_coord[2]]
# #   
# #   matchup_df[i,] <- c(ref_ind, compare_ind)
# #   
# #   
# #   if (i == total_length) break
# #     
# #   #remove those inds from consideration
# #   dist_mat <- dist_mat[-min_coord[1],-min_coord[2],drop=FALSE]
# #   row_ind <- row_ind[-min_coord[1]]
# #   col_ind <- col_ind[-min_coord[2]]
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # #par(mfrow = c(1,1))
# # plot(true_locs$x, true_locs$y, col = 'blue', xlim = c(0, 52), ylim = c(0, 52))
# # points(estimate_locs$x, estimate_locs$y, col = 'gray')
# # text(true_locs$x, true_locs$y, labels = true_locs$group, col = 'blue', pos = 3, cex = 0.8)
# # text(estimate_locs$x, estimate_locs$y, labels = estimate_locs$group, col = 'gray', pos = 3, cex = 0.8)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ##########################################################################
# # ### INITIAL EXPLORATION OF REMOVING INTRA-TREE COVAR AND USING WRONG K ###
# # ##########################################################################
# # 
# # ### LOADING TREES AND OTHER INFO ###
# # base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME15000_DISP1.75_CC0.20_K2_Kreps5'
# # 
# # files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# # tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
# #                     function(x) ape::read.tree(x)) 
# # 
# # indiv_node_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq', 
# #                                    base_file_name,
# #                                    paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
# # indiv_node_info$n_node <- paste0('n', indiv_node_info$node)
# # 
# # 
# # 
# # offspring_info <- read.delim(here('simulation_output', 
# #                                   'postburnin_sims',
# #                                   paste0('offspring_info_', base_file_name, '.txt')), header = TRUE, sep = " ")
# # 
# # parent_locs <- read.delim(here('simulation_output', 
# #                                'postburnin_sims',
# #                                paste0('parents_locs_ped_', base_file_name, '.txt')), header = TRUE, sep = " ")
# # 
# # admix_node_info <- read.delim(here('simulation_output',
# #                                    'processed_treeseq',
# #                                    base_file_name,
# #                                    paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)
# # 
# # 
# # files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
# # tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
# #                     function(x) ape::read.tree(x)) 
# # 
# # 
# # 
# # ### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
# # admix_indivs <- unique(admix_node_info$indiv)
# # admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]
# # 
# # miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
# # rownames(miss_info) <- miss_info$n_node
# # miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# # 
# # #tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# # tree_list_subset1 <- tree_list[1:300]
# # root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
# #                                       miss_info_mat, 
# #                                       admix_nodes)
# # 
# # 
# # ### INFERRING LOCATION FOR A SINGLE ADMIXED INDIVIDUAL ###
# # MISSING_INDIV <- admix_indivs[5]
# # MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
# # REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
# # 
# # 
# # miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
# # rownames(miss_info) <- miss_info$n_node
# # miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# # 
# # tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
# # 
# # k2_process <- process_trees(
# #   tree = tree_list_subset, 
# #   trait = miss_info_mat, 
# #   brown_rate_mat = root_rate_info$mean_brown_rate, 
# #   root_coords_list = root_rate_info$root_list,
# #   missing_inds = MISSING_NODES
# # )
# # 
# # 
# # cond_mean_list <- lapply(k2_process, function(x) x$cond_distr$condmean)
# # cond_var_list <- lapply(k2_process, function(x) x$cond_distr$condvar)
# # cond_var_list_alt <- lapply(cond_var_list, function(x) diag(diag(x)) ) #remove covariance
# # 
# # 
# # ### LOCATION ESTIMATION ###
# # hardassign_k3 <- em_hardassign_wrapper(times = 5,
# #                                     cond_mean = cond_mean_list, 
# #                                cond_var = cond_var_list, 
# #                                k = 3, 
# #                                max_stp = 30, 
# #                                conv_thresh = 1e-5, 
# #                                sp_bounds = list(x = c(-5, 105), 
# #                                                 y = c(-5, 105)),
# #                                starting_location_seed = NULL
# # )
# # 
# # softassign_k3 <- em_softassign_wrapper(times = 2,
# #                                     cond_mean = cond_mean_list, 
# #                                cond_var = cond_var_list, 
# #                                k = 3, 
# #                                max_stp = 30, 
# #                                conv_thresh = 1e-5, 
# #                                sp_bounds = list(x = c(-5, 105), 
# #                                                 y = c(-5, 105)),
# #                                starting_location_seed = NULL
# # )
# # 
# # hardassign_k2 <- em_hardassign_wrapper(times = 5,
# #                                cond_mean = cond_mean_list, 
# #                                cond_var = cond_var_list, 
# #                                k = 2, 
# #                                max_stp = 60, 
# #                                conv_thresh = 1e-5, 
# #                                sp_bounds = list(x = c(-5, 105), 
# #                                                 y = c(-5, 105)),
# #                                starting_location_seed = NULL
# # )
# # 
# # softassign_k2 <- em_softassign_wrapper(times = 2,
# #                                     cond_mean = cond_mean_list, 
# #                                cond_var = cond_var_list, 
# #                                k = 2, 
# #                                max_stp = 30, 
# #                                conv_thresh = 1e-5, 
# #                                sp_bounds = list(x = c(-5, 105), 
# #                                                 y = c(-5, 105)),
# #                                starting_location_seed = NULL
# # )
# # 
# # softassign_k2_alt <- em_softassign_wrapper(times = 2,
# #                                         cond_mean = cond_mean_list, 
# #                                    cond_var = cond_var_list_alt, 
# #                                    k = 2, 
# #                                    max_stp = 30, 
# #                                    conv_thresh = 1e-5, 
# #                                    sp_bounds = list(x = c(-5, 105), 
# #                                                     y = c(-5, 105)),
# #                                    starting_location_seed = NULL
# # )
# # 
# # softassign_k3_alt <- em_softassign_wrapper(times = 2,
# #                                         cond_mean = cond_mean_list, 
# #                                    cond_var = cond_var_list_alt, 
# #                                    k = 3, 
# #                                    max_stp = 30, 
# #                                    conv_thresh = 1e-5, 
# #                                    sp_bounds = list(x = c(-5, 105), 
# #                                                     y = c(-5, 105)),
# #                                    starting_location_seed = NULL
# # )
# # 
# # 
# # 
# # ### PLOTTING LOCATION ESTIMATES ###
# # cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
# # 
# # 
# # k3_soft_plot <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = '#f48373', size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = softassign_k3$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
# #   geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   ggtitle('Soft assign with K = 3')
# # 
# # 
# # k3_hard_plot <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = '#f48373', size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = hardassign_k3$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
# #   geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   ggtitle('Hard assign with K = 3')
# # 
# # 
# # k3_soft_alt_plot <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = '#f48373', size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = softassign_k3_alt$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
# #   geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   ggtitle('Soft assign with K = 3 (ignore intra-tree covar)')
# # 
# # 
# # k2_soft_plot <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = '#f48373', size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = softassign_k2$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
# #   geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   ggtitle('Soft assign with K = 2')
# # 
# # k2_hard_plot <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = '#f48373', size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = hardassign_k2$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
# #   geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   ggtitle('Hard assign with K = 2')
# # 
# # k2_soft_alt_plot <- ggplot() +
# #   annotate("rect", 
# #            xmin = 0, xmax = 100,       # Start and end on X axis
# #            ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
# #            fill = "#e5e5e5", #"#fac9c3",       # Background color
# #            alpha = 0.4) +
# #   geom_point(data = cond_mean_df,
# #              aes(x = x, y = y), shape = 8,
# #              color = '#f48373', size = 0.75, alpha = 0.4) +
# #   geom_point(data = as.data.frame(miss_info_mat),
# #              aes(x = x, y = y),
# #              color = '#999999', size = 0.7, alpha = 0.3) +
# #   geom_point(data = data.frame(matrix(data = softassign_k2_alt$top_x, ncol = 2, byrow = TRUE)),
# #              aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
# #   geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
# #              aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
# #   theme_bw() +
# #   theme(panel.grid.minor = element_blank(),
# #         #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
# #         panel.grid.major = element_blank(),
# #         panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
# #         plot.title = element_text(hjust = 0.5, size = 15)) +
# #   ggtitle('Soft assign with K = 2 (ignore intra-tree covar)')
# # 
# # 
# # #combine into multipanel and export
# # k_assignment_type_explore <- cowplot::plot_grid(
# #   k2_hard_plot, k2_soft_plot, k2_soft_alt_plot,
# #   k3_hard_plot, k3_soft_plot, k3_soft_alt_plot,
# #   nrow = 2)
# # 
# # ggsave(k_assignment_type_explore,
# #        filename = here('plots', 'tmp', 'k_assignment_type_explore.png'),
# #        width = 10*1.5, height = 6*1.5,
# #        bg = 'white')
