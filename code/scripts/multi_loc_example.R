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

source(here('code', 'scripts', 'location_est_code.R'))
source(here('code', 'scripts', 'custom_project_funcs.R'))



#################################
### PROCESSING INFO FOR K = 2 ###
#################################

### LOADING TREES AND OTHER INFO ###
base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME15000_DISP1.75_CC0.20_K2_Kreps5'

files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                         function(x) ape::read.tree(x)) 

indiv_node_info <- read.delim(here('simulation_output',
                                   'processed_treeseq', 
                                   base_file_name,
                                   paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
indiv_node_info$n_node <- paste0('n', indiv_node_info$node)



offspring_info <- read.delim(here('simulation_output', 
                                  'postburnin_sims',
                                  paste0('offspring_info_', base_file_name, '.txt')), header = TRUE, sep = " ")

parent_locs <- read.delim(here('simulation_output', 
                               'postburnin_sims',
                               paste0('parents_locs_ped_', base_file_name, '.txt')), header = TRUE, sep = " ")

admix_node_info <- read.delim(here('simulation_output',
                                   'processed_treeseq',
                                   base_file_name,
                                   paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)


files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 



### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
admix_indivs <- unique(admix_node_info$indiv)
admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]

miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
rownames(miss_info) <- miss_info$n_node
miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])

#tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
tree_list_subset1 <- tree_list[1:500]
root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
                                      miss_info_mat, 
                                      admix_nodes)


### LOCATION ESTIMATION AND PLOTTING ###

k2_est_list <- list()
k2_plot_list <- list()

color_vec <- c("#8cb369", "#f4e285", "#f4a259", "#0096c7", "#bc4b51")

for (MISSING_INDIV in admix_indivs) {
  index <- which(admix_indivs == MISSING_INDIV)
    
  MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
  REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
  
  
  miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
  rownames(miss_info) <- miss_info$n_node
  miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
  
  tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
  
  k2_process <- process_trees(
    tree = tree_list_subset, 
    trait = miss_info_mat, 
    brown_rate_mat = root_rate_info$mean_brown_rate, 
    root_coords_list = root_rate_info$root_list,
    missing_inds = MISSING_NODES
  )
  
  
  cond_mean_list <- lapply(k2_process, function(x) x$cond_distr$condmean)
  cond_var_list <- lapply(k2_process, function(x) x$cond_distr$condvar)
  
  k2_est_list[[paste0('IND', MISSING_INDIV)]] <- em_softassign_wrapper(times = 2,
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
                        admix_node_df = admix_node_info, 
                        offspring_df = offspring_info, 
                        parent_loc_df = parent_locs)
  
  k2_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
    annotate("rect", 
             xmin = 0, xmax = 100,       # Start and end on X axis
             ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
             fill = "#e5e5e5", #"#fac9c3",       # Background color
             alpha = 0.4) +
    geom_point(data = cond_mean_df,
               aes(x = x, y = y), shape = 8,
               color = plot_color, size = 0.75, alpha = 0.4) +
    geom_point(data = as.data.frame(miss_info_mat),
               aes(x = x, y = y),
               color = '#999999', size = 0.7, alpha = 0.3) +
    geom_point(data = data.frame(matrix(data = k2_est_list[[paste0('IND', MISSING_INDIV)]]$top_x, ncol = 2, byrow = TRUE)),
               aes(X1, X2), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
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


k2_plot_title <- title <- ggdraw() + draw_label("K2 (soft node assignment)", fontface = 'bold', x = 0.5, hjust = 0.5)
k2_multipan <- cowplot::plot_grid(plotlist = k2_plot_list, nrow = 1)

k2_multipan_withtitle <- plot_grid(k2_plot_title, k2_multipan, ncol = 1, rel_heights = c(0.1, 1))



#################################
### PROCESSING INFO FOR K = 3 ###
#################################

### LOADING TREES AND OTHER INFO ###
base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME15000_DISP1.75_CC0.20_K3_Kreps5'

files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 

indiv_node_info <- read.delim(here('simulation_output',
                                   'processed_treeseq', 
                                   base_file_name,
                                   paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
indiv_node_info$n_node <- paste0('n', indiv_node_info$node)



offspring_info <- read.delim(here('simulation_output', 
                                  'postburnin_sims',
                                  paste0('offspring_info_', base_file_name, '.txt')), header = TRUE, sep = " ")

parent_locs <- read.delim(here('simulation_output', 
                               'postburnin_sims',
                               paste0('parents_locs_ped_', base_file_name, '.txt')), header = TRUE, sep = " ")

admix_node_info <- read.delim(here('simulation_output',
                                   'processed_treeseq',
                                   base_file_name,
                                   paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)


files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 



### ESTIMATION ###
admix_indivs <- unique(admix_node_info$indiv)
admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]

miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
rownames(miss_info) <- miss_info$n_node
miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])

#tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
tree_list_subset1 <- tree_list
root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
                                      miss_info_mat, 
                                      admix_nodes)


### LOCATION ESTIMATION AND PLOTTING ###

k3_est_list <- list()
k3_plot_list <- list()

color_vec <- c("#b266b2", "#a78a7f", "#c8b6ff", "#4abcde", "#a69cac")

for (MISSING_INDIV in admix_indivs) {
  index <- which(admix_indivs == MISSING_INDIV)
  
  MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
  REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]
  
  
  miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
  rownames(miss_info) <- miss_info$n_node
  miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
  
  tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
  
  k3_process <- process_trees(
    tree = tree_list_subset, 
    trait = miss_info_mat, 
    brown_rate_mat = root_rate_info$mean_brown_rate, 
    root_coords_list = root_rate_info$root_list,
    missing_inds = MISSING_NODES
  )
  
  
  cond_mean_list <- lapply(k3_process, function(x) x$cond_distr$condmean)
  cond_var_list <- lapply(k3_process, function(x) x$cond_distr$condvar)
  
  k3_est_list[[paste0('IND', MISSING_INDIV)]] <- em_hardassign_wrapper(times = 5,
                                                                       cond_mean = cond_mean_list, 
                                                                       cond_var = cond_var_list, 
                                                                       k = 3, 
                                                                       max_stp = 15, 
                                                                       conv_thresh = 1e-5, 
                                                                       sp_bounds = list(x = c(-5, 105), 
                                                                                        y = c(-5, 105)),
                                                                       starting_location_seed = NULL)
  
  cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))
  
  plot_color <- color_vec[index]
  
  ancestor_locs <- pedigree_founder_info(indiv_id = MISSING_INDIV, 
                                         admix_node_df = admix_node_info, 
                                         offspring_df = offspring_info, 
                                         parent_loc_df = parent_locs)
  
  
  k3_plot_list[[paste0('IND', MISSING_INDIV)]] <- ggplot() +
    annotate("rect", 
             xmin = 0, xmax = 100,       # Start and end on X axis
             ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
             fill = "#e5e5e5", #"#fac9c3",       # Background color
             alpha = 0.4) +
    geom_point(data = cond_mean_df,
               aes(x = x, y = y), shape = 8,
               color = plot_color, size = 0.75, alpha = 0.4) +
    geom_point(data = as.data.frame(miss_info_mat),
               aes(x = x, y = y),
               color = '#999999', size = 0.7, alpha = 0.3) +
    geom_point(data = data.frame(matrix(data = k3_est_list[[paste0('IND', MISSING_INDIV)]]$top_x, ncol = 2, byrow = TRUE)),
               aes(X1, X2), shape = 21, fill = plot_color, color = 'white', size = 5, alpha = 1) +
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

k3_plot_title <- title <- ggdraw() + draw_label("K3 (hard node assignment)", fontface = 'bold', x = 0.5, hjust = 0.5)
k3_multipan <- cowplot::plot_grid(plotlist = k3_plot_list, nrow = 1)

k3_multipan_withtitle <- plot_grid(k3_plot_title, k3_multipan, ncol = 1, rel_heights = c(0.1, 1))



### COMBINING PLOTS ###
infer_examples_multipan <- cowplot::plot_grid(k2_multipan_withtitle, k3_multipan_withtitle, nrow = 2)

ggsave(infer_examples_multipan,
       filename = here('plots', 'tmp', 'infer_examples_multipan.png'),
       width = 12*1.4, height = 5*1.4,
       bg = 'white')



##########################################################################
### INITIAL EXPLORATION OF REMOVING INTRA-TREE COVAR AND USING WRONG K ###
##########################################################################

### LOADING TREES AND OTHER INFO ###
base_file_name <- 'POSTBURNIN_BURNINREP1_RUNTIME15000_DISP1.75_CC0.20_K2_Kreps5'

files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 

indiv_node_info <- read.delim(here('simulation_output',
                                   'processed_treeseq', 
                                   base_file_name,
                                   paste0(base_file_name, '_indiv_node_info.txt')), header = TRUE)
indiv_node_info$n_node <- paste0('n', indiv_node_info$node)



offspring_info <- read.delim(here('simulation_output', 
                                  'postburnin_sims',
                                  paste0('offspring_info_', base_file_name, '.txt')), header = TRUE, sep = " ")

parent_locs <- read.delim(here('simulation_output', 
                               'postburnin_sims',
                               paste0('parents_locs_ped_', base_file_name, '.txt')), header = TRUE, sep = " ")

admix_node_info <- read.delim(here('simulation_output',
                                   'processed_treeseq',
                                   base_file_name,
                                   paste0(base_file_name, '_admix_node_info.txt')), header = TRUE)


files_vec <- list.files(here('simulation_output', 'processed_treeseq', base_file_name), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 



### ESTIMATION OF ROOT COORDINATES AND BROWNIAN RATE MATRIX ###
admix_indivs <- unique(admix_node_info$indiv)
admix_nodes <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_indivs]

miss_info <- indiv_node_info[!indiv_node_info$n_node %in% admix_nodes,]
rownames(miss_info) <- miss_info$n_node
miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])

#tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
tree_list_subset1 <- tree_list[1:300]
root_rate_info <- root_rate_multitree(tree = tree_list_subset1, 
                                      miss_info_mat, 
                                      admix_nodes)


### INFERRING LOCATION FOR A SINGLE ADMIXED INDIVIDUAL ###
MISSING_INDIV <- admix_indivs[5]
MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% MISSING_INDIV]#[i]#[2]
REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & !indiv_node_info$n_node %in% MISSING_NODES]


miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, MISSING_NODES),]
rownames(miss_info) <- miss_info$n_node
miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])

tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))

k2_process <- process_trees(
  tree = tree_list_subset, 
  trait = miss_info_mat, 
  brown_rate_mat = root_rate_info$mean_brown_rate, 
  root_coords_list = root_rate_info$root_list,
  missing_inds = MISSING_NODES
)


cond_mean_list <- lapply(k2_process, function(x) x$cond_distr$condmean)
cond_var_list <- lapply(k2_process, function(x) x$cond_distr$condvar)
cond_var_list_alt <- lapply(cond_var_list, function(x) diag(diag(x)) ) #remove covariance


### LOCATION ESTIMATION ###
hardassign_k3 <- em_hardassign_wrapper(times = 5,
                                    cond_mean = cond_mean_list, 
                               cond_var = cond_var_list, 
                               k = 3, 
                               max_stp = 30, 
                               conv_thresh = 1e-5, 
                               sp_bounds = list(x = c(-5, 105), 
                                                y = c(-5, 105)),
                               starting_location_seed = NULL
)

softassign_k3 <- em_softassign_wrapper(times = 2,
                                    cond_mean = cond_mean_list, 
                               cond_var = cond_var_list, 
                               k = 3, 
                               max_stp = 30, 
                               conv_thresh = 1e-5, 
                               sp_bounds = list(x = c(-5, 105), 
                                                y = c(-5, 105)),
                               starting_location_seed = NULL
)

hardassign_k2 <- em_hardassign_wrapper(times = 5,
                               cond_mean = cond_mean_list, 
                               cond_var = cond_var_list, 
                               k = 2, 
                               max_stp = 60, 
                               conv_thresh = 1e-5, 
                               sp_bounds = list(x = c(-5, 105), 
                                                y = c(-5, 105)),
                               starting_location_seed = NULL
)

softassign_k2 <- em_softassign_wrapper(times = 2,
                                    cond_mean = cond_mean_list, 
                               cond_var = cond_var_list, 
                               k = 2, 
                               max_stp = 30, 
                               conv_thresh = 1e-5, 
                               sp_bounds = list(x = c(-5, 105), 
                                                y = c(-5, 105)),
                               starting_location_seed = NULL
)

softassign_k2_alt <- em_softassign_wrapper(times = 2,
                                        cond_mean = cond_mean_list, 
                                   cond_var = cond_var_list_alt, 
                                   k = 2, 
                                   max_stp = 30, 
                                   conv_thresh = 1e-5, 
                                   sp_bounds = list(x = c(-5, 105), 
                                                    y = c(-5, 105)),
                                   starting_location_seed = NULL
)

softassign_k3_alt <- em_softassign_wrapper(times = 2,
                                        cond_mean = cond_mean_list, 
                                   cond_var = cond_var_list_alt, 
                                   k = 3, 
                                   max_stp = 30, 
                                   conv_thresh = 1e-5, 
                                   sp_bounds = list(x = c(-5, 105), 
                                                    y = c(-5, 105)),
                                   starting_location_seed = NULL
)



### PLOTTING LOCATION ESTIMATES ###
cond_mean_df <- do.call(rbind,lapply(cond_mean_list, function(df) data.frame(x = df[c(1,3)], y = df[c(2,4)])))


k3_soft_plot <- ggplot() +
  annotate("rect", 
           xmin = 0, xmax = 100,       # Start and end on X axis
           ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
           fill = "#e5e5e5", #"#fac9c3",       # Background color
           alpha = 0.4) +
  geom_point(data = cond_mean_df,
             aes(x = x, y = y), shape = 8,
             color = '#f48373', size = 0.75, alpha = 0.4) +
  geom_point(data = as.data.frame(miss_info_mat),
             aes(x = x, y = y),
             color = '#999999', size = 0.7, alpha = 0.3) +
  geom_point(data = data.frame(matrix(data = softassign_k3$top_x, ncol = 2, byrow = TRUE)),
             aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
  geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
             aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  ggtitle('Soft assign with K = 3')


k3_hard_plot <- ggplot() +
  annotate("rect", 
           xmin = 0, xmax = 100,       # Start and end on X axis
           ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
           fill = "#e5e5e5", #"#fac9c3",       # Background color
           alpha = 0.4) +
  geom_point(data = cond_mean_df,
             aes(x = x, y = y), shape = 8,
             color = '#f48373', size = 0.75, alpha = 0.4) +
  geom_point(data = as.data.frame(miss_info_mat),
             aes(x = x, y = y),
             color = '#999999', size = 0.7, alpha = 0.3) +
  geom_point(data = data.frame(matrix(data = hardassign_k3$top_x, ncol = 2, byrow = TRUE)),
             aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
  geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
             aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  ggtitle('Hard assign with K = 3')


k3_soft_alt_plot <- ggplot() +
  annotate("rect", 
           xmin = 0, xmax = 100,       # Start and end on X axis
           ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
           fill = "#e5e5e5", #"#fac9c3",       # Background color
           alpha = 0.4) +
  geom_point(data = cond_mean_df,
             aes(x = x, y = y), shape = 8,
             color = '#f48373', size = 0.75, alpha = 0.4) +
  geom_point(data = as.data.frame(miss_info_mat),
             aes(x = x, y = y),
             color = '#999999', size = 0.7, alpha = 0.3) +
  geom_point(data = data.frame(matrix(data = softassign_k3_alt$top_x, ncol = 2, byrow = TRUE)),
             aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
  geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
             aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  ggtitle('Soft assign with K = 3 (ignore intra-tree covar)')


k2_soft_plot <- ggplot() +
  annotate("rect", 
           xmin = 0, xmax = 100,       # Start and end on X axis
           ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
           fill = "#e5e5e5", #"#fac9c3",       # Background color
           alpha = 0.4) +
  geom_point(data = cond_mean_df,
             aes(x = x, y = y), shape = 8,
             color = '#f48373', size = 0.75, alpha = 0.4) +
  geom_point(data = as.data.frame(miss_info_mat),
             aes(x = x, y = y),
             color = '#999999', size = 0.7, alpha = 0.3) +
  geom_point(data = data.frame(matrix(data = softassign_k2$top_x, ncol = 2, byrow = TRUE)),
             aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
  geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
             aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  ggtitle('Soft assign with K = 2')

k2_hard_plot <- ggplot() +
  annotate("rect", 
           xmin = 0, xmax = 100,       # Start and end on X axis
           ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
           fill = "#e5e5e5", #"#fac9c3",       # Background color
           alpha = 0.4) +
  geom_point(data = cond_mean_df,
             aes(x = x, y = y), shape = 8,
             color = '#f48373', size = 0.75, alpha = 0.4) +
  geom_point(data = as.data.frame(miss_info_mat),
             aes(x = x, y = y),
             color = '#999999', size = 0.7, alpha = 0.3) +
  geom_point(data = data.frame(matrix(data = hardassign_k2$top_x, ncol = 2, byrow = TRUE)),
             aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
  geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
             aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  ggtitle('Hard assign with K = 2')

k2_soft_alt_plot <- ggplot() +
  annotate("rect", 
           xmin = 0, xmax = 100,       # Start and end on X axis
           ymin = 0, ymax = 100,  # Stretch completely from bottom to top on Y axis
           fill = "#e5e5e5", #"#fac9c3",       # Background color
           alpha = 0.4) +
  geom_point(data = cond_mean_df,
             aes(x = x, y = y), shape = 8,
             color = '#f48373', size = 0.75, alpha = 0.4) +
  geom_point(data = as.data.frame(miss_info_mat),
             aes(x = x, y = y),
             color = '#999999', size = 0.7, alpha = 0.3) +
  geom_point(data = data.frame(matrix(data = softassign_k2_alt$top_x, ncol = 2, byrow = TRUE)),
             aes(X1, X2), color = '#f48373', size = 5, alpha = 1) +
  geom_point(data = as.data.frame(parent_locs[parent_locs$subpop == 6,c('x', 'y')]),
             aes(x, y), shape = 21, size = 1.5, fill = 'black', color = '#f48373', stroke = 1.15) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(color = '#e5e5e5', linetype = 'dashed'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = '#e5e5e5', fill = NA, linewidth = 0.75),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  ggtitle('Soft assign with K = 2 (ignore intra-tree covar)')


#combine into multipanel and export
k_assignment_type_explore <- cowplot::plot_grid(
  k2_hard_plot, k2_soft_plot, k2_soft_alt_plot,
  k3_hard_plot, k3_soft_plot, k3_soft_alt_plot,
  nrow = 2)

ggsave(k_assignment_type_explore,
       filename = here('plots', 'tmp', 'k_assignment_type_explore.png'),
       width = 10*1.5, height = 6*1.5,
       bg = 'white')



#################################
### CODE NOT CURRENTLY IN USE ###
#################################
# library(phytools)
# 
# x <- setNames(miss_info_mat[,2], rownames(miss_info_mat))
# 
# dotTree(drop.tip(tree_list_subset[[1]], MISSING_NODES), x)
# 
# plotTree.barplot(
#   drop.tip(tree_list_subset[[1]], MISSING_NODES),
#   x,
#   args.barplot = list(
#     col = "steelblue",
#     border = NA
#   )
# )

