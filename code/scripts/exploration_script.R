###############################################
###############################################
### SCRIPT FOR TESTING NEW PROJECT ELEMENTS ###
###############################################
###############################################

#####################
### SCRIPT SET-UP ###
#####################

### PACKAGES ###
library(here)
library(dplyr)
library(ggplot2)
library(ape)
library(phytools)
library(treesliceR)


### LOADING CUSTOM FUNCTIONS ###
source(here('code', 'scripts', 'location_est_code.R'))


### LOADING DATA (GENEALOGIES, INDIVIDUAL INFO, ETC...) ###
#here('simulation_output', 'tree_dir')
files_vec <- list.files(here('simulation_output'), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 


indiv_node_info <- read.delim(here('simulation_output', 'indiv_node_info.txt'), header = TRUE)
indiv_node_info$n_node <- paste0('n', indiv_node_info$node)

tree_position_info <- read.delim(here('simulation_output', 'tree_position_info.txt'), header = TRUE)



######################################################
### PROCESSING GENEALOGIES AND SPATIAL INFORMATION ###
######################################################

#the evolutionary rate matrix to use in location estimation
true_rate_mat <- diag(4, nrow = 2)

#removing the missing individual's nodes from the spatial information
MISSING_INDIV <- '72' #19 is very wrong
MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]#[2] #4 is weird


### PROCESSING WITH FULL TREES ###
miss_info <- indiv_node_info[!indiv_node_info$n_node %in% MISSING_NODES,]
rownames(miss_info) <- miss_info$n_node
miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])

vcv_list <- calc_vcv_multitree(tree_list = tree_list)

data_prep_list <- data_prep_multitree(vcv_list,
                                      trait = miss_info_mat, 
                                      rate_mat = true_rate_mat, 
                                      missing_inds = MISSING_NODES)

cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list)

tree_mean_vec <- c(t(tree_position_info[,c('x', 'y')]))

cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
                                      prep_list = data_prep_list,
                                      tree_means = tree_mean_vec)


### CHOPPED TREES VERSION ###
chop_time <- 2500

tree_list_chopped <- lapply(tree_list, function(x, time) {
  if (max(phytools::nodeHeights(x)) > time) {
    return(treesliceR::squeeze_root(tree = x, time = time))
  } else {
    return(x)
  }
}, time = chop_time)

vcv_list_chopped <- calc_vcv_multitree(tree_list = tree_list_chopped)

data_prep_list_chopped <- data_prep_multitree(vcv_list_chopped,
                                              trait = miss_info_mat, 
                                              rate_mat = true_rate_mat, 
                                              missing_inds = MISSING_NODES)

cond_var_elements_list_chopped <- calc_cond_var_elements_multitree(prep_list = data_prep_list_chopped)



cond_mean_list_chopped <- calc_cond_mean_list(cond_param_list = cond_var_elements_list_chopped,
                                              prep_list = data_prep_list_chopped,
                                              tree_means = tree_mean_vec)



####################################
### ESTIMATING UNKNOWN LOCATIONS ###
####################################

#using all the trees currently
tree_count <- length(cond_var_elements_list)

loc_est <- custom_metrop(iters = 1500, #number of MCMC iterations after the initial step
                         init_params = rep(50, 2), #c(50, 50), #initial parameter values
                         cond_param_list = cond_var_elements_list[1:tree_count],
                         prep_list = data_prep_list[1:tree_count],
                         tree_means = cond_mean_list[1:tree_count],
                         weights = 1,
                         proposal_scale = c(0.1, 0.1), #first proposal is for locations, 2nd is for root states
                         miss_n = 2,
                         type = c('fixed_mean', 'mean_estimate')[1],
                         report_progress = TRUE)


loc_est_weights <- custom_metrop(iters = 2500, #number of MCMC iterations after the initial step
                                 init_params = rep(50, 2), #c(50, 50), #initial parameter values
                                 cond_param_list = cond_var_elements_list[1:tree_count],
                                 prep_list = data_prep_list[1:tree_count],
                                 tree_means = cond_mean_list[1:tree_count],
                                 weights = tree_position_info$span[1:tree_count]/sum(tree_position_info$span[1:tree_count]),
                                 proposal_scale = c(4, 4), #first proposal is for locations, 2nd is for root states
                                 miss_n = 2,
                                 type = c('fixed_mean', 'mean_estimate')[1],
                                 report_progress = TRUE)

loc_est_chopped <- custom_metrop(iters = 2500, #number of MCMC iterations after the initial step
                                 init_params = rep(50, 2), #c(50, 50), #initial parameter values
                                 cond_param_list = cond_var_elements_list_chopped[1:tree_count],
                                 prep_list = data_prep_list_chopped[1:tree_count],
                                 tree_means = cond_mean_list_chopped[1:tree_count],
                                 weights = tree_position_info$span[1:tree_count]/sum(tree_position_info$span[1:tree_count]),
                                 proposal_scale = c(4.5, 4.5), #first proposal is for locations, 2nd is for root states
                                 miss_n = 2,
                                 type = c('fixed_mean', 'mean_estimate')[1],
                                 report_progress = TRUE)

loc_est_chopped_meanest <- custom_metrop(iters = 1500, #number of MCMC iterations after the initial step
                                 init_params = rep(50, 2 + 200*2), #c(50, 50), #initial parameter values
                                 cond_param_list = cond_var_elements_list_chopped[1:200],
                                 prep_list = data_prep_list_chopped[1:200],
                                 #tree_means = cond_mean_list_chopped[1:tree_count],
                                 weights = tree_position_info$span[1:200]/sum(tree_position_info$span[1:200]),
                                 proposal_scale = c(4.5, 0.3), #first proposal is for locations, 2nd is for root states
                                 miss_n = 2,
                                 type = c('fixed_mean', 'mean_estimate')[2],
                                 report_progress = TRUE)


different_posteriors_plot <- loc_est_chopped$post_mat %>% 
  as.data.frame() %>% 
  mutate(order = 1:n()) %>% 
  ggplot() +
 geom_density_2d(data = . %>%
                   slice_tail(n = 1000),
                 aes(x= V1, y = V2),
                 color = '#01cdfe', 
                 linewidth = 0.6,
                 bins = 4) +
 geom_density_2d(data = loc_est_weights$post_mat %>%
                   as.data.frame() %>%
                   mutate(order = 1:n()) %>%
                   slice_tail(n = 1000),
                 aes(x= V1, y = V2),
                 color = '#419873',
                 linewidth = 0.6,
                 bins = 4) +
  geom_density_2d(data = loc_est_chopped_meanest$post_mat %>%
                    as.data.frame() %>% 
                    mutate(order = 1:n()) %>% 
                    slice_tail(n = 1000),
                  aes(x= V1, y = V2),
                  color = 'gray', 
                  linewidth = 0.6,
                  bins = 4) +
  geom_density_2d(data = loc_est$post_mat %>%
                    as.data.frame() %>% 
                    mutate(order = 1:n()) %>% 
                    slice_tail(n = 1000),
                  aes(x= V1, y = V2),
                  color = '#ff71ce',
                  linewidth = 0.6,
                  bins = 4) +
  geom_point(data = indiv_node_info %>% 
               filter(n_node %in% MISSING_NODES),
             aes(x, y), shape = 21, color = 'purple', size = 3, stroke = 1.25) +
  geom_point(data = indiv_node_info %>% 
               filter(!n_node %in% MISSING_NODES),
             aes(x, y), color = 'black') +
  theme_bw() +
  #scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = '#d0fced', linetype = 'dashed'),
        plot.title = element_text(size = 10),
        legend.position = 'none') +
  #ggtitle(paste0('INDIV = ', 1, '; Green = span weighted; Blue = chopped at 2500; Pink = full unweighted')) +
  xlab('x') + ylab('y')


ggsave(here('plots', 'different_posteriors_plot.png'), 
       different_posteriors_plot, 
       width = 5, height = 4, units = "in")



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# loc_est_chopped$post_mat %>% 
#   as.data.frame() %>% 
#   mutate(order = 1:n()) %>% 
#   ggplot() +
#   geom_line(aes(x = order, y = V1)) +
#   theme_bw()
# 
# loc_est_weights$post_mat %>% 
#   as.data.frame() %>% 
#   mutate(order = 1:n()) %>% 
#   ggplot() +
#   geom_line(aes(x = order, y = V1)) +
#   theme_bw()
# 
# loc_est$post_mat %>% 
#   as.data.frame() %>% 
#   mutate(order = 1:n()) %>% 
#   ggplot() +
#   geom_line(aes(x = order, y = V1)) +
#   theme_bw()
# 
# loc_est_chopped_meanest$post_mat %>%
#   as.data.frame() %>%
#   mutate(order = 1:n()) %>%
#   ggplot() +
#   geom_line(aes(x = order, y = V1)) +
#   theme_bw()


# for (INDIV in c('19', '60', '48'))
# 
# MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == '30']#[2] #4 is weird
# miss_info <- indiv_node_info[!indiv_node_info$n_node %in% MISSING_NODES,]
# rownames(miss_info) <- miss_info$n_node
# miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# 
# 
# 
# vcv_list <- calc_vcv_multitree(tree_list = tree_list)
# 
# data_prep_list <- data_prep_multitree(vcv_list,
#                                       trait = miss_info_mat, 
#                                       rate_mat = true_rate_mat, 
#                                       missing_inds = MISSING_NODES)
# 
# cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list)
# 
# tree_mean_vec <- c(t(tree_position_info[,c('x', 'y')]))
# 
# cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                 prep_list = data_prep_list,
#                                 tree_means = tree_mean_vec)
# 
# ### chopping ###
# 
# tree_list_chopped <- lapply(tree_list, function(x) {
#   if (max(phytools::nodeHeights(x)) > 2500) {
#     return(treesliceR::squeeze_root(tree = x, time = 2500))
#   } else {
#     return(x)
#   }
# })
# 
# vcv_list_chopped <- calc_vcv_multitree(tree_list = tree_list_chopped)
# 
# data_prep_list_chopped <- data_prep_multitree(vcv_list_chopped,
#                                       trait = miss_info_mat, 
#                                       rate_mat = true_rate_mat, 
#                                       missing_inds = MISSING_NODES)
# 
# cond_var_elements_list_chopped <- calc_cond_var_elements_multitree(prep_list = data_prep_list_chopped)
# 
# 
# 
# cond_mean_list_chopped <- calc_cond_mean_list(cond_param_list = cond_var_elements_list_chopped,
#                                       prep_list = data_prep_list_chopped,
#                                       tree_means = tree_mean_vec)
# 
# ###
# 
# 
# tree_count <- 450 #567 #150
# 
# loc_est <- custom_metrop(iters = 1500, #number of MCMC iterations after the initial step
#                           init_params = rep(50, 2), #c(50, 50), #initial parameter values
#                           cond_param_list = cond_var_elements_list[1:tree_count],
#                           prep_list = data_prep_list[1:tree_count],
#                           tree_means = cond_mean_list[1:tree_count],
#                          weights = 1,
#                           proposal_scale = c(0.1, 0.1), #first proposal is for locations, 2nd is for root states
#                           miss_n = 2,
#                           type = c('fixed_mean', 'mean_estimate')[1],
#                           report_progress = TRUE)
# 
# 
# loc_est_weights <- custom_metrop(iters = 1500, #number of MCMC iterations after the initial step
#                          init_params = rep(50, 2), #c(50, 50), #initial parameter values
#                          cond_param_list = cond_var_elements_list[1:tree_count],
#                          prep_list = data_prep_list[1:tree_count],
#                          tree_means = cond_mean_list[1:tree_count],
#                          weights = tree_position_info$span[1:tree_count]/sum(tree_position_info$span[1:tree_count]),
#                          proposal_scale = c(3, 3), #first proposal is for locations, 2nd is for root states
#                          miss_n = 2,
#                          type = c('fixed_mean', 'mean_estimate')[1],
#                          report_progress = TRUE)
# 
# loc_est_chopped <- custom_metrop(iters = 1500, #number of MCMC iterations after the initial step
#                                  init_params = rep(50, 2), #c(50, 50), #initial parameter values
#                                  cond_param_list = cond_var_elements_list_chopped[1:tree_count],
#                                  prep_list = data_prep_list_chopped[1:tree_count],
#                                  tree_means = cond_mean_list_chopped[1:tree_count],
#                                  weights = tree_position_info$span[1:tree_count]/sum(tree_position_info$span[1:tree_count]),
#                                  proposal_scale = c(3.5, 3.5), #first proposal is for locations, 2nd is for root states
#                                  miss_n = 2,
#                                  type = c('fixed_mean', 'mean_estimate')[1],
#                                  report_progress = TRUE)
# 
# 
# # plot_no_weights <- loc_est$post_mat %>% 
# #   as.data.frame() %>% 
# #   mutate(order = 1:n()) %>% 
# #   ggplot() +
# #   ##geom_line(aes(x = V1, y = V2)) +
# #   # geom_path(aes(x=V1, y=V2, color = order),linewidth = 0.85) +
# #   # scale_colour_gradient(
# #   #   low = "#ffe5ea",
# #   #   high = "#b20000"
# #   # ) +
# #   # geom_point(data = cond_mean_df,
# #   #            aes(x, y, fill = node), 
# #   #            shape = 21, color = 'white', size = 2.5, stroke = 0.25) +
# #   geom_density_2d(data = . %>%
# #                     slice_tail(n = 500),
# #                   aes(x= V1, y = V2),
# #                   color = 'green', linewidth = 0.4) +
# #   geom_point(data = indiv_node_info %>% 
# #                filter(n_node %in% MISSING_NODES),
# #              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
# #   geom_point(data = indiv_node_info %>% 
# #                filter(!n_node %in% MISSING_NODES),
# #              aes(x, y), color = 'black') +
# #   theme_bw() +
# #   #scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
# #   theme(panel.grid.minor = element_blank(),
# #         panel.grid.major = element_line(color = '#d0fced', linetype = 'dashed'),
# #         plot.title = element_text(size = 10),
# #         legend.position = 'none') +
# #   ggtitle('Equal tree weights (full trees)') +
# #   xlab('x') + ylab('y')
# 
# 
# # plot_weights <- loc_est_weights$post_mat %>% 
# #   as.data.frame() %>% 
# #   mutate(order = 1:n()) %>% 
# #   ggplot() +
# #   ##geom_line(aes(x = V1, y = V2)) +
# #   # geom_path(aes(x=V1, y=V2, color = order),linewidth = 0.85) +
# #   # scale_colour_gradient(
# #   #   low = "#ffe5ea",
# #   #   high = "#b20000"
# #   # ) +
# #   # geom_point(data = cond_mean_df,
# #   #            aes(x, y, fill = node), 
# #   #            shape = 21, color = 'white', size = 2.5, stroke = 0.25) +
# #   geom_density_2d(data = . %>%
# #                     slice_tail(n = 500),
# #                   aes(x= V1, y = V2),
# #                   color = 'green', linewidth = 0.4) +
# #   geom_point(data = indiv_node_info %>% 
# #                filter(n_node %in% MISSING_NODES),
# #              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
# #   geom_point(data = indiv_node_info %>% 
# #                filter(!n_node %in% MISSING_NODES),
# #              aes(x, y), color = 'black') +
# #   theme_bw() +
# #   #scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
# #   theme(panel.grid.minor = element_blank(),
# #         panel.grid.major = element_line(color = '#d0fced', linetype = 'dashed'),
# #         plot.title = element_text(size = 10),
# #         legend.position = 'none') +
# #   ggtitle('Span weighted (full trees)') +
# #   xlab('x') + ylab('y')
# 
# different_posteriors_plot <- loc_est_chopped$post_mat %>% 
#   as.data.frame() %>% 
#   mutate(order = 1:n()) %>% 
#   ggplot() +
#   ##geom_line(aes(x = V1, y = V2)) +
#   # geom_path(aes(x=V1, y=V2, color = order),linewidth = 0.85) +
#   # scale_colour_gradient(
#   #   low = "#ffe5ea",
#   #   high = "#b20000"
#   # ) +
#   # geom_point(data = cond_mean_df,
#   #            aes(x, y, fill = node), 
#   #            shape = 21, color = 'white', size = 2.5, stroke = 0.25) +
#   geom_density_2d(data = . %>%
#                     slice_tail(n = 1000),
#                   aes(x= V1, y = V2),
#                   color = '#01cdfe', linewidth = 0.4) +
#   geom_density_2d(data = loc_est_weights$post_mat %>%
#                     as.data.frame() %>% 
#                     mutate(order = 1:n()) %>% 
#                     slice_tail(n = 1000),
#                   aes(x= V1, y = V2),
#                   color = '#419873', linewidth = 0.4) +
#   geom_density_2d(data = loc_est$post_mat %>%
#                     as.data.frame() %>% 
#                     mutate(order = 1:n()) %>% 
#                     slice_tail(n = 1000),
#                   aes(x= V1, y = V2),
#                   color = '#ff71ce', linewidth = 0.4) +
#   geom_point(data = indiv_node_info %>% 
#                filter(n_node %in% MISSING_NODES),
#              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
#   geom_point(data = indiv_node_info %>% 
#                filter(!n_node %in% MISSING_NODES),
#              aes(x, y), color = 'black') +
#   theme_bw() +
#   #scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = '#d0fced', linetype = 'dashed'),
#         plot.title = element_text(size = 10),
#         legend.position = 'none') +
#   ggtitle(paste0('INDIV = ', 1, '; Green = span weighted; Blue = chopped at 2500; Pink = full unweighted')) +
#   xlab('x') + ylab('y')
# 
# # loc_est_weights$post_mat %>% 
# #   as.data.frame() %>% 
# #   mutate(order = 1:n()) %>% 
# #   ggplot() +
# #   geom_line(aes(x = order, y = V1))
# # 
# # loc_est_chopped$post_mat %>% 
# #   as.data.frame() %>% 
# #   mutate(order = 1:n()) %>% 
# #   ggplot() +
# #   geom_line(aes(x = order, y = V1))
# 
# 
# 
# 
# cond_mean_list_chopped_df <- list()
# for (i in 1:length(cond_mean_list_chopped)) {
#   
#   cond_mean_list_chopped_df[[i]] <- data.frame(tree_ind = i,
#                                             node = gsub(":[xy]", "", names(cond_mean_list_chopped[[i]][grep('x',rownames(cond_mean_list_chopped[[i]])),])),
#                                             x = cond_mean_list_chopped[[i]][grep('x',rownames(cond_mean_list_chopped[[i]])),],
#                                             y = cond_mean_list_chopped[[i]][grep('y',rownames(cond_mean_list_chopped[[i]])),])
# }
# 
# cond_mean_df_chopped <- do.call(rbind, cond_mean_list_chopped_df)
# 
# cond_mean_list_df <- list()
# for (i in 1:length(cond_mean_list)) {
#   
#   cond_mean_list_df[[i]] <- data.frame(tree_ind = i,
#                                                node = gsub(":[xy]", "", names(cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),])),
#                                                x = cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),],
#                                                y = cond_mean_list[[i]][grep('y',rownames(cond_mean_list[[i]])),])
# }
# 
# cond_mean_df <- do.call(rbind, cond_mean_list_df)
# 
# dist_compare_boxplot <- data.frame(dist = c(apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# }),apply(cond_mean_df_chopped[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })),
# type = c(rep('full', nrow(cond_mean_df) ), rep('chopped', nrow(cond_mean_df) ))
# 
# ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = type, y = dist)) +
#   ylab('Dist to true location') +
#   theme_bw() +
#   ggtitle(paste0("INDIV = ", 1, "; Dist(true location, conditional mean)"))
# 
# # cond_mean_df_combined %>% 
# #   ggplot() +
# #   geom_point(aes(x, y, color = span)) +
# #   geom_point(data = indiv_node_info %>% 
# #                filter(!n_node %in% MISSING_NODES),
# #              aes(x, y), color = 'blue', alpha = 0.25, size = 0.8) +
# #   scale_color_gradient(low = '#b2d8b2', high = '#004000') +
# #   geom_point(data = indiv_node_info %>% 
# #                filter(n_node %in% MISSING_NODES),
# #              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
# #   theme_bw()
# 
# 
# 
# cond_mean_df <- do.call(rbind, cond_mean_list)
# 
# cond_mean_df_list <- list()
# for (i in 1:length(cond_mean_list)) {
#   cond_mean_df_list[[i]] <- data.frame(tree_index = i - 1,
#                                        node = gsub(":[xy]", "", names(cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),])),
#                                        x = cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),],
#                                        y = cond_mean_list[[i]][grep('y',rownames(cond_mean_list[[i]])),])
# }
# 
# cond_mean_df <- do.call(rbind, cond_mean_df_list)
# 
# cond_mean_df$dist_to_true <- apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })
# 
# cond_mean_df_combined <- cond_mean_df %>% 
#   #rename('tree_index' = tree_ind) %>% 
#   left_join(., tree_position_info[,c('tree_index', 'span')], 
#             by = 'tree_index')
# 
# cor_val <- cor.test(cond_mean_df_combined$dist_to_true, cond_mean_df_combined$span, method = 'spearman')
# 
# 
# span_dist_plot <- cond_mean_df_combined %>% 
#   ggplot() +
#   geom_point(aes(x = span, y = dist_to_true), alpha = 0.5) +
#   theme_bw() +
#   xlab('Span length') +
#   ylab('Distance to true location') +
#   ggtitle(paste0('INDIV = ',  1, '; Spearman cor: ', round(cor_val$estimate, 3)))
# 
# 
# span_color_map_plot <- cond_mean_df_combined %>% 
#   ggplot() +
#   geom_point(aes(x, y, color = span)) +
#   geom_point(data = indiv_node_info %>% 
#                filter(!n_node %in% MISSING_NODES),
#              aes(x, y), color = 'blue', alpha = 0.25, size = 0.8) +
#   scale_color_gradient(low = '#b2d8b2', high = '#004000') +
#   geom_point(data = indiv_node_info %>% 
#                filter(n_node %in% MISSING_NODES),
#              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
#   theme_bw() +
#   ggtitle(paste0('INDIV = ', 1) )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mean_full <- apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })
# 
# mean_chopped <- apply(cond_mean_df_chopped[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })
# 
# data.frame(distance_to_true = c(mean_full, mean_chopped),
#            type = c(rep('full', length(mean_full)), rep('chopped', length(mean_chopped)))) %>% 
#   ggplot() +
#   geom_density(aes(x = distance_to_true, group = type, color = type, fill = type), alpha = 0.5) +
#   theme_bw()
# 
# 
# #cowplot::plot_grid(plot_no_weights, plot_weights, plot_chopped, ncol = 3)
# 
# 
# # cond_mean_output <- calc_cond_mean(
# #   covar_12 = cond_var_elements_list[[2]]$covar_12,
# #   inv_covar22 = cond_var_elements_list[[2]]$inv_covar_22,
# #   trait_vec = data_prep_list[[2]]$vector_trait,
# #   means = unlist(tree_position_info[2,c('x', 'y')]),
# #   no_ind = data_prep_list[[2]]$no_ind,
# #   yes_ind =  data_prep_list[[2]]$yes_ind
# # )
# 
# 
# cond_mean_df <- do.call(rbind, cond_mean_list)
# 
# cond_mean_df_list <- list()
# for (i in 1:length(cond_mean_list)) {
#   cond_mean_df_list[[i]] <- data.frame(tree_index = i - 1,
#                                     node = gsub(":[xy]", "", names(cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),])),
#                                     x = cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),],
#                                     y = cond_mean_list[[i]][grep('y',rownames(cond_mean_list[[i]])),])
# }
# 
# cond_mean_df <- do.call(rbind, cond_mean_df_list)
# 
# cond_mean_df$dist_to_true <- apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })
# 
# cond_mean_df_combined <- cond_mean_df %>% 
#   #rename('tree_index' = tree_ind) %>% 
#   left_join(., tree_position_info[,c('tree_index', 'span')], 
#             by = 'tree_index')
# 
# cor_val <- cor.test(cond_mean_df_combined$dist_to_true, cond_mean_df_combined$span, method = 'spearman')
# 
# 
# span_dist_plot <- cond_mean_df_combined %>% 
#   ggplot() +
#   geom_point(aes(x = span, y = dist_to_true), alpha = 0.5) +
#   theme_bw() +
#   xlab('Span length') +
#   ylab('Distance to true location') +
#   ggtitle(paste0('Spearman cor: ', round(cor_val$estimate, 3)))
# 
# 
# span_color_map_plot <- cond_mean_df_combined %>% 
#   ggplot() +
#   geom_point(aes(x, y, color = span)) +
#   geom_point(data = indiv_node_info %>% 
#                filter(!n_node %in% MISSING_NODES),
#              aes(x, y), color = 'blue', alpha = 0.25, size = 0.8) +
#   scale_color_gradient(low = '#b2d8b2', high = '#004000') +
#   geom_point(data = indiv_node_info %>% 
#                filter(n_node %in% MISSING_NODES),
#              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
#   theme_bw()
# 
# 
# 
# 
# span_dist_plot
# span_color_map_plot
# 
# 
# 
# 
# 
# 
# 
# 
# loc_est$post_mat %>% 
#   as.data.frame() %>% 
#   mutate(order = 1:n()) %>% 
#   ggplot() +
#   ##geom_line(aes(x = V1, y = V2)) +
#   # geom_path(aes(x=V1, y=V2, color = order),linewidth = 0.85) +
#   # scale_colour_gradient(
#   #   low = "#ffe5ea",
#   #   high = "#b20000"
#   # ) +
#   # geom_point(data = cond_mean_df,
#   #            aes(x, y, fill = node), 
#   #            shape = 21, color = 'white', size = 2.5, stroke = 0.25) +
#   geom_density_2d(data = . %>%
#                     slice_tail(n = 500),
#                   aes(x= V1, y = V2),
#                   color = 'green', linewidth = 0.4) +
#   geom_point(data = indiv_node_info %>% 
#                filter(n_node %in% MISSING_NODES),
#              aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
#   geom_point(data = indiv_node_info %>% 
#                filter(!n_node %in% MISSING_NODES),
#              aes(x, y), color = 'black') +
#   theme_bw() +
#   #scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = '#d0fced', linetype = 'dashed'),
#         plot.title = element_text(size = 10),
#         legend.position = 'none') +
#   ggtitle(paste0('Using both nodes (using true root locations)')) +
#   xlab('x') + ylab('y')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# tree_list_chopped <- lapply(tree_list, function(x) {
#   if (max(phytools::nodeHeights(x)) > 800) {
#     return(treesliceR::squeeze_root(tree = x, time = 800))
#   } else {
#     return(x)
#   }
# })
# 
# 
# tip_order_chopped <- tree_list_chopped[[1]]$tip.label
# vcv_list_chopped <- lapply(tree_list_chopped, function(x, tip_order) {
#   vcv_mat <- vcv.phylo(x)
#   return(vcv_mat[match(tip_order, rownames(vcv_mat)),
#                  match(tip_order, colnames(vcv_mat))]
#   )
# }, tip_order = tip_order_chopped)
# 
# nfm_prep_list_chopped <- lapply(vcv_list_chopped, function(x, trait, rate_mat, miss_inds) {
#   missing_trait_lik_prep_cond_nfm(trait = trait,
#                                   phylo_covar = x,
#                                   rate_mat = rate_mat,
#                                   missing_inds = miss_inds)
# }, trait = miss_info_mat, rate_mat = true_rate_mat, miss_inds = MISSING_NODES)
# 
# 
# calc_cond_var_elements_list_chopped <- lapply(nfm_prep_list_chopped, function(x) {
#   calc_cond_var_elements(x$kron_varcovar_ratemat, x$no_ind, x$yes_ind)
# })
# 
# 
# cond_mean_list_chopped <- list()
# for (i in 1:length(calc_cond_var_elements_list_chopped)) {
#   
#   cond_mean_output <- calc_cond_mean(
#     covar_12 = calc_cond_var_elements_list_chopped[[i]]$covar_12,
#     inv_covar22 = calc_cond_var_elements_list_chopped[[i]]$inv_covar_22,
#     trait_vec = nfm_prep_list_chopped[[i]]$vector_trait,
#     means = unlist(tree_position_info[i,c('x', 'y')]),
#     no_ind = nfm_prep_list_chopped[[i]]$no_ind,
#     yes_ind =  nfm_prep_list_chopped[[i]]$yes_ind)
#   
#   cond_mean_list_chopped[[i]] <- data.frame(tree_ind = i,
#                                             node = gsub(":[xy]", "", names(cond_mean_output[grep('x',rownames(cond_mean_output)),])),
#                                             x = cond_mean_output[grep('x',rownames(cond_mean_output)),],
#                                             y = cond_mean_output[grep('y',rownames(cond_mean_output)),])
# }
# 
# cond_mean_df_chopped <- do.call(rbind, cond_mean_list_chopped)
# 
# cond_mean_list <- list()
# for (i in 1:length(calc_cond_var_elements_list)) {
#   
#   cond_mean_output <- calc_cond_mean(
#     covar_12 = calc_cond_var_elements_list[[i]]$covar_12,
#     inv_covar22 = calc_cond_var_elements_list[[i]]$inv_covar_22,
#     trait_vec = nfm_prep_list[[i]]$vector_trait,
#     means = unlist(tree_position_info[i,c('x', 'y')]),
#     no_ind = nfm_prep_list[[i]]$no_ind,
#     yes_ind =  nfm_prep_list[[i]]$yes_ind)
#   
#   cond_mean_list[[i]] <- data.frame(tree_ind = i,
#                                     node = gsub(":[xy]", "", names(cond_mean_output[grep('x',rownames(cond_mean_output)),])),
#                                     x = cond_mean_output[grep('x',rownames(cond_mean_output)),],
#                                     y = cond_mean_output[grep('y',rownames(cond_mean_output)),])
# }
# 
# cond_mean_df <- do.call(rbind, cond_mean_list)
# 
# # 
# # data.frame(dist = c(apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
# #   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# # }),apply(cond_mean_df_chopped[,c('x', 'y')], 1, function(x) {
# #   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# # })),
# # type = c(rep('full', nrow(cond_mean_df) ), rep('chopped', nrow(cond_mean_df) ))
# # 
# # ) %>%
# #   ggplot() +
# #   geom_boxplot(aes(x = type, y = dist))
# 
# 
# 
# mean_full <- apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })
# 
# mean_chopped <- apply(cond_mean_df_chopped[,c('x', 'y')], 1, function(x) {
#   dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
# })
# 
# data.frame(distance_to_true = c(mean_full, mean_chopped),
#            type = c(rep('full', length(mean_full)), rep('chopped', length(mean_chopped)))) %>% 
#   ggplot() +
#   geom_density(aes(x = distance_to_true, group = type, color = type, fill = type), alpha = 0.5) +
#   theme_bw()
# 
# 
# 
# 
# 
# 
# dist_list[[INDIV_INDEX]] <- mean_full - mean_chopped
# message(INDIV_INDEX)
# }
# 
# data.frame(a = unlist(dist_list)) %>% 
#   ggplot() +
#   geom_violin(aes(x = 'm', y = a)) +
#   geom_hline(aes(yintercept = 0))
# 
# 
# chopped_plot <- ggplot(data = cond_mean_df_chopped) + 
#   geom_point(aes(x, y, fill = node), 
#              shape = 21, color = 'white', size = 2.5, stroke = 0.25) +
#   geom_density_2d(aes(x= x, y = y)) +
#   scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
#   geom_point(data = indiv_node_info %>% 
#                filter(!n_node %in% MISSING_NODES),
#              aes(x, y), color = 'black', size = 1) +
#   geom_point(data = indiv_node_info %>% 
#                filter(n_node %in% MISSING_NODES),
#              aes(x, y), shape = 21, color = 'purple', size = 4, stroke = 2) +
#   theme_bw() +
#   theme(legend.position = 'none') +
#   ggtitle('chopped')
# 
# 
# 
# 
# 
# 
# # locs <- c(50, 50)
# # means <- c(t(tree_position_info[tree_indices,c('x', 'y')]))
# # log_lik_vec <- vector(mode = 'numeric', length = length(cond_var_elements_list))
# # miss_n <- 2
# # for (i in seq_along(log_lik_vec)) {
# #   #message(i)
# #   cond_mean <- calc_cond_mean(covar_12 = cond_var_elements_list[[i]]$covar_12, 
# #                               inv_covar22 = cond_var_elements_list[[i]]$inv_covar_22, 
# #                               trait_vec = data_prep_list[[i]]$vector_trait, 
# #                               means = means[c(i*2 - 1, i*2)], 
# #                               no_ind = data_prep_list[[i]]$no_ind, 
# #                               yes_ind = data_prep_list[[i]]$yes_ind)
# #   
# #   
# #   log_lik_vec[i] <- mvtnorm::dmvnorm(rep(locs, times = miss_n), 
# #                                      cond_mean, 
# #                                      sigma = as.matrix(Matrix::forceSymmetric(cond_var_elements_list[[i]]$cond_var)),
# #                                      log = TRUE)
# # }
# #   
# #   
# #   
# #   
# # custom_metrop_cond_nfm_fixedmean(iters = 1500, 
# #                                  means = c(t(tree_position_info[tree_indices,c('x', 'y')])),
# #                                  init_params = c(50, 50), #sample(0:13, tree_count*2 + 2, replace = TRUE),#sample((100-13):100, tree_count*2 + 2, replace = TRUE),
# #                                  cond_param_list = calc_cond_var_elements_list[tree_indices],
# #                                  nfm_prep_list = nfm_prep_list[tree_indices],
# #                                  proposal_scale = 1,
# #                                  miss_n = 2)
# # 
# # 
# # log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
# #   
# #   for (i in seq_along(log_lik_vec)) {
# #     #message(i)
# #     cond_mean <- calc_cond_mean(covar_12 = cond_list[[i]]$covar_12, 
# #                                 inv_covar22 = cond_list[[i]]$inv_covar_22, 
# #                                 trait_vec = prep_list[[i]]$vector_trait, 
# #                                 means = means[c(i*2 - 1, i*2)], 
# #                                 no_ind = prep_list[[i]]$no_ind, 
# #                                 yes_ind = prep_list[[i]]$yes_ind)
# #     
# #     
# #     log_lik_vec[i] <- mvtnorm::dmvnorm(rep(locs, times = miss_n), 
# #                                        cond_mean, 
# #                                        sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
# #                                        log = TRUE)
# # 
# # 
# # 
# # tip_order <- tree_list[[1]]$tip.label
# # vcv_list <- lapply(tree_list, function(x, tip_order) {
# #   vcv_mat <- vcv.phylo(x)
# #   return(vcv_mat[match(tip_order, rownames(vcv_mat)),
# #                  match(tip_order, colnames(vcv_mat))]
# #   )
# # }, tip_order = tip_order)
# # 
# # 
# # true_rate_mat <- matrix(c(4, 0, 0, 4), nrow = 2)
# # miss_info <- indiv_node_info[!indiv_node_info$n_node %in% MISSING_NODES,]
# # #miss_info <- indiv_node_info[indiv_node_info$indiv != MISSING_INDIV,]
# # rownames(miss_info) <- miss_info$n_node
# # miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
# # 
# # nfm_prep_list <- lapply(vcv_list, function(x, trait, rate_mat, miss_inds) {
# #   missing_trait_lik_prep_cond_nfm(trait = trait,
# #                                   phylo_covar = x,
# #                                   rate_mat = rate_mat,
# #                                   missing_inds = miss_inds)
# # }, trait = miss_info_mat, rate_mat = true_rate_mat, miss_inds = MISSING_NODES)
# # 
# # 
# # calc_cond_var_elements_list <- lapply(nfm_prep_list, function(x) {
# #   calc_cond_var_elements(x$kron_varcovar_ratemat, x$no_ind, x$yes_ind)
# # })
# 
# 
# 
# 
# posterior_plot_list <- list()
# dist_compare_boxplot_list <- list()
# span_dist_plot_list <- list()
# span_color_map_plot_list <- list()
# 
# 
# true_rate_mat <- diag(4, nrow = 2)
# 
# #35 --> really different
# 
# 
# #indivs --> 19
# 
# for (INDIV in c('19', '10', '48')) {
#   INDIV <- '10'
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == INDIV]#[2] #4 is weird
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% MISSING_NODES,]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
#   
#   
#   
#   vcv_list <- calc_vcv_multitree(tree_list = tree_list)
#   
#   data_prep_list <- data_prep_multitree(vcv_list,
#                                         trait = miss_info_mat, 
#                                         rate_mat = true_rate_mat, 
#                                         missing_inds = MISSING_NODES)
#   
#   cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list)
#   
#   tree_mean_vec <- c(t(tree_position_info[,c('x', 'y')]))
#   
#   cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                         prep_list = data_prep_list,
#                                         tree_means = tree_mean_vec)
#   
#   ### chopping ###
#   
#   tree_list_chopped <- lapply(tree_list, function(x) {
#     if (max(phytools::nodeHeights(x)) > 2500) {
#       return(treesliceR::squeeze_root(tree = x, time = 2500))
#     } else {
#       return(x)
#     }
#   })
#   
#   vcv_list_chopped <- calc_vcv_multitree(tree_list = tree_list_chopped)
#   
#   data_prep_list_chopped <- data_prep_multitree(vcv_list_chopped,
#                                                 trait = miss_info_mat, 
#                                                 rate_mat = true_rate_mat, 
#                                                 missing_inds = MISSING_NODES)
#   
#   cond_var_elements_list_chopped <- calc_cond_var_elements_multitree(prep_list = data_prep_list_chopped)
#   
#   
#   
#   cond_mean_list_chopped <- calc_cond_mean_list(cond_param_list = cond_var_elements_list_chopped,
#                                                 prep_list = data_prep_list_chopped,
#                                                 tree_means = tree_mean_vec)
#   
#   ###
#   
#   
#   tree_count <- 460 #567 #150
#   
#   loc_est <- custom_metrop(iters = 2000, #number of MCMC iterations after the initial step
#                            init_params = rep(50, 2), #c(50, 50), #initial parameter values
#                            cond_param_list = cond_var_elements_list[1:tree_count],
#                            prep_list = data_prep_list[1:tree_count],
#                            tree_means = cond_mean_list[1:tree_count],
#                            weights = 1,
#                            proposal_scale = c(0.1, 0.1), #first proposal is for locations, 2nd is for root states
#                            miss_n = 2,
#                            type = c('fixed_mean', 'mean_estimate')[1],
#                            report_progress = TRUE)
#   
#   
#   loc_est_weights <- custom_metrop(iters = 2000, #number of MCMC iterations after the initial step
#                                    init_params = rep(50, 2), #c(50, 50), #initial parameter values
#                                    cond_param_list = cond_var_elements_list[1:tree_count],
#                                    prep_list = data_prep_list[1:tree_count],
#                                    tree_means = cond_mean_list[1:tree_count],
#                                    weights = tree_position_info$span[1:tree_count]/sum(tree_position_info$span[1:tree_count]),
#                                    proposal_scale = c(3, 3), #first proposal is for locations, 2nd is for root states
#                                    miss_n = 2,
#                                    type = c('fixed_mean', 'mean_estimate')[1],
#                                    report_progress = TRUE)
#   
#   loc_est_chopped <- custom_metrop(iters = 2000, #number of MCMC iterations after the initial step
#                                    init_params = rep(50, 2), #c(50, 50), #initial parameter values
#                                    cond_param_list = cond_var_elements_list_chopped[1:tree_count],
#                                    prep_list = data_prep_list_chopped[1:tree_count],
#                                    tree_means = cond_mean_list_chopped[1:tree_count],
#                                    weights = tree_position_info$span[1:tree_count]/sum(tree_position_info$span[1:tree_count]),
#                                    proposal_scale = c(3.5, 3.5), #first proposal is for locations, 2nd is for root states
#                                    miss_n = 2,
#                                    type = c('fixed_mean', 'mean_estimate')[1],
#                                    report_progress = TRUE)
#   
#   
#   posterior_plot_list[[paste0('INDIV', INDIV)]] <- loc_est_chopped$post_mat %>% 
#     as.data.frame() %>% 
#     mutate(order = 1:n()) %>% 
#     ggplot() +
#     geom_density_2d(data = loc_est_weights$post_mat %>%
#                       as.data.frame() %>% 
#                       mutate(order = 1:n()) %>% 
#                       slice_tail(n = 1500),
#                     aes(x= V1, y = V2),
#                     color = '#419873', linewidth = 0.4) +
#     geom_density_2d(data = . %>%
#                       slice_tail(n = 1500),
#                     aes(x= V1, y = V2),
#                     color = '#01cdfe', linewidth = 0.4) +
#     geom_density_2d(data = loc_est$post_mat %>%
#                       as.data.frame() %>% 
#                       mutate(order = 1:n()) %>% 
#                       slice_tail(n = 1500),
#                     aes(x= V1, y = V2),
#                     color = '#ff71ce', linewidth = 0.4) +
#     geom_point(data = indiv_node_info %>% 
#                  filter(n_node %in% MISSING_NODES),
#                aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
#     geom_point(data = indiv_node_info %>% 
#                  filter(!n_node %in% MISSING_NODES),
#                aes(x, y), color = 'black') +
#     theme_bw() +
#     #scale_fill_manual(values = c("#FB8542", '#1E90FF')) +
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(color = '#d0fced', linetype = 'dashed'),
#           plot.title = element_text(size = 10),
#           legend.position = 'none') +
#     ggtitle(paste0('INDIV = ', INDIV, '; Green = span weighted; Blue = chopped at 2500; Pink = full unweighted')) +
#     xlab('x') + ylab('y')
#   
#   
#   cond_mean_list_chopped_df <- list()
#   for (i in 1:length(cond_mean_list_chopped)) {
#     
#     cond_mean_list_chopped_df[[i]] <- data.frame(tree_ind = i,
#                                                  node = gsub(":[xy]", "", names(cond_mean_list_chopped[[i]][grep('x',rownames(cond_mean_list_chopped[[i]])),])),
#                                                  x = cond_mean_list_chopped[[i]][grep('x',rownames(cond_mean_list_chopped[[i]])),],
#                                                  y = cond_mean_list_chopped[[i]][grep('y',rownames(cond_mean_list_chopped[[i]])),])
#   }
#   
#   cond_mean_df_chopped <- do.call(rbind, cond_mean_list_chopped_df)
#   
#   cond_mean_list_df <- list()
#   for (i in 1:length(cond_mean_list)) {
#     
#     cond_mean_list_df[[i]] <- data.frame(tree_ind = i,
#                                          node = gsub(":[xy]", "", names(cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),])),
#                                          x = cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),],
#                                          y = cond_mean_list[[i]][grep('y',rownames(cond_mean_list[[i]])),])
#   }
#   
#   cond_mean_df <- do.call(rbind, cond_mean_list_df)
#   
#   dist_compare_boxplot_list[[paste0('INDIV', INDIV)]] <- data.frame(dist = c(apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#     dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
#   }),apply(cond_mean_df_chopped[,c('x', 'y')], 1, function(x) {
#     dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
#   })),
#   type = c(rep('full', nrow(cond_mean_df) ), rep('chopped', nrow(cond_mean_df) ))
#   
#   ) %>%
#     ggplot() +
#     geom_boxplot(aes(x = type, y = dist)) +
#     ylab('Dist to true location') +
#     theme_bw() +
#     ggtitle(paste0("INDIV = ", INDIV, "; Dist(true location, conditional mean)"))
#   
#   
#   cond_mean_df <- do.call(rbind, cond_mean_list)
#   
#   cond_mean_df_list <- list()
#   for (i in 1:length(cond_mean_list)) {
#     cond_mean_df_list[[i]] <- data.frame(tree_index = i - 1,
#                                          node = gsub(":[xy]", "", names(cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),])),
#                                          x = cond_mean_list[[i]][grep('x',rownames(cond_mean_list[[i]])),],
#                                          y = cond_mean_list[[i]][grep('y',rownames(cond_mean_list[[i]])),])
#   }
#   
#   cond_mean_df <- do.call(rbind, cond_mean_df_list)
#   
#   cond_mean_df$dist_to_true <- apply(cond_mean_df[,c('x', 'y')], 1, function(x) {
#     dist(rbind(x, indiv_node_info[indiv_node_info$n_node %in% MISSING_NODES,c('x', 'y')][1,]))
#   })
#   
#   cond_mean_df_combined <- cond_mean_df %>% 
#     #rename('tree_index' = tree_ind) %>% 
#     left_join(., tree_position_info[,c('tree_index', 'span')], 
#               by = 'tree_index')
#   
#   cor_val <- cor.test(cond_mean_df_combined$dist_to_true, cond_mean_df_combined$span, method = 'spearman')
#   
#   
#   span_dist_plot_list[[paste0('INDIV', INDIV)]] <- cond_mean_df_combined %>% 
#     ggplot() +
#     geom_point(aes(x = span, y = dist_to_true), alpha = 0.5) +
#     theme_bw() +
#     xlab('Span length') +
#     ylab('Distance to true location') +
#     ggtitle(paste0('INDIV = ',  INDIV, '; Spearman cor: ', round(cor_val$estimate, 3)))
#   
#   
#   span_color_map_plot_list[[paste0('INDIV', INDIV)]] <- cond_mean_df_combined %>% 
#     ggplot() +
#     geom_point(aes(x, y, color = span)) +
#     geom_point(data = indiv_node_info %>% 
#                  filter(!n_node %in% MISSING_NODES),
#                aes(x, y), color = 'blue', alpha = 0.25, size = 0.8) +
#     scale_color_gradient(low = '#b2d8b2', high = '#004000') +
#     geom_point(data = indiv_node_info %>% 
#                  filter(n_node %in% MISSING_NODES),
#                aes(x, y), shape = 21, color = 'purple', size = 2.5, stroke = 1.15) +
#     theme_bw() +
#     ggtitle(paste0('INDIV = ', INDIV) )
# }
# 
# multipan_different_posteriors <- cowplot::plot_grid(posterior_plot_list[[1]], posterior_plot_list[[2]], posterior_plot_list[[3]],
#                    dist_compare_boxplot_list[[1]], dist_compare_boxplot_list[[2]], dist_compare_boxplot_list[[3]],
#                    span_dist_plot_list[[1]], span_dist_plot_list[[2]], span_dist_plot_list[[3]],
#                    span_color_map_plot_list[[1]], span_color_map_plot_list[[2]], span_color_map_plot_list[[3]],
#                    ncol = 3)
# 
# ggsave('/Users/alexlewanski/Documents/michigan_state/research/tmp/multipan_different_posteriors.png', 
#        multipan_different_posteriors, 
#        width = 12*1.1, height = 11*1.1, units = "in")


