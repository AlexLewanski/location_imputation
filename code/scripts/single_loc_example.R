#################################################################
#################################################################
### EXAMPLE OF ESTIMATING LOCATIONS IN A NON-ADMIXED SCENARIO ###
#################################################################
#################################################################

########################
### LOADING PACKAGES ###
########################
library(here)
library(dplyr)
library(ggplot2)
library(ape)
#library(phytools)

source(here('code', 'scripts', 'location_est_code.R'))


#######################
### PROCESSING INFO ###
#######################

#true_rate_mat <- diag(4, nrow = 2)


files_vec <- list.files(here('simulation_output', 'test_notracking'), full.names = TRUE)
tree_list <- lapply(files_vec[grep(".*nwk$", files_vec)],
                    function(x) ape::read.tree(x)) 


indiv_node_info <- read.delim(here('simulation_output', 'test_notracking', 'indiv_node_info.txt'), header = TRUE)
indiv_node_info$n_node <- paste0('n', indiv_node_info$node)

#let's just use 75 trees for the example
tree_list_subset1 <- tree_list[1:75]



####################################################################
### DEMONSTRATION OF LOCATION ESTIMATION FOR A SINGLE INDIVIDUAL ###
####################################################################

### EXAMPLE FOR ONE INDIVIDUAL ###
#let's pick an individual to use as the unknown individual 
MISSING_INDIV <- 29
MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]

#the location information for the georeferenced individuals (unknown individual is removed)
miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES),]
rownames(miss_info) <- miss_info$n_node
miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])


### DERIVING DISTRIBUTION INFORMATION FROM THE TREES AND COORDINATE INFO ###
root_rate_info <- root_rate_multitree(tree = tree_list_subset1, miss_info_mat, MISSING_NODES)

multi_tree_test <- process_trees(
  tree = tree_list_subset1, 
  trait = miss_info_mat, 
  brown_rate_mat = root_rate_info$mean_brown_rate, 
  root_coords_list = root_rate_info$root_list,
  missing_inds = MISSING_NODES
)

# library(phytools)
# x <- setNames(miss_info_mat[,2], rownames(miss_info_mat))
# 
# plotTree.barplot(
#   tree_list_subset1[[1]],
#   x,
#   args.barplot = list(
#     col = "steelblue",
#     border = NA
#   )
# )



### EXAMPLE FOR ONE INDIVIDUAL: PRODUCT OF GAUSSIANS ###

#APPROACH 1: MULTIPLYING THE DISTRIBUTIONS AND THEN CONDITIONING
A <- matrix(c(
  1, 0, -1, 0,
  0, 1, 0, -1
), 2, 4, byrow = TRUE)

mvnorm_product <- multiply_mvnorms(lapply(multi_tree_test, function(x) x$cond_distr$condmean), 
                                   lapply(multi_tree_test, function(x) x$cond_distr$condvar))

mvnorm_condition <- condition_mvn(mu = mvnorm_product$mu, 
                                  V = mvnorm_product$Sigma, 
                                  A = A)



# APPROACH 2: CONDITIONING EACH DISTRIBUTION AND THEN MULTIPLYING ###

conditioned_tree_dist_list <- list()

for (i in seq_len(length(multi_tree_test))) {
  
  conditioned_tree_dist_list[[i]] <- condition_mvn(mu = multi_tree_test[[i]]$cond_distr$condmean, 
                                                   V = multi_tree_test[[i]]$cond_distr$condvar, 
                                                   A = A)
}

#I'm currently just taking a single x and y dimension from a single node
conditioned_mvnorm_product <- multiply_mvnorms(lapply(conditioned_tree_dist_list, function(x) x$mean[1:2,1]), 
                                   lapply(conditioned_tree_dist_list, function(x) x$covariance[1:2, 1:2]))



### EXAMPLE FOR ONE INDIVIDUAL: GENERALIZED LEAST SQUARES ###

# weighted_least_squares(create_design(mean_vec), 
#                        Sigma_list = multi_tree_test$cond_var, 
#                        mu = mean_vec
# )

weighted_least_squares(create_design(unlist(lapply(multi_tree_test, function(x) x$cond_distr$condmean))), 
                       Sigma_list = lapply(multi_tree_test, function(x) x$cond_distr$condvar), 
                       mu = unlist(lapply(multi_tree_test, function(x) x$cond_distr$condmean))
)



### EXAMPLE FOR ONE INDIVIDUAL: ESTIMATION VIA THE EM FUNCTIONS ###

cond_var_list <- lapply(multi_tree_test, function(x) x$cond_distr$condvar)
cond_mean_list <- lapply(multi_tree_test, function(x) x$cond_distr$condmean)

#hard assign EM
em_singleloc_hardassign <- em_hardassign(cond_mean = cond_mean_list, 
                                 cond_var = cond_var_list, 
                                 k = 1, 
                                 max_stp = 5, 
                                 conv_thresh = 1e-10, 
                                 sp_bounds = list(x = c(-5, 105), 
                                                  y = c(-5, 105)),
                                 starting_location_seed = NULL
)

#soft assign EM
em_singleloc_softassign <- em_softassign(cond_mean = cond_mean_list, 
                                         cond_var = cond_var_list, 
                                         k = 1, 
                                         max_stp = 5, 
                                         conv_thresh = 1e-10, 
                                         sp_bounds = list(x = c(-5, 105), 
                                                          y = c(-5, 105)),
                                         starting_location_seed = NULL
)

#same estimates as the other approaches
em_singleloc_hardassign$X
em_singleloc_softassign$X



#######################################################################
### DEMONSTRATION OF LOCATION ESTIMATION FOR A BUNCH OF INDIVIDUALS ###
#######################################################################

tree_count <- 50
indiv_count <- 75

tree_list_subset2 <- tree_list[1:tree_count]

indiv_vec <- unique(indiv_node_info$indiv)[1:indiv_count]

prep_list <- list()
for (MISSING_INDIV in indiv_vec) {
  
  MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]
  
  #the location information for the georeferenced individuals (unknown individual is removed)
  miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES),]
  rownames(miss_info) <- miss_info$n_node
  miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
  
  
  ### DERIVING DISTRIBUTION INFORMATION FROM THE TREES AND COORDINATE INFO ###
  root_rate_info <- root_rate_multitree(tree = tree_list_subset2, miss_info_mat, MISSING_NODES)
  
  prep_list[[paste0('indiv', MISSING_INDIV)]] <- process_trees(tree = tree_list_subset2, 
                                                               trait = miss_info_mat, 
                                                               brown_rate_mat = root_rate_info$mean_brown_rate, 
                                                               root_coords_list = root_rate_info$root_list,
                                                               missing_inds = MISSING_NODES
                                                               )
  message('finished ', MISSING_INDIV)
}


loc_est_df <- data.frame(x_est = rep(NA, length(prep_list)),
                         y_est = rep(NA, length(prep_list)))

for (i in seq_len(length(prep_list))) {
  
  loc_est <- weighted_least_squares(create_design(unlist(lapply(prep_list[[i]], function(x) x$cond_distr$condmean))), 
                                    Sigma_list = lapply(prep_list[[i]], function(x) x$cond_distr$condvar), 
                                    mu = unlist(lapply(prep_list[[i]], function(x) x$cond_distr$condmean))
  )
  
  loc_est_df[i,] <- t(loc_est)
  
}


loc_est_df$indiv <- as.numeric(gsub(pattern = 'indiv', "", names(prep_list)))

loc_est_df_new <- merge(loc_est_df, indiv_node_info[!duplicated(indiv_node_info$indiv),c('indiv', 'x', 'y')],
                        by = c('indiv'))

loc_est_df_new %>%
  ggplot() +
  geom_segment(aes(x = x_est, y = y_est, xend = x, yend = y), color = 'gray') +
  geom_point(data = tidyr::pivot_longer(loc_est_df_new,
                                        cols = c(x, x_est, y, y_est),
                                        names_to = c(".value", "type"),
                                        names_pattern = "(x|y)(_est)?") %>%
               mutate(type = if_else(type == '_est', 'estimate', 'true')),
  aes(x, y, color = type), size = 2.5) +
  geom_point(data = indiv_node_info[!duplicated(indiv_node_info$indiv),c('indiv', 'x', 'y')],
             aes(x = x, y =y), size = 0.5, color = 'black') +
  theme_bw()

                    

#######################
### CODE NOT IN USE ###
#######################

# prep_list <- list()
# for (MISSING_INDIV in indiv_vec) {
#   
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]
#   
#   #the location information for the georeferenced individuals (unknown individual is removed)
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES),]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
#   
#   
#   ### DERIVING DISTRIBUTION INFORMATION FROM THE TREES AND COORDINATE INFO ###
#   prep_list[[paste0('indiv', MISSING_INDIV)]] <- process_multitree(tree = tree_list_subset2,
#                                                                    trait_mat = miss_info_mat,
#                                                                    rate_mat = true_rate_mat,
#                                                                    missing_inds = MISSING_NODES,
#                                                                    extract_distribution_info = TRUE)
#   message('finished ', MISSING_INDIV)
# }



# 
# weighted_least_squares(create_design(mean_vec), 
#                        Sigma_list = multi_tree_test$cond_var, 
#                        mu = mean_vec
# )
# 
# 
# test_cond <- condition_mvn(mu = test$mu, Sigma = as.matrix(Matrix::forceSymmetric(test$Sigma)), 
#                            A, b)
# 
# 
# 
# 
# #preparation of material for 
# 
# 
# TREE <- 1
# MISSING_NODES
# true_rate_mat
# 
# 
# process_single_tree(tree = tree_list[[4]],
#                                 trait_mat = miss_info_mat,
#                                 rate_mat = true_rate_mat,
#                                 missing_inds = MISSING_NODES,
#                                 prior_variance = 500,
#                                 retain = c('tree_covar', 
#                                            'data_prep_output', 
#                                            'ls_root', 
#                                            'marginal_variance', 
#                                            'cond_var_elements',
#                                            'cond_mean',
#                                            'cond_var')[c(1, 3, 7)])
# 
# multi_tree_test <- process_multitree(tree = tree_list_subset[1:100],
#                    trait_mat = miss_info_mat,
#                    rate_mat = true_rate_mat,
#                    missing_inds = MISSING_NODES,
#                    prior_variance = 1000,
#                    extract_distribution_info = TRUE)
# 
# 
# condition_mvn_alt <- function(mu, Sigma, A, b) {
#   
#   ASAt <- A %*% Sigma %*% t(A)
#   K <- Sigma %*% t(A) %*% solve(ASAt)
#   mu_cond <- mu + K %*% (b - A %*% mu)
#   Sigma_cond <- Sigma - K %*% A %*% Sigma
#   
#   list(
#     mean = mu_cond,
#     covariance = Sigma_cond
#   )
# }
# 
# 
# 
# A <- matrix(c(
#   1, 0, -1, 0,
#   0, 1, 0, -1
# ), 2, 4, byrow = TRUE)
# 
# b <- c(0, 0)
# 
# 
# cond_var_list <- lapply(prep_list_missing_indiv[[1]]$cond_var_elements_list, function(x) x$cond_var)
# 
# test <- multiply_mvnorms(prep_list_missing_indiv[[1]]$cond_mean_list[seq(1, 300, 10)], cond_var_list[seq(1, 300, 10)])
# 
# test_cond <- condition_mvn(mu = test$mu, Sigma = as.matrix(Matrix::forceSymmetric(test$Sigma)), 
#                            A, b)
# 
# condition_mvn_alt(mu = test$mu, 
#               Sigma = as.matrix(Matrix::forceSymmetric(test$Sigma)), 
#               A = A,
#               b = 0)
# 
# 
# condition_mvn(mu = test$mu, 
#               V = as.matrix(Matrix::forceSymmetric(test$Sigma)), 
#               A = A)
# 
# 
# #calc_covar
# tree_covar <- calc_vcv_single_tree(tree_list_subset[[1]])
# 
# 
# #next processing steps
# data_prep_output <- data_prep(trait = miss_info_mat,
#                               phylo_covar = tree_covar,
#                               rate_mat = true_rate_mat,
#                               missing_inds = MISSING_NODES)
# 
# 
# #calculate least squares root value for prior (or just create independent prior)
# 
# ls_root <- calc_root_coords(missing_inds = MISSING_NODES, 
#                             kron_mat = data_prep_output$kron_varcovar_ratemat, 
#                             trait_vector = data_prep_output$vector_trait)
# 
# 
# marginal_variance <- add_variance(kron_mat = data_prep_output$kron_varcovar_ratemat, var = 10000)
# 
# 
# cond_var_elements <- calc_cond_var_elements(
#   kron_varcovar_ratemat = marginal_variance,
#   no_ind = data_prep_output$no_ind,
#   yes_ind = data_prep_output$yes_ind
#   )
# 
# 
# #calculate conditional mean
# cond_mean <- calc_cond_mean(covar_12 = cond_var_elements$covar_12,
#                inv_covar22 = cond_var_elements$inv_covar_22,
#                trait_vec = data_prep_output$vector_trait,
#                means = ls_root, #currently set at least squares estimate of root coordinates
#                no_ind = data_prep_output$no_ind,
#                yes_ind = data_prep_output$yes_ind)
# 
# 
# 
# cond_var_list <- lapply(prep_list_missing_indiv[[1]]$cond_var_elements_list, function(x) x$cond_var)
# 
# 
# 
# #calculate conditional mean
# calc_cond_mean(covar_12 = cond_var_elements$covar_12,
#                inv_covar22 = cond_var_elements$inv_covar_22,
#                trait_vec,
#                means = ls_root,
#                no_ind = data_prep_output$no_ind,
#                yes_ind = data_prep_output$yes_ind)
# 
# 
# cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                       prep_list = data_prep_list_subset,
#                                       tree_means = expected_root_vec)
# 
# 
# 
# 
# 
# #integrative out root value to get marginal distribution of unknown node coordinates 
# 
# 
# 
# #condition on known values
# 
# 
# 
# #condition on node location constraints
# 
# 
# 
# #####################
# #####################
# #####################
# 
# vcv_list_subset <- calc_vcv_multitree(tree_list = tree_list_subset)
# 
# data_prep_list_subset <- data_prep_multitree(vcv_list_subset,
#                                              trait = miss_info_mat,
#                                              rate_mat = true_rate_mat,
#                                              missing_inds = MISSING_NODES)
# 
# cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list_subset)
# 
# 
# true_rate_mat <- diag(3, nrow = 2)
# 
# 
# cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                       prep_list = data_prep_list_subset,
#                                       tree_means = tree_mean_vec)
# 
# 
# data_prep_info <- list()
# 
# for (i in 1:800) {
#   data_prep_info[[i]] <- data_prep(trait = miss_info_mat,
#                                    phylo_covar = vcv_list_subset[[i]],
#                                    rate_mat = true_rate_mat,
#                                    missing_inds = MISSING_NODES)
# }
# 
# data_prep_list <- data_prep(trait = miss_info_mat,
#                             phylo_covar = vcv_list_subset[[1]],
#                             rate_mat = true_rate_mat,
#                             missing_inds = MISSING_NODES)
# 
# 
# data_prep_list$kron_varcovar_ratemat + diag(700, nrow = nrow(data_prep_list$kron_varcovar_ratemat))
# 
# 
# rownames(data_prep_list$kron_varcovar_ratemat)
# names(data_prep_list$vector_trait)
# 
# 
# length(data_prep_list)
# 
# 
# 
# #plot(
# #  do.call(rbind, lapply(expected_root_list, function(x) t(x))),
# #  xlim = c(0, 100), ylim = c(0, 100)
# #)
# 
# 
# 
# 
# vcv_subset <- data_prep_list$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(data_prep_list$kron_varcovar_ratemat)) %in% MISSING_NODES),
#                                      !(gsub(":[xy]", "", colnames(data_prep_list$kron_varcovar_ratemat)) %in% MISSING_NODES)]
# 
# design_mat <- kronecker(rep(1, nrow(data_prep_list$trait_mat)), diag(2))
# 
# solve(
#   t(design_mat) %*% solve(vcv_subset) %*% design_mat
# ) %*%
#   t(design_mat) %*%
#   solve(vcv_subset) %*%
#   data_prep_list$vector_trait
# 
# 
# vcv_subset <- data_prep_list$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(data_prep_list$kron_varcovar_ratemat)) %in% MISSING_NODES),
#                                                    !(gsub(":[xy]", "", colnames(data_prep_list$kron_varcovar_ratemat)) %in% MISSING_NODES)]
# 
# 
# 
# design_mat <- kronecker(rep(1, nrow(Y)), diag(2))
# 
# 
# rate_mat <- diag(2)
# dimnames(rate_mat) <- list(c('x', 'y'), c('x', 'y'))
# phylo_tree <- vcv.phylo(tree)
# k_phylo_tree <- kronecker(phylo_tree,rate_mat)
# 
# solve(
#   t(design_mat) %*% solve(k_phylo_tree) %*% design_mat
# ) %*%
#   t(design_mat) %*%
#   solve(k_phylo_tree) %*%
#   c(t(Y))
# 
# 
# vcv_subset <- data_prep_list$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(data_prep_list$kron_varcovar_ratemat)) %in% MISSING_NODES),
#                                                    !(gsub(":[xy]", "", colnames(data_prep_list$kron_varcovar_ratemat)) %in% MISSING_NODES)]
# 
# design_mat <- kronecker(rep(1, nrow(data_prep_list$trait_mat)), diag(2))
# 
# vcv.phylo(tree)
# 
# solve(
#   t(design_mat) %*% solve(vcv_subset) %*% design_mat
# ) %*%
#   t(design_mat) %*%
#   solve(vcv_subset) %*%
#   data_prep_list$vector_trait
# 
# 
# 
# identical(rownames(vcv_subset), rownames(data_prep_list$trait_mat))
# 
# one <- rep(1, nrow(data_prep_list$trait_mat))
# 
# t(one) %*% solve(vcv_subset) %*% data_prep_list$trait_mat %*% solve(t(one) %*% solve(vcv_subset) %*% one)
# 
# 
# 
# estimate_root_BM <- function(tree, x) {
#   
#   # Ensure the trait vector is ordered like the tree tips
#   x <- x[tree$tip.label]
#   
#   # Brownian covariance matrix (shared branch lengths)
#   C <- vcv.phylo(tree)
#   
#   # Vector of ones
#   one <- rep(1, length(x))
#   
#   # Generalized least squares estimate
#   root <- as.numeric(
#     t(one) %*%
#       solve(C) %*%
#       x /
#       (t(one) %*%
#          solve(C) %*%
#          one)
#   )
#   
#   return(root)
# }
# 
# 
# 
# library(ape)
# 
# set.seed(1)
# 
# tree <- rtree(20)
# 
# x <- rTraitCont(tree, sigma = 1)
# estimate_root_BM(tree, x)
# 
# 
# estimate_root_BM(tree, x)
# as.numeric(
#   t(one) %*%
#     solve(vcv_subset) %*%
#     data_prep_list$trait_mat /
#     (t(one) %*%
#        solve(vcv_subset) %*%
#        one)
# )
# 
# 
# data_prep_list$kron_varcovar_ratemat
# 
# data_prep_list$vector_trait
# 
# 
# data_prep_list$trait_mat
# 
# as.numeric(
#   t(one) %*%
#     solve(C) %*%
#     x /
#     (t(one) %*%
#        solve(C) %*%
#        one)
# )
# 
# 
# estimate_root_BM <- function(tree, x) {
#   
#   # Ensure the trait vector is ordered like the tree tips
#   x <- x[tree$tip.label]
#   
#   # Brownian covariance matrix (shared branch lengths)
#   C <- vcv.phylo(tree)
#   
#   # Vector of ones
#   one <- rep(1, length(x))
#   
#   # Generalized least squares estimate
#   root <- as.numeric(
#     t(one) %*%
#       solve(C) %*%
#       x /
#       (t(one) %*%
#          solve(C) %*%
#          one)
#   )
#   
#   return(root)
# }
# 
# 
# for (i in 1:800) {
#   
#   vcv_subset <- data_prep_info[[i]]$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(data_prep_info[[i]]$kron_varcovar_ratemat)) %in% MISSING_NODES),
#                                                           !(gsub(":[xy]", "", colnames(data_prep_info[[i]]$kron_varcovar_ratemat)) %in% MISSING_NODES)]
#   
#   expected_root_list[[i]] <- expected_root_state(vcv_subset, 
#                                                  data_prep_info[[i]]$vector_trait)
# }
# 
# 
# tree_list_subset1 <- tree_list[1:500]
# 
# missing_vec <- 50
# MISSING_INDIV <- 50
# #inferring 80
# #removing the missing individual's nodes from the spatial information
# prep_list_missing_indiv <- list()
# for (MISSING_INDIV in missing_vec) {
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]
#   #REMOVE_NODES <- indiv_node_info$n_node[indiv_node_info$indiv %in% admix_node_info$indiv & indiv_node_info$n_node != MISSING_NODES]
#   REMOVE_NODES <- 'NONE'
#   
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, REMOVE_NODES),]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
#   
#   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
#   
#   # tree_list_chopped_subset <- lapply(tree_list_subset, function(x, time) {
#   #   if (max(phytools::nodeHeights(x)) > time) {
#   #     return(treesliceR::squeeze_root(tree = x, time = time))
#   #   } else {
#   #     return(x)
#   #   }
#   # }, time = chop_time)
#   
#   
#   #tree_mean_vec <- c(t(tree_position_info[,c('x', 'y')]))
#   
#   tree_mean_vec <- rep(50, length(tree_list_subset)*2)
#   
#   #test <- c(t(matrix(data = c(1, 2, 3, 4), ncol = 2)))
#   
#   
#   ### full ###
#   vcv_list_subset <- calc_vcv_multitree(tree_list = tree_list_subset)
#   
#   data_prep_list_subset <- data_prep_multitree(vcv_list_subset,
#                                                trait = miss_info_mat,
#                                                rate_mat = true_rate_mat,
#                                                missing_inds = MISSING_NODES)
#   
#   cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list_subset)
#   
#   
#   cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                         prep_list = data_prep_list_subset,
#                                         tree_means = tree_mean_vec)
#   
#   
#   ### chopped ###
#   #vcv_list_chopped_subset <- calc_vcv_multitree(tree_list = tree_list_chopped_subset)
#   
#   #data_prep_list_chopped_subset <- data_prep_multitree(vcv_list_chopped_subset,
#   #                                                     trait = miss_info_mat, 
#   #                                                     rate_mat = true_rate_mat, 
#   #                                                     missing_inds = MISSING_NODES)
#   
#   #cond_var_elements_list_chopped <- calc_cond_var_elements_multitree(prep_list = data_prep_list_chopped_subset)
#   
#   #cond_mean_list_chopped <- calc_cond_mean_list(cond_param_list = cond_var_elements_list_chopped,
#   #                                              prep_list = data_prep_list_chopped_subset,
#   #                                              tree_means = tree_mean_vec)
#   
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['MISSING_NODES']] <- MISSING_NODES
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['REMOVE_NODES']] <- REMOVE_NODES
#   #prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['vcv_list']] <- vcv_list_chopped_subset
#   #prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['data_prep_list']] <- data_prep_list_chopped_subset
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['miss_info_mat']] <- miss_info_mat
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['tree_mean_vec']] <- tree_mean_vec
#   #prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_var_elements_list_chopped']] <- cond_var_elements_list_chopped
#   #prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_mean_list_chopped']] <- cond_mean_list_chopped
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_var_elements_list']] <- cond_var_elements_list
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_mean_list']] <- cond_mean_list
# }
# 
# 
# 
# expected_root_list
# 
# 
# 
# tree_list_subset1 <- tree_list[1:300]
# 
# missing_vec <- 49
# MISSING_INDIV <- 49
# #inferring 80
# #removing the missing individual's nodes from the spatial information
# prep_list_missing_indiv <- list()
# for (MISSING_INDIV in missing_vec) {
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]
#   REMOVE_NODES <- 'NONE'
#   
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, REMOVE_NODES),]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
#   
#   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
#   
#   
#   #tree_mean_vec <- rep(50, length(tree_list_subset)*2)
#   
#   
#   ### full ###
#   vcv_list_subset <- calc_vcv_multitree(tree_list = tree_list_subset)
#   message('here')
#   data_prep_list_subset <- data_prep_multitree(vcv_list_subset,
#                                                trait = miss_info_mat,
#                                                rate_mat = true_rate_mat,
#                                                missing_inds = MISSING_NODES)
#   
#   #cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list_subset)
#   
#   cond_var_elements_list <- calc_cond_var_elements_multitree_add_variance(prep_list = data_prep_list_subset,
#                                                                           variance = 700)
#   message('here')
#   #cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#   #                                      prep_list = data_prep_list_subset,
#   #                                      tree_means = tree_mean_vec)
#   
#   for (i in 1:length(data_prep_list_subset)) {
#     vcv_subset <- data_prep_list_subset[[i]]$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(data_prep_list_subset[[i]]$kron_varcovar_ratemat)) %in% MISSING_NODES),
#                                                                    !(gsub(":[xy]", "", colnames(data_prep_list_subset[[i]]$kron_varcovar_ratemat)) %in% MISSING_NODES)]
#     
#     expected_root_list[[i]] <- t(expected_root_state(vcv_subset, 
#                                                      data_prep_list_subset[[i]]$vector_trait))
#   }
#   message('here')
#   expected_root_vec <- c(t(do.call(rbind,expected_root_list)))
#   
#   cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                         prep_list = data_prep_list_subset,
#                                         tree_means = expected_root_vec)
#   
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['MISSING_NODES']] <- MISSING_NODES
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['REMOVE_NODES']] <- REMOVE_NODES
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['miss_info_mat']] <- miss_info_mat
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['tree_mean_vec']] <- expected_root_vec
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_var_elements_list']] <- cond_var_elements_list
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_mean_list']] <- cond_mean_list
#   
# }
# 
# 
# 
# tree_list_subset1 <- tree_list[1:300]
# 
# missing_vec <- 8
# MISSING_INDIV <- 8
# #inferring 80
# #removing the missing individual's nodes from the spatial information
# prep_list_missing_indiv <- list()
# for (MISSING_INDIV in missing_vec) {
#   MISSING_NODES <- indiv_node_info$n_node[indiv_node_info$indiv == MISSING_INDIV]
#   REMOVE_NODES <- 'NONE'
#   
#   miss_info <- indiv_node_info[!indiv_node_info$n_node %in% c(MISSING_NODES, REMOVE_NODES),]
#   rownames(miss_info) <- miss_info$n_node
#   miss_info_mat <- as.matrix(miss_info[,c('x', 'y')])
#   
#   tree_list_subset <- lapply(tree_list_subset1, function(x) drop.tip(x, REMOVE_NODES))
#   
#   
#   tree_mean_vec <- rep(50, length(tree_list_subset)*2)
#   
#   
#   ### full ###
#   vcv_list_subset <- calc_vcv_multitree(tree_list = tree_list_subset)
#   
#   data_prep_list_subset <- data_prep_multitree(vcv_list_subset,
#                                                trait = miss_info_mat,
#                                                rate_mat = true_rate_mat,
#                                                missing_inds = MISSING_NODES)
#   
#   #cond_var_elements_list <- calc_cond_var_elements_multitree(prep_list = data_prep_list_subset)
#   
#   cond_var_elements_list <- calc_cond_var_elements_multitree_add_variance(prep_list = data_prep_list_subset,
#                                                                           variance = 0)
#   
#   cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                        prep_list = data_prep_list_subset,
#                                        tree_means = tree_mean_vec)
#   # 
#   # for (i in 1:length(data_prep_list_subset)) {
#   #   vcv_subset <- data_prep_list_subset[[i]]$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(data_prep_list_subset[[i]]$kron_varcovar_ratemat)) %in% MISSING_NODES),
#   #                                                                  !(gsub(":[xy]", "", colnames(data_prep_list_subset[[i]]$kron_varcovar_ratemat)) %in% MISSING_NODES)]
#   #   
#   #   expected_root_list[[i]] <- t(expected_root_state(vcv_subset, 
#   #                                                    data_prep_list_subset[[i]]$vector_trait))
#   # }
#   # 
#   # expected_root_vec <- c(t(do.call(rbind,expected_root_list)))
#   # 
#   # cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#   #                                       prep_list = data_prep_list_subset,
#   #                                       tree_means = expected_root_vec)
#   
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['MISSING_NODES']] <- MISSING_NODES
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['REMOVE_NODES']] <- REMOVE_NODES
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['miss_info_mat']] <- miss_info_mat
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['tree_mean_vec']] <- tree_mean_vec
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_var_elements_list']] <- cond_var_elements_list
#   prep_list_missing_indiv[[paste0('indiv', MISSING_INDIV)]][['cond_mean_list']] <- cond_mean_list
# }
# 
# 
# 
# 
# # cond_mean_list <- list()
# # for (INDIV in names(prep_list_missing_indiv)) {
# #   split_condmean_list <- lapply(prep_list_missing_indiv[[INDIV]]$cond_mean_list, function(x) {
# #     wide_df <- as.data.frame(x) %>% 
# #       tibble::rownames_to_column() %>%
# #       tidyr::separate(rowname, into = c("id", "var"), sep = ":") %>%
# #       tidyr::pivot_wider(names_from = var, values_from = V1)
# #     
# #     split(wide_df, wide_df$id)
# #   })
# #   
# #   cond_mean_vec_list <- list()
# #   for (NODE in names(split_condmean_list[[1]])) {
# #     cond_mean_vec_list[[NODE]] <- do.call(rbind, lapply(split_condmean_list, function(x, i) x[[i]], i = NODE)) %>% 
# #       select(-1) %>% 
# #       as.matrix()
# #   }
# #   
# #   cond_mean_list[[INDIV]] <- cond_mean_vec_list
# # }
# 
# 
# 
# # cond_var_list <- list()
# # for (INDIV in names(prep_list_missing_indiv)) {
# #   
# #   for (NODE in prep_list_missing_indiv[[INDIV]]$MISSING_NODES) {
# #     cond_var_vec <- vector(mode = 'numeric', length = length(prep_list_missing_indiv[[INDIV]]$cond_var_elements_list))
# #     for (tree_ind in 1:length(prep_list_missing_indiv[[INDIV]]$cond_var_elements_list)) {
# #       cond_var_mat <- prep_list_missing_indiv[[INDIV]]$cond_var_elements_list[[tree_ind]]$cond_var
# #       cond_var_vec[tree_ind] <- cond_var_mat[rownames(cond_var_mat) == paste0(NODE, ':x'),colnames(cond_var_mat) == paste0(NODE, ':x')]
# #     }
# #     cond_var_list[[INDIV]][[NODE]] <- cond_var_vec
# #   }
# # }
# 
# 
# # var_mat <- kronecker(rep(1, nrow(data_prep_list_subset[[1]]$trait_mat))%*%t(rep(1, nrow(data_prep_list_subset[[1]]$trait_mat))), diag(2)*700)
# # 
# # cond_var_elements_list <- calc_cond_var_elements_multitree_add_variance(prep_list = data_prep_list_subset,
# #                                                                         variance = 700)
# # 
# # 
# # cond_var_elements_list <- calc_cond_var_elements_multitree_add_variance(prep_list = data_prep_list_subset,
# #                                                                         variance = 700)
# # 
# 
# 
# A <- matrix(c(
#   1, 0, -1, 0,
#   0, 1, 0, -1
# ), 2, 4, byrow = TRUE)
# 
# b <- c(0, 0)
# 
# 
# cond_var_list <- lapply(prep_list_missing_indiv[[1]]$cond_var_elements_list, function(x) x$cond_var)
# 
# test <- multiply_mvnorms(prep_list_missing_indiv[[1]]$cond_mean_list[seq(1, 300, 10)], cond_var_list[seq(1, 300, 10)])
# 
# test_cond <- condition_mvn(mu = test$mu, Sigma = as.matrix(Matrix::forceSymmetric(test$Sigma)), 
#                            A, b)
# 
# draws_cond <- mvtnorm::rmvnorm(
#   n = 1000,
#   mean = test_cond$mean[1:2],
#   sigma = as.matrix(Matrix::forceSymmetric(test_cond$covariance[1:2, 1:2]))
# )
# 
# 
# plot(draws_cond, xlim = c(0, 100), ylim = c(0,100), cex = 0.2, col = 'purple')
# points(indiv_node_info[indiv_node_info$n_node %in% prep_list_missing_indiv[[1]]$MISSING_NODES,][1, c('x', 'y')], 
#        col = 'black', pch = 19)
# points(indiv_node_info[!indiv_node_info$n_node %in% prep_list_missing_indiv[[1]]$MISSING_NODES,c('x', 'y')], 
#        col = 'gray', pch = 8)
# 
# 
# 
# draws_cond <- mvtnorm::rmvnorm(
#   n = 1000,
#   mean = test$mu,
#   sigma = as.matrix(Matrix::forceSymmetric(test$Sigma))
# )
# 
# plot(draws_cond[,c(1,2)], xlim = c(0, 100), ylim = c(0,100), cex = 0.2, col = 'purple')
# points(draws_cond[,c(3,4)], xlim = c(0, 100), ylim = c(0,100), cex = 0.2, col = 'blue')
# points(indiv_node_info[!indiv_node_info$n_node %in% prep_list_missing_indiv[[1]]$MISSING_NODES,c('x', 'y')], 
#        col = 'gray', pch = 8)
# points(indiv_node_info[indiv_node_info$n_node %in% prep_list_missing_indiv[[1]]$MISSING_NODES,][1, c('x', 'y')], 
#        col = 'black', pch = 19)
# 
# 
# 
# fit_list <- list()
# 
# for (i in 1:1000) {
#   
#   tree_inds <- sample(1:300, 1, replace = TRUE)
#   test <- multiply_mvnorms(prep_list_missing_indiv[[1]]$cond_mean_list[tree_inds], cond_var_list[tree_inds])
#   test_cond <- condition_mvn(mu = test$mu, Sigma = as.matrix(Matrix::forceSymmetric(test$Sigma)), 
#                              A, b)
#   fit_list[[i]] <- test_cond$mean[1:2]
#   message(i)
# }
# 
# 
# plot(do.call(rbind,fit_list), xlim = c(0, 100), ylim = c(0,100), cex = 0.2, col = 'purple')
# points(indiv_node_info[indiv_node_info$n_node %in% prep_list_missing_indiv[[1]]$MISSING_NODES,][1, c('x', 'y')], 
#        col = 'black', pch = 19)
# points(indiv_node_info[!indiv_node_info$n_node %in% prep_list_missing_indiv[[1]]$MISSING_NODES,c('x', 'y')], 
#        col = 'gray', pch = 8)
# 
# 
# 
# 
# 
# # $mean
# # n193:x   n193:y   n192:x   n192:y 
# # 32.57253 81.34221 32.57253 81.34221 
# # 
# # $covariance
# # n193:x     n193:y     n192:x     n192:y
# # n193:x 0.05680363 0.00000000 0.05680363 0.00000000
# # n193:y 0.00000000 0.05680363 0.00000000 0.05680363
# # n192:x 0.05680363 0.00000000 0.05680363 0.00000000
# # n192:y 0.00000000 0.05680363 0.00000000 0.05680363
# 
# 
# prep_list_missing_indiv[[1]]$cond_mean_list[seq(1, 100, 5)]
# 
# 
# average_first_tmrca_age[average_first_tmrca_age$indiv == MISSING_INDIV,]
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
# cond_var_elements_list[[1]]$covar_11
# cond_var_elements_list0[[1]]$covar_12
# cond_var_elements_list700[[1]]$covar_12
# 
# data_prep_list_subset[[6]]$kron_varcovar_ratemat
# 
# 
# cond_var_elements_list[[9]]$cond_var
# 
# cond_var_elements_list0[[9]]$cond_var[3,3]
# cond_var_elements_list700[[9]]$cond_var[3,3]
# 
# #cond_var_elements_list1 <- calc_cond_var_elements_multitree(prep_list = data_prep_list_subset)
# 
# cond_var_elements_list[[5]]
# 
# cond_var_elements_list[[20]]$cond_var
# cond_var_elements_list1[[20]]$cond_var
#cond_mean_list <- calc_cond_mean_list(cond_param_list = cond_var_elements_list,
#                                      prep_list = data_prep_list_subset,
#                                      tree_means = tree_mean_vec)




# calc_cond_var_elements_multitree <- function(prep_list, variance = 0) {
#   lapply(prep_list, function(x, var) {
#     tip_count <- nrow(x$kron_varcovar_ratemat)/2
#     #var_mat <- kronecker(rep(1, tip_count)%*%t(rep(1, tip_count)), diag(2)*var)
#     #return(list(kron = x$kron_varcovar_ratemat,
#     #            var_mat =  var_mat))
#     calc_cond_var_elements(x$kron_varcovar_ratemat, 
#                            x$no_ind, x$yes_ind)
#   }, var = variance)
# }
# 
# 
# calc_vcv_single_tree <- function(tree) {
#   
#   tip_order <- tree$tip.label
#   vcv_mat <- vcv.phylo(tree)
#   
#   return(vcv_mat[match(tip_order, rownames(vcv_mat)),
#                  match(tip_order, colnames(vcv_mat))]
#   )
# }
# 
# 
# calc_vcv_multitree <- function(tree_list) {
#   tip_order <- tree_list[[1]]$tip.label
#   return(
#     lapply(tree_list, function(x, tip_order) {
#       vcv_mat <- vcv.phylo(x)
#       return(vcv_mat[match(tip_order, rownames(vcv_mat)),
#                      match(tip_order, colnames(vcv_mat))]
#       )
#     }, tip_order = tip_order)
#   )
# }


# data_prep_multitree <- function(vcv_list, trait, rate_mat, missing_inds) {
#   lapply(vcv_list, function(x, trait, rate_mat, miss_inds) {
#     data_prep(trait = trait,
#               phylo_covar = x,
#               rate_mat = rate_mat,
#               missing_inds = miss_inds)
#   }, trait = trait, rate_mat = rate_mat, miss_inds = missing_inds)
# } 


# data_prep <- function(trait,
#                       phylo_covar,
#                       rate_mat,
#                       missing_inds) {
#   
#   dimnames(rate_mat) <- list(c('x', 'y'), c('x', 'y'))
#   
#   trait_reorder <- trait[rank(match(rownames(phylo_covar), rownames(trait)), na.last = NA),]
#   
#   kron_varcovar_ratemat <- kronecker(phylo_covar, rate_mat, make.dimnames = TRUE)
#   
#   yes_ind <- which(!gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
#   no_ind <- which(gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
#   
#   
#   return(
#     list(kron_varcovar_ratemat = kron_varcovar_ratemat,
#          trait_mat = trait_reorder,
#          vector_trait = setNames(c(t(trait_reorder)), nm = rep(rownames(trait_reorder), each = 2)) ,
#          yes_ind = yes_ind,
#          no_ind = no_ind)
#   )
# }


# calc_cond_var_elements <- function(kron_varcovar_ratemat, 
#                                    no_ind, 
#                                    yes_ind) {
#   
#   covar_11 <- kron_varcovar_ratemat[no_ind,no_ind]
#   covar_12 <- kron_varcovar_ratemat[no_ind,yes_ind]
#   inv_covar_22 <- solve(kron_varcovar_ratemat[yes_ind,yes_ind])
#   covar_21 <- kron_varcovar_ratemat[yes_ind,no_ind]
#   
#   return(list(
#     covar_11 = covar_11,
#     covar_12 = covar_12,
#     inv_covar_22 = inv_covar_22,
#     covar_21 = covar_21,
#     cond_var = covar_11 - covar_12%*%inv_covar_22%*%covar_21
#   ))
#   
# }

# calc_cond_var_elements_multitree <- function(prep_list) {
#   lapply(prep_list, function(x) {
#     calc_cond_var_elements(x$kron_varcovar_ratemat, x$no_ind, x$yes_ind)
#   })
# }


# calc_cond_mean <- function(covar_12, inv_covar22, 
#                            trait_vec, means, 
#                            no_ind, yes_ind) {
#   mean1 <- rep(means, times = length(no_ind)/2)
#   mean2 <- rep(means, times = length(yes_ind)/2)
#   return(
#     mean1 + (covar_12 %*% inv_covar22 %*% (trait_vec -  mean2))
#   )
# }
# 
# calc_cond_mean_list <- function(cond_param_list,
#                                 prep_list,
#                                 tree_means) {
#   cond_mean_list <- list()
#   for (i in 1:length(cond_param_list)) {
#     cond_mean_list[[i]] <- calc_cond_mean(covar_12 = cond_param_list[[i]]$covar_12, 
#                                           inv_covar22 = cond_param_list[[i]]$inv_covar_22, 
#                                           trait_vec = prep_list[[i]]$vector_trait, 
#                                           means = tree_means[c(i*2 - 1, i*2)], 
#                                           no_ind = prep_list[[i]]$no_ind, 
#                                           yes_ind = prep_list[[i]]$yes_ind)
#   }
#   
#   return(cond_mean_list)
# }


# data_reorder <- function(trait,
#                          tree_covar#,
#                          #rate_mat,
#                          #missing_inds
# ) {
#   
#   
#   #dimnames(rate_mat) <- list(c('x', 'y'), c('x', 'y'))
#   
#   #trait_reorder <- trait[match(rownames(trait), rownames(phylo_covar)),]
#   trait_reorder <- trait[rank(match(rownames(tree_covar), rownames(trait)), na.last = NA),]
#   
#   #kron_varcovar_ratemat <- kronecker(phylo_covar, rate_mat, make.dimnames = TRUE)
#   #yes_ind <- which(!gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
#   #no_ind <- which(gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
#   
#   return(
#     list(trait_mat = trait_reorder,
#          vector_trait = c(t(trait_reorder))
#     )
#   )
# }

# expected_root_state <- function(vcv, trait_vec) {
#   #recover()
#   design_mat <- kronecker(rep(1, length(trait_vec)/2), diag(2))
#   vcv_inv <- solve(vcv)
#   
#   return(
#     solve(
#       t(design_mat) %*% vcv_inv %*% design_mat
#     ) %*% t(design_mat) %*% vcv_inv %*% trait_vec
#   )
#   
# }

# tree_list_subset1 <- tree_list[1:75]
# 
# 
# multi_tree_test <- process_multitree(tree = tree_list_subset1[1:10],
#                                      trait_mat = miss_info_mat,
#                                      rate_mat = true_rate_mat,
#                                      missing_inds = MISSING_NODES,
#                                      extract_distribution_info = FALSE)
# 
# 
# vcv_subset <- multi_tree_test[[1]]$data_prep_output$kron_varcovar_ratemat[!(gsub(":[xy]", "", rownames(multi_tree_test[[1]]$data_prep_output$kron_varcovar_ratemat)) %in% MISSING_NODES),
#                                                                           !(gsub(":[xy]", "", colnames(multi_tree_test[[1]]$data_prep_output$kron_varcovar_ratemat)) %in% MISSING_NODES)]
# 
# root_state <- expected_root_state(vcv_subset, 
#                                   multi_tree_test[[1]]$data_prep_output$vector_trait)
# 
# 
# covar_subset <- multi_tree_test[[1]]$tree_covar[!rownames(multi_tree_test[[1]]$tree_covar) %in% MISSING_NODES,
#                                                 !colnames(multi_tree_test[[1]]$tree_covar) %in% MISSING_NODES]
# 
# 
# 
# term1 <- multi_tree_test[[1]]$data_prep_output$trait_mat - one_col_vec%*%t(root_state)
# 
# t(term1)%*%solve(covar_subset)%*%term1/nrow( multi_tree_test[[1]]$data_prep_output$trait_mat)
# 
# 
# root_state <- calc_root_coords(missing_inds = MISSING_NODES, 
#                                tree_covar = multi_tree_test[[1]]$tree_covar, 
#                                trait_mat = multi_tree_test[[1]]$data_prep_output$trait_mat)


#est_brownian_rate_mat <- function(tree_covar, root_vals, trait_mat) {
#  term1 <- trait_mat - one_col_vec%*%t(root_vals)
#  return(t(term1)%*%solve(tree_covar)%*%term1/nrow(trait_mat))
#}
# multi_tree_test <- process_multitree(tree = tree_list_subset1,
#                                      trait_mat = miss_info_mat,
#                                      rate_mat = true_rate_mat,
#                                      missing_inds = MISSING_NODES,
#                                      extract_distribution_info = TRUE)