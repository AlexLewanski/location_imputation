#########################################################
#########################################################
### FUNCTIONS FOR SPATIAL IMPUTATION WITH GENEALOGIES ###
#########################################################
#########################################################
message('REQUIRED PACKAGES: MASS, Matrix, mvtnorm')



###################################################
### ROOT AND BROWNIAN RATE ESTIMATION FUNCTIONS ###
###################################################

expected_root_state <- function(vcv, trait_mat) {
  #recover()
  
  one_col_vec <- matrix(rep(1, nrow(trait_mat)), ncol = 1)
  
  vcv_inv <- solve(vcv)
  
  return(
    t(solve(t(one_col_vec)%*%vcv_inv%*%one_col_vec)%*%t(one_col_vec)%*%vcv_inv%*%trait_mat)
  )
  
}

calc_root_coords <- function(missing_inds, tree_covar, trait_mat) {
  vcv_subset <- tree_covar[!(rownames(tree_covar) %in% missing_inds),
                           !(colnames(tree_covar) %in% missing_inds)]
  
  root_state <- expected_root_state(vcv_subset, 
                                    trait_mat)
  
  return(root_state)
}



est_brownian_rate_mat <- function(missing_inds = NULL, tree_covar, root_vals, trait_mat) {
  #recover()
  vcv_subset <- tree_covar[!(rownames(tree_covar) %in% missing_inds),
                           !(colnames(tree_covar) %in% missing_inds)]
  one_col_vec <- matrix(rep(1, nrow(trait_mat)), 
                        ncol = 1)
  
  term1 <- trait_mat - one_col_vec%*%t(root_vals)
  return(t(term1)%*%solve(vcv_subset)%*%term1/nrow(trait_mat))
}


estimate_root_and_rate <- function(tree, missing_inds, trait_info) {
  #recover()
  tree_covar <- vcv.phylo(tree)
  vcv_subset <- tree_covar[!(rownames(tree_covar) %in% missing_inds),
                           !(colnames(tree_covar) %in% missing_inds)]
  trait_reorder <- trait_info[rank(match(rownames(tree_covar), rownames(trait_info)), na.last = NA),]
  
  
  root_coords <- calc_root_coords(missing_inds = NULL, 
                                  tree_covar = vcv_subset, 
                                  trait_mat = trait_reorder)
  
  brown_rate_mat <- est_brownian_rate_mat(missing_inds = NULL,
                                          tree_covar = vcv_subset, 
                                          root_vals = root_coords, 
                                          trait_mat = trait_reorder)
  
  return(
    list(
      root_coords = root_coords,
      brown_rate_mat = brown_rate_mat
    )
  )
}


root_rate_multitree <- function(tree_list, trait, missing_inds) {
  #recover()
  info_list <- lapply(tree_list, function(x, missing, trait) {
    estimate_root_and_rate(tree = x, 
                           missing_inds = missing, 
                           trait_info = trait)
    
  }, missing = missing_inds, trait = trait)
  
  
  return(
    list(
      root_list = lapply(info_list, function(x) x$root_coords),
      brown_rate_list = lapply(info_list, function(x) x$brown_rate_mat),
      mean_brown_rate = Reduce("+", lapply(info_list, function(x) x$brown_rate_mat))/length(info_list)
    )
  )
}



#########################################################
### PROCESSING AND DISTRIBUTION CALCULATION FUNCTIONS ###
#########################################################

process_tree <- function(tree, trait, brown_rate_mat, root_coords, missing_inds, 
                         retain = c('tree_covar', 'trait_reorder', 'root_coords', 'brown_rate_mat', 'cond_distr')) {
  
  element_list <- list()
  
  element_list[['tree_covar']] <- vcv.phylo(tree)
  
  element_list[['trait_reorder']] <- trait[rank(match(rownames(element_list[['tree_covar']]), rownames(trait)), na.last = NA),]
  
  element_list[['root_coords']] <- root_coords
  element_list[['brown_rate_mat']] <- brown_rate_mat
  
  element_list[['cond_distr']] <- calc_cond_distr(element_list[['tree_covar']],
                                                  brown_rate_mat = element_list[['brown_rate_mat']],
                                                  root_coord = element_list[['root_coords']],
                                                  trait_mat = element_list[['trait_reorder']],
                                                  missing_inds = missing_inds)
  
  return(element_list[names(element_list) %in% retain])
}


process_trees <- function(tree_list, trait, brown_rate_mat, root_coords_list, missing_inds, progress = TRUE) {
  
  if (isTRUE(progress)) prog_bar <- txtProgressBar(min = 1, max = length(tree_list), style = 3, char = "+")
  
  cond_info_list <- list()
  for (i in 1:length(tree_list)) {
    cond_info_list[[i]] <- process_tree(tree = tree_list[[i]], 
                                        trait = miss_info_mat, 
                                        brown_rate_mat = brown_rate_mat, 
                                        root_coords = root_coords_list[[i]],
                                        missing_inds = missing_inds,
                                        retain = 'cond_distr')
    
    if (isTRUE(progress)) setTxtProgressBar(prog_bar, i)
  }
  
  return(cond_info_list)
  
}


calc_kronmat <- function(tree_covar, rate_mat, missing_inds = NULL, include_inds = TRUE) {
  
  output_list <- list()
  
  output_list[['kron_mat']] <- kronecker(tree_covar, rate_mat, make.dimnames = TRUE)
  
  if (isTRUE(include_inds)) {
    output_list[['yes_ind']] <- which(!gsub(':[xy]', '', rownames(output_list[['kron_mat']])) %in% missing_inds)
    output_list[['no_ind']] <- which(gsub(':[xy]', '', rownames(output_list[['kron_mat']])) %in% missing_inds)
  }
  
  return(output_list)
}

calc_cond_mean <- function(covar_12, inv_covar22, 
                           trait_vec, means, 
                           no_ind, yes_ind) {
  mean1 <- rep(means, times = length(no_ind)/2)
  mean2 <- rep(means, times = length(yes_ind)/2)
  return(
    mean1 + (covar_12 %*% inv_covar22 %*% (trait_vec -  mean2))
  )
}

calc_cond_var_elements <- function(kron_varcovar_ratemat,
                                   no_ind,
                                   yes_ind) {

  covar_11 <- kron_varcovar_ratemat[no_ind,no_ind]
  covar_12 <- kron_varcovar_ratemat[no_ind,yes_ind]
  inv_covar_22 <- solve(kron_varcovar_ratemat[yes_ind,yes_ind])
  covar_21 <- kron_varcovar_ratemat[yes_ind,no_ind]

  #cond_mean <- mean1 + (covar_12 %*% solve(covar_22) %*% (trait_vec -  mean2))
  #cond_var <- covar_11 - covar_12%*%inv_covar_22%*%covar_21

  return(list(
    covar_11 = covar_11,
    covar_12 = covar_12,
    inv_covar_22 = inv_covar_22,
    covar_21 = covar_21,
    cond_var = covar_11 - covar_12%*%inv_covar_22%*%covar_21
  ))

  } 


calc_cond_distr <- function(tree_covar, brown_rate_mat, root_coords, trait_mat, missing_inds) {
  
  trait_reorder <- trait_mat[rank(match(rownames(tree_covar), rownames(trait_mat)), na.last = NA),]
  
  element_list <- list()
  
  kron <- calc_kronmat(tree_covar = tree_covar, 
                       rate_mat = brown_rate_mat,
                       missing_inds = missing_inds, 
                       include_inds = TRUE)
  
  element_list[['cond_var_elements']] <- calc_cond_var_elements(
    kron_varcovar_ratemat = kron$kron_mat,
    no_ind = kron$no_ind,
    yes_ind = kron$yes_ind
  )
  
  
  #calculate conditional mean
  element_list[['cond_mean']] <- calc_cond_mean(covar_12 = element_list[['cond_var_elements']]$covar_12,
                                                inv_covar22 = element_list[['cond_var_elements']]$inv_covar_22,
                                                trait_vec = c(t(trait_reorder)), #element_list[['data_prep_output']]$vector_trait,
                                                means = root_coords, #element_list[['ls_root']], #currently set at least squares estimate of root coordinates
                                                no_ind = kron$no_ind,
                                                yes_ind = kron$yes_ind)
  
  
  #element_list[['cond_var']] <- element_list[['cond_var_elements']]$cond_var
  
  #return(element_list[names(element_list) %in% retain])
  return(
    list(
      condmean = element_list[['cond_mean']],
      condvar = element_list[['cond_var_elements']]$cond_var
    )
  )
}



##################################
### SINGLE LOCATION ESTIMATION ###
##################################

multiply_mvnorms <- function(mu_list, Sigma_list) {
  
  d <- length(mu_list[[1]])
  
  precision_sum <- matrix(0, d, d)
  weighted_mean_sum <- rep(0, d)
  
  for(i in seq_along(mu_list)) {
    P <- solve(Sigma_list[[i]])
    
    precision_sum <- precision_sum + P
    weighted_mean_sum <- weighted_mean_sum + P %*% mu_list[[i]]
  }
  
  Sigma <- solve(precision_sum)
  mu <- Sigma %*% weighted_mean_sum
  
  list(mu = drop(mu),
       Sigma = Sigma)
}


calc_vcv_single_tree <- function(tree) {

  tip_order <- tree$tip.label
  vcv_mat <- vcv.phylo(tree)

  return(vcv_mat[match(tip_order, rownames(vcv_mat)),
                 match(tip_order, colnames(vcv_mat))]
  )
}


create_design <- function(mean_vec) {
  
  cbind(
    rep(c(1, 0), times = length(mean_vec)/2),
    rep(c(0, 1), times = length(mean_vec)/2)
  )
  
}


weighted_least_squares <- function(A, Sigma_list, mu) {
  sigma_inv <- as.matrix(Matrix::bdiag(lapply(Sigma_list, solve)))
  
  elementone <- solve(t(A)%*%sigma_inv%*%A)
  elementtwo <- t(A)%*%sigma_inv%*%mu
  
  return(elementone %*% elementtwo)
}



condition_mvn <- function(mu, V, A) {
  
  inv_AVAt <- solve(A %*% V %*% t(A))
  
  mu_c <- mu - V%*%t(A)%*%inv_AVAt%*%(A%*%mu)
  var_c <- V - V%*%t(A)%*%inv_AVAt%*%(A%*%V)
  
  list(
    mean = mu_c,
    covariance = var_c
  )
}


cond_log_lik_fixed_mean <- function(cond_list,
                                    prep_list,
                                    params,
                                    miss_n,
                                    weights = 1,
                                    means_list = NULL) {
  
  log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
  
  for (i in seq_along(log_lik_vec)) {
    log_lik_vec[i] <- mvtnorm::dmvnorm(rep(params, times = miss_n),
                                       means_list[[i]],
                                       sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
                                       log = TRUE)
  }
  return(sum(log_lik_vec*weights))
}



  
cond_log_lik_wrapper <- function(cond_list, 
                                 prep_list, 
                                 params, 
                                 miss_n, 
                                 means = NULL,
                                 weights = 1,
                                 type = c('fixed_mean', 'mean_estimate')) {
  
  return(
    switch(type,
           fixed_mean = {
             cond_log_lik_fixed_mean(cond_list = cond_list,
                                     prep_list = prep_list,
                                     params = params,
                                     miss_n = miss_n,
                                     weights = weights,
                                     means_list = means)
           },
           mean_estimate = {
             cond_log_lik_mean_est(
               cond_list = cond_list,
               prep_list = prep_list,
               params = params,
               weights = weights,
               miss_n = miss_n
             )
           })
  )
}



#####################
### OLD FUNCTIONS ###
#####################

# posterior_prob_locs <- function(pars,
#                                 cond_param_list,
#                                 prep_list,
#                                 miss_n,
#                                 means = NULL,
#                                 weights = 1,
#                                 type = c('fixed_mean', 'mean_estimate')) {
#   
#   locs_prior <- sum(dunif(pars, -10, 110, log = TRUE))
#   
#   #return(list(mu = mu, rate_mat = rate_mat, prior_mu = prior_mu, prior_rat_mat = prior_rat_mat))
#   output <- cond_log_lik_wrapper(cond_list = cond_param_list, 
#                                  prep_list = prep_list, 
#                                  params = pars, 
#                                  miss_n = miss_n, 
#                                  means = means,
#                                  weights = weights,
#                                  type = type) + locs_prior
#   
#   if (is.infinite(output)) output <- -50000
#   return(
#     output
#   )
# }


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
# 
# 
# data_prep_multitree <- function(vcv_list, trait, rate_mat, missing_inds) {
#   lapply(vcv_list, function(x, trait, rate_mat, miss_inds) {
#     data_prep(trait = trait,
#               phylo_covar = x,
#               rate_mat = rate_mat,
#               missing_inds = miss_inds)
#   }, trait = trait, rate_mat = rate_mat, miss_inds = missing_inds)
# } 

# expected_root_state_old <- function(vcv, trait_vec) {
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
# 
# 
# calc_root_coords_old <- function(missing_inds, kron_mat, trait_vector) {
#   vcv_subset <- kron_mat[!(gsub(":[xy]", "", rownames(kron_mat)) %in% missing_inds),
#                          !(gsub(":[xy]", "", colnames(kron_mat)) %in% missing_inds)]
#   
#   root_state <- expected_root_state(vcv_subset, 
#                                     trait_vector)
#   
#   return(t(root_state))
# }



#' process_single_tree_old <- function(tree,
#'                                     trait_mat,
#'                                     rate_mat,
#'                                     missing_inds,
#'                                     #prior_variance,
#'                                     retain = c('tree_covar', 
#'                                                'data_prep_output', 
#'                                                'ls_root', 
#'                                                #'marginal_variance', 
#'                                                'cond_var_elements',
#'                                                'cond_mean',
#'                                                'cond_var')) {
#'   
#'   if (identical(retain, 'identical')) retain <- c('tree_covar', 'data_prep_output', 'ls_root', 'marginal_variance', 'cond_var_elements','cond_mean','cond_var')
#'   
#'   element_list <- list()
#'   #calc_covar
#'   element_list[['tree_covar']] <- calc_vcv_single_tree(tree)
#'   
#'   
#'   #next processing steps
#'   element_list[['data_prep_output']] <- data_prep(trait = trait_mat,
#'                                                   phylo_covar = element_list[['tree_covar']],
#'                                                   rate_mat = rate_mat,
#'                                                   missing_inds = missing_inds)
#'   
#'   
#'   #calculate least squares root value for prior (or just create independent prior)
#'   
#'   element_list[['ls_root']] <- calc_root_coords(missing_inds = missing_inds, 
#'                                                 kron_mat = element_list[['data_prep_output']]$kron_varcovar_ratemat, 
#'                                                 trait_vector = element_list[['data_prep_output']]$vector_trait)
#'   
#'   
#'   #element_list[['marginal_variance']] <- add_variance(kron_mat = element_list[['data_prep_output']]$kron_varcovar_ratemat, 
#'   #                                                    var = prior_variance)
#'   
#'   element_list[['cond_var_elements']] <- calc_cond_var_elements(
#'     kron_varcovar_ratemat = element_list[['data_prep_output']]$kron_varcovar_ratemat,
#'     no_ind = element_list[['data_prep_output']]$no_ind,
#'     yes_ind = element_list[['data_prep_output']]$yes_ind
#'   )
#'   
#'   
#'   #calculate conditional mean
#'   element_list[['cond_mean']] <- calc_cond_mean(covar_12 = element_list[['cond_var_elements']]$covar_12,
#'                                                 inv_covar22 = element_list[['cond_var_elements']]$inv_covar_22,
#'                                                 trait_vec = element_list[['data_prep_output']]$vector_trait,
#'                                                 means = element_list[['ls_root']], #currently set at least squares estimate of root coordinates
#'                                                 no_ind = element_list[['data_prep_output']]$no_ind,
#'                                                 yes_ind = element_list[['data_prep_output']]$yes_ind)
#'   
#'   
#'   element_list[['cond_var']] <- element_list[['cond_var_elements']]$cond_var
#'   
#'   return(element_list[names(element_list) %in% retain])
#'   
#' }
#' 
#' 
#' 
#' process_single_tree <- function(tree,
#'                                 trait_mat,
#'                                 rate_mat,
#'                                 missing_inds,
#'                                 #prior_variance,
#'                                 retain = c('tree_covar', 
#'                                            'data_prep_output', 
#'                                            'ls_root', 
#'                                            #'marginal_variance', 
#'                                            'cond_var_elements',
#'                                            'cond_mean',
#'                                            'cond_var')) {
#'   
#'   if (identical(retain, 'identical')) retain <- c('tree_covar', 'data_prep_output', 'ls_root', 'marginal_variance', 'cond_var_elements','cond_mean','cond_var')
#'   
#'   element_list <- list()
#'   #calc_covar
#'   element_list[['tree_covar']] <- calc_vcv_single_tree(tree)
#'   
#'   
#'   #next processing steps
#'   element_list[['data_prep_output']] <- data_prep(trait = trait_mat,
#'                                                   phylo_covar = element_list[['tree_covar']],
#'                                                   rate_mat = rate_mat,
#'                                                   missing_inds = missing_inds)
#'   
#'   
#'   #calculate least squares root value for prior (or just create independent prior)
#'   
#'   element_list[['ls_root']] <- calc_root_coords(missing_inds = missing_inds, 
#'                                                 kron_mat = element_list[['data_prep_output']]$kron_varcovar_ratemat, 
#'                                                 trait_vector = element_list[['data_prep_output']]$vector_trait)
#'   
#'   
#'   #element_list[['marginal_variance']] <- add_variance(kron_mat = element_list[['data_prep_output']]$kron_varcovar_ratemat, 
#'   #                                                    var = prior_variance)
#'   
#'   element_list[['cond_var_elements']] <- calc_cond_var_elements(
#'     kron_varcovar_ratemat = element_list[['data_prep_output']]$kron_varcovar_ratemat,
#'     no_ind = element_list[['data_prep_output']]$no_ind,
#'     yes_ind = element_list[['data_prep_output']]$yes_ind
#'   )
#'   
#'   
#'   #calculate conditional mean
#'   element_list[['cond_mean']] <- calc_cond_mean(covar_12 = element_list[['cond_var_elements']]$covar_12,
#'                                                 inv_covar22 = element_list[['cond_var_elements']]$inv_covar_22,
#'                                                 trait_vec = element_list[['data_prep_output']]$vector_trait,
#'                                                 means = element_list[['ls_root']], #currently set at least squares estimate of root coordinates
#'                                                 no_ind = element_list[['data_prep_output']]$no_ind,
#'                                                 yes_ind = element_list[['data_prep_output']]$yes_ind)
#'   
#'   
#'   element_list[['cond_var']] <- element_list[['cond_var_elements']]$cond_var
#'   
#'   return(element_list[names(element_list) %in% retain])
#'   
#' }
#' 
#' process_multitree <- function(tree_list, 
#'                               trait_mat,
#'                               rate_mat,
#'                               missing_inds,
#'                               #prior_variance,
#'                               retain = c('tree_covar',
#'                                          'data_prep_output',
#'                                          'ls_root',
#'                                          #'marginal_variance', 
#'                                          'cond_var_elements',
#'                                          'cond_mean',
#'                                          'cond_var'),
#'                               extract_distribution_info = TRUE) {
#'   
#'   #if (isTRUE(extract_distribution_info) & !all(c('cond_mean', 'cond_var') %in% retain))
#'   #  stop('If you want to extract the conditional means and covariances, you must include those in the retain argument.')
#'   
#'   if (isTRUE(extract_distribution_info)) retain <- c('cond_mean','cond_var')
#'   tree_process_list <- list()
#'   
#'   prog_bar <- txtProgressBar(min = 1, max = length(tree_list), style = 3)
#'   
#'   for (TREE in seq_len(length(tree_list))) {
#'     tree_process_list[[TREE]] <- process_single_tree_old(tree_list[[TREE]],
#'                                                          trait_mat,
#'                                                          rate_mat,
#'                                                          missing_inds,
#'                                                          #prior_variance,
#'                                                          retain = retain)
#'     setTxtProgressBar(prog_bar, TREE)
#'   }
#'   
#'   
#'   if (isTRUE(extract_distribution_info)) {
#'     output_list <- list()
#'     output_list[['cond_mean']] <- lapply(tree_process_list, function(x) drop(x$cond_mean))
#'     output_list[['cond_var']] <- lapply(tree_process_list, function(x) x$cond_var)
#'     return(output_list)
#'   }
#'   
#'   return(tree_process_list)
#' }

# custom_metrop <- function(iters = 100, #number of MCMC iterations after the initial step
#                           init_params, #initial parameter values
#                           cond_param_list,
#                           prep_list,
#                           tree_means = NULL,
#                           weights = 1,
#                           proposal_scale = c(0.1, 0.1), #first proposal is for locations, 2nd is for root states
#                           miss_n,
#                           type = c('fixed_mean', 'mean_estimate'),
#                           report_progress = TRUE) {
#   
#   #=== PREP STEPS ===
#   param_count = length(init_params) #number of estimated parameters
#   
#   #covariance matrix of proposal distribution
#   #(off-diagonals are 0 so this is equivalent to a bunch of 1D normals)
#   sigma_prop_mat <- diag(param_count)*c(rep(proposal_scale[1], 2), rep(proposal_scale[2], param_count - 2))
#   
#   #matrix to hold posterior values
#   posterior_mat <- matrix(data = NA, nrow = iters + 1, 
#                           ncol = param_count)
#   #==================
#   
#   
#   #=== INITIAL MCMC STEP ===
#   if (isTRUE(report_progress)) prog_bar <- txtProgressBar(min = 1, max = iters, style = 3, char = "+")
#   
#   posterior_mat[1,] <- init_params
#   
#   log_prop <- vector(mode = 'numeric', length = iters + 1)
#   
#   # log_prop[1] <- posterior_prob_locs_center_prop_cond_nfm(pars = init_params, 
#   #                                                         prep_list = prep_list,
#   #                                                         cond_param_list = cond_param_list,
#   #                                                         miss_n = miss_n)
#   log_prop[1] <- posterior_prob_locs(pars = init_params,
#                                      cond_param_list = cond_param_list,
#                                      prep_list = prep_list,
#                                      miss_n = miss_n,
#                                      means = tree_means,
#                                      weights = weights,
#                                      type = type)
#   
#   if (isTRUE(report_progress)) setTxtProgressBar(prog_bar, 1)
#   #=========================
#   
#   
#   #=== ALL SUBSEQUENT MCMC STEPS ===
#   for (i in 2:(iters + 1)) {
#     #step 1: proposal
#     proposal_vals <- posterior_mat[i-1,] + MASS::mvrnorm(1, rep(0, param_count), sigma_prop_mat) #is the right sigma mat
#     
#     #step 2: accept/reject
#     # proposal_prop <- posterior_prob_locs_center_prop_cond_nfm(pars = proposal_vals, 
#     #                                                           prep_list = prep_list,
#     #                                                           cond_param_list = cond_param_list,
#     #                                                           miss_n = miss_n)
#     
#     proposal_prop <- posterior_prob_locs(pars = proposal_vals,
#                                          cond_param_list = cond_param_list,
#                                          prep_list = prep_list,
#                                          miss_n = miss_n,
#                                          means = tree_means,
#                                          weights = weights,
#                                          type = type)
#     
#     ### Currently using a symmetric proposal distribution so I'm not explicitly calculating
#     ### the Hastings ratio
#     #q_old_new <- dnorm(posterior_mat[i-1,], proposal_vals, proposal_sd, log = TRUE)
#     #q_new_old <- dnorm(proposal_vals, posterior_mat[i-1,], proposal_sd, log = TRUE)
#     
#     ##log_hastings_rat <- sum(dnorm(posterior_mat[i-1,], proposal_vals, c(0.1, 0.1, 0.08, 0.08, 0.08), log = TRUE) - dnorm(proposal_vals, posterior_mat[1,], c(0.1, 0.1, 0.08, 0.08, 0.08), log = TRUE))
#     #acceptance_ratio <- min((proposal_prop - log_prop[i-1]) + sum(q_old_new - q_new_old), 0)
#     
#     acceptance_ratio <- min(proposal_prop - log_prop[i-1], 0)
#     
#     #if accepted --> add proposed values to chain.
#     #if not --> stay at current values
#     if (log(runif(1)) <= acceptance_ratio) {
#       posterior_mat[i,] <- proposal_vals
#       log_prop[i] <- proposal_prop
#     } else {
#       posterior_mat[i,] <- posterior_mat[i-1,]
#       log_prop[i] <- log_prop[i - 1]
#     }
#     if (isTRUE(report_progress)) setTxtProgressBar(prog_bar, i)
#   }
#   
#   return(
#     list(logprop_vec = log_prop, #log posterior values
#          post_mat = posterior_mat) #parameter vlaues
#   )
# }
# 
# 
# 
# cond_log_lik_mean_est <- function(cond_list,
#                                   prep_list,
#                                   params,
#                                   weights = 1,
#                                   miss_n) {
#   
#   means <- params[3:length(params)]
#   locs <- params[c(1, 2)]
#   
#   log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
#   
#   for (i in seq_along(log_lik_vec)) {
#     #message(i)
#     cond_mean <- calc_cond_mean(covar_12 = cond_list[[i]]$covar_12,
#                                 inv_covar22 = cond_list[[i]]$inv_covar_22,
#                                 trait_vec = prep_list[[i]]$vector_trait,
#                                 means = means[c(i*2 - 1, i*2)],
#                                 no_ind = prep_list[[i]]$no_ind,
#                                 yes_ind = prep_list[[i]]$yes_ind)
#     
#     
#     log_lik_vec[i] <- mvtnorm::dmvnorm(rep(locs, times = miss_n),
#                                        cond_mean,
#                                        sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
#                                        log = TRUE)
#   }
#   return(sum(log_lik_vec*weights))
#}

# data_prep <- function(trait,
#                       phylo_covar,
#                       rate_mat,
#                       missing_inds) {
#   
#   dimnames(rate_mat) <- list(c('x', 'y'), c('x', 'y'))
#   
#   #trait_reorder <- trait[match(rownames(trait), rownames(phylo_covar)),]
#   trait_reorder <- trait[rank(match(rownames(phylo_covar), rownames(trait)), na.last = NA),]
#   
#   ### centering phylo covar and traits ###
#   # n_trait <- nrow(trait_reorder)
#   # n_phylo <- nrow(phylo_covar)
#   # mat1_trait <- matrix(-1/n_trait, nrow = n_trait - 1, ncol = n_trait)
#   # diag(mat1_trait) <- (n_trait - 1)/n_trait
#   # mat1_phylo <- matrix(-1/n_phylo, nrow = n_phylo - 1, ncol = n_phylo)
#   # diag(mat1_phylo) <- (n_phylo - 1)/n_phylo
#   
#   #centered_mat <- mat1_trait %*% trait_reorder
#   #rownames(centered_mat) <- rownames(trait_reorder)[-nrow(trait_reorder)]
#   
#   #center_phylo_covar <- mat1_phylo%*%phylo_covar%*%t(mat1_phylo)
#   #dimnames(center_phylo_covar) <- lapply(dimnames(phylo_covar), function(x) x[-length(x)])
#   
#   kron_varcovar_ratemat <- kronecker(phylo_covar, rate_mat, make.dimnames = TRUE)
#   #inverse_kron_varcovar_ratemat <- solve(kron_varcovar_ratemat, tol = 1e-18)
#   
#   #centered_trait_NA <- matrix(NA, n_phylo - 1, 2, 
#   #                            dimnames = list(rownames(center_phylo_covar), c('t1', 't2')))
#   #centered_trait_NA[match(rownames(centered_mat), rownames(center_phylo_covar)),] <- centered_mat
#   
#   yes_ind <- which(!gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
#   no_ind <- which(gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
#   
#   
#   return(
#     list(kron_varcovar_ratemat = kron_varcovar_ratemat,
#          trait_mat = trait_reorder,
#          vector_trait = c(t(trait_reorder)),
#          yes_ind = yes_ind,
#          no_ind = no_ind)
#   )
# }
# 
# 
#   
# }
# 
# calc_cond_var_elements_multitree <- function(prep_list) {
#   lapply(prep_list, function(x) {
#     calc_cond_var_elements(x$kron_varcovar_ratemat, x$no_ind, x$yes_ind)
#   })
# }
# 
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


#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# posterior_prob_locs_center_prop_cond_nfm_fixedmean <- function(pars,
#                                                                means,
#                                                                cond_param_list,
#                                                                prep_list,
#                                                                miss_n) {
#   
#   locs_prior <- sum(dunif(pars, -10, 110, log = TRUE))
#   
#   #return(list(mu = mu, rate_mat = rate_mat, prior_mu = prior_mu, prior_rat_mat = prior_rat_mat))
#   output <- cond_log_lik_nfm_fixedmean(cond_list = cond_param_list, 
#                                        prep_list = prep_list,
#                                        means = means,
#                                        params = pars, 
#                                        miss_n = miss_n) + locs_prior
#   
#   if (is.infinite(output)) output <- -50000
#   return(
#     output
#   )
# }


# cond_log_lik_nfm <- function(cond_list, prep_list, params, miss_n) {
#   
#   means <- params[3:length(params)]
#   locs <- params[c(1, 2)]
#   
#   log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
#   
#   for (i in seq_along(log_lik_vec)) {
#     #message(i)
#     cond_mean <- calc_cond_mean(covar_12 = cond_list[[i]]$covar_12, 
#                                 inv_covar22 = cond_list[[i]]$inv_covar_22, 
#                                 trait_vec = prep_list[[i]]$vector_trait, 
#                                 means = means[c(i*2 - 1, i*2)], 
#                                 no_ind = prep_list[[i]]$no_ind, 
#                                 yes_ind = prep_list[[i]]$yes_ind)
#     
#     
#     log_lik_vec[i] <- mvtnorm::dmvnorm(rep(locs, times = miss_n), 
#                                        cond_mean, 
#                                        sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
#                                        log = TRUE)
#   }
#   return(sum(log_lik_vec))
# }
# 
# 
# 
# cond_log_lik_fm <- function(cond_list, prep_list, means, params, miss_n) {
#   
#   log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
#   
#   for (i in seq_along(log_lik_vec)) {
#     cond_mean <- calc_cond_mean(covar_12 = cond_list[[i]]$covar_12, 
#                                 inv_covar22 = cond_list[[i]]$inv_covar_22, 
#                                 trait_vec = prep_list[[i]]$vector_trait, 
#                                 means = means[c(i*2 - 1, i*2)], 
#                                 no_ind = prep_list[[i]]$no_ind, 
#                                 yes_ind = prep_list[[i]]$yes_ind)
#     
#     
#     log_lik_vec[i] <- mvtnorm::dmvnorm(rep(params, times = miss_n), 
#                                        cond_mean, 
#                                        sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
#                                        log = TRUE)
#   }
#   return(sum(log_lik_vec))
# }

# cond_log_lik <- function(cond_list,
#                          prep_list,
#                          params,
#                          miss_n,
#                          means = NULL,
#                          type = c('fixed_mean', 'mean_estimate')) {
#   if (type == 'fixed_mean') {
#     locs <- params
#     if (is.null(means)) stop('means cannot be null when using the fixed mean option.')
#   } else if (type == 'mean_estimate') {
#     means <- params[3:length(params)]
#     locs <- params[c(1, 2)]
#   }
# 
#   log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
# 
#   for (i in seq_along(log_lik_vec)) {
#     #message(i)
#     cond_mean <- calc_cond_mean(covar_12 = cond_list[[i]]$covar_12,
#                                 inv_covar22 = cond_list[[i]]$inv_covar_22,
#                                 trait_vec = prep_list[[i]]$vector_trait,
#                                 means = means[c(i*2 - 1, i*2)],
#                                 no_ind = prep_list[[i]]$no_ind,
#                                 yes_ind = prep_list[[i]]$yes_ind)
# 
# 
#     log_lik_vec[i] <- mvtnorm::dmvnorm(rep(locs, times = miss_n),
#                                        cond_mean,
#                                        sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
#                                        log = TRUE)
#   }
#   return(sum(log_lik_vec))
# }


# add_variance <- function(kron_mat, var) {
#   tip_count <- nrow(kron_mat)/2
#   one_mat <- matrix(1, nrow = tip_count, ncol = tip_count) #equivalent to rep(1, tip_count)%*%t(rep(1, tip_count))
#   diagonal_var_mat <- diag(2)*var
#   var_mat <- kronecker(one_mat, diagonal_var_mat)
#   
#   return(kron_mat + var_mat)
# }

# calc_cond_var_elements_multitree_add_variance <- function(prep_list, variance = 0) {
#   lapply(prep_list, function(x, var) {
#     tip_count <- nrow(x$kron_varcovar_ratemat)/2
#     var_mat <- kronecker(rep(1, tip_count)%*%t(rep(1, tip_count)), diag(2)*var)
#     #return(list(kron = x$kron_varcovar_ratemat,
#     #            var_mat =  var_mat))
#     calc_cond_var_elements(x$kron_varcovar_ratemat + var_mat, 
#                            x$no_ind, x$yes_ind)
#   }, var = variance)
# }


#' process_single_tree <- function(tree,
#'                                 trait_mat,
#'                                 rate_mat = NULL,
#'                                 estimate_rate_mat = TRUE,
#'                                 missing_inds,
#'                                 #prior_variance,
#'                                 retain = c('tree_covar',
#'                                            'data_prep_output',
#'                                            'ls_root',
#'                                            #'marginal_variance',
#'                                            'cond_var_elements',
#'                                            'cond_mean',
#'                                            'cond_var')) {
#' 
#'   recover()
#'   if (identical(retain, 'identical')) retain <- c('tree_covar', 'data_prep_output', 'ls_root', 'marginal_variance', 'cond_var_elements','cond_mean','cond_var')
#' 
#'   element_list <- list()
#'   #calc_covar
#'   element_list[['tree_covar']] <- calc_vcv_single_tree(tree)
#' 
#' 
#'   #next processing steps
#'   #element_list[['data_prep_output']] <- data_prep(trait = trait_mat,
#'   #                                                phylo_covar = element_list[['tree_covar']],
#'   #                                                rate_mat = rate_mat,
#'   #                                                missing_inds = missing_inds)
#'   element_list[['data_prep_output']] <- data_reorder(trait_mat,
#'                                                      element_list[['tree_covar']])
#' 
#'   element_list[['ls_root']] <- calc_root_coords(missing_inds = missing_inds,
#'                                                 tree_covar = element_list[['tree_covar']],
#'                                                 trait_mat = element_list[['data_prep_output']]$trait_mat
#'   )
#' 
#'   if (isTRUE(estimate_rate_mat)) {
#'     element_list[['rate_mat']] <- est_brownian_rate_mat(tree_covar = element_list[['tree_covar']],
#'                                                         root_vals = element_list[['ls_root']],
#'                                                         trait_mat = element_list[['data_prep_output']]$trait_mat)
#'   } else {
#'     element_list[['rate_mat']] <- rate_mat
#'   }
#' 
#' 
#'   element_list[['kron']] <- calc_kronmat(tree_covar = element_list[['tree_covar']],
#'                                          rate_mat = element_list[['rate_mat']],
#'                                          missing_inds = missing_inds,
#'                                          include_inds = TRUE)
#' 
#'   element_list[['cond_var_elements']] <- calc_cond_var_elements(
#'     kron_varcovar_ratemat = element_list[['kron']],
#'     no_ind = element_list[['data_prep_output']]$no_ind,
#'     yes_ind = element_list[['data_prep_output']]$yes_ind
#'   )
#' 
#' 
#'   #calculate conditional mean
#'   element_list[['cond_mean']] <- calc_cond_mean(covar_12 = element_list[['cond_var_elements']]$covar_12,
#'                                                 inv_covar22 = element_list[['cond_var_elements']]$inv_covar_22,
#'                                                 trait_vec = element_list[['data_prep_output']]$vector_trait,
#'                                                 means = element_list[['ls_root']], #currently set at least squares estimate of root coordinates
#'                                                 no_ind = element_list[['data_prep_output']]$no_ind,
#'                                                 yes_ind = element_list[['data_prep_output']]$yes_ind)
#' 
#' 
#'   element_list[['cond_var']] <- element_list[['cond_var_elements']]$cond_var
#' 
#'   return(element_list[names(element_list) %in% retain])
#' 
#' }
#' 
#' 
#' 
#' process_single_tree <- function(tree,
#'                                 trait_mat,
#'                                 rate_mat = NULL,
#'                                 estimate_rate_mat = TRUE,
#'                                 missing_inds,
#'                                 #prior_variance,
#'                                 retain = c('tree_covar',
#'                                            'data_prep_output',
#'                                            'ls_root',
#'                                            #'marginal_variance',
#'                                            'cond_var_elements',
#'                                            'cond_mean',
#'                                            'cond_var')) {
#' 
#'   recover()
#'   if (identical(retain, 'identical')) retain <- c('tree_covar', 'data_prep_output', 'ls_root', 'marginal_variance', 'cond_var_elements','cond_mean','cond_var')
#' 
#'   element_list <- list()
#'   #calc_covar
#'   element_list[['tree_covar']] <- calc_vcv_single_tree(tree)
#' 
#' 
#'   #next processing steps
#'   #element_list[['data_prep_output']] <- data_prep(trait = trait_mat,
#'   #                                                phylo_covar = element_list[['tree_covar']],
#'   #                                                rate_mat = rate_mat,
#'   #                                                missing_inds = missing_inds)
#'   element_list[['data_prep_output']] <- data_reorder(trait,
#'                                                      tree_covar,
#'                                                      rate_mat,
#'                                                      missing_inds)
#' 
#'   element_list[['ls_root']] <- calc_root_coords(missing_inds = missing_inds,
#'                                                 tree_covar = element_list[['tree_covar']],
#'                                                 trait_mat = element_list[['data_prep_output']]$trait_mat
#'   )
#' 
#'   if (isTRUE(estimate_rate_mat)) {
#'     element_list[['rate_mat']] <- est_brownian_rate_mat(tree_covar = element_list[['tree_covar']],
#'                                                         root_vals = element_list[['ls_root']],
#'                                                         trait_mat = element_list[['data_prep_output']]$trait_mat)
#'   } else {
#'     element_list[['rate_mat']] <- rate_mat
#'   }
#' 
#' 
#'   element_list[['kron']] <- calc_kronmat(tree_covar = element_list[['tree_covar']],
#'                                          rate_mat = element_list[['rate_mat']],
#'                                          missing_inds = missing_inds,
#'                                          include_inds = TRUE)
#' 
#'   element_list[['cond_var_elements']] <- calc_cond_var_elements(
#'     kron_varcovar_ratemat = element_list[['kron']],
#'     no_ind = element_list[['data_prep_output']]$no_ind,
#'     yes_ind = element_list[['data_prep_output']]$yes_ind
#'   )
#' 
#' 
#'   #calculate conditional mean
#'   element_list[['cond_mean']] <- calc_cond_mean(covar_12 = element_list[['cond_var_elements']]$covar_12,
#'                                                 inv_covar22 = element_list[['cond_var_elements']]$inv_covar_22,
#'                                                 trait_vec = element_list[['data_prep_output']]$vector_trait,
#'                                                 means = element_list[['ls_root']], #currently set at least squares estimate of root coordinates
#'                                                 no_ind = element_list[['data_prep_output']]$no_ind,
#'                                                 yes_ind = element_list[['data_prep_output']]$yes_ind)
#' 
#' 
#'   element_list[['cond_var']] <- element_list[['cond_var_elements']]$cond_var
#' 
#'   return(element_list[names(element_list) %in% retain])
#' 
#' }

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
