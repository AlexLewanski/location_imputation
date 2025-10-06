#########################################################
#########################################################
### FUNCTIONS FOR SPATIAL IMPUTATION WITH GENEALOGIES ###
#########################################################
#########################################################
message('REQUIRED PACKAGES: MASS, Matrix, mvtnorm')


custom_metrop <- function(iters = 100, #number of MCMC iterations after the initial step
                          init_params, #initial parameter values
                          cond_param_list,
                          prep_list,
                          tree_means = NULL,
                          weights = 1,
                          proposal_scale = c(0.1, 0.1), #first proposal is for locations, 2nd is for root states
                          miss_n,
                          type = c('fixed_mean', 'mean_estimate'),
                          report_progress = TRUE) {
  
  #=== PREP STEPS ===
  param_count = length(init_params) #number of estimated parameters
  
  #covariance matrix of proposal distribution
  #(off-diagonals are 0 so this is equivalent to a bunch of 1D normals)
  sigma_prop_mat <- diag(param_count)*c(rep(proposal_scale[1], 2), rep(proposal_scale[2], param_count - 2))
  
  #matrix to hold posterior values
  posterior_mat <- matrix(data = NA, nrow = iters + 1, 
                          ncol = param_count)
  #==================
  
  
  #=== INITIAL MCMC STEP ===
  if (isTRUE(report_progress)) prog_bar <- txtProgressBar(min = 1, max = iters, style = 3, char = "+")
  
  posterior_mat[1,] <- init_params
  
  log_prop <- vector(mode = 'numeric', length = iters + 1)
  
  # log_prop[1] <- posterior_prob_locs_center_prop_cond_nfm(pars = init_params, 
  #                                                         prep_list = prep_list,
  #                                                         cond_param_list = cond_param_list,
  #                                                         miss_n = miss_n)
  log_prop[1] <- posterior_prob_locs(pars = init_params,
                                     cond_param_list = cond_param_list,
                                     prep_list = prep_list,
                                     miss_n = miss_n,
                                     means = tree_means,
                                     weights = weights,
                                     type = type)
  
  if (isTRUE(report_progress)) setTxtProgressBar(prog_bar, 1)
  #=========================
  
  
  #=== ALL SUBSEQUENT MCMC STEPS ===
  for (i in 2:(iters + 1)) {
    #step 1: proposal
    proposal_vals <- posterior_mat[i-1,] + MASS::mvrnorm(1, rep(0, param_count), sigma_prop_mat) #is the right sigma mat
    
    #step 2: accept/reject
    # proposal_prop <- posterior_prob_locs_center_prop_cond_nfm(pars = proposal_vals, 
    #                                                           prep_list = prep_list,
    #                                                           cond_param_list = cond_param_list,
    #                                                           miss_n = miss_n)
    
    proposal_prop <- posterior_prob_locs(pars = proposal_vals,
                                         cond_param_list = cond_param_list,
                                         prep_list = prep_list,
                                         miss_n = miss_n,
                                         means = tree_means,
                                         weights = weights,
                                         type = type)
    
    ### Currently using a symmetric proposal distribution so I'm not explicitly calculating
    ### the Hastings ratio
    #q_old_new <- dnorm(posterior_mat[i-1,], proposal_vals, proposal_sd, log = TRUE)
    #q_new_old <- dnorm(proposal_vals, posterior_mat[i-1,], proposal_sd, log = TRUE)
    
    ##log_hastings_rat <- sum(dnorm(posterior_mat[i-1,], proposal_vals, c(0.1, 0.1, 0.08, 0.08, 0.08), log = TRUE) - dnorm(proposal_vals, posterior_mat[1,], c(0.1, 0.1, 0.08, 0.08, 0.08), log = TRUE))
    #acceptance_ratio <- min((proposal_prop - log_prop[i-1]) + sum(q_old_new - q_new_old), 0)
    
    acceptance_ratio <- min(proposal_prop - log_prop[i-1], 0)
    
    #if accepted --> add proposed values to chain.
    #if not --> stay at current values
    if (log(runif(1)) <= acceptance_ratio) {
      posterior_mat[i,] <- proposal_vals
      log_prop[i] <- proposal_prop
    } else {
      posterior_mat[i,] <- posterior_mat[i-1,]
      log_prop[i] <- log_prop[i - 1]
    }
    if (isTRUE(report_progress)) setTxtProgressBar(prog_bar, i)
  }
  
  return(
    list(logprop_vec = log_prop, #log posterior values
         post_mat = posterior_mat) #parameter vlaues
  )
}



cond_log_lik_mean_est <- function(cond_list,
                                  prep_list,
                                  params,
                                  weights = 1,
                                  miss_n) {
  
  means <- params[3:length(params)]
  locs <- params[c(1, 2)]
  
  log_lik_vec <- vector(mode = 'numeric', length = length(cond_list))
  
  for (i in seq_along(log_lik_vec)) {
    #message(i)
    cond_mean <- calc_cond_mean(covar_12 = cond_list[[i]]$covar_12,
                                inv_covar22 = cond_list[[i]]$inv_covar_22,
                                trait_vec = prep_list[[i]]$vector_trait,
                                means = means[c(i*2 - 1, i*2)],
                                no_ind = prep_list[[i]]$no_ind,
                                yes_ind = prep_list[[i]]$yes_ind)
    
    
    log_lik_vec[i] <- mvtnorm::dmvnorm(rep(locs, times = miss_n),
                                       cond_mean,
                                       sigma = as.matrix(Matrix::forceSymmetric(cond_list[[i]]$cond_var)),
                                       log = TRUE)
  }
  return(sum(log_lik_vec*weights))
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






posterior_prob_locs <- function(pars,
                                cond_param_list,
                                prep_list,
                                miss_n,
                                means = NULL,
                                weights = 1,
                                type = c('fixed_mean', 'mean_estimate')) {
  
  locs_prior <- sum(dunif(pars, -10, 110, log = TRUE))
  
  #return(list(mu = mu, rate_mat = rate_mat, prior_mu = prior_mu, prior_rat_mat = prior_rat_mat))
  output <- cond_log_lik_wrapper(cond_list = cond_param_list, 
                                 prep_list = prep_list, 
                                 params = pars, 
                                 miss_n = miss_n, 
                                 means = means,
                                 weights = weights,
                                 type = type) + locs_prior
  
  if (is.infinite(output)) output <- -50000
  return(
    output
  )
}


calc_vcv_multitree <- function(tree_list) {
  tip_order <- tree_list[[1]]$tip.label
  return(
    lapply(tree_list, function(x, tip_order) {
      vcv_mat <- vcv.phylo(x)
      return(vcv_mat[match(tip_order, rownames(vcv_mat)),
                     match(tip_order, colnames(vcv_mat))]
      )
    }, tip_order = tip_order)
  )
}


data_prep_multitree <- function(vcv_list, trait, rate_mat, missing_inds) {
  lapply(vcv_list, function(x, trait, rate_mat, miss_inds) {
    data_prep(trait = trait,
              phylo_covar = x,
              rate_mat = rate_mat,
              missing_inds = miss_inds)
  }, trait = trait, rate_mat = rate_mat, miss_inds = missing_inds)
} 


data_prep <- function(trait,
                      phylo_covar,
                      rate_mat,
                      missing_inds) {
  
  dimnames(rate_mat) <- list(c('x', 'y'), c('x', 'y'))
  
  #trait_reorder <- trait[match(rownames(trait), rownames(phylo_covar)),]
  trait_reorder <- trait[rank(match(rownames(phylo_covar), rownames(trait)), na.last = NA),]
  
  ### centering phylo covar and traits ###
  # n_trait <- nrow(trait_reorder)
  # n_phylo <- nrow(phylo_covar)
  # mat1_trait <- matrix(-1/n_trait, nrow = n_trait - 1, ncol = n_trait)
  # diag(mat1_trait) <- (n_trait - 1)/n_trait
  # mat1_phylo <- matrix(-1/n_phylo, nrow = n_phylo - 1, ncol = n_phylo)
  # diag(mat1_phylo) <- (n_phylo - 1)/n_phylo
  
  #centered_mat <- mat1_trait %*% trait_reorder
  #rownames(centered_mat) <- rownames(trait_reorder)[-nrow(trait_reorder)]
  
  #center_phylo_covar <- mat1_phylo%*%phylo_covar%*%t(mat1_phylo)
  #dimnames(center_phylo_covar) <- lapply(dimnames(phylo_covar), function(x) x[-length(x)])
  
  kron_varcovar_ratemat <- kronecker(phylo_covar, rate_mat, make.dimnames = TRUE)
  #inverse_kron_varcovar_ratemat <- solve(kron_varcovar_ratemat, tol = 1e-18)
  
  #centered_trait_NA <- matrix(NA, n_phylo - 1, 2, 
  #                            dimnames = list(rownames(center_phylo_covar), c('t1', 't2')))
  #centered_trait_NA[match(rownames(centered_mat), rownames(center_phylo_covar)),] <- centered_mat
  
  yes_ind <- which(!gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
  no_ind <- which(gsub(':[xy]', '', rownames(kron_varcovar_ratemat)) %in% missing_inds)
  
  
  return(
    list(kron_varcovar_ratemat = kron_varcovar_ratemat,
         trait_mat = trait_reorder,
         vector_trait = c(t(trait_reorder)),
         yes_ind = yes_ind,
         no_ind = no_ind)
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

calc_cond_var_elements_multitree <- function(prep_list) {
  lapply(prep_list, function(x) {
    calc_cond_var_elements(x$kron_varcovar_ratemat, x$no_ind, x$yes_ind)
  })
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

calc_cond_mean_list <- function(cond_param_list,
                                prep_list,
                                tree_means) {
  cond_mean_list <- list()
  for (i in 1:length(cond_param_list)) {
    cond_mean_list[[i]] <- calc_cond_mean(covar_12 = cond_param_list[[i]]$covar_12, 
                                          inv_covar22 = cond_param_list[[i]]$inv_covar_22, 
                                          trait_vec = prep_list[[i]]$vector_trait, 
                                          means = tree_means[c(i*2 - 1, i*2)], 
                                          no_ind = prep_list[[i]]$no_ind, 
                                          yes_ind = prep_list[[i]]$yes_ind)
  }
  
  return(cond_mean_list)
}



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
