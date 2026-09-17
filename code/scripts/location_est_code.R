#########################################################
#########################################################
### FUNCTIONS FOR SPATIAL IMPUTATION WITH GENEALOGIES ###
#########################################################
#########################################################
message('REQUIRED PACKAGES: MASS, Matrix, mvtnorm')

################
### OVERVIEW ###
################
#This script contains functions for estimating locations from genome-wide genealogies.
#Other custom functions that are used for the project but are not directly involved
#in the location inference are housed in custom_project_funcs.R.


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

# weighted_least_squares_softassign <- function(A, weight_mat, Sigma_list, mu) {
#   #recover()
#   sigma_inv <- as.matrix(Matrix::bdiag(lapply(Sigma_list, solve)))
#   sigma_weight_mat <- sigma_inv*weight_mat #this should be elementwise multiplication I think
#   elementone <- solve(t(A)%*%sigma_weight_mat%*%A)
#   elementtwo <- t(A)%*%sigma_weight_mat%*%mu
#   
#   return(elementone %*% elementtwo)
# }



condition_mvn <- function(mu, V, A) {
  
  inv_AVAt <- solve(A %*% V %*% t(A))
  
  mu_c <- mu - V%*%t(A)%*%inv_AVAt%*%(A%*%mu)
  var_c <- V - V%*%t(A)%*%inv_AVAt%*%(A%*%V)
  
  list(
    mean = mu_c,
    covariance = var_c
  )
}


#can probably delete. Not currently in use. 
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



#can probably delete. Not currently in use.  
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



##############################################
### FUNCTIONS FOR EXPECTATION-MAXIMIZATION ###
##############################################

em_hardassign_wrapper <- function(times = 1,
                                  cond_mean = cond_mean_list, 
                                  cond_var = cond_var_list, 
                                  k, 
                                  max_stp = 20, 
                                  conv_thresh = 1e-5, 
                                  sp_bounds,
                                  starting_location_seed = NULL) {
  #need to fix seed
  est_list <- list()
  
  for (est in seq_len(times)) {
    est_list[[est]] <- em_hardassign(cond_mean = cond_mean, 
                                     cond_var = cond_var, 
                                     k = k, 
                                     max_stp = max_stp, 
                                     conv_thresh = conv_thresh, 
                                     sp_bounds = sp_bounds,
                                     starting_location_seed = starting_location_seed
    )
    message('finished ', est)
  }
  
  #extract the likelihood vectors from each EM run
  like_list <- lapply(est_list, function(x) x$like_vec)
  
  #find EM run with biggest likelihood
  top_like_run <- which.max(sapply(like_list, function(x) max(x, na.rm = TRUE)))
  
  #in the biggest likelihood EM run, find the max likelihood index
  max_iter <- which.max(like_list[[top_like_run]])
  
  return(
    list(
      top_run_index = top_like_run, #index of top EM run
      top_iter_index = max_iter, #index of iter with highest likelihood in top EM run
      top_lik = like_list[[top_like_run]][[max_iter]], #top likelihood
      top_x = est_list[[top_like_run]]$xvec_list[[max_iter]], #top coordinates estimate
      output_list = est_list #output from all the EM runs
    )
  )
  
}


em_softassign_wrapper <- function(times = 1,
                                  cond_mean = cond_mean_list, 
                                  cond_var = cond_var_list, 
                                  k, 
                                  max_stp = 20, 
                                  conv_thresh = 1e-5, 
                                  sp_bounds,
                                  starting_location_seed = NULL) {
  
  #recover()
  
  #need to fix seed
  est_list <- list()
  
  for (est in seq_len(times)) {
    est_list[[est]] <- em_softassign(cond_mean = cond_mean, 
                                     cond_var = cond_var, 
                                     k = k, 
                                     max_stp = max_stp, 
                                     conv_thresh = conv_thresh, 
                                     sp_bounds = sp_bounds,
                                     starting_location_seed = starting_location_seed
    )
    message('finished ', est)
  }
  
  #extract the likelihood vectors from each EM run
  like_list <- lapply(est_list, function(x) x$like_vec)
  
  #find EM run with biggest likelihood
  top_like_run <- which.max(sapply(like_list, function(x) max(x, na.rm = TRUE)))
  
  #in the biggest likelihood EM run, find the max likelihood index
  max_iter <- which.max(like_list[[top_like_run]])
  
  return(
    list(
      top_run_index = top_like_run, #index of top EM run
      top_iter_index = max_iter, #index of iter with highest likelihood in top EM run
      top_lik = like_list[[top_like_run]][[max_iter]], #top likelihood
      top_x = est_list[[top_like_run]]$xvec_list[[max_iter]], #top coordinates estimate
      output_list = est_list #output from all the EM runs
    )
  )
  
}


em_hardassign <- function(cond_mean, cond_var, k, max_stp, conv_thresh, sp_bounds, starting_location_seed = NULL) {
  
  #recover()
  
  #=======================
  #=== FUNCTION SET-UP ===
  #=======================
  
  tree_count <- length(cond_mean)
  
  #matrix of all pairwise location combos
  #group_pairs <- expand.grid(x1 = 1:K, x2 = 1:K)
  
  #initialize node membership vec
  group_vec <- rep(1, times = tree_count*2)
  
  #counter for location initialization
  init_counter <- 1
  
  #lists to store things
  zvec_list <- list() #node membership vectors
  xvec_list <- list() #location coordinate matrices
  
  #all_lik_list <- list()
  
  loglik_list <- list()
  loglik_list[[1]] <- NA
  
  
  #===============================
  #=== LOCATION INITIALIZATION ===
  #===============================
  
  repeat {
    set.seed(starting_location_seed)
    Xold <- as.vector(t(cbind(runif(k,min=sp_bounds$x[1],max=sp_bounds$x[2]),
                              runif(k,min=sp_bounds$y[1],max=sp_bounds$y[2]))))
    set.seed(NULL)
    if (k == 1) break
    
    group_info_list <- list()
    for (i in seq_len(tree_count)) {
      lik_info <- logL_groupmat(Xmat = matrix(Xold, ncol = 2, byrow = TRUE), 
                                group_count = k, 
                                cond_mean[[i]], 
                                cond_var[[i]],
                                normalize = FALSE)
      #find the node location pairing that maximizes the likelihood
      group_info_list[[i]] <- which(lik_info == max(lik_info), arr.ind = TRUE)[1,]
    }
    
    group_vec <- unlist(group_info_list) #vector version of node memberships
    
    init_counter <- init_counter + 1
    
    if (length(unique(group_vec)) == k) {
      break # Exit the loop if the condition is no longer met
    }
    
    if (init_counter >= 500) stop('Coordination initialization failed.')
  }
  
  message("Initialization attempts: ", init_counter)
  
  
  #initial param vals
  #zvec_list[[1]] <- group_vec
  zvec_list[[1]] <- NA
  xvec_list[[1]] <- Xold
  
  
  
  #================
  #=== EM ITERS ===
  #================
  stp <- 2
  
  while (stp < max_stp) {
    
    ### EXPECTATION: given location coordinates, what are the node location memberships? ###
    group_info_list <- list()
    # group_like <- list()
    for (i in seq_len(tree_count)) {
      lik_info <- logL_groupmat(Xmat = matrix(Xold, ncol = 2, byrow = TRUE),
                                group_count = k,
                                cond_mean[[i]],
                                cond_var[[i]],
                                normalize = FALSE)
      #find the node location pairing that maximizes the likelihood
      group_info_list[[i]] <- which(lik_info == max(lik_info), arr.ind = TRUE)[1,]
      
      # group_like[[i]] <- logL_groupmat(X = matrix(Xold, ncol = 2, byrow = TRUE), 
      #                                                        group_count = k, 
      #                                                        cond_mean[[i]], 
      #                                                        cond_var[[i]],
      #                                                      normalize = FALSE)
      # group_info_list[[i]] <- which(group_like[[i]] == max(group_like[[i]]), arr.ind = TRUE)[1,]
      
    }
    
    #all_lik_list[[stp]] <- group_like
    
    group_vec <- unlist(group_info_list) #vector version of node memberships
    
    
    ### MAXIMIZATION: DERIVE NEW GROUP COORDINATES USING THE WEIGHTS ###
    #I'm currently doing this by finding the generalized least squares solution
    #based on the node memberships found during the expectation step
    Xnew <- as.vector(weighted_least_squares(multigroup_design(max_group = k,  group_vec), 
                                             Sigma_list = cond_var, 
                                             mu = unlist(cond_mean)))
    
    
    #list of all log likelihoods using the new parameter values
    like_list <- logL_groupmat_multitree(Xvec = Xnew, 
                                         group_count = k, 
                                         cond_mean_list = cond_mean, 
                                         cond_var_list = cond_var, 
                                         normalize = FALSE)
    
    #add new parameter vals to x_vec and assign Xnew to Xold
    xvec_list[[stp]] <- Xnew
    Xold <- Xnew
    
    zvec_list[[stp]] <- group_vec
    
    #the observed data log-likelihood is calculated by 
    #(1) summing the likelihoods (not log likelihoods) within trees
    #(2) taking the log of the summed likelihood for each tree
    #(3) summing the log summed likelihood across trees
    loglik_list[[stp]] <- sum(unlist(lapply(like_list, function(x) {log(sum(exp(x)))}))) #sum(unlist(like_list)) #log likelihood
    
    message(paste0('iter ', stp, '; loglik = ', loglik_list[[stp]]))
    
    
    ### CHECKS FOR THE NEXT ITER ###
    # proceed to next iter if this was the first iter
    if (stp == 2) {
      stp <- stp + 1
      next
    }
    
    #if the current lik is less than the previous lik, throw an error
    if (loglik_list[[stp]] < loglik_list[[stp - 1]] ) {
      warning('The likelihood decreased!')
      
      # return(
      #   list(
      #     newlik = loglik_list[[stp]],
      #     oldlik = loglik_list[[stp - 1]],
      #     zvec_list = zvec_list,
      #     xvec_list = xvec_list
      #   )
      # )
      
    }
    
    #if the difference in log liks is below the stopping threshold, break out of loop
    loglik_dif <- loglik_list[[stp]] - loglik_list[[stp - 1]]
    if ( loglik_dif <= conv_thresh & loglik_dif > 0 ) {
      convergence <- TRUE
      break
    }
    #################################
    
    stp <- stp + 1
  }
  
  return(
    list(final_x = Xnew,
         final_loglik = loglik_list[[length(loglik_list)]],
         xvec_list = xvec_list,
         iter_count = stp,
         group_membership_list = zvec_list,
         like_vec = unlist(loglik_list))
  )
}


em_softassign <- function(cond_mean, 
                          cond_var, 
                          k, 
                          max_stp, 
                          conv_thresh, 
                          sp_bounds, 
                          starting_location_seed = NULL) {
  
  #recover()
  
  #=======================
  #=== FUNCTION SET-UP ===
  #=======================
  tree_count <- length(cond_mean) #number of trees
  
  #lists to store things
  weight_list <- list() #node membership vectors
  weight_list[[1]] <- NA
  
  xvec_list <- list() #location coordinate matrices
  loglik_list <- list()
  loglik_list[[1]] <- NA
  
  q_list <- list()
  
  convergence <- FALSE
  
  
  
  #===============================
  #=== LOCATION INITIALIZATION ===
  #===============================
  
  #set.seed(3543)
  #currently doing this in the simplest way possible: uniformly sampling across
  #the area defined by sp_bounds
  set.seed(starting_location_seed)
  Xold <- as.vector(t(cbind(runif(k, min = sp_bounds$x[1], max = sp_bounds$x[2]),
                            runif(k, min = sp_bounds$y[1], max = sp_bounds$y[2]))))
  set.seed(NULL)
  #initial param vals
  xvec_list[[1]] <- Xold
  
  
  
  #================
  #=== EM ITERS ===
  #================
  stp <- 2
  
  while(stp < max_stp){
    
    ### EXPECTATION: given location coordinates, what are the node location memberships? ###
    pair_w_list <- logL_groupmat_multitree(Xvec = Xold, #matrix(Xold, ncol = 2, byrow = TRUE), 
                                           group_count = k, 
                                           cond_mean_list = cond_mean, 
                                           cond_var_list = cond_var, 
                                           normalize = TRUE)
    
    
    ### MAXIMIZATION: DERIVE NEW GROUP COORDINATES USING THE WEIGHTS ###
    #I'm currently doing this by maximizing the product of likelihood calculated
    #across all possible latent variable (location membership) states and the
    #corresponding responsibilities
    #maximizing the expected complete-data log-likelihood
    m_val <- optim(
      par = Xold, #starting parameter values at Xold 
      fn = m_step,
      r_list= pair_w_list, 
      cond_mean_list = cond_mean, 
      cond_var_list = cond_var,
      group_count = k,
      method = "BFGS",
      control = c(maxit = 1000)
    )
    
    message('convergence: ', m_val$convergence)
    message('message: ', m_val$message)
    Xnew <- m_val$par
    
    qnew <- m_step(param = Xnew, 
                   r_list = pair_w_list, 
                   cond_mean_list = cond_mean, 
                   cond_var_list = cond_var, 
                   group_count = k)
    
    #The observed-data likelihood does the opposite
    #list of all log likelihoods using the new parameter values
    like_list <- logL_groupmat_multitree(Xvec = Xnew, #matrix(Xnew, ncol = 2, byrow = TRUE), 
                                         group_count = k, 
                                         cond_mean_list = cond_mean, 
                                         cond_var_list = cond_var, 
                                         normalize = FALSE)
    
    #add new parameter vals to x_vec and assign Xnew to Xold
    xvec_list[[stp]] <- Xnew
    q_list[[stp]] <- qnew
    Xold <- Xnew
    
    weight_list[[stp]] <- pair_w_list #list of pairwise responsibilities
    
    #the observed data log-likelihood is calculated by 
    #(1) summing the likelihoods (not log likelihoods) within trees
    #(2) taking the log of the summed likelihood for each tree
    #(3) summing the log summed likelihood across trees
    loglik_list[[stp]] <- sum(unlist(lapply(like_list, function(x) {log(sum(exp(x)))}))) #sum(unlist(like_list)) #log likelihood
    
    message(paste0('iter ', stp, '; loglik = ', loglik_list[[stp]]))
    
    
    ### CHECKS FOR THE NEXT ITER ###
    # proceed to next iter if this was the first iter
    if (stp == 2) {
      stp <- stp + 1
      next
    }
    
    #if the current lik is less than the previous lik, throw an error
    if (loglik_list[[stp]] < loglik_list[[stp - 1]] ) {
      message('The likelihood decreased. Halting estimation and outputting some info.')
      
      return(
        list(
          newlik = loglik_list[[stp]],
          oldlik = loglik_list[[stp - 1]],
          weight_list = weight_list,
          xvec_list = xvec_list,
          q_list = q_list
        )
      )
      
    }
    
    #if the difference in log liks is below the stopping threshold, break out of loop
    if ( (loglik_list[[stp]] - loglik_list[[stp - 1]]) <= conv_thresh) {
      convergence <- TRUE
      break
    }
    #################################
    
    stp <- stp + 1
  }
  
  return(
    list(final_x = Xnew, #final location estimate
         final_loglik = loglik_list[[length(loglik_list)]], #final loglik
         convergence = convergence,
         xvec_list = xvec_list, #location estimates from each iter
         weight_list = weight_list, #list of weights from each step
         like_vec = unlist(loglik_list), #log liks from each step of EM
         iter_count = stp, #number of iters,
         q_list = q_list
    )
  )
  
}


multigroup_design <- function(max_group, node_membership) {
  design_mat <- matrix(0, 
                       nrow = length(node_membership)*2, 
                       ncol = max_group*2)
  
  #based on the group membership
  #design matrix is organized as: x1, y1, x2, y2, x3, y3, x4, y4
  groupind <- unlist(lapply(node_membership, function(x) c((x - 1)*2 + 1, (x - 1)*2 + 2)))
  
  #design_mat[cbind(1:(length(node_membership)*2), groupind)] <- 1
  design_mat[cbind(seq_len(length(node_membership)*2), groupind)] <- 1
  
  return(design_mat)
}


m_step <- function(param, r_list, cond_mean_list, cond_var_list, group_count) {
  ##m_step for the soft assign EM algorithm where you multiply the likelihoods
  #by the corresponding responsibilities
  
  if (!is.vector(param) || length(param) != group_count*2) {
    stop("param must be a vector with length equal to two times the group_count (2 corresponds to x and y coordinates per group).")
  }
  
  Xmat <- matrix(param, ncol = 2, byrow = TRUE)
  
  log_like_list_r <- list()
  for (TREE in seq_along(cond_mean_list)) {
    log_like_list_r[[TREE]] <- logL_groupmat(Xmat = Xmat, 
                                             group_count = group_count, 
                                             cond_mean_list[[TREE]], 
                                             cond_var_list[[TREE]],
                                             normalize = FALSE)*r_list[[TREE]]
  }
  
  -sum(unlist(log_like_list_r))
}



normalize_loglik <- function(loglik_mat) {
  #NORMALIZE: normalize the likelihoods so that they sum to 1. Because we are
  #initially calculating log likelihoods, I am using the log-sum-exp approach to
  #avoid underflow issues
  #https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  m <- max(loglik_mat)
  m1 <- m + log(sum(exp(loglik_mat - m)))
  
  return(exp(loglik_mat - m1))
}



logL_groupmat_multitree <- function(Xvec, 
                                    group_count, 
                                    cond_mean_list, 
                                    cond_var_list, 
                                    normalize = FALSE) {
  
  if (!is.vector(Xvec) || length(Xvec) != group_count*2) {
    stop("Xvec must be a matrix with 2 columns and group_count rows.")
  }
  
  Xmat <- matrix(Xvec, ncol = 2, byrow = TRUE)
  output_list <- list()
  for (i in seq_along(cond_mean_list)) {
    
    #normalizing so that the responsibilities all sum to 1. This is simply done
    #by dividing each likelihood by the total sum of likelihoods
    output_list[[i]] <- logL_groupmat(Xmat = Xmat, 
                                      group_count = group_count, 
                                      cond_mean = cond_mean_list[[i]], 
                                      cond_var = cond_var_list[[i]],
                                      normalize = normalize)
  }
  return(output_list)
}




logL_groupmat <- function(Xmat, group_count, cond_mean, cond_var, normalize = FALSE) {
  #recover()
  if (!is.matrix(Xmat) || ncol(Xmat) != 2 || nrow(Xmat) != group_count) {
    stop("Xmat must be a matrix with 2 columns and group_count rows.")
  }
  
  #QUESTIONS:
  # - the conditional variance/covariance is slighly non-symmetrical. Why is this?
  #   Is this just a precision problem? Currently, I'm dealing with this by
  #   forcing the matrix to be symmetric, but this might be me ignoring a bigger
  #   issue with my conditional variance calculation.
  
  
  #matrix to hold loglik values. The matrix is organized with node 1 along the
  #rows and node 2 along the columns and the locations indexed by the column and row
  #indices (e.g., location 1 for node 1 is in row 1; location 1 for node 2 is column 1)
  #For k = 2, the matrix is organized as:
  ##  ________________________    ________________________
  ##| (node1-loc1; node2-loc1)    (node1-loc1; node2-loc2) |
  ##| (node1-loc2; node2-loc1)    (node1-loc2; node2-loc2) |
  ##  ________________________    ________________________
  
  like_mat <- matrix(data = 0, nrow = group_count, ncol = group_count)
  
  for (NODE1 in seq_len(group_count)) {
    for (NODE2 in seq_len(group_count)) {
      #I'm sacrificing concision here for the sake of readability. We are
      #evaluating the probability at c(x_node1, y_node1, x_node2, y_node2). A one
      #line version to sub into x argument of mvtnorm::dmvnorm is: as.vector(t(X[c(NODE1, NODE2),]))
      NODE1_coords <- Xmat[NODE1,]
      NODE2_coords <- Xmat[NODE2,]
      like_mat[NODE1, NODE2] <- mvtnorm::dmvnorm(x = c(NODE1_coords, NODE2_coords),
                                                 mean = cond_mean,
                                                 sigma = as.matrix(Matrix::forceSymmetric(cond_var)),
                                                 log=TRUE)
    }
  }
  
  #NORMALIZE: normalize the likelihoods so that they sum to 1. Because we are
  #initially calculating log likelihoods, I using the log-sum-exp approach to
  #avoid underflow issues
  #https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  if (isTRUE(normalize)) return(normalize_loglik(loglik_mat = like_mat))
  
  return(like_mat)
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

# em_hardassign <- function(cond_mean, cond_var, k, max_stp, conv_thresh, sp_bounds, starting_location_seed = NULL) {
#   
#   recover()
#   
#   #=================
#   #=== EM SET-UP ===
#   #=================
#   
#   tree_count <- length(cond_mean)
#   
#   #matrix of all pairwise location combos
#   #group_pairs <- expand.grid(x1 = 1:K, x2 = 1:K)
#   
#   #initialize node membership vec
#   group_vec <- rep(1, times = tree_count*2)
#   
#   #counter for location initialization
#   init_counter <- 1
#   
#   #lists to store things
#   zvec_list <- list() #node membership vectors
#   xvec_list <- list() #location coordinate matrices
#   
#   loglik_list <- list()
#   loglik_list[[1]] <- NA
#   
#   
#   #===============================
#   #=== LOCATION INITIALIZATION ===
#   #===============================
#   
#   repeat {
#     set.seed(starting_location_seed)
#     Xold <- as.vector(t(cbind(runif(k,min=sp_bounds$x[1],max=sp_bounds$x[2]),
#                               runif(k,min=sp_bounds$y[1],max=sp_bounds$y[2]))))
#     set.seed(NULL)
#     if (k == 1) break
#     
#     
#     #An "dummy" expectation step to make sure all groups are represented in the
#     #initial groupings. This initial coordinate generation step is repeated until
#     #this is achieved
#     group_info_list <- list()
#     for (i in seq_len(tree_count)) {
#       
#       lik_info <- logL_groupmat(Xmat = matrix(Xold, ncol = 2, byrow = TRUE), 
#                                 group_count = k, 
#                                 cond_mean[[i]], 
#                                 cond_var[[i]],
#                                 normalize = FALSE)
#       
#       group_info_list[[i]] <- which(lik_info == max(lik_info), arr.ind = TRUE)[1,]
#     }
#     
#     #vector version of node memberships
#     group_vec <- unlist(group_info_list)
#     
#     init_counter <- init_counter + 1
#     
#     if (length(unique(group_vec)) == k) {
#       break # Exit the loop if the condition is no longer met
#     }
#   }
#   
#   message("Initialization attempts: ", init_counter)
#   
#   
#   #initial param vals
#   #zvec_list[[1]] <- group_vec
#   zvec_list[[1]] <- NA
#   xvec_list[[1]] <- Xold
#   
#   
#   
#   #================
#   #=== EM ITERS ===
#   #================
#   stp <- 2
#   
#   while(stp < max_stp){
#     
#     ### EXPECTATION: given location coordinates, what are the node location memberships? ###
#     group_info_list <- list()
#     for (i in seq_len(tree_count)) {
#       
#       lik_info <- logL_groupmat(Xmat = matrix(Xold, ncol = 2, byrow = TRUE), 
#                                 group_count = k, 
#                                 cond_mean[[i]], 
#                                 cond_var[[i]],
#                                 normalize = FALSE)
#       
#       group_info_list[[i]] <- which(lik_info == max(lik_info), arr.ind = TRUE)[1,]
#     }
#     
#     group_vec <- unlist(group_info_list)
#     
#     
#     ### MAXIMIZATION: DERIVE NEW GROUP COORDINATES USING THE WEIGHTS ###
#     #I'm currently doing this by maximizing the product of likelihood calculated
#     #across all possible latent variable (location membership) states and the
#     #corresponding responsibilities
#     #maximizing the expected complete-data log-likelihood
#     Xnew <- matrix(weighted_least_squares(multigroup_design(max_group = k,  group_vec), 
#                                           Sigma_list = cond_var, 
#                                           mu = unlist(cond_mean)
#     ), ncol = 2, byrow = TRUE)
#     
#     #list of all log likelihoods using the new parameter values
#     like_list <- logL_groupmat_multitree(Xnew, 
#                                          group_count = k, 
#                                          cond_mean_list = cond_mean, 
#                                          cond_var_list = cond_var, 
#                                          normalize = FALSE)
#     
#     #add new parameter vals to x_vec and assign Xnew to Xold
#     xvec_list[[stp]] <- Xnew
#     Xold <- Xnew
#     
#     zvec_list[[stp]] <- group_vec
#     
#     #the observed data log-likelihood is calculated by 
#     #(1) summing the likelihoods (not log likelihoods) within trees
#     #(2) taking the log of the summed likelihood for each tree
#     #(3) summing the log summed likelihood across trees
#     loglik_list[[stp]] <- sum(unlist(lapply(like_list, function(x) {log(sum(exp(x)))}))) #sum(unlist(like_list)) #log likelihood
#     
#     message(paste0('iter ', stp, '; loglik = ', loglik_list[[stp]]))
#     
#     
#     ### CHECKS FOR THE NEXT ITER ###
#     # proceed to next iter if this was the first iter
#     if (stp == 2) {
#       stp <- stp + 1
#       next
#     }
#     
#     #if the current lik is less than the previous lik, throw an error
#     if (loglik_list[[stp]] < loglik_list[[stp - 1]] ) {
#       message('The likelihood decreased. Halting estimation and outputting some info.')
# 
#       return(
#         list(
#           newlik = loglik_list[[stp]],
#           oldlik = loglik_list[[stp - 1]],
#           zvec_list = zvec_list,
#           xvec_list = xvec_list
#         )
#       )
# 
#     }
# 
#     #if the difference in log liks is below the stopping threshold, break out of loop
#     if ( (loglik_list[[stp]] - loglik_list[[stp - 1]]) <= conv_thresh) {
#      convergence <- TRUE
#      break
#     }
#     #################################
#     
#     stp <- stp + 1
#   }
#   
#   return(
#     list(X = Xnew,
#          xvec_list = xvec_list,
#          iter_count = stp,
#          group_membership_list = zvec_list,
#          like_vec = unlist(loglik_list),
#          loglik = loglik_list[[length(loglik_list)]])
#   )
# }

# multigroup_design_softassign <- function(node_count, max_group) {
#   design_mat <- matrix(0, 
#                        nrow = length(node_count)*2, 
#                        ncol = max_group*2)
#   design_mat[seq(1, nrow(design_mat), by = 2),seq(1, ncol(design_mat), by = 2)] <- 1
#   design_mat[seq(2, nrow(design_mat), by = 2),seq(2, ncol(design_mat), by = 2)] <- 1
#   #design matrix is organized as: x1, y1, x2, y2, x3, y3, x4, y4
#   return(design_mat)
# }

# lnL <- function(X,nTrees,condMeans,condVars){
#   lnLs <- sapply(seq_len(nTrees),
#                  function(i){
#                    mvtnorm::dmvnorm(x=X,
#                                     mean=condMeans[i,],
#                                     sigma=diag(rep(condVars[i],2)),
#                                     log=TRUE)
#                  })
#   return(lnLs)
# }

# group_logL_single <- function(X,condMean,condVar){
#   return(
#     mvtnorm::dmvnorm(x=X,
#                      mean=condMean,
#                      sigma=as.matrix(Matrix::forceSymmetric(condVar)),
#                      log=TRUE)
#   )
# }


# logL_group <- function(X, group_pairs, condMean, condVar){
#   #recover()
#   like_vec <- vector(mode = 'numeric', nrow(group_pairs))
#   for (i in 1:nrow(group_pairs)) {
#     #as.vector(t(X[unlist(group_pairs[2,]),])) --> vector of coordinates from the locations
#     #indexed using the group_pairs table
#     like_vec[i] <- mvtnorm::dmvnorm(x=as.vector(t(X[unlist(group_pairs[i,]),])),
#                                     mean=condMean,
#                                     sigma=as.matrix(Matrix::forceSymmetric(condVar)),
#                                     log=TRUE)
#   }
#   
#   return(like_vec)
# }

# inferparams_kloc <- function(condMeans1,
#                              condMeans2,
#                              condVars1,
#                              condVars2,
#                              K,
#                              spBounds,
#                              maxStp=100,
#                              init = c('random', 'kpp'),
#                              conv_thresh = 1e-5){
#   #recover()
#   
#   mean_check <- sapply(list(condMeans1, condMeans2), nrow)
#   var_check <- sapply(list(condVars1, condVars2), length)
#   if ( length(unique(c(mean_check, var_check))) != 1)
#     stop('condMeans1, condMeans2, condVars1, and condVars2 all must be the same length.')
#   
#   nTrees <- nrow(condMeans1)
#   init <- match.arg(init)
#   
#   lik_list <- list()
#   Zvec_node1 <- list()
#   Zvec_node2 <- list()
#   
#   Xvec <- list() #vector("list", length=maxStp)
#   
#   #z_list <- replicate(K, rep(1,nTrees), simplify = FALSE)
#   Z_node1 <- rep(1,nTrees)
#   Z_node2 <- rep(1,nTrees)
#   
#   #initialize
#   init_iter <- 1
#   while( (length(unique(Z_node1)) != K | length(unique(Z_node2)) != K) | K == 1){
#     #while(any(sapply(z_list, function(x) length(unique(x)) == 1))) {
#     # choose random initial locations for X
#     
#     if (init == 'random') {
#       X <- cbind(runif(K,min=spBounds$x[1],max=spBounds$x[2]),
#                  runif(K,min=spBounds$y[1],max=spBounds$y[2]))
#     } else if (init == 'kpp') {
#       X <- init_coord_generator(condMeans1, condMeans2, condVars1,condVars2, K = K, power = 3)
#     }
#     
#     print(X)
#     # calculate the likelihood of X given the conditional means/variances
#     lik_mat1 <- apply(X, 1, function(x) lnL(x,nTrees,condMeans1,condVars1))
#     #Z1 <- ifelse((lnX1A-lnX1B) > 0,1,2)
#     Z_node1 <- apply(lik_mat1, 1, which.max)
#     
#     lik_mat2 <- apply(X, 1, function(x) lnL(x,nTrees,condMeans2,condVars2))
#     #Z2 <- ifelse((lnX2A-lnX2B) > 0,1,2)
#     Z_node2 <- apply(lik_mat2, 1, which.max)
#     
#     # assign each tree to one location or the other based on which has the higher likelihood
#     if (K == 1) break
#     init_iter <- init_iter + 1
#   }
#   message("Initialization attempts: ", init_iter)
#   print(X)
#   
#   #initial param vals
#   Xvec[[1]] <- X
#   Zvec_node1[[1]] <- Z_node1
#   Zvec_node2[[1]] <- Z_node2
#   
#   # the MLE of X_1 is the centroid of the conditional means that belong to group 1
#   #X[1,] <- colMeans(rbind(condMeans1[Z1==1,]*invVarWts1[Z1==1],
#   #                        condMeans2[Z2==1,]*invVarWts2[Z2==1]))
#   # X[1,] <- colMeans(rbind(condMeans1[Z1==1,]*condVars1[Z1==1],
#   #                         condMeans2[Z2==1,]*condVars2[Z2==1]))
#   
#   for (i in seq_len(K)) {
#     X[i,] <- colSums(rbind(condMeans1[Z_node1==i,], condMeans2[Z_node2==i,])*calc_inv_weight(c(condVars1[Z_node1==i], condVars2[Z_node2==i])))
#   }
#   
#   #param vals after first maximization step
#   Xvec[[2]] <- X
#   
#   Zvec_node1[[2]] <- Z_node1
#   Zvec_node2[[2]] <- Z_node2
#   
#   lik1 <- lnL_multigroup(X,nTrees,condMeans1,condVars1, Z_node1)
#   lik2 <- lnL_multigroup(X,nTrees,condMeans2,condVars2, Z_node2)
#   lik_list[[1]] <- NA
#   lik_list[[2]] <- sum(c(lik1, lik2))
#   
#   # the MLE of X_2 is the centroid of the conditional means that belong to group 2
#   ##X[2,] <- colMeans(rbind(condMeans1[Z1==2,]*invVarWts1[Z1==2],
#   ##                        condMeans2[Z2==2,]*invVarWts2[Z2==2]))
#   #X[2,] <- colSums(rbind(condMeans1[Z1==2,], condMeans2[Z2==2,])*calc_inv_weight(c(condVars1[Z1==2], condVars2[Z2==2])))
#   
#   stp <- 3
#   # repeat the above until you hit a max step limit
#   #	could also add something to diagnose convergence
#   while(stp < maxStp){
#     
#     lik_mat1 <- apply(X, 1, function(x) lnL(x,nTrees,condMeans1,condVars1))
#     #Z1 <- ifelse((lnX1A-lnX1B) > 0,1,2)
#     Z_node1 <- apply(lik_mat1, 1, which.max)
#     
#     lik_mat2 <- apply(X, 1, function(x) lnL(x,nTrees,condMeans2,condVars2))
#     #Z2 <- ifelse((lnX2A-lnX2B) > 0,1,2)
#     Z_node2 <- apply(lik_mat2, 1, which.max)
#     
#     
#     for (i in seq_len(K)) {
#       X[i,] <- colSums(rbind(condMeans1[Z_node1==i,], condMeans2[Z_node2==i,])*calc_inv_weight(c(condVars1[Z_node1==i], condVars2[Z_node2==i])))
#     }
#     
#     Zvec_node1[[stp]] <- Z_node1
#     Zvec_node2[[stp]] <- Z_node2
#     Xvec[[stp]] <- X
#     
#     lik1 <- lnL_multigroup(X,nTrees,condMeans1,condVars1, Z_node1)
#     lik2 <- lnL_multigroup(X,nTrees,condMeans2,condVars2, Z_node2)
#     lik_list[[stp]] <- sum(c(lik1, lik2))
#     
#     if ( (lik_list[[stp]] - lik_list[[stp - 1]]) <= conv_thresh)
#       break
#     
#     stp <- stp + 1
#   }
#   
#   
#   return(
#     list(X = X,
#          Z1 = Z_node1,
#          Z2 = Z_node2,
#          Xvec = Xvec,
#          Z1vec = Zvec_node1,
#          Z2vec = Zvec_node2,
#          iter_count = stp,
#          like_vec = unlist(lik_list),
#          loglik = lik_list[[length(lik_list)]])
#   )
# }
# 

# logL_multi <- function(X,nTrees = NULL,condMeans,condVars){
#   if (is.null(nTrees)) nTrees <- length(condMeans)
#   lnLs <- sapply(seq_len(nTrees),
#                  function(i, LOC, mean, var){
#                    mvtnorm::dmvnorm(x=LOC,
#                                     mean=mean[[i]],
#                                     sigma=as.matrix(Matrix::forceSymmetric(var[[i]])),
#                                     log=TRUE)
#                  }, LOC = X, mean = cond_mean_list, var = cond_var_list)
#   return(lnLs)
# }
# 
# logL_single <- function(X,condMean,condVar){
#   return(
#     mvtnorm::dmvnorm(x=X,
#                      mean=condMean,
#                      sigma=as.matrix(Matrix::forceSymmetric(condVar)),
#                      log=TRUE)
#   )
# }


