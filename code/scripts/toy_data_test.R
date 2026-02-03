###################################################################
###################################################################
### TOY DATASET :testing an EM algorithm for location inference ###
###################################################################
###################################################################

library(here)
library(dplyr)
library(ggplot2)
library(ape)
library(phytools)
#library(treesliceR)

source(here('code', 'scripts', 'location_est_code.R'))
#spBounds <- list("x"=c(5,95),"y"=c(5,95))


##########

toy_dat_k5 <- generate_toy_data(ntrees = 400,
                  spatial_bounds = list("x"=c(5,95),"y"=c(5,95)), # spatial bounds for X and Y coordinates
                  K = 4,
                  min_cluster_dist = 10,
                  condmean1_sigma = 60,
                  condmean2_sigma = 50,
                  cond_var_min = c(3),
                  cond_var_max = c(10)
)


#set.seed(891)
toy_knode_infer <- inferparams_kloc(condMeans1=toy_dat_k5$condMeans1,
                               condMeans2=toy_dat_k5$condMeans2,
                               condVars1=toy_dat_k5$condVars1,
                               condVars2=toy_dat_k5$condVars2,
                               K=4,
                               #spBounds=spBounds_infer,
                               maxStp=150,
                               init = c('random', 'kpp')[2])

#set.seed(22)
iter_df <- lapply(1:(length(toy_knode_infer$Xvec) - 1), function(x, df) {
  x_df <- as.data.frame(df[[x]])
  x_df$loc <- as.character(1:nrow(x_df))
  x_df$iter <- x
  return(x_df)
}, df = toy_knode_infer$Xvec[1:(length(toy_knode_infer$Xvec) - 1)]) %>% 
  do.call(rbind,.)


as.data.frame(rbind(toy_dat_k5$condMeans1,toy_dat_k5$condMeans2)) %>% 
  ggplot() +
  geom_point(aes(x, y)) +
  geom_path(data = iter_df %>% 
              group_by(loc) %>% 
              arrange(iter),
            aes(V1, V2, group = loc, color = loc), linewidth = 2) +
  geom_point(data = iter_df, aes(V1, V2), color = 'blue', size = 3) +
  geom_point(data = as.data.frame(toy_dat_k5$centers), #true centers
             aes(x, y), size = 4, color = 'red') +
  geom_point(data = as.data.frame(toy_knode_infer$X),
             aes(x = V1, y = V2), size = 3, shape = 21, stroke = 2,
             color = 'green') +
  theme_bw()



#################################
### CODE NOT CURRENTLY IN USE ###
#################################
#***this code has been packaged up into a function (generate_toy_data) that
#***is currently stored in location_est_code.R

# #This isn't simulating under the statistical model, it's just coming up
# #with reasonable values that imitate the structure of the data
# 
# # number of trees
# nTrees <- 500
# # spatial bounds for X and Y coordinates
# spBounds <- list("x"=c(5,95),"y"=c(5,95))
# # number of locations and individual is "from"
# #	***NOTE: code will only work for K=2 right now***
# K <- 5
# # cluster membership for each node across all trees
# #Z1 <- t(replicate(nTrees,sample(c(0,1),K,replace=FALSE)))
# Z1 <- t(rmultinom(nTrees, size = 1, prob = rep(1,K)/K))
# 
# #Z2 <- t(replicate(nTrees,sample(c(0,1),2,replace=FALSE)))
# Z2 <- t(rmultinom(nTrees, size = 1, prob = rep(1,K)/K))
# 
# # simulate geographic locations that are far enough
# #	away from each other that inference is likely to work
# xdist <- matrix(0,nrow=K,ncol=K)
# 
# attempt <- 0
# # rejects any locations that are less than 30 distance units away from each other
# while(any(xdist[upper.tri(xdist,diag=FALSE)] < 25) & attempt < 100){
#   X <- cbind(runif(K,min=spBounds$x[1],max=spBounds$x[2]),
#              runif(K,min=spBounds$y[1],max=spBounds$y[2]))
#   colnames(X) <- c("x","y")
#   row.names(X) <- paste0("X_",1:K)
#   xdist <- fields::rdist(X)
#   attempt <- attempt + 1
# }
# 
# # simulate some "conditional means"
# #	not a procedural simulation, just pulling some numbers 
# #	from some normals centered around the true locations (X)
# #	so that they're likely to give good results
# condMeans1 <- Reduce("rbind",lapply(1:nTrees,
#                                     function(i){
#                                       MASS::mvrnorm(n=1,mu=X[which(Z1[i,]==1),],Sigma=diag(rep(50,2)))
#                                     }),init=NULL)
# condMeans2 <- Reduce("rbind",lapply(1:nTrees,
#                                     function(i){
#                                       MASS::mvrnorm(n=1,mu=X[which(Z2[i,]==1),],Sigma=diag(rep(55,2)))
#                                     }),init=NULL)
# 
# # simulate some random uniform values to imitate conditional variance values
# condVars1 <- runif(nTrees,min=0.1,max=10)
# condVars2 <- runif(nTrees,min=0.1,max=10)
