#################################################
#################################################
### GENERATING PEDIGREES FOR SLIM SIMULATIONS ###
#################################################
#################################################

library(here)

generate_annual_simple_ped <- function(input_founders = setNames(c(5, 5), c("M", "F")),
                                       input_type = c('founder_info', 'function_generate'),
                                       consider_sex = TRUE, #not doing anything with this yet
                                       repro_cycles = 5,
                                       mean_repro = 3,
                                       generate_death_info = TRUE, #not doing anything with this yet
                                       mating_sucess_prob = 1,
                                       mating_option = 'F_choose_M', #not doing anything with this yet
                                       runaway = 100000,
                                       seed = NULL) {
  
  #column_names <- c('id', 'fid', 'mid', 'sex')
  set.seed(seed)
  if (input_type == 'founder_info') {
    input_founders$tick <- 0
    pedigree_list <- list(input_founders[,c('id', 'fid', 'mid', 'sex','tick')])
    
  } else if (input_type == 'function_generate') {
    pedigree_list <- list(
      data.frame(id = seq_len(sum(input_founders)), #paste0('ID', seq_len(sum(input_founders)), '_0'),
                 fid = 0,
                 mid = 0,
                 sex = as.integer(rep(names(input_founders), input_founders) == 'F'),
                 tick = 0)
      
    )
  }
  

  for (CYCLE in seq_len(repro_cycles)) {
    total_female_vec <- pedigree_list[[CYCLE]]$id[pedigree_list[[CYCLE]]$sex == 1]
    total_male_vec <- pedigree_list[[CYCLE]]$id[pedigree_list[[CYCLE]]$sex == 0]
      
    if (length(total_female_vec) == 0 | length(total_male_vec) == 0) {
        #attempt <- attempt + 1
        #next
      stop("Pedigree generation is halted. One or both sexes no longer exist.")
    }
      
    female_list <- list()
    male_list <- list()
      #this could be made more efficient by just generating the probabilities in one go
      #I'm keeping this purposefully inefficient for future, more complicated mating choices
    list_index <- 1
    for (i in seq_along(total_female_vec)) {
      if (runif(1, 0, 1) < mating_sucess_prob) {
        female_list[[list_index]] <- total_female_vec[i]
        male_list[[list_index]] <- sample(total_male_vec, 1)
        list_index <- list_index + 1
      }
    }
      
    repro_count_vec <- rpois(length(female_list), mean_repro)
      
    ### Some checks before proceeding with pedigree generation###
    #were there any mating events?
    if (length(female_list) == 0) {
      stop('No reproduction events occurred. Halting pedigree generation.')
    }
        
      
    #Are any offspring produced from the mating events?
    if (sum(repro_count_vec) == 0) {
      stop('No offspring were generated. Halting pedigree generation.')
    }
      
    #Does the number of offspring meet or exceed the maximum count threshold?
    if (sum(repro_count_vec) >= runaway) {
      stop('Offspring count hit the runaway threshold. Halting pedigree generation.')
    }
      
    female_mate_vec <- rep(unlist(female_list), times = repro_count_vec)
    male_mate_vec <- rep(unlist(male_list), times = repro_count_vec)
      
    pedigree_list[[CYCLE + 1]] <- data.frame(id = max(pedigree_list[[CYCLE]]$id) + seq_along(female_mate_vec),# paste0('ID', seq_along(female_mate_vec), '_', CYCLE),
                                              fid = male_mate_vec,
                                              mid = female_mate_vec,
                                              sex = sample(c(0, 1), 
                                                           size = length(female_mate_vec), 
                                                           replace = TRUE, 
                                                           prob = c(0.5, 0.5)),
                                               tick = CYCLE
    )
  }
  
  pedigree_final <- do.call(rbind, pedigree_list)
  pedigree_final[,c('id', 'fid', 'mid')] <- pedigree_final[,c('id', 'fid', 'mid')]*-1
  
  return(
    list(
      pedigree = pedigree_final,
      death_info = death_info <- data.frame(id = pedigree_final$id,
                                            tick = pedigree_final$tick + 1)
    )
  )
}


generate_annual_simple_ped_wrapper <- function(input_founders = setNames(c(15, 15), c("M", "F")),
                                               input_type = 'function_generate',
                                               repro_cycles = 15,
                                               mean_repro = 1.8,
                                               generate_death_info = TRUE,
                                               mating_sucess_prob = 0.95,
                                               generation_attempts = 1,
                                               seed = NULL) {
  
  attempts <- 1
  ped_output <- 'error'
  while(attempts <= generation_attempts & identical(ped_output, 'error')) {
     ped_output <- tryCatch(generate_annual_simple_ped(input_founders = input_founders,
                                                      input_type = input_type,
                                                      repro_cycles = repro_cycles,
                                                      mean_repro = mean_repro,
                                                      generate_death_info = generate_death_info,
                                                      mating_sucess_prob = mating_sucess_prob,
                                                      seed = seed),
                           error = function(c) "error"
     )
     attempts <- attempts + 1
   }

  if (identical(ped_output, 'error')) message('Pedigree unable to be generated')
  message("Number of pedigree generation attempts: ", attempts - 1)
  return(ped_output)
}


# test_ped_sim <- generate_annual_simple_ped_wrapper(input_founders = setNames(c(2, 2), c("M", "F")), 
#                                            input_type = 'function_generate',
#                                            repro_cycles = 2,
#                                            mean_repro = 2,
#                                            generate_death_info = TRUE,
#                                            mating_sucess_prob = 1,
#                                            generation_attempts = 100,
#                                            seed = 21)


#two full siblings in 3rd gen
ped_2gen <- data.frame(
  id =   c(1, 2, 3, 4, 5, 6, 7, 8),
  par1 = c(0, 0, 0, 0, 1, 3, 5, 5),
  par2 = c(0, 0, 0, 0, 2, 4, 6, 6),
  sex =  c(0, 1, 0, 1, 0, 1, 0, 0),
  tick = c(0, 0, 0, 0, 1, 1, 2, 2)
)

death_2gen <- data.frame(id = ped_2gen$id,
                         tick = ped_2gen$tick + 1)


write.table(ped_2gen[-nrow(ped_2gen),],
            file = here('code', 'slim', 'ped_sim_files', 'gen3_mate.txt'),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = '\t')

write.table(death_2gen[-nrow(death_2gen),],
            file = here('code', 'slim', 'ped_sim_files', 'gen3_death.txt'),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = '\t')




# library(pedtools)
# library(kinship2)
# 
# plot(
#   pedigree(ped_2gen$id, 
#            ped_2gen$par1, 
#            ped_2gen$par2,
#            sex = ped_2gen$sex)
# )

# write.table(test_ped_sim$pedigree, 
#             file = here('code', 'slim', 'ped_sim_files', 'test_mate.txt'), 
#             row.names = FALSE, 
#             col.names = FALSE, 
#             quote = FALSE, 
#             sep = '\t')
# 
# write.table(test_ped_sim$death_info, 
#             file = here('code', 'slim', 'ped_sim_files', 'test_death.txt'), 
#             row.names = FALSE, 
#             col.names = FALSE, 
#             quote = FALSE, 
#             sep = '\t')


#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# input_founders = setNames(c(5, 5), c("M", "F")) 
# consider_sex = TRUE #not doing anything with this yet
# repro_cycles = 5
# mean_repro = 3
# generate_death_info = TRUE
# seed = NULL
# mating_sucess_prob = 1
# 
# pedigree_list <- list(
#   data.frame(id = paste0('ID', seq_len(sum(input_founders)), '_0'),
#              fid = '0',
#              mid = '0',
#              sex = rep(names(input_founders), input_founders),
#              tick = 0)
#               
# )
# 
# for (CYCLE in seq_len(repro_cycles)) {
#   total_female_vec <- pedigree_list[[CYCLE]]$id[pedigree_list[[CYCLE]]$sex == 'F']
#   total_male_vec <- pedigree_list[[CYCLE]]$id[pedigree_list[[CYCLE]]$sex == 'M']
#   
#   if (length(total_female_vec) == 0 | length(total_male_vec) == 0) {
#     cat("Pedigree generation is halted. Either there are no males or females.")
#     pedigree_final <- do.call(rbind, pedigree_list)
#     return(
#       list(
#         pedigree = pedigree_final,
#         death_info = data.frame(id = pedigree_final$id,
#                                 tick = pedigree_final$tick + 1)
#       )
#     )
#   }
#   
#   female_list <- list()
#   male_list <- list()
#   #this could be made more efficient by just generating the probabilites in one go
#   #I'm keeping this purposefully inefficient for future, more complicated mating choices
#   list_index <- 1
#   for (i in seq_along(total_female_vec)) {
#     if (runif(1, 0, 1) < mating_sucess_prob) {
#       female_list[[list_index]] <- total_female_vec[i]
#       male_list[[list_index]] <- sample(total_male_vec, 1)
#       list_index <- list_index + 1
#     }
#   }
#   
#   repro_count_vec <- rpois(length(female_list), mean_repro)
#   female_mate_vec <- rep(unlist(female_list), times = repro_count_vec)
#   male_mate_vec <- rep(unlist(male_list), times = repro_count_vec)
#   
#   pedigree_list[[CYCLE + 1]] <- data.frame(id = paste0('ID', seq_along(female_mate_vec), '_', CYCLE),
#                                            fid = male_mate_vec,
#                                            mid = female_mate_vec,
#                                            sex = sample(c("M", "F"), 
#                                                         size = length(female_mate_vec), 
#                                                         replace = TRUE, 
#                                                         prob = c(0.5, 0.5)),
#                                            tick = CYCLE
#                                            )
# }
# 
# pedigree_final <- do.call(rbind, pedigree_list)
# return(
#   list(
#     pedigree = pedigree_final
#     death_info = death_info <- data.frame(id = pedigree_final$id,
#                                           tick = pedigree_final$tick + 1)
#   )
# )
# pedigree_final <- do.call(rbind, pedigree_list)
# return(
#   list(
#     pedigree = pedigree_final,
#     death_info = data.frame(id = pedigree_final$id,
#                             tick = pedigree_final$tick + 1)
#   )
# )

# test <- 1
# x <- 1
# 
# while(TRUE) {
#   new_var1 <- 1000
#   break
# }
# 
# print('HERE')
