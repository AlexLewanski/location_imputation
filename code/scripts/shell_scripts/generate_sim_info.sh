#! /bin/sh

#####################
### SCRIPT SET-UP ###
#####################

burnin_reps=1
k_reps=16 #number of reps per combo of k value and burnin
k_vals=$(seq 1 1)
sim_runtime=9000
disperse=(2.0)
carrying_capacity=(0.17)


top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
burnin_seed_file_name=burnin_seed_info.txt
postburnin_seed_file_name=postburnin_seed_info.txt


if [[ -f "$top_level_dir/$seed_file_name" ]]
then
    echo "Burnin seed file already exists. Delete and re-run if you actually want to create a new file."
    exit
fi

if [[ -f "$top_level_dir/$postburnin_seed_file_name" ]]
then
    echo "Postburnin seed file already exists. Delete and re-run if you actually want to create a new file."
    exit
fi



#############################
### SEEDS FOR BURNIN SIMS ###
#############################
printf REP"\t"SEED"\t"BURNIN_RUNTIME"\t"sigma_disperse"\t"carrying_capacity"\n" > $top_level_dir/$burnin_seed_file_name

RANDOM=506087

for DISPERSE in "${disperse[@]}"; do
  for CARRYING_CAPACITY in "${carrying_capacity[@]}"; do
    for BURNIN_REP in $(seq 1 $burnin_reps); do
      printf $BURNIN_REP"\t"$RANDOM"\t"${sim_runtime}"\t"${DISPERSE}"\t"${CARRYING_CAPACITY}"\n" >> $top_level_dir/$burnin_seed_file_name
    done
  done
done



##################################
### SEEDS FOR POST-BURNIN SIMS ###
##################################

printf BURNIN_REP"\t"K"\t"REP"\t"SEED"\t"BURNIN_RUNTIME"\t"sigma_disperse"\t"carrying_capacity"\n" > $top_level_dir/$postburnin_seed_file_name

RANDOM=98331
for DISPERSE in "${disperse[@]}"; do
  for CARRYING_CAPACITY in "${carrying_capacity[@]}"; do
    for K in $k_vals; do
      for BURN in $(seq 1 $burnin_reps); do
          printf $BURN"\t"$K"\t"$k_reps"\t"$RANDOM"\t"${sim_runtime}"\t"${DISPERSE}"\t"${CARRYING_CAPACITY}"\n" >> $top_level_dir/$postburnin_seed_file_name
      done
    done
  done
done



#################################
### CODE NOT CURRENTLY IN USE ###
#################################
##printf BURNIN_REP"\t"K"\t"REP"\t"SEED"\t"BURNIN_RUNTIME"\n" > $top_level_dir/$seed_file_name

#RANDOM=98331
#for BURN in $(seq 1 $burnin_count); do
#  for K in $k_vals; do
#    for i in $(seq 1 $k_reps); do
#      printf $BURN"\t"$K"\t"$i"\t"$RANDOM"\n" >> $top_level_dir/$seed_file_name
#    done
#  done
#done
