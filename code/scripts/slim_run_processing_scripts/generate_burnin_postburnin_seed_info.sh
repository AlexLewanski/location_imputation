#!/bin/bash

################
### OVERVIEW ###
################
#THIS SCRIPT GENERATES THE INFO TO RUN THE BURNIN AND POSTBURNIN SIMULATIONS. IT ALSO OUTPUTS A FILE
#THAT CONTAINS THE "BASENAME" OF THE POSTBURNIN TREE FILES, WHICH SUBSEQUENTLY IS USED FOR PROCESSING
#OF THE TREE FILES. 


#####################
### SCRIPT SET-UP ###
#####################

### PATHS AND FILE NAMES ####
top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"

burnin_seed_file_name=burnin_seed_info.txt
postburnin_seed_file_name=postburnin_seed_info.txt
base_file_name=basename_file.txt


### VARIABLES FOR BURNIN SIMS ###
burnin_reps=2

### VARIABLES FOR POSTBURNIN SIMS ###
k_reps=4
k_vals=($(seq 2 4))

### VARIABLES FOR BOTH BURNIN AND POSTBURNIN SIMS ###
BURNIN_RUNTIME=(5000)
SIGMA_DISPERSE=(1.5 2.0 2.5)
CARRYING_CAPACITY=(0.15)


### CHECKS TO MAKE SURE THE FILES DON'T ALREADY EXIST ###
if [[ -f "$top_level_dir/$burnin_seed_file_name" ]]
then
    echo "Burnin seed file already exists. Delete and re-run if you actually want to create a new file."
    exit
fi


if [[ -f "$top_level_dir/$postburnin_seed_file_name" ]]
then
    echo "Postburnin seed file already exists. Delete and re-run if you actually want to create a new file."
    exit
fi



#######################
### BURNIN SIM FILE ###
#######################

printf REP"\t"SEED"\t"BURNIN_RUNTIME"\t"SIGMA_DISPERSE"\t"CARRYING_CAPACITY"\n" > "$top_level_dir/$burnin_seed_file_name"

RANDOM=506087
for BURN in "${BURNIN_RUNTIME[@]}"; do
  for DISP in "${SIGMA_DISPERSE[@]}"; do
    for CAR in "${CARRYING_CAPACITY[@]}"; do
      for i in $(seq 1 "$burnin_reps"); do
        printf "%s\t%s\t%s\t%s\t%s\n" "$i" "$RANDOM" "$BURN" "$DISP" "$CAR" >> "$top_level_dir/$burnin_seed_file_name"

      done
    done
  done
done



###########################
### POSTBURNIN SIM FILE ###
###########################
#burnin_count=$(tail -n +2 $top_level_dir/$burnin_seed_file_name | wc -l)

printf BURNIN_REP"\t"K"\t"REP"\t"SEED"\t"BURNIN_RUNTIME"\t"SIGMA_DISPERSE"\t"CARRYING_CAPACITY"\n" > "$top_level_dir/$postburnin_seed_file_name"

RANDOM=98331

for BURN in "${BURNIN_RUNTIME[@]}"; do
  for DISP in "${SIGMA_DISPERSE[@]}"; do
    for CAR in "${CARRYING_CAPACITY[@]}"; do
      for i in $(seq 1 "$burnin_reps"); do
        for K in "${k_vals[@]}"; do
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$i" "$K" "$k_reps" "$RANDOM" "$BURN" "$DISP" "$CAR" >> "$top_level_dir/$postburnin_seed_file_name"
        done
      done
    done
  done
done



#####################
### SIM NAME FILE ###
#####################
for BURN in "${BURNIN_RUNTIME[@]}"; do
  for DISP in "${SIGMA_DISPERSE[@]}"; do
    for CAR in "${CARRYING_CAPACITY[@]}"; do
      for i in $(seq 1 "$burnin_reps"); do
        for K in "${k_vals[@]}"; do
          printf  "%s%s%s%s%s%s%s%s%s%s%s%s\n" "POSTBURNIN_BURNINREP" "$i" "_RUNTIME" "$BURN" "_DISP" "$DISP" "_CC" "$CAR" "_K" "$K" "_Kreps" "$k_reps" >> "$top_level_dir/$base_file_name"
        done
      done
    done
  done
done

#POSTBURNIN_BURNINREP1_RUNTIME15000_DISP1.75_CC0.20_K3_Kreps5.trees




#for BURN in $(seq 1 $burnin_count); do
#  for K in $k_vals; do
#      printf $BURN"\t"$K"\t"$k_reps"\t"$RANDOM"\t"$BURNIN_RUNTIME"\t"$SIGMA_DISPERSE"\t"$CARRYING_CAPACITY"\n" >> "$top_level_dir/$postburnin_seed_file_name"
#  done
#done
