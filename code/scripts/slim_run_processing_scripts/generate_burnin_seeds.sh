#! /bin/sh

# *** NOTE ***
#9/17/2026: THIS SCRIPT HAS BEEN SUPERSEDED BY A SCRIPT THAT GENERATES INFO FOR BOTH THE
#           BURNIN AND POSTBURNIN SIMS (generate_burnin_postburnin_seed_info.sh).

burnin_reps=1
top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
seed_file_name=burnin_seed_info.txt

if [[ -f "$top_level_dir/$seed_file_name" ]]
then
    echo "Seed file already exists. Delete and re-run if you actually want to create a new file."
    exit
fi

printf REP"\t"SEED"\t"BURNIN_RUNTIME"\t"SIGMA_DISPERSE"\t"CARRYING_CAPACITY"\n" > $top_level_dir/$seed_file_name

RANDOM=506087
BURNIN_RUNTIME=15000
SIGMA_DISPERSE=(0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0)
CARRYING_CAPACITY=(0.20)
for CAPAC in "${SIGMA_DISPERSE[@]}"; do
  for DISP in "${SIGMA_DISPERSE[@]}"; do
    for i in "$(seq 1 $burnin_reps)"; do
      printf $i"\t"$RANDOM"\t"$BURNIN_RUNTIME"\t"$SIGMA_DISPERSE"\t"$CARRYING_CAPACITY"\n" >> $top_level_dir/$seed_file_name
    done
  done
done


#printf REP"\t"SEED"\n" > $top_level_dir/$seed_file_name

#RANDOM=506087
#for i in $(seq 1 $burnin_reps); do
#  printf $i"\t"$RANDOM"\n" >> $top_level_dir/$seed_file_name
#done

#printf REP"\t"SEED"\n" > $top_level_dir/$seed_file_name

#RANDOM=506087
#for i in $(seq 1 $burnin_reps); do
#  printf $i"\t"$RANDOM"\n" >> $top_level_dir/$seed_file_name
#done

