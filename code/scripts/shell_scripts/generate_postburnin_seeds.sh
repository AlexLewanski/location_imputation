#! /bin/sh

top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
burnin_seed_file="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/burnin_seed_info.txt"
seed_file_name=postburnin_seed_info.txt

if [[ -f "$top_level_dir/$seed_file_name" ]]
then
    echo "Seed file already exists. Delete and re-run if you actually want to create a a new file."
    exit
fi

k_reps=4
k_vals=$(seq 2 4)
burnin_count=$(tail -n +2 $burnin_seed_file | wc -l)

printf BURNIN_REP"\t"K"\t"REP"\t"SEED"\n" > $top_level_dir/$seed_file_name

RANDOM=98331
for BURN in $(seq 1 $burnin_count); do
  for K in $k_vals; do
      printf $BURN"\t"$K"\t"$k_reps"\t"$RANDOM"\n" >> $top_level_dir/$seed_file_name
  done
done


##printf BURNIN_REP"\t"K"\t"REP"\t"SEED"\n" > $top_level_dir/$seed_file_name

#RANDOM=98331
#for BURN in $(seq 1 $burnin_count); do
#  for K in $k_vals; do
#    for i in $(seq 1 $k_reps); do
#      printf $BURN"\t"$K"\t"$i"\t"$RANDOM"\n" >> $top_level_dir/$seed_file_name
#    done
#  done
#done
