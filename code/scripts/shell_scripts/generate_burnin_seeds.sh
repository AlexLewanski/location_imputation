#! /bin/sh

burnin_reps=3
top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
seed_file_name=burnin_seed_info.txt

if [[ -f "$top_level_dir/$seed_file_name" ]]
then
    echo "Seed file already exists. Delete and re-run if you actually want to create a a new file."
    exit
fi

printf REP"\t"SEED"\n" > $top_level_dir/$seed_file_name

RANDOM=506087
for i in $(seq 1 $burnin_reps); do
  printf $i"\t"$RANDOM"\n" >> $top_level_dir/$seed_file_name
done

