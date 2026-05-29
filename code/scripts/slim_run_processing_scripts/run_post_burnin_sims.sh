#! /bin/sh

#top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
#slim_script_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/code/slim"
sim_info="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/postburnin_seed_info.txt"
burnin_tree_path="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/burnin_sims"
output_path="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"

output_dir="postburnin_sims"
WITHIN_LOCATION_SAMPLE_EXTENT=5
BETWEEN_LOCATION_MIN_DIST=15

reps=$(tail -n +2 $sim_info | wc -l)
rep_vec=$(seq 1 $reps)

for i in $rep_vec; do
  
  REP=$(( $i + 1 )) #index but avoiding the header
  burnin_rep=$(awk -v REP=$REP 'NR == REP {print $1}' $sim_info)
  k_val=$(awk -v REP=$REP 'NR == REP {print $2}' $sim_info)
  cross_count=$(awk -v REP=$REP 'NR == REP {print $3}' $sim_info)
  seed_val=$(awk -v REP=$REP 'NR == REP {print $4}' $sim_info)
  START_TICK=$(awk -v REP=$REP 'NR == REP {print $5}' $sim_info)
  SIGMA_DISPERSE=$(awk -v REP=$REP 'NR == REP {print $6}' $sim_info)
  CARRYING_CAPACITY=$(awk -v REP=$REP 'NR == REP {print $7}' $sim_info)

  name=BURNINREP${burnin_rep}_RUNTIME${START_TICK}_DISP${SIGMA_DISPERSE}_CC${CARRYING_CAPACITY}_K${k_val}_Kreps${cross_count}
  #burnin_tree_name=BURNIN${burnin_rep}.trees
  burnin_tree_name=BURNINREP${burnin_rep}_RUNTIME${START_TICK}_DISP${SIGMA_DISPERSE}_CC${CARRYING_CAPACITY}.trees
  

  #echo "burnin_rep $burnin_rep"
  #echo "k_val $k_val"
  #echo "k_val_rep $k_val_rep"
  #echo "seed_val $seed_val"
  #echo "name $name"

  if [[ $k_val -eq 1 ]]; then
    location_input_vec="4"
  elif [[ $k_val -eq 2 ]]; then
    location_input_vec="2 2"
  elif [[ $k_val -eq 3 ]]; then
    location_input_vec="2 1 1"
  elif [[ $k_val -eq 4 ]]; then
    location_input_vec="1 1 1 1"
  else
    echo "The script only works with K of 2, 3, or 4 right now."
    exit 1
  fi

  #echo "K ${k_val}"
  #echo "location_input_vec ${location_input_vec}"
  #echo "cross_count ${cross_count}"
  #echo ${WITHIN_LOCATION_SAMPLE_EXTENT}
  #echo ${BETWEEN_LOCATION_MIN_DIST}

 #feed in the correct file
  sh run_single_postburnin.sh \
     ${name} \
     ${output_dir} \
     ${output_path} \
     ${burnin_tree_name} \
     ${burnin_tree_path} \
     ${seed_val} \
     ${k_val} \
     "${location_input_vec}" \
     ${cross_count} \
     ${WITHIN_LOCATION_SAMPLE_EXTENT} \
     ${BETWEEN_LOCATION_MIN_DIST} \
     ${START_TICK} \
     ${SIGMA_DISPERSE} \
     ${CARRYING_CAPACITY}

done
 
