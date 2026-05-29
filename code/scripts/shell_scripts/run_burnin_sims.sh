#! /bin/sh

top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
seed_file="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/burnin_seed_info.txt"

total_burnin=$(tail -n +2 $seed_file | wc -l)
total_burnin_vec=$(seq 1 $total_burnin)

#echo total burnin vec: $total_burnin_vec

for a in $total_burnin_vec; do
  x=$(( $a + 1 )) #index but avoiding the header

  REP=$(awk -v LINE=$x 'NR == LINE {print $1}' $seed_file)
  SEED=$(awk -v LINE=$x 'NR == LINE {print $2}' $seed_file)
  BURNIN_RUNTIME=$(awk -v LINE=$x 'NR == LINE {print $3}' $seed_file)
  SIGMA_DISPERSE=$(awk -v LINE=$x 'NR == LINE {print $4}' $seed_file)
  CARRYING_CAPACITY=$(awk -v LINE=$x 'NR == LINE {print $5}' $seed_file)

  #echo "starting $x_final"

  #echo "SEED: ${SEED}"
  #echo "RUNTIME: ${RUNTIME}"
 #printf "SEED: %s\n" "$SEED"
 #echo "BURNIN_RUNTIME: ${BURNIN_RUNTIME}"
 #echo "SIGMA_DISPERSE ${SIGMA_DISPERSE}"
  suffix_name=REP${REP}_RUNTIME${BURNIN_RUNTIME}_DISP${SIGMA_DISPERSE}_CC${CARRYING_CAPACITY}

  sh run_single_burnin.sh \
     ${suffix_name} \
     burnin_sims \
     ${top_level_dir} \
     ${SEED} \
     ${BURNIN_RUNTIME} \
     ${SIGMA_DISPERSE} \
     ${CARRYING_CAPACITY}
done



#slim \
#-d "SUFFIX='${kval}_${ped_gen}_rep${rep}'" \
#-d "ouput_path='${top_level_dir}/$i'"
#-d SEED=${seed} \
#-d RUNTIME=${runtime} \
#${slim_script_dir}/burnin_sim.slim

#page 665 of manual
