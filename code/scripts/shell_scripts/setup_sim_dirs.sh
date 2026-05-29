#! /bin/sh

top_level_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output"
slim_script_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/code/slim"

#create top level dir if it doesn't exist
if [[ ! -d ${top_level_dir} ]]
    then
        mkdir ${top_level_dir}
    fi

#create subdirs if they don't exist
for i in burnin_sims postburnin_sims
do
    if [[ ! -d ${top_level_dir}/$i ]]
    then
        mkdir ${top_level_dir}/$i
    fi

done
