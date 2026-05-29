#!/bin/bash

#make sure to have slim_work conda enviro activated


raw_tree_path="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/postburnin_sims"
#raw_tree_file="POSTBURNIN_BURNINREP1_RUNTIME10000_DISP1.80_CC0.15_K3_Kreps15.trees"
#raw_tree_file="POSTBURNIN_BURNINREP1_RUNTIME10000_DISP1.25_CC0.15_K3_Kreps16.trees"
#raw_tree_file="POSTBURNIN_BURNINREP1_RUNTIME8000_DISP1.25_CC0.11_K3_Kreps16.trees"
#raw_tree_file="POSTBURNIN_BURNINREP1_RUNTIME8000_DISP1.50_CC0.13_K2_Kreps16.trees"
#raw_tree_file="POSTBURNIN_BURNINREP1_RUNTIME8000_DISP1.50_CC0.15_K2_Kreps16.trees"
raw_tree_file="POSTBURNIN_BURNINREP1_RUNTIME9000_DISP2.0_CC0.17_K1_Kreps16.trees"

top_level_dir=/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/processed_treeseq

output_file_name=$(basename ${raw_tree_file} .trees)
mkdir ${top_level_dir}/${output_file_name}

python3 process_tree_file_commandline.py \
--tree_sep 30 \
--sample_size 250 \
--ts_upload_path ${raw_tree_path}/ \
--ts_name ${raw_tree_file} \
--output_file_path ${top_level_dir}/${output_file_name}/ \
--output_tree_path ${top_level_dir}/${output_file_name}/ \
--output_prefix $output_file_name \
--recomb_rate 1e-8 \
--ne 500 \
--seed 342342

#python3 test.py \
#-n 3 \
#-u 50 \
#--file_name "green.txt"


#python3 process_tree_file_commandline.py \
#--tree_sep 50 \
#--ts_upload_path "/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_3pop/" \
#--ts_name "example_sim_admix_threepop.trees" \
#--output_file_path "/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_3pop/" \
#--output_tree_path "/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_3pop/" \
#--output_prefix "example_sim_admix_threepop" \
#--recomb_rate 1e-8 \
#--ne 500
