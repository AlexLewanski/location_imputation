#!/bin/bash

#make sure to have slim_work conda enviro activated

basename_file="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/basename_file.txt"
raw_tree_path="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/postburnin_sims"
output_dir=/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/processed_treeseq

sample_size=(50 150)

while IFS= read -r FILE_NAME; do
    raw_tree_file=${FILE_NAME}.trees

    if [ ! -f "${raw_tree_path}/$raw_tree_file" ]; then
        echo "The tree file does not exist!"
        exit 1
    else
        echo "Tree file exists! Proceeding with processing."
    fi

    for SAMP_SIZE in "${sample_size[@]}"; do

        #echo $raw_tree_file

        output_file_name=${FILE_NAME}_sampsize${SAMP_SIZE}
        #output_file_name=$(basename ${raw_tree_file} .trees)_sampsize${SAMP_SIZE}
        mkdir ${output_dir}/${output_file_name}

        python3 process_tree_file_commandline_alt.py \
        --tree_sep 50 \
        --sample_size $SAMP_SIZE \
        --ts_upload_path ${raw_tree_path}/ \
        --ts_name ${raw_tree_file} \
        --output_file_path ${output_dir}/${output_file_name}/ \
        --output_tree_path ${output_dir}/${output_file_name}/ \
        --output_prefix $output_file_name \
        --recomb_rate 1e-8 \
        --ne 500 \
        --seed $RANDOM
    done
done < ${basename_file}




#for FILE in core_names; do
#    for SAMP_SIZE in "${sample_size[@]}"; do

#        output_file_name=$(basename ${raw_tree_file} .trees)_sampsize${SAMP_SIZE}
#        mkdir ${top_level_dir}/${output_file_name}

#        python3 process_tree_file_commandline_alt.py \
#        --tree_sep 50 \
#        --sample_size 150 \
#        --ts_upload_path ${raw_tree_path}/ \
#        --ts_name ${raw_tree_file} \
#        --output_file_path ${top_level_dir}/${output_file_name}/ \
#        --output_tree_path ${top_level_dir}/${output_file_name}/ \
#        --output_prefix $output_file_name \
#        --recomb_rate 1e-8 \
#        --ne 500 \
#        --seed 342322
#    done
#done

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
