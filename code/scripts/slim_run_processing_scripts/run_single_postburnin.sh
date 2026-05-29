#! /bin/sh

slim_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/code/slim"

#test_vec="1 3"
#slim \
#-d "SUFFIX='TEST'" \
#-d "OUTDIR='postburnin_sims'" \
#-d "OUTPATH='/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output'" \
#-d SEED=251 \
#-d LOCATION_COUNT=2 \
#-d "FOUNDER_COUNT_VEC='${test_vec}'" \
#-d CROSS_COUNT=8 \
#${slim_dir}/run_from_burnin_multicross.slim

SUFFIX=$1
OUTDIR=$2
OUTPATH=$3
BURNIN_TREE_NAME=$4
BURNIN_TREE_PATH=$5
SEED=$6
LOCATION_COUNT=$7
FOUNDER_COUNT_VEC=$8
CROSS_COUNT=$9
WITHIN_LOCATION_SAMPLE_EXTENT=${10}
BETWEEN_LOCATION_MIN_DIST=${11}
START_TICK=${12}
SIGMA_DISPERSE=${13}
CARRYING_CAPACITY=${14}

#echo "++++++++"
#echo "run_single_postburnin"
#echo "${FOUNDER_COUNT_VEC}"
#echo "++++++++"

slim \
-d "SUFFIX='${SUFFIX}'" \
-d "OUTDIR='${OUTDIR}'" \
-d "OUTPATH='${OUTPATH}'" \
-d "BURNIN_TREE_NAME='${BURNIN_TREE_NAME}'" \
-d "BURNIN_TREE_PATH='${BURNIN_TREE_PATH}'" \
-d SEED=${SEED} \
-d LOCATION_COUNT=${LOCATION_COUNT} \
-d "FOUNDER_COUNT_VEC='${FOUNDER_COUNT_VEC}'" \
-d CROSS_COUNT=${CROSS_COUNT} \
-d WITHIN_LOCATION_SAMPLE_EXTENT=${WITHIN_LOCATION_SAMPLE_EXTENT} \
-d BETWEEN_LOCATION_MIN_DIST=${BETWEEN_LOCATION_MIN_DIST} \
-d START_TICK=${START_TICK} \
-d sigma_disp=${SIGMA_DISPERSE} \
-d K=${CARRYING_CAPACITY} \
${slim_dir}/run_from_burnin_multicross.slim
