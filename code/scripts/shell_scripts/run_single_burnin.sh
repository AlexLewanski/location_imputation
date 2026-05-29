#! /bin/sh

slim_dir="/Users/alexlewanski/Documents/michigan_state/research/location_imputation/code/slim"

SUFFIX=$1
OUTDIR=$2
OUTPATH=$3
SEED=$4
RUNTIME=$5
SIGMA_DISPERSE=$6
CARRYING_CAPACITY=$7

echo "CARRYING_CAPACITY ${CARRYING_CAPACITY}"

slim \
-d "SUFFIX='${SUFFIX}'" \
-d "OUTDIR='${OUTDIR}'" \
-d "OUTPATH='${OUTPATH}'" \
-d SEED=${SEED} \
-d RUNTIME=${RUNTIME} \
-d sigma_disp=${SIGMA_DISPERSE} \
-d K=${CARRYING_CAPACITY} \
${slim_dir}/burnin_sim.slim

#-d SUFFIX="${kval}_${ped_gen}_rep${rep}" \
#https://github.com/MesserLab/SLiM/issues/32
