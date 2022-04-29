#!/bin/bash
#SBATCH --partition=scavenger # Queue selection
#SBATCH --job-name=stack # Job name
#SBATCH --ntasks=36 # Run on 36 CPUs
#SBATCH --cpus-per-task=1
#SBATCH --mem=187gb # Job memory request
#SBATCH --time=1-12:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/only_ensemble_fitting_%A-%a.log# Standard output/error
#SBATCH --array=6,7,8,9,10,13,16,17,23,24,30,31,32,33,34,35,36,37,38,40,41,43,44,46,47,48,49,50,51,52,214,215,216,218,219,221,222,223,224,225,226,227,228,229,230,231,232,235,237,252,253,254,266,267,268,269,270,272,273,275,276,277,278,279,280,281,282,283,284,285,286,289,291,301,304,311,315,316,317,319,320
#SBATCH --qos=scavenger
#
cd /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/
#
module load anaconda/5.1
source activate newr
#
for ITER in $SLURM_ARRAY_TASK_ID
#
do NNET=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_nnet_tunefile_grid50_compute.rda)
XGBOOST=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_xgboost_tunefile_grid50_compute.rda)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/train_ensemble.R -n $NNET -x $XGBOOST
#
done

