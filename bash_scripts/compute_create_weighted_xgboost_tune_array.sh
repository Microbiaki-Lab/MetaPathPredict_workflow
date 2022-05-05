#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=c_xgb # Job name
#SBATCH --ntasks=36 # Run on 36 CPUs
#SBATCH --cpus-per-task=1
#SBATCH --mem=187gb # Job memory request
#SBATCH --time=1-12:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/create_weighted_xgboost_tune_%A-%a.log# Standard output/error
#SBATCH --array=1-80
#SBATCH --qos=unlim
#
cd /vortexfs1/home/dgellermcgrath/ensemble/rdata/old_splits
#
module load anaconda/5.1
source activate newr
#
VARS=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/*_"$SLURM_ARRAY_TASK_ID"_ensemble_env_vars.rdata)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/with_weights_create_xgboost_tune.R -f $VARS -n $SLURM_ARRAY_TASK_ID
