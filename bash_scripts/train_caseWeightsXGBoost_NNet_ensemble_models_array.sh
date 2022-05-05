#!/bin/bash
#SBATCH --partition=bigmem # Queue selection
#SBATCH --job-name=cwxg_nn_ens # Job name
#SBATCH --ntasks=5 # Run on 5 CPUs
#SBATCH --cpus-per-task=1
#SBATCH --mem=185gb # Job memory request
#SBATCH --time=1-12:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/LLP_weighted_XGBoost_NNet_ensemble_training_%A-%a.log# Standard output/error
#SBATCH --array=1-424
#SBATCH --qos=unlim
#
cd /vortexfs1/home/dgellermcgrath/ensemble/rdata/old_splits
#
module load anaconda/5.1
source activate newr
#
NNET=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_nnet_tunefile_grid50_compute.rda)
XGBOOST=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_tune_files/*_"$SLURM_ARRAY_TASK_ID"_xgboost_tunefile_grid50_compute.rda)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/retrain_poor_ensemble_with_less_lasso_penalty.R -n $NNET -x $XGBOOST