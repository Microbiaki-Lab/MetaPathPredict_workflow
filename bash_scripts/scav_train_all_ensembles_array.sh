#!/bin/bash
#SBATCH --partition=scavenger # Queue selection
#SBATCH --job-name=ensmbl # Job name
#SBATCH --ntasks=36 # Run on 36 CPUs
#SBATCH --cpus-per-task=1
#SBATCH --mem=187gb # Job memory request
#SBATCH --time=1-12:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/ensemble_training_%A-%a.log# Standard output/error
#SBATCH --array=29,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,138,139,150,151,152,153,155,162,163,174,185,192,193,227,228,230,232,237,253,257,266,267,268,269,272,275,305,308,322,326,327,328,331,334,424
#SBATCH --qos=scavenger
#
cd /vortexfs1/home/dgellermcgrath/ensemble/rdata/old_splits
#
module load anaconda/5.1
source activate newr
#
for ITER in $SLURM_ARRAY_TASK_ID
do TRAIN=$(ls *train*_$ITER.RDA)
TEST=$(ls *test*_$ITER.RDA)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/create_model_env_vars.R --train $TRAIN --test $TEST -n $SLURM_ARRAY_TASK_ID
#
VARS=$(realpath /vortexfs1/home/dgellermcgrath/ensemble/rdata/ensemble_env_vars/*_"$SLURM_ARRAY_TASK_ID"_ensemble_env_vars.rdata)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/create_nnet_tune.R -f $VARS -n $SLURM_ARRAY_TASK_ID
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/create_xgboost_tune.R -f $VARS -n $SLURM_ARRAY_TASK_ID
#
NNET=$(realpath /vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_nnet_tunefile_grid50_compute.rda)
XGBOOST=$(realpath /vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_xgboost_tunefile_grid50_compute.rda)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/train_ensemble.R -n $NNET -x $XGBOOST
#
#
#
#
mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/ensemble_env_vars/*_"$SLURM_ARRAY_TASK_ID"_ensemble_env_vars.rdata /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars
#
mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_nnet_tunefile_grid50_compute.rda /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files
#
mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/leftover_training_files/*_"$SLURM_ARRAY_TASK_ID"_nnet_leftover_env_vars_grid50.rdata /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/leftover_training_files
#
mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_xgboost_tunefile_grid50_compute.rda /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files
#
mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/leftover_training_files/*_"$SLURM_ARRAY_TASK_ID"_xgboost_leftover_env_vars_grid50.rdata /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/leftover_training_files
#
#mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/ensemble_models/*_"$SLURM_ARRAY_TASK_ID"_ensemble_model.rds /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models
#
done

