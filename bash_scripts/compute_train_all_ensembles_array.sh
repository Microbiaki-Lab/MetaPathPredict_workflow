#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=ensmbl # Job name
#SBATCH --ntasks=36 # Run on 36 CPUs
#SBATCH --cpus-per-task=1
#SBATCH --mem=187gb # Job memory request
#SBATCH --time=1-12:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/ensemble_training_%A-%a.log# Standard output/error
#SBATCH --array=282,283,284,285,286,289,291,292,293,294,299,300,301,302,303,304,306,307,309,310,311,312,313,314,315,316,317,318,319,320,323,325,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,400,403,405,406,407,408,409,410,411,412,413,414,415,416,419,421
#SBATCH --qos=unlim
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
mv /vortexfs1/home/dgellermcgrath/ensemble/rdata/ensemble_models/*_"$SLURM_ARRAY_TASK_ID"_ensemble_model.rda /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models
#
done

