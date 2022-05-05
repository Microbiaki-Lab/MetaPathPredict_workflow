#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=rensmbl # Job name
#SBATCH --ntasks=36 # Run on 36 CPUs
#SBATCH --cpus-per-task=1
#SBATCH --mem=187gb # Job memory request
#SBATCH --time=1-12:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/redo_ensemble_training_%A-%a.log# Standard output/error
#SBATCH --array=105,107,108,111,127,135,148,156,159,160,161,170,171,182,186,187
#SBATCH --qos=unlim
#
cd /vortexfs1/home/dgellermcgrath/ensemble/rdata/old_splits
#
module load anaconda/5.1
source activate newr
#
NNET=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_nnet_tunefile_grid50_compute.rda)
XGBOOST=$(realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/*_"$SLURM_ARRAY_TASK_ID"_xgboost_tunefile_grid50_compute.rda)
#
Rscript /vortexfs1/home/dgellermcgrath/ensemble/rscripts/retrain_poor_ensemble_with_less_lasso_penalty.R -n $NNET -x $XGBOOST



#


#105,107,108,111,127,135,148,156,159,160,161,170,171,182,186,187,
#206,212,222,225,226,233,248,250,256,260,263,27,271,28,290,295,298,299,300,310,324,326,327,337,345,352,355,367,374,380,388,39,393,401,405,414,417,419,421,423,50,55,69,77,88,91,92,93