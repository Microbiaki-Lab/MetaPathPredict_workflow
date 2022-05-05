library(tidymodels)
library(stacks, lib.loc = '/vortexfs1/home/dgellermcgrath/manually_dl_R_libs')
library(tidyverse)
library(recipeselectors)
library(furrr)
library(butcher)
library(doParallel)
library(optparse)


#set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-n', '--nnetFile'), action = 'store',
              help = 'Path to nnet tune file', type = 'character'),
  make_option(c('-x', '--xgboostFile'), action = 'store',
              help = 'Path to xgboost tune file', type = 'character')
)

#parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

module = str_replace(argv$nnetFile, '\\/.*\\/+(M\\d{5}_\\d+)_nnet.*', '\\1') #FILL THIS IN

# load req'd tune files
xgboost_tune = readRDS(argv$xgboostFile)
nnet_tune = readRDS(argv$nnetFile)

data_stack <-
  stacks() %>%
  add_candidates(xgboost_tune) %>%
  add_candidates(nnet_tune)

rm(xgboost_tune, nnet_tune)

ensemble_model =
  data_stack %>%
  blend_predictions(
    control = control_grid(allow_par = TRUE),
    penalty = 10^(-6:-3)
  )

rm(data_stack)

# register a parallel backend to run blend_predictions and fit_ensemble w/ parallel processing
options(cores = 5)
registerDoParallel()

ensemble_model =
  ensemble_model %>%
  fit_members()

#reduce model size
ensemble_model =
  ensemble_model %>%
  axe_call() %>%
  axe_data() %>%
  axe_env() %>%
  axe_fitted()

if (length(ensemble_model$member_fits) == 0) {
  message('No ensemble members were retained in the ensemble model.')
}

#saveRDS(ensemble_model, file = paste0('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models_lower_lasso_penalty/', module, '_ensemble_model.rds'))

saveRDS(ensemble_model, file = paste0('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models_LLP_caseWeightsXGBoost_NNet/', module, '_ensemble_model.rds'))
