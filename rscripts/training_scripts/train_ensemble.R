library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors)#, lib.loc='/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)
library(butcher)
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

## worked running it like this!!! perhaps allow_par is defaulted as TRUE and was causing the node to run out of memory...
ensemble_model =
  data_stack %>%
  blend_predictions(control = control_grid(allow_par = FALSE))

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

saveRDS(ensemble_model, file = paste0('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models/', module, '_ensemble_model.rds'))

