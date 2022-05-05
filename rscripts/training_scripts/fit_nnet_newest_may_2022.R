library(tidyverse)
library(tidymodels)
library(recipeselectors)
library(furrr)
library(butcher)

#set parallel processing parameters
plan(multicore, workers = 50) #49000 * 1024 ^2 - 49Gb per core, 50 cores total
options(future.globals.maxSize = 51380224000) #set maximum memory per parallal worker

#tuning/hyperparameter values & training performance metrics (auc_roc, kappa)
tune_files = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/*nnet*', intern = TRUE)

#model names, extracted
model_names = str_replace(tune_files, '\\/.*\\/+(M\\d{5})_\\d+.*', '\\1')

#other various important environmental variables from model training, excluding tune tibble
leftover_vars = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/leftover_training_files/*nnet*',
                       intern = TRUE)

env_vars = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/*', intern = TRUE)

#filter just to keep var files for which a tune file exists
leftover_vars = leftover_vars[str_detect(leftover_vars, paste0(model_names, '.*nnet*', collapse = '|'))]
env_vars = env_vars[str_detect(env_vars, paste0(model_names, collapse = '|'))]

############ check obj sizes before starting fittting
# get_size = function(.data) {
#   print(paste0(.data, ': ',
#                (get(.data) %>%
#                   object.size() %>%
#                   format(units = "Mb"))
#   ))
# }
#
# walk(ls(), ~ get_size(.x))
###################################################

#function to fit model
fit_model = function(leftoverVars, tuneFile, modelName, envVars) {

  load(leftoverVars) # load the env variables
  load(envVars)

  nnet_tune = readRDS(tuneFile) # load tune file

  # select best hyperparameters based on highest Kohen's Kappa result
  best_kap <- select_best(nnet_tune, "kap")

  #finalize the tidymodels workflow
  final_nnet <- finalize_workflow(
    nnet_workflow,
    best_kap
  )

  #fit the final model
  model <- fit(final_nnet, training_data_with_downsampled_observations_incl)

  model <- butcher::axe_call(model)
  model <- butcher::axe_ctrl(model)
  model <- butcher::axe_data(model)
  model <- butcher::axe_env(model)

  #save the model
  saveRDS(model, file = paste0(
    '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/butchered_nnet_models/',
    modelName, '_nnet_model.rds'
  ))

  #saveRDS(model, file = paste0(
  #  '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/nnet_models/',
  #  modelName, '_nnet_model.rds'
  #))
}

#wrap function in trycatch to record any errors that occur in any loop iteration
safely_fit = safely(fit_model)

#run model fits in parallel
results = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    leftover_vars,
    tune_files,
    model_names,
    env_vars
  ), ~
    safely_fit(..1, ..2, ..3, ..4)
)

