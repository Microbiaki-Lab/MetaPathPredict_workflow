library(tidyverse)
library(tidymodels)
library(recipeselectors)
library(furrr)
library(butcher)

#set parallel processing parameters
plan(multicore, workers = 50) #49000 * 1024 ^2 - 49Gb per core, 50 cores total
options(future.globals.maxSize = 51380224000) #set maximum memory per parallal worker

#tuning/hyperparameter values & training performance metrics (auc_roc, kappa)
tune_files = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_tune_files/*xgboost*', intern = TRUE)

#model names, extracted
model_names = str_replace(tune_files, '\\/.*\\/+(M\\d{5})_\\d+.*', '\\1')

#other various important environmental variables from model training, excluding tune tibble
leftover_vars = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_leftover_training_files/*xgboost*', intern = TRUE)

# contains training data to fit the model
ensemble_env_vars = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/*', intern = TRUE)

#filter just to keep ensemble var files for which a tune file exists
ensemble_env_vars = ensemble_env_vars[str_detect(
  ensemble_env_vars,
  paste0(model_names, collapse = '|'))]


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
fit_model = function(leftoverVars, ensembleEnvVars, tuneFile, modelName) {

  load(ensembleEnvVars) # load ensemble env vars - contains needed training data
  load(leftoverVars) # load the env variables
  xgboost_tune = readRDS(tuneFile) # load tune file


  ####################################################################################
  class_counts = table(training_data_with_downsampled_observations_incl$y)

  # create named vector for the majority class; its name will be either '1' or '0'
  maj_class = max(class_counts)
  names(maj_class) = names(which(class_counts == max(class_counts)))

  # create named vector for the minority class; its name will be again either '1' or '0'
  min_class = min(class_counts)
  names(min_class) = names(which(class_counts == min(class_counts)))

  #
  maj_class = (maj_class + min_class) / (2 * maj_class) # majority class weight
  min_class = (maj_class + min_class) / (2 * min_class) # minority class weight

  # add a importance weights column to the training data
  training_data_with_downsampled_observations_incl =
    training_data_with_downsampled_observations_incl %>%
    mutate(
      case_wts = case_when(
        names(maj_class) == '0' ~ if_else(y == 0, maj_class, min_class),
        names(maj_class) == '1' ~ if_else(y == 1, maj_class, min_class)
      ),
      case_wts = importance_weights(case_wts)
    )
  ####################################################################################


  # select best hyperparameters based on highest Kohen's Kappa result
  best_kap <- select_best(xgboost_tune, "kap")

  #finalize the tidymodels workflow
  final_xgb <- finalize_workflow(
    xgboost_workflow,
    best_kap
  )

  #fit the final model
  model <- fit(final_xgb, training_data_with_downsampled_observations_incl)

  model <- butcher::axe_call(model)
  model <- butcher::axe_ctrl(model)
  model <- butcher::axe_data(model)
  model <- butcher::axe_env(model)

  #save the model
  saveRDS(model, file = paste0(
    '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_butchered_xgboost_models/',
    modelName, '_weighted_xgboost_model.rds'
  ))
}

#wrap function in trycatch to record any errors that occur in any loop
safely_fit = safely(fit_model)

#run model fits in parallel
results = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    leftover_vars,
    ensemble_env_vars,
    tune_files,
    model_names
  ), ~
    safely_fit(..1, ..2, ..3, ..4)
)



# fit models that finished later on
finished_models = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_butchered_xgboost_models/*', intern = TRUE)
finished_models = str_replace(finished_models, '\\/.*\\/+(M\\d{5})_.*', '\\1')

#remaining tune files
tune_files = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_tune_files/*xgboost*', intern = TRUE)

tune_files = tune_files[!str_detect(
  tune_files,
  paste0(finished_models, collapse = '|'))]

#remaining model names
model_names = str_replace(tune_files, '\\/.*\\/+(M\\d{5})_\\d+.*', '\\1')

model_names = model_names[!str_detect(
  model_names,
  paste0(finished_models, collapse = '|'))]

#remaining leftover_vars
leftover_vars = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_leftover_training_files/*xgboost*', intern = TRUE)

leftover_vars = leftover_vars[!str_detect(
  leftover_vars,
  paste0(finished_models, collapse = '|'))]

# remaining ensemble_env_vars
ensemble_env_vars = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/*', intern = TRUE)

ensemble_env_vars = ensemble_env_vars[str_detect(
  ensemble_env_vars,
  paste0(model_names, collapse = '|'))]

# fit the remaining models that require fitting (this step can be repeated)
results = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    leftover_vars,
    ensemble_env_vars,
    tune_files,
    model_names
  ), ~
    safely_fit(..1, ..2, ..3, ..4)
)











# evaluate xgboost models (fit with case weights) on held-out test --------
#int
#conda activate newr
#R

library(tidyverse)
library(tidymodels)
library(recipeselectors)
library(furrr)

#set parallel processing parameters
plan(multicore, workers = 50) #49000 * 1024 ^2 - 49Gb per core, 50 cores total
options(future.globals.maxSize = 51380224000) #set maximum memory per parallal worker

#testing fn before applying
#x = readRDS('weighted_butchered_xgboost_models/M00511_weighted_xgboost_model.rds')
#y = readRDS('test_datasets/M00511_251_test_data_incl_downsampled_versions.rds')

model_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_butchered_xgboost_models/*', intern = TRUE)
testdata_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/test_datasets/*', intern = TRUE)

# appear to be only 423/424 testdata files; see which one is missing later..
model_names = str_replace(
  testdata_paths,
  '\\/.*\\/+(M\\d{5})_.*', '\\1'
)

# keep the 423 model paths that we have testdata paths for
model_paths = model_paths[str_detect(
  model_paths, paste0(model_names, collapse = '|')
)]

test_model = function(.model, .testdata, .model_name) {
  x = readRDS(.model)
  y = readRDS(.testdata)

  performance_metrics = y %>%
    group_by(prop_counts_retained) %>%
    group_split() %>%
    map_dfr(~ .x %>%
              mutate(
                preds = predict(x, .) %>%
                  pull(1)
              ) %>%
              select(y, preds, prop_counts_retained) %>%
              conf_mat(., truth = y, preds) %>%
              summary(event_level = 'second') %>%
              mutate(prop_counts_ret = .x %>%
                       pull(prop_counts_retained) %>%
                       unique())
    ) %>%
    select(-.estimator) %>%
    pivot_wider(
      names_from = .metric,
      values_from = .estimate
    ) %>%
    relocate(c(f_meas, precision, recall, spec), .after = 1) %>%
    mutate(model_name = .model_name, .before = 1)

  return(performance_metrics)
}

weighted_xgboost_performance = future_pmap_dfr(
  .options = furrr_options(seed = TRUE),
  .progress = TRUE,
  list(
    model_paths,
    testdata_paths,
    model_names
  ), ~
    test_model(..1, ..2, ..3)
)

saveRDS(weighted_xgboost_performance, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/weighted_butchered_xgboost_model_performances_on_heldout_downsampled_protein_annotation_testdata_051622.rda')








































