library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors, lib.loc='/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)
library(optparse)

plan(multicore, workers = 10) #16000 * 1024 ^2 - 16Gb per core, 10 cores total
options(future.globals.maxSize = 16777216000)


#set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-f', '--file'), action = 'store',
              help = 'Path to environmental variables rdata file', type = 'character'),
  make_option(c('-n', '--number'), action = 'store',
              help = 'Path to input dataset', type = 'integer')
)

#parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

#load environmental variables
load(argv$file)

nnet_rec <-
  recipe(vctrs::vec_slice(training_data_with_downsampled_observations_incl, 0)) %>%
  update_role(everything()) %>%
  update_role(y, new_role = 'outcome') %>%
  step_zv(all_predictors()) %>% # will remove features with zero variance
  step_select_infgain(all_predictors(), outcome = "y", top_p = tune()) %>%
  step_normalize(all_predictors()) # center and scale all features based on training data

nnet_spec <- # tuning specifications for the model
  mlp(
    mode = "classification",
    hidden_units = tune(),
    penalty = tune(),
    dropout = tune(),
    epochs = 1000,
    activation = tune()
  ) %>%
  set_engine("nnet",
             MaxNWts = 5000)

nnet_workflow <- #create a workflow for this model
  workflow() %>%
  add_recipe(nnet_rec) %>%
  add_model(nnet_spec)

#changed from control_grid()
ctrl <- control_grid(verbose = TRUE, # be chatty about the fitting
                     event_level = "second",  # the '1' level, which indicates an event for this model, is the second factor level
                     parallel_over = "resamples",
                     save_pred = TRUE,
                     save_workflow = TRUE
)

nnet_grid <- grid_latin_hypercube( # create a space-covering tuning grid that should adequately cover the hyperparamter tuning space
  hidden_units(range = c(1, 20)),
  penalty(),
  dropout(),
  activation(),
  top_p(range = c(10, 300)),
  size = 50
)

#get_size = function(.data) {
#  print(paste0(.data, ': ',
#               (get(.data) %>%
#                  object.size() %>%
#                  format(units = "Mb"))
#  ))
#}

rm(training_data_with_downsampled_observations_incl)

#walk(ls(), ~ get_size(.x))

library(doMC)
options(cores = 5)
cl = makeForkCluster(outfile = '')
registerDoMC()

set.seed(1234)

nnet_tune <- # tune the model ! should get notifications from parallel workers as tuning takes place
  tune_grid(nnet_workflow,
            resamples = rsamp_method,
            grid = nnet_grid,
            metrics = metric_set(kap, roc_auc),
            control = ctrl)

stopCluster(cl) # close the cluster

saveRDS(nnet_tune, file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/', module, '_', argv$number, '_nnet_tunefile_grid50_compute.rda'))
rm(nnet_tune)
save(list = ls(), file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/leftover_training_files/', module, '_', argv$number, '_nnet_leftover_env_vars_grid50.rdata'))

