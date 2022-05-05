# trying to improve models w/ severely imbalanced class distributions

# load reqd libraries
library(tidymodels)
library(stacks, lib.loc = '/vortexfs1/home/dgellermcgrath/manually_dl_R_libs')
library(tidyverse)
library(recipeselectors)
library(furrr)
library(optparse)

plan(multicore, workers = 10) #16000 * 1024 ^2 - 16Gb per core, 10 cores total
options(future.globals.maxSize = 16777216000)

#set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-f', '--file'), action = 'store',
              help = 'Path to input dataset', type = 'character'),
  make_option(c('-n', '--number'), action = 'store',
              help = 'Task array ID number', type = 'integer')
)

#parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

#load environmental variables
load(argv$file)

# create a table with class counts ['KEGG module present' class = 1, 'module absent' class = 0]
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


# xgboost receipe specification; includes downsampling (is it in the appropriate order?)
xgboost_recipe <-
  recipe(y ~ ., data = vctrs::vec_slice(training_data_with_downsampled_observations_incl, 0)) %>%
  #update_role(everything()) %>%
  #update_role(y, new_role = 'outcome') %>%
  step_zv(all_predictors()) %>% # will remove features with zero variance
  step_select_infgain(all_predictors(), outcome = "y", top_p = tune()) #%>%
  #themis::step_downsample(y, under_ratio = underRatio)


xgboost_spec <- # tuning specifications for the model
  boost_tree(trees = tune(), # list parameters that will be tuned via repeated cross-validation
             mtry = tune(),
             min_n = tune(),
             tree_depth = tune(),
             learn_rate = tune(),
             loss_reduction = tune(),
             sample_size = tune()) %>%
  set_mode("classification") %>% # set method to classification
  set_engine("xgboost")

xgboost_workflow <- #create a workflow for this model
  workflow() %>%
  add_recipe(xgboost_recipe) %>%
  add_model(xgboost_spec) %>%
  add_case_weights(case_wts) # ADD case weights column to workflow

#changed from control_grid()
ctrl <- control_grid(verbose = TRUE, # be chatty about the fitting
                     event_level = "second",  # the '1' level, which indicates an event for this model, is the second factor level
                     parallel_over = "resamples", # CHANGE BACK TO RESAMPLES for compute nodes!!!
                     save_pred = TRUE,
                     save_workflow = TRUE
)

xgb_grid <- grid_latin_hypercube( # create a space-covering tuning grid that should adequately cover the hyperparamter tuning space
  trees(),
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), training_data_with_downsampled_observations_incl),
  learn_rate(),
  top_p(range = c(10, 300)),
  size = 50
)

# get_size = function(.data) {
#   print(paste0(.data, ': ',
#                (get(.data) %>%
#                   object.size() %>%
#                   format(units = "Mb"))
#   ))
# }

rm(training_data_with_downsampled_observations_incl)

## add weight column to the rsamp_method object; it will be needed
for (j in 1:5) {
  rsamp_method$splits[[j]]$data = rsamp_method$splits[[j]]$data %>%
    mutate(
      case_wts = case_when(
        names(maj_class) == '0' ~ if_else(y == 0, maj_class, min_class),
        names(maj_class) == '1' ~ if_else(y == 1, maj_class, min_class)
      ),
      case_wts = importance_weights(case_wts)
    )
}
rm(j)

#walk(ls(), ~ get_size(.x))

library(doMC)
options(cores = 5)
cl = makeForkCluster(outfile = '')
registerDoMC()

set.seed(1234)

xgboost_tune <- # tune the model ! should get notifications from parallel workers as tuning takes place
  tune_grid(xgboost_workflow,
            resamples = rsamp_method,
            grid = xgb_grid,
            metrics = metric_set(kap, roc_auc),
            control = ctrl)

stopCluster(cl) # close the cluster

# save tune file then env without the tune file
saveRDS(xgboost_tune, file = paste0('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_tune_files/', module, '_', argv$number, '_xgboost_tunefile_grid50_compute.rda'))
rm(xgboost_tune)
save(list = ls(), file = paste0('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/weighted_leftover_training_files/', module, '_', argv$number, '_xgboost_leftover_env_vars_grid50.rdata'))
