

# trying to improve models w/ severely imbalanced class distributions


## to determine the weight for the model classes:

 # the majority class weight will be: whatever integer number
 # gets you from the true class ratio to the class ratio you get
 # by downsampling the majority class

 # the minority class weight will remain 1 (I think??) because
 # it is not being up/downsampled or otherwise modified



# an example of how to add case_weights to a model workflow


#training_sim <-
#  training_sim %>%
#  mutate(
#    case_wts = ifelse(class == "class_1", 60, 1),
#    case_wts = importance_weights(case_wts)
#  )

## specify the model recipe, tuning specifications, etc.

## then add the case_weights arg to the model's workflow

# lr_wflow <-
#   workflow() %>%
#   add_model(lr_spec) %>%
#   add_recipe(sim_rec) %>%
#   add_case_weights(case_wts)   # <<- added here


#### full workflow - working on sample model M00038 - one that is highly imbalanced and coudl not be stacked w/ stacks

# load reqd libraries
library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors)
library(furrr)

#set furrr parameters
#plan(multicore, workers = 10) #16000 * 1024 ^2 - 16Gb per core, 10 cores total
#options(future.globals.maxSize = 16777216000)


set.seed(224)

# load model environmental variables that were pre-computed
load('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/M00038_29_ensemble_env_vars.rdata')



# create a table with class counts ['KEGG module present' class = 1, 'module absent' class = 0]
class_counts = table(training_data_with_downsampled_observations_incl$y)

# create named vector for the majority class; its name will be either '1' or '0'
maj_class = max(class_counts)
names(maj_class) = names(which(class_counts == max(class_counts)))

# create named vector for the minority class; its name will be again either '1' or '0'
min_class = min(class_counts)
names(min_class) = names(which(class_counts == min(class_counts)))

# the true majority class:minority class ratio before any downsampling
true_majority_to_minority_class_ratio = maj_class / min_class # 102/1

# encode the majority:minority class ratio we would like to implement on the training data
  # via random downsampling of majority class observations
underRatio = 10

# the majority class:minority class ratio we would like to achieve with random downsampling
  # of the majority class in the training data
downsampled_majority_to_minority_class_ratio = underRatio

# calculate the weight that will be given to the majority class after downsampling
majority_class_weight = true_majority_to_minority_class_ratio / downsampled_majority_to_minority_class_ratio


# add a importance weights column to the training data
training_data_with_downsampled_observations_incl = training_data_with_downsampled_observations_incl %>%
  mutate(
    case_wts = case_when(
      names(maj_class) == '0' ~ if_else(y == 0, majority_class_weight, 1),
      names(maj_class) == '1' ~ if_else(y == 1, majority_class_weight, 1)
    ),
    case_wts = importance_weights(case_wts)
  )


# neural network receipe specification; includes downsampling (is it in the appropriate order?)
nnet_rec <-
  recipe(y ~ ., data = vctrs::vec_slice(training_data_with_downsampled_observations_incl, 0)) %>%
  #update_role(everything()) %>%
  #update_role(y, new_role = 'outcome') %>%
  step_zv(all_predictors()) %>% # will remove features with zero variance
  step_select_infgain(all_predictors(), outcome = "y", top_p = tune()) %>%
  themis::step_downsample(y, under_ratio = underRatio) %>%
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
             MaxNWts = 8000,
             weights = case_wts)

nnet_workflow <- #create a workflow for this model
  workflow() %>%
  add_recipe(nnet_rec) %>%
  add_model(nnet_spec) %>%
  add_case_weights(case_wts) # ADD case weights column to workflow

#changed from control_grid()
ctrl <- control_grid(verbose = TRUE, # be chatty about the fitting
                     event_level = "second",  # the '1' level, which indicates an event for this model, is the second factor level
                     parallel_over = "everything", # CHANGE BACK TO RESAMPLES for compute nodes!!!
                     save_pred = TRUE,
                     save_workflow = TRUE
)

nnet_grid <- grid_latin_hypercube( # create a space-covering tuning grid that should adequately cover the hyperparamter tuning space
  hidden_units(range = c(1, 20)),
  penalty(),
  dropout(),
  activation(),
  top_p(range = c(10, 300)),
  size = 5 # CHANGE BACK TO 50 AFTER TESTING
)

# get_size = function(.data) {
#  print(paste0(.data, ': ',
#               (get(.data) %>%
#                  object.size() %>%
#                  format(units = "Mb"))
#  ))
# }

rm(training_data_with_downsampled_observations_incl)

## add weight column to the rsamp_method object; it will be needed

for (j in 1:10) {
  rsamp_method$splits[[j]]$data = rsamp_method$splits[[j]]$data %>%
    mutate(
      case_wts = case_when(
        names(maj_class) == '0' ~ if_else(y == 0, majority_class_weight, 1),
        names(maj_class) == '1' ~ if_else(y == 1, majority_class_weight, 1)
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

nnet_tune <- # tune the model ! should get notifications from parallel workers as tuning takes place
  tune_grid(nnet_workflow,
            resamples = rsamp_method,
            grid = nnet_grid,
            metrics = metric_set(kap, roc_auc),
            control = ctrl)

stopCluster(cl) # close the cluster


## errors/warnings from doing running the above code - same error for all folds:

# > nnet_tune$.notes[[1]]$note[1]
# [1] "`terms_select()` was deprecated in recipes 0.1.17.\nPlease use `recipes_eval_select()` instead."
# > nnet_tune$.notes[[1]]$note[2]
# [1] "\033[1m\033[33mError\033[39m in \033[1m\033[1m`check_case_weights()`:\033[22m\n\033[33m!\033[39m Case weights are not enabled by the underlying model implementation."































