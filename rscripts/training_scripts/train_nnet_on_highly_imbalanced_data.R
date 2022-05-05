

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
plan(multicore, workers = 10) #16000 * 1024 ^2 - 16Gb per core, 10 cores total
options(future.globals.maxSize = 16777216000)



x = readRDS('/vortexfs1/home/dgellermcgrath/ensemble/rdata/old_splits/M00038-train_data_29.RDA')
y = readRDS('/vortexfs1/home/dgellermcgrath/ensemble/rdata/old_splits/M00038-test_data_29.RDA')


# format the dataset
df = x %>%
  bind_rows(y) %>%
  select(-source) %>%
  rename(y = 1) %>%
  mutate(y = if_else(y == '-1', '0', y)) %>% #check that  this step works...
  mutate(y = factor(y, levels = c('0', '1')))

#rm unneeded objects
rm(x, y)

#import rcpp fn for downsampling train and test splits later
Rcpp::sourceCpp('/vortexfs1/home/dgellermcgrath/cpp_fns_files/rarefy.cpp')

set.seed(224)

split <- initial_split(df, strata = "y", prop = 0.80)
train <- training(split)
test <- testing(split)

rm(df, split)

train_metadata = select(train, c(y, genome_name))
train = select(train, -c(y, genome_name))


if (all(map_chr(train, ~ class(.x)) == 'integer')) {
  train = as.matrix(train)
} else {
  stop('Error: Not all columns were integer class.')
}

create_downsampled_train_from = function(.train, .metadata) {
  future_map(
    .options = furrr_options(seed = TRUE),
    seq(1, 0.1, by = -0.1), ~
      rarefy(.train, sample_rate = .x) %>%
      as_tibble() %>%
      bind_cols(.metadata) %>%
      relocate(c(y, genome_name), .before = 1) %>%
      mutate(prop_counts_retained = .x, .after = genome_name)
  )
}

training_data_with_downsampled_observations_incl =
  create_downsampled_train_from(train, train_metadata)

rm(train)

training_data_with_downsampled_observations_incl =
  training_data_with_downsampled_observations_incl %>%
  map_dfr(~ .x) #%>%
#slice_sample(prop = 0.10)

# the metadata columns of the training data including downsampled observations
multi_train_metadata = training_data_with_downsampled_observations_incl %>%
  select(c(y, genome_name, prop_counts_retained))

# remove some cols
training_data_with_downsampled_observations_incl =
  training_data_with_downsampled_observations_incl %>%
  select(-c(genome_name, prop_counts_retained)) # keep the y column - the dependent variable - kegg module pres/abs

# 10-fold cross-validation
#stratifying the cv folds - balancing the amt of data from each downsampling increment in each fold
rsamp_method = vfold_cv(
  training_data_with_downsampled_observations_incl %>%
    mutate(ds = multi_train_metadata$prop_counts_retained, .after = y),
  v = 5, strata = "ds")

# remove the ds column from the rsamp object folds
for (j in 1:5) {
  rsamp_method$splits[[j]]$data = select(rsamp_method$splits[[j]]$data, -ds)
}
rm(j)



set.seed(224)

# load model environmental variables that were pre-computed
#load('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/M00038_29_ensemble_env_vars.rdata')




# neural network receipe specification; includes downsampling (is it in the appropriate order?)
nnet_rec <-
  recipe(y ~ ., data = vctrs::vec_slice(training_data_with_downsampled_observations_incl, 0)) %>%
  #update_role(everything()) %>%
  #update_role(y, new_role = 'outcome') %>%
  step_zv(all_predictors()) %>% # will remove features with zero variance
  step_select_infgain(all_predictors(), outcome = "y", top_p = tune()) %>%
  #themis::step_downsample(y, under_ratio = underRatio) %>%
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
             MaxNWts = 8000)#,
             #weights = case_wts)

nnet_workflow <- #create a workflow for this model
  workflow() %>%
  add_recipe(nnet_rec) %>%
  add_model(nnet_spec) #%>%
#add_case_weights(case_wts) # ADD case weights column to workflow

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
  size = 50 # CHANGE BACK TO 50 AFTER TESTING
)

get_size = function(.data) {
  print(paste0(.data, ': ',
               (get(.data) %>%
                  object.size() %>%
                  format(units = "Mb"))
  ))
}

#rm(training_data_with_downsampled_observations_incl)

## add weight column to the rsamp_method object; it will be needed

walk(ls(), ~ get_size(.x))

library(doMC)
options(cores = 35)
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


## it worked


# reload env vars
#load('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_env_vars/M00038_29_ensemble_env_vars.rdata')

# training_data_with_downsampled_observations_incl = training_data_with_downsampled_observations_incl %>%
#   mutate(
#     case_wts = case_when(
#       names(maj_class) == '0' ~ if_else(y == 0, majority_class_weight, 1),
#       names(maj_class) == '1' ~ if_else(y == 1, majority_class_weight, 1)
#     ),
#     case_wts = importance_weights(case_wts)
#   )


# try fitting with best kappa; try auroc later
best_kap <- select_best(nnet_tune, "kap")

#finalize the tidymodels workflow


final_nnet <- finalize_workflow(
  nnet_workflow,
  best_kap
)

#fit the final model
model <- fit(final_nnet, training_data_with_downsampled_observations_incl)




test_and_get_conf_matrices_for = function(.test_data_path,
                                          .module_name,
                                          .model) {
  # load test_data and model
  .test_data = readRDS(.test_data_path)
  #.model = readRDS(.model_path)

  # add column of predictions to the test_data
  .test_data = .test_data %>%
    relocate(y, .after = genome_name) %>%
    mutate(preds = predict(.model, .) %>% pull(1), .after = y) %>%
    select(c(prop_counts_retained, genome_name, y, preds))

  # create a confusion matrix for predictions
  cm = caret::confusionMatrix(.test_data$y, .test_data$preds, positive = '1')

  return(list(pred_md = .test_data, conf_mat = cm))
}

res = test_and_get_conf_matrices_for(
  '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/test_datasets/M00038_29_test_data_incl_downsampled_versions.rds',
  'M00038',
  model
)

strat = res$pred_md %>%
  group_by(prop_counts_retained) %>%
  group_split() %>%
  set_names(
    map_chr(., ~ paste0('prop_',
                        unique(.x$prop_counts_retained)
    ))) %>%
  map2_dfr(names(.), function(.x, .name) {
    caret::confusionMatrix(.x$y, .x$preds, positive = '1') %>%
      keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
      flatten() %>%
      as_tibble() %>%
      mutate(prop_counts_retained = .name, module_name = 'M00038') %>%
      relocate(c(prop_counts_retained, module_name,
                 F1, Precision, Recall,
                 Specificity, `Balanced Accuracy`),
               .before = 1)
  })



cms = res$pred_md %>%
  group_by(prop_counts_retained) %>%
  group_split() %>%
  set_names(
    map_chr(., ~ paste0('prop_',
                        unique(.x$prop_counts_retained)
    ))) %>%
  map2(names(.), function(.x, .name) {
    caret::confusionMatrix(.x$y, .x$preds, positive = '1')
  })
