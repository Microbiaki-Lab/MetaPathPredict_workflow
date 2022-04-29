#load libraries
library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors)#, lib.loc = '/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)
library(optparse)

#set furrr parallel process parameters
plan(multicore, workers = 10)
options(future.globals.maxSize = 16777216000) #16000 * 1024 ^2 - 16Gb max per worker, 10 cores total

#set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-t', '--train'), action = 'store',
              help = 'Path to input dataset', type = 'character'),
  make_option(c('-p', '--test'), action = 'store',
              help = 'Path to input dataset', type = 'character'),
  make_option(c('-n', '--number'), action = 'store',
              help = 'Slurm array task ID', type = 'integer')
)

#parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

#extract module name from input file name
module = str_replace(argv$train, '(M\\d{5})-.*', '\\1')

# creating XGBoost/nnet ensemble model stack ------------------------------

#read in dataset flatfile
#df = vroom::vroom(
#  file = argv$file,
#  col_types = paste0(c('ic'), paste0(rep('i', 10137), collapse = ''))
#)

x = readRDS(argv$train)
y = readRDS(argv$test)


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
  map_dfr(~ .x) %>%
  slice_sample(prop = 0.10)

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
  v = 10, strata = "ds")

# remove the ds column from the rsamp object folds
for (j in 1:10) {
  rsamp_method$splits[[j]]$data = select(rsamp_method$splits[[j]]$data, -ds)
}
rm(j)

save(list = ls(), file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/ensemble_env_vars/', module, '_', argv$number, '_ensemble_env_vars.rdata'))




# tune neural network hyperparameters on training data --------------------

library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors)#, lib.loc='/vortexfs1/home/dgellermcgrath/rpacks/')
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

saveRDS(nnet_tune, file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/', module, '_', argv$number, '_nnet_tunefile_grid50_compute.rds'))
rm(nnet_tune)
save(list = ls(), file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/leftover_training_files/', module, '_', argv$number, '_nnet_leftover_env_vars_grid50.rdata'))







# tune xgboost model hyperparameters --------------------------------------

#### create xgboost models on same information
library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors)#, lib.loc='/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)
library(optparse)

plan(multicore, workers = 10) #16000 * 1024 ^2 - 16Gb per core, 10 cores total
options(future.globals.maxSize = 16777216000)


#set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-f', '--file'), action = 'store',
              help = 'Path to input dataset', type = 'character'),
  make_option(c('-n', '--number'), action = 'store',
              help = 'Path to input dataset', type = 'integer')
)

#parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

#load environmental variables
load(argv$file)


xgboost_recipe <-
  recipe(vctrs::vec_slice(training_data_with_downsampled_observations_incl, 0)) %>%
  update_role(everything()) %>%
  update_role(y, new_role = 'outcome') %>%
  step_zv(all_predictors()) %>% # will remove features with zero variance
  step_select_infgain(all_predictors(), outcome = "y", top_p = tune())

xgboost_spec <- # tuning specifications for the model
  boost_tree(trees = tune(), # list parameters that will be tuned via repeated cross-validation using tune()
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
  add_model(xgboost_spec)

#changed from control_grid()
ctrl <- control_grid(verbose = TRUE, # be chatty about the fitting
                     event_level = "second",  # the '1' level, which indicates an event for this model, is the second factor level
                     parallel_over = "resamples",
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

############ check obj sizes before starting tuning
#get_size = function(.data) {
#  print(paste0(.data, ': ',
#               (get(.data) %>%
#                  object.size() %>%
#                  format(units = "Mb"))
#  ))
#}

#walk(ls(), ~ get_size(.x))
###################################################

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
saveRDS(xgboost_tune, file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/tune_files/', module, '_', argv$number, '_xgboost_tunefile_grid50_compute.rds'))
rm(xgboost_tune)
save(list = ls(), file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/leftover_training_files/', module, '_', argv$number, '_xgboost_leftover_env_vars_grid50.rdata'))



# train stacked ensemble model --------------------------------------------

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

























## testing the ensemble model
confMats = x10_test %>%
  map(~ caret::confusionMatrix(
    factor(predict(ensemble_model, .x)[[1]], levels = c('-1', '1')),
    factor(.x$y, levels = c('-1', '1')),
    positive = '1')
  )

confMats = set_names(confMats, paste0('sim', seq(100, 10, by = -10)))

perf_metrics <- map(
  list(confMats), ~
    .x %>%
    map2_dfr(., names(.), function(.x, .y) {
      .x %>%
        purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
        purrr::flatten() %>%
        as_tibble() %>%
        mutate(sim_comp = .y, .before = 1) %>%
        relocate(c(F1, Precision, Recall, `Balanced Accuracy`), .after = 2)
    })
)



# analyze individual ensemble model performances compared to full ensemble model

#member_preds = test %>%
#  select(y) %>%
#  bind_cols(predict(ensemble_model, test, members = TRUE))

#indiv_preds = member_preds %>%
#  map_dfr(~ f_meas(truth = y, data = member_preds, estimate = .x)) %>%
#  mutate(member = colnames(member_preds))

#indiv_preds = member_preds %>%
#  map_dfr(~ f_meas(truth = y, data = member_preds, estimate = .x)) %>%
#  mutate(member = colnames(member_preds))


























































# to test nnet the model -------------------------------------------------------

# M00086

library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors, lib.loc='/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)

Rcpp::sourceCpp('/vortexfs1/home/dgellermcgrath/cpp_fns_files/rarefy.cpp')

plan(multicore, workers = 10) #16000 * 1024 ^2 - 16Gb per core, 10 cores total
options(future.globals.maxSize = 16777216000)

load('nnet_env_50_ensemble_test.rdata')
nnet_tune = readRDS('nnet_grid_50_compute_ensemble_test.rds')

load('xgb_env_50_ensemble_test.rdata')



best_kap <- select_best(nnet_tune, "kap")

#kap_info <- nnet_tune %>%
#  show_best()

final_nnet <- finalize_workflow(
  nnet_workflow,
  best_kap
)

#fit the final model
model <- fit(final_nnet, training_data_with_downsampled_observations_incl)

#### test the model on simulated incomplete genomes.. 10 - 100 %


test_metadata = select(test, c(y, genome_name))
test = select(test, -c(y, genome_name))


if (all(map_chr(test, ~ class(.x)) == 'integer')) {
  test = as.matrix(test)
} else {
  stop('Error: Not all columns were integer class.')
}


rm(training_data_with_downsampled_observations_incl, rsamp_method, nnet_tune, nnet_workflow)
rm(ctrl, nnet_grid, nnet_rec, nnet_spec)
rm(xgb_grid, xgboost_recipe, xgboost_spec, xgboost_workflow)
rm(final_nnet)

#test_data_with_downsampled_observations_incl =
#  create_downsampled_train_from(test, test_metadata)

# trying this 10x, take average of results

#test_data_with_downsampled_observations_incl =

create_downsampled_seed = function(.train, .metadata) {
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

######
walk(ls(), ~ get_size(.x))


######


res =
  map(1:10, function(.x) {
    x = create_downsampled_seed(test, test_metadata)

    res = map(x, ~
                caret::confusionMatrix(
                  factor(predict(model, .x)[[1]], levels = c('-1', '1')),
                  factor(test_metadata$y, levels = c('-1', '1')),
                  positive = '1'))
    return(res)
  })




res = map(res, ~ .x %>% set_names(., paste0('prop_', seq(100, 10, by = -10))))


perf_metrics <- map(
  res, ~
    .x %>%
    map2_dfr(., names(.), function(.x, .y) {
      .x %>%
        purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
        purrr::flatten() %>%
        as_tibble() %>%
        mutate(sim_comp = .y, .before = 1) %>%
        relocate(c(F1, Precision, Recall, `Balanced Accuracy`), .after = 2)
    })
)


avg_perf = reduce(map(perf_metrics, ~ select(.x, -sim_comp)), `+`)
avg_perf = avg_perf / 10
avg_perf = as_tibble(avg_perf) %>%
  mutate(sim_comp = paste0('prop_', seq(100, 10, by = -10)), .before = 1)
avg_perf



saveRDS(nnet_tune, file = 'M00086_grid_400_nnet_tune.rda')






























































































xgb = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/svm_models_november_2021/cross_validated/xgboost_tidymodels/models/M00035_intermediate_data_27_model.rda')


xgb_res =
  map(1:10, function(.x) {
    x = create_downsampled_seed(test, test_metadata, .seed = .x)

    res = map(x, ~
                caret::confusionMatrix(
                  factor(predict(xgb, .x)[[1]], levels = c('-1', '1')),
                  factor(test_metadata$y, levels = c('-1', '1')),
                  positive = '1'))
    return(res)
  })


xgb_res = map(xgb_res, ~ .x %>% set_names(., paste0('prop_', seq(100, 10, by = -10))))


xgb_perf_metrics <- map(
  xgb_res, ~
    .x %>%
    map2_dfr(., names(.), function(.x, .y) {
      .x %>%
        purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
        purrr::flatten() %>%
        as_tibble() %>%
        mutate(sim_comp = .y, .before = 1) %>%
        relocate(c(F1, Precision, Recall, `Balanced Accuracy`), .after = 2)
    })
)


xgb_avg_perf = reduce(map(xgb_perf_metrics, ~ select(.x, -sim_comp)), `+`)
xgb_avg_perf = xgb_avg_perf / 10
xgb_avg_perf = as_tibble(xgb_avg_perf) %>%
  mutate(sim_comp = paste0('prop_', seq(100, 10, by = -10)), .before = 1)
xgb_avg_perf














































