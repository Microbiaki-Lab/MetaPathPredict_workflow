# evaluate ensemble models on held-out test data --------------------------

#run interactively on a Poseidon HPC bigmem2 node

# load required libraries
library(tidyverse)
library(tidymodels)
library(furrr)
library(stacks)
library(recipeselectors)

# set furrr parallel process parameters for a Poseidon compute node
plan(multicore, workers = 50)
options(future.globals.maxSize = 51380224000) #49000 * 1024 ^2 - 49Gb max per worker, 50 cores total


# set the working directory to the one containing all saved ensemble models
setwd('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models/')

# create vector of simulated incomplete test datasets
sim_data_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/test_datasets/*', intern = TRUE)

# create vector model_paths of all full paths to saved model files
model_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/svm_models_november_2021/cross_validated/xgboost_tidymodels/models/*.rda', intern = TRUE)

#extract KEGG module identifiers from model_paths
module_names = str_replace(model_paths, '\\/.*\\/+(M\\d{5})_.*', '\\1')

## make sure there's an ensemble model for each test dataset
sim_data_paths = sim_data_paths[str_detect(sim_data_paths, paste0(module_names, collapse = '|'))]

## making doubly sure
module_names = str_replace(sim_data_paths, '\\/.*\\/+(M\\d{5})_.*', '\\1')
model_paths = model_paths[str_detect(model_paths, paste0(module_names, collapse = '|'))]


# modified just for old models - changed prediction factor levels to 0,1 from -1,1
test_and_get_conf_matrices_for = function(.test_data_path,
                                          .module_name,
                                          .model_path) {
  # load test_data and model
  .test_data = readRDS(.test_data_path)
  .model = readRDS(.model_path)

  # add column of predictions to the test_data
  .test_data = .test_data %>%
    relocate(y, .after = genome_name) %>%
    mutate(preds = predict(.model, .) %>% pull(1), .after = y) %>%
    select(c(prop_counts_retained, genome_name, y, preds)) %>%
    mutate(preds = as.character(preds),
           preds = if_else(preds == '-1', '0', preds),
           preds = factor(preds, levels = c('0', '1')))

  # create a confusion matrix for predictions
  cm = caret::confusionMatrix(.test_data$y, .test_data$preds, positive = '1')

  return(list(pred_md = .test_data, conf_mat = cm))
}

# wrap function in a try-catch fn; fn will either run successfully and save result,
# or fail and save error; importantly it will not break a loop if one or more iterations fail
safely_test_and_get_conf_matrices_for = safely(test_and_get_conf_matrices_for)

old_xgb_testdata_conf_matrices = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    sim_data_paths,
    module_names,
    model_paths), ~
    safely_test_and_get_conf_matrices_for(..1, ..2, ..3)
)

old_xgb_testdata_conf_matrices = old_xgb_testdata_conf_matrices %>%
  map(~ .x$result) %>%
  set_names(module_names) %>%
  keep(~ !is.null(.x))

old_xgb_testdata_gen_perf_metrics <- purrr::imap_dfr(
  old_xgb_testdata_conf_matrices, ~
    .x$conf_mat %>%
    purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
    purrr::flatten() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(model_name = .y, .before = 1) %>%
    dplyr::relocate(c(F1, Precision, Recall, Specificity, `Balanced Accuracy`), .after = 1))

old_xgb_testdata_strat =
  old_xgb_testdata_conf_matrices %>%
  imap(~ .x$pred_md %>%
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
             mutate(prop_counts_retained = .name, module_name = .y) %>%
             relocate(c(prop_counts_retained, module_name,
                        F1, Precision, Recall,
                        Specificity, `Balanced Accuracy`),
                      .before = 1)
         })
  )

save(old_xgb_testdata_conf_matrices, old_xgb_testdata_gen_perf_metrics, old_xgb_testdata_strat,
     file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/old_xgb_models_perf_on_testdata.rdata')
