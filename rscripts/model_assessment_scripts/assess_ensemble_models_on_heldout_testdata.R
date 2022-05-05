# evaluate ensemble models on held-out test data --------------------------

#run interactively on a Poseidon HPC compute node

# load required libraries
library(tidyverse)
library(tidymodels)
library(furrr)
library(stacks)
library(recipeselectors)

# set furrr parallel process parameters for a Poseidon BIGMEM node
plan(multicore, workers = 50)
options(future.globals.maxSize = 51380224000) #49000 * 1024 ^2 - 49Gb max per worker, 20 cores total


# set the working directory to the one containing all saved ensemble models
setwd('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models/')

# create vector of simulated incomplete test datasets
sim_data_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/test_datasets/*', intern = TRUE)

# create vector model_paths of all full paths to saved model files
model_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models_lower_lasso_penalty/*', intern = TRUE)

#extract KEGG module identifiers from model_paths
module_names = str_replace(model_paths, '\\/.*\\/+(M\\d{5})_.*', '\\1')

## make sure there's an ensemble model for each test dataset
sim_data_paths = sim_data_paths[str_detect(sim_data_paths, paste0(module_names, collapse = '|'))]

## making doubly sure
module_names = str_replace(sim_data_paths, '\\/.*\\/+(M\\d{5})_.*', '\\1')
model_paths = model_paths[str_detect(model_paths, paste0(module_names, collapse = '|'))]


#
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
    select(c(prop_counts_retained, genome_name, y, preds))

  # create a confusion matrix for predictions
  cm = caret::confusionMatrix(.test_data$y, .test_data$preds, positive = '1')

  return(list(pred_md = .test_data, conf_mat = cm))
}

# wrap function in a try-catch fn; fn will either run successfully and save result,
# or fail and save error; importantly it will not break a loop if one or more iterations fail
safely_test_and_get_conf_matrices_for = safely(test_and_get_conf_matrices_for)

ens_testdata_conf_matrices = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    sim_data_paths,
    module_names,
    model_paths), ~
    safely_test_and_get_conf_matrices_for(..1, ..2, ..3)
)

ens_testdata_conf_matrices = ens_testdata_conf_matrices %>%
  map(~ .x$result) %>%
  set_names(module_names) %>%
  keep(~ !is.null(.x))

ens_testdata_gen_perf_metrics <- purrr::imap_dfr(
  ens_testdata_conf_matrices, ~
    .x$conf_mat %>%
    purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
    purrr::flatten() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(model_name = .y, .before = 1) %>%
    dplyr::relocate(c(F1, Precision, Recall, Specificity, `Balanced Accuracy`), .after = 1))

ens_testdata_strat =
  ens_testdata_conf_matrices %>%
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


save(ens_testdata_strat, ens_testdata_gen_perf_metrics, ens_testdata_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/UPDATED_051222_ensemble_models_DOWNSAMPLED_PROTEIN_ANNO_TESTDATA_test_results.rdata')


#save(ens_testdata_strat, ens_testdata_gen_perf_metrics, ens_testdata_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/ensemble_models_protein_annotation_sparsification_test_results_329_with_redone_poor_perf_models.rdata')


# avg metrics of 185 ensembles; 17 produced NAs...
#reduce(
#  ens_testdata_strat %>%
#    map(~ .x %>% select(F1, Precision, Recall)) %>% #, `Pos Pred Value`, `Neg Pred Value`)) %>%
#    keep(~ all(!is.na(.x$F1))),
# `+`) %>%
#  mutate(across(everything(), ~ .x / 228))

          # F1 Precision    Recall
# 1  0.7211380 0.6254326 0.9474751
# 2  0.8458874 0.7891998 0.9482835
# 3  0.8902927 0.8501335 0.9574711
# 4  0.9165961 0.8862457 0.9653334
# 5  0.9342556 0.9106865 0.9709803
# 6  0.9466202 0.9285114 0.9749507
# 7  0.9581390 0.9452643 0.9785949
# 8  0.9666373 0.9580508 0.9816356
# 9  0.9740329 0.9700535 0.9833462
# 10 0.9805292 0.9809025 0.9850131


# xgboost perfomrance just on the 228 models so far analyzed for ensemble models above
#        F1     Precision  Recall    Pos Pred Value  Neg Pred Value
# 1  0.7480872 0.6485532 0.9256438      0.6485532      0.9623316
# 2  0.8676022 0.8251719 0.9233507      0.8251719      0.9576838
# 3  0.9057893 0.8821835 0.9347301      0.8821835      0.9615938
# 4  0.9296705 0.9139561 0.9481874      0.9139561      0.9678094
# 5  0.9465582 0.9363647 0.9581823      0.9363647      0.9730716
# 6  0.9575453 0.9505100 0.9655137      0.9505100      0.9772620
# 7  0.9669173 0.9634564 0.9708969      0.9634564      0.9806221
# 8  0.9754973 0.9754585 0.9757874      0.9754585      0.9833908
# 9  0.9816402 0.9850872 0.9783425      0.9850872      0.9851241
# 10 0.9874152 0.9943998 0.9806370      0.9943998      0.9864741


## nnet performance on the 228 models above
#       F1     Precision   Recall
# 1  0.7570249 0.6704492 0.8888201
# 2  0.8525659 0.8195047 0.8934318
# 3  0.8911519 0.8760006 0.9092156
# 4  0.9160748 0.9092841 0.9240611
# 5  0.9316498 0.9297935 0.9340481
# 6  0.9436789 0.9445514 0.9432020
# 7  0.9537256 0.9565546 0.9511792
# 8  0.9621427 0.9668233 0.9577175
# 9  0.9685506 0.9749900 0.9623490
# 10 0.9745782 0.9825823 0.9668648




#save(ens_testdata_strat, ens_testdata_gen_perf_metrics, ens_testdata_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/ensemble_models_sparsification_test_results_so_far_323.rdata')












































