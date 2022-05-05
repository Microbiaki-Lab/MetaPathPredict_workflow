# evaluate ensemble models on held-out test data --------------------------

#run interactively on a Poseidon HPC compute node

# load required libraries
library(tidyverse)
library(tidymodels)
library(furrr)
library(stacks)
library(recipeselectors)

# set furrr parallel process parameters for a Poseidon compute node
plan(multicore, workers = 20)
options(future.globals.maxSize = 8388608000) #8000 * 1024 ^2 - 8Gb max per worker, 20 cores total


# set the working directory to the one containing all saved ensemble models
setwd('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models/')

# create vector of simulated incomplete test datasets
sim_data_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/test_datasets/*', intern = TRUE)

# create vector model_paths of all full paths to saved model files
model_paths = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/nnet_models/*', intern = TRUE)

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

nnet_testdata_conf_matrices = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    sim_data_paths,
    module_names,
    model_paths), ~
    safely_test_and_get_conf_matrices_for(..1, ..2, ..3)
)

nnet_testdata_conf_matrices = nnet_testdata_conf_matrices %>%
  map(~ .x$result) %>%
  set_names(module_names) %>%
  keep(~ !is.null(.x))

nnet_testdata_gen_perf_metrics <- purrr::imap_dfr(
  nnet_testdata_conf_matrices, ~
    .x$conf_mat %>%
    purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
    purrr::flatten() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(model_name = .y, .before = 1) %>%
    dplyr::relocate(c(F1, Precision, Recall, Specificity, `Balanced Accuracy`), .after = 1))

nnet_testdata_strat =
  nnet_testdata_conf_matrices %>%
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

# avg metrics of 185 ensembles; 17 produced NAs...
reduce(
  nnet_testdata_strat %>%
    map(~ .x %>% select(F1, Precision, Recall)) %>% #397
    keep(~ all(!is.na(.x$F1))),  #389
  `+`) %>%
  mutate(across(everything(), ~ .x / 389))


# F1 Precision    Recall
# 1  0.6930658 0.5989652 0.8726299
# 2  0.8106271 0.7660687 0.8760055
# 3  0.8595839 0.8356148 0.8930932
# 4  0.8917209 0.8780510 0.9105927
# 5  0.9114872 0.9052176 0.9209126
# 6  0.9269644 0.9244926 0.9310474
# 7  0.9403547 0.9419593 0.9403027
# 8  0.9508297 0.9545702 0.9481932
# 9  0.9600315 0.9676296 0.9531959
# 10 0.9671509 0.9770125 0.9582250




compl_ens_model_names = "M00001,M00002,M00003,M00004,M00006,M00007,M00008,M00010,M00012,M00015,M00016,M00017,M00018,M00019,M00020,M00021,M00022,M00023,M00024,M00025,M00026,M00027,M00028,M00033,M00035,M00036,M00045,M00049,M00050,M00052,M00060,M00061,M00063,M00064,M00093,M00096,M00098,M00115,M00116,M00117,M00119,M00120,M00122,M00123,M00127,M00133,M00134,M00135,M00140,M00149,M00153,M00155,M00156,M00157,M00159,M00166,M00167,M00168,M00169,M00175,M00178,M00183,M00185,M00186,M00188,M00189,M00190,M00192,M00193,M00196,M00200,M00201,M00207,M00208,M00209,M00210,M00211,M00212,M00215,M00217,M00279,M00280,M00283,M00287,M00300,M00301,M00302,M00306,M00317,M00319,M00320,M00323,M00325,M00328,M00331,M00332,M00333,M00334,M00335,M00336,M00339,M00345,M00362,M00364,M00365,M00394,M00416,M00417,M00432,M00434,M00438,M00443,M00444,M00445,M00449,M00450,M00451,M00452,M00453,M00454,M00455,M00456,M00457,M00458,M00459,M00461,M00470,M00471,M00472,M00473,M00474,M00475,M00477,M00479,M00481,M00482,M00483,M00486,M00488,M00491,M00492,M00493,M00495,M00497,M00498,M00499,M00500,M00501,M00502,M00503,M00504,M00505,M00506,M00510,M00512,M00514,M00515,M00520,M00521,M00523,M00527,M00533,M00545,M00549,M00550,M00551,M00552,M00554,M00555,M00565,M00568,M00570,M00572,M00573,M00579,M00580,M00581,M00582,M00593,M00600,M00605,M00610,M00622,M00627,M00632,M00642,M00648,M00654,M00700,M00701,M00702,M00705,M00706,M00707,M00708,M00712,M00715,M00716,M00720,M00722,M00724,M00725,M00727,M00728,M00729,M00739,M00742,M00743,M00761,M00762,M00764,M00766,M00792,M00793,M00807,M00816,M00817,M00818,M00844,M00855,M00866,M00880,M00896,M00899,M00909,M00916,M00918,M00926" %>%
  str_split(pattern = ',') %>%
  unlist()



reduce(
  nnet_testdata_strat %>%
    map(~ .x %>% select(F1, Precision, Recall)) %>%
    {function(.data) {.data[names(.data) %in% compl_ens_model_names]}}() %>%
    keep(~ all(!is.na(.x$F1))),
  `+`) %>%
  mutate(across(everything(), ~ .x / 228))














save(nnet_testdata_strat, nnet_testdata_gen_perf_metrics, nnet_testdata_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/nnet_models_sparsification_test_results_so_far.rdata')












































