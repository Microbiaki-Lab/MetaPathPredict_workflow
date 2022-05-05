#on poseidon, bigmem node
#in dir: /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/leftover_training_files
library(tidyverse)
library(furrr)

#SETTINGS FOR BIGMEM NODE
plan(multicore, workers = 50) #49000 * 1024 ^2 - 30Gb per core, 10 cores total
options(future.globals.maxSize = 51380224000)

# xgboost
xgb_files = system('realpath *xgb*', intern = TRUE)
xgb_model_names = str_replace(xgb_files, '\\/.*\\/+(M\\d{5})_.*', '\\1')

xgb_genomes_list = future_map(
  .progress = TRUE,
  xgb_files, function(.rdata) {
    load(.rdata)
    multi_train_metadata %>%
      pull(genome_name) %>%
      unique()
  }
) %>%
  set_names(xgb_model_names)

saveRDS(xgb_genomes_list, file = '../auxiliary_r_data/xgboost_model_training_genome_names_list.rda')


# nnet
nnet_files = system('realpath *nnet*', intern = TRUE)
nnet_model_names = str_replace(nnet_files, '\\/.*\\/+(M\\d{5})_.*', '\\1')

nnet_genomes_list = future_map(
  .progress = TRUE,
  nnet_files, function(.rdata) {
    load(.rdata)
    multi_train_metadata %>%
      pull(genome_name) %>%
      unique()
  }
) %>%
  set_names(nnet_model_names)

saveRDS(nnet_genomes_list, file = '../auxiliary_r_data/nnet_model_training_genome_names_list.rda')

# run on new R session on bigmem (could run on compute too)
library(tidyverse)
library(furrr)

#SETTINGS FOR BIGMEM NODE
plan(multicore, workers = 50) #49000 * 1024 ^2 - 30Gb per core, 10 cores total
options(future.globals.maxSize = 51380224000)


load('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/UPDATED_051222_xgboost_models_SPARSIFICATION_test_results_so_far.rdata')

ensTrainingList = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/xgboost_model_training_genome_names_list.rda')

sparse_metadata = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/genomes_used_forsparsification_50_metadata.rda')

all(names(xgb_conf_matrices) == names(ensTrainingList)) #true


get_true_test_results = function(confMat,
                                 modelName,
                                 sparseMetadata,
                                 ensTrainingList) {

  filteredTestResults = confMat$pred_md %>%
    filter(genome_name != 'DRR161213_100-ko.tsv') %>%
    mutate(
      sample_name = str_replace(
        genome_name,
        '(.*)_\\d+.*', '\\1')
    ) %>%
    left_join(
      sparseMetadata %>%
        select(run, genome_id) %>%
        rename(sample_name = run),
      by = 'sample_name'
    ) %>%
    filter(
      !is.na(genome_id),
      !(genome_id %in% ensTrainingList)
    ) %>%
    mutate(
      perc_reads_ret = str_replace(
        genome_name, '.*_(\\d+)-contigs.*', '\\1'),
      perc_reads_ret = if_else(
        str_detect(
          perc_reads_ret,
          '^101$|^104$|^107$'),
        str_replace(
          perc_reads_ret,
          '\\d{1}(\\d{1})(\\d{1})', '\\1.\\2'),
        perc_reads_ret),
      perc_reads_ret = as.numeric(perc_reads_ret)
    )

  perfMetrics = filteredTestResults %>%
    group_by(perc_reads_ret) %>%
    group_map(~
                yardstick::conf_mat(
                  data = .x,
                  truth = y,
                  estimate = preds
                ) %>%
                {bind_cols(
                  pivot_wider(rename(tidy(.),
                                     .metric = 1,
                                     .estimate = 2),
                              names_from = .metric,
                              values_from = .estimate
                  ),
                  summary(.)
                )} %>%
                rename(
                  true_neg = 1,
                  false_pos = 2,
                  false_neg = 3,
                  true_pos = 4
                ) %>%
                #             summary() %>%
                mutate(perc_reads_ret = .y[[1]]) %>%
                mutate(model_name = modelName, .before = 1)
    ) %>%
    map_dfr(~ .x)

  return(perfMetrics)
}


safely_get_true_test_results = safely(get_true_test_results)

trueTestResults = future_pmap(
  list(
    xgb_conf_matrices,
    names(xgb_conf_matrices),
    list(sparse_metadata),
    ensTrainingList),
  ~ safely_get_true_test_results(..1, ..2, ..3, ..4)
)

map_lgl(trueTestResults, ~ is.null(.x$error)) %>% table()
#TRUE
#424


trueTestResults = map_dfr(trueTestResults, ~ .x$result)

saveRDS(trueTestResults, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/xgb_ACTUAL_SPARSIFICATION_TEST_RESULTS_TIBBLE.rda')


## run locally
trueTestResults = readRDS('rdata_files/model_comparisons/xgb_ACTUAL_SPARSIFICATION_TEST_RESULTS_TIBBLE.rda')

trueTestResultsMod = trueTestResults %>%
  mutate(specificity = if_else(.metric == 'spec', .estimate, NA_real_)) %>%
  group_by(model_name, perc_reads_ret) %>%
  fill(specificity, .direction = 'updown') %>%
  ungroup() %>%
  mutate(.estimate = case_when(
    (str_detect(.metric, 'recall|precision|f_meas') & specificity == 1 & is.na(.estimate)) ~ 1,
    TRUE ~ .estimate
  ))




trueTestResults %>%
  mutate(specificity = if_else(.metric == 'spec', .estimate, NA_real_)) %>%
  group_by(model_name, perc_reads_ret) %>%
  fill(specificity, .direction = 'updown') %>%
  ungroup() %>%
  filter(str_detect(.metric, 'f_meas|precision|recall')) %>%
  filter(is.na(.estimate)) %>%
  View()


#
trueTestResults %>%
  ggplot(aes(as_factor(perc_reads_ret), .estimate)) + #, fill = `Model type`)) +
  geom_boxplot() + #outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  #geom_smooth(method = "loess", se = FALSE, aes(group = `Model type`, color = `Model type`)) +
  facet_wrap(~ as_factor(.metric), scales = 'free')

###























































###
  b = xgb_conf_matrices$M00001$pred_md %>%
    filter(genome_name != 'DRR161213_100-ko.tsv') %>%
    mutate(
      sample_name = str_replace(
        genome_name,
        '(.*)_\\d+.*', '\\1')
    ) %>%
    left_join(
      sparse_metadata %>%
        select(run, genome_id) %>%
        rename(sample_name = run),
      by = 'sample_name'
    ) %>%
    filter(
      !is.na(genome_id),
      !(genome_id %in% ensTrainingList$M00001)
    ) %>%
    mutate(
      perc_reads_ret = str_replace(
        genome_name, '.*_(\\d+)-contigs.*', '\\1'),
      perc_reads_ret = if_else(
        str_detect(
          perc_reads_ret,
          '^101$|^104$|^107$'),
        str_replace(
          perc_reads_ret,
          '\\d{1}(\\d{1})(\\d{1})', '\\1.\\2'),
        perc_reads_ret),
      perc_reads_ret = as.numeric(perc_reads_ret)
    )

  perfMetrics = b %>%
    group_by(perc_reads_ret) %>%
    group_map(~
                yardstick::conf_mat(
                  data = .x,
                  truth = y,
                  estimate = preds
                 ) %>%
                {bind_cols(
                  pivot_wider(rename(tidy(.),
                         .metric = 1,
                         .estimate = 2),
                        names_from = .metric,
                        values_from = .estimate
                         ),
                  summary(.)
                )} %>%
                rename(
                  true_neg = 1,
                  false_pos = 2,
                  false_neg = 3,
                  true_pos = 4
                ) %>%
    #             summary() %>%
                 mutate(perc_reads_ret = .y[[1]]) %>%
                 mutate(model_name = modelName, .before = 1)
     ) %>%
    map_dfr(~ .x)






