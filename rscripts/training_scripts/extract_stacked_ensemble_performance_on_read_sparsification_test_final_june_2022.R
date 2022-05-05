library(tidyverse)
library(tidymodels)

load('~/Documents/MetaPredict_workflow/rdata_files/model_comparisons/UPDATED_051222_ensemble_models_SPARSIFICATION_test_results_so_far.rdata')

ensTrainingList = readRDS('~/Documents/MetaPredict_workflow/rdata_files/model_metadata/xgboost_model_training_genome_names_list.rda')

sparse_metadata = readRDS('~/Documents/MetaPredict_workflow/rdata_files/model_metadata/genomes_used_forsparsification_50_metadata.rda')

checkm_sparse_reads_res = read_tsv(
  '~/Documents/MetaPredict_workflow/flat_files/checkm_sparse_genomes.txt',
  skip = 3,
  col_names = FALSE,
  show_col_types = FALSE
  ) %>% 
  mutate(group = rep(1:(nrow(.)/2), each = 2)) %>% 
  group_by(group) %>% 
  summarize(X1 = paste0(X1, collapse = '    ')) %>% 
  select(-group) %>% 
  separate(
    X1,
    into = c('Bin','Id','num_genomes',
             'num_markers','marker sets','0','1','2','3','4','5+',
             'completeness','contamination','strain_heterogeneity'),
    sep = '\\s{2,}') %>% 
  separate(
    Id,
    into = c('marker_lineage', 'id'),
    sep = ' ') %>% 
  rename(bin = Bin) %>% 
  select(c(bin, marker_lineage, completeness, contamination, strain_heterogeneity)) %>% 
  filter(str_detect(bin, pattern = paste0(sparse_metadata$run, collapse = '|'))) %>% 
  mutate(perc_reads_ret = str_replace(bin, '.*_(\\d+)-contigs.*', '\\1'))




#
ensTrainingList = ensTrainingList[
  names(ensTrainingList) %in% names(ens_conf_matrices)
]

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

trueTestResults = pmap(
  list(
    ens_conf_matrices,
    names(ens_conf_matrices),
    list(sparse_metadata),
    ensTrainingList),
  ~ safely_get_true_test_results(..1, ..2, ..3, ..4)
)

map_lgl(trueTestResults, ~ is.null(.x$error)) %>% table()
# should all be TRUE

ensemble_sparse_reads_testResults_df = trueTestResults %>% 
  map_dfr(~ .x$result)



####
ensemble_sparse_reads_testResults_df




ensemble_sparse_reads_testResults_df = ensemble_sparse_reads_testResults_df %>%
  mutate(specificity = if_else(.metric == 'spec', .estimate, NA_real_)) %>% 
  group_by(model_name, perc_reads_ret) %>% 
  fill(specificity, .direction = 'updown') %>%
  ungroup() %>% 
  mutate(.estimate = case_when(
    (str_detect(.metric, 'recall|precision|f_meas') & specificity == 1 & is.na(.estimate)) ~ 1,
    TRUE ~ .estimate
  ))




ensemble_sparse_reads_testResults_df %>% 
  filter(str_detect(.metric, 'f_meas|precision|recall')) %>% 
  filter(is.na(.estimate)) %>% 
  View()


#
ens_sparse_plot_df = ensemble_sparse_reads_testResults_df %>%
  filter(
    str_detect(
      .metric, 'f_meas|precision|recall|spec|ppv|npv'),
    perc_reads_ret <= 10) %>%
  mutate(.estimate =
    case_when(
      is.na(.estimate) & (false_pos == 0 & false_neg == 0) ~ 1, 
      is.na(.estimate) & (false_pos == 0 | false_neg == 0) ~ 0,
      TRUE ~ .estimate
    ))
  
  
  
  
ens_sparse_plot_df



  
  
  
  
  
  
  


  ggplot(aes(as_factor(perc_reads_ret), .estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    aes(fill = 'black'),
      #fill = p_pos,
    #    shape = is_stacked),
    alpha = 0.3,
    pch = 21,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.3),
    size = 1
  ) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  facet_wrap(~ as_factor(.metric)) +
  # scale_fill_gradient2(low = '#ABD9E9',
  #                      mid = '#FEE090',
  #                      high = 'orange',
  #                      midpoint = 0.50,
  #                      breaks = seq(0, 1, by = 0.2),
  #                      limits = c(0, 1)) +
  theme(
    strip.background = element_rect(fill = 'wheat1'),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 19),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)) +
  labs(
    x = '\nProportion of protein annotation counts retained',
    title = 'MetaPredict model performances on held-out test data (number of models = 415)'
#     fill = 'Proportion of genomes
# that contain the module'
  )



























