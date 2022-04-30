library(tidyverse)
library(tidymodels)

load('rdata_files/model_metrics/ensemble_79_models_sparsification_test_results.rdata')


xgbGenPerfMetrics = read_rds('rdata_files/model_metrics/xgboost_overall_test_results_on_sparsified_genomes.rds')
xgbConfMats = read_rds('rdata_files/model_metrics/xgboost_conf_matrices_from_model_results_on_sparsified_test_genomes.rds')
oldEnsGenPerfMetrics = ens_gen_perf_metrics
oldEnsConfMats = ens_conf_matrices

rm(ens_gen_perf_metrics, ens_conf_matrices)

load('rdata_files/model_metrics/ensemble_models_sparsification_test_results_so_far.rdata')
newEnsGenPerfMetrics = ens_gen_perf_metrics
newEnsConfMats = ens_conf_matrices

rm(ens_gen_perf_metrics, ens_conf_matrices)


oldEnsGenPerfMetrics %>%
  left_join(
      filter(xgbGenPerfMetrics, model_name %in% .$model_name),
      by = 'model_name',
      suffix = c('_oldEns', '_xgb')
      ) %>%
  left_join(
    filter(newEnsGenPerfMetrics, model_name %in% .$model_name),
    by = 'model_name'
  ) %>%
  pivot_longer(
    cols = 2:last_col(),
    values_to = 'value',
    names_to = 'metric'
  ) %>%
  mutate(model_type = case_when(
    str_detect(metric, '_oldEns') ~ 'old_ensemble',
    str_detect(metric, '_xgb') ~ 'xgboost',
    TRUE ~ 'new_ensemble'),
    metric = str_replace(
      metric, '(.*)_.*', '\\1'
  )) %>%
  ggplot(aes(model_name, value,
             shape = model_type,
             color = model_type)) +
  geom_point(
    alpha = 0.4#,
    #color = 'blue',
    #shape = 4
    ) +
  # geom_point(
  #   aes(model_name, xgb_F1),
  #   alpha = 0.4
  #   ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ metric)







































