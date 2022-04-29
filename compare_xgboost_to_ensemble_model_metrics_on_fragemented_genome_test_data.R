library(tidyverse)
library(tidymodels)

load('rdata_files/model_metrics/ensemble_79_models_sparsification_test_results.rdata')

xgbGenPerfMetrics = read_rds('rdata_files/model_metrics/xgboost_overall_test_results_on_sparsified_genomes.rds')
xgbConfMats = read_rds('rdata_files/model_metrics/xgboost_conf_matrices_from_model_results_on_sparsified_test_genomes.rds')
ensGenPerfMetrics = ens_gen_perf_metrics
ensConfMats = ens_conf_matrices

rm(ens_gen_perf_metrics, ens_conf_matrices)

ensGenPerfMetrics %>%
  left_join(
      filter(xgbGenPerfMetrics, model_name %in% .$model_name),
      by = 'model_name',
      suffix = c('_ens', '_xgb')
      ) %>%
  pivot_longer(
    cols = 2:last_col(),
    values_to = 'values',
    names_to = 'names'
  ) %>%
  mutate(model_type = if_else(
    str_detect(names, '_ens'), 'ensemble', 'xgboost'),

    ) %>%
  ggplot(aes(model_name, F1)) +
  geom_point(
    alpha = 0.4,
    color = 'blue',
    shape = 4
    ) +
  geom_point(
    aes(model_name, xgb_F1),
    alpha = 0.4
    ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
  ) #+
  #scale_color_manual(values = c('M00001' = 'seagreen'))







































