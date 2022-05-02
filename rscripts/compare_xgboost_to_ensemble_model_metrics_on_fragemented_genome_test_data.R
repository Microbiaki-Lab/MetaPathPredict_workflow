library(tidyverse)
library(tidymodels)

load('rdata_files/model_metrics/ensemble_79_models_sparsification_test_results.rdata')


xgbGenPerfMetrics = read_rds('rdata_files/model_metrics/xgboost_overall_test_results_on_sparsified_genomes.rds')
xgbConfMats = read_rds('rdata_files/model_metrics/xgboost_conf_matrices_from_model_results_on_sparsified_test_genomes.rds')
oldEnsGenPerfMetrics = ens_gen_perf_metrics
oldEnsConfMats = ens_conf_matrices

rm(ens_gen_perf_metrics, ens_conf_matrices)

load('rdata_files/model_comparisons/ensemble_models_sparsification_test_results_so_far.rdata')
newEnsGenPerfMetrics = ens_gen_perf_metrics
newEnsConfMats = ens_conf_matrices

rm(ens_gen_perf_metrics, ens_conf_matrices)


load('rdata_files/model_comparisons/xgboost_models_sparsification_test_results_so_far.rdata')
xgbGenPerfMetrics = xgb_gen_perf_metrics
xgbConfMats = xgb_conf_matrices

rm(xgb_gen_perf_metrics, xgb_conf_matrices)

newEnsGenPerfMetrics %>%
  inner_join(
      xgbGenPerfMetrics,
      by = 'model_name',
      suffix = c('_ens', '_xgb')
      ) %>%
#  left_join(
#    filter(newEnsGenPerfMetrics, model_name %in% .$model_name),
#    by = 'model_name'
#  ) %>%
  pivot_longer(
    cols = 2:last_col(),
    values_to = 'value',
    names_to = 'metric'
  ) %>%
  mutate(model_type = case_when(
    #str_detect(metric, '_oldEns') ~ 'old_ensemble',
    str_detect(metric, '_xgb') ~ 'xgboost',
    TRUE ~ 'ensemble'),
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



# some side-by-side comparisons of specific model metrics; ensemble vs xgboost

mod_compare = newEnsGenPerfMetrics %>%
  inner_join(
    xgbGenPerfMetrics,
    by = 'model_name',
    suffix = c('_ens', '_xgb')
  )

#F1 score
table(mod_compare$F1_ens > mod_compare$F1_xgb | mod_compare$F1_ens == mod_compare$F1_xgb )

#precision
table(mod_compare$Precision_ens > mod_compare$Precision_xgb | mod_compare$Precision_ens == mod_compare$Precision_xgb )

#recall
table(mod_compare$Recall_ens > mod_compare$Recall_xgb | mod_compare$Recall_ens == mod_compare$Recall_xgb)


table(mod_compare$`Pos Pred Value_ens` > mod_compare$`Pos Pred Value_xgb` |
        mod_compare$`Pos Pred Value_ens` == mod_compare$`Pos Pred Value_xgb`)


table(mod_compare$`Neg Pred Value_ens` > mod_compare$`Neg Pred Value_xgb` |
        mod_compare$`Neg Pred Value_ens` == mod_compare$`Neg Pred Value_xgb` )

mod_compare %>%
  select(model_name, Precision_ens, Precision_xgb) %>%
  mutate(diff = Precision_xgb - Precision_ens) %>%
  arrange(desc(diff)) %>%
  View()


mod_compare %>%
  select(model_name, Recall_ens, Recall_xgb) %>%
  mutate(diff = Recall_xgb - Recall_ens) %>%
  arrange(desc(diff)) %>%
  View()




























