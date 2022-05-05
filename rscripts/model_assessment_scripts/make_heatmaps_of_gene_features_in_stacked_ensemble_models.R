#extract features selected during training

#for xgboost member models
#x$member_fits$xgboost_tune_47_1$fit$fit$fit$feature_names



xgboost_members = x$member_fits[str_detect(names(x$member_fits), 'xgboost')]

xgboost_model_features =
  xgboost_members %>%
  map(~ .x$fit$fit$fit$feature_names)


nnet_members = x$member_fits[str_detect(names(x$member_fits), 'nnet')]

nnet_model_features =
  nnet_members %>%
  map(~ .x$fit$fit$fit$coefnames)


all_features_list = c(
  xgboost_model_features,
  nnet_model_features
)

saveRDS(all_features_list, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/all_ensemble_model_coef_names_M00001_list.rds')


all_features_list = readRDS('~/Downloads/all_ensemble_model_coef_names_M00001_list.rds')

View(all_features_list)


tibble(
  cols = all_features_list[[1]],
  values = 1L
) %>%
  pivot_wider(names_from = cols, values_from = values)



feat_tbl =
all_features_list %>%
  map(~
    tibble(
      cols = .x,
      values = 1L
    ) #%>%
      #pivot_wider(names_from = cols, values_from = values)
  ) %>%
  reduce(full_join, by = 'cols') %>%
  column_to_rownames(var = 'cols') %>%
  mutate(across(everything(), ~ case_when(is.na(.x) ~ 0L,
                                          TRUE ~ .x))) %>%
  rename_with(~ names(all_features_list))



feat_tbl_heat_M00001_plot = feat_tbl %>%
  pheatmap::pheatmap(
    color = c('lightblue', 'white'),
    clustering_method = 'ward.D',
    show_rownames = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    legend = FALSE#,

    #cluster_rows = FALSE,
    #cluster_cols = FALSE,
  )

tiff(filename = '~/Downloads/metapredict_plots/heat_of_model_coefs_M00001.tiff',
     width = 10, height = 10, units = 'in', res = 100)
feat_tbl_heat_M00001_plot
dev.off()
#






















































