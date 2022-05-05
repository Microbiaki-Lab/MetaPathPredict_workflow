library(tidyverse)

poor_perf_modnames = model_comparisons_data %>%
  filter(
    model_type == 'ensemble',
    metric == 'Precision',
    is.na(value) | (prop_counts_retained < 0.5 & value < 0.3)
  ) %>% #View()
  pull(module_name) %>%
  unique()

poor_perf_modnames %>%
  paste0('*', .,'*_') %>%
  paste0(collapse = '|') %>%
  clipr::write_clip()

length(poor_perf_modnames)

## on command line:
## in poseidon dir: /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models

# for FILE in $(ls | egrep "*M00035*_|*M00036*_|*M00087*_|*M00120*_|*M00134*_|*M00167*_|*M00186*_|*M00199*_|*M00204*_|*M00205*_|*M00206*_|*M00218*_|*M00220*_|*M00221*_|*M00224*_|*M00253*_|*M00265*_|*M00302*_|*M00319*_|*M00325*_|*M00326*_|*M00328*_|*M00339*_|*M00342*_|*M00429*_|*M00436*_|*M00438*_|*M00460*_|*M00468*_|*M00478*_|*M00482*_|*M00483*_|*M00490*_|*M00506*_|*M00510*_|*M00518*_|*M00522*_|*M00525*_|*M00535*_|*M00577*_|*M00583*_|*M00592*_|*M00593*_|*M00595*_|*M00622*_|*M00649*_|*M00655*_|*M00656*_|*M00699*_|*M00708*_|*M00717*_|*M00720*_|*M00733*_|*M00747*_|*M00766*_|*M00792*_|*M00818*_|*M00851*_|*M00866*_|*M00896*_|*M00912*_|*M00918*_|*M00925*_|*M00930*_"); do cp $FILE ../mods_to_dl; done

# scp dgellermcgrath@poseidon.whoi.edu:/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/mods_to_dl/* ~/Documents/MetaPredict_workflow/rdata_files/poor_performance_ensembles/






#### looking into ensemble models that failed to stack

# vector of all model ensemble array job batch names
ens_model_array_names = read_lines('flat_files/ensemble_model_array_names.txt') %>%
  str_split(pattern = ',') %>%
  unlist()

# vector of those ensembles that successfully stacked
ens_arrays_successful = read_lines('flat_files/ensemble_model_arrays_successfully_fitted_as_of_may_9_2022.txt') %>%
  str_split(pattern = ',') %>%
  unlist()

# vector of unsuccessfully stacked models
ens_arrays_unsuccessful = ens_model_array_names[!(ens_model_array_names %in% ens_arrays_successful)]

# negative:positive class ratios of unsuccessful model stacks
view(class_ratios %>%
  filter(module_name %in%
           (ens_arrays_unsuccessful %>%
              str_replace('_\\d+', '')
           ))

)
### see ratios for successful models, looking for those with the highest imbalance to see what kinds of
# ensembles with imbalanced data were successfully stacked

class_ratios %>%
  filter(!(module_name %in%
           (ens_arrays_unsuccessful %>%
              str_replace('_\\d+', '')
           ))) %>%
  view()

# they all have extremely imbalanced negative:positive class ratios


# but the issue may be that when boostrapping prediction results
# of the various ensemble model candidates, not enough minority class
# prediction outcomes are making it one or more of the bootstraps
# which causes glmnet to error out



# possible solutions:

# 1. we will try stacking M00038 with looser lasso penalty values - this fixed some weird models with few ensemble members that just predicted '0' for all observations. using M00038 as a sample model that didnt stack. if successful, will apply whatever changes worked for M00038 to all the other models that failed to stack

# 2. we will creating tune files used in stacking with more observations in them

# 3. we will look at the stacks code to see if it can be manually adjusted to work

# 4. we will use whichever NN/XGBoost model did the best on the remaining unstackable
    # models for now, until a different solution is found

library(tidyverse)
library(tidymodels)
library(stacks)
library(recipeselectors)

setwd('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models')

# load req'd tune files
xgboost_tune = readRDS('../tune_files/M00038_29_xgboost_tunefile_grid50_compute.rda')
nnet_tune = readRDS('../tune_files/M00038_29_nnet_tunefile_grid50_compute.rda')

data_stack <-
  stacks() %>%
  add_candidates(xgboost_tune) %>%
  add_candidates(nnet_tune)

rm(xgboost_tune, nnet_tune)

ensemble_model =
  data_stack %>%
  blend_predictions(
    control = control_grid(allow_par = TRUE),
    penalty = 10^(-20:-0.01) #10^(-6:-3)
  )


## all boostraps continue to fail even with a broad range of lasso penalty values

rm(data_stack)

# register a parallel backend to run blend_predictions and fit_ensemble w/ parallel processing
options(cores = 30)
registerDoParallel()

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




































