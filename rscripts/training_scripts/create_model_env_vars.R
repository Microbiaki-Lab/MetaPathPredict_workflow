#load libraries
library(tidymodels)
library(stacks)
library(tidyverse)
library(recipeselectors, lib.loc = '/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)
library(optparse)

#set furrr parallel process parameters
plan(multicore, workers = 10)
options(future.globals.maxSize = 16777216000) #16000 * 1024 ^2 - 16Gb max per worker, 10 cores total

#set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-t', '--train'), action = 'store',
              help = 'Path to input dataset', type = 'character'),
  make_option(c('-p', '--test'), action = 'store',
              help = 'Path to input dataset', type = 'character'),
  make_option(c('-n', '--number'), action = 'store',
              help = 'Slurm array task ID', type = 'integer')
)

#parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

#extract module name from input file name
module = str_replace(argv$train, '(M\\d{5})-.*', '\\1')

# creating XGBoost/nnet ensemble model stack ------------------------------

#read in dataset flatfile
#df = vroom::vroom(
#  file = argv$file,
#  col_types = paste0(c('ic'), paste0(rep('i', 10137), collapse = ''))
#)

x = readRDS(argv$train)
y = readRDS(argv$test)


# format the dataset
df = x %>%
  bind_rows(y) %>%
  select(-source) %>%
  rename(y = 1) %>%
  mutate(y = if_else(y == '-1', '0', y)) %>% #check that  this step works...
  mutate(y = factor(y, levels = c('0', '1')))

#rm unneeded objects
rm(x, y)

#import rcpp fn for downsampling train and test splits later
Rcpp::sourceCpp('/vortexfs1/home/dgellermcgrath/cpp_fns_files/rarefy.cpp')

set.seed(224)

split <- initial_split(df, strata = "y", prop = 0.80)
train <- training(split)
test <- testing(split)

rm(df, split)

train_metadata = select(train, c(y, genome_name))
train = select(train, -c(y, genome_name))


if (all(map_chr(train, ~ class(.x)) == 'integer')) {
  train = as.matrix(train)
} else {
  stop('Error: Not all columns were integer class.')
}

create_downsampled_train_from = function(.train, .metadata) {
  future_map(
    .options = furrr_options(seed = TRUE),
    seq(1, 0.1, by = -0.1), ~
      rarefy(.train, sample_rate = .x) %>%
      as_tibble() %>%
      bind_cols(.metadata) %>%
      relocate(c(y, genome_name), .before = 1) %>%
      mutate(prop_counts_retained = .x, .after = genome_name)
  )
}

training_data_with_downsampled_observations_incl =
  create_downsampled_train_from(train, train_metadata)

rm(train)

training_data_with_downsampled_observations_incl =
  training_data_with_downsampled_observations_incl %>%
  map_dfr(~ .x) %>%
  slice_sample(prop = 0.10)

# the metadata columns of the training data including downsampled observations
multi_train_metadata = training_data_with_downsampled_observations_incl %>%
  select(c(y, genome_name, prop_counts_retained))

# remove some cols
training_data_with_downsampled_observations_incl =
  training_data_with_downsampled_observations_incl %>%
  select(-c(genome_name, prop_counts_retained)) # keep the y column - the dependent variable - kegg module pres/abs

# 10-fold cross-validation
#stratifying the cv folds - balancing the amt of data from each downsampling increment in each fold
rsamp_method = vfold_cv(
  training_data_with_downsampled_observations_incl %>%
    mutate(ds = multi_train_metadata$prop_counts_retained, .after = y),
  v = 10, strata = "ds")

# remove the ds column from the rsamp object folds
for (j in 1:10) {
  rsamp_method$splits[[j]]$data = select(rsamp_method$splits[[j]]$data, -ds)
}
rm(j)

save(list = ls(), file = paste0('/vortexfs1/home/dgellermcgrath/ensemble/rdata/ensemble_env_vars/', module, '_', argv$number, '_ensemble_env_vars.rdata'))

