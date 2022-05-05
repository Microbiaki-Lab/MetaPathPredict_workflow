# generate simulated incomplete protein annotation datasets from held-out test data --------

# load required libraries
library(tidyverse)
library(furrr)
library(optparse)

# set furrr parallel process parameters for a Poseidon compute node
plan(multicore, workers = 10)
options(future.globals.maxSize = 16777216000) #16000 * 1024 ^2 - 16Gb max per worker, 10 cores total

# set command line flags to be called in slurm script
option_list <- list(
  make_option(c('-e', '--env-vars'), action = 'store',
              help = 'Path to ensemble environmental variable rdata file',
              type = 'character')
)

# parse command line arguments
argv <- parse_args(OptionParser(option_list = option_list))

# extract module identifier plus the slurm_array_task_id
module_info = str_replace(
  argv$`env-vars`, '\\/.*\\/+(M\\d{5}_\\d+)_.*', '\\1'
)

# load environmental variables rdata file
load(argv$`env-vars`)

# remove unneeded environmental variables
rm(multi_train_metadata, rsamp_method, train_metadata,
   training_data_with_downsampled_observations_incl,
   create_downsampled_train_from, rarefy)

# compile Rcpp function for fast test data downsampling routine
Rcpp::sourceCpp('/vortexfs1/home/dgellermcgrath/cpp_fns_files/rarefy.cpp')

# create a test_metadata tibble; remove response var/genome name columns from test
test_metadata = select(test, c(y, genome_name))
test = select(test, -c(y, genome_name))

# convert test tibble into a matrix for fast sparsification routine via C++ (Rcpp) fn
if (all(map_chr(test, ~ class(.x)) == 'integer')) {
  test = as.matrix(test)
} else {
  stop('Error: Not all columns were integer class.')
}

# function to create sparsified/downsampled versions of the test dataset observations
# this is to simulate incomplete genome annotations, and then to test models on this
create_downsampled_test_from = function(.test, .metadata) {
  future_map(
    .options = furrr_options(seed = TRUE),
    seq(1, 0.1, by = -0.1), ~
      rarefy(.test, sample_rate = .x) %>%
      as_tibble() %>%
      bind_cols(.metadata) %>%
      relocate(c(y, genome_name), .before = 1) %>%
      mutate(prop_counts_retained = .x, .after = genome_name)
  )
}

# create list of 10 datasets - the original full test data annotations, and
# incrementally downsampled protein annotations, in 10% increments from
# 10% of protein annotation counts randomly removed, to 90% randomly removed
test = create_downsampled_test_from(test, test_metadata)

# merge all 10 list component tibbles into one tibble
test = map_dfr(test, ~ .x)

# save compressed binary .rds file of the final test dataset tibble
saveRDS(test, file = paste0(
  '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/test_datasets/',
  module_info, '_test_data_incl_downsampled_versions.rds'
))
