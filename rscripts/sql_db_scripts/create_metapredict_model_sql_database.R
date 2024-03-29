# save handful of ensemble models,  on Poseidon ---------------------------

setwd('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models')
models = list.files()
model_names = str_replace(models, '(M\\d{5})_.*', '\\1')

for (x in 1:10) {
  assign(model_names[x], readRDS(models[x]))
}
rm(x)

save(list = ls(), file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/test_10_ens_model_db.rdata')




# testing creating sql database full of metapredict models, locally --------
library(tidyverse)
library(RSQL)
library(RSQLite)

load('rdata_files/test_db/test_10_ens_model_db.rdata')
rm(M00005, M00006, M00007, M00008)
rm(M00010, M00012, M00015)

# create char vector containing names of all model objects in global env
model_names = ls()

# create empty SQLite database called "model_sqlite.db"
model_db <- dbConnect(drv = RSQLite::SQLite(),
                 dbname = "model_sqlite.db")

# create a table in SQL database: col 1 = model name, col 2 = raw model
  # raw model is basically the whole model saved as a binary "blob"
dbSendQuery(conn = model_db,
            "CREATE TABLE models (model_name CHARACTER, raw_model BLOB)")

# iteratively input each model name and serialized model object into the
  #SQL table
for (x in model_names) {
  RSQLite::dbExecute(model_db,
                     'INSERT INTO models VALUES (:model_name, :raw_model)',
                     params = list(
                       model_name = x,
                       raw_model = list(serialize(get(x), NULL))
                     ))
}

# create a reference obj to the database using dbplyr
model_db_ref = tbl(model_db, 'models')

# function to properly index and then unserialize a raw model "blob"
unserialize_model = function(.data) {
  unserialize(.data[[1]])
}

# code to extract a specified raw model from the SQL database and
  # then unserialize it/load it into RAM
modelx = model_db_ref %>%
  filter(model_name == 'M00001') %>%
  collect() %>%
  pull(raw_model) %>%
  unserialize_model()





# connect to a pre-made sql database,  load a model -----------------------
library(tidyverse)
#library(RSQL)
#library(RSQLite)

size = function(.object) {
  object.size(.object) %>%
    format(units = "Mb")
}


#connect to a pre-made SQL database
model_db <- RSQLite::dbConnect(drv = RSQLite::SQLite(),
                      dbname = "model_sqlite.db")

# create a reference obj to the database using dbplyr
model_db_ref = tbl(model_db, 'models')

# function to properly index and then unserialize a raw model "blob"
unserialize_model = function(.data) {
  unserialize(.data[[1]])
}

# code to extract a specified raw model from the SQL database and
# then unserialize it/load it into RAM
modelx = model_db_ref %>%
  filter(model_name == 'M00001') %>%
  collect() %>%
  pull(raw_model) %>%
  unserialize_model()






# experimenting with remove redunant ensemble model info ------------------
library(tidyverse)
library(tidymodels)
library(stacks)
library(recipeselectors)

mod_analysis = butcher::weigh(modelx)

map(modelx$member_fits$xgboost_tune_09_1$pre$actions$recipe$recipe, ~
      object.size(.x) %>%
      format(units = "Mb"))

#fn to remove xgboost redundant information
shrink_xgb_recipe = function(.data) {
  #names of features + response var included in final model
  varsToRetain = names(.data$pre$mold$blueprint$recipe$levels)

  #remove bloated xgboost recipe sections of ensemble models
  .data$pre$actions$recipe$recipe$var_info =
    .data$pre$actions$recipe$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$term_info =
    .data$pre$actions$recipe$recipe$term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$template =
    .data$pre$actions$recipe$recipe$template %>%
    select(all_of(varsToRetain))

  .data$pre$mold$blueprint$recipe$var_info =
    .data$pre$mold$blueprint$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$last_term_info =
    .data$pre$mold$blueprint$recipe$last_term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$orig_lvls =
    .data$pre$mold$blueprint$recipe$orig_lvls[varsToRetain]

  .data$pre$mold$blueprint$ptypes$predictors =
    .data$pre$mold$blueprint$ptypes$predictors %>%
    select(varsToRetain[2:length(varsToRetain)])

  return(.data)
}

#fn to remove neural network redundant information
shrink_nnet_recipe = function(.data) {
  #names of features + response var included in final model
  varsToRetain = .data$pre$mold$blueprint$recipe$term_info$variable

  .data$pre$actions$recipe$recipe$var_info =
    .data$pre$actions$recipe$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$term_info =
    .data$pre$actions$recipe$recipe$term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$template =
    .data$pre$actions$recipe$recipe$template %>%
    select(all_of(varsToRetain))

  .data$pre$mold$blueprint$recipe$var_info =
    .data$pre$mold$blueprint$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$last_term_info =
    .data$pre$mold$blueprint$recipe$last_term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$orig_lvls =
    .data$pre$mold$blueprint$recipe$orig_lvls[varsToRetain]

  .data$pre$mold$blueprint$ptypes$predictors =
    .data$pre$mold$blueprint$ptypes$predictors %>%
    select(varsToRetain[2:length(varsToRetain)])

  return(.data)
}

shrink_modeldefs = function(.data, varsToRetain) {

  ## xgboost
  .data$model_defs$xgboost_tune$pre$actions$recipe$recipe$var_info = NULL
    #.data$model_defs$xgboost_tune$pre$actions$recipe$recipe$var_info %>%
    #filter(variable %in% varsToRetain)

  .data$model_defs$xgboost_tune$pre$actions$recipe$recipe$term_info = NULL
    #.data$model_defs$xgboost_tune$pre$actions$recipe$recipe$term_info %>%
    #filter(variable %in% varsToRetain)

  .data$model_defs$xgboost_tune$pre$actions$recipe$recipe$template = NULL
    #.data$model_defs$xgboost_tune$pre$actions$recipe$recipe$template %>%
    #select(varsToRetain)


  ## nnet
  .data$model_defs$nnet_tune$pre$actions$recipe$recipe$var_info = NULL
    #.data$model_defs$nnet_tune$pre$actions$recipe$recipe$var_info %>%
    #filter(variable %in% varsToRetain)

  .data$model_defs$nnet_tune$pre$actions$recipe$recipe$term_info = NULL
    #.data$model_defs$nnet_tune$pre$actions$recipe$recipe$term_info %>%
    #filter(variable %in% varsToRetain)

  .data$model_defs$nnet_tune$pre$actions$recipe$recipe$template = NULL
    #.data$model_defs$nnet_tune$pre$actions$recipe$recipe$template %>%
    #select(varsToRetain)

  return(.data)
}

b = modelx

b$member_fits[str_detect(names(b$member_fits), 'xgboost')] =
  b$member_fits[str_detect(names(b$member_fits), 'xgboost')] %>%
  map(~ shrink_xgb_recipe(.x))

b$member_fits[str_detect(names(b$member_fits), 'nnet')] =
  b$member_fits[str_detect(names(b$member_fits), 'nnet')] %>%
  map(~ shrink_nnet_recipe(.x))

b = shrink_modeldefs(b)

#
#b = b[!(names(b) %in% c('model_metrics'))]




#see how much this shrinks the size of the ensemble model
size(b) # 96.8 Mb
size(modelx) # 352.3 Mb
# the model with redundant info removed is less than 1/3 the size

# does it still predict properly?
test_data = readRDS('rdata_files/simulated_test_data/M00001_intermediate_data_1_sim_incomp_data.rda')

pred_full = predict(modelx, test_data$prop.10)
pred_reduced = predict(b, test_data$prop.10)

identical(pred_full, pred_reduced) # [1] TRUE

# the predictions appear to be identical.
# note: need to run more thorough test to confirm



# create official MetaPredict SQL database --------------------------------

# open up an interactive node on Poseidon

#load reqd libraries
library(tidyverse)
library(furrr)
library(RSQL)
library(RSQLite)

#set parallel processing parameters
plan(multicore, workers = 35) #2700 * 1024 ^2 - 2.7Gb per core, 35 cores total
options(future.globals.maxSize = 2831155200) #set maximum memory per parallal worker

#
models = list.files()
model_names = str_replace(models, '(M\\d{5})_.*', '\\1')


# fn to remove xgboost redundant information
shrink_xgb_recipe = function(.data) {
  #names of features + response var included in final model
  varsToRetain = names(.data$pre$mold$blueprint$recipe$levels)

  #remove bloated xgboost recipe sections of ensemble models
  .data$pre$actions$recipe$recipe$var_info =
    .data$pre$actions$recipe$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$term_info =
    .data$pre$actions$recipe$recipe$term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$template =
    .data$pre$actions$recipe$recipe$template %>%
    select(all_of(varsToRetain))

  .data$pre$mold$blueprint$recipe$var_info =
    .data$pre$mold$blueprint$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$last_term_info =
    .data$pre$mold$blueprint$recipe$last_term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$orig_lvls =
    .data$pre$mold$blueprint$recipe$orig_lvls[varsToRetain]

  .data$pre$mold$blueprint$ptypes$predictors =
    .data$pre$mold$blueprint$ptypes$predictors %>%
    select(varsToRetain[2:length(varsToRetain)])

  return(.data)
}


# fn to remove neural network redundant information
shrink_nnet_recipe = function(.data) {
  #names of features + response var included in final model
  varsToRetain = .data$pre$mold$blueprint$recipe$term_info$variable

  .data$pre$actions$recipe$recipe$var_info =
    .data$pre$actions$recipe$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$term_info =
    .data$pre$actions$recipe$recipe$term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$actions$recipe$recipe$template =
    .data$pre$actions$recipe$recipe$template %>%
    select(all_of(varsToRetain))

  .data$pre$mold$blueprint$recipe$var_info =
    .data$pre$mold$blueprint$recipe$var_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$last_term_info =
    .data$pre$mold$blueprint$recipe$last_term_info %>%
    filter(variable %in% varsToRetain)

  .data$pre$mold$blueprint$recipe$orig_lvls =
    .data$pre$mold$blueprint$recipe$orig_lvls[varsToRetain]

  .data$pre$mold$blueprint$ptypes$predictors =
    .data$pre$mold$blueprint$ptypes$predictors %>%
    select(varsToRetain[2:length(varsToRetain)])

  return(.data)
}


# remove unneeded model definitions from ensemble model objects
shrink_modeldefs = function(.data) {

  ## xgboost
  .data$model_defs$xgboost_tune$pre$actions$recipe$recipe$var_info = NULL
  .data$model_defs$xgboost_tune$pre$actions$recipe$recipe$term_info = NULL
  .data$model_defs$xgboost_tune$pre$actions$recipe$recipe$template = NULL

  ## nnet
  .data$model_defs$nnet_tune$pre$actions$recipe$recipe$var_info = NULL
  .data$model_defs$nnet_tune$pre$actions$recipe$recipe$term_info = NULL
  .data$model_defs$nnet_tune$pre$actions$recipe$recipe$template = NULL

  return(.data)
}


#fn to load a model, and remove info not reqd for predictions on new data
load_and_trim = function(model) {

  mod = readRDS(model)

  mod$member_fits[str_detect(names(mod$member_fits), 'xgboost')] =
    mod$member_fits[str_detect(names(mod$member_fits), 'xgboost')] %>%
    map(~ shrink_xgb_recipe(.x))

  mod$member_fits[str_detect(names(mod$member_fits), 'nnet')] =
    mod$member_fits[str_detect(names(mod$member_fits), 'nnet')] %>%
    map(~ shrink_nnet_recipe(.x))

  mod = shrink_modeldefs(mod)

  #assign(model_name, mod, envir = .GlobalEnv)
  #rm(mod)

  #return(get(model_name))
  #return(assign(model_name, mod))

  return(mod)
}

# load and trim down models in parallel, 35 at a time
# use 'walk' instead of 'map' bc we are just loading & trimming model objects
#future_walk2(models, model_names, ~ assign(.y, load_and_trim(.x), envir = .GlobalEnv))
#the above command doesn't load anything; not sure how to assign from fn

for (x in seq_along(models)) {
  assign(model_names[x], load_and_trim(models[x]))
}


# create empty SQLite database called "model_sqlite.db"
model_db <- dbConnect(drv = RSQLite::SQLite(),
                      dbname = "/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/sql_db/model_sqlite.db")

# create a table in SQL database: col 1 = model name, col 2 = raw model
# raw model is basically the whole model saved as a binary "blob"
dbSendQuery(conn = model_db,
            "CREATE TABLE models (model_name CHARACTER, raw_model BLOB)")

# iteratively input each model name and serialized model object into the
#SQL table
for (x in model_names) {
  RSQLite::dbExecute(model_db,
                     'INSERT INTO models VALUES (:model_name, :raw_model)',
                     params = list(
                       model_name = x,
                       raw_model = list(serialize(get(x), NULL))
                     ))
}

# create a reference obj to the database using dbplyr
model_db_ref = tbl(model_db, 'models')





# restart R session; then connect to SQL db and try loading a model -------

library(tidyverse)
library(tidymodels)
library(stacks)
library(recipeselectors)


#check object size
size = function(.data) {
  .data %>%
    object.size() %>%
    format(units = 'Mb')
}

# function to properly index and then unserialize a raw model "blob"
unserialize_model = function(.data) {
  unserialize(.data[[1]])
}

#connect to the SQL database
model_db <- RSQLite::dbConnect(drv = RSQLite::SQLite(),
                      dbname = "/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/sql_db/model_sqlite.db")

#set up a reference to the database
model_db_ref = tbl(model_db, 'models')

# code to extract a specified raw model from the SQL database and
# then unserialize it/load it into RAM
modelx = model_db_ref %>%
  filter(model_name == 'M00002') %>%
  collect() %>%
  pull(raw_model) %>%
  unserialize_model()

# works












































