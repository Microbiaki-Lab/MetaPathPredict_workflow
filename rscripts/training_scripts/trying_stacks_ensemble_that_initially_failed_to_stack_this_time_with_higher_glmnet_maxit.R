#withr::with_libpaths(new = '/vortexfs1/home/dgellermcgrath/manually_dl_R_libs',
#                     devtools::install_github("tidymodels/stacks", ref = "main",
#                                              force = TRUE))



custom_blend_predictions <- function(data_stack,
                                     penalty = 10 ^ (-6:-1),
                                     mixture = 1,
                                     non_negative = TRUE,
                                     metric = NULL,
                                     control = tune::control_grid(),
                                     ...) {
  check_inherits(data_stack, "data_stack")
  check_blend_data_stack(data_stack)
  check_regularization(penalty, "penalty")
  check_regularization(mixture, "mixture")
  check_inherits(non_negative, "logical")
  if (!is.null(metric)) {
    check_inherits(metric, "metric_set")
  }
  check_inherits(control, "control_grid")
  check_empty_ellipses(...)

  outcome <- attr(data_stack, "outcome")

  preds_formula <-
    rlang::new_formula(as.name(outcome), as.name("."), env = rlang::base_env()) # creates "label ~ ." formula

  lvls <- levels(data_stack[[outcome]]) # the levels of the outcome variable, 'y'; either 0 or 1

  dat <- process_data_stack(data_stack) # creates a tibble with .pred_0_* cols and y

  ll <- if (non_negative) {0} else {-Inf}

  tune_quo <- rlang::new_quosure(tune::tune(), env = rlang::empty_env())

  if (attr(data_stack, "mode") == "regression") {
    model_spec <-
      parsnip::linear_reg(penalty = !!tune_quo, mixture = !!tune_quo) %>%
      parsnip::set_engine("glmnet", lower.limits = !!ll, lambda.min.ratio = 0, maxit = 10^6)

    preds_wf <-
      workflows::workflow() %>%
      workflows::add_model(model_spec) %>%
      workflows::add_formula(preds_formula)
  } else {
    # The class probabilities add up to one so we remove the probability columns
    # associated with the first level of the outcome.
    col_filter <- paste0(".pred_", lvls[1])
    dat <- dat %>% dplyr::select(-dplyr::starts_with(!!col_filter))
    if (length(lvls) == 2) {
      model_spec <-
        parsnip::logistic_reg(penalty = !!tune_quo, mixture = !!tune_quo) %>%
        parsnip::set_engine("glmnet", lower.limits = !!ll, lambda.min.ratio = 0, maxit = 10^6) %>%
        parsnip::set_mode("classification")
    } else {
      model_spec <-
        parsnip::multinom_reg(penalty = !!tune_quo, mixture = !!tune_quo) %>%
        parsnip::set_engine("glmnet", lower.limits = !!ll, lambda.min.ratio = 0, maxit = 10^6) %>%
        parsnip::set_mode("classification")
    }

    preds_wf <-
      workflows::workflow() %>%
      workflows::add_recipe(
        recipes::recipe(
          preds_formula,
          data = dat
        )
      ) %>%
      workflows::add_model(model_spec)
  }

  get_models <- function(x) {
    x %>%
      workflows::extract_fit_parsnip() %>%
      purrr::pluck("fit")
  }

  control$extract <- get_models

  candidates <-
    preds_wf %>%
    tune::tune_grid(
      resamples = rsample::bootstraps(dat), #, strata = y),
      grid = purrr::cross_df(
        list(
          penalty = penalty,
          mixture = mixture
        )
      ),
      metrics = metric,
      control = control
    )

  metric <- tune::.get_tune_metric_names(candidates)[1]
  best_param <- tune::select_best(candidates, metric = metric)
  coefs <-
    model_spec %>%
    tune::finalize_model(best_param) %>%
    parsnip::fit(formula = preds_formula, data = dat)

  model_stack <-
    structure(
      list(model_defs = attr(data_stack, "model_defs"),
           coefs = coefs,
           penalty = list(
             penalty = best_param$penalty,
             mixture = best_param$mixture,
             metric = metric
           ),
           metrics = glmnet_metrics(candidates),
           equations = get_expressions(coefs),
           cols_map = attr(data_stack, "cols_map"),
           model_metrics = attr(data_stack, "model_metrics"),
           train = attr(data_stack, "train"),
           mode = attr(data_stack, "mode"),
           outcome = attr(data_stack, "outcome"),
           data_stack = dat,
           splits = attr(data_stack, "splits")),
      class = c("linear_stack", "model_stack", "list")
    )

  if (model_stack_constr(model_stack)) {model_stack}
}


environment(custom_blend_predictions) <- asNamespace('stacks')
assignInNamespace("blend_predictions", custom_blend_predictions, ns = "stacks")


library(tidymodels)
library(stacks) #, lib.loc = '/vortexfs1/home/dgellermcgrath/rpacks/')
library(tidyverse)
library(recipeselectors)#, lib.loc='/vortexfs1/home/dgellermcgrath/rpacks/')
library(furrr)
library(butcher)
library(doParallel)


# load req'd tune files
xgboost_tune = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/M00038_29_xgboost_tunefile_grid50_compute.rda')
nnet_tune = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/tune_files/M00038_29_nnet_tunefile_grid50_compute.rda')

data_stack <-
  stacks() %>%
  add_candidates(xgboost_tune) %>%
  add_candidates(nnet_tune)

rm(xgboost_tune, nnet_tune)

## worked running it like this!!! perhaps allow_par is defaulted as TRUE and was causing the node to run out of memory...
#get_ensemble_model = function(data_st) {
  data_st %>%
    blend_predictions(control = control_grid(allow_par = TRUE))
}

#safely_get_ensemble_model = safely(get_ensemble_model)

#rm_rows = sample(sample(which(data_stack$y == 0)), 19000)

#data_stack2 = data_stack %>%
#  slice(-rm_rows)

#ensemble_model = data_stack2 %>%
#  blend_predictions(control = control_grid(allow_par = TRUE))



# preds_wf <-
#   workflows::workflow() %>%
#   workflows::add_recipe(
#     recipes::recipe(
#       y ~ .,
#       data = b %>% select(-starts_with('.pred_0'))
#     )
#   ) #%>%
# #  workflows::add_model(model_spec)
#
#
# get_models <- function(x) {
#   x %>%
#     workflows::extract_fit_parsnip() %>%
#     purrr::pluck("fit")
# }
#
# control$extract <- get_models
#
# candidates <-
#   preds_wf %>%
#   tune::tune_grid(
#     resamples = rsample::bootstraps(b %>% select(-starts_with('.pred_0')), strata = 'y'),
#     grid = 5,
#     metrics = 'roc_auc')




# rm_rows = sample(sample(which(data_stack$y == 0)), 15000)
#
# data_stack2 = data_stack %>%
#   slice(-rm_rows)

### data stack from xgboost tuned on full M00038 dataset on bigmem- including simulated observations
data_stack = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/data_stack_M00038.rds')
###

ensemble_model = data_stack %>%
  blend_predictions(
    penalty = 10^(-6:-3),
    control = control_grid(allow_par = TRUE)
)
   # control = control_grid(allow_par = TRUE),
    #test_strata = y)

options(cores = 5)
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

saveRDS(ensemble_model, file = paste0('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models/', module, '_ensemble_model.rds'))

