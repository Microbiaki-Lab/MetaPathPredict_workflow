# for xgboost models --------------------------------------------------

library(tidyverse)
library(dtplyr)
library(tidymodels)
library(recipeselectors)
library(furrr)
library(stacks)

#BIGMEM NODE settings
plan(multicore, workers = 50) #49000 * 1024 ^2 - 49Gb per core, 50 cores total
options(future.globals.maxSize = 51380224000)

#files = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/SRA_FILES/feb_2022_sparsification_redo/annotations/*.tsv', intern = TRUE)

#genome_names = str_replace(files, '\\/.*\\/+(.*.tsv)', '\\1')


load('/vortexfs1/home/dgellermcgrath/FRAGMENTED_GENOMES_REQD_OBJECTS_MAY_2022.rds')


models = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/butchered_nnet_models/*', intern = TRUE)
#/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemba/ensemble_models

model_names = str_replace(models, '\\/.*\\/+(M\\d{5})_.*', '\\1')


## NEXT: GO DOWN TO test_and_get_ens_conf_matrices_for AND RUN ALL OF THAT CODE + BELOW IT





# load all_kegg_modules, patt.kegg_modules, filler (0 row tibble of all columns new test data needs)
#load('/vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/allKeggModules-PattKeggModules-080321.RData')
#filler = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/SRA_FILES/NEW_ATTEMPT_DECEMBER_03_21/assemblies/up-to-date-filler-feb-2022.rda')


get_size = function(.data) {
  .data %>%
    object.size() %>%
    format(units = "Mb")
}


read_kofam = function(.data, .genome_name, cutoff = 1e-7) {
  vroom::vroom(.data, show_col_types = FALSE, delim = '\t') %>%
    rename(adaptive_threshold = 1,
           gene_name = `gene name`,
           k_number = KO,
           e_value = `E-value`) %>%
    filter(!str_detect(e_value, '---')) %>%
    select(adaptive_threshold, gene_name, k_number, e_value,
           score, thrshld, `KO definition`) %>%
    mutate(e_value = as.numeric(e_value),
           score = as.numeric(score),
           thrshld = as.numeric(thrshld)) %>%
    mutate(final_score = (score / thrshld), .after = thrshld) %>%
    #filter(`E-value` <= 1e-7 | `E-value` == 0) %>%
    filter(e_value <= cutoff | !is.na(adaptive_threshold)) %>%
    group_by(gene_name) %>%
    #filter(`E-value` == min(`E-value`)) %>%
    mutate(best_score =
             if (any(!is.na(final_score))) {
               max(na.omit(final_score)) #| is.na(final_score)
             } else {
               NA_real_
             }
    ) %>%
    filter(is.na(best_score) | best_score == final_score) %>%
    ungroup() %>%
    dplyr::mutate(genome_name = .genome_name, .before = 1) %>%
    select(c(genome_name, k_number))
}


# read genome kofam tsv gene annotation files in, one by one
kofam_data = future_map2(files, genome_names, ~ read_kofam(.x, .y))

# name the resulting list
kofam_data = kofam_data %>%
  set_names(nm = map_chr(kofam_data, ~ unique(.x$genome_name)))

# check there are no NA values in any columns
all(map_lgl(kofam_data, ~ anyNA(.x)) == FALSE)
#[1] TRUE

kofam_data = kofam_data[str_detect(names(kofam_data), 'SRR13267350_|SRR1955549_|SRR11947723_|SRR5891495_|SRR5070925_|SRR2075005_|SRR3707433_|SRR4156088_|SRR9006113_|SRR3112265_|SRR1955597_|SRR12936662_|SRR2819198_|SRR2819878_|SRR4272053_|SRR1955924_|SRR9641500_|SRR1956102_|SRR3321824_|SRR3666125_|SRR2587483_|SRR10253821_|SRR11909886_|SRR11909882_|SRR5070990_|SRR5070707_|SRR10253785_|SRR10253891_|SRR10258082_|SRR9006053_|SRR4995472_|SRR12296540_|SRR7417595_|SRR3980215_|SRR1956333_|SRR5054237_|SRR11194319_|SRR5003786_|SRR13414004_|SRR13296353_|DRR041183_|DRR140917_|DRR161213_|DRR161233_|DRR255734_|ERR1023744_|ERR1015333_|SRR7441440_|SRR10436013_|SRR1609236_')]

# detect modules only in complete versions of the input genomes
#complete_genomes = kofam_data[str_detect(names(kofam_data), '100')]


formatted_test_data = future_imap_dfr(kofam_data, ~ {
  .x %>%
    dtplyr::lazy_dt() %>%
    mutate(count = 1L) %>%
    group_by(genome_name, k_number) %>%
    summarize(count = sum(count), .groups = 'drop') %>%
    pivot_wider(names_from = k_number, values_from = count) %>%
    as_tibble() %>%
    dplyr::bind_cols(dplyr::select(filler, -c(colnames(filler)[colnames(filler) %in% colnames(.)]))) %>%
    dplyr::select(colnames(filler)) %>%
    dplyr::relocate(colnames(filler)) %>%
    mutate(genome_name = .y, .before = 1)
})


test_and_get_nnet_conf_matrices_for = function(formatted_test_data,
                                              module_name,
                                              model) {

  resVars = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/res_vars.rds')

  model = readRDS(model)

  bb = formatted_test_data %>%
    mutate(temp = str_replace(
      genome_name, '^(.*)_\\d+-contigs.*', '\\1'), .before = 1) %>%
    left_join(resVars %>%
                select(run, all_of(module_name)) %>%
                rename(y = module_name,
                       temp = run),
              by = 'temp') %>%
    relocate(y, .after = genome_name) %>%
    select(-temp) %>%
    mutate(preds = predict(model, .) %>% pull(1), .after = y) %>%
    select(c(genome_name, y, preds))

  # create a confusion matrix for predictions
  cm = caret::confusionMatrix(bb$y, bb$preds, positive = '1')

  return(list(pred_md = bb, conf_mat = cm))
}


safely_test_and_get_nnet_conf_matrices_for = safely(test_and_get_nnet_conf_matrices_for)


#set.seed(123)
nnet_conf_matrices = future_pmap(
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    list(formatted_test_data),
    model_names,
    models), ~
    safely_test_and_get_nnet_conf_matrices_for(..1, ..2, ..3)
)

nnet_conf_matrices = nnet_conf_matrices %>%
  map(~ .x$result) %>%
  set_names(model_names) %>%
  keep(~ !is.null(.x))

nnet_gen_perf_metrics <- purrr::imap_dfr(
  nnet_conf_matrices, ~
    .x$conf_mat %>%
    purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
    purrr::flatten() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(model_name = .y, .before = 1) %>%
    dplyr::relocate(c(F1, Precision, Recall, Specificity, `Balanced Accuracy`), .after = 1))

nnet_strat =
  nnet_conf_matrices %>%
  imap(~ .x$pred_md %>%
         mutate(grouping_col = str_replace(
           genome_name, '^.*_(\\d+)-contigs.*', '\\1')) %>%
         group_by(grouping_col) %>%
         group_split() %>%
         keep(~ unique(.x$grouping_col != 'DRR161213_100-ko.tsv')) %>%
         set_names(
           map_chr(., ~ paste0(
             unique(.x$grouping_col), '_', nrow(.x))
           )) %>%
         map2_dfr(names(.), function(.x, .name) {
           caret::confusionMatrix(.x$y, .x$preds, positive = '1') %>%
             keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
             flatten() %>%
             as_tibble() %>%
             mutate(prr = .name,
                    module_name = .y) %>%
             relocate(c(prr, module_name,
                        F1, Precision, Recall,
                        Specificity, `Balanced Accuracy`),
                      .before = 1) %>%
             separate(prr,
                      into = c('prr', 'n'),
                      sep = '_') %>%
             mutate(prr = case_when(
               prr == '101' ~ '0.1',
               prr == '104' ~ '0.4',
               prr == '107' ~ '0.7',
               TRUE ~ prr),
               prr = as.numeric(prr))
         })
  ) %>%
  map(~ .x %>% arrange(prr))



save(nnet_strat, nnet_gen_perf_metrics, nnet_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/UPDATED_051222_nnet_models_SPARSIFICATION_test_results_so_far')


#save(nnet_strat, nnet_gen_perf_metrics, nnet_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/xgboost_models_sparsification_test_results_so_far.rdata')
#
