# for ensemble models --------------------------------------------------

library(tidyverse)
library(dtplyr)
library(tidymodels)
library(recipeselectors)
library(furrr)
library(stacks)

#SETTINGS FOR BIGMEM NODE
plan(multicore, workers = 30) #81000 * 1024 ^2 - 30Gb per core, 10 cores total
options(future.globals.maxSize = 84934656000)

#files = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/SRA_FILES/feb_2022_sparsification_redo/annotations/*.tsv', intern = TRUE)

#genome_names = str_replace(files, '\\/.*\\/+(.*.tsv)', '\\1')

load('/vortexfs1/home/dgellermcgrath/FRAGMENTED_GENOMES_REQD_OBJECTS_MAY_2022.rds')

models = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/ensemble_models_lower_lasso_penalty/*', intern = TRUE)
#/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemba/ensemble_models

model_names = str_replace(models, '\\/.*\\/+(M\\d{5})_\\d.*', '\\1')


## NEXT: GO DOWN TO test_and_get_ens_conf_matrices_for AND RUN ALL OF THAT CODE + BELOW IT

# load all_kegg_modules, patt.kegg_modules, filler (0 row tibble of all columns new test data needs)
#load('/vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/allKeggModules-PattKeggModules-080321.RData')
#filler = readRDS('/vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/SRA_FILES/NEW_ATTEMPT_DECEMBER_03_21/assemblies/up-to-date-filler-feb-2022.rda')


size = function(.data) {
  print(paste0(.data, ': ',
               (get(.data) %>%
                  object.size() %>%
                  format(units = "Mb"))
  ))
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


#kofam_info = map2_dfr(kofam_data, names(kofam_data), ~
#                   tibble(filename = .y,
#                          n_proteins = nrow(.x))
#                   )

#saveRDS(kofam_info, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/SRA_FILES/feb_2022_sparsification_redo/rdata_files/full_kofam_info.rda')


test_and_get_ens_conf_matrices_for = function(formatted_test_data,
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


safely_test_and_get_ens_conf_matrices_for = safely(test_and_get_ens_conf_matrices_for)

#plan(multicore, workers = 50) #49000 * 1024 ^2 - 30Gb per core, 10 cores total
#options(future.globals.maxSize = 51380224000)

#set.seed(123)
ens_conf_matrices = future_pmap(.progress = TRUE,
  .options = furrr_options(seed = TRUE),
  list(
    list(formatted_test_data),
    model_names,
    models), ~
    safely_test_and_get_ens_conf_matrices_for(..1, ..2, ..3)
)

ens_conf_matrices = ens_conf_matrices %>%
  map(~ .x$result) %>%
  set_names(model_names) %>%
  keep(~ !is.null(.x))

ens_gen_perf_metrics <- purrr::imap_dfr(
  ens_conf_matrices, ~
    .x$conf_mat %>%
    purrr::keep(~ 'Kappa' %in% names(.x) | 'F1' %in% names(.x)) %>%
    purrr::flatten() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(model_name = .y, .before = 1) %>%
    dplyr::relocate(c(F1, Precision, Recall, Specificity, `Balanced Accuracy`), .after = 1))

ens_strat =
  ens_conf_matrices %>%
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



save(ens_strat, ens_gen_perf_metrics, ens_conf_matrices, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/UPDATED_051222_ensemble_models_SPARSIFICATION_test_results_so_far.rdata')
#























# sparsified genomes re-work ----------------------------------------------

## matching sparsified genomes from sra to ncbi/gtdb db genome names
xgbConfMats$M00001$pred_md$genome_name %>%
  str_replace('(.*)_\\d+-.*', '\\1') %>%
  unique() %>%
  paste0(collapse = '|') %>%
  clipr::write_clip()

# load results from sra explorer from using above info
sra_metadata = vroom::vroom('~/Downloads/sra_explorer_metadata.tsv')
genome_metadata = readRDS('~/Documents/lasso-regression-metapredict/lasso_metadata.rda')


#on https://www.ncbi.nlm.nih.gov/search/
xgbConfMats$M00001$pred_md$genome_name %>%
  str_replace('(.*)_\\d+-.*', '\\1') %>%
  unique() %>%
  paste0(collapse = '|') %>%
  clipr::write_clip()



genome_metadata %>%
  slice(1:100) %>%
  pull(genome_id) %>%
  paste0('|') %>%
  clipr::write_clip()


load('~/Documents/import-metagenome-databases/bac-ncbi-genome-names-refseq-genbank-02232021.RData')

load('~/Documents/import-metagenome-databases/training-genome-information-tibbles-gtdb-refseq.RData')

sparse_genome_metadata = read_csv('~/Downloads/SraRunInfo.csv')


table(refseq_cur$gbrs_paired_asm %in% genome_metadata$genome_id)
table(new_fnb$Assembly %in% genome_metadata$genome_id)

metadata_key = refseq_cur %>%
  select(gbrs_paired_asm, biosample) %>%
  rename(genome_id = gbrs_paired_asm) %>%
  bind_rows(
    new_fnb %>%
      select(Assembly, BioSample) %>%
      rename(genome_id = Assembly,
             biosample = BioSample)
  ) %>%
  filter(genome_id %in% genome_metadata$genome_id)


table(sparse_genome_metadata$BioSample %in% metadata_key$biosample)

new_sparse = sparse_genome_metadata %>%
  select(Run, BioSample, LibrarySource) %>%
  rename_with(~ str_to_lower(.x)) %>%
  left_join(metadata_key, by = 'biosample') %>%
  distinct() %>%
  group_by(genome_id) %>%
  add_tally() %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

rm(refseq_cur, new_fnb, gtdb_95.tidy_trimmed, genome_metadata)

saveRDS(new_sparse, file = 'rdata_files/model_metadata/genomes_used_forsparsification_50_metadata.rda')

#
new_sparse$run %>%
  paste0(collapse = '_|') %>%
  str_replace('(.*)', '\\1_') %>%
  clipr::write_clip()



paste0('"', new_sparse$genome_id, '"') %>%
  paste0(collapse = ',') %>%
  clipr::write_clip()


## create response variables
trains = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/imba_train_data_dec_8_2021/*', intern = TRUE)

tests = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/imba_test_data_dec_8_2021/*', intern = TRUE)


train = readRDS(trains[1])
test = readRDS(tests[1])

module = colnames(train)[1]

# df = train %>%
#   select(c(1, genome_name)) %>%
#   bind_rows(
#     test %>%
#       select(c(1, genome_name))
#   ) %>%
#   mutate({{module}} := if_else(!!sym(module) == '-1', '0', !!sym(module))) %>%
#   mutate({{module}} := factor(!!sym(module), levels = c('0', '1'))) %>%
#   relocate(genome_name) %>%
#   filter(
#     genome_name %in% c(
#       "GCA_015999605.1","GCA_001629755.1","GCA_012689525.1","GCA_002761575.1","GCA_001988115.1","GCA_002128385.1","GCA_003130605.1","GCA_000191485.1","GCA_004008595.1","GCA_002787315.1","GCA_001018725.2","GCA_016728205.1","GCA_001559135.2","GCA_002951455.1","GCA_002073535.2","GCA_001019575.2","GCA_013267575.1","GCA_001718535.1","GCA_001606315.1","GCA_001677215.1","GCA_001448985.1","GCA_009626915.1","GCA_014170055.1","GCA_014170015.1","GCA_001987875.1","GCA_001987915.1","GCA_009626355.1","GCA_009626135.1","GCA_006770325.1","GCA_004636025.1","GCA_003382795.1","GCA_014058685.1","GCA_003952485.1","GCA_001697545.1","GCA_001718975.1","GCA_005537135.1","GCA_009756605.1","GCA_001913215.1","GCA_016694755.1","GCA_015654185.1","GCA_002355155.1","GCA_003925895.1","GCA_010730195.1","GCA_010731775.1","GCA_016755975.1","GCA_904863205.1","GCA_904866495.1","GCA_003491385.1","GCA_009676865.1","GCA_000567945.1"
#     )
#   )

# verify the above works
# df2 = train %>%
#   select(c(1, genome_name)) %>%
#   bind_rows(
#     test %>%
#       select(c(1, genome_name))
#   ) %>%
#   rename(y = 1) %>%
#   mutate(y = if_else(y == '-1', '0', y)) %>% #check that  this step works...
#   mutate(y = factor(y, levels = c('0', '1'))) %>%
#   relocate(genome_name) %>%
#   filter(
#     genome_name %in% c(
#       "GCA_015999605.1","GCA_001629755.1","GCA_012689525.1","GCA_002761575.1","GCA_001988115.1",
#       "GCA_002128385.1","GCA_003130605.1","GCA_000191485.1","GCA_004008595.1","GCA_002787315.1",
#       "GCA_001018725.2","GCA_016728205.1","GCA_001559135.2","GCA_002951455.1","GCA_002073535.2",
#       "GCA_001019575.2","GCA_013267575.1","GCA_001718535.1","GCA_001606315.1","GCA_001677215.1",
#       "GCA_001448985.1","GCA_009626915.1","GCA_014170055.1","GCA_014170015.1","GCA_001987875.1",
#       "GCA_001987915.1","GCA_009626355.1","GCA_009626135.1","GCA_006770325.1","GCA_004636025.1",
#       "GCA_003382795.1","GCA_014058685.1","GCA_003952485.1","GCA_001697545.1","GCA_001718975.1",
#       "GCA_005537135.1","GCA_009756605.1","GCA_001913215.1","GCA_016694755.1","GCA_015654185.1",
#       "GCA_002355155.1","GCA_003925895.1","GCA_010730195.1","GCA_010731775.1","GCA_016755975.1",
#       "GCA_904863205.1","GCA_904866495.1","GCA_003491385.1","GCA_009676865.1","GCA_000567945.1"
#     )
#   )
#
# df3 = df %>% rename(y = 2)
#
# identical(df2, df3) #TRUE


#on bigmem
library(tidyverse)
library(furrr)

#set furrr p-processing params
plan(multicore, workers = 50) #49000 * 1024 ^2 - 70Gb per core, 35 cores total
options(future.globals.maxSize = 51380224000)

get_rv = function(trainFile, testFile) {
  train = readRDS(trainFile)
  test = readRDS(testFile)

  module = colnames(train)[1]

  df = train %>%
    select(c(1, genome_name)) %>%
    bind_rows(
      test %>%
        select(c(1, genome_name))
    ) %>%
    mutate({{module}} := if_else(!!sym(module) == '-1', '0', !!sym(module))) %>%
    mutate({{module}} := factor(!!sym(module), levels = c('0', '1'))) %>%
    relocate(genome_name) %>%
    filter(
      genome_name %in% c(
        "GCA_015999605.1","GCA_001629755.1","GCA_012689525.1","GCA_002761575.1","GCA_001988115.1",
        "GCA_002128385.1","GCA_003130605.1","GCA_000191485.1","GCA_004008595.1","GCA_002787315.1",
        "GCA_001018725.2","GCA_016728205.1","GCA_001559135.2","GCA_002951455.1","GCA_002073535.2",
        "GCA_001019575.2","GCA_013267575.1","GCA_001718535.1","GCA_001606315.1","GCA_001677215.1",
        "GCA_001448985.1","GCA_009626915.1","GCA_014170055.1","GCA_014170015.1","GCA_001987875.1",
        "GCA_001987915.1","GCA_009626355.1","GCA_009626135.1","GCA_006770325.1","GCA_004636025.1",
        "GCA_003382795.1","GCA_014058685.1","GCA_003952485.1","GCA_001697545.1","GCA_001718975.1",
        "GCA_005537135.1","GCA_009756605.1","GCA_001913215.1","GCA_016694755.1","GCA_015654185.1",
        "GCA_002355155.1","GCA_003925895.1","GCA_010730195.1","GCA_010731775.1","GCA_016755975.1",
        "GCA_904863205.1","GCA_904866495.1","GCA_003491385.1","GCA_009676865.1","GCA_000567945.1"
      )
    )
}

trains = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/imba_train_data_dec_8_2021/*', intern = TRUE)
tests = system('realpath /vortexfs1/omics/pachiadaki/dgellermcgrath/lasso/feature-tables/training-data-072821/hq_bact_ko_analysis_November_2021/imba_test_data_dec_8_2021/*', intern = TRUE)

resVars = future_map2(trains, tests, ~ get_rv(.x, .y))

x = resVars[[1]] %>%
  reduce(
    resVars[2:length(resVars)], ~
      left_join(.x, .y, by = 'genome_name'),
    .init = .)

##on poseidon
#new_sparse %>% clipr::write_clip()
#new_sparse = read_tsv('...paste the clip here...')

newSparse = read_tsv('/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_flat_files/resVars_april_2022.tsv')

x = x %>%
  left_join(
    newSparse %>%
      select(run, genome_id) %>%
      rename(genome_name = genome_id),
    by = 'genome_name'
  ) %>%
  relocate(run, .before = 1)

saveRDS(x, file = '/vortexfs1/omics/pachiadaki/dgellermcgrath/ensemble/rdata/auxiliary_r_data/res_vars.rds')



