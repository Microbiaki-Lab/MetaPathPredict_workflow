library(reticulate)
py_install('pandas')

pd <- import('pandas')

reg_modules = pd$read_pickle('~/Downloads/01.KEGG_Regular_Module_Information.pickle') 
bifur_modules = pd$read_pickle('~/Downloads/02.KEGG_Bifurcating_Module_Information.pickle')
struct_modules = pd$read_pickle('~/Downloads/03.KEGG_Structural_Module_Information.pickle')


all_microbeAnnotator_modules = c(reg_modules, bifur_modules, struct_modules)
current_metapredict_modules = nnet_gen_perf_metrics$model_name
all_kegg_modules = readRDS('~/Documents/construct_yl_042121/allKeggModules-052221.RDA')


table(names(all_microbeAnnotator_modules) %in% current_metapredict_modules)
table(names(all_microbeAnnotator_modules) %in% all_kegg_modules$module)


pd$
