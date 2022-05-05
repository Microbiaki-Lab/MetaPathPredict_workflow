library(tidyverse)
library(tidymodels)
library(stacks)
library(cowplot)

setwd('~/Documents/MetaPredict_workflow/rdata_files/poor_performance_ensembles')

model_file_paths = system('find ~/Documents/MetaPredict_workflow/rdata_files/poor_performance_ensembles -name "*" -type f', intern = TRUE)
model_names = str_replace(model_file_paths, '\\/.*\\/+(M\\d{5}).*', '\\1')

#fn to get size of an environmental object
size = function(.obj) {
  .obj %>% 
    object.size() %>% 
    format(units = 'Mb')
}

# read in the first 20 models to the global env
for (x in 1:20) {
  assign(model_names[x], readRDS(model_file_paths[x]))
}

# get size of the models in Mb;
# commented out section calculates total size in memory (Mb) of all 20 models
map_chr(model_names[1:20], ~ size(get(.x))) %>%
  set_names(model_names[1:20]) #%>% 
  #str_replace(' Mb', '') %>% 
  #as.numeric() %>% 
  #reduce(`+`)

# see how many ensemble members are in each of the first 20 "poor performing" models
map_int(model_names[1:20], ~ length(get(.x)$member_fits))%>%
  set_names(model_names[1:20])
        
# overwrite the model env objects with their autoplots showing useful info about each model ensemble
for (model in model_names[1:20]) {
  assign(
    model,
    autoplot(get(model)) +
      labs(title = model)
  )
}


# write model plot names to clipboard
ls() %>%
  {function(.data) {.data[str_detect(.data, 'M\\d{5}$')]}}() %>% 
  paste0(collapse = ', ') %>% 
  clipr::write_clip()


#plot them together
plot_grid(M00120, M00206, M00221, M00265, M00326, M00342, M00436, M00478, M00490, M00518, M00522, M00525, M00592, M00656, M00699, M00708, M00818, M00896, M00918, M00930)


# read in each model, then overwrite the written object with a stacks autoplot 
# that gives important visual information about the model's ensemble fit
for (x in seq_along(model_file_paths)) {
  assign(
    model_names[x], 
    readRDS(model_file_paths[x])
  )
  
  assign(
    model_names[x],
    autoplot(get(model_names[x])) +
      labs(title = model_names[x])
  )
}

#1-20
model_names[1:20] %>% 
  paste0(collapse = ', ') %>% 
  clipr::write_clip()

#plot them together
plot_grid(M00918, M00930, M00326, M00699, M00522, M00478, M00221, M00592, M00818, M00708, M00120, M00436, M00518, M00490, M00265, M00896, M00206, M00342, M00525, M00656)



#21-40
model_names[21:40] %>% 
  paste0(collapse = ', ') %>% 
  clipr::write_clip()

#plot them together
plot_grid(M00510, M00851, M00925, M00912, M00535, M00319, M00460, M00205, M00218, M00595, M00429, M00167, M00717, M00134, M00792, M00325, M00747, M00720, M00593, M00482)



#41-60
model_names[41:60] %>% 
  paste0(collapse = ', ') %>% 
  clipr::write_clip()

#plot them together
plot_grid(M00220, M00468, M00577, M00302, M00506, M00035, M00583, M00328, M00036, M00483, M00438, M00087, M00199, M00622, M00733, M00253, M00186, M00204, M00339, M00766)



#61-64
model_names[61:64] %>% 
  paste0(collapse = ', ') %>% 
  clipr::write_clip()

#plot them together
plot_grid(M00866, M00649, M00655, M00224)


# based on above plots: we can try to re-fit these 64 models, setting the lasso penalty upper bound to 1e-03; 1e-01 seems to create 
# problematic ensembles with few members that have horrible performance; while lowering the penalty and having more members
# creates good models (based on one detailed look into M00930)









































