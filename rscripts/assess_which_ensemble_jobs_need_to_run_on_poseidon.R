library(tidyverse)

## find which jobs still need to be run

scav_models = 150:424

failed_models = c(150,151,155,162,174,185,227,228,230,29,38,40,58,65,232,237,253,266,267,268,269,272,275,301,302,66,286,304,307,311,313,315,316,317,319,320)

completed_ensembles = c(154,156,157,158,159,160,165,167,172,175,176,178,179,181,182,183,184,214,215,216,218,219,221,222,223,224,225,226,229,231,235,252,254,270,273,276,277,278,279,280,281,282,283,284,285,289,291,292,293,299,300,303,306,309,310,312,314,318,323,325)

mods_to_run_still = scav_models[!(scav_models %in% c(failed_models, completed_ensembles))]

mods_to_run_still %>%
  paste0(collapse = ',') %>%
  clipr::write_clip()




