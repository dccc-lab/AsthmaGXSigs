library("ezknitr")
library("here")

ezspin(file = here("code", "data", "explore_metadata.R"),
       out_dir = here("reports", "diagnostics"), 
       verbose = TRUE)
