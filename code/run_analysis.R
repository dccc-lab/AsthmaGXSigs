library("ezknitr")
library("here")

ezspin(file = here("code", "visualization", "figure01.R"),
       out_dir = here("reports", "publication"), 
       verbose = TRUE)

ezspin(file = here("code", "visualization", "table01.R"),
       out_dir = here("reports", "publication"), 
       verbose = TRUE)

save.image(file = here("data", "final", paste0("GSEA_Run_", rundate_gsea,".rda" )))