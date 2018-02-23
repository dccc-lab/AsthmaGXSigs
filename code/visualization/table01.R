## ---- setup_env
library("here")

if(!exists("es_ubiopred_wb")) source(here("code", "data", "create_ubiopred_eset.R"), echo = TRUE)

library("tableone")

output_dir <- "reports/publication"


## ---- make_table01
cat_pheno_list <- c("sex", "race")
all_pheno_list <- c(cat_pheno_list, "RIN", "WBC", "PCTEOS", "PCTLYMPH", "PCTMONO", "PCTNEUT")
ubiopred_phenos <- es_ubiopred_wb %>% pData %>% select("cohort", all_pheno_list) %>% mutate(sex = fct_rev(sex))

table01 <- CreateTableOne(data = ubiopred_phenos, 
                          vars = all_pheno_list, 
                          factorVars = cat_pheno_list,
                          strata = "cohort")

table01_summary <- summary(table01)
table01_print <- print(table01, nonnormal = TRUE)

write.csv(table01_print, file.path(here(), output_dir, "table01.txt"))


## ---- clean_up_env
session_info()