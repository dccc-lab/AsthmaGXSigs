## ---- setup_env
options(rgl.useNULL = TRUE) # prevent X11 display warning

library("iCheck")
library("CellMix")
library("dplyr")
library("stringr")
library("broom")
library("forcats")
library("magrittr")
library("genefilter")
library("devtools")
library("readr")

# custom annotation package:
# see "make_hthgu133pluspm_db.R"
library("hthgu133pluspm.db") 


## ---- load_data
base_data_dir <- "data/external"

load(file.path(base_data_dir, "GSE69683/data/GSE69683.rda"))

es_ubiopred_wb <- GSE69683[[1]]

# annotation(es_ubiopred_wb) <- "GPL13158" # aka Affymetrix HT HG-U133+ PM Array
annotation(es_ubiopred_wb) <- "hthgu133pluspm"


## ---- setup_pheno_data
pData(es_ubiopred_wb)$cohort <- es_ubiopred_wb %>% pData %>% use_series(characteristics_ch1) %>% as.character %>% str_replace("cohort: ", "")
pData(es_ubiopred_wb)$sex <- es_ubiopred_wb %>% pData %>% use_series(characteristics_ch1.1) %>% as.character %>% str_replace("gender: ", "")
pData(es_ubiopred_wb)$race <- es_ubiopred_wb %>% pData %>% use_series(characteristics_ch1.2) %>% as.character %>% str_replace("race: ", "") %>% factor %>% fct_inorder
pData(es_ubiopred_wb)$title <- es_ubiopred_wb %>% pData %>% use_series(title) %>% as.character
pData(es_ubiopred_wb)$geo_accession <- es_ubiopred_wb %>% pData %>% use_series(geo_accession) %>% as.character

pData(es_ubiopred_wb)$asthma <- es_ubiopred_wb %>% pData %>% use_series(cohort) %in% c("Severe asthma, non-smoking", "Moderate asthma, non-smoking", "Severe asthma, smoking")
pData(es_ubiopred_wb)$severe <- es_ubiopred_wb %>% pData %>% use_series(cohort) %in% c("Severe asthma, non-smoking", "Severe asthma, smoking")
pData(es_ubiopred_wb)$smoker <- es_ubiopred_wb %>% pData %>% use_series(cohort) %in% c("Severe asthma, smoking")


## ---- load_more_metadata
# incorporate important metadata available not through GEO but through direct data access request to U-BIOPRED
full_metadata <- read_tsv(file.path(base_data_dir, "ubiopred_bigler_metadata.txt"))

full_metadata %<>%
    transmute(
        title = Patient_ID,
        geo_accession = GEO_ID,
        site = Site_Code %>% as.factor,
        WBC = `Wbcs_(xE03_/uL)` %>% na_if("NULL") %>% as.numeric,
        PCTEOS = Eosinophils_Pct %>% na_if("NULL") %>% as.numeric,
        PCTLYMPH = Lymphocytes_Pct %>% na_if("NULL") %>% as.numeric,
        PCTMONO = Monocytes_Pct %>% na_if("NULL") %>% as.numeric,
        PCTNEUT = Neutrophils_Pct %>% na_if("NULL") %>% as.numeric,
        RIN = RIN_RNA_QC %>% na_if("NA") %>% na_if(".") %>% as.numeric
    ) %>% 
    # replace missing metadata values among 36 subjects with median values
    mutate_at(vars(WBC:RIN), ~ifelse(is.na(.), median(., na.rm = TRUE), .))

pData(es_ubiopred_wb) %<>% left_join(full_metadata)


## ---- setup_feature_data
fData(es_ubiopred_wb)$ID <- fData(es_ubiopred_wb)$ID %>% as.character

fData(es_ubiopred_wb)$symbol <- es_ubiopred_wb %>% fData %>% use_series(`Gene Symbol`) %>% as.character # nb: 'DEFA1 /// DEFA1B /// DEFA3'
fData(es_ubiopred_wb)$symbol2 <- unlist(mget(fData(es_ubiopred_wb)$ID, hthgu133pluspmSYMBOL)) # somewhat diff from original GEO anno data

fData(es_ubiopred_wb)$entrez <- unlist(mget(fData(es_ubiopred_wb)$ID, hthgu133pluspmENTREZID))

temp_chr_df <- mget(fData(es_ubiopred_wb)$ID, hthgu133pluspmCHR) %>% unlist %>% data_frame(ID = names(.), chromosome = .) %>% as.data.frame
fData(es_ubiopred_wb) %<>% left_join(temp_chr_df)
rownames(fData(es_ubiopred_wb)) <- fData(es_ubiopred_wb)$ID # put rownames back!
rm(temp_chr_df)


## ---- estimate_wb_cbcs
em1 <- ExpressionMix(es_ubiopred_wb)
res1 <- gedBlood(em1, verbose = TRUE, normalize = "none")
ecbc1 <- asCBC(res1)

cbc_phenos1 <- ecbc1 %>% 
    coef %>% 
    t %>%
    as_data_frame %>% 
    mutate(geo_accession = ecbc1 %>% coef %>% colnames,
           ESTLYMPH = Lymphocytes * 100, 
           ESTMONO = Monocytes * 100, 
           ESTNEUT = Neutrophils * 100) %>% 
    select(geo_accession, starts_with("EST"))

pData(es_ubiopred_wb) %<>% left_join(cbc_phenos1)


## ---- make_gx_pcs_wb
res_pca_ubiopred_wb <- getPCAFunc(es = es_ubiopred_wb, labelVariable = "title", requireLog2 = FALSE, corFlag = TRUE)

# merge principal components (PC) of gene expression data into phenotype data
pDat_wb <- es_ubiopred_wb %>% pData %>% tibble::rownames_to_column()
pcDat_wb <- res_pca_ubiopred_wb %>% use_series(pcs) %>% use_series(x) %>% as.data.frame %>%
    set_colnames(paste0("gx", colnames(.))) %>% tibble::rownames_to_column() %>% dplyr::select(rowname, gxPC1:gxPC10)

pData(es_ubiopred_wb) <- pDat_wb %>% left_join(pcDat_wb, by = c("title" = "rowname")) %>%
    set_rownames(pDat_wb %>% use_series(rowname)) %>% dplyr::select(-rowname) %>% as.data.frame


## ---- clean_up_env
rm(list = ls() %>% extract(grepl("^es_ubiopred|^res_", .) %>% not))