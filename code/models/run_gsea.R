## ---- setup_env
library("here")

source(here("code", "data", "create_ubiopred_eset.R"), echo = TRUE)

library("fgsea")

my.seed <- 40799039 # a fixed random number seed to aid reproducibility
set.seed(my.seed)

rundate_gsea <- format(Sys.Date(), format = "%Y%m%d")


## ---- run_dge
# cases vs. ctrls. (all subjects)
dge1 <- lmFitWrapper(es_ubiopred_wb, 
                     formula = ~asthma + severe + smoker + sex + gxPC1 + gxPC2 + ESTNEUT,
                     probeID.var = "ID",
                     gene.var = "symbol",
                     chr.var = "chromosome")

# cases vs. ctrls. (no smokers)
dge2 <- lmFitWrapper(es_ubiopred_wb[, !es_ubiopred_wb$smoker], 
                     formula = ~asthma + severe + sex + gxPC1 + gxPC2 + ESTNEUT,
                     probeID.var = "ID",
                     gene.var = "symbol",
                     chr.var = "chromosome")

# severe cases vs. moderate cases (no smokers)
dge3 <- lmFitWrapper(es_ubiopred_wb[, !es_ubiopred_wb$smoker & es_ubiopred_wb$asthma], 
                     formula = ~severe + sex + gxPC1 + gxPC2 + ESTNEUT,
                     probeID.var = "ID",
                     gene.var = "symbol",
                     chr.var = "chromosome")

# severe cases vs. ctrls. (no smokers)
dge4 <- lmFitWrapper(es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Severe asthma, non-smoking")], 
                     formula = ~asthma + sex + gxPC1 + gxPC2 + ESTNEUT,
                     probeID.var = "ID",
                     gene.var = "symbol",
                     chr.var = "chromosome")

# moderate cases vs. ctrls. (no smokers)
dge5 <- lmFitWrapper(es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Moderate asthma, non-smoking")], 
                     formula = ~asthma + sex + gxPC1 + gxPC2 + ESTNEUT,
                     probeID.var = "ID",
                     gene.var = "symbol",
                     chr.var = "chromosome")


## ---- extract_dge_results
dge1_rnks <- dge1 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
dge1_rnks <- setNames(dge1_rnks$stats, dge1_rnks$geneSymbols)

dge2_rnks <- dge2 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
dge2_rnks <- setNames(dge2_rnks$stats, dge2_rnks$geneSymbols)

dge3_rnks <- dge3 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
dge3_rnks <- setNames(dge3_rnks$stats, dge3_rnks$geneSymbols)

dge4_rnks <- dge4 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
dge4_rnks <- setNames(dge4_rnks$stats, dge4_rnks$geneSymbols)

dge5_rnks <- dge5 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
dge5_rnks <- setNames(dge5_rnks$stats, dge5_rnks$geneSymbols)


## ---- load_input_gene_sets
# Focus on TREM1-related leading edge gene sets of asthma control
# See Supp. Table 5 from Croteau-Chonka et al. (2017) [https://dx.doi.org/10.1164/rccm.201601-0107OC]
gsc_c7_control <- gmtPathways("data/external/c7.asthmactrl.leading.edges.v4.0.symbols.hs.gmt")
gsc_c7_control_trem1 <- gsc_c7_control[grepl("^GSE9988", names(gsc_c7_control))]


## ---- run_gsea
fgseaRes1 <- fgsea(
    pathways = gsc_c7_control_trem1,
    stats = dge1_rnks,
    minSize = 15,
    maxSize = 500,
    nperm = 100000
)

fgseaRes2 <- fgsea(
    pathways = gsc_c7_control_trem1,
    stats = dge2_rnks,
    minSize = 15,
    maxSize = 500,
    nperm = 100000
)

fgseaRes3 <- fgsea(
    pathways = gsc_c7_control_trem1,
    stats = dge3_rnks,
    minSize = 15,
    maxSize = 500,
    nperm = 100000
)

fgseaRes4 <- fgsea(
    pathways = gsc_c7_control_trem1,
    stats = dge4_rnks,
    minSize = 15,
    maxSize = 500,
    nperm = 100000
)

fgseaRes5 <- fgsea(
    pathways = gsc_c7_control_trem1,
    stats = dge5_rnks,
    minSize = 15,
    maxSize = 500,
    nperm = 100000
)
