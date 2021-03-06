

```r
library("here")

source(here("code", "models", "run_gsea.R"), echo = TRUE)
```

```
## 
## > ## ---- setup_env
## > library("here")
## 
## > source(here("code", "data", "create_ubiopred_eset.R"), echo = TRUE)
## 
## > ## ---- setup_env
## > options(rgl.useNULL = TRUE) # prevent X11 display warning
## 
## > library("iCheck")
```

```
## Loading required package: Biobase
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: lumi
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```
## Loading required package: gplots
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```

```
## 
## > library("CellMix")
```

```
## Loading required package: pkgmaker
```

```
## Loading required package: registry
```

```
## 
## Attaching package: 'pkgmaker'
```

```
## The following object is masked from 'package:base':
## 
##     isNamespaceLoaded
```

```
## Loading required package: NMF
```

```
## Loading required package: rngtools
```

```
## Loading required package: cluster
```

```
## NMF - BioConductor layer [OK] | Shared memory capabilities [NO: synchronicity] | Cores 2/3
```

```
##   To enable shared memory capabilities, try: install.extras('
## NMF
## ')
```

```
## Loading required package: csSAM
```

```
## Loading required package: compiler
```

```
## Loading required package: stringr
```

```
## Loading required package: GSEABase
```

```
## Loading required package: annotate
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:NMF':
## 
##     nrun
```

```
## The following object is masked from 'package:pkgmaker':
## 
##     new2
```

```
## The following object is masked from 'package:gplots':
## 
##     space
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: XML
```

```
## Loading required package: graph
```

```
## 
## Attaching package: 'graph'
```

```
## The following object is masked from 'package:XML':
## 
##     addNode
```

```
## The following object is masked from 'package:stringr':
## 
##     boundary
```

```
## Warning: replacing previous import 'graph::boundary' by 'stringr::boundary'
## when loading 'CellMix'
```

```
## 
## Attaching package: 'CellMix'
```

```
## The following object is masked from 'package:IRanges':
## 
##     reverse
```

```
## 
## > library("dplyr")
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:CellMix':
## 
##     combine, intersect
```

```
## The following objects are masked from 'package:GSEABase':
## 
##     intersect, setdiff, union
```

```
## The following object is masked from 'package:graph':
## 
##     union
```

```
## The following object is masked from 'package:AnnotationDbi':
## 
##     select
```

```
## The following objects are masked from 'package:IRanges':
## 
##     collapse, desc, intersect, setdiff, slice, union
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     first, intersect, rename, setdiff, setequal, union
```

```
## The following object is masked from 'package:lumi':
## 
##     combine
```

```
## The following object is masked from 'package:Biobase':
## 
##     combine
```

```
## The following objects are masked from 'package:BiocGenerics':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following object is masked from '.env':
## 
##     n
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## 
## > library("stringr")
## 
## > library("broom")
## 
## > library("forcats")
## 
## > library("magrittr")
## 
## > library("genefilter")
## 
## > library("devtools")
## 
## > library("readr")
```

```
## 
## Attaching package: 'readr'
```

```
## The following object is masked from 'package:genefilter':
## 
##     spec
```

```
## 
## > # custom annotation package:
## > # see "make_hthgu133pluspm_db.R"
## > library("hthgu133pluspm.db")
```

```
## Loading required package: org.Hs.eg.db
```

```
## 
```

```
## 
```

```
## 
## > ## ---- load_data
## > base_data_dir <- "data/external"
## 
## > load(file.path(base_data_dir, "GSE69683/data/GSE69683.rda"))
## 
## > es_ubiopred_wb <- GSE69683[[1]]
## 
## > # annotation(es_ubiopred_wb) <- "GPL13158" # aka Affymetrix HT HG-U133+ PM Array
## > annotation(es_ubiopred_wb) <- "hthgu133pluspm"
## 
## > ## ---- setup_pheno_data
## > pData(es_ubiopred_wb)$cohort <- es_ubiopred_wb %>% pData %>% use_series(characteristics_ch1) %>% as.character %>% str_rep .... [TRUNCATED] 
## 
## > pData(es_ubiopred_wb)$sex <- es_ubiopred_wb %>% pData %>% use_series(characteristics_ch1.1) %>% as.character %>% str_replace("gender: ", "")
## 
## > pData(es_ubiopred_wb)$race <- es_ubiopred_wb %>% pData %>% use_series(characteristics_ch1.2) %>% as.character %>% str_replace("race: ", "") %>% fact .... [TRUNCATED] 
## 
## > pData(es_ubiopred_wb)$title <- es_ubiopred_wb %>% pData %>% use_series(title) %>% as.character
## 
## > pData(es_ubiopred_wb)$geo_accession <- es_ubiopred_wb %>% pData %>% use_series(geo_accession) %>% as.character
## 
## > pData(es_ubiopred_wb)$asthma <- es_ubiopred_wb %>% pData %>% use_series(cohort) %in% c("Severe asthma, non-smoking", "Moderate asthma, non-smoking", .... [TRUNCATED] 
## 
## > pData(es_ubiopred_wb)$severe <- es_ubiopred_wb %>% pData %>% use_series(cohort) %in% c("Severe asthma, non-smoking", "Severe asthma, smoking")
## 
## > pData(es_ubiopred_wb)$smoker <- es_ubiopred_wb %>% pData %>% use_series(cohort) %in% c("Severe asthma, smoking")
## 
## > ## ---- load_more_metadata
## > # incorporate important metadata available not through GEO but through direct data access request to U-BIOPRED
## > full_m .... [TRUNCATED]
```

```
## Parsed with column specification:
## cols(
##   GEO_ID = col_character(),
##   Patient_ID = col_character(),
##   Site_Code = col_character(),
##   `Wbcs_(xE03_/uL)` = col_character(),
##   Monocytes_Pct = col_character(),
##   Lymphocytes_Pct = col_character(),
##   Neutrophils_Pct = col_character(),
##   Eosinophils_Pct = col_character(),
##   RIN_RNA_QC = col_character()
## )
```

```
## 
## > full_metadata %<>%
## +     transmute(
## +         title = Patient_ID,
## +         geo_accession = GEO_ID,
## +         site = Site_Code %>% as.factor,
## +      .... [TRUNCATED] 
## 
## > pData(es_ubiopred_wb) %<>% left_join(full_metadata)
```

```
## Joining, by = c("title", "geo_accession")
```

```
## 
## > ## ---- setup_feature_data
## > fData(es_ubiopred_wb)$ID <- fData(es_ubiopred_wb)$ID %>% as.character
## 
## > fData(es_ubiopred_wb)$symbol <- es_ubiopred_wb %>% fData %>% use_series(`Gene Symbol`) %>% as.character # nb: 'DEFA1 /// DEFA1B /// DEFA3'
## 
## > fData(es_ubiopred_wb)$symbol2 <- unlist(mget(fData(es_ubiopred_wb)$ID, hthgu133pluspmSYMBOL)) # somewhat diff from original GEO anno data
## 
## > fData(es_ubiopred_wb)$entrez <- unlist(mget(fData(es_ubiopred_wb)$ID, hthgu133pluspmENTREZID))
## 
## > temp_chr_df <- mget(fData(es_ubiopred_wb)$ID, hthgu133pluspmCHR) %>% unlist %>% data_frame(ID = names(.), chromosome = .) %>% as.data.frame
## 
## > fData(es_ubiopred_wb) %<>% left_join(temp_chr_df)
```

```
## Joining, by = "ID"
```

```
## 
## > rownames(fData(es_ubiopred_wb)) <- fData(es_ubiopred_wb)$ID # put rownames back!
## 
## > rm(temp_chr_df)
## 
## > ## ---- estimate_wb_cbcs
## > em1 <- ExpressionMix(es_ubiopred_wb)
## 
## > res1 <- gedBlood(em1, verbose = TRUE, normalize = "none")
## Loading basis signature from Abbas et al. (2009) ... OK [359 features x 17 cell types]
## Estimating proportions for blood cell subset(s): WB
## Mapping signature ids onto target ids (method: auto) ... OK [218 features x 17 cell types]
## Limit/reorder to common set of features ... OK [218 features x 17 cell types]
## Checking data dimension compatibility ... OK [218 features x 17 cell types]
## Using cell type signatures: 'Th', 'Th act', ..., 'neutro' [17 total]
## Checking log-scale ... data:YES - signatures:NO
## Applying log-transform to signatures (base 2) ... OK
## Normalizing signatures and target together (method: none) ... SKIP
##  Using ged algorithm: "lsfit"
##   Estimating cell proportions from cell-specific signatures [lsfit: ls]
##  Timing:
##    user  system elapsed 
##  33.723   0.226  34.085 
##  GED final wrap up ...  OK
## 
## > ecbc1 <- asCBC(res1)
## 
## > cbc_phenos1 <- ecbc1 %>% 
## +     coef %>% 
## +     t %>%
## +     as_data_frame %>% 
## +     mutate(geo_accession = ecbc1 %>% coef %>% colnames,
## +           .... [TRUNCATED] 
## 
## > pData(es_ubiopred_wb) %<>% left_join(cbc_phenos1)
```

```
## Joining, by = "geo_accession"
```

```
## 
## > ## ---- make_gx_pcs_wb
## > res_pca_ubiopred_wb <- getPCAFunc(es = es_ubiopred_wb, labelVariable = "title", requireLog2 = FALSE, corFlag = TRUE)
## 
## > # merge principal components (PC) of gene expression data into phenotype data
## > pDat_wb <- es_ubiopred_wb %>% pData %>% tibble::rownames_to_column()
## 
## > pcDat_wb <- res_pca_ubiopred_wb %>% use_series(pcs) %>% use_series(x) %>% as.data.frame %>%
## +     set_colnames(paste0("gx", colnames(.))) %>% tibble .... [TRUNCATED] 
## 
## > pData(es_ubiopred_wb) <- pDat_wb %>% left_join(pcDat_wb, by = c("title" = "rowname")) %>%
## +     set_rownames(pDat_wb %>% use_series(rowname)) %>% dp .... [TRUNCATED] 
## 
## > ## ---- clean_up_env
## > rm(list = ls() %>% extract(grepl("^es_ubiopred|^res_", .) %>% not))
## 
## > library("fgsea")
```

```
## Loading required package: Rcpp
```

```
## 
## > my.seed <- 40799039 # a fixed random number seed to aid reproducibility
## 
## > set.seed(my.seed)
## 
## > rundate_gsea <- format(Sys.Date(), format = "%Y%m%d")
## 
## > ## ---- run_dge
## > # cases vs. ctrls. (all subjects)
## > dge1 <- lmFitWrapper(es_ubiopred_wb, 
## +                      formula = ~asthma + severe + smok .... [TRUNCATED] 
## dim(dat)>>
## [1] 54715   498
## 
## Running lmFit...
## Running eBayes...
## Preparing output...
##           probeIDs geneSymbols  chr     stats         pval       p.adj
## 1     226999_PM_at       RNPC3    1 -5.345353 1.402865e-07 0.002317314
## 2     209124_PM_at       MYD88    3  5.340417 1.439383e-07 0.002317314
## 3     232441_PM_at        KRR1   12 -5.254879 2.239805e-07 0.002317314
## 4   227280_PM_s_at      CCNYL1    2  5.236354 2.462980e-07 0.002317314
## 5     231698_PM_at    FLJ36848 <NA> -5.231373 2.526577e-07 0.002317314
## 6     201705_PM_at       PSMD7   16  5.215879 2.734722e-07 0.002317314
## 7     244822_PM_at               21 -5.189181 3.132979e-07 0.002317314
## 8   224635_PM_s_at       BIRC6    2  5.173751 3.388195e-07 0.002317314
## 9  1554890_PM_a_at        TIA1    2 -5.102278 4.857451e-07 0.002682127
## 10   1566342_PM_at                6  5.048957 6.337757e-07 0.002682127
## 11  204485_PM_s_at      TOM1L1   17  5.046407 6.418537e-07 0.002682127
## 12    241320_PM_at             <NA> -5.043169 6.522530e-07 0.002682127
## 13   1555845_PM_at             <NA>  5.039312 6.648524e-07 0.002682127
## 14    226219_PM_at    ARHGAP30    1  5.020441 7.299736e-07 0.002682127
## 15    229410_PM_at     SLC35E1   19  5.018970 7.352993e-07 0.002682127
## 16    239160_PM_at             <NA> -4.985051 8.690954e-07 0.002972035
## 17    215316_PM_at             <NA> -4.939373 1.086907e-06 0.003318577
## 18    206249_PM_at     MAP3K13    3 -4.914963 1.224021e-06 0.003318577
## 19    226673_PM_at      SH2D3C    9  4.913194 1.234581e-06 0.003318577
## 20   1570352_PM_at         ATM   11 -4.910937 1.248181e-06 0.003318577
##      pos
## 1  36255
## 2  18539
## 3  41696
## 4  36536
## 5  40953
## 6  11154
## 7  54073
## 8  33895
## 9   1911
## 10  8251
## 11 13933
## 12 50570
## 13  2617
## 14 35476
## 15 38665
## 16 48410
## 17 24611
## 18 15696
## 19 35929
## 20  9801
## 
## pvalue quantiles for intercept and covariates>>
##          (Intercept)   asthmaTRUE   severeTRUE   smokerTRUE   sexmale
## min    1.876594e-178 1.402865e-07 8.163286e-07 1.608934e-12 0.0000000
## 25%     3.308695e-07 1.070850e-01 1.672244e-01 2.486899e-01 0.1026222
## median  2.745956e-04 3.527247e-01 4.222741e-01 4.965248e-01 0.3512727
## 75%     9.504010e-03 6.591279e-01 7.061102e-01 7.457661e-01 0.6626072
## max     9.995002e-01 9.999963e-01 9.999972e-01 9.999655e-01 0.9999922
##                gxPC1        gxPC2          RIN        siteb        sitec
## min    1.917397e-192 3.231879e-90 2.560358e-24 2.220438e-08 1.833634e-10
## 25%     2.849438e-21 6.547070e-12 8.635667e-03 1.622850e-01 1.106517e-01
## median  1.052482e-08 8.241990e-05 1.416665e-01 4.164055e-01 3.525695e-01
## 75%     4.911457e-03 6.296949e-02 4.978617e-01 7.009676e-01 6.595139e-01
## max     9.999733e-01 9.997822e-01 9.999680e-01 9.999809e-01 9.999973e-01
##               sited        sitee        sitef        siteg        siteh
## min    4.357972e-09 5.752869e-19 1.510799e-11 6.437463e-13 4.184373e-13
## 25%    1.526343e-01 6.911670e-02 1.365643e-01 7.029154e-02 1.056053e-01
## median 4.002489e-01 2.963214e-01 3.862277e-01 2.900238e-01 3.425045e-01
## 75%    6.904556e-01 6.226426e-01 6.819283e-01 6.139466e-01 6.543753e-01
## max    9.999727e-01 9.999888e-01 9.999954e-01 9.999837e-01 9.999736e-01
##               sitei        sitej        sitek        sitel        sitem
## min    7.081163e-10 1.409959e-10 2.267729e-09 3.077715e-09 3.228545e-07
## 25%    1.312245e-01 1.340553e-01 1.320509e-01 1.575869e-01 1.904151e-01
## median 3.736414e-01 3.837946e-01 3.773356e-01 4.109016e-01 4.388570e-01
## 75%    6.745540e-01 6.827344e-01 6.804147e-01 6.992246e-01 7.172319e-01
## max    9.999745e-01 9.999731e-01 9.999846e-01 9.999975e-01 9.999907e-01
##               siten        siteo       PCTEOS     PCTLYMPH      PCTMONO
## min    5.672569e-07 6.522614e-10 1.028078e-07 1.814048e-06 1.265702e-11
## 25%    2.078192e-01 1.892494e-01 2.844465e-01 2.088326e-01 1.932834e-01
## median 4.592334e-01 4.449223e-01 5.311246e-01 4.608083e-01 4.512684e-01
## 75%    7.273473e-01 7.192073e-01 7.679577e-01 7.311528e-01 7.222551e-01
## max    9.999999e-01 9.999657e-01 9.999991e-01 9.999798e-01 9.999479e-01
##             PCTNEUT          WBC
## min    1.507858e-05 5.836582e-19
## 25%    2.893692e-01 3.879432e-02
## median 5.348120e-01 2.480734e-01
## 75%    7.695300e-01 5.939337e-01
## max    9.999837e-01 9.999999e-01
## 
## formula>>
## ~asthma + severe + smoker + sex + gxPC1 + gxPC2 + RIN + site + 
##     PCTEOS + PCTLYMPH + PCTMONO + PCTNEUT + WBC
## 
## covariate of interest is  asthma 
## Number of tests= 54715 
## Number of arrays= 498 
## Number of significant tests (raw p-value <  0.05 )= 9031 
## Number of significant tests after p-value adjustments= 1177 
## 
## 
## **********************************************
## 
## > # cases vs. ctrls. (no smokers)
## > dge2 <- lmFitWrapper(es_ubiopred_wb[, !es_ubiopred_wb$smoker], 
## +                      formula = ~asthma + severe  .... [TRUNCATED] 
## dim(dat)>>
## [1] 54715   410
## 
## Running lmFit...
## Running eBayes...
## Preparing output...
##           probeIDs geneSymbols  chr     stats         pval       p.adj
## 1     231698_PM_at    FLJ36848 <NA> -5.563610 4.926384e-08 0.001052156
## 2     226999_PM_at       RNPC3    1 -5.529819 5.892018e-08 0.001052156
## 3    1554364_PM_at     PPP2R5C   14 -5.457328 8.624801e-08 0.001052156
## 4     209124_PM_at       MYD88    3  5.435561 9.662755e-08 0.001052156
## 5   224635_PM_s_at       BIRC6    2  5.345799 1.537924e-07 0.001052156
## 6   227280_PM_s_at      CCNYL1    2  5.307695 1.869815e-07 0.001052156
## 7   204485_PM_s_at      TOM1L1   17  5.306788 1.878505e-07 0.001052156
## 8     232441_PM_at        KRR1   12 -5.297357 1.971233e-07 0.001052156
## 9     226219_PM_at    ARHGAP30    1  5.292626 2.019407e-07 0.001052156
## 10   1566342_PM_at                6  5.265591 2.317399e-07 0.001052156
## 11    201705_PM_at       PSMD7   16  5.250365 2.503563e-07 0.001052156
## 12    229410_PM_at     SLC35E1   19  5.242877 2.600362e-07 0.001052156
## 13   1555845_PM_at             <NA>  5.224020 2.860469e-07 0.001052156
## 14 1554890_PM_a_at        TIA1    2 -5.216546 2.970400e-07 0.001052156
## 15    241320_PM_at             <NA> -5.214096 3.007311e-07 0.001052156
## 16    226673_PM_at      SH2D3C    9  5.209565 3.076760e-07 0.001052156
## 17    202940_PM_at        WNK1   12 -5.127527 4.638939e-07 0.001432373
## 18   1570352_PM_at         ATM   11 -5.124376 4.712185e-07 0.001432373
## 19    226976_PM_at       KPNA6    1  5.085176 5.722403e-07 0.001551392
## 20    239160_PM_at             <NA> -5.083443 5.771576e-07 0.001551392
##      pos
## 1  40953
## 2  36255
## 3   1515
## 4  18539
## 5  33895
## 6  36536
## 7  13933
## 8  41696
## 9  35476
## 10  8251
## 11 11154
## 12 38665
## 13  2617
## 14  1911
## 15 50570
## 16 35929
## 17 12390
## 18  9801
## 19 36232
## 20 48410
## 
## pvalue quantiles for intercept and covariates>>
##          (Intercept)   asthmaTRUE   severeTRUE   sexmale         gxPC1
## min    7.686753e-151 4.926384e-08 9.646046e-07 0.0000000 6.735356e-156
## 25%     9.951420e-07 9.783880e-02 1.699683e-01 0.1199723  5.215142e-18
## median  4.906777e-04 3.411243e-01 4.230084e-01 0.3756570  1.597667e-07
## 75%     1.266049e-02 6.502384e-01 7.072330e-01 0.6783356  9.720785e-03
## max     9.974904e-01 9.999893e-01 9.999989e-01 0.9999887  9.999960e-01
##               gxPC2          RIN        siteb        sitec        sited
## min    9.224360e-78 3.705286e-23 4.856296e-09 3.795600e-10 7.449322e-09
## 25%    8.343417e-10 1.542207e-02 1.587523e-01 1.140202e-01 1.334662e-01
## median 4.466917e-04 1.720762e-01 4.089722e-01 3.608276e-01 3.770632e-01
## 75%    1.000411e-01 5.211883e-01 7.006320e-01 6.673619e-01 6.754138e-01
## max    9.997788e-01 9.999687e-01 9.998527e-01 9.999597e-01 9.999775e-01
##               sitee        sitef        siteg        siteh        sitei
## min    5.252842e-18 4.152717e-10 3.132139e-12 2.884281e-10 6.705243e-10
## 25%    8.800256e-02 1.408952e-01 9.168549e-02 1.192017e-01 1.380445e-01
## median 3.298749e-01 3.873301e-01 3.271199e-01 3.624866e-01 3.836885e-01
## 75%    6.496997e-01 6.856562e-01 6.432604e-01 6.711982e-01 6.813581e-01
## max    9.999544e-01 9.999557e-01 9.998571e-01 9.999909e-01 9.999912e-01
##               sitej        sitek        sitel        sitem        siten
## min    5.031622e-10 5.526412e-10 8.603977e-09 1.006798e-08 3.673556e-07
## 25%    1.283551e-01 1.370766e-01 1.454662e-01 1.755915e-01 2.009133e-01
## median 3.751173e-01 3.855738e-01 3.992734e-01 4.273501e-01 4.524106e-01
## 75%    6.783809e-01 6.819015e-01 6.894302e-01 7.118065e-01 7.237052e-01
## max    9.999930e-01 9.999763e-01 9.999417e-01 9.999654e-01 9.999841e-01
##               siteo       PCTEOS     PCTLYMPH      PCTMONO      PCTNEUT
## min    9.929998e-09 2.027448e-08 2.150033e-05 2.983990e-12 2.332702e-05
## 25%    1.945914e-01 2.824090e-01 2.177349e-01 1.817328e-01 2.949015e-01
## median 4.499866e-01 5.330661e-01 4.675438e-01 4.389241e-01 5.423018e-01
## 75%    7.232283e-01 7.683422e-01 7.310959e-01 7.158968e-01 7.765997e-01
## max    9.999719e-01 9.999854e-01 9.999823e-01 9.999977e-01 9.999906e-01
##                 WBC
## min    1.613369e-13
## 25%    4.959764e-02
## median 2.643108e-01
## 75%    6.019941e-01
## max    9.999091e-01
## 
## formula>>
## ~asthma + severe + sex + gxPC1 + gxPC2 + RIN + site + PCTEOS + 
##     PCTLYMPH + PCTMONO + PCTNEUT + WBC
## 
## covariate of interest is  asthma 
## Number of tests= 54715 
## Number of arrays= 410 
## Number of significant tests (raw p-value <  0.05 )= 9682 
## Number of significant tests after p-value adjustments= 1748 
## 
## 
## **********************************************
## 
## > # severe cases vs. moderate cases (no smokers)
## > dge3 <- lmFitWrapper(es_ubiopred_wb[, !es_ubiopred_wb$smoker & es_ubiopred_wb$asthma], 
## +           .... [TRUNCATED] 
## dim(dat)>>
## [1] 54715   323
## 
## Running lmFit...
## Running eBayes...
## Preparing output...
##          probeIDs geneSymbols  chr     stats         pval      p.adj   pos
## 1    236131_PM_at             <NA> -4.908480 1.502375e-06 0.06943988 45381
## 2    231819_PM_at                2  4.696637 4.015808e-06 0.06943988 41074
## 3    214038_PM_at        CCL8   17 -4.651314 4.933963e-06 0.06943988 23338
## 4  216598_PM_s_at        CCL2   17 -4.645018 5.076479e-06 0.06943988 25889
## 5    225316_PM_at      MFSD2A    1 -4.494357 9.943767e-06 0.10007843 34574
## 6  216262_PM_s_at       TGIF2   20 -4.399363 1.505617e-05 0.10007843 25555
## 7    232379_PM_at        SKIL    3  4.350590 1.857908e-05 0.10007843 41634
## 8     35436_PM_at      GOLGA2    9  4.292392 2.381848e-05 0.10007843 54225
## 9    226773_PM_at                4 -4.261693 2.712423e-05 0.10007843 36029
## 10   200989_PM_at       HIF1A   14 -4.222680 3.196098e-05 0.10007843 10438
## 11   209684_PM_at        RIN2   20 -4.208750 3.387955e-05 0.10007843 19093
## 12   218856_PM_at    TNFRSF21    6 -4.206723 3.416774e-05 0.10007843 28141
## 13 200732_PM_s_at      PTP4A1    6 -4.205076 3.440360e-05 0.10007843 10181
## 14 226571_PM_s_at       PTPRS   19 -4.194921 3.589266e-05 0.10007843 35828
## 15   214370_PM_at      S100A8    1  4.188940 3.679822e-05 0.10007843 23670
## 16 229465_PM_s_at               19 -4.184672 3.745775e-05 0.10007843 38720
## 17   225136_PM_at     PLEKHA2    8 -4.180339 3.813875e-05 0.10007843 34394
## 18   215029_PM_at                1  4.164009 4.081334e-05 0.10007843 24324
## 19 202917_PM_s_at      S100A8    1  4.159797 4.153149e-05 0.10007843 12367
## 20 200733_PM_s_at      PTP4A1    6 -4.152425 4.281796e-05 0.10007843 10182
## 
## pvalue quantiles for intercept and covariates>>
##          (Intercept)   severeTRUE   sexmale         gxPC1        gxPC2
## min    2.521140e-110 1.502375e-06 0.0000000 3.058848e-129 4.346663e-61
## 25%     2.249605e-05 1.797369e-01 0.1618508  5.757109e-15 4.788662e-08
## median  2.558141e-03 4.355219e-01 0.4190819  2.784846e-06 1.598463e-03
## 75%     3.880083e-02 7.124054e-01 0.7050558  2.227387e-02 1.407206e-01
## max     9.978293e-01 9.999999e-01 0.9999845  9.999185e-01 9.999097e-01
##                 RIN        siteb        sitec        sited        sitee
## min    4.281228e-22 7.236509e-08 9.613267e-10 7.281822e-08 2.988004e-15
## 25%    1.730347e-02 1.693761e-01 1.159589e-01 1.507079e-01 1.043983e-01
## median 1.784885e-01 4.188090e-01 3.600134e-01 3.980649e-01 3.459978e-01
## 75%    5.328131e-01 7.037766e-01 6.704034e-01 6.886177e-01 6.535828e-01
## max    9.999506e-01 9.999974e-01 9.999467e-01 9.999845e-01 9.999248e-01
##               sitef        siteg        siteh        sitei        sitej
## min    1.784615e-08 1.574913e-10 4.378004e-09 7.774404e-08 5.794149e-09
## 25%    1.478148e-01 1.042761e-01 1.282602e-01 2.109283e-01 1.409048e-01
## median 3.953336e-01 3.442866e-01 3.782510e-01 4.631356e-01 3.929307e-01
## 75%    6.875856e-01 6.518550e-01 6.786865e-01 7.300777e-01 6.867948e-01
## max    9.999978e-01 9.999652e-01 9.999888e-01 9.999288e-01 9.999770e-01
##               sitek        sitel        sitem        siten        siteo
## min    3.296176e-09 2.955485e-07 1.309945e-08 2.871356e-06 6.857258e-11
## 25%    1.396823e-01 1.786168e-01 1.846225e-01 2.229158e-01 2.149713e-01
## median 3.845806e-01 4.299968e-01 4.376276e-01 4.729022e-01 4.689018e-01
## 75%    6.827955e-01 7.095458e-01 7.142048e-01 7.338887e-01 7.309303e-01
## max    9.999676e-01 9.999948e-01 9.999756e-01 9.999830e-01 9.999611e-01
##              PCTEOS     PCTLYMPH      PCTMONO      PCTNEUT          WBC
## min    8.434901e-07 1.899119e-05 1.350499e-10 1.358092e-05 2.793841e-13
## 25%    2.940593e-01 2.265007e-01 1.826490e-01 3.064206e-01 6.029020e-02
## median 5.385072e-01 4.750638e-01 4.405649e-01 5.501708e-01 2.842416e-01
## 75%    7.730992e-01 7.358124e-01 7.162832e-01 7.792685e-01 6.139853e-01
## max    9.999950e-01 9.999578e-01 9.999828e-01 9.999428e-01 9.999410e-01
## 
## formula>>
## ~severe + sex + gxPC1 + gxPC2 + RIN + site + PCTEOS + PCTLYMPH + 
##     PCTMONO + PCTNEUT + WBC
## 
## covariate of interest is  severe 
## Number of tests= 54715 
## Number of arrays= 323 
## Number of significant tests (raw p-value <  0.05 )= 5157 
## Number of significant tests after p-value adjustments= 0 
## 
## 
## **********************************************
## No genes are differentially expressed!
## 
## > # severe cases vs. ctrls. (no smokers)
## > dge4 <- lmFitWrapper(es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Severe asthma,  ..." ... [TRUNCATED] 
## dim(dat)>>
## [1] 54715   333
## 
## Running lmFit...
## Running eBayes...
## Preparing output...
##          probeIDs geneSymbols  chr     stats         pval        p.adj
## 1  212637_PM_s_at        WWP1    8  5.851709 1.226890e-08 0.0003957693
## 2    232034_PM_at   LOC203274    9  5.816176 1.486803e-08 0.0003957693
## 3    238606_PM_at      ZNF747   16 -5.745772 2.169986e-08 0.0003957693
## 4    231283_PM_at      MGAT4A    2  5.570722 5.471942e-08 0.0007484932
## 5    242431_PM_at             <NA> -5.453456 1.004411e-07 0.0008279190
## 6  205383_PM_s_at      ZBTB20    3 -5.427044 1.150076e-07 0.0008279190
## 7    242920_PM_at             <NA> -5.424339 1.166107e-07 0.0008279190
## 8    235204_PM_at      ENTPD7 <NA> -5.393650 1.363872e-07 0.0008279190
## 9  213795_PM_s_at       PTPRA   20  5.377510 1.480589e-07 0.0008279190
## 10 203743_PM_s_at         TDG   12 -5.373229 1.513148e-07 0.0008279190
## 11 224635_PM_s_at       BIRC6    2  5.326405 1.917830e-07 0.0009539462
## 12   243249_PM_at             <NA> -5.269228 2.555942e-07 0.0010979516
## 13   236696_PM_at       SR140    3 -5.265144 2.608676e-07 0.0010979516
## 14   227290_PM_at               16  5.131951 5.043886e-07 0.0019712588
## 15   232140_PM_at                9  5.113806 5.512272e-07 0.0020106931
## 16 202693_PM_s_at      STK17A    7  5.032407 8.185093e-07 0.0027990460
## 17   214814_PM_at      YTHDC1    4 -5.012763 8.997766e-07 0.0028959575
## 18 211733_PM_x_at        SCP2    1  4.991757 9.953047e-07 0.0029453645
## 19   239294_PM_at                7  4.973604 1.085697e-06 0.0029453645
## 20   217979_PM_at     TSPAN13    7 -4.966537 1.122995e-06 0.0029453645
##      pos
## 1  21943
## 2  41289
## 3  47856
## 4  40538
## 5  51681
## 6  14831
## 7  52170
## 8  44454
## 9  23096
## 10 13191
## 11 33895
## 12 52499
## 13 45946
## 14 36546
## 15 41395
## 16 12142
## 17 24110
## 18 21054
## 19 48544
## 20 27265
## 
## pvalue quantiles for intercept and covariates>>
##          (Intercept)   asthmaTRUE   sexmale         gxPC1        gxPC2
## min    3.916189e-116 1.226890e-08 0.0000000 9.837992e-128 4.959007e-64
## 25%     8.378811e-06 1.182730e-01 0.1206524  1.526137e-15 2.016072e-08
## median  1.689672e-03 3.660940e-01 0.3756914  1.305358e-06 1.271668e-03
## 75%     2.970501e-02 6.710450e-01 0.6770830  1.746576e-02 1.309926e-01
## max     9.994832e-01 9.999737e-01 0.9999952  9.999525e-01 9.999544e-01
##                 RIN        siteb        sitec        sited        sitee
## min    1.519086e-20 3.626020e-07 9.524663e-10 2.420565e-07 3.893431e-15
## 25%    1.889595e-02 1.924424e-01 1.450170e-01 1.518222e-01 1.044831e-01
## median 1.846825e-01 4.424725e-01 3.990977e-01 4.040504e-01 3.544894e-01
## 75%    5.354916e-01 7.177623e-01 6.916693e-01 6.939684e-01 6.657629e-01
## max    9.999678e-01 9.999933e-01 9.999856e-01 9.999741e-01 9.999806e-01
##               sitef        siteg        siteh        sitei        sitej
## min    5.245888e-09 8.526174e-11 1.713772e-08 2.617989e-08 7.345943e-09
## 25%    1.562233e-01 1.183542e-01 1.465190e-01 1.527889e-01 1.713839e-01
## median 4.040378e-01 3.634895e-01 3.936999e-01 4.031482e-01 4.248593e-01
## 75%    6.935538e-01 6.680690e-01 6.878217e-01 6.945879e-01 7.056420e-01
## max    9.999950e-01 9.999734e-01 9.999911e-01 9.999698e-01 9.999871e-01
##               sitek        sitel        sitem        siten        siteo
## min    7.696817e-10 1.715770e-06 4.153912e-07 1.259767e-07 5.478469e-08
## 25%    1.596156e-01 1.923806e-01 1.993316e-01 2.107683e-01 1.961301e-01
## median 4.077163e-01 4.461822e-01 4.539855e-01 4.669852e-01 4.539715e-01
## 75%    6.962227e-01 7.195689e-01 7.240546e-01 7.314150e-01 7.262822e-01
## max    9.999935e-01 9.999865e-01 9.999917e-01 9.999941e-01 9.999841e-01
##              PCTEOS     PCTLYMPH      PCTMONO      PCTNEUT          WBC
## min    1.629514e-07 1.645911e-05 5.671318e-12 2.591793e-05 3.815529e-12
## 25%    2.572100e-01 2.089153e-01 1.806122e-01 2.721224e-01 6.239996e-02
## median 5.081725e-01 4.606191e-01 4.353422e-01 5.223342e-01 2.839399e-01
## 75%    7.553662e-01 7.287652e-01 7.139634e-01 7.611739e-01 6.158444e-01
## max    9.999044e-01 9.999953e-01 9.999598e-01 9.999886e-01 9.999928e-01
## 
## formula>>
## ~asthma + sex + gxPC1 + gxPC2 + RIN + site + PCTEOS + PCTLYMPH + 
##     PCTMONO + PCTNEUT + WBC
## 
## covariate of interest is  asthma 
## Number of tests= 54715 
## Number of arrays= 333 
## Number of significant tests (raw p-value <  0.05 )= 8304 
## Number of significant tests after p-value adjustments= 731 
## 
## 
## **********************************************
## 
## > # moderate cases vs. ctrls. (no smokers)
## > dge5 <- lmFitWrapper(es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Moderate asth ..." ... [TRUNCATED] 
## dim(dat)>>
## [1] 54715   164
## 
## Running lmFit...
## Running eBayes...
## Preparing output...
##           probeIDs geneSymbols  chr     stats         pval       p.adj
## 1     231698_PM_at    FLJ36848 <NA> -5.872801 2.816778e-08 0.001541200
## 2     229410_PM_at     SLC35E1   19  5.395368 2.727541e-07 0.005640470
## 3     217368_PM_at             <NA>  5.368224 3.092645e-07 0.005640470
## 4    1552348_PM_at      PRSS33   16  5.197862 6.743521e-07 0.008073822
## 5   224365_PM_s_at       TIGD7   16  5.177987 7.378070e-07 0.008073822
## 6   221542_PM_s_at      ERLIN2    8  5.092468 1.083736e-06 0.009418126
## 7     226673_PM_at      SH2D3C    9  5.068734 1.204914e-06 0.009418126
## 8     209169_PM_at       GPM6B    X -4.992904 1.686976e-06 0.010901974
## 9     212384_PM_at        BAT1 <NA> -4.979059 1.793252e-06 0.010901974
## 10  234902_PM_s_at      ZNF416   19  4.924296 2.281004e-06 0.012480511
## 11 1554890_PM_a_at        TIA1    2 -4.857872 3.046894e-06 0.013362356
## 12  221705_PM_s_at       SIKE1    1 -4.853700 3.102543e-06 0.013362356
## 13  218471_PM_s_at        BBS1   11  4.833296 3.389184e-06 0.013362356
## 14  221543_PM_s_at      ERLIN2    8  4.817230 3.632788e-06 0.013362356
## 15    230742_PM_at                3 -4.813288 3.695093e-06 0.013362356
## 16    213517_PM_at       PCBP2   12 -4.796196 3.977459e-06 0.013362356
## 17    226976_PM_at       KPNA6    1  4.782778 4.213659e-06 0.013362356
## 18    232441_PM_at        KRR1   12 -4.769693 4.457057e-06 0.013362356
## 19   1554769_PM_at      ZNF785   16 -4.758247 4.681079e-06 0.013362356
## 20   1554364_PM_at     PPP2R5C   14 -4.739356 5.074830e-06 0.013362356
##      pos
## 1  40953
## 2  38665
## 3  26654
## 4     75
## 5  33632
## 6  30826
## 7  35929
## 8  18584
## 9  21690
## 10 44152
## 11  1911
## 12 30986
## 13 27757
## 14 30827
## 15 39997
## 16 22819
## 17 36232
## 18 41696
## 19  1816
## 20  1515
## 
## pvalue quantiles for intercept and covariates>>
##         (Intercept)   asthmaTRUE       sexmale        gxPC1        gxPC2
## min    4.548118e-62 2.816778e-08 1.743427e-172 2.611468e-60 4.775337e-32
## 25%    4.556073e-04 1.170046e-01  1.927255e-01 1.651540e-06 2.552346e-04
## median 1.591470e-02 3.614972e-01  4.479254e-01 3.824653e-03 4.188633e-02
## 75%    1.232032e-01 6.648375e-01  7.206858e-01 1.453659e-01 3.408624e-01
## max    9.999841e-01 9.999933e-01  9.999568e-01 9.999185e-01 9.998952e-01
##                 RIN        siteb        sitec        sited        sitee
## min    4.964818e-09 1.559741e-07 2.870973e-07 8.623151e-08 2.652541e-09
## 25%    1.423972e-01 1.670616e-01 1.799586e-01 1.430584e-01 1.121166e-01
## median 3.950481e-01 4.175070e-01 4.310888e-01 3.967284e-01 3.659954e-01
## 75%    6.895926e-01 7.023244e-01 7.087992e-01 6.870643e-01 6.695603e-01
## max    9.999909e-01 9.999853e-01 9.999768e-01 9.999343e-01 9.999719e-01
##               sitef        siteg        siteh        sitei        sitej
## min    1.842574e-06 1.115511e-06 2.380054e-08 1.891479e-07 1.202880e-07
## 25%    1.624466e-01 1.702546e-01 1.461369e-01 1.959878e-01 1.757475e-01
## median 4.102621e-01 4.210903e-01 3.892383e-01 4.486947e-01 4.307166e-01
## 75%    6.981967e-01 7.024633e-01 6.843766e-01 7.190233e-01 7.104289e-01
## max    9.999873e-01 9.999847e-01 9.999962e-01 9.999842e-01 9.999655e-01
##               sitek        sitel        siten        siteo       PCTEOS
## min    3.666615e-10 2.151553e-07 2.573177e-05 4.216041e-08 4.695002e-06
## 25%    1.759791e-01 1.369254e-01 2.471475e-01 2.618065e-01 2.961084e-01
## median 4.290468e-01 3.858078e-01 4.928819e-01 5.089509e-01 5.436444e-01
## 75%    7.106522e-01 6.825835e-01 7.441235e-01 7.530559e-01 7.734195e-01
## max    9.999895e-01 9.999400e-01 9.999755e-01 9.999998e-01 9.999622e-01
##            PCTLYMPH      PCTMONO      PCTNEUT          WBC
## min    1.143520e-06 4.730521e-06 5.294631e-06 4.980111e-07
## 25%    2.653755e-01 2.418420e-01 2.907179e-01 1.440029e-01
## median 5.135958e-01 4.893101e-01 5.367872e-01 3.920409e-01
## 75%    7.584840e-01 7.449413e-01 7.716170e-01 6.881198e-01
## max    9.999850e-01 9.999747e-01 9.999922e-01 9.999731e-01
## 
## formula>>
## ~asthma + sex + gxPC1 + gxPC2 + RIN + site + PCTEOS + PCTLYMPH + 
##     PCTMONO + PCTNEUT + WBC
## 
## covariate of interest is  asthma 
## Number of tests= 54715 
## Number of arrays= 164 
## Number of significant tests (raw p-value <  0.05 )= 8275 
## Number of significant tests after p-value adjustments= 493 
## 
## 
## **********************************************
## 
## > ## ---- extract_dge_results
## > dge1_rnks <- dge1 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
## 
## > dge1_rnks <- setNames(dge1_rnks$stats, dge1_rnks$geneSymbols)
## 
## > dge2_rnks <- dge2 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
## 
## > dge2_rnks <- setNames(dge2_rnks$stats, dge2_rnks$geneSymbols)
## 
## > dge3_rnks <- dge3 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
## 
## > dge3_rnks <- setNames(dge3_rnks$stats, dge3_rnks$geneSymbols)
## 
## > dge4_rnks <- dge4 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
## 
## > dge4_rnks <- setNames(dge4_rnks$stats, dge4_rnks$geneSymbols)
## 
## > dge5_rnks <- dge5 %>% use_series(frame) %>% select(geneSymbols, stats) %>% arrange(stats)
## 
## > dge5_rnks <- setNames(dge5_rnks$stats, dge5_rnks$geneSymbols)
## 
## > ## ---- load_input_gene_sets
## > # Focus on TREM1-related leading edge gene sets of asthma control
## > # See Supp. Table 5 from Croteau-Chonka et al. (2 .... [TRUNCATED] 
## 
## > gsc_c7_control_trem1 <- gsc_c7_control[grepl("^GSE9988", names(gsc_c7_control))]
## 
## > ## ---- run_gsea
## > fgseaRes1 <- fgsea(
## +     pathways = gsc_c7_control_trem1,
## +     stats = dge1_rnks,
## +     minSize = 15,
## +     maxSize = 500,
## +    .... [TRUNCATED] 
## 
## > fgseaRes2 <- fgsea(
## +     pathways = gsc_c7_control_trem1,
## +     stats = dge2_rnks,
## +     minSize = 15,
## +     maxSize = 500,
## +     nperm = 100000
## +  .... [TRUNCATED] 
## 
## > fgseaRes3 <- fgsea(
## +     pathways = gsc_c7_control_trem1,
## +     stats = dge3_rnks,
## +     minSize = 15,
## +     maxSize = 500,
## +     nperm = 100000
## +  .... [TRUNCATED] 
## 
## > fgseaRes4 <- fgsea(
## +     pathways = gsc_c7_control_trem1,
## +     stats = dge4_rnks,
## +     minSize = 15,
## +     maxSize = 500,
## +     nperm = 100000
## +  .... [TRUNCATED] 
## 
## > fgseaRes5 <- fgsea(
## +     pathways = gsc_c7_control_trem1,
## +     stats = dge5_rnks,
## +     minSize = 15,
## +     maxSize = 500,
## +     nperm = 100000
## +  .... [TRUNCATED]
```

```r
trem1_pathways <- fgseaRes1 %>% arrange(desc(NES)) %>% filter(padj < 0.05) %>% use_series(pathway)

pdf(file = here("reports", "publication", "Figure01.pdf"), width = 8, height = 7, pointsize = 8)
plotGseaTable(pathways = gsc_c7_control[trem1_pathways],
              stats = dge1_rnks,
              fgseaRes = fgseaRes1,
              colwidths = c(8.0, 3.0, 0.8, 1.0, 1.0),
              gseaParam = 1)
dev.off()
```

```
## png 
##   2
```

```r
# How many cases and controls?
es_ubiopred_wb %>% pData %$% table(asthma)
```

```
## asthma
## FALSE  TRUE 
##    87   411
```

```r
es_ubiopred_wb %>% pData %$% table(cohort, asthma)
```

```
##                               asthma
## cohort                         FALSE TRUE
##   Healthy, non-smoking            87    0
##   Moderate asthma, non-smoking     0   77
##   Severe asthma, non-smoking       0  246
##   Severe asthma, smoking           0   88
```

```r
es_ubiopred_wb %>% pData %$% table(severe, asthma)
```

```
##        asthma
## severe  FALSE TRUE
##   FALSE    87   77
##   TRUE      0  334
```

```r
es_ubiopred_wb %>% pData %$% table(smoker, asthma)
```

```
##        asthma
## smoker  FALSE TRUE
##   FALSE    87  323
##   TRUE      0   88
```

```r
es_ubiopred_wb %>% pData %>% filter(asthma == TRUE) %$% table(severe)
```

```
## severe
## FALSE  TRUE 
##    77   334
```

```r
es_ubiopred_wb %>% pData %>% filter(asthma == TRUE) %$% table(smoker, severe)
```

```
##        severe
## smoker  FALSE TRUE
##   FALSE    77  246
##   TRUE      0   88
```

```r
es_ubiopred_wb %>% pData %$% table(cohort)
```

```
## cohort
##         Healthy, non-smoking Moderate asthma, non-smoking 
##                           87                           77 
##   Severe asthma, non-smoking       Severe asthma, smoking 
##                          246                           88
```

```r
es_ubiopred_wb[, !es_ubiopred_wb$smoker] %>% pData %$% table(cohort)
```

```
## cohort
##         Healthy, non-smoking Moderate asthma, non-smoking 
##                           87                           77 
##   Severe asthma, non-smoking 
##                          246
```

```r
es_ubiopred_wb[, !es_ubiopred_wb$smoker & es_ubiopred_wb$asthma] %>% pData %$% table(cohort)
```

```
## cohort
## Moderate asthma, non-smoking   Severe asthma, non-smoking 
##                           77                          246
```

```r
es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Severe asthma, non-smoking")] %>% pData %$% table(cohort)
```

```
## cohort
##       Healthy, non-smoking Severe asthma, non-smoking 
##                         87                        246
```

```r
es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Severe asthma, non-smoking")] %>% pData %$% table(cohort)
```

```
## cohort
##       Healthy, non-smoking Severe asthma, non-smoking 
##                         87                        246
```

```r
es_ubiopred_wb[, es_ubiopred_wb$cohort %in% c("Healthy, non-smoking", "Moderate asthma, non-smoking")] %>% pData %$% table(cohort)
```

```
## cohort
##         Healthy, non-smoking Moderate asthma, non-smoking 
##                           87                           77
```

```r
# How many genes tested for differential expression?
dge1 %>% use_series(frame) %>% nrow
```

```
## [1] 54715
```

```r
# How many differentially expressed genes in each comparison?
dge1 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # cases vs. ctrls. (all subjects)
```

```
## [1] 1177
```

```r
dge2 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # cases vs. ctrls. (no smokers)
```

```
## [1] 1748
```

```r
dge3 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # severe cases vs. moderate cases (no smokers)
```

```
## [1] 0
```

```r
dge4 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # severe cases vs. ctrls. (no smokers)
```

```
## [1] 731
```

```r
dge5 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # moderate cases vs. ctrls. (no smokers)
```

```
## [1] 493
```

```r
# How many TREM1 gene sets were tested?
fgseaRes1 %>% use_series(pathway) %>% length
```

```
## [1] 24
```

```r
# How large were the TREM-1 gene sets tested?
gsc_c7_control_trem1 %>% sapply(length) %>% summary
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   10.00   42.00   44.00   45.92   50.00   62.00
```

```r
# How many TREM1 gene sets were enriched?
trem1_pathways %>% length
```

```
## [1] 24
```

```r
# What are the core enriched TREM1 pathway genes?
core_genes <- fgseaRes1 %>% arrange(NES) %>% filter(padj < 0.05) %>% use_series(leadingEdge) %>% unlist %>% unique %>% sort
core_genes %>% print
```

```
##   [1] "ACOX1"        "ACSL1"        "ACSL3"        "ADO"         
##   [5] "AKIRIN2"      "ANKRD57"      "AQP9"         "ARAP3"       
##   [9] "ARID5A"       "ARL5B"        "ARRB1"        "ATF3"        
##  [13] "BCAR3"        "BID"          "C10orf128"    "C14orf102"   
##  [17] "C17orf79"     "C1orf122"     "C1orf38"      "C6orf1"      
##  [21] "CALHM2"       "CALM2"        "CASP5"        "CBL"         
##  [25] "CCL23"        "CCNG2"        "CCR2"         "CD300LB"     
##  [29] "CD68"         "CDKN1A"       "CEP350"       "CFLAR"       
##  [33] "CHST15"       "CRBN"         "CSNK1D"       "CTSB"        
##  [37] "CYTH4"        "DAB2"         "DCP2"         "DDIT3"       
##  [41] "DDX17"        "DOK2"         "DR1"          "DRAM1"       
##  [45] "DUSP18"       "DYRK3"        "EFR3A"        "EHD1"        
##  [49] "ELF2"         "ENG"          "ETS2"         "EVI2A"       
##  [53] "EVI2B"        "FAM105A"      "FAM177A1"     "FKBP15"      
##  [57] "FOXN2"        "FRAT1"        "FRAT2"        "FUCA1"       
##  [61] "GAPT"         "GBP1"         "GCLM"         "GDI2"        
##  [65] "GFOD1"        "GLA"          "GLUL"         "GRB2"        
##  [69] "GRN"          "HBEGF"        "HCK"          "HECA"        
##  [73] "HHEX"         "HSD3B7"       "IFI30"        "IFIH1"       
##  [77] "IFNGR2"       "IL1B"         "IL8"          "ITGAX"       
##  [81] "JMJD1C"       "KCNJ2"        "KDM6B"        "KLF6"        
##  [85] "KMO"          "LAIR1"        "LDLR"         "LFNG"        
##  [89] "LOC100129034" "LONRF3"       "LRRC25"       "LRRC33"      
##  [93] "LYSMD2"       "MAP3K14"      "MAP3K8"       "MAPK13"      
##  [97] "MARCKS"       "MBD2"         "MKL1"         "MXD1"        
## [101] "MYD88"        "N4BP1"        "NACC2"        "NCKAP1L"     
## [105] "NCOR2"        "NFATC1"       "NFIC"         "NFKBIZ"      
## [109] "NIN"          "NINJ1"        "NPEPPS"       "OLIG1"       
## [113] "OLIG2"        "PARP10"       "PDK4"         "PELI1"       
## [117] "PHACTR2"      "PHC2"         "PIK3AP1"      "PILRA"       
## [121] "PLAUR"        "PLEK"         "PLEKHO2"      "PLSCR1"      
## [125] "PLXNB2"       "PMAIP1"       "PPP2R5A"      "PSTPIP2"     
## [129] "PTAFR"        "PTGER4"       "PTGS2"        "PTPN6"       
## [133] "PXN"          "RAB11FIP4"    "RAB31"        "RAPGEF2"     
## [137] "RARA"         "RASSF2"       "RCOR3"        "RGS1"        
## [141] "RGS12"        "RGS19"        "RHOB"         "RHOH"        
## [145] "RHOU"         "RIN2"         "RIN3"         "RIT1"        
## [149] "RNF144B"      "RNF19B"       "RRAGC"        "RRAGD"       
## [153] "RYBP"         "S1PR3"        "SCO2"         "SERPINB1"    
## [157] "SESTD1"       "SH2D3C"       "SIPA1"        "SIRPA"       
## [161] "SLA"          "SLC16A3"      "SLC2A3"       "SLC2A6"      
## [165] "SLC38A2"      "SLC43A2"      "SNX20"        "SOD2"        
## [169] "SORL1"        "SPAG9"        "ST3GAL6"      "ST8SIA4"     
## [173] "STK17B"       "STX3"         "STX7"         "SYN2"        
## [177] "TAF5"         "TGFBR1"       "THBS1"        "TIFA"        
## [181] "TIMP2"        "TLE3"         "TLR1"         "TLR4"        
## [185] "TLR5"         "TLR7"         "TMBIM1"       "TMEM123"     
## [189] "TMEM88"       "TNF"          "TNFAIP2"      "TNFSF13B"    
## [193] "TNRC6A"       "TRIB1"        "TRIM33"       "TRIM38"      
## [197] "TSC22D3"      "UBE2W"        "USP22"        "VAV1"        
## [201] "VMO1"         "WDFY1"        "ZBTB2"        "ZC3HAV1"     
## [205] "ZFP36"        "ZFP36L2"      "ZFYVE16"      "ZNF227"      
## [209] "ZNF295"       "ZNF398"
```

```r
core_genes %>% length
```

```
## [1] 210
```

```r
# What is the overlap of this core with highlighted TREM1 genes from Table 3 of Croteau-Chonka et al. (2017)?
table3_genes <- c("CCL23", "OLIG1", "OLIG2", "GFOD1", "RHOBTB3", "HSD3B7")
intersect(core_genes, table3_genes)
```

```
## [1] "CCL23"  "GFOD1"  "HSD3B7" "OLIG1"  "OLIG2"
```

```r
dge1 %>% use_series(frame) %>% filter(geneSymbols %in% table3_genes)
```

```
##          probeIDs geneSymbols chr       stats         pval      p.adj
## 1    213825_PM_at       OLIG2  21  3.58553293 0.0003711393 0.02986388
## 2  210549_PM_s_at       CCL23  17  3.31400806 0.0009898511 0.04740000
## 3    210548_PM_at       CCL23  17  2.97668243 0.0030623328 0.08023879
## 4  202975_PM_s_at     RHOBTB3   5 -1.77133457 0.0771450373 0.36887099
## 5    228170_PM_at       OLIG1  21  1.73756507 0.0829344160 0.38126000
## 6    213824_PM_at       OLIG2  21  1.53356266 0.1258017999 0.46002223
## 7    222817_PM_at      HSD3B7  16  1.49018064 0.1368392202 0.47734510
## 8    216049_PM_at     RHOBTB3   5 -1.29303350 0.1966268353 0.55622155
## 9    240111_PM_at     RHOBTB3   5 -1.28483524 0.1994748878 0.55920533
## 10 219821_PM_s_at       GFOD1   6  0.97971027 0.3277269469 0.68522441
## 11 202976_PM_s_at     RHOBTB3   5  0.63151103 0.5280096617 0.81674615
## 12   225202_PM_at     RHOBTB3   5  0.38270546 0.7021090359 0.89758915
## 13 216048_PM_s_at     RHOBTB3   5 -0.07930557 0.9368229124 0.98233549
##      pos
## 1  23126
## 2  19938
## 3  19937
## 4  12425
## 5  37425
## 6  23125
## 7  32097
## 8  25342
## 9  49361
## 10 29106
## 11 12426
## 12 34460
## 13 25341
```

```r
# Which primary analysis gene sets were enriched in secondary analyses?
fgseaRes1 %>% arrange(desc(NES)) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
##                                                       pathway     pval
## 1               GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 2.06e-05
## 2                 GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 2.07e-05
## 3             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 2.07e-05
## 4           GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 2.07e-05
## 5     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 2.06e-05
## 6                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.07e-05
## 7                  GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.07e-05
## 8                       GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 2.07e-05
## 9                       GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 2.07e-05
## 10             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.06e-05
## 11 GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 2.07e-05
## 12         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 2.06e-05
## 13                  GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 2.07e-05
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.07e-05
## 15       GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.06e-05
## 16                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 4.13e-05
## 17              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 6.20e-05
## 18            GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 1.24e-04
## 19                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 1.86e-04
## 20    GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 5.58e-04
## 21             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 6.82e-04
## 22 GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 2.32e-03
## 23                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 2.00e-03
## 24         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 2.94e-03
##        padj    ES  NES nMoreExtreme size
## 1  3.31e-05 0.706 2.87            0   59
## 2  3.31e-05 0.724 2.75            0   44
## 3  3.31e-05 0.713 2.71            0   44
## 4  3.31e-05 0.654 2.61            0   55
## 5  3.31e-05 0.695 2.60            0   41
## 6  3.31e-05 0.657 2.55            0   48
## 7  3.31e-05 0.664 2.50            0   42
## 8  3.31e-05 0.630 2.49            0   52
## 9  3.31e-05 0.601 2.46            0   62
## 10 3.31e-05 0.649 2.46            0   43
## 11 3.31e-05 0.619 2.35            0   44
## 12 3.31e-05 0.627 2.34            0   40
## 13 3.31e-05 0.567 2.32            0   62
## 14 3.31e-05 0.589 2.21            0   42
## 15 3.31e-05 0.526 2.13            0   58
## 16 6.19e-05 0.543 2.12            1   50
## 17 8.75e-05 0.552 2.10            2   44
## 18 1.65e-04 0.539 2.05            5   44
## 19 2.35e-04 0.518 1.99            8   46
## 20 6.70e-04 0.485 1.89           26   49
## 21 7.79e-04 0.493 1.88           32   45
## 22 2.43e-03 0.486 1.79          112   38
## 23 2.19e-03 0.470 1.79           96   44
## 24 2.94e-03 0.468 1.76          141   42
```

```r
fgseaRes2 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
##                                                       pathway     pval
## 1     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 2.10e-05
## 2     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 2.10e-04
## 3  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 2.09e-05
## 4  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 1.40e-03
## 5        GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.11e-05
## 6             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 2.09e-05
## 7             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 8.36e-05
## 8                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 2.11e-05
## 9                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.09e-05
## 10                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 2.11e-05
## 11                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 2.10e-05
## 12         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 2.09e-05
## 13         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 1.15e-03
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.09e-05
## 15          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 2.10e-05
## 16                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 2.09e-05
## 17                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.30e-03
## 18             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.09e-05
## 19             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 5.02e-04
## 20                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 2.10e-05
## 21              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 4.18e-05
## 22              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 2.11e-05
## 23                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.09e-05
## 24                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 1.26e-04
##        padj    ES  NES nMoreExtreme size
## 1  3.17e-05 0.719 2.70            0   41
## 2  2.52e-04 0.504 1.97            9   49
## 3  3.17e-05 0.650 2.48            0   44
## 4  1.40e-03 0.501 1.85           66   38
## 5  3.17e-05 0.532 2.15            0   58
## 6  3.17e-05 0.719 2.74            0   44
## 7  1.12e-04 0.554 2.11            3   44
## 8  3.17e-05 0.578 2.37            0   62
## 9  3.17e-05 0.649 2.52            0   48
## 10 3.17e-05 0.589 2.41            0   62
## 11 3.17e-05 0.638 2.52            0   52
## 12 3.17e-05 0.648 2.42            0   40
## 13 1.26e-03 0.492 1.86           54   42
## 14 3.17e-05 0.592 2.23            0   42
## 15 3.17e-05 0.665 2.66            0   55
## 16 3.17e-05 0.741 2.83            0   44
## 17 1.35e-03 0.484 1.84           61   44
## 18 3.17e-05 0.665 2.52            0   43
## 19 5.74e-04 0.501 1.92           23   45
## 20 3.17e-05 0.549 2.15            0   50
## 21 5.90e-05 0.563 2.15            1   44
## 22 3.17e-05 0.702 2.85            0   59
## 23 3.17e-05 0.685 2.58            0   42
## 24 1.59e-04 0.525 2.02            5   46
```

```r
fgseaRes3 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
## [1] pathway      pval         padj         ES           NES         
## [6] nMoreExtreme size        
## <0 rows> (or 0-length row.names)
```

```r
fgseaRes4 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
##                                                       pathway     pval
## 1     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 2.33e-05
## 2     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 4.88e-03
## 3  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 2.35e-05
## 4  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 5.83e-03
## 5        GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.85e-02
## 6             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 2.35e-05
## 7             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 9.15e-04
## 8                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 3.38e-04
## 9                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.36e-05
## 10                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 6.03e-04
## 11                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 2.37e-05
## 12         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 2.33e-05
## 13          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 4.67e-05
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 2.38e-05
## 15                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 2.35e-05
## 16                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 7.08e-03
## 17             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.34e-05
## 18             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 2.33e-02
## 19                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.66e-04
## 20              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 4.69e-05
## 21              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 2.39e-05
## 22                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.33e-05
## 23                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 6.28e-03
##        padj    ES  NES nMoreExtreme size
## 1  5.23e-05 0.671 2.56            0   41
## 2  6.51e-03 0.420 1.67          205   49
## 3  5.23e-05 0.601 2.33            0   44
## 4  7.37e-03 0.450 1.69          250   38
## 5  2.98e-02 0.352 1.45         1192   58
## 6  5.23e-05 0.698 2.71            0   44
## 7  1.29e-03 0.476 1.85           38   44
## 8  5.40e-04 0.439 1.83           13   62
## 9  5.23e-05 0.597 2.36            0   48
## 10 9.05e-04 0.433 1.80           24   62
## 11 5.23e-05 0.540 2.17            0   52
## 12 5.23e-05 0.648 2.46            0   40
## 13 8.66e-05 0.556 2.13            1   42
## 14 5.23e-05 0.529 2.15            0   55
## 15 5.23e-05 0.654 2.54            0   44
## 16 8.09e-03 0.426 1.65          301   44
## 17 5.23e-05 0.673 2.59            0   43
## 18 2.54e-02 0.386 1.51          988   45
## 19 2.85e-04 0.491 1.96            6   50
## 20 8.66e-05 0.551 2.14            1   44
## 21 5.23e-05 0.564 2.33            0   59
## 22 5.23e-05 0.697 2.67            0   42
## 23 7.53e-03 0.424 1.66          266   46
```

```r
fgseaRes5 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
##                                                       pathway     pval
## 1     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 2.04e-05
## 2     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 7.16e-04
## 3  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 2.04e-05
## 4  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 2.55e-03
## 5        GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 4.09e-05
## 6             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 2.04e-05
## 7             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 1.64e-04
## 8                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 2.05e-05
## 9                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.05e-05
## 10                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 2.05e-05
## 11                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 2.04e-05
## 12         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 2.04e-05
## 13         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 3.94e-03
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 1.02e-04
## 15          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 2.04e-05
## 16                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 2.04e-05
## 17                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.23e-03
## 18             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.04e-05
## 19             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 2.66e-04
## 20                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 2.04e-05
## 21              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 1.43e-04
## 22              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 2.04e-05
## 23                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.04e-05
## 24                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 6.14e-05
##        padj    ES  NES nMoreExtreme size
## 1  3.51e-05 0.718 2.70            0   41
## 2  8.18e-04 0.468 1.83           34   49
## 3  3.51e-05 0.653 2.49            0   44
## 4  2.66e-03 0.481 1.78          124   38
## 5  6.55e-05 0.550 2.23            1   58
## 6  3.51e-05 0.720 2.75            0   44
## 7  2.07e-04 0.521 1.99            7   44
## 8  3.51e-05 0.611 2.51            0   62
## 9  3.51e-05 0.582 2.27            0   48
## 10 3.51e-05 0.607 2.49            0   62
## 11 3.51e-05 0.577 2.29            0   52
## 12 3.51e-05 0.655 2.45            0   40
## 13 3.94e-03 0.453 1.71          192   42
## 14 1.44e-04 0.539 2.04            4   42
## 15 3.51e-05 0.657 2.63            0   55
## 16 3.51e-05 0.736 2.81            0   44
## 17 1.34e-03 0.476 1.82           59   44
## 18 3.51e-05 0.673 2.56            0   43
## 19 3.19e-04 0.510 1.96           12   45
## 20 3.51e-05 0.548 2.15            0   50
## 21 1.91e-04 0.524 2.00            6   44
## 22 3.51e-05 0.696 2.83            0   59
## 23 3.51e-05 0.654 2.47            0   42
## 24 9.21e-05 0.544 2.10            2   46
```

```r
# Which primary analysis gene sets were NOT enriched in secondary analyses?
fgseaRes3 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj > 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
##                                                       pathway   pval  padj
## 1     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 0.5260 0.752
## 2     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 0.0163 0.110
## 3  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 0.7750 0.775
## 4  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 0.0150 0.110
## 5        GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 0.6260 0.752
## 6             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 0.1280 0.321
## 7             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 0.2970 0.548
## 8                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 0.7060 0.775
## 9                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 0.2690 0.542
## 10                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 0.6240 0.752
## 11                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 0.7670 0.775
## 12         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 0.4210 0.721
## 13         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 0.0253 0.122
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 0.1340 0.321
## 15          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 0.5970 0.752
## 16                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 0.5590 0.752
## 17                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 0.0580 0.232
## 18             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 0.2710 0.542
## 19             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 0.0184 0.110
## 20                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 0.0754 0.259
## 21              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 0.1260 0.321
## 22              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 0.4770 0.752
## 23                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 0.7410 0.775
## 24                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 0.0158 0.110
##        ES    NES nMoreExtreme size
## 1   0.250  0.956        22949   41
## 2   0.388  1.540          701   49
## 3  -0.220 -0.819        43827   44
## 4   0.422  1.580          659   38
## 5   0.222  0.914        26654   58
## 6   0.324  1.260         5568   44
## 7   0.282  1.100        12894   44
## 8   0.210  0.879        29802   62
## 9   0.281  1.110        11619   48
## 10  0.219  0.916        26336   62
## 11  0.208  0.839        32985   52
## 12  0.267  1.020        18435   40
## 13  0.393  1.510         1101   42
## 14  0.326  1.250         5819   42
## 15  0.227  0.926        25476   55
## 16 -0.252 -0.939        31589   44
## 17  0.357  1.390         2519   44
## 18  0.289  1.110        11799   43
## 19  0.396  1.550          798   45
## 20  0.334  1.330         3254   50
## 21  0.325  1.260         5474   44
## 22 -0.249 -0.987        27478   59
## 23 -0.227 -0.837        41845   42
## 24  0.397  1.560          682   46
```

```r
fgseaRes4 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj > 0.05) %>% mutate_at(c("pval", "padj", "ES", "NES"), funs(signif(., 3)))
```

```
##                                              pathway   pval   padj    ES
## 1 GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 0.0951 0.0951 0.342
##    NES nMoreExtreme size
## 1 1.31         4076   42
```

```r
session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.0 (2017-04-21)
##  system   x86_64, linux-gnu           
##  ui       RStudio (1.1.423)           
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       <NA>                        
##  date     2018-02-27
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version    date      
##  affy                   1.56.0     2017-12-12
##  affyio                 1.48.0     2017-12-12
##  annotate             * 1.56.1     2017-12-11
##  AnnotationDbi        * 1.40.0     2017-11-23
##  assertthat             0.2.0      2017-04-11
##  backports              1.1.1      2017-09-25
##  base                 * 3.4.0      2017-05-02
##  base64                 2.0        2016-05-10
##  beanplot               1.2        2014-09-19
##  beeswarm               0.2.3      2016-04-25
##  bibtex                 0.4.2      2017-06-30
##  bigmemory            * 4.5.31     2017-11-20
##  bigmemory.sri        * 0.1.3      2014-08-18
##  bindr                  0.1        2016-11-13
##  bindrcpp             * 0.2        2017-06-17
##  Biobase              * 2.38.0     2017-11-23
##  BiocGenerics         * 0.24.0     2017-11-23
##  BiocInstaller          1.28.0     2017-11-23
##  BiocParallel           1.12.0     2017-11-23
##  biomaRt                2.34.0     2017-11-23
##  Biostrings             2.46.0     2017-11-24
##  bit                    1.1-12     2014-04-09
##  bit64                  0.9-7      2017-05-08
##  bitops                 1.0-6      2013-08-17
##  blob                   1.1.0      2017-06-17
##  broom                * 0.4.3      2017-11-20
##  bumphunter             1.20.0     2017-12-12
##  caTools                1.17.1     2014-09-10
##  CellMix              * 1.6.2      2017-07-25
##  cluster              * 2.0.6      2017-03-16
##  codetools              0.2-15     2016-10-05
##  colorspace             1.3-2      2016-12-14
##  compiler             * 3.4.0      2017-05-02
##  corpcor                1.6.9      2017-04-01
##  csSAM                * 1.2.4      2013-05-13
##  data.table             1.10.4-3   2017-10-27
##  datasets             * 3.4.0      2017-05-02
##  DBI                    0.7        2017-06-18
##  DelayedArray           0.4.1      2017-11-23
##  devtools             * 1.13.4     2017-11-09
##  digest                 0.6.15     2018-01-28
##  doParallel             1.0.11     2017-09-28
##  doRNG                  1.6.6      2017-04-10
##  dplyr                * 0.7.4      2017-09-28
##  evaluate               0.10.1     2017-06-24
##  ezknitr              * 0.6        2016-09-16
##  fastmatch              1.1-0      2017-01-28
##  fgsea                * 1.4.0      2017-11-23
##  forcats              * 0.2.0      2017-01-23
##  foreach                1.4.3      2015-10-13
##  foreign                0.8-69     2017-06-21
##  gdata                  2.18.0     2017-06-06
##  genefilter           * 1.60.0     2017-11-23
##  GeneSelectMMD          2.22.0     2017-11-23
##  GenomeInfoDb           1.14.0     2017-11-23
##  GenomeInfoDbData       0.99.1     2017-12-11
##  GenomicAlignments      1.14.1     2017-11-23
##  GenomicFeatures        1.30.0     2017-11-23
##  GenomicRanges          1.30.0     2017-11-23
##  GEOquery               2.46.11    2017-12-11
##  ggplot2                2.2.1.9000 2018-02-05
##  glue                   1.2.0      2017-10-29
##  gplots               * 3.0.1      2016-03-30
##  graph                * 1.56.0     2017-12-11
##  graphics             * 3.4.0      2017-05-02
##  grDevices            * 3.4.0      2017-05-02
##  grid                   3.4.0      2017-05-02
##  gridBase               0.4-7      2014-02-24
##  gridExtra              2.3        2017-09-09
##  GSEABase             * 1.40.1     2017-11-23
##  gtable                 0.2.0      2016-02-26
##  gtools                 3.5.0      2015-05-29
##  here                 * 0.1        2017-09-26
##  hgu133a.db           * 3.2.3      2017-07-18
##  hgu133b.db           * 3.2.3      2017-07-25
##  hms                    0.4.0      2017-11-23
##  hthgu133pluspm.db    * 0.1        2017-07-20
##  htmltools              0.3.6      2017-04-28
##  htmlwidgets            0.9        2017-07-10
##  httpuv                 1.3.5      2017-07-04
##  iCheck               * 1.8.0      2017-11-23
##  illuminaio             0.20.0     2017-12-11
##  IRanges              * 2.12.0     2017-11-23
##  iterators              1.0.9      2017-12-12
##  jsonlite               1.5        2017-06-01
##  KernSmooth             2.23-15    2015-06-29
##  knitr                  1.17       2017-08-10
##  labeling               0.3        2014-08-23
##  lattice                0.20-35    2017-03-25
##  lazyeval               0.2.1      2017-10-29
##  limma                  3.34.3     2017-12-12
##  limSolve               1.5.5.3    2017-08-14
##  lmtest                 0.9-35     2017-02-11
##  locfit                 1.5-9.1    2013-04-20
##  lpSolve                5.6.13     2015-09-19
##  lumi                 * 2.30.0     2017-12-12
##  magrittr             * 1.5        2014-11-22
##  MASS                   7.3-47     2017-04-21
##  Matrix                 1.2-12     2017-11-16
##  matrixStats            0.52.2     2017-04-14
##  mclust                 5.4        2017-11-22
##  memoise                1.1.0      2017-04-21
##  methods              * 3.4.0      2017-05-02
##  methylumi              2.24.1     2017-12-12
##  mgcv                   1.8-22     2017-09-19
##  mime                   0.5        2016-07-07
##  minfi                  1.24.0     2017-11-23
##  mnormt                 1.5-5      2016-10-15
##  multtest               2.34.0     2017-12-11
##  munsell                0.4.3      2016-02-13
##  nleqslv                3.3.1      2017-07-06
##  nlme                   3.1-131    2017-02-06
##  NMF                  * 0.20.6     2015-05-26
##  nor1mix                1.2-3      2017-08-30
##  openssl                0.9.9      2017-11-10
##  org.Hs.eg.db         * 3.5.0      2017-12-11
##  parallel             * 3.4.0      2017-05-02
##  pillar                 1.1.0      2018-01-14
##  pkgconfig              2.0.1      2017-03-21
##  pkgmaker             * 0.22       2014-05-14
##  plyr                   1.8.4      2016-06-08
##  preprocessCore         1.40.0     2017-11-23
##  prettyunits            1.0.2      2015-07-13
##  progress               1.1.2      2016-12-14
##  psych                  1.7.8      2017-09-09
##  purrr                  0.2.4      2017-10-18
##  quadprog               1.5-5      2013-04-17
##  R.methodsS3            1.7.1      2016-02-16
##  R.oo                   1.21.0     2016-11-01
##  R.utils                2.6.0      2017-11-05
##  R6                     2.2.2      2017-06-17
##  randomForest           4.6-12     2015-10-07
##  RColorBrewer           1.1-2      2014-12-07
##  Rcpp                 * 0.12.14    2017-11-23
##  RCurl                  1.95-4.8   2016-03-01
##  readr                * 1.1.1      2017-05-16
##  registry             * 0.5        2017-12-03
##  reshape                0.8.7      2017-08-06
##  reshape2               1.4.3      2017-12-11
##  rgl                    0.98.1     2017-03-08
##  rlang                  0.1.6.9003 2018-02-05
##  RMySQL                 0.10.13    2017-08-14
##  rngtools             * 1.2.4      2014-03-06
##  rprojroot              1.2        2017-01-16
##  Rsamtools              1.30.0     2017-12-11
##  RSQLite                2.0        2017-06-19
##  rtracklayer            1.38.2     2017-12-12
##  S4Vectors            * 0.16.0     2017-11-23
##  scales                 0.5.0.9000 2018-02-05
##  scatterplot3d          0.3-40     2017-04-22
##  shiny                  1.0.5      2017-08-23
##  siggenes               1.52.0     2017-12-11
##  splines                3.4.0      2017-05-02
##  stats                * 3.4.0      2017-05-02
##  stats4               * 3.4.0      2017-05-02
##  stringi                1.1.6      2017-11-17
##  stringr              * 1.2.0      2017-02-18
##  SummarizedExperiment   1.8.0      2017-11-23
##  survival               2.41-3     2017-04-04
##  tibble                 1.4.2      2018-01-22
##  tidyr                  0.7.2      2017-10-16
##  tools                  3.4.0      2017-05-02
##  utils                * 3.4.0      2017-05-02
##  withr                  2.1.1.9000 2018-02-05
##  XML                  * 3.98-1.9   2017-06-19
##  xml2                   1.1.1      2017-01-24
##  xtable                 1.8-2      2016-02-05
##  XVector                0.18.0     2017-12-11
##  yaml                   2.1.15     2017-12-01
##  zlibbioc               1.24.0     2017-12-12
##  zoo                    1.8-0      2017-04-12
##  source                            
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  local                             
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  local                             
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  local                             
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  local                             
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  cran (@0.6.15)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  Github (thomasp85/ggplot2@f53b99f)
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  local                             
##  local                             
##  local                             
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Github (krlmlr/here@93593ee)      
##  Bioconductor                      
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  Bioconductor                      
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  local                             
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  local                             
##  cran (@1.1.0)                     
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Github (tidyverse/rlang@c6747f9)  
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  Bioconductor                      
##  Github (hadley/scales@d767915)    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  local                             
##  local                             
##  local                             
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  cran (@1.4.2)                     
##  CRAN (R 3.4.0)                    
##  local                             
##  local                             
##  Github (jimhester/withr@df18523)  
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)                    
##  Bioconductor                      
##  CRAN (R 3.4.0)
```

