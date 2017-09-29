

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
## > rm(list = ls() %>% extract(!grepl("^es_", .))) # clean up
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
##  32.979   0.008  33.014 
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
##          probeIDs geneSymbols  chr     stats         pval        p.adj
## 1    209124_PM_at       MYD88    3  6.137541 1.719762e-09 9.409679e-05
## 2  227280_PM_s_at      CCNYL1    2  5.931601 5.639866e-09 1.542926e-04
## 3    232441_PM_at        KRR1   12 -5.603745 3.487982e-08 4.507697e-04
## 4    226999_PM_at       RNPC3    1 -5.579693 3.973425e-08 4.507697e-04
## 5    201705_PM_at       PSMD7   16  5.573025 4.119252e-08 4.507697e-04
## 6    235645_PM_at       ESCO1   18 -5.471018 7.117843e-08 5.830729e-04
## 7   1566342_PM_at                6  5.462199 7.459582e-08 5.830729e-04
## 8    239629_PM_at       CFLAR    2 -5.429584 8.867353e-08 6.064715e-04
## 9    244822_PM_at               21 -5.392511 1.078158e-07 6.554601e-04
## 10   202026_PM_at        SDHD   11  5.319076 1.582718e-07 8.276226e-04
## 11   235216_PM_at       ESCO1   18 -5.309449 1.663867e-07 8.276226e-04
## 12   214814_PM_at      YTHDC1    4 -5.239576 2.386433e-07 1.026810e-03
## 13   219969_PM_at       TXLNG    X -5.227417 2.539965e-07 1.026810e-03
## 14   241320_PM_at             <NA> -5.209546 2.783132e-07 1.026810e-03
## 15   226219_PM_at    ARHGAP30    1  5.207319 2.814978e-07 1.026810e-03
## 16    40189_PM_at         SET    9  5.180730 3.223432e-07 1.039913e-03
## 17   229410_PM_at     SLC35E1   19  5.180268 3.231020e-07 1.039913e-03
## 18   211074_PM_at       FOLR1   11 -5.158299 3.612160e-07 1.051075e-03
## 19   215018_PM_at    KIAA1731   11 -5.145831 3.847447e-07 1.051075e-03
## 20 221705_PM_s_at       SIKE1    1 -5.125883 4.255087e-07 1.051075e-03
##      pos
## 1  18539
## 2  36536
## 3  41696
## 4  36255
## 5  11154
## 6  44895
## 7   8251
## 8  48879
## 9  54073
## 10 11475
## 11 44466
## 12 24110
## 13 29254
## 14 50570
## 15 35476
## 16 54349
## 17 38665
## 18 20439
## 19 24313
## 20 30986
## 
## pvalue quantiles for intercept and covariates>>
##         (Intercept)   asthmaTRUE   severeTRUE   smokerTRUE    sexmale
## min    0.000000e+00 1.719762e-09 7.107525e-10 2.604129e-14 0.00000000
## 25%    2.546901e-97 9.045413e-02 1.240844e-01 1.727951e-01 0.05351697
## median 5.360849e-67 3.275531e-01 3.767533e-01 4.250300e-01 0.27780405
## 75%    5.225016e-47 6.457142e-01 6.825190e-01 7.088523e-01 0.61493780
## max    9.861531e-01 9.999890e-01 9.999738e-01 9.999963e-01 0.99999435
##                gxPC1         gxPC2      ESTNEUT
## min    4.152111e-228 1.386093e-102 2.644747e-81
## 25%     1.116799e-32  4.692088e-13 9.342563e-07
## median  2.742412e-13  2.163824e-05 7.572454e-03
## 75%     2.202629e-04  4.130509e-02 2.173101e-01
## max     9.998273e-01  9.999473e-01 9.997582e-01
## 
## formula>>
## ~asthma + severe + smoker + sex + gxPC1 + gxPC2 + ESTNEUT
## 
## covariate of interest is  asthma 
## Number of tests= 54715 
## Number of arrays= 498 
## Number of significant tests (raw p-value <  0.05 )= 10128 
## Number of significant tests after p-value adjustments= 2305 
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
##           probeIDs geneSymbols  chr     stats         pval        p.adj
## 1     209124_PM_at       MYD88    3  6.167768 1.669291e-09 9.133527e-05
## 2   227280_PM_s_at      CCNYL1    2  5.876302 8.726027e-09 2.387223e-04
## 3     232441_PM_at        KRR1   12 -5.706957 2.215427e-08 3.656172e-04
## 4     235645_PM_at       ESCO1   18 -5.672353 2.672885e-08 3.656172e-04
## 5     226999_PM_at       RNPC3    1 -5.572903 4.560983e-08 4.464407e-04
## 6    1566342_PM_at                6  5.559620 4.895630e-08 4.464407e-04
## 7     201705_PM_at       PSMD7   16  5.472279 7.771758e-08 4.694030e-04
## 8     239629_PM_at       CFLAR    2 -5.459875 8.294961e-08 4.694030e-04
## 9     235216_PM_at       ESCO1   18 -5.456719 8.433485e-08 4.694030e-04
## 10    244822_PM_at               21 -5.453456 8.579054e-08 4.694030e-04
## 11    241320_PM_at             <NA> -5.393997 1.170126e-07 5.820312e-04
## 12 1564190_PM_x_at      ZNF519   18 -5.226719 2.760750e-07 9.397599e-04
## 13  229501_PM_s_at        USP8   15  5.224467 2.792424e-07 9.397599e-04
## 14   1552348_PM_at      PRSS33   16  5.222501 2.820364e-07 9.397599e-04
## 15    219969_PM_at       TXLNG    X -5.219820 2.858894e-07 9.397599e-04
## 16    211074_PM_at       FOLR1   11 -5.209339 3.014495e-07 9.397599e-04
## 17    226219_PM_at    ARHGAP30    1  5.196301 3.219525e-07 9.397599e-04
## 18    229410_PM_at     SLC35E1   19  5.195339 3.235188e-07 9.397599e-04
## 19    202026_PM_at        SDHD   11  5.192892 3.275325e-07 9.397599e-04
## 20   1570352_PM_at         ATM   11 -5.171333 3.650549e-07 9.397599e-04
##      pos
## 1  18539
## 2  36536
## 3  41696
## 4  44895
## 5  36255
## 6   8251
## 7  11154
## 8  48879
## 9  44466
## 10 54073
## 11 50570
## 12  7509
## 13 38756
## 14    75
## 15 29254
## 16 20439
## 17 35476
## 18 38665
## 19 11475
## 20  9801
## 
## pvalue quantiles for intercept and covariates>>
##         (Intercept)   asthmaTRUE   severeTRUE    sexmale         gxPC1
## min    0.000000e+00 1.669291e-09 5.373778e-10 0.00000000 7.963565e-190
## 25%    8.007602e-84 8.711253e-02 1.253311e-01 0.06792468  1.237271e-28
## median 7.430805e-58 3.223858e-01 3.764067e-01 0.30368424  8.133367e-12
## 75%    9.385885e-41 6.412002e-01 6.824862e-01 0.63487767  5.389431e-04
## max    9.867523e-01 9.999485e-01 9.999636e-01 0.99998131  9.993542e-01
##               gxPC2      ESTNEUT
## min    3.807316e-90 4.130943e-67
## 25%    1.641704e-11 5.134678e-06
## median 7.226723e-05 1.216234e-02
## 75%    5.610325e-02 2.481584e-01
## max    9.999893e-01 9.997975e-01
## 
## formula>>
## ~asthma + severe + sex + gxPC1 + gxPC2 + ESTNEUT
## 
## covariate of interest is  asthma 
## Number of tests= 54715 
## Number of arrays= 410 
## Number of significant tests (raw p-value <  0.05 )= 10377 
## Number of significant tests after p-value adjustments= 2479 
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
##          probeIDs            geneSymbols chr     stats         pval
## 1    214370_PM_at                 S100A8   1  6.245505 1.335494e-09
## 2    201094_PM_at                  RPS29  14  6.211049 1.624991e-09
## 3    224972_PM_at                  ROMO1  20  6.182344 1.912365e-09
## 4  202917_PM_s_at                 S100A8   1  6.104139 2.971706e-09
## 5    218856_PM_at               TNFRSF21   6 -6.050945 4.001110e-09
## 6  217249_PM_x_at                          4  6.041401 4.219546e-09
## 7  200834_PM_s_at LOC100291837 /// RPS21  20  6.033136 4.418113e-09
## 8    201597_PM_at                 COX7A2   6  5.973954 6.132549e-09
## 9    210244_PM_at                   CAMP   3  5.961576 6.565868e-09
## 10 216262_PM_s_at                  TGIF2  20 -5.951425 6.943411e-09
## 11   226773_PM_at                          4 -5.901115 9.150742e-09
## 12 201492_PM_s_at                  RPL41  12  5.826253 1.375377e-08
## 13 217773_PM_s_at                 NDUFA4   7  5.821030 1.414831e-08
## 14 224841_PM_x_at                   GAS5   1  5.775051 1.813293e-08
## 15 214334_PM_x_at                 DAZAP2  12 -5.749909 2.075503e-08
## 16 207573_PM_x_at                  ATP5L  11  5.724923 2.372564e-08
## 17 224741_PM_x_at                   GAS5   1  5.707949 2.597618e-08
## 18 229465_PM_s_at                         19 -5.699246 2.720954e-08
## 19 226845_PM_s_at                 MYEOV2   2  5.673586 3.118739e-08
## 20   221474_PM_at                 MYL12B  18  5.650673 3.521448e-08
##           p.adj   pos
## 1  3.453387e-05 23670
## 2  3.453387e-05 10543
## 3  3.453387e-05 34231
## 4  3.453387e-05 12367
## 5  3.453387e-05 28141
## 6  3.453387e-05 26535
## 7  3.453387e-05 10283
## 8  3.799088e-05 11046
## 9  3.799088e-05 19645
## 10 3.799088e-05 25555
## 11 4.551662e-05 36029
## 12 5.954806e-05 10941
## 13 5.954806e-05 27059
## 14 7.086739e-05 34100
## 15 7.570744e-05 23634
## 16 8.113429e-05 17015
## 17 8.270944e-05 34001
## 18 8.270944e-05 38720
## 19 8.981149e-05 36101
## 20 9.633801e-05 30758
## 
## pvalue quantiles for intercept and covariates>>
##         (Intercept)   severeTRUE    sexmale         gxPC1        gxPC2
## min    0.000000e+00 1.335494e-09 0.00000000 3.521421e-149 8.929113e-70
## 25%    2.386617e-62 1.258925e-01 0.09499583  3.318426e-22 4.402135e-09
## median 1.825830e-42 3.809366e-01 0.34563960  3.140942e-09 4.550965e-04
## 75%    1.519464e-29 6.851606e-01 0.65991725  3.350454e-03 8.994666e-02
## max    9.856730e-01 9.999947e-01 0.99993437  9.999073e-01 9.999171e-01
##             ESTNEUT
## min    2.212534e-49
## 25%    4.328954e-05
## median 2.203923e-02
## 75%    2.904202e-01
## max    9.999777e-01
## 
## formula>>
## ~severe + sex + gxPC1 + gxPC2 + ESTNEUT
## 
## covariate of interest is  severe 
## Number of tests= 54715 
## Number of arrays= 323 
## Number of significant tests (raw p-value <  0.05 )= 8178 
## Number of significant tests after p-value adjustments= 1447 
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
## [1] 2305
```

```r
dge2 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # cases vs. ctrls. (no smokers)
```

```
## [1] 2479
```

```r
dge3 %>% use_series(frame) %>% filter(p.adj < 0.05) %>% nrow # severe cases vs. moderate cases (no smokers)
```

```
## [1] 1447
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
##   [5] "ADRB2"        "AKIRIN2"      "ANKRD57"      "AQP9"        
##   [9] "ARAP3"        "ARID5A"       "ARL5B"        "ARRB1"       
##  [13] "ATF3"         "BCAR3"        "BID"          "C10orf128"   
##  [17] "C11orf75"     "C14orf102"    "C17orf79"     "C1orf122"    
##  [21] "C1orf38"      "C6orf1"       "CALHM2"       "CASP5"       
##  [25] "CBL"          "CCL23"        "CCNG2"        "CCR2"        
##  [29] "CD300LB"      "CD68"         "CDKN1A"       "CEP350"      
##  [33] "CHST15"       "CLCF1"        "CLIC4"        "CMTM3"       
##  [37] "CRBN"         "CSNK1D"       "CTSB"         "CYBB"        
##  [41] "CYTH4"        "DAB2"         "DCP2"         "DDX17"       
##  [45] "DOK2"         "DR1"          "DRAM1"        "DUSP18"      
##  [49] "DUSP3"        "DYRK3"        "EFR3A"        "EHD1"        
##  [53] "ELF2"         "ENG"          "EPB41L3"      "ETS2"        
##  [57] "EVI2B"        "FAM105A"      "FKBP15"       "FOXN2"       
##  [61] "FRAT1"        "FUCA1"        "GAPT"         "GBP1"        
##  [65] "GDI2"         "GFOD1"        "GK"           "GLA"         
##  [69] "GRB2"         "GRN"          "HBEGF"        "HCK"         
##  [73] "HECA"         "HHEX"         "HSD3B7"       "IER5"        
##  [77] "IFI30"        "IFIH1"        "IFNGR2"       "IL18"        
##  [81] "IL1B"         "ITGAX"        "JMJD1C"       "KCNJ2"       
##  [85] "KLF6"         "KMO"          "LAIR1"        "LDLR"        
##  [89] "LFNG"         "LOC100129034" "LRRC25"       "LRRC33"      
##  [93] "LYSMD2"       "MAP3K14"      "MAP3K8"       "MAPK13"      
##  [97] "MAPKAPK3"     "MARCKS"       "MBD2"         "MKL1"        
## [101] "MR1"          "MXD1"         "MYD88"        "N4BP1"       
## [105] "NACC2"        "NCKAP1L"      "NCOR2"        "NFATC1"      
## [109] "NFIC"         "NFKBIZ"       "NIN"          "NINJ1"       
## [113] "NPEPPS"       "OASL"         "OLIG1"        "OLIG2"       
## [117] "PARP10"       "PDK4"         "PELI1"        "PHACTR2"     
## [121] "PIK3AP1"      "PILRA"        "PLAUR"        "PLEK"        
## [125] "PLSCR1"       "PLXNB2"       "PMAIP1"       "PPP2R5A"     
## [129] "PSTPIP2"      "PTAFR"        "PTGER4"       "PTGS2"       
## [133] "PTPN6"        "RAB11FIP4"    "RAB31"        "RAPGEF2"     
## [137] "RARA"         "RASSF2"       "RCOR3"        "RGS1"        
## [141] "RGS12"        "RGS19"        "RHOB"         "RHOBTB3"     
## [145] "RHOH"         "RHOU"         "RIN2"         "RIN3"        
## [149] "RIT1"         "RNF144B"      "RNF19B"       "RRAGC"       
## [153] "RRAGD"        "RYBP"         "S1PR3"        "SCARB2"      
## [157] "SCO2"         "SERPINB1"     "SERPINB2"     "SESTD1"      
## [161] "SH2D3C"       "SIRPA"        "SLA"          "SLC16A3"     
## [165] "SLC2A3"       "SLC2A6"       "SLC38A2"      "SLC43A2"     
## [169] "SNX20"        "SOD2"         "SORL1"        "SPAG9"       
## [173] "ST3GAL6"      "ST8SIA4"      "STK17B"       "STX11"       
## [177] "STX3"         "SYN2"         "TAF5"         "TGFBR1"      
## [181] "THBS1"        "TIFA"         "TIMP2"        "TLR1"        
## [185] "TLR4"         "TLR7"         "TMBIM1"       "TMEM123"     
## [189] "TMEM88"       "TNF"          "TNFAIP2"      "TNFSF13B"    
## [193] "TNRC6A"       "TOR1B"        "TRIB1"        "TRIM33"      
## [197] "TRIM38"       "TSC22D3"      "UBE2W"        "USP22"       
## [201] "VAV1"         "WDFY1"        "ZBTB2"        "ZC3HAV1"     
## [205] "ZFP36"        "ZFYVE16"      "ZNF227"       "ZNF398"
```

```r
core_genes %>% length
```

```
## [1] 208
```

```r
# What is the overlap of this core with highlighted TREM1 genes from Table 3 of Croteau-Chonka et al. (2017)?
table3_genes <- c("CCL23", "OLIG1", "OLIG2", "GFOD1", "RHOBTB3", "HSD3B7")
intersect(core_genes, table3_genes)
```

```
## [1] "CCL23"   "GFOD1"   "HSD3B7"  "OLIG1"   "OLIG2"   "RHOBTB3"
```

```r
dge1 %>% use_series(frame) %>% filter(geneSymbols %in% table3_genes)
```

```
##          probeIDs geneSymbols chr      stats         pval       p.adj
## 1  210549_PM_s_at       CCL23  17  4.2643240 2.402434e-05 0.005134733
## 2    213825_PM_at       OLIG2  21  4.1486384 3.936892e-05 0.006483616
## 3    210548_PM_at       CCL23  17  4.1070741 4.688221e-05 0.007001643
## 4    228170_PM_at       OLIG1  21  2.7698083 5.819737e-03 0.085571990
## 5    222817_PM_at      HSD3B7  16  2.7328801 6.503029e-03 0.090838205
## 6  219821_PM_s_at       GFOD1   6  2.1847211 2.937779e-02 0.204857411
## 7  202976_PM_s_at     RHOBTB3   5  2.1198299 3.451761e-02 0.222545129
## 8    225202_PM_at     RHOBTB3   5  1.9700836 4.938602e-02 0.268121678
## 9    216049_PM_at     RHOBTB3   5 -1.5757498 1.157226e-01 0.407934595
## 10   213824_PM_at       OLIG2  21  1.3853653 1.665650e-01 0.483765151
## 11 216048_PM_s_at     RHOBTB3   5  0.3561856 7.218534e-01 0.895618405
## 12   240111_PM_at     RHOBTB3   5 -0.2917926 7.705675e-01 0.915899020
## 13 202975_PM_s_at     RHOBTB3   5 -0.1752128 8.609841e-01 0.951746043
##      pos
## 1  19938
## 2  23126
## 3  19937
## 4  37425
## 5  32097
## 6  29106
## 7  12426
## 8  34460
## 9  25342
## 10 23125
## 11 25341
## 12 49361
## 13 12425
```

```r
# Which primary analysis gene sets were enriched in secondary analyses?
fgseaRes1 %>% arrange(desc(NES)) %>% select(-leadingEdge) %>% filter(padj < 0.05)
```

```
##                                                       pathway         pval
## 1               GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 2.143577e-05
## 2     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 2.127931e-05
## 3             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 2.133333e-05
## 4                 GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 2.133333e-05
## 5           GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 2.145278e-05
## 6              GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.134381e-05
## 7  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 2.133333e-05
## 8                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.138397e-05
## 9                       GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 2.142704e-05
## 10                  GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 2.154708e-05
## 11         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 2.123999e-05
## 12                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 2.154708e-05
## 13                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.135201e-05
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.135201e-05
## 15              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 2.133333e-05
## 16       GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.144082e-05
## 17            GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 6.400000e-05
## 18                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.495886e-04
## 19    GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 3.632634e-04
## 20                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 8.115496e-04
## 21             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 1.067395e-03
## 22                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.856000e-03
## 23         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 3.822010e-03
## 24 GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 5.796424e-03
##            padj        ES      NES nMoreExtreme size
## 1  3.232062e-05 0.7462695 3.037452            0   59
## 2  3.232062e-05 0.7706867 2.892690            0   41
## 3  3.232062e-05 0.7525535 2.869805            0   44
## 4  3.232062e-05 0.7380543 2.814513            0   44
## 5  3.232062e-05 0.6914557 2.771448            0   55
## 6  3.232062e-05 0.7025062 2.664984            0   43
## 7  3.232062e-05 0.6843808 2.609833            0   44
## 8  3.232062e-05 0.6660411 2.590442            0   48
## 9  3.232062e-05 0.6458536 2.557266            0   52
## 10 3.232062e-05 0.6216615 2.553596            0   62
## 11 3.232062e-05 0.6833943 2.551336            0   40
## 12 3.232062e-05 0.6160032 2.530353            0   62
## 13 3.232062e-05 0.6679562 2.519729            0   42
## 14 3.232062e-05 0.6294299 2.374397            0   42
## 15 3.232062e-05 0.6052969 2.308253            0   44
## 16 3.232062e-05 0.5334658 2.162636            0   58
## 17 9.035294e-05 0.5517037 2.103880            2   44
## 18 1.994515e-04 0.5247718 2.060198            6   50
## 19 4.588590e-04 0.4925250 1.925896           16   49
## 20 9.738596e-04 0.4833068 1.861876           37   46
## 21 1.219880e-03 0.4796964 1.839242           49   45
## 22 2.024727e-03 0.4715886 1.798367           86   44
## 23 3.988184e-03 0.4599951 1.735238          178   42
## 24 5.796424e-03 0.4582767 1.691799          272   38
```

```r
fgseaRes2 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05)
```

```
##                                                       pathway         pval
## 1     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_DN 2.136250e-05
## 2     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 2.361579e-04
## 3  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_DN 2.134472e-05
## 4  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 2.960218e-03
## 5        GSE9988_ANTI_TREM1_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.159081e-05
## 6             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 2.134472e-05
## 7             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 8.537887e-05
## 8                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN 2.164783e-05
## 9                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.138626e-05
## 10                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN 2.164783e-05
## 11                      GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 2.151556e-05
## 12         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_DN 2.131469e-05
## 13         GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 2.408817e-03
## 14          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 2.131696e-05
## 15          GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_UP 2.156892e-05
## 16                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_DN 2.134472e-05
## 17                GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.664888e-03
## 18             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.130379e-05
## 19             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 8.330663e-04
## 20                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 4.287245e-05
## 21              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 2.134472e-05
## 22              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_UP 2.164081e-05
## 23                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_DN 2.131696e-05
## 24                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 4.277617e-04
##            padj        ES      NES nMoreExtreme size
## 1  3.247175e-05 0.7779490 2.925557            0   41
## 2  2.983047e-04 0.5003731 1.957904           10   49
## 3  3.247175e-05 0.6895618 2.635430            0   44
## 4  2.960218e-03 0.4743669 1.752870          138   38
## 5  3.247175e-05 0.5455014 2.211442            0   58
## 6  3.247175e-05 0.7461772 2.851808            0   44
## 7  1.138385e-04 0.5545600 2.119468            3   44
## 8  3.247175e-05 0.6244654 2.566773            0   62
## 9  3.247175e-05 0.6680436 2.603946            0   48
## 10 3.247175e-05 0.6184207 2.541927            0   62
## 11 3.247175e-05 0.6462016 2.559111            0   52
## 12 3.247175e-05 0.6911018 2.584515            0   40
## 13 2.513548e-03 0.4693161 1.775373          112   42
## 14 3.247175e-05 0.6235836 2.358951            0   42
## 15 3.247175e-05 0.6955706 2.788095            0   55
## 16 3.247175e-05 0.7472163 2.855779            0   44
## 17 1.816241e-03 0.4748922 1.814986           77   44
## 18 3.247175e-05 0.7127222 2.710032            0   43
## 19 9.520758e-04 0.4854217 1.864894           38   45
## 20 6.052582e-05 0.5356944 2.105468            1   50
## 21 3.247175e-05 0.6012046 2.297739            0   44
## 22 3.247175e-05 0.7475125 3.039678            0   59
## 23 3.247175e-05 0.6841896 2.588217            0   42
## 24 5.133141e-04 0.4960656 1.914845           19   46
```

```r
fgseaRes3 %>% filter(pathway %in% trem1_pathways) %>% select(-leadingEdge) %>% filter(padj < 0.05)
```

```
##                                                       pathway         pval
## 1     GSE9988_ANTI_TREM1_AND_LPS_VS_CTRL_TREATED_MONOCYTES_UP 1.290123e-04
## 2  GSE9988_ANTI_TREM1_AND_LPS_VS_VEHICLE_TREATED_MONOCYTES_UP 2.508466e-05
## 3             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_DN 5.755616e-03
## 4             GSE9988_ANTI_TREM1_VS_CTRL_TREATED_MONOCYTES_UP 1.375236e-03
## 5                   GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_UP 2.405578e-02
## 6                       GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_UP 1.684282e-02
## 7          GSE9988_ANTI_TREM1_VS_VEHICLE_TREATED_MONOCYTES_UP 5.062522e-05
## 8           GSE9988_LOW_LPS_VS_ANTI_TREM1_AND_LPS_MONOCYTE_DN 1.769351e-02
## 9                 GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP 1.298834e-02
## 10             GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 2.814307e-03
## 11                    GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP 8.170653e-03
## 12              GSE9988_LPS_VS_LPS_AND_ANTI_TREM1_MONOCYTE_DN 5.857485e-04
## 13                 GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP 4.254556e-03
##            padj        ES      NES nMoreExtreme size
## 1  0.0010320983 0.5160203 2.029216            4   49
## 2  0.0006020319 0.5956238 2.211400            0   38
## 3  0.0172668466 0.4357837 1.671790          225   44
## 4  0.0066011307 0.4735213 1.816563           53   44
## 5  0.0444106681 0.3813409 1.492764          934   48
## 6  0.0353870298 0.3835667 1.528109          646   52
## 7  0.0006075027 0.5565813 2.113514            1   42
## 8  0.0353870298 0.4083608 1.550674          698   42
## 9  0.0311720063 0.4104216 1.574494          509   44
## 10 0.0112572277 0.4505652 1.736675          109   45
## 11 0.0217884077 0.4098414 1.619508          315   50
## 12 0.0035144909 0.4990881 1.914644           22   44
## 13 0.0145870481 0.4385599 1.698821          165   46
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
##  ui       RStudio (1.0.143)           
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       <NA>                        
##  date     2017-09-26
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version  date       source                      
##  affy                   1.54.0   2017-05-02 Bioconductor                
##  affyio                 1.46.0   2017-05-02 Bioconductor                
##  annotate             * 1.54.0   2017-05-02 Bioconductor                
##  AnnotationDbi        * 1.38.2   2017-08-16 Bioconductor                
##  assertthat             0.2.0    2017-04-11 CRAN (R 3.4.0)              
##  backports              1.1.0    2017-05-22 CRAN (R 3.4.0)              
##  base                 * 3.4.0    2017-05-02 local                       
##  base64                 2.0      2016-05-10 CRAN (R 3.4.0)              
##  beanplot               1.2      2014-09-19 CRAN (R 3.4.0)              
##  beeswarm               0.2.3    2016-04-25 CRAN (R 3.4.0)              
##  bibtex                 0.4.2    2017-06-30 CRAN (R 3.4.0)              
##  bigmemory            * 4.5.19   2016-03-28 CRAN (R 3.4.0)              
##  bigmemory.sri        * 0.1.3    2014-08-18 CRAN (R 3.4.0)              
##  bindr                  0.1      2016-11-13 CRAN (R 3.4.0)              
##  bindrcpp             * 0.2      2017-06-17 CRAN (R 3.4.0)              
##  Biobase              * 2.36.2   2017-07-19 Bioconductor                
##  BiocGenerics         * 0.22.0   2017-05-02 Bioconductor                
##  BiocInstaller          1.26.0   2017-05-02 Bioconductor                
##  BiocParallel           1.10.1   2017-07-19 Bioconductor                
##  biomaRt                2.32.1   2017-07-19 Bioconductor                
##  Biostrings             2.44.2   2017-07-25 Bioconductor                
##  bit                    1.1-12   2014-04-09 CRAN (R 3.4.0)              
##  bit64                  0.9-7    2017-05-08 CRAN (R 3.4.0)              
##  bitops                 1.0-6    2013-08-17 CRAN (R 3.4.0)              
##  blob                   1.1.0    2017-06-17 CRAN (R 3.4.0)              
##  broom                * 0.4.2    2017-02-13 CRAN (R 3.4.0)              
##  bumphunter             1.16.0   2017-05-02 Bioconductor                
##  caTools                1.17.1   2014-09-10 CRAN (R 3.4.0)              
##  CellMix              * 1.6.2    2017-07-25 local                       
##  cluster              * 2.0.6    2017-03-16 CRAN (R 3.4.0)              
##  codetools              0.2-15   2016-10-05 CRAN (R 3.4.0)              
##  colorspace             1.3-2    2016-12-14 CRAN (R 3.4.0)              
##  compiler             * 3.4.0    2017-05-02 local                       
##  corpcor                1.6.9    2017-04-01 CRAN (R 3.4.0)              
##  csSAM                * 1.2.4    2013-05-13 CRAN (R 3.4.0)              
##  data.table             1.10.4   2017-02-01 CRAN (R 3.4.0)              
##  datasets             * 3.4.0    2017-05-02 local                       
##  DBI                    0.7      2017-06-18 CRAN (R 3.4.0)              
##  DelayedArray           0.2.7    2017-07-19 Bioconductor                
##  devtools             * 1.13.3   2017-08-02 CRAN (R 3.4.0)              
##  digest                 0.6.12   2017-01-27 CRAN (R 3.4.0)              
##  doParallel             1.0.10   2015-10-14 CRAN (R 3.4.0)              
##  doRNG                  1.6.6    2017-04-10 CRAN (R 3.4.0)              
##  dplyr                * 0.7.2    2017-07-20 CRAN (R 3.4.0)              
##  evaluate               0.10.1   2017-06-24 CRAN (R 3.4.0)              
##  ezknitr              * 0.6      2016-09-16 CRAN (R 3.4.0)              
##  fastmatch              1.1-0    2017-01-28 CRAN (R 3.4.0)              
##  fgsea                * 1.2.1    2017-07-18 Bioconductor                
##  forcats              * 0.2.0    2017-01-23 CRAN (R 3.4.0)              
##  foreach                1.4.3    2015-10-13 CRAN (R 3.4.0)              
##  foreign                0.8-69   2017-06-21 CRAN (R 3.4.0)              
##  gdata                  2.18.0   2017-06-06 CRAN (R 3.4.0)              
##  genefilter           * 1.58.1   2017-07-19 Bioconductor                
##  GeneSelectMMD          2.20.1   2017-07-07 Bioconductor                
##  GenomeInfoDb           1.12.2   2017-07-19 Bioconductor                
##  GenomeInfoDbData       0.99.0   2017-05-02 Bioconductor                
##  GenomicAlignments      1.12.1   2017-07-19 Bioconductor                
##  GenomicFeatures        1.28.4   2017-07-19 Bioconductor                
##  GenomicRanges          1.28.4   2017-07-19 Bioconductor                
##  GEOquery               2.42.0   2017-05-02 Bioconductor                
##  ggplot2                2.2.1    2016-12-30 CRAN (R 3.4.0)              
##  glue                   1.1.1    2017-06-21 CRAN (R 3.4.0)              
##  gplots               * 3.0.1    2016-03-30 CRAN (R 3.4.0)              
##  graph                * 1.54.0   2017-05-02 Bioconductor                
##  graphics             * 3.4.0    2017-05-02 local                       
##  grDevices            * 3.4.0    2017-05-02 local                       
##  grid                   3.4.0    2017-05-02 local                       
##  gridBase               0.4-7    2014-02-24 CRAN (R 3.4.0)              
##  gridExtra              2.2.1    2016-02-29 CRAN (R 3.4.0)              
##  GSEABase             * 1.38.0   2017-07-18 Bioconductor                
##  gtable                 0.2.0    2016-02-26 CRAN (R 3.4.0)              
##  gtools                 3.5.0    2015-05-29 CRAN (R 3.4.0)              
##  here                 * 0.1      2017-09-26 Github (krlmlr/here@93593ee)
##  hgu133a.db           * 3.2.3    2017-07-18 Bioconductor                
##  hgu133b.db           * 3.2.3    2017-07-25 Bioconductor                
##  hthgu133pluspm.db    * 0.1      2017-07-20 Bioconductor                
##  htmltools              0.3.6    2017-04-28 CRAN (R 3.4.0)              
##  htmlwidgets            0.9      2017-07-10 CRAN (R 3.4.0)              
##  httpuv                 1.3.5    2017-07-04 CRAN (R 3.4.0)              
##  httr                   1.3.0    2017-08-16 CRAN (R 3.4.0)              
##  iCheck               * 1.6.0    2017-07-18 Bioconductor                
##  illuminaio             0.18.0   2017-05-02 Bioconductor                
##  IRanges              * 2.10.2   2017-07-19 Bioconductor                
##  iterators              1.0.8    2015-10-13 CRAN (R 3.4.0)              
##  jsonlite               1.5      2017-06-01 CRAN (R 3.4.0)              
##  KernSmooth             2.23-15  2015-06-29 CRAN (R 3.4.0)              
##  knitr                  1.17     2017-08-10 CRAN (R 3.4.0)              
##  labeling               0.3      2014-08-23 CRAN (R 3.4.0)              
##  lattice                0.20-35  2017-03-25 CRAN (R 3.4.0)              
##  lazyeval               0.2.0    2016-06-12 CRAN (R 3.4.0)              
##  limma                  3.32.5   2017-08-16 Bioconductor                
##  limSolve               1.5.5.3  2017-08-14 CRAN (R 3.4.0)              
##  lmtest                 0.9-35   2017-02-11 CRAN (R 3.4.0)              
##  locfit                 1.5-9.1  2013-04-20 CRAN (R 3.4.0)              
##  lpSolve                5.6.13   2015-09-19 CRAN (R 3.4.0)              
##  lumi                 * 2.28.0   2017-05-02 Bioconductor                
##  magrittr             * 1.5      2014-11-22 CRAN (R 3.4.0)              
##  MASS                   7.3-47   2017-04-21 CRAN (R 3.4.0)              
##  Matrix                 1.2-11   2017-08-16 CRAN (R 3.4.0)              
##  matrixStats            0.52.2   2017-04-14 CRAN (R 3.4.0)              
##  mclust                 5.3      2017-05-21 CRAN (R 3.4.0)              
##  memoise                1.1.0    2017-04-21 CRAN (R 3.4.0)              
##  methods              * 3.4.0    2017-05-02 local                       
##  methylumi              2.22.0   2017-05-02 Bioconductor                
##  mgcv                   1.8-18   2017-07-28 CRAN (R 3.4.0)              
##  mime                   0.5      2016-07-07 CRAN (R 3.4.0)              
##  minfi                  1.22.1   2017-07-19 Bioconductor                
##  mnormt                 1.5-5    2016-10-15 CRAN (R 3.4.0)              
##  multtest               2.32.0   2017-05-03 Bioconductor                
##  munsell                0.4.3    2016-02-13 CRAN (R 3.4.0)              
##  nleqslv                3.3.1    2017-07-06 CRAN (R 3.4.0)              
##  nlme                   3.1-131  2017-02-06 CRAN (R 3.4.0)              
##  NMF                  * 0.20.6   2015-05-26 CRAN (R 3.4.0)              
##  nor1mix                1.2-2    2016-08-25 CRAN (R 3.4.0)              
##  openssl                0.9.6    2016-12-31 CRAN (R 3.4.0)              
##  org.Hs.eg.db         * 3.4.1    2017-05-03 Bioconductor                
##  parallel             * 3.4.0    2017-05-02 local                       
##  pkgconfig              2.0.1    2017-03-21 CRAN (R 3.4.0)              
##  pkgmaker             * 0.22     2014-05-14 CRAN (R 3.4.0)              
##  plyr                   1.8.4    2016-06-08 CRAN (R 3.4.0)              
##  preprocessCore         1.38.1   2017-07-19 Bioconductor                
##  psych                  1.7.5    2017-05-03 CRAN (R 3.4.0)              
##  purrr                  0.2.3    2017-08-02 CRAN (R 3.4.0)              
##  quadprog               1.5-5    2013-04-17 CRAN (R 3.4.0)              
##  R.methodsS3            1.7.1    2016-02-16 CRAN (R 3.4.0)              
##  R.oo                   1.21.0   2016-11-01 CRAN (R 3.4.0)              
##  R.utils                2.5.0    2016-11-07 CRAN (R 3.4.0)              
##  R6                     2.2.2    2017-06-17 CRAN (R 3.4.0)              
##  randomForest           4.6-12   2015-10-07 CRAN (R 3.4.0)              
##  RColorBrewer           1.1-2    2014-12-07 CRAN (R 3.4.0)              
##  Rcpp                 * 0.12.12  2017-07-15 CRAN (R 3.4.0)              
##  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.4.0)              
##  registry             * 0.3      2015-07-08 CRAN (R 3.4.0)              
##  reshape                0.8.7    2017-08-06 CRAN (R 3.4.0)              
##  reshape2               1.4.2    2016-10-22 CRAN (R 3.4.0)              
##  rgl                    0.98.1   2017-03-08 CRAN (R 3.4.0)              
##  rlang                  0.1.2    2017-08-09 CRAN (R 3.4.0)              
##  rngtools             * 1.2.4    2014-03-06 CRAN (R 3.4.0)              
##  rprojroot              1.2      2017-01-16 CRAN (R 3.4.0)              
##  Rsamtools              1.28.0   2017-05-03 Bioconductor                
##  RSQLite                2.0      2017-06-19 CRAN (R 3.4.0)              
##  rtracklayer            1.36.4   2017-07-19 Bioconductor                
##  S4Vectors            * 0.14.3   2017-07-19 Bioconductor                
##  scales                 0.4.1    2016-11-09 CRAN (R 3.4.0)              
##  scatterplot3d          0.3-40   2017-04-22 CRAN (R 3.4.0)              
##  shiny                  1.0.4    2017-08-14 CRAN (R 3.4.0)              
##  siggenes               1.50.0   2017-05-03 Bioconductor                
##  splines                3.4.0    2017-05-02 local                       
##  stats                * 3.4.0    2017-05-02 local                       
##  stats4               * 3.4.0    2017-05-02 local                       
##  stringi                1.1.5    2017-04-07 CRAN (R 3.4.0)              
##  stringr              * 1.2.0    2017-02-18 CRAN (R 3.4.0)              
##  SummarizedExperiment   1.6.3    2017-07-19 Bioconductor                
##  survival               2.41-3   2017-04-04 CRAN (R 3.4.0)              
##  tibble                 1.3.3    2017-05-28 CRAN (R 3.4.0)              
##  tidyr                  0.7.0    2017-08-16 CRAN (R 3.4.0)              
##  tools                  3.4.0    2017-05-02 local                       
##  utils                * 3.4.0    2017-05-02 local                       
##  withr                  2.0.0    2017-07-28 CRAN (R 3.4.0)              
##  XML                  * 3.98-1.9 2017-06-19 CRAN (R 3.4.0)              
##  xtable                 1.8-2    2016-02-05 CRAN (R 3.4.0)              
##  XVector                0.16.0   2017-05-03 Bioconductor                
##  zlibbioc               1.22.0   2017-05-03 Bioconductor                
##  zoo                    1.8-0    2017-04-12 CRAN (R 3.4.0)
```

