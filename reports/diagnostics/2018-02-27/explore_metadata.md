

```r
library("here")

source(here("code", "data", "create_ubiopred_eset.R"), echo = TRUE)
```

```
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
##  33.228   0.173  33.479 
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
```

```r
library("reshape2")
library("ggplot2")
library("corrr")
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ cohort) %>% tidy
```

```
##                                 term  estimate std.error  statistic
## 1                        (Intercept) -8.794578  10.17976 -0.8639276
## 2 cohortModerate asthma, non-smoking -6.319154  14.85642 -0.4253483
## 3   cohortSevere asthma, non-smoking 21.685718  11.84383  1.8309716
## 4       cohortSevere asthma, smoking -5.322863  14.35540 -0.3707917
##      p.value
## 1 0.38804693
## 2 0.67076798
## 3 0.06770695
## 4 0.71095166
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ cohort) %>% tidy
```

```
##                                 term   estimate std.error statistic
## 1                        (Intercept)  27.877999  7.139687  3.904653
## 2 cohortModerate asthma, non-smoking  -6.249118 10.419713 -0.599740
## 3   cohortSevere asthma, non-smoking -42.494174  8.306800 -5.115589
## 4       cohortSevere asthma, smoking -33.505621 10.068317 -3.327828
##        p.value
## 1 1.074759e-04
## 2 5.489545e-01
## 3 4.483931e-07
## 4 9.405590e-04
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ sex) %>% tidy
```

```
##          term   estimate std.error  statistic   p.value
## 1 (Intercept) -0.5304473  5.766898 -0.0919814 0.9267500
## 2     sexmale  1.1845864  8.617960  0.1374555 0.8907265
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ sex) %>% tidy
```

```
##          term  estimate std.error statistic   p.value
## 1 (Intercept) -4.308919  4.140480 -1.040681 0.2985306
## 2     sexmale  9.622608  6.187468  1.555177 0.1205417
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ race) %>% tidy
```

```
##                        term     estimate std.error   statistic   p.value
## 1               (Intercept)  -0.08134523  4.516159 -0.01801204 0.9856366
## 2      racesouth_east_asian -61.39679027 55.556710 -1.10511926 0.2696514
## 3                 raceother   0.62881331 43.128727  0.01457992 0.9883733
## 4         raceblack_african  -5.54472989 25.171932 -0.22027431 0.8257495
## 5           racesouth_asian  33.94130272 30.663356  1.10690112 0.2688809
## 6 racearabic_north_heritage -44.48699142 43.128727 -1.03149327 0.3028195
## 7        racemultiple_races  16.35469014 36.530275  0.44770236 0.6545663
## 8         racecentral_asian  54.74878765 96.014856  0.57021163 0.5687961
## 9            raceeast_asian  18.51748926 96.014856  0.19286067 0.8471482
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ race) %>% tidy
```

```
##                        term  estimate std.error  statistic   p.value
## 1               (Intercept)  -1.32917  3.244163 -0.4097114 0.6821971
## 2      racesouth_east_asian  36.31465 39.908915  0.9099382 0.3633034
## 3                 raceother  30.52513 30.981328  0.9852751 0.3249762
## 4         raceblack_african  20.92979 18.082145  1.1574835 0.2476399
## 5           racesouth_asian -17.36906 22.026885 -0.7885392 0.4307635
## 6 racearabic_north_heritage  10.07185 30.981328  0.3250943 0.7452487
## 7        racemultiple_races  22.59660 26.241360  0.8611063 0.3896015
## 8         racecentral_asian  80.41403 68.971843  1.1658964 0.2442247
## 9            raceeast_asian -28.84846 68.971843 -0.4182644 0.6759375
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC1 ~ site) %>% tidy
```

```
##           term   estimate std.error  statistic      p.value
## 1  (Intercept)  24.095129  18.84276  1.2787473 2.016001e-01
## 2        siteb  -3.041407  26.64769 -0.1141340 9.091790e-01
## 3        sitec -30.277732  20.99728 -1.4419835 1.499551e-01
## 4        sited -57.267814  23.16685 -2.4719725 1.378055e-02
## 5        sitee  67.697023  25.38404  2.6669133 7.912062e-03
## 6        sitef -29.048843  25.38404 -1.1443745 2.530352e-01
## 7        siteg -91.723418  23.16685 -3.9592530 8.649424e-05
## 8        siteh  15.935292  20.97505  0.7597260 4.477890e-01
## 9        sitei -15.547105  27.30575 -0.5693710 5.693691e-01
## 10       sitej -40.281985  23.91703 -1.6842385 9.278163e-02
## 11       sitek -71.956410  24.32590 -2.9580165 3.248088e-03
## 12       sitel -38.856926  31.71709 -1.2251100 2.211309e-01
## 13       sitem -48.312531  31.71709 -1.5232333 1.283548e-01
## 14       siten -26.372244  31.71709 -0.8314837 4.061113e-01
## 15       siteo  36.340198  40.70501  0.8927696 3.724252e-01
```

```r
es_ubiopred_wb %>% pData %$% lm(gxPC2 ~ site) %>% tidy
```

```
##           term   estimate std.error   statistic      p.value
## 1  (Intercept) -14.022701  13.44642 -1.04285752 2.975360e-01
## 2        siteb  25.661923  19.01611  1.34948325 1.778140e-01
## 3        sitec  21.384331  14.98391  1.42715268 1.541820e-01
## 4        sited  -9.627869  16.53214 -0.58237269 5.605874e-01
## 5        sitee  78.066975  18.11435  4.30967466 1.982233e-05
## 6        sitef -18.864004  18.11435 -1.04138426 2.982182e-01
## 7        siteg -26.046820  16.53214 -1.57552593 1.157898e-01
## 8        siteh  21.757178  14.96805  1.45357469 1.467137e-01
## 9        sitei  48.078455  19.48572  2.46736923 1.395672e-02
## 10       sitej  31.355313  17.06748  1.83713761 6.680378e-02
## 11       sitek -36.699083  17.35925 -2.11409322 3.501920e-02
## 12       sitel  41.953042  22.63370  1.85356538 6.441091e-02
## 13       sitem   1.737400  22.63370  0.07676164 9.388450e-01
## 14       siten  34.400104  22.63370  1.51986219 1.291999e-01
## 15       siteo  98.377110  29.04759  3.38675641 7.649783e-04
```

```r
# pdf(file = here("reports", "diagnostics", "gxPCA.pdf"), width = 8, height = 8, pointsize = 8)

screeplot(res_pca_ubiopred_wb$pcs, type = "lines", main = "U-BIOPRED WB gxPCs")
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-1.png)

```r
summary(res_pca_ubiopred_wb$pcs) %>% use_series(importance) %>% extract(, 1:10)
```

```
##                             PC1      PC2      PC3      PC4      PC5
## Standard deviation     95.53874 68.76001 52.56790 47.00382 39.50203
## Proportion of Variance  0.16682  0.08641  0.05051  0.04038  0.02852
## Cumulative Proportion   0.16682  0.25323  0.30374  0.34412  0.37264
##                             PC6      PC7      PC8      PC9     PC10
## Standard deviation     34.03766 31.47322 28.80299 27.00969 25.41572
## Proportion of Variance  0.02117  0.01810  0.01516  0.01333  0.01181
## Cumulative Proportion   0.39381  0.41191  0.42708  0.44041  0.45222
```

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "cohort",
          title = "U-BIOPRED WB [Cohort]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-2.png)

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(3, 2),
          labelVariable = "cohort",
          title = "U-BIOPRED WB [Cohort]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC3",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-3.png)

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(3, 4),
          labelVariable = "cohort",
          title = "U-BIOPRED WB [Cohort]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC3",
          ylab = "PC4",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-4.png)

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(5, 4),
          labelVariable = "cohort",
          title = "U-BIOPRED WB [Cohort]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC5",
          ylab = "PC4",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-5.png)

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "sex",
          title = "U-BIOPRED WB [Sex]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-6.png)

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "race",
          title = "U-BIOPRED WB [Race]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-7.png)

```r
pca2DPlot(pcaObj = res_pca_ubiopred_wb,
          plot.dim = c(1, 2),
          labelVariable = "site",
          title = "U-BIOPRED WB [Site]",
          plotOutPutFlag = FALSE,
          equalRange = FALSE,
          xlab = "PC1",
          ylab = "PC2",
          cex.legend = 0.6,
          cex = 1,
          cex.lab = 1,
          cex.axis = 1,
          legendPosition = "topleft"
)
```

![plot of chunk explore_gx_pca](explore_metadata//explore_gx_pca-8.png)

```r
# dev.off()
```

```r
# pdf(file = here("reports", "diagnostics", "CBC.pdf"), width = 8, height = 8, pointsize = 8)

# plot RNA quality by asthma severity
pData(es_ubiopred_wb) %>%
    select(cohort, RIN) %>% 
    melt(id.vars = "cohort") %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = cohort), notch = TRUE) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "RNA Integrity Number",
        x = "",
        title = "U-BIOPRED WB") +
    ylim(5, 10)
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-1.png)

```r
# plot total WBC by asthma severity
pData(es_ubiopred_wb) %>%
    select(cohort, WBC) %>% 
    melt(id.vars = "cohort") %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = cohort), notch = TRUE) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Total White Blood Cells (10^3/uL)",
        x = "",
        title = "U-BIOPRED WB")
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-2.png)

```r
# plot estimated CBCs by asthma severity
pData(es_ubiopred_wb) %>%
    select(cohort, starts_with("EST")) %>% 
    melt(id.vars = "cohort") %>%
    mutate(variable = str_replace(variable, "EST", "")) %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = cohort), notch = TRUE) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Estimated Proportion",
        x = "Blood Cell Type",
        title = "U-BIOPRED WB") +
    ylim(0, 100)
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-3.png)

```r
# plot measured CBCs by asthma severity
pData(es_ubiopred_wb) %>%
    select(cohort, starts_with("PCT")) %>% 
    melt(id.vars = "cohort") %>%
    mutate(variable = str_replace(variable, "PCT", "")) %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = cohort), notch = TRUE) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Blood Cell Type",
        title = "U-BIOPRED WB") +
    ylim(0, 100)
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-4.png)

```r
# plot measured vs. estimated CBCs
pData(es_ubiopred_wb) %>% 
    select(cohort, starts_with("PCT"), starts_with("EST")) %>% 
    ggplot(aes(x = ESTLYMPH, y = PCTLYMPH)) + 
    geom_point() +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Estimated Proportion",
        title = "LYMPH") +
    ylim(0, 100) +
    xlim(0, 100) +
    geom_abline(intercept = 0, slope = 1, colour = "gray60", size = 0.25, linetype = "dashed") + 
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, size = 0.25) 
```

```
## Warning: Removed 3 rows containing missing values (geom_smooth).
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-5.png)

```r
pData(es_ubiopred_wb) %>% 
    select(cohort, starts_with("PCT"), starts_with("EST")) %>% 
    ggplot(aes(x = ESTMONO, y = PCTMONO)) + 
    geom_point() +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Estimated Proportion",
        title = "MONO") +
    ylim(0, 100) +
    xlim(0, 100) +
    geom_abline(intercept = 0, slope = 1, colour = "gray60", size = 0.25, linetype = "dashed") + 
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, size = 0.25) 
```

```
## Warning: Removed 42 rows containing missing values (geom_smooth).
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-6.png)

```r
pData(es_ubiopred_wb) %>% 
    select(cohort, starts_with("PCT"), starts_with("EST")) %>% 
    ggplot(aes(x = ESTNEUT, y = PCTNEUT)) + 
    geom_point() +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = "")) +
    labs(
        y = "Measured Proportion",
        x = "Estimated Proportion",
        title = "NEUT") +
    ylim(0, 100) +
    xlim(0, 100) +
    geom_abline(intercept = 0, slope = 1, colour = "gray60", size = 0.25, linetype = "dashed") + 
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, size = 0.25) 
```

```
## Warning: Removed 18 rows containing missing values (geom_smooth).
```

![plot of chunk explore_cbcs](explore_metadata//explore_cbcs-7.png)

```r
# dev.off()

# measure correlations of estimated and measured CBCs
pData(es_ubiopred_wb) %$% cor.test(ESTLYMPH, PCTLYMPH) %>% tidy
```

```
##    estimate statistic      p.value parameter  conf.low conf.high
## 1 0.7319707  23.92617 1.086425e-84       496 0.6883776 0.7702952
##                                 method alternative
## 1 Pearson's product-moment correlation   two.sided
```

```r
pData(es_ubiopred_wb) %$% cor.test(ESTMONO, PCTMONO) %>% tidy
```

```
##   estimate statistic    p.value parameter   conf.low conf.high
## 1 0.134012  3.011756 0.00272997       496 0.04669521 0.2192964
##                                 method alternative
## 1 Pearson's product-moment correlation   two.sided
```

```r
pData(es_ubiopred_wb) %$% cor.test(ESTNEUT, PCTNEUT) %>% tidy
```

```
##    estimate statistic      p.value parameter  conf.low conf.high
## 1 0.7117941   22.5693 4.009251e-78       496 0.6655531 0.7525915
##                                 method alternative
## 1 Pearson's product-moment correlation   two.sided
```

```r
pData(es_ubiopred_wb) %>% use_series(ESTMONO) %>% summary
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00000 0.00000 0.00000 0.01844 0.00000 2.02648
```

```r
pData(es_ubiopred_wb) %>% use_series(PCTMONO) %>% summary
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.740   5.513   7.000   7.065   8.640  15.970
```

```r
# measure all pair-wise CBC correlations
pData(es_ubiopred_wb) %>% select(starts_with("PCT"), starts_with("EST")) %>% correlate %>% fashion
```

```
##    rowname PCTEOS PCTLYMPH PCTMONO PCTNEUT ESTLYMPH ESTMONO ESTNEUT
## 1   PCTEOS             .07     .15    -.43      .13     .01    -.13
## 2 PCTLYMPH    .07              .19    -.90      .73    -.10    -.73
## 3  PCTMONO    .15      .19            -.41      .24     .13    -.25
## 4  PCTNEUT   -.43     -.90    -.41             -.71     .05     .71
## 5 ESTLYMPH    .13      .73     .24    -.71             -.10   -1.00
## 6  ESTMONO    .01     -.10     .13     .05     -.10             .08
## 7  ESTNEUT   -.13     -.73    -.25     .71    -1.00     .08
```

```r
pData(es_ubiopred_wb) %>% select(starts_with("PCT"), starts_with("EST")) %>% correlate %>% stretch %>% fashion %>% filter(x != y, grepl("^PCT", x))
```

```
##           x        y     r
## 1    PCTEOS PCTLYMPH   .07
## 2    PCTEOS  PCTMONO   .15
## 3    PCTEOS  PCTNEUT  -.43
## 4    PCTEOS ESTLYMPH   .13
## 5    PCTEOS  ESTMONO   .01
## 6    PCTEOS  ESTNEUT  -.13
## 7  PCTLYMPH   PCTEOS   .07
## 8  PCTLYMPH  PCTMONO   .19
## 9  PCTLYMPH  PCTNEUT  -.90
## 10 PCTLYMPH ESTLYMPH   .73
## 11 PCTLYMPH  ESTMONO  -.10
## 12 PCTLYMPH  ESTNEUT  -.73
## 13  PCTMONO   PCTEOS   .15
## 14  PCTMONO PCTLYMPH   .19
## 15  PCTMONO  PCTNEUT  -.41
## 16  PCTMONO ESTLYMPH   .24
## 17  PCTMONO  ESTMONO   .13
## 18  PCTMONO  ESTNEUT  -.25
## 19  PCTNEUT   PCTEOS  -.43
## 20  PCTNEUT PCTLYMPH  -.90
## 21  PCTNEUT  PCTMONO  -.41
## 22  PCTNEUT ESTLYMPH  -.71
## 23  PCTNEUT  ESTMONO   .05
## 24  PCTNEUT  ESTNEUT   .71
```

```r
# measure correlations of estimated CBCs
pData(es_ubiopred_wb) %>% select(starts_with("PCT")) %>% correlate %>% stretch %>% fashion %>% filter(x != y) %>% arrange(r) %>% filter(duplicated(r))
```

```
##          x        y    r
## 1 PCTLYMPH   PCTEOS  .07
## 2  PCTMONO   PCTEOS  .15
## 3  PCTMONO PCTLYMPH  .19
## 4  PCTNEUT  PCTMONO -.41
## 5  PCTNEUT   PCTEOS -.43
## 6  PCTNEUT PCTLYMPH -.90
```

```r
pData(es_ubiopred_wb) %$% cor.test(PCTNEUT, PCTEOS) %>% tidy
```

```
##     estimate statistic      p.value parameter   conf.low  conf.high
## 1 -0.4258415 -10.48184 2.346559e-23       496 -0.4951799 -0.3511126
##                                 method alternative
## 1 Pearson's product-moment correlation   two.sided
```

```r
pData(es_ubiopred_wb) %$% cor.test(PCTNEUT, PCTLYMPH) %>% tidy
```

```
##     estimate statistic       p.value parameter   conf.low  conf.high
## 1 -0.9034434 -46.93327 1.411671e-184       496 -0.9184047 -0.8859019
##                                 method alternative
## 1 Pearson's product-moment correlation   two.sided
```

```r
pData(es_ubiopred_wb) %$% cor.test(PCTNEUT, PCTMONO) %>% tidy
```

```
##     estimate statistic      p.value parameter   conf.low  conf.high
## 1 -0.4146805 -10.14913 4.063543e-22       496 -0.4848798 -0.3391721
##                                 method alternative
## 1 Pearson's product-moment correlation   two.sided
```

```r
# check relationships of CBCs with disease severity
pData(es_ubiopred_wb) %$% lm(WBC ~ cohort) %>% tidy
```

```
##                                 term estimate std.error  statistic
## 1                        (Intercept) 5.861494 0.2329215 25.1651049
## 2 cohortModerate asthma, non-smoking 0.278246 0.3399274  0.8185454
## 3   cohortSevere asthma, non-smoking 2.226351 0.2709968  8.2154154
## 4       cohortSevere asthma, smoking 2.557369 0.3284636  7.7858526
##        p.value
## 1 1.506111e-90
## 2 4.134405e-01
## 3 1.867256e-15
## 4 4.096769e-14
```

```r
pData(es_ubiopred_wb) %$% lm(PCTEOS ~ cohort) %>% tidy
```

```
##                                 term estimate std.error statistic
## 1                        (Intercept) 2.564598 0.3720788  6.892620
## 2 cohortModerate asthma, non-smoking 1.222935 0.5430145  2.252122
## 3   cohortSevere asthma, non-smoking 1.513898 0.4329019  3.497093
## 4       cohortSevere asthma, smoking 1.233811 0.5247018  2.351452
##        p.value
## 1 1.677690e-11
## 2 2.475291e-02
## 3 5.127991e-04
## 4 1.909215e-02
```

```r
pData(es_ubiopred_wb) %$% lm(PCTLYMPH ~ cohort) %>% tidy
```

```
##                                 term  estimate std.error statistic
## 1                        (Intercept) 32.037241 0.9327517 34.347022
## 2 cohortModerate asthma, non-smoking -1.612826 1.3612647 -1.184800
## 3   cohortSevere asthma, non-smoking -6.417810 1.0852270 -5.913795
## 4       cohortSevere asthma, smoking -4.921560 1.3153572 -3.741614
##         p.value
## 1 5.388653e-133
## 2  2.366664e-01
## 3  6.245951e-09
## 4  2.043079e-04
```

```r
pData(es_ubiopred_wb) %$% lm(ESTLYMPH ~ cohort) %>% tidy
```

```
##                                 term  estimate std.error statistic
## 1                        (Intercept) 29.887620 0.6702632 44.590870
## 2 cohortModerate asthma, non-smoking -1.427781 0.9781872 -1.459619
## 3   cohortSevere asthma, non-smoking -4.170467 0.7798300 -5.347918
## 4       cohortSevere asthma, smoking -3.707005 0.9451986 -3.921932
##         p.value
## 1 2.644872e-175
## 2  1.450304e-01
## 3  1.362793e-07
## 4  1.002677e-04
```

```r
pData(es_ubiopred_wb) %$% lm(PCTMONO ~ cohort) %>% tidy
```

```
##                                 term    estimate std.error  statistic
## 1                        (Intercept)  7.48264368 0.2442259 30.6382044
## 2 cohortModerate asthma, non-smoking -0.04355277 0.3564251 -0.1221933
## 3   cohortSevere asthma, non-smoking -0.60296888 0.2841491 -2.1220156
## 4       cohortSevere asthma, smoking -0.64025731 0.3444050 -1.8590245
##         p.value
## 1 2.669629e-116
## 2  9.027956e-01
## 3  3.433347e-02
## 4  6.361808e-02
```

```r
pData(es_ubiopred_wb) %$% lm(ESTMONO ~ cohort) %>% tidy
```

```
##                                 term     estimate  std.error   statistic
## 1                        (Intercept)  0.009062691 0.01668222  0.54325436
## 2 cohortModerate asthma, non-smoking -0.001736853 0.02434616 -0.07133991
## 3   cohortSevere asthma, non-smoking  0.018175802 0.01940924  0.93645102
## 4       cohortSevere asthma, smoking  0.003788878 0.02352511  0.16105678
##     p.value
## 1 0.5871997
## 2 0.9431561
## 3 0.3494986
## 4 0.8721145
```

```r
pData(es_ubiopred_wb) %$% lm(PCTNEUT ~ cohort) %>% tidy
```

```
##                                 term   estimate std.error  statistic
## 1                        (Intercept) 57.2800000  1.139988 50.2461620
## 2 cohortModerate asthma, non-smoking  0.2548052  1.663706  0.1531551
## 3   cohortSevere asthma, non-smoking  5.2228049  1.326339  3.9377589
## 4       cohortSevere asthma, smoking  4.0952273  1.607599  2.5474180
##         p.value
## 1 2.687467e-196
## 2  8.783385e-01
## 3  9.406908e-05
## 4  1.115429e-02
```

```r
pData(es_ubiopred_wb) %$% lm(ESTNEUT ~ cohort) %>% tidy
```

```
##                                 term  estimate std.error  statistic
## 1                        (Intercept) 70.103317 0.6689029 104.803435
## 2 cohortModerate asthma, non-smoking  1.429517 0.9762018   1.464367
## 3   cohortSevere asthma, non-smoking  4.152291 0.7782473   5.335439
## 4       cohortSevere asthma, smoking  3.703216 0.9432802   3.925892
##        p.value
## 1 0.000000e+00
## 2 1.437298e-01
## 3 1.454440e-07
## 4 9.868163e-05
```

