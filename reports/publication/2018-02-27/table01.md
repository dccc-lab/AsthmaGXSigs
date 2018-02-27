

```r
library("here")

if(!exists("es_ubiopred_wb")) source(here("code", "data", "create_ubiopred_eset.R"), echo = TRUE)

library("tableone")

output_dir <- "reports/publication"
```

```r
cat_pheno_list <- c("sex", "race")
all_pheno_list <- c(cat_pheno_list, "RIN", "WBC", "PCTEOS", "PCTLYMPH", "PCTMONO", "PCTNEUT")
ubiopred_phenos <- es_ubiopred_wb %>% pData %>% select("cohort", all_pheno_list) %>% mutate(sex = fct_rev(sex))

table01 <- CreateTableOne(data = ubiopred_phenos, 
                          vars = all_pheno_list, 
                          factorVars = cat_pheno_list,
                          strata = "cohort")

table01_summary <- summary(table01)
```

```
## 
##      ### Summary of continuous variables ###
## 
## cohort: Healthy, non-smoking
##           n miss p.miss mean   sd median p25 p75 min max skew kurt
## RIN      87    0      0    9  0.7      9   8   9   6  10 -1.5  1.9
## WBC      87    0      0    6  2.0      5   5   6   3  16  1.9  6.3
## PCTEOS   87    0      0    3  1.9      2   1   3   0  13  2.2  8.5
## PCTLYMPH 87    0      0   32  9.0     31  26  38   7  59  0.4  1.0
## PCTMONO  87    0      0    7  2.1      7   6   9   3  13 -0.1 -0.4
## PCTNEUT  87    0      0   57  9.5     57  51  64  32  81 -0.3  0.5
## -------------------------------------------------------- 
## cohort: Moderate asthma, non-smoking
##           n miss p.miss mean  sd median p25 p75  min max  skew kurt
## RIN      77    0      0    9 0.5      9   9   9  6.9  10 -1.51  2.6
## WBC      77    0      0    6 1.5      6   5   7  2.8  10  0.24 -0.4
## PCTEOS   77    0      0    4 2.8      3   2   5  0.5  14  1.45  1.8
## PCTLYMPH 77    0      0   30 7.2     31  26  33 15.9  53  0.26  0.4
## PCTMONO  77    0      0    7 2.0      7   6   9  3.5  14  0.36  0.3
## PCTNEUT  77    0      0   58 8.2     57  52  63 38.7  76  0.03 -0.2
## -------------------------------------------------------- 
## cohort: Severe asthma, non-smoking
##            n miss p.miss mean   sd median p25 p75  min max skew  kurt
## RIN      246    0      0    9  0.7      9   8   9  6.5  10 -0.9  0.21
## WBC      246    0      0    8  2.4      8   6   9  3.7  17  0.9  1.03
## PCTEOS   246    0      0    4  4.2      3   1   6  0.0  43  4.1 32.54
## PCTLYMPH 246    0      0   26  8.9     27  19  31  4.6  53  0.1 -0.04
## PCTMONO  246    0      0    7  2.4      7   5   8  0.7  16  0.5  0.79
## PCTNEUT  246    0      0   63 11.3     61  55  69 31.9  93  0.3  0.16
## -------------------------------------------------------- 
## cohort: Severe asthma, smoking
##           n miss p.miss mean   sd median p25 p75 min max skew  kurt
## RIN      88    0      0    9  0.8      9   8   9   6  10 -1.3  1.20
## WBC      88    0      0    8  2.2      8   7  10   4  16  0.7  0.80
## PCTEOS   88    0      0    4  3.0      3   2   5   0  14  1.2  0.94
## PCTLYMPH 88    0      0   27  8.9     28  21  33   5  49 -0.1 -0.02
## PCTMONO  88    0      0    7  2.4      7   5   9   1  14  0.2  0.19
## PCTNEUT  88    0      0   61 11.6     60  55  68  37  93  0.3  0.07
## 
## p-values
##               pNormal   pNonNormal
## RIN      6.575854e-04 4.873021e-05
## WBC      5.657091e-22 1.007559e-23
## PCTEOS   6.721371e-03 1.205716e-02
## PCTLYMPH 4.327264e-09 1.243556e-08
## PCTMONO  5.962229e-02 2.644836e-02
## PCTNEUT  4.625806e-05 1.767656e-04
## 
## Standardize mean differences
##            average     1 vs 2    1 vs 3    1 vs 4     2 vs 3      2 vs 4
## RIN      0.3060675 0.26170160 0.2646772 0.2373018 0.56072547 0.502627874
## WBC      0.7856133 0.15807803 1.0139340 1.2147002 0.97819251 1.204996526
## PCTEOS   0.2712039 0.51118565 0.4681362 0.4852255 0.08213391 0.003718146
## PCTLYMPH 0.4385181 0.19764932 0.7142820 0.5495534 0.59226380 0.409480600
## PCTMONO  0.1854432 0.02092016 0.2678260 0.2840410 0.25377633 0.270346919
## PCTNEUT  0.3161670 0.02857879 0.5002149 0.3852678 0.50315209 0.381266201
##               3 vs 4
## RIN      0.009371157
## WBC      0.143778626
## PCTEOS   0.076824274
## PCTLYMPH 0.167879567
## PCTMONO  0.015748692
## PCTNEUT  0.098522349
## 
## =======================================================================================
## 
##      ### Summary of categorical variables ### 
## 
## cohort: Healthy, non-smoking
##   var  n miss p.miss                 level freq percent cum.percent
##   sex 87    0    0.0                  male   53    60.9        60.9
##                                     female   34    39.1       100.0
##                                                                    
##  race 87    0    0.0       white_caucasian   80    92.0        92.0
##                           south_east_asian    0     0.0        92.0
##                                      other    0     0.0        92.0
##                              black_african    3     3.4        95.4
##                                south_asian    1     1.1        96.6
##                      arabic_north_heritage    0     0.0        96.6
##                             multiple_races    2     2.3        98.9
##                              central_asian    1     1.1       100.0
##                                 east_asian    0     0.0       100.0
##                                                                    
## -------------------------------------------------------- 
## cohort: Moderate asthma, non-smoking
##   var  n miss p.miss                 level freq percent cum.percent
##   sex 77    0    0.0                  male   40    51.9        51.9
##                                     female   37    48.1       100.0
##                                                                    
##  race 77    0    0.0       white_caucasian   72    93.5        93.5
##                           south_east_asian    2     2.6        96.1
##                                      other    0     0.0        96.1
##                              black_african    3     3.9       100.0
##                                south_asian    0     0.0       100.0
##                      arabic_north_heritage    0     0.0       100.0
##                             multiple_races    0     0.0       100.0
##                              central_asian    0     0.0       100.0
##                                 east_asian    0     0.0       100.0
##                                                                    
## -------------------------------------------------------- 
## cohort: Severe asthma, non-smoking
##   var   n miss p.miss                 level freq percent cum.percent
##   sex 246    0    0.0                  male   85    34.6        34.6
##                                      female  161    65.4       100.0
##                                                                     
##  race 246    0    0.0       white_caucasian  215    87.4        87.4
##                            south_east_asian    1     0.4        87.8
##                                       other    5     2.0        89.8
##                               black_african    8     3.3        93.1
##                                 south_asian    7     2.8        95.9
##                       arabic_north_heritage    5     2.0        98.0
##                              multiple_races    4     1.6        99.6
##                               central_asian    0     0.0        99.6
##                                  east_asian    1     0.4       100.0
##                                                                     
## -------------------------------------------------------- 
## cohort: Severe asthma, smoking
##   var  n miss p.miss                 level freq percent cum.percent
##   sex 88    0    0.0                  male   45    51.1        51.1
##                                     female   43    48.9       100.0
##                                                                    
##  race 88    0    0.0       white_caucasian   84    95.5        95.5
##                           south_east_asian    0     0.0        95.5
##                                      other    0     0.0        95.5
##                              black_african    1     1.1        96.6
##                                south_asian    2     2.3        98.9
##                      arabic_north_heritage    0     0.0        98.9
##                             multiple_races    1     1.1       100.0
##                              central_asian    0     0.0       100.0
##                                 east_asian    0     0.0       100.0
##                                                                    
## 
## p-values
##           pApprox       pExact
## sex  4.873413e-05 4.480032e-05
## race 2.324385e-01           NA
## 
## Standardize mean differences
##        average    1 vs 2    1 vs 3    1 vs 4    2 vs 3     2 vs 4
## sex  0.2733085 0.1816789 0.5472820 0.1980653 0.3566572 0.01624203
## race 0.3763677 0.3874787 0.3804501 0.2527024 0.4779624 0.39586299
##         3 vs 4
## sex  0.3399256
## race 0.3637497
```

```r
table01_print <- print(table01, nonnormal = TRUE)
```

```
##                           Stratified by cohort
##                            Healthy, non-smoking
##   n                           87               
##   sex = female (%)            34 (39.1)        
##   race (%)                                     
##      white_caucasian          80 (92.0)        
##      south_east_asian          0 ( 0.0)        
##      other                     0 ( 0.0)        
##      black_african             3 ( 3.4)        
##      south_asian               1 ( 1.1)        
##      arabic_north_heritage     0 ( 0.0)        
##      multiple_races            2 ( 2.3)        
##      central_asian             1 ( 1.1)        
##      east_asian                0 ( 0.0)        
##   RIN (median [IQR])        9.00 [8.50, 9.20]  
##   WBC (median [IQR])        5.33 [4.70, 6.40]  
##   PCTEOS (median [IQR])     2.17 [1.50, 3.12]  
##   PCTLYMPH (median [IQR])  31.20 [26.12, 37.91]
##   PCTMONO (median [IQR])    7.40 [6.12, 9.15]  
##   PCTNEUT (median [IQR])   57.14 [50.96, 63.52]
##                           Stratified by cohort
##                            Moderate asthma, non-smoking
##   n                           77                       
##   sex = female (%)            37 (48.1)                
##   race (%)                                             
##      white_caucasian          72 (93.5)                
##      south_east_asian          2 ( 2.6)                
##      other                     0 ( 0.0)                
##      black_african             3 ( 3.9)                
##      south_asian               0 ( 0.0)                
##      arabic_north_heritage     0 ( 0.0)                
##      multiple_races            0 ( 0.0)                
##      central_asian             0 ( 0.0)                
##      east_asian                0 ( 0.0)                
##   RIN (median [IQR])        9.00 [8.80, 9.30]          
##   WBC (median [IQR])        5.80 [5.10, 7.20]          
##   PCTEOS (median [IQR])     2.90 [1.88, 4.55]          
##   PCTLYMPH (median [IQR])  30.53 [26.40, 33.33]        
##   PCTMONO (median [IQR])    7.27 [5.80, 8.95]          
##   PCTNEUT (median [IQR])   56.86 [52.20, 63.40]        
##                           Stratified by cohort
##                            Severe asthma, non-smoking
##   n                          246                     
##   sex = female (%)           161 (65.4)              
##   race (%)                                           
##      white_caucasian         215 (87.4)              
##      south_east_asian          1 ( 0.4)              
##      other                     5 ( 2.0)              
##      black_african             8 ( 3.3)              
##      south_asian               7 ( 2.8)              
##      arabic_north_heritage     5 ( 2.0)              
##      multiple_races            4 ( 1.6)              
##      central_asian             0 ( 0.0)              
##      east_asian                1 ( 0.4)              
##   RIN (median [IQR])        8.75 [8.20, 9.10]        
##   WBC (median [IQR])        7.79 [6.30, 9.40]        
##   PCTEOS (median [IQR])     2.90 [1.40, 6.17]        
##   PCTLYMPH (median [IQR])  26.51 [18.70, 31.09]      
##   PCTMONO (median [IQR])    6.78 [5.39, 8.36]        
##   PCTNEUT (median [IQR])   61.11 [54.90, 69.03]      
##                           Stratified by cohort
##                            Severe asthma, smoking p      test   
##   n                           88                                
##   sex = female (%)            43 (48.9)           <0.001        
##   race (%)                                         0.232        
##      white_caucasian          84 (95.5)                         
##      south_east_asian          0 ( 0.0)                         
##      other                     0 ( 0.0)                         
##      black_african             1 ( 1.1)                         
##      south_asian               2 ( 2.3)                         
##      arabic_north_heritage     0 ( 0.0)                         
##      multiple_races            1 ( 1.1)                         
##      central_asian             0 ( 0.0)                         
##      east_asian                0 ( 0.0)                         
##   RIN (median [IQR])        8.75 [8.28, 9.10]     <0.001 nonnorm
##   WBC (median [IQR])        7.94 [6.74, 9.83]     <0.001 nonnorm
##   PCTEOS (median [IQR])     2.88 [1.58, 5.14]      0.012 nonnorm
##   PCTLYMPH (median [IQR])  27.55 [21.27, 32.98]   <0.001 nonnorm
##   PCTMONO (median [IQR])    6.79 [5.08, 8.53]      0.026 nonnorm
##   PCTNEUT (median [IQR])   60.39 [54.64, 68.01]   <0.001 nonnorm
```

```r
write.csv(table01_print, file.path(here(), output_dir, "table01.txt"))
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
##  class                  7.3-14     2015-08-30
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
##  e1071                  1.6-8      2017-02-02
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
##  markdown               0.8        2017-04-20
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
##  rstudioapi             0.7        2017-09-07
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
##  survey                 3.32-1     2017-06-22
##  survival               2.41-3     2017-04-04
##  tableone             * 0.8.1      2017-06-17
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
##  CRAN (R 3.4.0)                    
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

