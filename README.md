# TREM-1 response signatures common to expression profiles of both asthma affection and asthma control

Croteau-Chonka, D.C. and Raby, B.A. (2018) [DOI TBD]

****

### Order of operations for analysis pipeline:

* Run `code/data/make_hthgu133pluspm_db.R`, which will create the custom annotation package you will then need to install to be able to process the U-BIOPRED GEO data (GSE69683).

* Run `code/run_analysis.R`, which will run the following scripts:

    + `code/data/create_ubiopred_eset.R`, which will generate a processed ExpressionSet with gene expression principal components and estimated blood cell counts.
    
    + `code/models/run_gsea.R`, which will run the differential gene expression (DGE) analysis and then the gene set enrichment analysis (GSEA).
    
    + `code/visualization/figure01.R`, which will make a summary figure for the primary GSEA results and provide other additional analysis details.
    
    + `code/visualization/table01.R`, which will make a summary table describing the characteristics of the U-BIOPRED study cohort.
    
The data objects the from analysis environment will end up in the `data/final/` directory and the summary figure and other analysis output will end up in the  `reports/publication/` directory.
