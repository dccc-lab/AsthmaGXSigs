library("AnnotationForge")
library("org.Hs.eg.db")
library("human.db0")

# "GPL13158" == "HT_HG_U133_Plus_PM" (Affymetrix HT HG-U133+ PM Array)
# Example: https://support.bioconductor.org/p/68567/#76261
# Affy Support Page: http://www.affymetrix.com/support/technical/byproduct.affx?product=ht_hg-u133_pm_ap

makeDBPackage(
    schema = "HUMANCHIP_DB",
    affy = TRUE,
    prefix = "hthgu133pluspm",
    chipName = "hthgu133pluspm",
    # http://www.affymetrix.com/analysis/downloads/na36/ivt/HT_HG-U133_Plus_PM.na36.annot.csv.zip
    fileName = "data/external/HT_HG-U133_Plus_PM.na36.annot.csv",
    outputDir = "data/interim/",
    baseMapType = "gbNRef",
    version = "0.1",
    manufacturer = "Affymetrix",
    manufacturerUrl = "http://www.affymetrix.com/support/technical/annotationfilesmain.affx"
)

# Not run:
# install.packages(pkgs = "data/processed/hthgu133pluspm.db", repos = NULL)

# To see all probe mapping DBs associated with this customCDF
# ls("package:hthgu133pluspm.db")