## code to prepare `DATASET` dataset goes here
library(GEOquery)

GSE_36150 <- getGEO("GSE36150", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = tempdir())
example_data <- GSE_36150$GSE36150_series_matrix.txt.gz

# library(GEOmirror)
# eSet=geoChina('GSE36150')
# example_data = eSet$GSE36150_series_matrix.txt.gz

usethis::use_data(example_data, overwrite = TRUE)

load("data-raw/merged_geneList.RData")
example_gsets <- merged_geneList
usethis::use_data(example_gsets)
