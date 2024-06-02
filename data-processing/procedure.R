# Step-by-step instructions (after installing Bioconductor and GEOquery):
library("GEOquery")

# Reading whole TRN into table
TRRUST <- read.csv("data/trrust_rawdata.human.tsv", sep = "\t", header=FALSE)

# Reading series matrices .txt files into ExpressionSet objects:
gse_15852 <- getGEO(filename="data/GSE15852_series_matrix.txt.gz")
gse_3268 <- getGEO(filename="data/GSE3268_series_matrix.txt.gz")

# At this point, if you have a look at the large ExpressionSets then you can see that
# you will get the information about normal vs tumour tissue in the 'phenoData/data' section (the authors for different datasets have it under different names, so you need to do this manually),
# whilst you can get the actual gene expression data in 'assayData/exprs'.


# The procedure to get a matrix of data ordered according to tissue category is shown below, here for gse_3268.
# Extract expression data and phenotype data
exprs_3268 <- exprs(gse_3268)
pheno_3268 <- pData(gse_3268)

# Merge the expression data with phenotype data, so that we know which samples are from tumour tissue and which aren't.
merged_3268 <- data.frame(exprs_3268)
merged_3268["TissueType", ] <- pheno_3268[["description"]]
merged_3268 <- merged_3268[c(nrow(merged_3268), 1:(nrow(merged_3268)-1)), ]
# N.B. 'description' name is used because that is what the authors named it

# At this point in time, this gives us the gene expression values for all of our samples, and tells us which ones come from normal tissue and which ones from cancerous.

# TO DO LIST:
# Average out the data for each type of tissue, and check uniqueness etc. of IDs
# Translate name to name that would match TRRUST
# Remove the irrelevant data points from TRRUST.