# Step-by-step instructions (after installing Bioconductor and GEOquery):
library("GEOquery")

# (this line downloads GPL96.soft.gz, only needs to run once; move this file to data/)
# gpl96 <- getGEO('GPL96', destdir=".")

# Reading whole TRN into table
TRRUST <- read.csv("data/trrust_rawdata.human.tsv", sep = "\t", header=FALSE)

# Reading series matrix .txt file into ExpressionSet object:
gse_3268 <- getGEO(filename="data/GSE3268_series_matrix.txt.gz")

# Reading GPL file, used to map GEO gene identifiers to actual gene names:
gpl96 <- getGEO(filename='data/GPL96.soft.gz')

# At this point, if you have a look at the large ExpressionSet then you can see that
# you will get the information about normal vs tumour tissue in the 'phenoData/data' section
# (the authors for different datasets have it under different names, so you need to do this manually),
# whilst you can get the actual gene expression data in 'assayData/exprs'.

# The procedure to get a matrix of data ordered according to tissue category is shown below.
# Extract expression data and phenotype data
exprs_3268 <- exprs(gse_3268)
pheno_3268 <- pData(gse_3268)

# Merge the expression data with phenotype data, so that we know which samples are from tumour tissue and which aren't.
merged_3268 <- data.frame(exprs_3268)
merged_3268["TissueType", ] <- pheno_3268[["description"]]
merged_3268 <- merged_3268[c(nrow(merged_3268), 1:(nrow(merged_3268)-1)), ]
type_3268 <- merged_3268["TissueType", ]
# N.B. 'description' name is used because that is what the authors named it in this dataset.

# Calculating the mean expression values for normal vs tumour tissue
row_means_by_tissue <- function(data, tissue_types, tissue_type) {
  selected_columns <- which(tissue_types == tissue_type)
  selected_data <- data[, selected_columns]
  rowMeans(selected_data)
}

rmeans_nexpr <- row_means_by_tissue(exprs_3268, type_3268, "Normal cells")
rmeans_texpr <- row_means_by_tissue(exprs_3268, type_3268, "Tumor cells")

# Put rmeans double values into dataframe
# rmeans_expr_difference <- abs(rmeans_nexpr - rmeans_texpr)
expression_df <- data.frame(rmeans_nexpr, rmeans_texpr)

# Check uniqueness of IDs
n_repeats <- length(unique(names(rmeans_nexpr))) - length(names(rmeans_nexpr))
t_repeats <- length(unique(names(rmeans_texpr))) - length(names(rmeans_texpr))
cat("Number of repeats (N, T):", n_repeats, t_repeats)


# Translate name to name that would match TRRUST
# Ref: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/
gene_symbol_table <- Table(gpl96)[c("ID","Gene Symbol")]
ids <- gene_symbol_table[["ID"]]
symbols <- gene_symbol_table[["Gene Symbol"]]
id_symbol_map <- setNames(symbols, ids)

expression_df$gene = id_symbol_map[as.character(rownames(expression_df))]


# For convenience, to process the expression values for duplicate genes,
# and to take care of paired up genes by splitting on ///,
# we now dump everything to .csv files for processing in Python ~ Kevin
write.csv(TRRUST, "network-data/TRRUST.csv")
write.csv(expression_df, "network-data/expressions.csv")
