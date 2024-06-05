# Step-by-step instructions (after installing Bioconductor and GEOquery):
library("GEOquery")

# (this line downloads GPL96.soft.gz, only needs to run once)
# gpl96 <- getGEO('GPL96', destdir="data/")

# Reading entire TRN information into table
#TRRUST <- read.csv("data/trrust_rawdata.human.tsv", sep = "\t", header=FALSE)

# Download TFLink from https://tflink.net/download/ - Homo sapiens small and large-scale interaction table
TFLink <- read.csv("data/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz", sep = "\t", header=TRUE) 
TFLink <- TFLink[c("Name.TF", "Name.Target")]
colnames(TFLink) <- c("TF", "Target")

# Reading STRING database of proteins into table, for weighting edges
# Also reading STRING translation 'table' for proteins
STRING <- read.csv("data/9606.protein.links.v12.0.txt.gz", sep = " ", header=TRUE)
STRINGtrans <- read.csv("data/9606.protein.info.v12.0.txt.gz", sep = "\t", header=TRUE)

# Reading series matrix .txt file into ExpressionSet object:
gse_15852 <- getGEO(filename="data/GSE15852_series_matrix.txt.gz")

# Reading GPL file, used to map GEO gene identifiers to actual gene names:
gpl96 <- getGEO(filename='data/GPL96.soft.gz')

# At this point, if you have a look at the large ExpressionSet then you can see that
# you will get the information about normal vs tumour tissue in the 'phenoData/data' section
# (the authors for different datasets have it under different names, so you need to do this manually),
# whilst you can get the actual gene expression data in 'assayData/exprs'.

# The procedure to get a matrix of data ordered according to tissue category is shown below.
# Extract expression data and phenotype data
exprs_15852 <- exprs(gse_15852)
pheno_15852 <- pData(gse_15852)

# Merge the expression data with phenotype data, so that we know which samples are from tumour tissue and which aren't.
merged_15852 <- data.frame(exprs_15852)
merged_15852["TissueType", ] <- pheno_15852[["description"]]
merged_15852 <- merged_15852[c(nrow(merged_15852), 1:(nrow(merged_15852)-1)), ]
type_15852 <- merged_15852["TissueType", ]
# N.B. 'description' name is used because that is what the authors named it in this dataset.

# Calculating the mean expression values for normal vs tumour tissue
row_means_by_tissue <- function(data, tissue_types, tissue_type) {
  selected_columns <- which(tissue_types == tissue_type)
  selected_data <- data[, selected_columns]
  rowMeans(selected_data)
}

rmeans_nexpr <- row_means_by_tissue(exprs_15852, type_15852, "gene expression data from normal breast tissue")
rmeans_texpr <- row_means_by_tissue(exprs_15852, type_15852, "gene expression data from breast tumor")

# ~ Winner

# Put rmeans double values into dataframe
# rmeans_expr_difference <- abs(rmeans_nexpr - rmeans_texpr)
expression_df <- data.frame(rmeans_nexpr, rmeans_texpr)


# Check uniqueness of IDs
n_repeats <- length(unique(names(rmeans_nexpr))) - length(names(rmeans_nexpr))
t_repeats <- length(unique(names(rmeans_texpr))) - length(names(rmeans_texpr))
cat("Number of repeats (N, T):", n_repeats, t_repeats)


# Translate GEO ID to name that would match standard as in TRRUST
# Ref: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/
gene_symbol_table <- Table(gpl96)[c("ID","Gene Symbol")]
ids <- gene_symbol_table[["ID"]]
symbols <- gene_symbol_table[["Gene Symbol"]]
id_symbol_map <- setNames(symbols, ids)

expression_df$gene = id_symbol_map[as.character(rownames(expression_df))]

# Rename TRRUST columns
# https://pubmed.ncbi.nlm.nih.gov/29087512/#&gid=article-figures&pid=figure-3-uid-2
#colnames(TRRUST) <- c("name1","name2","mode","references")

# Translate STRING ID to name that would match standard as in TRRUST
ids <- STRINGtrans$X.string_protein_id
symbols <- STRINGtrans$preferred_name
STRING_symbol_map <- setNames(symbols, ids)

# To save memory, drop the STRING IDs from the database
STRING$name1 = STRING_symbol_map[STRING$protein1]
STRING$name2 = STRING_symbol_map[STRING$protein2]

STRING <- STRING[c('TF', 'Target', 'combined_score')]

# For convenience, to process the expression values for duplicate genes,
# and to take care of paired up genes by splitting on ///,
# we now dump everything to .csv files for processing in Python ~ Kevin
#write.csv(TRRUST, "network-data/TRRUST.csv")
write.csv(TFLink, "network-data/TFLink.csv")
write.csv(expression_df, "network-data/expressions.csv")
write.csv(STRING, "network-data/STRING_by_gene.csv")
