# Install and load GEOquery package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)
library(affy)
library(tidyverse)
install.packages("gtsummary")
install.packages("gt")
library(gtsummary)
library(gt)
# Fetch the dataset information
geo_data <- getGEO("GSE21942", GSEMatrix = TRUE)
expression_set <- geo_data[[1]]

# Sample annotations
metadata <- pData(expression_set)
head(metadata)            
table(metadata$group) 


# Expression matrix
expression_matrix <- exprs(expression_set)
dim(expression_matrix)            
head(expression_matrix[, 1:5])    

# Gene metadata
gene_metadata <- fData(expression_set)
head(gene_metadata)

#Saving expression matrix and associated metadata
write.csv(expression_matrix, "raw_data/expression_matrix.csv", row.names=FALSE)
write.csv(metadata, "raw_data/metadata.csv",row.names=FALSE)
write.csv(gene_metadata, "raw_data/gene_metadata.csv",row.names=FALSE)

#Initial Exploration
View(expression_matrix)
View(metadata)
View(gene_metadata)

# Check for disease state info
head(metadata$characteristics_ch1.1)

head(metadata$'disease state:ch1') 




# Summary report of normal people and patients
summary_report <- metadata |> 
  select(`disease state:ch1`, characteristics_ch1.1)|>
  tbl_summary(by = `disease state:ch1`)|>
  add_overall() |>
  modify_header(label = "**Disease State**") |>
  modify_footnote(all_stat_cols() ~ "Percentages are based on the total population")|>
  modify_caption("**Table 1: Summary of Disease States and Characteristics**") |>
  as_gt()|>
  gtsave("tables/Table2.docx")



summary_report






 
 