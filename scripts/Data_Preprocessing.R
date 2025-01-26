BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
install.packages("gridExtra")
#checking expression matrix data type and converting it into integer for dseq obj
class(expression_matrix)
summary(expression_matrix)
expression_matrix_int <- ceiling(expression_matrix)
summary(expression_matrix_int)
head(expression_matrix_int)
dim(expression_matrix)
dim(gene_metadata)
expression_matrix_int[1:10, 1:10] 

#prepare metadata to align with expression_matrix
colData_subset <- metadata[, c("disease state:ch1", "geo_accession")]

# Rename the column "disease state:ch1" to "disease state"
colnames(colData_subset)[colnames(colData_subset) == "disease state:ch1"] <- "disease_state"

# Check the updated column names
colnames(colData_subset)
view(colData_subset)
colnames(expression_matrix_int)
rownames(expression_matrix_int)
colnames(gene_metadata)
rownames(gene_metadata)

all(colnames(expression_matrix_int)== rownames(colData_subset))
all(colnames(expression_matrix_int)%in% rownames(colData_subset))
# 1. filtering out
keep <- rowSums(expression_matrix_int >= 10) > (0.5 * ncol(expression_matrix_int))
filtered_expression_matrix <- expression_matrix_int[keep, ]
all(colnames(filtered_expression_matrix)%in% rownames(colData_subset))
all(colnames(filtered_expression_matrix)== rownames(colData_subset))

# 2. Create DESeqDataSet object using the filtered count data
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_matrix,  
  colData = colData_subset,               
  design = ~ disease_state               
)
dds
# 3. Perform VST transformation
vst_expression_matrix <- varianceStabilizingTransformation(dds)

# 4. Access the transformed expression matrix
vst_expression_matrix <- assay(vst_expression_matrix)

# 5. Check the transformed data
head(vst_expression_matrix[, 1:5])

summary(vst_expression_matrix)

# install.packages("gridExtra")
library(gridExtra)

# Boxplot for raw  data
raw_data <- data.frame(t((expression_matrix_int + 1)))
colnames(raw_data) <- colnames(expression_matrix_int)

raw_plot <- ggplot(raw_data, aes(x = factor(1), y = raw_data[,1])) +
  geom_boxplot() +
  xlab("Samples") +
  ylab("Raw Data Expression") +
  ggtitle("Boxplot of Raw Data")

# Boxplot for normalized (VST-transformed) data
normalized_data <- data.frame(t(vst_expression_matrix))
colnames(normalized_data) <- colnames(vst_expression_matrix)

normalized_plot <- ggplot(normalized_data, aes(x = factor(1), y = normalized_data[,1])) +
  geom_boxplot() +
  xlab("Samples") +
  ylab("Normalized Expression") +
  ggtitle("Boxplot of Normalized Data (VST)")

# Arrange both plots side by side
normalization <- grid.arrange(raw_plot, normalized_plot, ncol = 2)

ggsave(
  filename = "figures/box_plot.jpg",  
  plot =  normalization,              
  width = 10,                             
  height = 6,                             
  dpi = 300,                              
  device = "jpeg"                        
)



