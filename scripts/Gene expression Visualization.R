#Gene expression Visualization
library(reshape2)
library(ggplot2)
# Top 10 upregulated genes
upregulated_genes <- head(significant_DEGs[order(-significant_DEGs$log2FoldChange), ], 10)
#Top 10 downregulated genes
downregulated_genes <- head(significant_DEGs[order(significant_DEGs$log2FoldChange), ], 10)

#View the top 10 upregulated and downregulated genes
upregulated_genes
downregulated_genes
dim(significant_DEGs)

# Categorizing genes based on log2FoldChange

significant_DEGs$Category <- ifelse(significant_DEGs$log2FoldChange > 0, "Significant", "Not Significant")

#VISULAIZATION OF SIGNIFICANT GENE EXPRESSION
volcano <- ggplot(significant_DEGs, aes(x = log2FoldChange, y = -log10(padj), color =Category)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted p-value") +
  theme(legend.position = "top")
# Save the volcano plot as a PNG file
ggsave("figures/volcano_plot.jpeg", plot = volcano, width = 8, height = 6, dpi = 300)




# Combine PROBEIDs of interest
key_genes <- c(upregulated_genes$PROBEID, downregulated_genes$PROBEID)

# Filter the expression matrix for the selected genes
filtered_expression_matrix <- expression_matrix[rownames(expression_matrix) %in% key_genes, ]
filtered_expression_matrix <- as.data.frame(filtered_expression_matrix)
filtered_expression_matrix$PROBEID <- rownames(filtered_expression_matrix)
# Merge expression_interest with significant_DEGs to add the SYMBOL column
filtered_expression_matrix <- merge(filtered_expression_matrix, 
                                    significant_DEGs[, c("PROBEID", "SYMBOL")], 
                                    by = "PROBEID")
# Each row in the resulting data frame corresponds to a single gene-sample combination
expression_long <- melt(filtered_expression_matrix, 
                        id.vars = c("PROBEID", "SYMBOL"), 
                        variable.name = "geo_accession", 
                        value.name = "Expression")
expression_long <- merge(expression_long, colData_subset, by = "geo_accession")


head(expression_long)
view(expression_long)

#visualization violin plots for key genes of interest across MS and healthy
VIOLIN_GG <- ggplot(expression_long, aes(x = disease_state, y = Expression, fill = disease_state)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ SYMBOL, scales = "free_y") +
  theme_minimal() +
  labs(title = "Expression of Key Genes in MS and Control Groups",
       x = "Disease State", y = "Expression Level")
ggsave("figures/Violin plots for Key Genes of Interest.jpeg",plot=VIOLIN_GG, width = 12,height = 8, dpi = 300 )


