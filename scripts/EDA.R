install.packages("pheatmap")

install.packages("Rtsne")
library(pheatmap)
library(Rtsne)
# Transpose the VST expression matrix 
pca_data <- t(vst_expression_matrix)

# 2. Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# 3. Extract PCA results for visualization
pca_df <- as.data.frame(pca_result$x)  
#Grouping to explore clustering between MS and control groups.
pca_df$disease_state <- colData_subset$disease_state

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = disease_state)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCA: Healthy vs Multiple Sclerosis",
    x = "Principal Component 1 (PC1)",
    y = "Principal Component 2 (PC2)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("healthy" = "green", "multiple sclerosis" = "blue"))  # Customize colors

# Display the plot
print(pca_plot)

ggsave(
  filename = "figures/PCA_plot_healthy_vs_MS.jpg",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300,
  device = "jpeg"
)


#Generating the heatmap
# 1. Calculate row variances for each gene
gene_variances <- apply(vst_expression_matrix, 1, var)

# 2. Sort genes by variance and select the top 50
top_variable_genes <- names(sort(gene_variances, decreasing = TRUE))[1:50]

# 3. Subset the vst expression matrix for the top 50 most variable genes
top_genes_matrix <- vst_expression_matrix[top_variable_genes, ]

# Generate the heatmap
heatmap_gg <- pheatmap(
  mat = top_genes_matrix,
  scale = "row",  
  clustering_distance_rows = "euclidean",  
  clustering_distance_cols = "euclidean",  
  clustering_method = "complete",  
  color = colorRampPalette(c("white", "blue"))(50),  
  main = "Heatmap of Top 50 Most Variable Genes",
  fontsize = 10,
  show_rownames = TRUE,
  show_colnames = FALSE
)
ggsave("figures/heatmap.jpeg", plot = heatmap_gg, width = 10, height = 80,limitsize = FALSE)

#T_Sne
top_pcs <- as.data.frame(pca_result$x[, 1:15]) 

tsne_results <- Rtsne(top_pcs, dims = 2, perplexity = 5, verbose = TRUE)

# 5. Extract t-SNE coordinates
# Add disease_state to tsne_coords
tsne_coords$disease_state <- vst_data$disease_state
tsne_coords <- as.data.frame(tsne_results$Y)
colnames(tsne_coords) <- c("tSNE1", "tSNE2")

# 6. Add metadata for plotting
tsne_coords$disease_state <- colData_subset$disease_state



# Plot t-SNE results with custom colors
tSNE_gg <- ggplot(tsne_coords, aes(x = tSNE1, y = tSNE2, color = disease_state)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "t-SNE Visualization of Groups",
    x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("healthy" = "blue", "multiple sclerosis" = "red"),
    labels = c("Healthy", "Multiple Sclerosis") # Custom labels for the legend
  ) +
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "right"       # Adjust legend position (optional)
  )
tSNE_gg
ggsave("figures/tSNE.jpeg", plot =tSNE_gg, width = 10, height = 12)

