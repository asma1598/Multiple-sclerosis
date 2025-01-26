BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot"))
BiocManager::install("hgu133plus2.db")
library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)    
library(AnnotationDbi)
library(hgu133plus2.db) 

# Convert probe IDs to Entrez IDs
gene_list <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(significant_DEGs),  # Probe IDs as rownames
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
sum(is.na(gene_list$ENTREZID))  # Count of probe IDs with no ENTREZID mapping
sum(is.na(gene_list$SYMBOL))   # Count of probe IDs with no SYMBOL mapping
gene_list_filtered <- gene_list[!is.na(gene_list$ENTREZID) & !is.na(gene_list$SYMBOL), ]

# Merge with DEGs data
significant_DEGs$PROBEID <- rownames(significant_DEGs)

significant_DEGs <- as.data.frame(significant_DEGs)

significant_DEGs <- merge(significant_DEGs, gene_list, by = "PROBEID")
head(significant_DEGs)

significant_DEGs <- na.omit(significant_DEGs)
dim(significant_DEGs)
sum(is.na(significant_DEGs))

go_enrichment <- enrichGO(
  gene =significant_DEGs$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",  # Perform analysis for BP, CC, and MF
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_enrichment

head(as.data.frame(go_enrichment))

dim(go_enrichment)

#GO enrichment for biological process
go_enrichment_bp <- enrichGO(
  gene =significant_DEGs$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Perform analysis for BP
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
# Convert GO enrichment results to data frame
go_enrichment_bp_df <- as.data.frame(go_enrichment_bp)
colnames(go_enrichment_bp_df)

#visualizing GO with bar plot
fit <-plot(barplot(go_enrichment_bp, showCategory=20))
png("figures/Bar plot for biological process using GO", res=250, width=1200, height=1000)
# Save the plot to a file (e.g., PNG, PDF, etc.)
ggsave("figures/Bar plot for biological process using GO.jpeg", plot = fit, width = 10, height = 18, dpi = 450)



#KEGG Pathway analysis
mapped_genes <- bitr(significant_DEGs$ENTREZID, 
                     fromType = "ENTREZID", 
                     toType = "PATH",  # Map to KEGG pathway IDs
                     OrgDb = org.Hs.eg.db)
head(mapped_genes)
dim(mapped_genes)
# Check the number of mapped genes
length(mapped_genes$ENTREZID)





kegg_enrichment <- enrichKEGG(
  gene = mapped_genes$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  minGSSize = 5,  # Minimum number of genes required to consider a pathway
  maxGSSize = 500  # Maximum number of genes in a pathway
)

# Check the enrichment results in detail
head(kegg_enrichment)
summary(kegg_enrichment)  
sum(is.na(mapped_genes$ENTREZID)) 
view(kegg_enrichment)

# Visualize the KEGG results
KEGG <- dotplot(kegg_enrichment, showCategory = 10)
ggsave("figures/Dot plot for enriched pathways using KEGG.jpeg", plot = KEGG, width = 10, height = 15, dpi = 450)



#SAVING FILES
write.csv(as.data.frame(go_enrichment), "clean_data/GO_Enrichment_ALL.csv")
write.csv(as.data.frame(go_enrichment_bp), "clean_data/GO_Enrichment_BP.csv")
write.csv(as.data.frame(kegg_enrichment), "clean_data/KEGG_Enrichment.csv")



