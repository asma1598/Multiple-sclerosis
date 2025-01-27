BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot"))
BiocManager::install("hgu133plus2.db")
library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)    
library(AnnotationDbi)
library(hgu133plus2.db) 

# Annotation Mapping: Entrez IDs and Gene Symbols to Probe IDs 
gene_list <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(significant_DEGs),  
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
dim(gene_list)
#Checking Null values
sum(is.na(gene_list$ENTREZID))  
sum(is.na(gene_list$SYMBOL)) 
#filter out null values
gene_list_filtered <- gene_list[!is.na(gene_list$ENTREZID) & !is.na(gene_list$SYMBOL), ]

# preparing DEGS for Merging with annotated gene_list data
significant_DEGs$PROBEID <- rownames(significant_DEGs)

significant_DEGs <- as.data.frame(significant_DEGs)

significant_DEGs<- merge(significant_DEGs, gene_list_filtered, by = "PROBEID")

dim(significant_DEGs)

#GO enrichment  analysis for BP, CC, and MF
go_enrichment <- enrichGO(
  gene =significant_DEGs$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",  
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
  ont = "BP",  
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
# Convert GO enrichment results to data frame
go_enrichment_bp_df <- as.data.frame(go_enrichment_bp)
colnames(go_enrichment_bp_df)

#visualizing GO with bar plot
fit <-plot(barplot(go_enrichment_bp, showCategory=20))

# Save the plot to a file 
ggsave("figures/Bar plot for biological process using GO.jpeg", plot = fit, width = 10, height = 18, dpi = 450)



#KEGG Pathway analysis
mapped_genes <- bitr(significant_DEGs$ENTREZID, 
                     fromType = "ENTREZID", 
                     toType = "PATH",  
                     OrgDb = org.Hs.eg.db)
head(mapped_genes)
dim(mapped_genes)
# Check the number of mapped genes
length(mapped_genes$ENTREZID)





kegg_enrichment <- enrichKEGG(
  gene = mapped_genes$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  minGSSize = 5,  
  maxGSSize = 500  
)

# Check the enrichment results in detail
head(kegg_enrichment)
summary(kegg_enrichment)  
sum(is.na(mapped_genes$ENTREZID)) 
view(kegg_enrichment)

# Visualize the KEGG results
KEGG <- dotplot(kegg_enrichment, showCategory = 10)
ggsave("figures/Dot plot for enriched pathways using KEGG.jpeg", plot = KEGG, width = 10, height = 15, dpi = 450)
KEGG


#SAVING FILES
write.csv(as.data.frame(go_enrichment), "clean_data/GO_Enrichment_ALL.csv")
write.csv(as.data.frame(go_enrichment_bp), "clean_data/GO_Enrichment_BP.csv")
write.csv(as.data.frame(kegg_enrichment), "clean_data/KEGG_Enrichment.csv")



