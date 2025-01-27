
#DESeq2 for differential gene expression analysis (DEG)
filtered_dds$disease_state  <- relevel(filtered_dds$disease_state, ref="multiple sclerosis")
filtered_dds$disease_state 
filtered_dds<-DESeq(filtered_dds)

#Results of DEGs
res<- results(filtered_dds)
summary(res)
View(res)
#Identifying significantly differentially expressed genes
significant_DEGs <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
significant_DEGs
summary(significant_DEGs)

#Exporting significant DEGs to csv
write.csv(as.data.frame(significant_DEGs), "clean_data/significant_DEGs_DESeq2.csv")

dim(significant_DEGs)
