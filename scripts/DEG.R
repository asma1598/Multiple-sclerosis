

dds$disease_state  <- relevel(dds$disease_state, ref="multiple sclerosis")
dds$disease_state 
dds <-DESeq(dds)

res<- results(dds)
summary(res)

significant_DEGs <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
significant_DEGs
summary(significant_DEGs)

write.csv(as.data.frame(significant_DEGs), "clean_data/significant_DEGs_DESeq2.csv")

view(significant_DEGs)
dim(dds)