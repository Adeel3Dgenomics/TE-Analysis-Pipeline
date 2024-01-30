data <- read.table("TEs_only_count.cnTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",3),rep("CGroup",3)))
min_read <- 10 ##consider only the TEs higher then 10 read counts
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
write.table(data,"clean_TE_count.txt", sep = "\t", quote = F)
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds)
write.table(res, file="expressed_TE.txt", sep="\t",quote=F)
##resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) & (abs(res$log2FoldChange)> 2.000000)), ] 
##write.table(resSig, file="new_25hpf_control_sigdiff_gene_TE_Only.txt",sep="\t", quote=F)
