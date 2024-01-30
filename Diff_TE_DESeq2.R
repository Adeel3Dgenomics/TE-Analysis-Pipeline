data <- read.table("test_TE_only.cnTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",3),rep("CGroup",3)))
min_read <- 0
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds)
write.table(res, file="1_new_25hpf_control_gene_TE.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 2.000000)), ]
write.table(resSig, file="new_25hpf_control_sigdiff_gene_TE_Only.txt",sep="\t", quote=F)
