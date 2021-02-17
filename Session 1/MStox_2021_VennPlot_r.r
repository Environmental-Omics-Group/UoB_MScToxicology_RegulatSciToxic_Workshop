##### source for installing packages #####
if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

# install packages eg:-
# BiocManager::install("plotly") 
# install.packages("VennDiagram")

##### load installed libraries to use in project #####

library(plotly)
library(VennDiagram)


df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)
dim(df)
View(df)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])

##### Count and subset of Differential Expressed Genes (DEGs) #####
## Number of significant DEG IDs by P-value and P-adjusted (Q) cutoffs

cat(paste0('significant p-values = ', sum(AluIndVsCon$pvalue <= 0.05, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(AluIndVsCon$padj <= 0.05, na.rm = T), '\n'))

# Subset significant DEG

padjCutoff <- 0.05
AluIndVsCon.sig <- subset(AluIndVsCon, AluIndVsCon$padj <= padjCutoff)


### Subset up / down regulated genes ###

AluIndVsConUp <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange > 0 & AluIndVsCon$padj <= padjCutoff)
mean(AluIndVsConUp$log2FoldChange)

AluIndVsConDown <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange < 0 & AluIndVsCon$padj <= padjCutoff)





##### Volcano plot ######
valid_index <- which(!is.na(AluIndVsCon$padj))
df <- data.frame(pvals = AluIndVsCon$pvalue[valid_index], lfcs = AluIndVsCon$log2FoldChange[valid_index], gene = rownames(AluIndVsCon)[valid_index])
df$significant = as.factor(df$pvals <= 0.05 & abs(df$lfcs) >= 2)

p <- ggplot(df, aes(x = lfcs, y = -log10(pvals), colour = significant, gene = gene)) +
 geom_point(alpha = 0.8, size = 1) +
 labs(title = 'Volcano Plot',
 x = 'log2 fold change',
 y = '-log10 p-value') +
 scale_color_manual(values = c('#619CFF', '#F564E3'))

ggplotly(p)



##### Venn Diagram #####
# Subset significant genes and draw Venn diagram #
AluIndVsCon.sig <- subset(AluIndVsCon, AluIndVsCon$padj <= padjCutoff)
AluIndVsCon.genes <- row.names(AluIndVsCon.sig)

AluVsCon.sig <- subset(AluVsCon, AluVsCon$padj <= padjCutoff)
AluVsCon.genes <- row.names(AluVsCon.sig)

IndVsCon.sig <- subset(IndVsCon, IndVsCon$padj <= padjCutoff)
IndVsCon.genes <- row.names(IndVsCon.sig)

pdf("vennPlot_significantGenes.pdf")

vennPlot <- venn.diagram(list(AluIndVsCon.genes, AluVsCon.genes, IndVsCon.genes), NULL, fill=c("red", "green", "blue"), alpha=c(0.5,0.5,0.5), cex=3, cat.fonface=4, category.names = c("AluIndVsCon","AluVsCon", "IndVsCon"))
grid.draw(vennPlot)

dev.off()

