##### Set working direcotry #####
getwd()
setwd("/Users/vigneshdh/Desktop/MStox2021") # Change the path to your woring directory


##### source for installing packages #####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install packages eg:-
# BiocManager::install("plotly") 
# install.packages("VennDiagram")

##### load intalled libraries to use in project #####
library(plotly)
library(VennDiagram)


df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)
dim(df)
View(df)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])


##### Count and subset DEGs #####
## Number of Differentially Expresed  Genes (DEGs)
cat(paste0('significant p-values = ', sum(AluIndVsCon$pvalue <= 0.05, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(AluIndVsCon$padj <= 0.05, na.rm = T), '\n'))

# Subset significant DEGs gene IDs by P-value and P-adjacent cutoffs
padjCutoff <- 0.05
pvalCutoff <- 0.05
lfcCutoff <- 2

AluIndVsCon.sig <- row.names(AluIndVsCon)[which(AluIndVsCon$padj <= padjCutoff & AluIndVsCon$pvalue <= pvalCutoff)]
length(AluIndVsCon.sig)

### Subset up / down regulated genes ###

AluIndVsConUp.genes <- row.names(AluIndVsCon)[which(AluIndVsCon$padj <= padjCutoff & AluIndVsCon$log2FoldChange >= lfcCutoff)]

AluIndVsConDown.genes <- row.names(AluIndVsCon)[which(AluIndVsCon$padj <= padjCutoff & AluIndVsCon$log2FoldChange >= -lfcCutoff)]


##### Volcano plot ######
valid_index <- which(!is.na(AluIndVsCon$padj))
df <- data.frame(pvals = AluIndVsCon$padj[valid_index], lfcs = AluIndVsCon$log2FoldChange[valid_index], gene = rownames(AluIndVsCon)[valid_index])
df$significant = as.factor(df$pvals <= 0.05 & abs(df$lfcs) >= 2)

p <- ggplot(df, aes(x = lfcs, y = -log10(pvals), colour = significant, gene = gene)) + 
  geom_point(alpha = 0.8, size = 1) +
  labs(title = 'Volcano Plot', 
       x = 'log2 fold change', 
       y = '-log10 p-value') +
  scale_color_manual(values = c('#619CFF', '#F564E3'))

ggplotly(p)


##### Venn Diagram #####

# Subset significant genes and draw venn diagram #

AluIndVsCon.sig <- subset(AluIndVsCon, AluIndVsCon$padj <= padjCutoff)
AluIndVsCon.genes <- row.names(AluIndVsCon.sig)

AluVsCon.sig <- subset(AluVsCon, AluVsCon$padj <= padjCutoff)
AluVsCon.genes <- row.names(AluVsCon.sig)

IndVsCon.sig <- subset(IndVsCon, IndVsCon$padj <= padjCutoff)
IndVsCon.genes <- row.names(IndVsCon.sig)

#pdf("vennPlot_significantGenes.pdf")

vennPlot <- venn.diagram(list(AluIndVsCon.genes, AluVsCon.genes, IndVsCon.genes), NULL, fill=c("red", "green", "blue"), alpha=c(0.5,0.5,0.5), cex=3, cat.fonface=4, category.names = c("AluIndVsCon","AluVsCon", "IndVsCon"))
grid.draw(vennPlot)

dev.off()





