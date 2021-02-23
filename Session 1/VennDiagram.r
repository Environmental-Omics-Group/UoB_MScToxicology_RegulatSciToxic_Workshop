##### Venn Diagram #####
# Subset significant genes and draw Venn diagram #
library(plotly)
library(VennDiagram)

df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])

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

# AluIndVsConUp <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange > 0 & AluIndVsCon$padj <= padjCutoff)
# AluIndVsConUp.genes <- row.names(AluIndVsConUp)

# AluIndVsConDown <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange < 0 & AluIndVsCon$padj <= padjCutoff)
# AluIndVsConDown.genes <- row.names(AluIndVsConDown)