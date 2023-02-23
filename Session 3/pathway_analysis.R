# read DESeq2 results
df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)

dim(df)
head(df)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])


# number of significant DEG IDs by P-value and P-adjusted (Q) cutoff
pcutoff <- 0.05

AluVsConSig <- subset(AluVsCon, AluVsCon$padj <= pcutoff)
cat(paste0('significant p-adjust = ', nrow(AluVsConSig), '\n'))

IndVsConSig <- subset(IndVsCon, IndVsCon$padj <= pcutoff)
cat(paste0('significant p-adjust = ', nrow(IndVsConSig), '\n'))

AluIndVsConSig <- subset(AluIndVsCon, AluIndVsCon$padj <= pcutoff)
cat(paste0('significant p-adjust = ', nrow(AluIndVsConSig), '\n'))


# subset up (down) regulated genes
AluVsConSigUp <- subset(AluVsConSig, AluVsConSig$log2FoldChange > 0)
cat(paste0('significant p-adjust and up-regulated = ', nrow(AluVsConSigUp), '\n'))

IndVsConSigUp <- subset(IndVsConSig, IndVsConSig$log2FoldChange > 0)
cat(paste0('significant p-adjust and up-regulated = ', nrow(IndVsConSigUp), '\n'))

AluIndVsConSigUp <- subset(AluIndVsConSig, AluIndVsConSig$log2FoldChange > 0)
cat(paste0('significant p-adjust and up-regulated = ', nrow(AluIndVsConSigUp), '\n'))


# mapping genes
orthoMap <- read.csv("odb10v1_dma_to_hsa_level_33208_gn.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

RefOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(df))
RefOrtho <- unique(unlist(strsplit(RefOrtho$V2, split = ';')))
length(RefOrtho)
writeLines(RefOrtho, 'reference_orthologs_hsa.txt')

AluVsConSigOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(AluVsConSig))
AluVsConSigOrtho <- unique(unlist(strsplit(AluVsConSigOrtho$V2, split = ';')))
length(AluVsConSigOrtho)
writeLines(AluVsConSigOrtho, 'AluVsConSig_orthologs_hsa.txt')

IndVsConSigOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(IndVsConSig))
IndVsConSigOrtho <- unique(unlist(strsplit(IndVsConSigOrtho$V2, split = ';')))
length(IndVsConSigOrtho)
writeLines(IndVsConSigOrtho, 'IndVsConSigOrtho_orthologs_hsa.txt')

AluIndVsConSigOrtho <- subset(orthoMap, row.names(orthoMap) %in% row.names(AluIndVsConSig))
AluIndVsConSigOrtho <- unique(unlist(strsplit(AluIndVsConSigOrtho$V2, split = ';')))
length(AluIndVsConSigOrtho)
writeLines(AluIndVsConSigOrtho, 'AluIndVsConSigOrtho_orthologs_hsa.txt')

