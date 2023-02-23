# load libraries
library(ggplot2)
library(plotly)

# Enable plotting in BearPortal
options(bitmapType='cairo')


# -----------------------------------------------------
# 1. load DESeq2 results
# -----------------------------------------------------
aluminium_indium_deseq_file <- 'aluminium_indium_vs_control_DESeq2_results.csv'
aluminium_indium_dres <- read.table(aluminium_indium_deseq_file, header = T, sep = ',',  row.names = 1, check.names = F)

aluminium_deseq_file <- 'aluminium_vs_control_DESeq2_results.csv'
aluminium_dres <- read.table(aluminium_deseq_file, header = T, sep = ',',  row.names = 1, check.names = F)

indium_deseq_file <- 'indium_vs_control_DESeq2_results.csv'
indium_dres <- read.table(indium_deseq_file, header = T, sep = ',',  row.names = 1, check.names = F)

# check the loaded values
head(aluminium_indium_dres)
dim(aluminium_indium_dres)

head(aluminium_dres)
dim(aluminium_dres)

head(indium_dres)
dim(indium_dres)


# -----------------------------------------------------
# 2. calculate observed and estimated log2FC for DEGs 
# -----------------------------------------------------
# =====================================
# QUESTION 1: Please try to used original log2FoldChange values rather than the absolute values
# How does the additive plot change?
# =====================================
aluminium_indium_lfc <- aluminium_indium_dres$log2FoldChange
aluminium_lfc <- aluminium_dres$log2FoldChange
indium_lfc <- indium_dres$log2FoldChange

aluminium_indium_lfc <- abs(aluminium_indium_lfc)
aluminium_lfc <- abs(aluminium_lfc)
indium_lfc <- abs(indium_lfc)

observed_lfc <- aluminium_indium_lfc
expected_lfc <- aluminium_lfc + indium_lfc


# -----------------------------------------------------
# 3. find the overlapping DEGs
# -----------------------------------------------------
# =====================================
# QUESTION 2: Please try to used DEGs in just two comparisons:
#   a. aluminium_indium and aluminium
#   b. aluminium_indium and indium
#   c. aluminium and indium
# How does the additive plot change?
# =====================================
# =====================================
# QUESTION 3: Please try to used DEGs with different criteria:
#   a. padj < 0.01
#   b. padj < 0.1
#   c. all genes without padj cutoff
# How does the additive plot change?
# If we use pvalue rather than padj, how does the additive plot change?
# =====================================
padj.cutoff <- 0.05

# DEGs on all three comparisons
degs <- (aluminium_indium_dres$padj < padj.cutoff) & (aluminium_dres$padj < padj.cutoff) & (indium_dres$padj < padj.cutoff)
# OR, DEGs on any two comparisons
degs <- ((aluminium_indium_dres$padj < padj.cutoff) & (aluminium_dres$padj < padj.cutoff)) | 
        ((aluminium_indium_dres$padj < padj.cutoff) & (indium_dres$padj < padj.cutoff)) |
        ((aluminium_dres$padj < padj.cutoff) & (indium_dres$padj < padj.cutoff))
# OR, all genes
degs <- !(is.na(aluminium_indium_dres$padj) | is.na(aluminium_dres$padj) | is.na(indium_dres$padj))

# remove na and outliers
degs[is.na(degs)] <- F
degs[(aluminium_lfc < quantile(aluminium_lfc, 0.001)) | (aluminium_lfc > quantile(aluminium_lfc, 0.999))] <- F
degs[(indium_lfc < quantile(indium_lfc, 0.001)) | (indium_lfc > quantile(indium_lfc, 0.999))] <- F

# show the number of DEGs
sum(degs)


# -----------------------------------------------------
# 4. plot additivity and fit linear model 
# -----------------------------------------------------
plot.data <- data.frame(Observed = observed_lfc[degs], Expected = expected_lfc[degs])

p <- ggplot(plot.data, aes(Observed,Expected)) +
     geom_point() + 
     geom_smooth(method = 'lm', formula = y~x) + 
     geom_abline(intercept = -0.5, slope = 2, color = 'red', linetype = 'dashed', linewidth = 1.5)
p
#ggplotly(p)
