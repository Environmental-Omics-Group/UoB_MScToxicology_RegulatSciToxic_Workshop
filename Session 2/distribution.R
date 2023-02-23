# load libraries
library(reshape2)
library(ggplot2)
library(plotly)

# Enable plotting in BearPortal
options(bitmapType='cairo')


# -----------------------------------------------------
# 1. load sample sheet metadata & read counts
# -----------------------------------------------------
sample_sheet_file <- 'sample_sheet.csv'
sample_sheet <- read.table(sample_sheet_file, header = T, sep = ',', row.names = 1)

read_counts_file <- 'gene_norm_counts.csv'
read_counts <- read.table(read_counts_file, header = T, sep = ',',  row.names = 1, check.names = F)
read_counts <- as.matrix(read_counts)


# -----------------------------------------------------
# 2. match the order of sample sheet and read counts
# -----------------------------------------------------
match_index <- match(sample_sheet$Sample_Name, colnames(read_counts))
read_counts <- read_counts[,match_index]

stopifnot(sample_sheet$Sample_Name == colnames(read_counts))
stopifnot(!any(is.na(sample_sheet$Sample_Name)))
stopifnot(!any(is.na(colnames(read_counts))))

head(read_counts)
dim(read_counts)


# -----------------------------------------------------
# 3. sample distribution
# -----------------------------------------------------
df <- melt(log10(read_counts+1), varnames = c('gene', 'sample'))
df <- data.frame(df, condition = sample_sheet[df$sample,]$Condition)

p <- ggplot(df, aes(x = sample, y = value, fill = condition)) + 
  geom_boxplot() +
  ylab('log10(count + 1)')
p
# ggplotly(p)


# -----------------------------------------------------
# 4. counts distribution
# -----------------------------------------------------
df <- melt(read_counts, varnames = c('gene', 'sample'))
df <- data.frame(df, condition = sample_sheet[df$sample,]$Condition)

# df <- df[df$condition == 'Control',]
# df <- df[df$value < 1e+5,]

p <- ggplot(df, aes(x = value, colour = sample, fill = sample)) +
  geom_density(alpha = 0.2, linewidth = 1.25) + 
  facet_wrap(~ condition) +
  theme(legend.position = 'top') + xlab('read count')
p
# ggplotly(p)
