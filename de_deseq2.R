###A base editing platform for the correction of cancer driver mutations unmasks conserved p53 transcription programs###

###Differential expression analysis with DESeq2###
library(magrittr)
library(ggplot2)
library(DESeq2)
library(gridExtra)
library(ComplexHeatmap)
library(vsn)
library(pheatmap)
library(EnhancedVolcano)
library(writexl)
library(tidyverse)

setwd("../results/")

###########
#####input counts####
#reading in the counts after quantification with featureCounts
temp <- list.files("./count_matrices/", pattern = "*.txt$")
col <- gsub("counts_", "", gsub(".txt", "", temp))

filenames <- list.files(path = "./count_matrices/", pattern = "*.txt$", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, 
                                                     sep = "\t", stringsAsFactors = FALSE)[,6]})

#counts
count_matrix <- sapply(datalist, cbind)
nrow(count_matrix)
#63140
ncol(count_matrix)
#54

#sample names
colnames(count_matrix) <- col
rownames(count_matrix) <- (rownames(read.table("./count_matrices/counts_A_273_36_1.txt", row.names = 1, header = TRUE)))

write.csv(count_matrix, file = "de/raw_counts.csv")

colData <- data.frame(col)
names(colData) <- "col"

colData$time <- NULL
colData[grep("48", colData$col), "time"] <- "48"
colData[grep("36", colData$col), "time"] <- "36"
colData[grep("72", colData$col), "time"] <- "72"
colData[grep("WT", colData$col), "time"] <- "WT"
colData[45:47, "time"] <- "72"
colData$time <- factor(colData$time,levels = c("36", "48", "72", "WT"))


######extract gene names from gtf########
x <- data.frame(rtracklayer::import("../genome/gencode_v46_annotation.gtf"))
u <- unique(x[,c("gene_id", "gene_name")])
write.csv(u,"gene_names.csv")

#assigning gene names to gene ids
gene_names <- read.csv("gene_names.csv")
gene_names <- gene_names[, !(names(gene_names) == "X")]
nrow(gene_names) == nrow(count_matrix)
#TRUE


#add column for gene names matching w/ rownames in count matrix
count_matrix$rows <- gene_names[gene_names$gene_id == rownames(count_matrix), "gene_name"]
#OR match
count_matrix$rows <- gene_names$gene_name[match(rownames(count_matrix), gene_names$gene_id)]
sum(!is.null(count_matrix$rows))

#cleaning up rownames
rownames(count_matrix) <- make.names(count_matrix$rows, unique = T)
count_matrix <- count_matrix[, !(names(count_matrix) == "rows")]




####exploration####
######DESeq2######
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData[1:2],
                              design = ~ time)
dds


head(colData(dds))
head(rowData(dds))
assay(dds, "counts") %>% head


#getting rid of genes with 0 counts
dds <- dds[rowSums(counts(dds)) > 0,]
colSums(counts(dds))

#mean variance plot
meanSdPlot(as.matrix(count_matrix), ranks = F)

#size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

colData(dds)
colData <- colData(dds)[1:2]

###
##transformation of raw counts
vsd <- vst(dds, blind = F)
#n of samples = 54
meanSdPlot(assay(vsd), ranks = F)



#####PCA#####
var_genes <- order(rowVars(as.matrix(count_matrix)), decreasing = T)
#OR
#var_genes <- order(rowVars(counts(dds, normalized = T)), decreasing = T)
pca <- prcomp(t(log(count_matrix[var_genes,] + 1)), scale = FALSE)
variance_prop <- round(as.data.frame(summary(pca)[[6]])[2,],3)
variance_prop <- c(variance_prop[1], variance_prop[2])

colData$colf <- substring(colData$col, 1, 5)
coldata_pca <- colData
coldata_pca <- microeco::dropallfactors(coldata_pca)

##if a PCA needs to be done without the SMAD4 mutation
coldata_pca[coldata_pca$colf == "H_VEG" | coldata_pca$colf == "A_VEG", "time" ] <- "CTRL"
coldata_pca <- coldata_pca[!coldata_pca$colf == "H_311",]
coldata_pca[coldata_pca$time == "WT", "time"] <- "CTRL"
coldata_pca[coldata_pca$colf == "A_273" | coldata_pca$colf == "A_VEG", "colf" ] <- "A431"
coldata_pca[coldata_pca$colf == "H_273" | coldata_pca$colf == "H_VEG", "colf" ] <- "HT-29"

pc <- ggplot(data.frame(PC1 = pca$x[,"PC1"], PC2 = pca$x[,"PC2"], col = factor(coldata_pca$colf), time = coldata_pca$time), 
             aes(x = PC1, y = PC2, color = col, shape = time)) +
  geom_point(size = 3.5)+
  scale_colour_manual(values = c(
    "firebrick",  "darkgreen", "grey50"))+
  labs(
    x = paste("PC1 (",variance_prop[[1]]*100, "% )", sep = ""),
    y = paste("PC2 (",variance_prop[[2]]*100, "% )", sep = ""),
    color="cell lines",
    shape="h.p.i"
  )+
  theme_minimal()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size=14))+
  guides(color = guide_legend(override.aes = list(size = 8)))+
  theme(legend.key.size = unit(0.8, 'cm'))
pc
ggsave("plots/pca/pca__hpi.svg", plot = pc,
       dpi = 300, device = grDevices::svg, width = 12, height = 8)



####p53 target list####
target_list <- readLines("gene_list_p53target_genes.list")


#####comparison for each cell line#####
####for each timepoint, edited vs CTRL/WT###

#####PANC1 at 36h####

colData_pair <- colData[colData$colf == "PANC1" & (colData$time == 36 | colData$time == "WT"),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~time)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#31311

res_PANC1_WT_vs_36 <- results(dds_pair, contrast = c("time", "36", "WT"))
res_PANC1_WT_vs_36 <- data.frame(res_PANC1_WT_vs_36)
res_PANC1_WT_vs_36 <- res_PANC1_WT_vs_36[complete.cases(res_PANC1_WT_vs_36),]


res_PANC1_WT_vs_36_sig <- res_PANC1_WT_vs_36[res_PANC1_WT_vs_36$padj < 0.05,]
positive_PANC1_WT_vs_36 <- res_PANC1_WT_vs_36_sig[res_PANC1_WT_vs_36_sig$log2FoldChange > 1,]

negative_PANC1_WT_vs_36 <- res_PANC1_WT_vs_36_sig[res_PANC1_WT_vs_36_sig$log2FoldChange < -1,]

nrow(positive_PANC1_WT_vs_36) + nrow(negative_PANC1_WT_vs_36)
#2233

#upregulated
top_genes <- c(rownames(positive_PANC1_WT_vs_36[order(positive_PANC1_WT_vs_36$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_PANC1_WT_vs_36))
#downregulated
top_genes <- c(rownames(negative_PANC1_WT_vs_36[order(negative_PANC1_WT_vs_36$log2FoldChange, decreasing = F),]),
               rownames(positive_PANC1_WT_vs_36))

select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(, 
  log10(select[,order(df$time)] + 1)[1:30,], 
  cluster_rows=F, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T, main = "PANC1 cell line, WT vs 36h, up@36") 


###p53 target genes###
#upregulated at 36
p53_PANC1_WT_36_at_36 <- intersect(rownames(positive_PANC1_WT_vs_36), target_list)
#5
#up at WT
p53_PANC1_WT_36_at_WT <- intersect(rownames(negative_PANC1_WT_vs_36), target_list)
#7

######modifed volcano####
p53_PANC1_36 <- union(p53_PANC1_WT_36_at_36, p53_PANC1_WT_36_at_WT)

keyvals <- ifelse(p53_PANC1_36 %in% rownames(res_PANC1_WT_vs_36), 'gold1',
                  ifelse(res_PANC1_WT_vs_36$log2FoldChange < -1 & res_PANC1_WT_vs_36$padj < 0.05, 'skyblue3', 
                         ifelse(res_PANC1_WT_vs_36$log2FoldChange > 1 & res_PANC1_WT_vs_36$padj < 0.05, 'firebrick2','gray65')))
keyvals[is.na(keyvals)] <- 'gray65'

res_PANC1_WT_vs_36$keyvals <- ""
res_PANC1_WT_vs_36$keyvals <- ifelse(res_PANC1_WT_vs_36$log2FoldChange < -1 & res_PANC1_WT_vs_36$padj < 0.05, 'skyblue3', res_PANC1_WT_vs_36$keyvals)
res_PANC1_WT_vs_36$keyvals <- ifelse(res_PANC1_WT_vs_36$log2FoldChange > 1 & res_PANC1_WT_vs_36$padj < 0.05, 'lightcoral',res_PANC1_WT_vs_36$keyvals)
res_PANC1_WT_vs_36$keyvals <- ifelse(rownames(res_PANC1_WT_vs_36) %in% p53_PANC1_36, 'springgreen4', res_PANC1_WT_vs_36$keyvals)
res_PANC1_WT_vs_36[res_PANC1_WT_vs_36$keyvals == "", "keyvals"] <- 'gray80'


keyvals <- res_PANC1_WT_vs_36$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_PANC1_WT_vs_36[order(res_PANC1_WT_vs_36$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_PANC1_WT_vs_36[order(res_PANC1_WT_vs_36$keyvals, decreasing = F),],
                     lab = rownames(res_PANC1_WT_vs_36[order(res_PANC1_WT_vs_36$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 12,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'PANC1 36 vs WT',
                     caption = bquote("Total DE genes: 2233 (overlapping 12)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)

p <- p + coord_cartesian(xlim = c(-8, 10)) +
  scale_x_continuous(breaks = seq(-8,10,1))
ggsave("plots/volcanoes_overlapping_116_targets/PANC1_@_36.svg", plot = p, 
       width = 10, height = 8, dpi = 300, device = grDevices::svg)


####PANC1 at 48h####
colData_pair <- colData[colData$colf == "PANC1" & (colData$time == 48 | colData$time == "WT"),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~time)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#31585


res_PANC1_WT_vs_48 <- results(dds_pair, contrast = c("time", "48", "WT"))
res_PANC1_WT_vs_48 <- data.frame(res_PANC1_WT_vs_48)
res_PANC1_WT_vs_48 <- res_PANC1_WT_vs_48[complete.cases(res_PANC1_WT_vs_48),]

res_PANC1_WT_vs_48_sig <- res_PANC1_WT_vs_48[res_PANC1_WT_vs_48$padj < 0.05,]
positive_PANC1_WT_vs_48 <- res_PANC1_WT_vs_48_sig[res_PANC1_WT_vs_48_sig$log2FoldChange > 1,]

negative_PANC1_WT_vs_48 <- res_PANC1_WT_vs_48_sig[res_PANC1_WT_vs_48_sig$log2FoldChange < -1,]

nrow(positive_PANC1_WT_vs_48) + nrow(negative_PANC1_WT_vs_48)
#2409 total DE

#upregulated
top_genes <- c(rownames(positive_PANC1_WT_vs_48[order(positive_PANC1_WT_vs_48$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_PANC1_WT_vs_48))
#downregulated
top_genes <- c(rownames(negative_PANC1_WT_vs_48[order(negative_PANC1_WT_vs_48$log2FoldChange, decreasing = F),]),
               rownames(positive_PANC1_WT_vs_48))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(, 
  log10(select[,order(df$time)] + 1)[1:30,], 
  cluster_rows=F, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T, main = "PANC1 cell line, WT vs 48h, up@WT") 


###p53 target genes###
#upregulated at 48
p53_PANC1_WT_48_at_48 <- intersect(rownames(positive_PANC1_WT_vs_48), target_list)
#51
#up at WT
p53_PANC1_WT_48_at_WT <- intersect(rownames(negative_PANC1_WT_vs_48), target_list)
#2


######modifed volcano####

p53_PANC1_48 <- union(p53_PANC1_WT_48_at_48, p53_PANC1_WT_48_at_WT)

res_PANC1_WT_vs_48$keyvals <- ""
res_PANC1_WT_vs_48$keyvals <- ifelse(res_PANC1_WT_vs_48$log2FoldChange < -1 & res_PANC1_WT_vs_48$padj < 0.05, 'skyblue3', res_PANC1_WT_vs_48$keyvals)
res_PANC1_WT_vs_48$keyvals <- ifelse(res_PANC1_WT_vs_48$log2FoldChange > 1 & res_PANC1_WT_vs_48$padj < 0.05, 'lightcoral',res_PANC1_WT_vs_48$keyvals)
res_PANC1_WT_vs_48$keyvals <- ifelse(rownames(res_PANC1_WT_vs_48) %in% p53_PANC1_48, 'springgreen4', res_PANC1_WT_vs_48$keyvals)
res_PANC1_WT_vs_48[res_PANC1_WT_vs_48$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_PANC1_WT_vs_48$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_PANC1_WT_vs_48[order(res_PANC1_WT_vs_48$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_PANC1_WT_vs_48[order(res_PANC1_WT_vs_48$keyvals, decreasing = F),],
                     lab = rownames(res_PANC1_WT_vs_48[order(res_PANC1_WT_vs_48$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'PANC1 48 vs WT',
                     caption = bquote("Total DE genes: 2409 (overlapping 53)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-7, 9)) +
  scale_x_continuous(breaks = seq(-7,9,1))

ggsave("plots/volcanoes_overlapping_116_targets/PANC1_@_48.svg", plot = p, 
       width = 10, height = 8, dpi = 300, device = grDevices::svg)


#####PANC1 at 72h####
colData_pair <- colData[colData$colf == "PANC1" & (colData$time == 72 | colData$time == "WT"),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~time)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#32330

res_PANC1_WT_vs_72 <- results(dds_pair, contrast = c("time", "72", "WT"))
res_PANC1_WT_vs_72 <- data.frame(res_PANC1_WT_vs_72)
res_PANC1_WT_vs_72 <- res_PANC1_WT_vs_72[complete.cases(res_PANC1_WT_vs_72),]

res_PANC1_WT_vs_72_sig <- res_PANC1_WT_vs_72[res_PANC1_WT_vs_72$padj < 0.05,]
positive_PANC1_WT_vs_72 <- res_PANC1_WT_vs_72_sig[res_PANC1_WT_vs_72_sig$log2FoldChange > 1,]

negative_PANC1_WT_vs_72 <- res_PANC1_WT_vs_72_sig[res_PANC1_WT_vs_72_sig$log2FoldChange < -1,]

nrow(positive_PANC1_WT_vs_72) + nrow(negative_PANC1_WT_vs_72)
#1472

#upregulated
top_genes <- c(rownames(positive_PANC1_WT_vs_72[order(positive_PANC1_WT_vs_72$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_PANC1_WT_vs_72))
#downregulated
top_genes <- c(rownames(negative_PANC1_WT_vs_72[order(negative_PANC1_WT_vs_72$log2FoldChange, decreasing = F),]),
               rownames(positive_PANC1_WT_vs_72))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(, 
  log10(select[,order(df$time)] + 1)[1:30,], 
  cluster_rows=F, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T, main = "PANC1 cell line, WT vs 72h, up@72") 


###p53 target genes###
#upregulated at 72
p53_PANC1_WT_72_at_72 <- intersect(rownames(positive_PANC1_WT_vs_72), target_list)
#79
#up at WT
p53_PANC1_WT_72_at_WT <- intersect(rownames(negative_PANC1_WT_vs_72), target_list)
#none

####modifed volcano####
p53_PANC1_72 <- union(p53_PANC1_WT_72_at_72, p53_PANC1_WT_72_at_WT)
#79

res_PANC1_WT_vs_72$keyvals <- ""
res_PANC1_WT_vs_72$keyvals <- ifelse(res_PANC1_WT_vs_72$log2FoldChange < -1 & res_PANC1_WT_vs_72$padj < 0.05, 'skyblue3', res_PANC1_WT_vs_72$keyvals)
res_PANC1_WT_vs_72$keyvals <- ifelse(res_PANC1_WT_vs_72$log2FoldChange > 1 & res_PANC1_WT_vs_72$padj < 0.05, 'lightcoral',res_PANC1_WT_vs_72$keyvals)
res_PANC1_WT_vs_72$keyvals <- ifelse(rownames(res_PANC1_WT_vs_72) %in% p53_PANC1_72, 'springgreen4', res_PANC1_WT_vs_72$keyvals)
res_PANC1_WT_vs_72[res_PANC1_WT_vs_72$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_PANC1_WT_vs_72$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_PANC1_WT_vs_72[order(res_PANC1_WT_vs_72$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_PANC1_WT_vs_72[order(res_PANC1_WT_vs_72$keyvals, decreasing = F),],
                     lab = rownames(res_PANC1_WT_vs_72[order(res_PANC1_WT_vs_72$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     axisLabSize = 30,
                     captionLabSize = 24,
                     legendLabSize = 20,
                     legendIconSize = 8,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 10,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'PANC-1, 72h p.i.',
                     subtitle = '',
                     caption = bquote("Total DE genes: 1472 (overlapping 79)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-5, 10)) +
  scale_x_continuous(breaks = seq(-5,10,1))
ggsave("plots/volcanoes_overlapping_116_targets/PANC1_@_72.svg", plot = p, 
       width = 14, height = 10, dpi = 300, device = grDevices::svg)




####a heatmap counts across all timepoints####
#heatmap for RNASeq data
tmp <- union(rownames(head(positive_PANC1_WT_vs_36[order(positive_PANC1_WT_vs_36$padj),], 20)),
             rownames(head(positive_PANC1_WT_vs_48[order(positive_PANC1_WT_vs_48$padj),]),20))
PANC1_up <- as.data.frame(counts(dds, normalized=T))
PANC1_up <- PANC1_up[union(tmp,
                           rownames(head(positive_PANC1_WT_vs_72[order(positive_PANC1_WT_vs_72$padj),], 20))),
                     colData[colData$colf == "PANC1" & colData$time != "WT", "col"]]
df <- data.frame(time = colData[colData$colf == "PANC1" & colData$time != "WT", "time"])
rownames(df) <- colnames(PANC1_up)

##z-scores
PANC1_up_z <-t(apply(PANC1_up, 1, scale))
colnames(PANC1_up_z) <- colnames(PANC1_up)

rb <- data.frame(c("36" = rownames(head(positive_PANC1_WT_vs_36[order(positive_PANC1_WT_vs_36$padj),], 20)),
                   "48" = rownames(head(positive_PANC1_WT_vs_48[order(positive_PANC1_WT_vs_48$padj),]),20),
                   "72" = rownames(head(positive_PANC1_WT_vs_72[order(positive_PANC1_WT_vs_72$padj),], 20))))
##

PANC1_up$PANC1_WT <- NA
df <- microeco::dropallfactors(df)
df[10,] <- c("0")
df$time <- df[order(df$time),]
rownames(df)[rownames(df) == 10] <- "PANC1_WT"

ann_colors <- list(time = c("36" = "gray70", "48" = "gray40", "72" = "gray5", "0" = "white"))
p <- pheatmap( 
  log10(PANC1_up[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F,
  cellwidth = 90,
  legend_breaks = c(0:4),
  gaps_col = c(1,4,7),
  fontsize = 28,
  annotation_colors = ann_colors,
  labels_col = c("WT", paste0("36_hpi_", 1:3), paste0("48_hpi_", 1:3), paste0("72_hpi_", 1:3)),
  main = "PANC1",
  na_col = "steelblue4")
ggsave("plots/heatmaps/PANC1_upregulated.svg", plot = p, 
       width = 17, height = 18, dpi = 300, device = grDevices::svg, units = "in")


#z-score heatmap
PANC1_up_z <- as.data.frame(PANC1_up_z)
PANC1_up_z$PANC1_WT <- 0
pheatmap( 
  log10(PANC1_up_z[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T,
  main = "upregulated genes in PANC1, z scores",
  na_col = "steelblue4")




######A431 at 36h#####

colData_pair <- colData[(colData$colf == "A_273" | colData$colf == "A_VEG") & (colData$time == 36),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#27490

res_A431_36 <- results(dds_pair, contrast = c("colf", "A_273", "A_VEG"))
res_A431_36 <- data.frame(res_A431_36)
res_A431_36 <- res_A431_36[complete.cases(res_A431_36),]

res_A431_36_sig <- res_A431_36[res_A431_36$padj < 0.05,]
positive_A431_36 <- res_A431_36_sig[res_A431_36_sig$log2FoldChange > 1,]

negative_A431_36 <- res_A431_36_sig[res_A431_36_sig$log2FoldChange < -1,]

nrow(negative_A431_36) + nrow(positive_A431_36)
#11 total de

#upregulated
top_genes <- c(rownames(positive_A431_36[order(positive_A431_36$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_A431_36))
#downregulated
top_genes <- c(rownames(negative_A431_36[order(negative_A431_36$log2FoldChange, decreasing = F),]),
               rownames(positive_A431_36))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap( 
  log10(select[,order(df$time)] + 1), 
  cluster_rows=F, show_rownames=T,
  #cluster_rows=T, show_rownames=F,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "A431 cell line, 36h") 


###p53 target genes###
#upregulated at corrected
p53_A431_36 <- intersect(rownames(positive_A431_36), target_list)
#7

#up at WT
p53_A431_VEGF <- intersect(rownames(negative_A431_36), target_list)
#none

######modifed volcano####
p53_A431_36_c <- union(p53_A431_36, p53_A431_VEGF)
#7

res_A431_36$keyvals <- ""
res_A431_36$keyvals <- ifelse(res_A431_36$log2FoldChange < -1 & res_A431_36$padj < 0.05, 'skyblue3', res_A431_36$keyvals)
res_A431_36$keyvals <- ifelse(res_A431_36$log2FoldChange > 1 & res_A431_36$padj < 0.05, 'lightcoral',res_A431_36$keyvals)
res_A431_36$keyvals <- ifelse(rownames(res_A431_36) %in% p53_A431_36_c, 'springgreen4', res_A431_36$keyvals)
res_A431_36[res_A431_36$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_A431_36$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_A431_36[order(res_A431_36$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_A431_36[order(res_A431_36$keyvals, decreasing = F),],
                     lab = rownames(res_A431_36[order(res_A431_36$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'A431 at 36h',
                     caption = bquote("Total DE genes: 11 (overlapping 7)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)

p <- p + coord_cartesian(xlim = c(-5, 6)) +
  scale_x_continuous(breaks = seq(-5,6,1))

ggsave("plots/volcanoes_overlapping_116_targets/A431_@_36.svg", plot = p, 
       width = 10, height = 8, dpi = 300, device = grDevices::svg)


#####A431 at 48h#####

colData_pair <- colData[(colData$colf == "A_273" | colData$colf == "A_VEG") & (colData$time == 48),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#27665

res_A431_48 <- results(dds_pair, contrast = c("colf", "A_273", "A_VEG"))
res_A431_48 <- data.frame(res_A431_48)
res_A431_48 <- res_A431_48[complete.cases(res_A431_48),]

res_A431_48_sig <- res_A431_48[res_A431_48$padj < 0.05,]

positive_A431_48 <- res_A431_48_sig[res_A431_48_sig$log2FoldChange > 1,]

negative_A431_48 <- res_A431_48_sig[res_A431_48_sig$log2FoldChange < -1,]

nrow(positive_A431_48) + nrow(negative_A431_48)
#55 total de

#upregulated
top_genes <- c(rownames(positive_A431_48[order(positive_A431_48$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_A431_48))
#downregulated
top_genes <- c(rownames(negative_A431_48[order(negative_A431_48$log2FoldChange, decreasing = F),]),
               rownames(positive_A431_48))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap( 
  log10(select[,order(df$time)] + 1), 
  cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "A431 cell line, 48h") 


###p53 target genes###
#upregulated at corrected
p53_A431_48 <- intersect(rownames(positive_A431_48), target_list)
#22

#up at WT
p53_A431_VEGF <- intersect(rownames(negative_A431_48), target_list)
#none

####modifed volcano####
p53_A431_48_c <- union(p53_A431_48, p53_A431_VEGF)
#22


res_A431_48$keyvals <- ""
res_A431_48$keyvals <- ifelse(res_A431_48$log2FoldChange < -1 & res_A431_48$padj < 0.05, 'skyblue3', res_A431_48$keyvals)
res_A431_48$keyvals <- ifelse(res_A431_48$log2FoldChange > 1 & res_A431_48$padj < 0.05, 'lightcoral',res_A431_48$keyvals)
res_A431_48$keyvals <- ifelse(rownames(res_A431_48) %in% p53_A431_48_c, 'springgreen4', res_A431_48$keyvals)
res_A431_48[res_A431_48$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_A431_48$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_A431_48[order(res_A431_48$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_A431_48[order(res_A431_48$keyvals, decreasing = F),],
                     lab = rownames(res_A431_48[order(res_A431_48$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 10,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'A431 at 48h',
                     caption = bquote("Total DE genes: 55 (overlapping 22)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)

p <- p + coord_cartesian(xlim = c(-1.5, 6.5)) +
  scale_x_continuous(breaks = seq(-2,6.5,1))

ggsave("plots/volcanoes_overlapping_116_targets/A431_@_48.svg", plot = p, 
       width = 10, height = 8, dpi = 300, device = grDevices::svg)


#####A431 at 72h####

colData_pair <- colData[(colData$colf == "A_273" | colData$colf == "A_VEG") & (colData$time == 72),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#29278

res_A431_72 <- results(dds_pair, contrast = c("colf", "A_273", "A_VEG"))
res_A431_72 <- data.frame(res_A431_72)
res_A431_72 <- res_A431_72[complete.cases(res_A431_72),]

res_A431_72_sig <- res_A431_72[res_A431_72$padj < 0.05,]

positive_A431_72 <- res_A431_72_sig[res_A431_72_sig$log2FoldChange > 1,]

negative_A431_72 <- res_A431_72_sig[res_A431_72_sig$log2FoldChange < -1,]

nrow(positive_A431_72) + nrow(negative_A431_72)
#327

#upregulated
top_genes <- c(rownames(positive_A431_72[order(positive_A431_72$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_A431_72))
#downregulated
top_genes <- c(rownames(negative_A431_72[order(negative_A431_72$log2FoldChange, decreasing = F),]),
               rownames(positive_A431_72))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#jannotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap( 
  log10(select[,order(df$time)] + 1)[1:30,], 
  cluster_rows=T, show_rownames=T,
  #cluster_rows=T, show_rownames=F,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "A431 cell line, 72h, up @ 72h") 


###p53 target genes###
#upregulated at corrected
p53_A431_72 <- intersect(rownames(positive_A431_72), target_list)
#56

#up at corrected
p53_A431_72_VEGF <- intersect(rownames(negative_A431_72), target_list)
#none

######modifed volcano####
p53_A431_72_c <- union(p53_A431_72, p53_A431_72_VEGF)
#56

res_A431_72$keyvals <- ""
res_A431_72$keyvals <- ifelse(res_A431_72$log2FoldChange < -1 & res_A431_72$padj < 0.05, 'skyblue3', res_A431_72$keyvals)
res_A431_72$keyvals <- ifelse(res_A431_72$log2FoldChange > 1 & res_A431_72$padj < 0.05, 'lightcoral',res_A431_72$keyvals)
res_A431_72$keyvals <- ifelse(rownames(res_A431_72) %in% p53_A431_72_c, 'springgreen4', res_A431_72$keyvals)
res_A431_72[res_A431_72$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_A431_72$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_A431_72[order(res_A431_72$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_A431_72[order(res_A431_72$keyvals, decreasing = F),],
                     lab = rownames(res_A431_72[order(res_A431_72$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 10,
                     labFace = 'italic',
                     axisLabSize = 30,
                     captionLabSize = 24,
                     legendLabSize = 20,
                     legendIconSize = 8,
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'A431, 72h p.i.',
                     subtitle = '',
                     caption = bquote("Total DE genes: 327 (overlapping 56)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-2.5, 10)) +
  scale_x_continuous(breaks = seq(-3,10,1))

ggsave("plots/volcanoes_overlapping_116_targets/A431_@_72.svg", plot = p, 
       width = 14, height = 10, dpi = 300, device = grDevices::svg)


#heatmap for p53 target genes across all the timepoints
p53_up_at_test <- as.data.frame(counts(dds, normalized=T))
tmp <- union(p53_A431_36, p53_A431_48)
p53_up_at_test <- p53_up_at_test[union(p53_A431_72, tmp),
                                 colData[colData$colf == "A_273" & colData$time != "WT", "col"]]
df <- data.frame(time = colData[colData$colf == "A_273" & colData$time != "WT", "time"])
rownames(df) <- colnames(p53_up_at_test)

#z score
p53_up_at_test_z <-t(apply(p53_up_at_test, 1, scale))
colnames(p53_up_at_test_z) <- colnames(p53_up_at_test)

p53_up_at_test$A431_VEGF <- NA
df <- microeco::dropallfactors(df)
df[7,] <- c("VEGF")
rownames(df)[rownames(df) == 7] <- "A431_VEGF"


ann_colors <- list(time = c("36" = "gray70", "48" = "gray40", "72" = "gray5", "VEGF" = "white"))
h <- pheatmap(#a[,order(df$time)], 
  log10(p53_up_at_test[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T,
  na_col = "steelblue4",
  cellwidth = 120,
  #cellheight = 25,
  fontsize = 20,
  legend_breaks = 0:4,
  annotation_colors = ann_colors,
  labels_col = c(paste0("36_hpi_", 1:3), paste0("48_hpi_", 1:3), paste0("72_hpi_", 1:3), "VEGF"),
  main = "A431")

ggsave("plots/heatmaps/A431__p53_1.svg", plot = h,
       dpi = 300, width = 20, height = 26, units = "in", grDevices::svg)


#heatmap for RNASeq data
tmp <- union(rownames(head(positive_A431_36[order(positive_A431_36$padj),], 20)),
             rownames(head(positive_A431_48[order(positive_A431_48$padj),]),20))
A431_up <- as.data.frame(counts(dds, normalized=T))
A431_up <- A431_up[union(tmp,
                         rownames(head(positive_A431_72[order(positive_A431_72$padj),], 40))),
                   colData[colData$colf == "A_273", "col"]]
df <- data.frame(time = colData[colData$colf == "A_273", "time"])
rownames(df) <- colnames(A431_up)

##z-scores
A431_up_z <-t(apply(A431_up, 1, scale))
colnames(A431_up_z) <- colnames(A431_up)

rb_A431 <- data.frame(c("36" = rownames(head(positive_A431_36[order(positive_A431_36$padj),], 20)),
                        "48" = rownames(head(positive_A431_48[order(positive_A431_48$padj),]),20),
                        "72" = rownames(head(positive_A431_72[order(positive_A431_72$padj),], 40))))
##

A431_up$A431_VEGF <- NA
df$time <- droplevels(df$time)
df <- microeco::dropallfactors(df)
df[7,] <- c("0")
rownames(df)[rownames(df) == 7] <- "A431_VEGF"

ann_colors <- list(time = c("36" = "gray70", "48" = "gray40", "72" = "gray5", "0" = "white"))
p <- pheatmap( 
  log10(A431_up[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F,
  main = "A431",
  gaps_col = c(1,3,5),
  fontsize = 28,
  cellwidth = 120,
  legend_breaks = 1:4,
  annotation_colors = ann_colors,
  labels_col = c("VEGF", paste0("36_hpi_", 1:2), paste0("48_hpi_", 1:2), paste0("72_hpi_", 1:2)),
  na_col = "steelblue4")
ggsave("plots/heatmaps/A431_upregulated.svg", plot = p, 
       width = 17, height = 18, dpi = 300, device = grDevices::svg, units = "in")

#z-score heatmap
A431_up_z <- as.data.frame(A431_up_z)
A431_up_z$A431_VEGF <- 0
pheatmap( 
  log10(A431_up_z[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T,
  main = "upregulated genes in A431, z scores",
  na_col = "steelblue4")



#####HT29 at 36h####

colData_pair <- colData[(colData$colf == "H_273" | colData$colf == "H_VEG") & (colData$time == 36),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#28619

res_HT29_36 <- results(dds_pair, contrast = c("colf", "H_273", "H_VEG"))
res_HT29_36 <- data.frame(res_HT29_36)
res_HT29_36 <- res_HT29_36[complete.cases(res_HT29_36),]

res_HT29_36_sig <- res_HT29_36[res_HT29_36$padj < 0.05,]

positive_HT29_36 <- res_HT29_36_sig[res_HT29_36_sig$log2FoldChange > 1,]

negative_HT29_36 <- res_HT29_36_sig[res_HT29_36_sig$log2FoldChange < -1,]

nrow(positive_HT29_36) + nrow(negative_HT29_36)
#1

#upregulated
top_genes <- c(rownames(positive_HT29_36[order(positive_HT29_36$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_HT29_36))
#downregulated
top_genes <- c(rownames(negative_HT29_36[order(negative_HT29_36$log2FoldChange, decreasing = F),]),
               rownames(positive_HT29_36))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(
  log10(select[,order(df$time)] + 1), 
  cluster_rows=F, show_rownames=T,
  #cluster_rows=T, show_rownames=F,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "HT29 cell line, 36h") 

###p53 target genes###
#upregulated at corrected
p53_HT29_36 <- intersect(rownames(positive_HT29_36), target_list)
#1
#up at WT
p53_HT29_36_VEGF <- intersect(rownames(negative_HT29_36), target_list)
#none

######modifed volcano####
p53_HT29_36_c <- union(p53_HT29_36, p53_HT29_36_VEGF)
#1

res_HT29_36$keyvals <- ""
res_HT29_36$keyvals <- ifelse(res_HT29_36$log2FoldChange < -1 & res_HT29_36$padj < 0.05, 'skyblue3', res_HT29_36$keyvals)
res_HT29_36$keyvals <- ifelse(res_HT29_36$log2FoldChange > 1 & res_HT29_36$padj < 0.05, 'lightcoral',res_HT29_36$keyvals)
res_HT29_36$keyvals <- ifelse(rownames(res_HT29_36) %in% p53_HT29_36_c, 'springgreen4', res_HT29_36$keyvals)
res_HT29_36[res_HT29_36$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_HT29_36$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_HT29_36[order(res_HT29_36$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_HT29_36[order(res_HT29_36$keyvals, decreasing = F),],
                     lab = rownames(res_HT29_36[order(res_HT29_36$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 20,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'HT29 at 36h',
                     caption = bquote("Total DE genes: 1 (overlapping 1)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-5, 5)) +
  scale_x_continuous(breaks = seq(-5,5,1))

ggsave("plots/volcanoes_overlapping_116_targets/HT29_@_36.svg", plot = p, 
       width = 8, height = 8, dpi = 300, device = grDevices::svg)


#####HT29 at 48h####

colData_pair <- colData[(colData$colf == "H_273" | colData$colf == "H_VEG") & (colData$time == 48),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#28720

res_HT29_48 <- results(dds_pair, contrast = c("colf", "H_273", "H_VEG"))
res_HT29_48 <- data.frame(res_HT29_48)
res_HT29_48 <- res_HT29_48[complete.cases(res_HT29_48),]

res_HT29_48_sig <- res_HT29_48[res_HT29_48$padj < 0.05,]

positive_HT29_48 <- res_HT29_48_sig[res_HT29_48_sig$log2FoldChange > 1,]

negative_HT29_48 <- res_HT29_48_sig[res_HT29_48_sig$log2FoldChange < -1,]

nrow(positive_HT29_48) + nrow(negative_HT29_48)
#54

#upregulated
top_genes <- c(rownames(positive_HT29_48[order(positive_HT29_48$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_HT29_48))
#downregulated
top_genes <- c(rownames(negative_HT29_48[order(negative_HT29_48$log2FoldChange, decreasing = F),]),
               rownames(positive_HT29_48))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#just having time as col annotation
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(, 
  log10(select[,order(df$time)] + 1), 
  cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "HT29 cell line, 48h, up @ corrected") 


###p53 target genes###
#upregulated at corrected
p53_HT29_48 <- intersect(rownames(positive_HT29_48), target_list)
#31
#up at VEGF
p53_HT29_48_VEGF <- intersect(rownames(negative_HT29_48), target_list)
#none

######modifed volcano####
p53_HT29_48_c <- union(p53_HT29_48, p53_HT29_48_VEGF)
#31

res_HT29_48$keyvals <- ""
res_HT29_48$keyvals <- ifelse(res_HT29_48$log2FoldChange < -1 & res_HT29_48$padj < 0.05, 'skyblue3', res_HT29_48$keyvals)
res_HT29_48$keyvals <- ifelse(res_HT29_48$log2FoldChange > 1 & res_HT29_48$padj < 0.05, 'lightcoral',res_HT29_48$keyvals)
res_HT29_48$keyvals <- ifelse(rownames(res_HT29_48) %in% p53_HT29_48_c, 'springgreen4', res_HT29_48$keyvals)
res_HT29_48[res_HT29_48$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_HT29_48$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_HT29_48[order(res_HT29_48$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_HT29_48[order(res_HT29_48$keyvals, decreasing = F),],
                     lab = rownames(res_HT29_48[order(res_HT29_48$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'HT29 at 48h',
                     caption = bquote("Total DE genes: 54 (overlapping 31)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)

p <- p + coord_cartesian(xlim = c(-4, 8.5)) +
  scale_x_continuous(breaks = seq(-4,9,1))
ggsave("plots/volcanoes_overlapping_116_targets/HT29_@_48.svg", plot = p, 
       width = 10, height = 8, dpi = 300, device = grDevices::svg)



#####HT29 at 72h####

colData_pair <- colData[(colData$colf == "H_273" | colData$colf == "H_VEG") & (colData$time == 72),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#30309

res_HT29_72 <- results(dds_pair, contrast = c("colf", "H_273", "H_VEG"))
res_HT29_72 <- data.frame(res_HT29_72)
res_HT29_72 <- res_HT29_72[complete.cases(res_HT29_72),]

res_HT29_72_sig <- res_HT29_72[res_HT29_72$padj < 0.05,]

positive_HT29_72 <- res_HT29_72_sig[res_HT29_72_sig$log2FoldChange > 1,]

negative_HT29_72 <- res_HT29_72_sig[res_HT29_72_sig$log2FoldChange < -1,]

nrow(positive_HT29_72) + nrow(negative_HT29_72)
#260

#upregulated
top_genes <- c(rownames(positive_HT29_72[order(positive_HT29_72$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_HT29_72))
#downregulated
top_genes <- c(rownames(negative_HT29_72[order(negative_HT29_72$log2FoldChange, decreasing = F),]),
               rownames(positive_HT29_72))

select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(, 
  log10(select[,order(df$time)] + 1)[1:30,], 
  cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "HT29 cell line, 72h, up @ corrected") 


###p53 target genes###
#upregulated at corrected
p53_HT29_72 <- intersect(rownames(positive_HT29_72), target_list)
#70
#up at WT
p53_HT29_72_VEGF <- intersect(rownames(negative_HT29_72), target_list)
#none

######modifed volcano####
p53_HT29_72_c <- union(p53_HT29_72, p53_HT29_72_VEGF)
#70

res_HT29_72$keyvals <- ""
res_HT29_72$keyvals <- ifelse(res_HT29_72$log2FoldChange < -1 & res_HT29_72$padj < 0.05, 'skyblue3', res_HT29_72$keyvals)
res_HT29_72$keyvals <- ifelse(res_HT29_72$log2FoldChange > 1 & res_HT29_72$padj < 0.05, 'lightcoral',res_HT29_72$keyvals)
res_HT29_72$keyvals <- ifelse(rownames(res_HT29_72) %in% p53_HT29_72_c, 'springgreen4', res_HT29_72$keyvals)
res_HT29_72[res_HT29_72$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_HT29_72$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
names(keyvals)[keyvals == 'springgreen4'] <- 'High confidence p53 targets'
tmp <- res_HT29_72[order(res_HT29_72$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])

p <- EnhancedVolcano(res_HT29_72[order(res_HT29_72$keyvals, decreasing = F),],
                     lab = rownames(res_HT29_72[order(res_HT29_72$keyvals, decreasing = F),]),
                     selectLab = rownames(overlap_labels),
                     drawConnectors = T,
                     maxoverlapsConnectors = 10,
                     colCustom = (keyvals[order(keyvals, decreasing = F)]),
                     colAlpha = 0.8,
                     labSize = 10,
                     labFace = 'italic',
                     axisLabSize = 30,
                     captionLabSize = 24,
                     legendLabSize = 20,
                     legendIconSize = 8,
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'HT-29, 72h p.i.',
                     subtitle = '',
                     caption = bquote("Total DE genes: 260 (overlapping 70)"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-3, 11.5)) +
  scale_x_continuous(breaks = seq(-3,12,1))
ggsave("plots/volcanoes_overlapping_116_targets/HT29_@_72.svg", plot = p, 
       width = 14, height = 10, dpi = 300, device = grDevices::svg, units = "in")


#heatmap for RNASeq data
tmp <- union(rownames(head(positive_HT29_36[order(positive_HT29_36$padj),], 20)),
             rownames(head(positive_HT29_48[order(positive_HT29_48$padj),]),20))
HT29_up <- as.data.frame(counts(dds, normalized=T))
HT29_up <- HT29_up[union(tmp,
                         rownames(head(positive_HT29_72[order(positive_HT29_72$padj),], 40))),
                   colData[colData$colf == "H_273", "col"]]
df <- data.frame(time = colData[colData$colf == "H_273", "time"])
rownames(df) <- colnames(HT29_up)

##z-scores
HT29_up_z <-t(apply(HT29_up, 1, scale))
colnames(HT29_up_z) <- colnames(HT29_up)

rb_ht29 <- data.frame(c("36" = rownames(head(positive_HT29_36[order(positive_HT29_36$padj),], 20)),
                        "48" = rownames(head(positive_HT29_48[order(positive_HT29_48$padj),]),20),
                        "72" = rownames(head(positive_HT29_72[order(positive_HT29_72$padj),], 40))))
##

HT29_up$HT29_VEGF <- NA
df$time <- droplevels(df$time)
df <- microeco::dropallfactors(df)
df[7,] <- c("0")
rownames(df)[rownames(df) == 7] <- "HT29_VEGF"

ann_colors <- list(time = c("36" = "gray70", "48" = "gray40", "72" = "gray5", "0" = "white"))
p <- pheatmap(
  log10(HT29_up[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F,
  main = "HT-29",
  annotation_colors = ann_colors,
  gaps_col = c(1,3,5),
  fontsize = 28,
  cellwidth = 115,
  labels_col = c("VEGF", paste0("36_hpi_", 1:2), paste0("48_hpi_", 1:2), paste0("72_hpi_", 1:2)),
  na_col = "steelblue4")
ggsave("plots/heatmaps/HT29_upregulated.svg", plot = p, 
       width = 17, height = 18, dpi = 300, device = grDevices::svg, units = "in")

#z-score heatmap
HT29_up_z <- as.data.frame(HT29_up_z)
HT29_up_z$HT29_VEGF <- 0
pheatmap(#a[,order(df$time)], 
  log10(HT29_up_z[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T,
  main = "upregulated genes in HT29, z scores",
  na_col = "steelblue4")


####H311 (SMAD4) at 36h####
colData_pair <- colData[(colData$colf == "H_311" | colData$colf == "H_VEG") & (colData$time == 36),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#28761

res_H311_36 <- results(dds_pair, contrast = c("colf", "H_311", "H_VEG"))
res_H311_36 <- data.frame(res_H311_36)
res_H311_36 <- res_H311_36[complete.cases(res_H311_36),]

res_H311_36_sig <- res_H311_36[res_H311_36$padj < 0.05,]
positive_H311_36 <- res_H311_36_sig[res_H311_36_sig$log2FoldChange > 1,]

negative_H311_36 <- res_H311_36_sig[res_H311_36_sig$log2FoldChange < -1,]

nrow(positive_H311_36) + nrow(negative_H311_36)
#5

#upregulated
top_genes <- c(rownames(positive_H311_36[order(positive_H311_36$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_H311_36))
#downregulated
top_genes <- c(rownames(negative_H311_36[order(negative_H311_36$log2FoldChange, decreasing = F),]),
               rownames(positive_H311_36))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap( 
  log10(select[,order(df$time)] + 1), 
  cluster_rows=F, show_rownames=T,
  #cluster_rows=T, show_rownames=F,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "SMAD4 mutation, 36h") 

######modified volcano####
res_H311_36$keyvals <- ""
res_H311_36$keyvals <- ifelse(res_H311_36$log2FoldChange < -1 & res_H311_36$padj < 0.05, 'skyblue3', res_H311_36$keyvals)
res_H311_36$keyvals <- ifelse(res_H311_36$log2FoldChange > 1 & res_H311_36$padj < 0.05, 'lightcoral',res_H311_36$keyvals)
res_H311_36[res_H311_36$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_H311_36$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
tmp <- res_H311_36[order(res_H311_36$keyvals, decreasing = F),]
overlap_labels <- (tmp[ tmp$keyvals == 'springgreen4', ])


p <- EnhancedVolcano(res_H311_36,
                     lab = rownames(res_H311_36),
                     x = 'log2FoldChange',
                     y = 'padj',
                     ylab = bquote('p-value'),
                     drawConnectors=T,
                     maxoverlapsConnectors = 10,
                     colCustom = keyvals,
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     title = 'SMAD4 at 36h',
                     caption = bquote("Total DE genes: 5"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-4, 5)) +
  scale_x_continuous(breaks = seq(-4,5,1))
ggsave("plots/volcanoes/SMAD4_@_36.svg", plot = p, 
       width = 12, height = 8, dpi = 300, device = grDevices::svg, units = "in")


#####H311 at 48h####
colData_pair <- colData[(colData$colf == "H_311" | colData$colf == "H_VEG") & (colData$time == 48),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#29133

res_H311_48 <- results(dds_pair, contrast = c("colf", "H_311", "H_VEG"))
res_H311_48 <- data.frame(res_H311_48)
res_H311_48 <- res_H311_48[complete.cases(res_H311_48),]

res_H311_48_sig <- res_H311_48[res_H311_48$padj < 0.05,]
positive_H311_48 <- res_H311_48_sig[res_H311_48_sig$log2FoldChange > 1,]

negative_H311_48 <- res_H311_48_sig[res_H311_48_sig$log2FoldChange < -1,]

nrow(positive_H311_48) + nrow(negative_H311_48)
#14

#upregulated
top_genes <- c(rownames(positive_H311_48[order(positive_H311_48$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_H311_48))
#downregulated
top_genes <- c(rownames(negative_H311_48[order(negative_H311_48$log2FoldChange, decreasing = F),]),
               rownames(positive_H311_48))


select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap(
  log10(select[,order(df$time)] + 1), 
  cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "SMAD4 mutation, 48h") 


######modified volcano####
res_H311_48$keyvals <- ""
res_H311_48$keyvals <- ifelse(res_H311_48$log2FoldChange < -1 & res_H311_48$padj < 0.05, 'skyblue3', res_H311_48$keyvals)
res_H311_48$keyvals <- ifelse(res_H311_48$log2FoldChange > 1 & res_H311_48$padj < 0.05, 'lightcoral',res_H311_48$keyvals)
res_H311_48[res_H311_48$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_H311_48$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
tmp <- res_H311_48[order(res_H311_48$keyvals, decreasing = F),]


p <- EnhancedVolcano(res_H311_48,
                     lab = rownames(res_H311_48),
                     x = 'log2FoldChange',
                     y = 'padj',
                     ylab = bquote('p-value'),
                     drawConnectors=T,
                     maxoverlapsConnectors = 10,
                     colCustom = keyvals,
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     title = 'SMAD4 at 48h',
                     caption = bquote("Total DE genes: 14"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p + coord_cartesian(xlim = c(-2, 7)) +
  scale_x_continuous(breaks = seq(-2,7,1))
ggsave("plots/volcanoes/SMAD4_@_48.svg", plot = p, 
       width = 10, height = 8, dpi = 300, device = grDevices::svg, units = "in")


#####H311 at 72h####
colData_pair <- colData[(colData$colf == "H_311" | colData$colf == "H_VEG") & (colData$time == 72),]
rownames(colData_pair) <- colData_pair$col
dds_pair <- DESeqDataSetFromMatrix(countData = count_matrix[,rownames(colData_pair)],
                                   colData = colData_pair,
                                   design= ~colf)
dds_pair <- DESeq(dds_pair)
nrow(counts(dds_pair))
#63140

dds_pair <- dds_pair[rowSums(counts(dds_pair)) > 0,]
nrow(counts(dds_pair))
#30894

res_H311_72 <- results(dds_pair, contrast = c("colf", "H_311", "H_VEG"))
res_H311_72 <- data.frame(res_H311_72)
res_H311_72 <- res_H311_72[complete.cases(res_H311_72),]

res_H311_72_sig <- res_H311_72[res_H311_72$padj < 0.05,]
positive_H311_72 <- res_H311_72_sig[res_H311_72_sig$log2FoldChange > 1,]

negative_H311_72 <- res_H311_72_sig[res_H311_72_sig$log2FoldChange < -1,]

nrow(positive_H311_72) + nrow(negative_H311_72)
#107

#upregulated
top_genes <- c(rownames(positive_H311_72[order(positive_H311_72$log2FoldChange, decreasing = TRUE),]), 
               rownames(negative_H311_72))
#downregulated
top_genes <- c(rownames(negative_H311_72[order(negative_H311_72$log2FoldChange, decreasing = F),]),
               rownames(positive_H311_72))

select <- counts(dds_pair, normalized = T)
select <- select[top_genes, rownames(colData_pair)]

#annotation for heatmap
df <- data.frame(time = colData(dds_pair)$time)
rownames(df) <- colnames(select)

#sort according to timepoints
pheatmap( 
  log10(select[,order(df$time)] + 1)[1:40,], 
  cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = F, main = "SMAD4 mutation, 72h, up @ VEGF") 

######modified volcano####
res_H311_72$keyvals <- ""
res_H311_72$keyvals <- ifelse(res_H311_72$log2FoldChange < -1 & res_H311_72$padj < 0.05, 'skyblue3', res_H311_72$keyvals)
res_H311_72$keyvals <- ifelse(res_H311_72$log2FoldChange > 1 & res_H311_72$padj < 0.05, 'lightcoral',res_H311_72$keyvals)
res_H311_72[res_H311_72$keyvals == "", "keyvals"] <- 'gray80'

keyvals <- res_H311_72$keyvals
names(keyvals)[keyvals == 'lightcoral'] <- 'sig. upregulated'
names(keyvals)[keyvals == 'skyblue3'] <- 'sig. downregulated'
names(keyvals)[keyvals == 'gray80'] <- 'not significant'
tmp <- res_H311_48[order(res_H311_48$keyvals, decreasing = F),]


p <- EnhancedVolcano(res_H311_72,
                     lab = rownames(res_H311_72),
                     x = 'log2FoldChange',
                     y = 'padj',
                     ylab = bquote('p-value'),
                     drawConnectors=T,
                     maxoverlapsConnectors = 20,
                     colCustom = keyvals,
                     colAlpha = 0.8,
                     labSize = 4,
                     labFace = 'italic',
                     title = 'SMAD4 at 72h',
                     caption = bquote("Total DE genes: 107"),
                     pCutoff = 0.05,
                     FCcutoff = 1)
p <- p + coord_cartesian(xlim = c(-3, 10)) +
  scale_x_continuous(breaks = seq(-3,10,1))
ggsave("plots/volcanoes/SMAD4_@_72.svg", plot = p, 
       width = 14, height = 10, dpi = 300, device = grDevices::svg, units = "in")

#heatmap for RNASeq data
tmp <- union(rownames(head(positive_H311_36[order(positive_H311_36$padj),], 20)),
             rownames(head(positive_H311_48[order(positive_H311_48$padj),]),20))
H311_up <- as.data.frame(counts(dds, normalized=T))
H311_up <- H311_up[union(rownames(head(positive_H311_48[order(positive_H311_48$padj),], 20)),
                         rownames(head(positive_H311_72[order(positive_H311_72$padj),], 40))),
                   colData[colData$colf == "H_311" & colData$time != 36, "col"]]
df <- data.frame(time = colData[colData$colf == "H_311" & colData$time != 36, "time"])
rownames(df) <- colnames(H311_up)

##z-scores
H311_up_z <-t(apply(H311_up, 1, scale))
colnames(H311_up_z) <- colnames(H311_up)

rb_H311 <- data.frame(c("36" = rownames(head(positive_H311_36[order(positive_H311_36$padj),], 20)),
                        "48" = rownames(head(positive_H311_48[order(positive_H311_48$padj),]),20),
                        "72" = rownames(head(positive_H311_72[order(positive_H311_72$padj),], 40))))
##

H311_up$H311_VEGF <- NA
df$time <- droplevels(df$time)
df <- microeco::dropallfactors(df)
df[5,] <- c("VEGF")
rownames(df)[rownames(df) == 5] <- "H311_VEGF"

ann_colors <- list(time = c("48" = "gray60", "72" = "gray5", "VEGF" = "white"))
p <- pheatmap( 
  log10(H311_up[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T,
  fontsize = 20,
  cellwidth = 110,
  main = "",
  #labels_col = c(paste0("36_hpi_", 1:2), paste0("48_hpi_", 1:2), paste0("72_hpi_", 1:2), "VEGF"),
  labels_col = c(paste0("48_hpi_", 1:2), paste0("72_hpi_", 1:2), "VEGF"),
  annotation_colors = ann_colors,
  legend_breaks = 1:4,
  na_col = "steelblue4")

ggsave("plots/heatmaps/SMAD4__48_72_upregulated.svg", plot = p, 
       width = 1, height = 18, dpi = 300, device = grDevices::svg, units = "in")

#z-score heatmap
H311_up_z <- as.data.frame(H311_up_z)
H311_up_z$H311_VEGF <- 0
pheatmap(#a[,order(df$time)], 
  log10(H311_up_z[,order(df$time)] + 1), cluster_rows=T, show_rownames=T,
  cluster_cols=F, annotation_col=df, annotation_legend = T,
  main = "upregulated genes in H311, z scores",
  na_col = "steelblue4")



####unique genes across cell lines####
##differentially regulated at 72h##

#unique genes in upregulated in A431, not present in PANC1
uniq_up_A431_PANC1 <- setdiff(rownames(positive_A431_72), rownames(positive_PANC1_WT_vs_72))
#unique genes in upregulated in A431, not present in HT29
uniq_up_A431 <- setdiff(uniq_up_A431_PANC1, rownames(positive_HT29_72))
#87

#unique genes in upregulated in HT29, not present in PANC1
uniq_up_HT29 <- setdiff(rownames(positive_HT29_72), rownames(positive_PANC1_WT_vs_72))
#unique genes in upregulated in HT29, not present in A431
uniq_up_HT29 <- setdiff(uniq_up_HT29, rownames(positive_A431_72))
#101

#unique genes in upregulated in PANC1, not present in A431
uniq_up_PANC1 <- setdiff(rownames(positive_PANC1_WT_vs_72), rownames(positive_A431_72))
#unique genes in upregulated in PANC1, not present in HT29
uniq_up_PANC1 <- setdiff(uniq_up_PANC1, rownames(positive_HT29_72))
#753

write_xlsx(data.frame(uniq_up_A431) ,"unique_up_A431.xlsx")
write_xlsx(data.frame(uniq_up_HT29) ,"unique_up_HT29.xlsx")
write_xlsx(data.frame(uniq_up_PANC1) ,"unique_up_PANC1.xlsx")

#unique genes in downregulated in A431, not present in PANC1
uniq_down_A431_PANC1 <- setdiff(rownames(negative_A431_72), rownames(negative_PANC1_WT_vs_72))
#unique genes in upregulated in A431, not present in HT29
uniq_down_A431 <- setdiff(uniq_down_A431_PANC1, rownames(negative_HT29_72))
#23

#unique genes in downregulated in HT29, not present in PANC1
uniq_down_HT29 <- setdiff(rownames(negative_HT29_72), rownames(negative_PANC1_WT_vs_72))
#unique genes in downregulated in HT29, not present in A431
uniq_down_HT29 <- setdiff(uniq_down_HT29, rownames(negative_A431_72))
#3

#unique genes in downregulated in PANC1, not present in A431
uniq_down_PANC1 <- setdiff(rownames(negative_PANC1_WT_vs_72), rownames(negative_A431_72))
#unique genes in downregulated in PANC1, not present in HT29
uniq_down_PANC1 <- setdiff(uniq_down_PANC1, rownames(negative_HT29_72))
#474

write_xlsx(data.frame(uniq_down_A431) ,"unique_down_A431.xlsx")
write_xlsx(data.frame(uniq_down_HT29) ,"unique_down_HT29.xlsx")
write_xlsx(data.frame(uniq_down_PANC1) ,"unique_down_PANC1.xlsx")

#####common genes across all the cell lines####
###now comparing the corrected vs the wild types at each timepoint###

#A273 vs all
#all DE genes in A273

up_72_A431_PANC1 <- intersect(rownames(positive_A431_72), rownames(positive_PANC1_WT_vs_72))
#94
write.csv(up_72_A431_PANC1, "common_up_A431_PANC1.csv")
write_xlsx(as.data.frame(up_72_A431_PANC1), "common_up_A431_PANC1.xlsx")

up_72_A431_HT29 <- intersect(rownames(positive_A431_72), rownames(positive_HT29_72))
#95
write.csv(up_72_A431_HT29, "common_up_A431_HT29.csv")
write_xlsx(as.data.frame(up_72_A431_HT29), "common_up_A431_HT29.xlsx")

up_72_A431_all <- intersect(up_72_A431_PANC1, up_72_A431_HT29)
#63
write.csv(up_72_A431_all, "common_up_A431_all.csv")
write_xlsx(as.data.frame(up_72_A431_all), "common_up_A431_all.xlsx")

down_72_A431_PANC1 <- intersect(rownames(negative_A431_72), rownames(negative_PANC1_WT_vs_72))
#91
write.csv(down_72_A431_PANC1, "common_down_A431_PANC1.csv")
write_xlsx(as.data.frame(down_72_A431_PANC1), "common_down_A431_PANC1.xlsx")

down_72_A431_HT29 <- intersect(rownames(negative_A431_72), rownames(negative_HT29_72))
#1
write.csv(down_72_A431_HT29, "common_down_A431_HT29.csv")
write_xlsx(as.data.frame(down_72_A431_HT29), "common_down_A431_HT29.xlsx")

down_72_A431_all <- intersect(down_72_A431_HT29, down_72_A431_PANC1)
#1!
write.csv(down_72_A431_all, "common_down_A431_all.csv")
write_xlsx(as.data.frame(down_72_A431_all), "common_down_A431_all.xlsx")
all_72_A431 <- union(up_72_A431_all, down_72_A431_all)
#64
length(all_72_A431)
#64


#HT29 vs all
up_72_HT29_PANC1 <- intersect(rownames(positive_HT29_72), rownames(positive_PANC1_WT_vs_72))
#121
write.csv(up_72_HT29_PANC1, "common_up_HT29_PANC1.csv")
write_xlsx(as.data.frame(up_72_HT29_PANC1), "common_up_HT29_PANC1.xlsx")

up_72_HT29_A431 <- intersect(rownames(positive_HT29_72), rownames(positive_A431_72))
#95

up_72_HT29_all <- intersect(up_72_HT29_A431, up_72_HT29_PANC1)
#63
write.csv(up_72_HT29_all, "common_up_HT29_all.csv")
write_xlsx(as.data.frame(up_72_HT29_all), "common_up_HT29_all.xlsx")


down_72_HT29_PANC1 <- intersect(rownames(negative_HT29_72), rownames(negative_PANC1_WT_vs_72))
#3
write.csv(down_72_HT29_PANC1, "common_down_HT29_PANC1.csv")
write_xlsx(as.data.frame(down_72_HT29_PANC1), "common_down_HT29_PANC1.xlsx")

down_72_HT29_A431 <- intersect(rownames(negative_HT29_72), rownames(negative_A431_72))
#1!
down_72_HT29_all <- intersect(down_72_HT29_A431, down_72_A431_PANC1)
#1
write.csv(down_72_HT29_all, "common_down_HT29_all.csv")
write_xlsx(as.data.frame(down_72_HT29_all), "common_down_HT29_all.xlsx")

all_72_HT29 <- union(up_72_HT29_all, down_72_HT29_all)
#64
length(all_72_HT29)
#64


#PANC1 vs all
up_72_PANC1_all <- intersect(up_72_HT29_PANC1, up_72_A431_PANC1)
write.csv(up_72_PANC1_all, "common_up_PANC1_all.csv")
write_xlsx(as.data.frame(up_72_PANC1_all), "common_up_PANC1_all.xlsx")

#63
down_72_PANC1_all <- intersect(down_72_HT29_PANC1, down_72_A431_PANC1)
#1
write.csv(down_72_PANC1_all, "common_down_PANC1_all.csv")
write_xlsx(as.data.frame(down_72_PANC1_all), "common_down_PANC1_all.xlsx")
all_72_PANC1 <- union(up_72_PANC1_all, down_72_PANC1_all)
length(all_72_PANC1)
#64


common_72 <- intersect(all_72_A431, all_72_HT29)
length(common_72)
#64
common_72 <- intersect(common_72, all_72_PANC1)
length(common_72)
#64
intersect(target_list, common_72)

down_72 <- intersect(down_72_A431_all, down_72_HT29_all)
down_72 <- intersect(down_72, down_72_PANC1_all)
#1
write.csv(down_72, "common_down__all.csv")
write_xlsx(as.data.frame(down_72), "common_down__all.xlsx")

up_72 <- intersect(up_72_A431_all, up_72_HT29_all)
up_72 <- intersect(up_72, up_72_PANC1_all)
#63
write.csv(up_72, "common_up__all.csv")
write_xlsx(as.data.frame(up_72), "common_up__all.xlsx")


write.csv(intersect(target_list, common_72), "common_target_fischer_genes.csv")
write.csv((common_72), "common_DE_genes__across_cell_lines.csv")
