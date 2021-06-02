
# Load functions
more_3 <- function(row){
  return(length(which(row >= 3)))
} 

more_12 <- function(row){
  return(length(which(row >= 12)))
}

less_4 <- function(row){
  return(length(which(row < 4)))
}

library(dplyr)
library(ggplot2)

# import processed count file from Data.prep.R
counts <- read.csv("pp_counts.txt", sep = "\t")
counts_name_type <- counts[,c("gene_name","Type")]

# select low count threshold
mean <- rowMeans(counts[,-1][,-1])
mean <- as.data.frame(mean)
log_mean <- log2(mean + 0.5)

ggplot(log_mean, aes(y=mean)) + geom_boxplot() + theme_linedraw() + ylab("log2 mean of counts") + theme(text = element_text(size = 20))

# value 3 selected as low value
2^1.5

p <- mean %>%
  filter(mean < 100) %>%
  ggplot(aes(x=mean)) + geom_histogram(binwidth = 1) + theme_classic() + xlab( "mean counts") + ylab("frequency") + theme(text = element_text(size = 20))

p


# Plot distribution of human and chimpanzee types
# subset human and chimp groups
human <- counts[,c("Type", "DA028", "DA035", "DA036", "DA050", "DA056", "DA057")]
chimp <- counts[,c("Type", "DA026", "DA027", "DA058", "DA059", "DA061")]

# filter out low counts from human and chimp groups
human <- human[which(apply(human[,-1], 1, FUN = more_3) > 3),]
chimp <- chimp[which(apply(chimp[,-1], 1, FUN = more_3) > 3),]

all_Type <- counts

# count types for human columns
human_type <- (table(human$Type))
human_type <- as.data.frame(human_type)

# count types for chimp columns
chimp_type <- (table(chimp$Type))
chimp_type <- as.data.frame(chimp_type)

# count types for both
all_type <- (table(all_Type$Type))
all_type <- as.data.frame(all_type)

chimp <- chimp_type
human <- human_type
type <- all_type

chimp$Var1 <- as.character(chimp$Var1)
human$Var1 <- as.character(human$Var1)
type$Var1 <- as.character(type$Var1)

zhuman <- as.data.frame(human)
h <- sum(zhuman$Freq)
zhuman$percentage <- (zhuman$Freq/h)*100
zchimp <- as.data.frame(chimp)
c <- sum(zchimp$Freq)
zchimp$percentage <- (zchimp$Freq/c)*100
ztype <- as.data.frame(type)

# plot distribution of small RNA types
# human group
ggplot(zhuman, aes(x=Var1, y=percentage, fill=Var1, label=Freq)) + geom_bar(stat="identity") + scale_fill_manual(values = c("red", "blue", "orange", "green", "purple", "grey")) + xlab("") + ylab("percentage %") + geom_text(size = 6, position = position_stack(vjust = 0.5)) +
  guides(fill=FALSE, color=FALSE) + theme_classic() + ylim(0,60) + theme(text = element_text(size = 20))
# chimpanzee group
ggplot(zchimp, aes(x=Var1, y=percentage, fill=Var1, label=Freq)) + geom_bar(stat="identity") + scale_fill_manual(values = c("red", "blue", "orange", "green", "purple", "grey")) + xlab("") + ylab("percentage %") + geom_text(size = 6, position = position_stack(vjust = 0.5)) +
  guides(fill=FALSE, color=FALSE) + theme_classic() + ylim(0,60) + theme(text = element_text(size = 20))


# remove low counts from countfile before input into DESeq2
countdata <- counts
countdata <- countdata[which(apply(countdata[,-1][,-1], 1, FUN = more_3) > 3),]
# make ids rownames
rownames(countdata) <- countdata$gene_name
#rownames(countdata) <- countdata$V1
# delete first 2 columns
countdata <- countdata[-1]
countdata <- countdata[-1]

library(DESeq2)
# DESEQ2
# Metadata
(condition <- factor(c("Chimp", "Chimp", "Human", "Human", "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp")))
(coldata <- data.frame(row.names=colnames(countdata),sample=colnames(countdata), condition))

# DESeq object
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)

# Set factor and run DESeq
dds$condition <- factor(dds$condition, levels = c('Chimp', 'Human'))
dds <- DESeq(dds)

# DESeq results object
res <- results(dds)
summary(res)
res

# Regularized log transformation for visualization
# Normalization
norm <- counts(dds, normalized=TRUE)
norm_dds <- norm
rld <- rlogTransformation(dds)

# Adjusted p-value, up and down-regulated genes
table(res$padj < 0.05)
table(res$padj < 0.05 & res$log2FoldChange > 1)
table(res$padj < 0.05 & res$log2FoldChange < -1)
table(res$padj < 1)

# Print to table expressed, up and downregulated genes
up_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
down_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
down_regulated <- as.data.frame(down_regulated)
all_genes <- subset(res, res$padj < 0.05)
all_genes_significant <- subset(all_genes, all_genes$log2FoldChange >= -1 & all_genes$log2FoldChange <= 1)
All <- subset(res)

write.table(up_regulated, file="up_regulated.txt", sep = "\t", quote = FALSE)
write.table(down_regulated, file= "down_regulated.txt", sep = "\t", quote = FALSE)
write.table(all_genes, file="all_regulated.txt", sep = "\t", quote = FALSE)


library(data.table)
## Visualization
# Principal Component Analysis
plotPCA(rld) + theme_classic()
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd) + theme_classic()

# MA plot
plotMA(res, alpha=0.05, ylim=c(-15,15))

# Mean plot
norm_dds <- as.data.frame(norm_dds)
norm_dds_test <- norm_dds
# Mean of human and chimp
human_mean <- rowMeans(norm_dds_test[,subset(coldata, coldata$condition == "Human")$sample])
h <- as.data.frame(human_mean)
chimp_mean <- rowMeans(norm_dds_test[,subset(coldata, coldata$condition == "Chimp")$sample])
c <- as.data.frame(chimp_mean)
plot_gg <- cbind(c,h)

# Genes found in humans and chimp
human_genes_expressed <- norm_dds_test[,subset(coldata, coldata$condition == "Human")$sample]
human_genes_expressed <- human_genes_expressed[rowSums(human_genes_expressed)>0,]
chimp_genes_expressed <- norm_dds_test[,subset(coldata, coldata$condition == "Chimp")$sample]
chimp_genes_expressed <- chimp_genes_expressed[rowSums(chimp_genes_expressed)>0,]

# Mark up-regulate genes red and down-regulated genes blue
norm_dds_test$cex <- rep("300", nrow(norm_dds_test))
norm_dds_test[rownames(up_regulated),] <- "1"
norm_dds_test[rownames(down_regulated),] <- "4000"

ggplot(plot_gg, aes(x=log2(chimp_mean + 0.5), y=log2(human_mean+0.5), color = norm_dds_test$cex)) + 
  (geom_point(size=0.5)) + scale_colour_manual(values = c("#FC4E07","#999999","#0072B2"), name = "", labels = c("up-regulated (423)", "not significant (6830)", "down-regulated (278)")) +
  theme_classic() + ylab("log2(mean Human)") + xlab("log2(mean Chimpanzee") + theme(legend.text = element_text(size = 12))


# Plot up and down-regulated types after running DEA
# merge upregulated with type 
up <- read.csv("up_regulated.txt", sep = "\t")
setDT(up, keep.rownames="gene_name")[] 
up_1 <- merge(up, counts, by = "gene_name")
up_1 <- up_1[, c(1,8,2,3,4,5,6,7)]
up_1 <- as.data.frame(up_1)
write.table(up_1, file="up-regulated_type.txt", sep = "\t", quote = FALSE)

# merge down regulated with type
down <- read.csv("down_regulated.txt", sep = "\t")
setDT(down, keep.rownames="gene_name")[] 
down_1 <- merge(down, counts, by = "gene_name")
down_1 <- down_1[, c(1,8,2,3,4,5,6,7)]
down_1 <- as.data.frame(down_1)
write.table(down_1, file="down-regulated_type.txt", sep = "\t", quote = FALSE)

# merge regulated with type
all_1 <- read.csv("all_regulated.txt", sep= "\t")
setDT(all_1, keep.rownames="gene_name")[] 
all_1 <- merge(all_1, counts, by = "gene_name")
all_1 <- all_1[, c(1,8,2,3,4,5,6,7)]

library(multcomp)
# bar-plot of smallRNAs and normalized counts from all samples
signdiff_up <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange > 1)),]
sum_columns <- as.data.frame(rowSums(signdiff_up))
names(sum_columns)[1] <- "sum_counts"
setDT(sum_columns, keep.rownames = "gene_name")[]
normalized_counts_up <- as.data.frame(signdiff_up)
setDT(normalized_counts_up, keep.rownames="gene_name")[] 
normalized_counts_up <- merge(normalized_counts_up, sum_columns, by="gene_name")
normalized_counts_up<- merge(normalized_counts_up, up_1, by = "gene_name")

normal_counts <- normalized_counts_up[, c(14,15)]
SE <- function(x) (sum(x))
se.norm <- aggregate(list(SE=normal_counts$baseMean), by=list(type=normal_counts$Type), FUN=SE)

# plot up-regulated 
up_type <- (table(normal_counts$Type))
up_type <- as.data.frame(up_type)
up_perc <- sum(up_type$Freq)
up_type$percentage <- (up_type$Freq/up_perc)*100
up_type <- up_type[order(up_type$Freq, decreasing = TRUE),]

ggplot(up_type, aes(x=Var1, y=percentage, fill=Var1, label=Freq)) + geom_bar(stat="identity") + scale_fill_manual(values = c("red", "blue", "orange", "green", "purple", "grey")) + xlab("") + ylab("percentage %") + geom_text(size = 6, position = position_stack(vjust = 0.5)) +
  guides(fill=FALSE, color=FALSE) + theme_classic() + ylim(0,60) + theme(text = element_text(size = 20))

# down regulated
signdiff_low <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange < -1)),]
sum_columns <- as.data.frame(rowSums(signdiff_low))
names(sum_columns)[1] <- "sum_counts"
setDT(sum_columns, keep.rownames = "gene_name")[]
normalized_counts_down <- as.data.frame(signdiff_low)
setDT(normalized_counts_down, keep.rownames="gene_name")[] 
normalized_counts_down <- merge(normalized_counts_down, sum_columns, by="gene_name")
normalized_counts_down <- merge(normalized_counts_down, down_1, by = "gene_name")

normal_counts <- normalized_counts_down[, c(14,15)]
SE <- function(x) (sum(x))
se.norm <- aggregate(list(SE=normal_counts$baseMean), by=list(type=normal_counts$Type), FUN=SE)

# plot down-regulated
down_type <- (table(normal_counts$Type))
down_type <- as.data.frame(down_type)
down_perc <- sum(down_type$Freq)
down_type$percentage <- (down_type$Freq/down_perc)*100
down_type <- down_type[order(down_type$Freq, decreasing = TRUE),]

ggplot(down_type, aes(x=Var1, y=percentage, fill=Var1, label=Freq)) + geom_bar(stat="identity") + scale_fill_manual(values = c("red", "blue", "orange", "green", "purple", "grey")) + xlab("") + ylab("percentage %") + geom_text(size = 6, position = position_stack(vjust = 0.5)) +
  guides(fill=FALSE, color=FALSE) + theme_classic() + ylim(0,40) + theme(text = element_text(size = 20))


# Plot different heatmaps to explore the DEA
# Heatmaps
library(gplots)           
library(RColorBrewer)
library(pheatmap)

# select condition to be compared in heatmaps
coldata <- coldata[order(coldata$condition),]
df <- as.data.frame(colData(rld)[,c("condition"),drop=FALSE])

# Heatmap for all regulated small RNAs
nt <- norm[rownames(subset(res, res$padj < 0.05)),]
nt <- as.data.frame(nt)
ntb <- order(rowMeans(nt), decreasing = TRUE)

counts_scaled <- log2(nt + 0.5)
counts_scaled <- as.data.frame(counts_scaled)
mycol <- colorpanel(100, "yellow", "black", "red")

plot_heatmaps <- function(x, x_dataFrame) {
  plot <- pheatmap(as.matrix(counts_scaled[rownames(x),rownames(x_dataFrame)]),
                   show_colnames = F, annotation_col=x_dataFrame,
                   labels_col = x_dataFrame[,1], annotation_names_row = F,
                   show_rownames = F, cluster_cols = F,
                   clustering_distance_rows = "correlation", clustering_method = "ward.D2",
                   annotation_names_col=F, col=mycol, scale="row", breaks = breaks, fontsize = 15)
  return(plot)
}

breaks <- seq(-2,2, length.out = 101)
breaks[length(breaks)] <- max(max(counts_scaled), max(breaks))
breaks[1] <- min(-max(counts_scaled), min(breaks))

# Heatmap for all the small RNAs regulated
human_chimp_dataframe <- as.data.frame(as.character(coldata$condition), row.names = coldata$sample)
colnames(human_chimp_dataframe) <- "condition"
plot_heatmaps(nt, human_chimp_dataframe)


# Heatmap up regulated small RNAs
signdiff_up <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange > 1)),]
signdiff_overlaps <- as.data.frame(signdiff_up)
write.table(signdiff_overlaps, file="signdiff_up.txt", sep = "\t", quote = FALSE, col.names = NA)

upp <- as.data.frame(signdiff_up)
uppp <- upp[,c("DA028", "DA035", "DA036", "DA050", "DA056", "DA057")]
signdiff_uppp <- order(rowMeans(uppp), decreasing = TRUE) [1:40]

pheatmap(log2(signdiff_up[signdiff_uppp, as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = F, annotation_col = df, fontsize = 14)

# Heatmap of down regulated small RNAs
signdiff_low <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange < -1)),]

down <- as.data.frame(signdiff_low)
down_1 <- down[,c("DA026", "DA027", "DA058", "DA059", "DA061")]
signdiff_lowww <- order(rowMeans(down_1), decreasing = TRUE) [1:40]

pheatmap(log2(signdiff_low[signdiff_lowww, as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = F, annotation_col = df, fontsize = 14)


# Heatmap for piRNAs & miRNAs
# extract piRNAs and miRNAs
expressed <- norm[rownames(subset(res, res$padj < 0.05)),]
expressed <- setDT(as.data.frame(expressed), keep.rownames = "gene_name")[]
expressed <- merge(expressed, counts_name_type, by="gene_name")
expressed <- as.data.frame(expressed)
rownames(expressed) <- expressed$gene_name
expressed <- expressed[-1]

piRNA <-  subset(expressed, Type == "piRNA")
piRNA <- as.data.frame(piRNA[,c("DA026", "DA027","DA028", "DA035", "DA036", "DA050", "DA056", "DA057","DA058", "DA059", "DA061")])
piRNA_select <- order(rowMeans(piRNA), decreasing = TRUE)

miRNA <- subset(expressed, Type == "miRNA")
miRNA <- miRNA[,c("DA026", "DA027","DA028", "DA035", "DA036", "DA050", "DA056", "DA057","DA058", "DA059", "DA061")]
miRNA_select <- order(rowMeans(miRNA), decreasing = TRUE)

# plot heatmap piRNA
counts_scaled <- log2(piRNA + 0.5)
plot_heatmaps(piRNA, human_chimp_dataframe)

# plot heatmap miRNA
counts_scaled <- log2(miRNA+0.5)
plot_heatmaps(miRNA, human_chimp_dataframe)



# DEA for high counts HUMAN vs low counts CHIMP

# extract human and chimp columns with Type
human <- counts[,c("Type", "DA028", "DA035", "DA036", "DA050", "DA056", "DA057")]
chimp <- counts[,c("Type", "DA026", "DA027", "DA058", "DA059", "DA061")]
# subset human and chimp with gene name 
human_id <- counts[,c("gene_name", "DA028", "DA035", "DA036", "DA050", "DA056", "DA057")]
chimp_id <- counts[,c("gene_name", "DA026", "DA027", "DA058", "DA059", "DA061")]

# find out threshold high counts for human and low counts for chimp
# human
mean <- rowMeans(human[,-1])
mean <- as.data.frame(mean)
log_mean <- log2(mean + 0.5)
log2_human <- log_mean
log2_human$condition <- rep("human", nrow(log2_human))

ggplot(log_mean, aes(y=mean)) + geom_boxplot() + theme_linedraw() + theme(text = element_text(size = 20))

# threshold for high counts
2^3.5

p <- log_mean %>%
  ggplot(aes(x=mean)) + 
  geom_histogram(binwidth = 1) + theme_classic() + xlab("log2 mean counts") + ylab("frequency")

p
p <- mean %>%
  filter(mean < 100) %>%
  ggplot(aes(x=mean)) + geom_histogram(binwidth = 1) + theme_classic() + xlab( "mean counts") + ylab("frequency") + theme(text = element_text(size = 20))

p

# chimp
mean <- rowMeans(chimp[,-1])
mean <- as.data.frame(mean)
log_mean <- log2(mean + 0.5)
log2_chimp <- log_mean
log2_chimp$condition <- rep("chimpanzee", nrow(log2_chimp))

ggplot(log_mean, aes(y=mean)) + geom_boxplot() + theme_linedraw() + theme(text = element_text(size = 20))

# threshold for low counts
2^1.5

p <- log_mean %>%
  ggplot(aes(x=mean)) + 
  geom_histogram(binwidth = 1) + theme_classic() + xlab("log2 mean counts") + ylab("frequency")

p
p <- mean %>%
  filter(mean < 100) %>%
  ggplot(aes(x=mean)) + geom_histogram(binwidth = 1) + theme_classic() + xlab( "mean counts") + ylab("frequency") + theme(text = element_text(size = 20))

p

# boxplot with both human and chimp
log2_ch <- rbind(log2_human, log2_chimp)
ggplot(log2_ch, aes(x=condition, y=mean, fill=condition)) + geom_boxplot() + theme_linedraw() + ylab("log2(mean)")


# Before running the DEA again, subset human group with high counts and subset chimp with low counts
# filter out high counts from chimp and low counts from human
human_id <- human_id[which(apply(human_id[,-1], 1, FUN = more_12) > 3),]
nrow(human_id)
chimp_id <- chimp_id[which(apply(chimp_id[,-1], 1, FUN = less_4) > 3),]
nrow(chimp_id)
human_chimp_significant <- merge(human_id, chimp_id, by="gene_name")
nrow(human_chimp_significant)

human_chimp_significant <- human_chimp_significant[,c("gene_name","DA026", "DA027", "DA028", "DA035", "DA036", "DA050", "DA056", "DA057", "DA058", "DA059", "DA061")]


# DESEQ2 
countdata <- human_chimp_significant
# make ids rownames
rownames(countdata) <- countdata$gene_name
# delete first 2 columns
countdata <- countdata[-1]

# Metadata
(condition <- factor(c("Chimp", "Chimp", "Human", "Human", "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp")))
(coldata <- data.frame(row.names=colnames(countdata),sample=colnames(countdata), condition))

# DESeq object
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)

# Set factor and run DESeq
dds$condition <- factor(dds$condition, levels = c('Chimp', 'Human'))
dds <- DESeq(dds)

# DESeq results object
res <- results(dds)
summary(res)
res

# Regularized log transformation for visualization
# Normalization
norm <- counts(dds, normalized=TRUE)
norm_dds <- norm

# Adjusted p-value, up and down-regulated genes
table(res$padj < 0.05)
table(res$padj < 0.05 & res$log2FoldChange > 1)
table(res$padj < 0.05 & res$log2FoldChange < -1)

# Print to table expressed, up and downregulated genes
up_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
down_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
all_genes <- subset(res, res$padj < 0.05)

write.table(up_regulated, file="up_regulated_humanhigh.txt", sep = "\t", quote = FALSE)
write.table(down_regulated, file= "down_regulated_humanhigh.txt", sep = "\t", quote = FALSE)
write.table(all_genes, file="all_regulated_humanhigh.txt", sep = "\t", quote = FALSE)

# Heatmaps
# select condition to be compared in heatmaps
coldata <- coldata[order(coldata$condition),]
# up and down regulated small RNAs
signdiff_up <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange > 1)),]
signdiff_upppp <- order(rowMeans(signdiff_up), decreasing = TRUE)

# Heatmap showing highly expressed small-RNAs in humans vs chimpanzee
pheatmap(log2(signdiff_up[signdiff_upppp, as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = F, annotation_col = df, fontsize = 14)

# barplot showing up-regulated types
up_plot <- as.data.frame(up_regulated)
setDT(up_plot, keep.rownames = "gene_name")[]
up_plot <- merge(up_plot, counts, by="gene_name")

up_human <- (table(up_plot$Type))
up_human <- as.data.frame(up_human)
up_perc <- sum(up_human$Freq)
up_human$percentage <- (up_human$Freq/up_perc)*100
up_human <- up_human[order(up_human$Freq, decreasing = TRUE),]

ggplot(up_human, aes(x=Var1, y=percentage, fill=Var1, label=Freq)) + geom_bar(stat="identity") + scale_fill_manual(values = c("red", "blue", "orange", "green", "grey")) + xlab("") + ylab("percentage %") + geom_text(size = 6, position = position_stack(vjust = 0.5)) +
  guides(fill=FALSE, color=FALSE) + theme_classic() + ylim(0,80) + theme(text = element_text(size = 20))

