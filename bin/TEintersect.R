
# Run first DEA.R script to load in data.

# load files
overlapps <- read.csv("TE_overlapp.txt", sep = "\t", header = FALSE)
up_overlapps <- read.csv("signdiff_up.txt", sep = "\t", header = TRUE)
counts <- read.csv("pp_counts.txt", sep = "\t")
counts_name_type <- counts[,c("gene_name","Type")]

# rename columns and merge files together
colnames(overlapps)[10] <- "gene_name"
colnames(up_overlapps)[1] <- "gene_name"
heatmap_overlapps <- merge(overlapps, up_overlapps, by="gene_name")
barplot_overlapps <- merge(heatmap_overlapps, counts, by="gene_name")
#select columns for heatmap
barplot_overlapps <- barplot_overlapps[, c(7, 23)]
heatmap_overlapps <- heatmap_overlapps[, c(1, 12, 13,14,15,16,17,18,19,20,21,22)]
#make gene_name rownames and delete first column
rownames(heatmap_overlapps) <- heatmap_overlapps$gene_name
heatmap_overlapps_1 <- heatmap_overlapps[,-1]
heatmap_overlapps_2 <- as.vector(heatmap_overlapps_1)
heatmap_overlapps_2 <- order(rowMeans(heatmap_overlapps_2), decreasing = TRUE)

# add small RNA types to data
types_heatmap <- merge(heatmap_overlapps, counts_name_type, by="gene_name")
types_heatmap <- types_heatmap[, c(1,13)]
rownames(types_heatmap) <- types_heatmap$gene_name
types_heat <- types_heatmap[-1]

#plot heatmap of small RNAs overlapping TEs
library(pheatmap)
pheatmap(log2(heatmap_overlapps[heatmap_overlapps_2, as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = F, annotation_col = df, fontsize = 14, annotation_row = types_heat)

# plot barplot of TEs overlapping small RNAs
library(dplyr)
library(ggplot2)
library(hrbrthemes)
barplot_overlapps <- barplot_overlapps %>% group_by(V6, Type) %>% summarize(small = n())
names(barplot_overlapps)[2] <- "type"
type <- barplot_overlapps$type

ggplot(barplot_overlapps, aes(fill=type, y=small, x=V6)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  theme_classic() +
  xlab("") + ylab("frequency") + labs(fill = "small-RNA types") + theme(text = element_text(size = 15))


