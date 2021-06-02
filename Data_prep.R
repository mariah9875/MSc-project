
# load functions
more_0 <- function(row){
  return(length(which(row > 0)))
}

'%ni%' <- Negate('%in%')


# import count files and add a column with type
miRNA <- read.csv("COMPSRA_MERGE_0_miRNA.txt", sep = "\t")
miRNA$Type <- "miRNA"
piRNA <- read.csv("COMPSRA_MERGE_0_piRNA.txt", sep = "\t")
piRNA$Type <- "piRNA"
snoRNA <- read.csv("COMPSRA_MERGE_0_snoRNA.txt", sep = "\t")
snoRNA$Type <- "snoRNA"
snRNA <- read.csv("COMPSRA_MERGE_0_snRNA.txt", sep = "\t")
snRNA$Type <- "snRNA"
tRNA <- read.csv("COMPSRA_MERGE_0_tRNA.txt", sep = "\t")
tRNA$Type <- "tRNA"
circRNA <- read.csv("COMPSRA_MERGE_0_circRNA.txt", sep = "\t")
circRNA$Type <- "circRNA"

library(dplyr)
# rename columns of count files
miRNA <- miRNA %>% rename("gene_name"=miRNA, "DA026"= DA026_SRNA_S7.trimmed_1_sorted_miRNA.txt, "DA027"=DA027_SRNA_S1.trimmed_1_sorted_miRNA.txt, "DA028"=DA028_SRNA_S2.trimmed_1_sorted_miRNA.txt, "DA035"=DA035_SRNA_S9.trimmed_1_sorted_miRNA.txt, "DA036"=DA036_SRNA_S3.trimmed_1_sorted_miRNA.txt, "DA050"=DA050_SRNA_S8.trimmed_1_sorted_miRNA.txt, "DA056"=DA056_SRNA_S5.trimmed_1_sorted_miRNA.txt, "DA057"=DA057_SRNA_S10.trimmed_1_sorted_miRNA.txt, "DA058"=DA058_SRNA_S11.trimmed_1_sorted_miRNA.txt, "DA059"=DA059_SRNA_S6.trimmed_1_sorted_miRNA.txt, "DA061"=DA061_SRNA_S4.trimmed_1_sorted_miRNA.txt)
piRNA <- piRNA %>% rename("gene_name"=piRNA, "DA026"=DA026_SRNA_S7.trimmed_1_sorted_piRNA.txt, "DA027"=DA027_SRNA_S1.trimmed_1_sorted_piRNA.txt, "DA028"=DA028_SRNA_S2.trimmed_1_sorted_piRNA.txt, "DA035"=DA035_SRNA_S9.trimmed_1_sorted_piRNA.txt, "DA036"=DA036_SRNA_S3.trimmed_1_sorted_piRNA.txt, "DA050"=DA050_SRNA_S8.trimmed_1_sorted_piRNA.txt, "DA056"=DA056_SRNA_S5.trimmed_1_sorted_piRNA.txt, "DA057"=DA057_SRNA_S10.trimmed_1_sorted_piRNA.txt, "DA058"=DA058_SRNA_S11.trimmed_1_sorted_piRNA.txt, "DA059"=DA059_SRNA_S6.trimmed_1_sorted_piRNA.txt, "DA061"=DA061_SRNA_S4.trimmed_1_sorted_piRNA.txt)
snoRNA <- snoRNA %>% rename("gene_name"=snoRNA, "DA026"=DA026_SRNA_S7.trimmed_1_sorted_snoRNA.txt, "DA027"=DA027_SRNA_S1.trimmed_1_sorted_snoRNA.txt, "DA028"=DA028_SRNA_S2.trimmed_1_sorted_snoRNA.txt, "DA035"=DA035_SRNA_S9.trimmed_1_sorted_snoRNA.txt, "DA036"=DA036_SRNA_S3.trimmed_1_sorted_snoRNA.txt, "DA050"=DA050_SRNA_S8.trimmed_1_sorted_snoRNA.txt, "DA056"=DA056_SRNA_S5.trimmed_1_sorted_snoRNA.txt, "DA057"=DA057_SRNA_S10.trimmed_1_sorted_snoRNA.txt, "DA058"=DA058_SRNA_S11.trimmed_1_sorted_snoRNA.txt, "DA059"=DA059_SRNA_S6.trimmed_1_sorted_snoRNA.txt, "DA061"=DA061_SRNA_S4.trimmed_1_sorted_snoRNA.txt)
snRNA <- snRNA %>% rename("gene_name"=snRNA, "DA026"=DA026_SRNA_S7.trimmed_1_sorted_snRNA.txt, "DA027"=DA027_SRNA_S1.trimmed_1_sorted_snRNA.txt, "DA028"=DA028_SRNA_S2.trimmed_1_sorted_snRNA.txt, "DA035"=DA035_SRNA_S9.trimmed_1_sorted_snRNA.txt, "DA036"=DA036_SRNA_S3.trimmed_1_sorted_snRNA.txt, "DA050"=DA050_SRNA_S8.trimmed_1_sorted_snRNA.txt, "DA056"=DA056_SRNA_S5.trimmed_1_sorted_snRNA.txt, "DA057"=DA057_SRNA_S10.trimmed_1_sorted_snRNA.txt, "DA058"=DA058_SRNA_S11.trimmed_1_sorted_snRNA.txt, "DA059"=DA059_SRNA_S6.trimmed_1_sorted_snRNA.txt, "DA061"=DA061_SRNA_S4.trimmed_1_sorted_snRNA.txt)
tRNA <- tRNA %>% rename("gene_name"=tRNA, "DA026"=DA026_SRNA_S7.trimmed_1_sorted_tRNA.txt, "DA027"=DA027_SRNA_S1.trimmed_1_sorted_tRNA.txt, "DA028"=DA028_SRNA_S2.trimmed_1_sorted_tRNA.txt, "DA035"=DA035_SRNA_S9.trimmed_1_sorted_tRNA.txt, "DA036"=DA036_SRNA_S3.trimmed_1_sorted_tRNA.txt, "DA050"=DA050_SRNA_S8.trimmed_1_sorted_tRNA.txt, "DA056"=DA056_SRNA_S5.trimmed_1_sorted_tRNA.txt, "DA057"=DA057_SRNA_S10.trimmed_1_sorted_tRNA.txt, "DA058"=DA058_SRNA_S11.trimmed_1_sorted_tRNA.txt, "DA059"=DA059_SRNA_S6.trimmed_1_sorted_tRNA.txt, "DA061"=DA061_SRNA_S4.trimmed_1_sorted_tRNA.txt)
circRNA <- circRNA %>% rename("gene_name"=circRNA, "DA026"=DA026_SRNA_S7.trimmed_1_sorted_circRNA.txt, "DA027"=DA027_SRNA_S1.trimmed_1_sorted_circRNA.txt, "DA028"=DA028_SRNA_S2.trimmed_1_sorted_circRNA.txt, "DA035"=DA035_SRNA_S9.trimmed_1_sorted_circRNA.txt, "DA036"=DA036_SRNA_S3.trimmed_1_sorted_circRNA.txt, "DA050"=DA050_SRNA_S8.trimmed_1_sorted_circRNA.txt, "DA056"=DA056_SRNA_S5.trimmed_1_sorted_circRNA.txt, "DA057"=DA057_SRNA_S10.trimmed_1_sorted_circRNA.txt, "DA058"=DA058_SRNA_S11.trimmed_1_sorted_circRNA.txt, "DA059"=DA059_SRNA_S6.trimmed_1_sorted_circRNA.txt, "DA061"=DA061_SRNA_S4.trimmed_1_sorted_circRNA.txt)

# stack data together into 1 countmatrix
merge_all <- rbind(miRNA, piRNA, snoRNA, snRNA, tRNA, circRNA)
with_id_type <- merge_all[, c(1,14, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]

# delete rows containing 0 counts
id_type <- with_id_type[which(apply(with_id_type[,-1][,-1], 1, FUN = more_0) > 3),]


# add genomic locations of small RNAs to countmatrix
# import annotation files
miRNA_db <- as.data.frame(rtracklayer::import('hsa_miRNAbase.gff3'))
circRNA_db <- read.csv("circRNA_liftoverhg38.txt", sep = "\t", header = FALSE)
circRNA_19 <- read.csv("circRNAs_hg19.txt", sep = "\t", header = FALSE)
piRNA_db <- read.csv("hsa_piRbase.txt", sep= "\t", header = FALSE)
piRNA_bank <- read.csv("piRNAbank_hg38.txt", sep= "\t", header = FALSE)
gencode <- as.data.frame(rtracklayer::import('gencode.v37.annotation.gtf'))
tRN_db <- read.csv("tRNA_liftoverhg38.txt", sep= "\t", header = FALSE)

# reformat gencode file and subset for snoRNA and snRNA
gtf_genes <- subset(gencode, gencode$type == "gene")
snoRNA <- subset(gtf_genes, gtf_genes$gene_type == "snoRNA")
snRNA <- subset(gtf_genes, gtf_genes$gene_type == "snRNA")


# rename columns of annotation files
miRNA_db <- miRNA_db %>% rename("gene_name"=Name, "chr"=seqnames, "size"=width)
circRNA_db <- circRNA_db %>% rename("chr"=V1, "start"=V2, "end"=V3, "strand"=V4, "id"=V5, "gene_name"=V6)
circRNA_db$size <- circRNA_db$end - circRNA_db$start
piRNA_db <- piRNA_db %>% rename("chr"=V1, "start"=V2, "end"=V3, "strand"=V6, "gene_name"=V4)
piRNA_db$size <- piRNA_db$end - piRNA_db$start
piRNA_bank <- piRNA_bank %>% rename("chr"=V1, "start"=V2, "end"=V3, "strand"=V4, "gene_name"=V5, "id"=V6)
piRNA_bank$size <- piRNA_bank$end - piRNA_bank$start
tRN_db <- tRN_db %>% rename("chr"=V1, "start"=V2, "end"=V3, "strand"=V4, "gene_name"=V5)
tRN_db$size <- tRN_db$end - tRN_db$start
circRNA_19 <- circRNA_19 %>% rename("id"=V1, "best_transcript"=V6, "gene_name"=V7)
snoRNA <- snoRNA %>% rename("chr"=seqnames, "size"=width)
snRNA <- snRNA %>% rename("chr"=seqnames, "size"=width)


# extract columns from circRNA annotation files
circRNA_data <- select(circRNA_db, chr, start, end, size, id, gene_name, strand)
circRNA_19data <- select(circRNA_19, id, best_transcript, gene_name)
# extract the transcript_id for circRNAs
cirncDC <- merge(circRNA_data, circRNA_19data, by="id")
cirncDC <- cirncDC %>% rename("gene_name"=best_transcript)


# select columns to look at small RNA location
miRNA_data <- select(miRNA_db, chr, size, gene_name)
circRNA_data <- cirncDC
piRNA_data <- select(piRNA_db, size, gene_name)
snoRNA_data <- select(snoRNA, chr, start, end, size, gene_name, gene_id, strand)
snRNA_data <- select(snRNA, chr, start, end, size, gene_name, gene_id, strand)
tRNA_data <- select(tRN_db, chr, start, end, size, gene_name, strand)
piRNA_bank_one <- select(piRNA_bank, size, gene_name)

# check if duplicates and delete
piRNAbank_one <- piRNA_bank_one %>% distinct()
miRNA_one <- miRNA_data %>% distinct()
piRNA_one <- piRNA_data %>% distinct()
snoRNA_one <- snoRNA_data %>% distinct()
snRNA_one <- snRNA_data %>% distinct()
tRNA_one <- tRNA_data %>% distinct()
circRNA_one <- circRNA_data %>% distinct()

# concatinate all annotation files with countmatrix and extract columns after deleting duplicates
miRNA_location <- merge(id_type,miRNA_one, by="gene_name")
miRNA_location <- select(miRNA_location, chr, gene_name, size)
piRNA_location <- merge(id_type,piRNA_one, by="gene_name")
piRNA_location <- select(piRNA_location, gene_name, size)
piRNAbank_location <- merge(id_type,piRNAbank_one, by="gene_name")
piRNAbank_location <- select(piRNAbank_location, gene_name,size)
snoRNA_location <- merge(id_type,snoRNA_one, by="gene_name")
snoRNA_location <- select(snoRNA_location, chr, start, end, gene_name, gene_id, size)
snRNA_location <- merge(id_type,snRNA_one, by="gene_name")
snRNA_location <- select(snRNA_location, chr, start, end, gene_name, gene_id, size)
tRNA_location <- merge(id_type,tRNA_one, by="gene_name")
tRNA_location <- select(tRNA_location, chr, start, end, gene_name, size)
circRNA_location <- merge(id_type,circRNA_one, by="gene_name")


# remove circRNAs, snoRNA, snRNA larger than 200 from countfile
circRNA_truesize <- circRNA_location[circRNA_location$size < 201, ]
snoRNA_truesize <- snoRNA_location[snoRNA_location$size < 201, ]
snRNA_truesize <- snRNA_location[snRNA_location$size > 201, ]

# countmatrix
id_type_true <- id_type

# remove >200 nt for circRNAs and merge with countmatrix
circO <- subset(id_type_true, Type == "circRNA")
circ <- circO[which(circO$gene_name %in% circRNA_truesize$gene_name), "gene_name"]
circ <- as.data.frame(circ)
colnames(circ)[1] <- "gene_name"
circ <- merge(circ, circO, by="gene_name")
circ_size <- merge(circ, circRNA_truesize, by="gene_name")

# remove >200 nt for snoRNAs and merge with countmatrix
snoR <- subset(id_type_true, Type == "snoRNA")
sno <- snoR[which(snoR$gene_name %in% snoRNA_truesize$gene_name), "gene_name"]
sno <- as.data.frame(sno)
colnames(sno)[1] <- "gene_name"
sno <- merge(sno, snoR, by="gene_name")
s <- snoR[which(snoR$gene_name %ni% snoRNA_location$gene_name), "gene_name"]
s <- as.data.frame(s)
colnames(s)[1] <- "gene_name"
s <- merge(s, snoR, by="gene_name")
snO <- rbind(sno, s)
sno_size <- merge(snO, snoRNA_truesize, by="gene_name")

# remove >200 nt for snRNAs and merge with countmatrix
sn0 <- subset(id_type_true, Type == "snRNA")
sn <- sn0[which(sn0$gene_name %ni% snRNA_truesize$gene_name), "gene_name"]
sn <- as.data.frame(sn)
colnames(sn)[1] <- "gene_name"
sn <- merge(sn, sn0, by="gene_name")
sn_size <- merge(sn, snRNA_truesize, by="gene_name")

# merge countmatrix with size for miRNAs, tRNAs and piRNAs
t <- subset(id_type_true, Type == "tRNA")
t_size <- merge(t, tRNA_location, by="gene_name")
mi <- subset(id_type_true, Type == "miRNA")
mi_size <- merge(mi, miRNA_location, by="gene_name")
pi <- subset(id_type_true, Type == "piRNA")
pi_size <- merge(pi, piRNA_location, by="gene_name")
pi_z <- merge(pi, piRNAbank_location, by="gene_name")
p_size <- rbind(pi_size, pi_z)

id_type_true <- rbind(circ, snO, sn, t, mi, pi)

library(ggplot2)
# size distribution plots to verify that the small RNAs are <200 nt.
# merge size of all small RNA types in 1 file
circ_size <- select(circ_size, gene_name, size)
circ_size <- circ_size %>% distinct()
p_size <- select(p_size, gene_name, size)
mi_size <- select(mi_size, gene_name, size)
sno_size <- select(sno_size, gene_name, size)
t_size <- select(t_size, gene_name, size)
sn_size <- select(sn_size, gene_name, size)
all_size <- rbind(circ_size, p_size, mi_size, sno_size, t_size, sn_size)

# plot size distribution
size <- select(all_size, size)

size %>%
  filter( size<300 ) %>%
  ggplot( aes(x=size)) +
  stat_bin(breaks=seq(0,200,2), fill="#69b3a2", color="#e9ecef", alpha=0.9)  + theme_linedraw() + labs(x="length", y="frequency")


# plot piRNA and miRNA size 
plot_piRNA <- (table(p_size$size))
plot_piRNA <- as.data.frame(plot_piRNA)

plot_miRNA <- (table(mi_size$size))
plot_miRNA <- as.data.frame(plot_miRNA)

ggplot(p_size, aes(x="piRNAs", y=size, fill="transcripts")) + geom_violin() + ylim(c(10, 35)) + theme_classic() + labs(y="length", x="") + theme(legend.position =  "none") + theme(text = element_text(size = 20))

ggplot(mi_size, aes(x="miRNAs", y=size, fill="transcripts")) + geom_violin() + ylim(c(0, 200)) + theme_classic() + labs(y="length", x="") + theme(legend.position =  "none") + theme(text = element_text(size = 20))

# Remove small RNAs positioned on chr M
# concatenate all the annotation files to check for chr names
miRNA_chrM <- select(miRNA_db, chr, gene_name)
piRNA_chrM <- select(piRNA_db, chr, gene_name)
snoRNA_chrM <- select(snoRNA, chr, gene_name)
snRNA_chrM <- select(snRNA, chr, gene_name)
tRNA_chrM <- select(tRN_db, chr, gene_name)
piRNA_bank_chrM <- select(piRNA_bank, chr, gene_name)
circRNA_chrM <- select(circRNA_data, chr, gene_name)

all <- rbind(miRNA_chrM, piRNA_chrM, snoRNA_chrM, snRNA_chrM, tRNA_chrM, piRNA_bank_chrM, circRNA_chrM)

# delete small RNAs positioned on chrM
small_RNAs <- id_type_true

small <- merge(small_RNAs, all, by="gene_name")
small <- select(small, gene_name, chr)
small <- small %>% distinct()

# subset all small RNAs located on chr M
small_subset <- subset(small, chr == "chrM")

# manually check if these small RNAs are located on other chr
# if located on other chr than M then delete them from subset
small_subset <- small_subset[-c(2, 3, 5, 6, 23, 24, 29), ]

# if small RNAs not in chr M subset then print to new file
small_chrM <- small_RNAs[which(small_RNAs$gene_name %ni% small_subset$gene_name), "gene_name"]
small_chrM <- as.data.frame(small_chrM)
colnames(small_chrM)[1] <- "gene_name"

# new countmatrix ready for import in the DEA script, all small RNAs > 200 nt removed.
small_no_chrM <- merge(small_chrM, small_RNAs, by="gene_name")
small_no_chrM <- small_no_chrM[-c(1363),]

write.table(small_no_chrM, file="pp_counts.txt", sep = "\t", quote = FALSE)
