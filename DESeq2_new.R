
# This workflow excludes sg1

library(DESeq2)
library(dplyr)

# Set working directory
setwd("/Users/antioxidant233/Documents/Zhen_Lab/Workflow")

# Read count matrix
cts <- read.csv("geneCounts_wc-rev_s-un_exc-sg1.csv")
# View(cts)
#cg1-cg5, sg2-sg6

# Set first column (genes) as row names
dds2 <- cts[, -1]
rownames(dds2) <- cts$X
cts <- dds2


# Read condition file
coldata <- read.csv("coldata.csv", row.names = 1)
# View(coldata)

# Factor each column
coldata$conditions <- factor(coldata$conditions)
coldata$genotype <- factor(coldata$genotype)
coldata$combined <- factor(coldata$combined)

# Make dds
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ combined)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)


# Factor each column for dds
dds$conditions <- factor(dds$conditions, levels = c("sorted","whole"))
dds$genotype <- factor(dds$genotype, levels = c("gain", "wild", "loss"))
dds$combined <- relevel(dds$combined, ref = "whole_wild")


# Run differential expression analysis
dds <- DESeq(dds)


# Plot counts for the gene with minimum P-adj
plotCounts(dds, gene=which.min(res$padj), intgroup="combined")
# ?plotCounts


# View results for significant genes
# P-adj threshold: 0.05
# Compared: sorted_genotype vs whole_genotype

# ###### Compare sorted gf v.s. whole gf ########
# res <- results(dds, contrast=c("combined","sorted_gain","whole_gain"))
# 
# resultsNames(dds)
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj < 0.05, na.rm=TRUE)
# 
# resSig <- subset(resOrdered, padj < 0.05)
# 
# up <- subset(resSig, log2FoldChange > 2) # Up-regulated genes
# down <- subset(resSig, log2FoldChange < -2) # Down-regulated genes
# nrow(up)
# nrow(down)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# write(rownames(up), "sig_sg_vs_cg/up.txt")
# write(rownames(down), "sig_sg_vs_cg/down.txt")
# 
# 
# ###### Compare sorted lf v.s. whole lf ########
# 
# 
# res <- results(dds, contrast=c("combined","sorted_loss","whole_loss"))
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj < 0.05, na.rm=TRUE)
# 
# resSig <- subset(resOrdered, padj < 0.05)
# 
# up <- subset(resSig, log2FoldChange > 2) # Up-regulated genes
# down <- subset(resSig, log2FoldChange < -2) # Down-regulated genes
# nrow(up)
# nrow(down)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# write(rownames(up), "sig_sl_vs_cl/up.txt")
# write(rownames(down), "sig_sl_vs_cl/down.txt")
# 
# ###### Compare sorted wt v.s. whole wt ########
# 
# res <- results(dds, contrast=c("combined","sorted_wild","whole_wild"))
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj < 0.05, na.rm=TRUE)
# 
# 
# resSig <- subset(resOrdered, padj < 0.05) # Significant genes
# 
# up <- subset(resSig, log2FoldChange > 2) # Up-regulated genes
# down <- subset(resSig, log2FoldChange < -2) # Down-regulated genes
# nrow(up)
# nrow(down)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# write(rownames(up), "sig_sw_vs_cw/up.txt")
# write(rownames(down), "sig_sw_vs_cw/down.txt")
# 
#######################################


# View results for non-significant genes
# P-adj threshold: 0.05
# Compared: sorted_genotype vs whole_genotype

###### Compare sorted gf v.s. whole gf ########
# res <- results(dds, contrast=c("combined","sorted_gain","whole_gain"))
# 
# resultsNames(dds)
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj > 0.05, na.rm=TRUE)
# 
# resnonSig <- subset(resOrdered, padj > 0.05)
# 
# up <- subset(resnonSig, log2FoldChange > 0) # Up-regulated genes
# down <- subset(resnonSig, log2FoldChange < 0) # Down-regulated genes
# nrow(up)   
# nrow(down)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# write(rownames(up), "nonsig_sg_vs_cg/up.txt")
# write(rownames(down), "nonsig_sg_vs_cg/down.txt")
# 
# ###### Compare sorted lf v.s. whole lf ########
# 
# 
# res <- results(dds, contrast=c("combined","sorted_loss","whole_loss"))
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj > 0.05, na.rm=TRUE)
# 
# resnonSig <- subset(resOrdered, padj > 0.05)
# 
# up <- subset(resnonSig, log2FoldChange > 0) # Up-regulated genes
# down <- subset(resnonSig, log2FoldChange < 0) # Down-regulated genes
# nrow(up)
# nrow(down)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# write(rownames(up), "nonsig_sl_vs_cl/up.txt")
# write(rownames(down), "nonsig_sl_vs_cl/down.txt")
# 
# ###### Compare sorted wt v.s. whole wt ########
# 
# res <- results(dds, contrast=c("combined","sorted_wild","whole_wild"))
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj > 0.05, na.rm=TRUE)
# 
# resnonSig <- subset(resOrdered, padj > 0.05) # Significant genes
# 
# up <- subset(resnonSig, log2FoldChange > 0) # Up-regulated genes
# down <- subset(resnonSig, log2FoldChange < 0) # Down-regulated genes
# nrow(up)
# nrow(down)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# write(rownames(up), "nonsig_sw_vs_cw/up.txt")
# write(rownames(down), "nonsig_sw_vs_cw/down.txt")
# 
# 
# 
###############################################




# res <- results(dds, contrast=c("conditions","sorted","whole"), 
#                alpha = 0.05, lfcThreshold = 0.01)
# 
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj < 0.05, na.rm=TRUE)
# 
# resSig <- subset(resOrdered, padj < 0.05)
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# nrow(up)
# nrow(down)




###### Check sorted_gain vs sorted_wild #####
res <- results(dds, contrast=c("combined","sorted_gain","sorted_wild"))


resOrdered <- res[order(res$pvalue),]

sum(res$padj < 0.05, na.rm=TRUE)

resSig <- subset(resOrdered, padj < 0.05)

up <- subset(resSig, log2FoldChange > 0) # Up-regulated genes
down <- subset(resSig, log2FoldChange < -0) # Down-regulated genes

remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)

up <- up[ !(row.names(up) %in% remove_up), ]
down <- down[ !(row.names(down) %in% remove_down), ]

nrow(up)
nrow(down)

write(rownames(up), "new/sorted_wt_vs_sorted_gf/up.txt")
write(rownames(down), "new/sorted_wt_vs_sorted_gf/down.txt")


# ###### Check sorted_wild vs sorted_loss #####
# res <- results(dds, contrast=c("combined","sorted_wild","sorted_loss"))
# 
# 
# resOrdered <- res[order(res$pvalue),]
# 
# sum(res$padj < 0.05, na.rm=TRUE)
# 
# resSig <- subset(resOrdered, padj < 0.05)
# 
# up <- subset(resSig, log2FoldChange > 0) # Up-regulated genes
# down <- subset(resSig, log2FoldChange < 0) # Down-regulated genes
# 
# 
# remove_up <- grep("^[^WBGene].*$", rownames(up), value=TRUE)
# remove_down <- grep("^[^WBGene].*$", rownames(down), value=TRUE)
# 
# up <- up[ !(row.names(up) %in% remove_up), ]
# down <- down[ !(row.names(down) %in% remove_down), ]
# 
# nrow(up)
# nrow(down)

# write(rownames(up), "sorted_wt_vs_sorted_lf/up.txt")
# write(rownames(down), "sorted_wt_vs_sorted_lf/down.txt")


