# BiocManager::install("edge")
# library(edge)
sessionInfo()
#BiocManager::install("ggfortify")
library(ggfortify)
#BiocManager::install("limma")
library(limma)
# BiocManager::install("AnnotationHub")
library(AnnotationHub)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
# BiocManager::install("edgeR")
library(edgeR)
library(GenomicRanges)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

library(ggplot2)
feature_table = read.table(file="./feature_counts.txt", stringsAsFactors=F, sep='\t')#, row.names=TRUE, col.names=TRUE)
names(feature_table) <- as.character(unlist(feature_table[1,]))
feature_table = feature_table[-1,]
# remove low expression data
feature_table = feature_table[rowMeans(feature_table) > 10, ]

# create a SummarizedExperiment data
phenotype_table = read.table(file="phenotype.txt", stringsAsFactors=F)#col.names=TRUE, row.names=TRUE)
col_data = phenotype_table
row_data = relist(GRanges(), vector("list", length=nrow(feature_table)))
se = SummarizedExperiment(assays = list(counts = feature_table), rowRanges = row_data, colData = col_data)
#se = assays(list(counts = feature_table), rowRanges = row_data, colData = col_data)
print(se)

# make a boxplot of the expression levels for each sample
dge <- DGEList(counts = assay(se, "counts"), group = phenotype_table$age.group )
dge$samples <- merge(dge$samples, as.data.frame(colData(se)), by = 0)
png("dgecount.png", width = 350, height = 350)
boxplot(dge$counts)
dev.off()

png("dgecount_log2.png", width = 350, height = 350)
log2_dge_count = log2(dge$counts + 1)
boxplot(log2_dge_count)
dev.off()

# library(ggfortify)
# perform PCA
count_pca = prcomp(log2_dge_count, center=TRUE, scale=TRUE)
dat = data.frame(X=count_pca$rotation[,1], Y=count_pca$rotation[,2], age_group=phenotype_table$age.group, RIN=phenotype_table$RIN)

# scatterplot using PC1 and PC2, colored by RIN, shaped by age.group
ggplot(dat, aes(x=X, y=Y, shape=age_group, color=RIN)) + geom_point(size=5) + xlab("PC1") + ylab("PC2")
ggsave("pca1.png")


#--------------------------------------------------------------------------------------------------------------------------
# library(limma)
# library(edge)

# make log2 transformation and remove low expression
edata = assay(se)
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# fit the model by age.group and write results to a tab-delimited file with gene name, log2 fold-change, p-value and adjusted p-value
mod = model.matrix(~ se$age.group)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
limma_toptable = topTable(ebayes_limma,number=dim(edata)[1])
limma_table_output = limma_toptable[,c(1,4,5)]
write.table(limma_table_output, file="dif_exp_genes.txt", sep='\t', row.names=TRUE, col.names=TRUE)
head(limma_table_output)

# make the volcano plot, mark those gene with p-value less than 0.05 as red
png("volcano.png", width = 350, height = 350)
with(limma_toptable, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot"))
with(subset(limma_toptable, adj.P.Val < 0.05), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
dev.off()

print(paste0("Genes differentially expressed: ", sum(limma_table_output$adj.P.Val < 0.05)))
print(paste0("Genes differentially expressed and down-regulated from fetal to adult: ", sum(limma_table_output$adj.P.Val < 0.05 & limma_table_output$logFC > 1)))
print(paste0("Genes differentially expressed and up-regulated from fetal to adult: ", sum(limma_table_output$adj.P.Val < 0.05 & limma_table_output$logFC < -1)))


#-------------------------------------------------------------------------------------------------------------
# library(AnnotationHub)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# BiocManager::install("httr")
library(httr)
HEAD("https://annotationhub.bioconductor.org/metadata/annotationhub.sqlite3")
ah <- AnnotationHub()
ah <- subset(ah, species == "Homo sapiens")
print("Done1")
HEAD("https://annotationhub.bioconductor.org/fetch/35911")
ah_fetal = query(ah, c("EpigenomeRoadMap", "H3K4me3", "E081"))
ah_adult <- query(ah, c("EpigenomeRoadMap", "H3K4me3", "E073"))
#ah_liver <- query(ah, c("EpigenomeRoadMap", "H3K4me3", "E066"))
ah_liver <- query (ah, c ("EpigenomeRoadMap", "H3K4me3", "Liver"))
# download narrowPeak datasets
fetal_gr <- ah_fetal[[2]]
adult_gr <- ah_adult[[2]]
liver_gr <- ah_liver[[2]]
print("Done2")
# BiocManager::install("mygene")
library(mygene)
dif_exp_genes = row.names(limma_toptable[limma_toptable$adj.P.Val < 0.05,])
print("Done3")
dif_exp_gene_ids = queryMany(dif_exp_genes, scopes = "symbol", fields = "entrezgene", species = "human" )
print("Done4")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
print("Done5")
txdb_genes <- genes(txdb)
print("Done6")
dif_exp_promoters <- promoters(txdb_genes[dif_exp_gene_ids$entrezgene %in% txdb_genes$gene_id])
print("Done7")
adult_perc_peak = length(subsetByOverlaps(adult_gr, dif_exp_promoters, ignore.strand=TRUE)) / length(adult_gr)
fetal_perc_peak = length(subsetByOverlaps(fetal_gr, dif_exp_promoters, ignore.strand=TRUE)) / length(fetal_gr)
liver_perc_peak = length(subsetByOverlaps(liver_gr, dif_exp_promoters, ignore.strand=TRUE)) / length(liver_gr)
print("Done8")
print(paste0("Percentage of differentially expressed gene in adult narrowpeaks: ", round(adult_perc_peak, 3)))
print(paste0("Percentage of differentially expressed gene in fetal narrowpeaks: ", round(fetal_perc_peak, 3)))
print(paste0("Percentage of differentially expressed gene in adult liver narrowpeaks: ", round(liver_perc_peak, 3)))
odds_ratio = function(prom_counts, peak_counts, print=TRUE){
overlapMat <- matrix(0,, ncol = 2, nrow = 2)
colnames(overlapMat) <- c("in.peaks", "out.peaks")
rownames(overlapMat) <- c("in.promoters", "out.promoter")

prom <- reduce(prom_counts, ignore.strand = TRUE)
peaks <- reduce(peak_counts)
both <- intersect(prom, peaks)
only.prom <- setdiff(prom, both)
only.peaks <- setdiff(peaks, both)

overlapMat[1,1] <- sum(width(both))
overlapMat[1,2] <- sum(width(only.prom))
overlapMat[2,1] <- sum(width(only.peaks))
overlapMat[2,2] <- 1.5*10^9 - sum(overlapMat)

oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
return(oddsRatio)
}
print(paste0("Odds ratio in adult brain: ", round(odds_ratio(dif_exp_promoters, adult_gr), 2)))
print(paste0("Odds ratio in fetal brain: ", round(odds_ratio(dif_exp_promoters, fetal_gr), 2)))
print(paste0("Odds ratio in adult liver: ", round(odds_ratio(dif_exp_promoters, liver_gr), 2)))