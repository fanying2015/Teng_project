library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(rGREAT)
library(soGGi)
library(limma)
library(magrittr)
library(Rsamtools)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(Biostrings)
library(clusterProfiler)
library(ChIPQC)
library(GenomicAlignments)
library(tracktables)
library(Rsubread)
library(DESeq2)
library(Rtsne)
library(ggrepel)


sample_dir <- "/Users/fanyingtang/Documents/Teng_project/"

setwd(sample_dir)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## step1. Get the peaks, which are overlap between biological replicates## step1. Get the peaks, which are the IDR peaks between biological replicates (I copied the IDR peaks and re-named them for replicates of the same sample/genotype)
ATAC_peaks <- dir("/Users/fanyingtang/Documents/Teng_project/ATAC_macs_mm10/", pattern = "*.narrowPeak", 
                      full.names = TRUE)


ATAC_peaks_all <- as.list(ATAC_peaks)
ATAC_peaks_all
names(ATAC_peaks_all) <- c("KRPS1", "KRPS1R", "KRPS3", "KRPS3R")

sample_all <- ATAC_peaks_all
peaks <- ATAC_peaks

# Do peak annotation, and get a distribution
peakAnnoList <- lapply(sample_all, annotatePeak, 
                       tssRegion=c(-3000, 3000), verbose=FALSE, annoDb="org.Mm.eg.db", TxDb=txdb)
saveRDS(peakAnnoList, paste0(sample_dir, "ATAC_cellline_peakannolist.RDS"))
pdf(paste0(sample_dir, "reads_distribution_overlap_ATACseq.pdf"))
plotAnnoBar(peakAnnoList)

dev.off()

# get the peak coordinates and readcounts
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
peak_names <- c(names(ATAC_peaks_all))
names(myPeaks) <- peak_names
peak_group <- c("KRPS", "KRPSR", "KRPS", "KRPSR")
Group <- factor(peak_group)

consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")
consensusToCount_df <- as.data.frame(consensusToCount)
write.table(consensusToCount_df, paste0(sample_dir, "consensus.bed"), quote = F, col.names = F, row.names = F, sep = "\t")


## step2. Do readcounts with bam files
library(Rsubread)
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
  rowSums
table(occurrences) %>% rev %>% cumsum
consensusToCount <- consensusToCount[occurrences >= 2, ]

ATAC_bam <- dir("/Users/fanyingtang/Documents/Teng_project/ATAC_bam_mm10/", full.names = TRUE, pattern = "*.\\.bam$")
bamsToCount <- c(ATAC_bam)

# indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
                                            start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))

#[Generate the Count table] The following step is time-consuming!!!Do count of reads
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE, 
                           countMultiMappingReads = FALSE, maxFragLength=1000, nthreads=4)

myCounts <- fcResults$counts
dim(myCounts) #130867      6
length(peak_names) #6
colnames(myCounts) <- peak_names

saveRDS(myCounts, file = paste0(sample_dir, "countsFromATAC_cellline_all_region.RDS"))
myCounts <- readRDS(paste0(sample_dir, "countsFromATAC_cellline_all_region.RDS"))



## step3. Use DESeq2 to do normalization
library(DESeq2)

# [generate the metadata] to describe/group the samples
metaData <- data.frame(Group, row.names = colnames(myCounts))
metaData$pair <- c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3")
saveRDS(metaData, "metaData_6samples.RDS")
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
saveRDS(atacDDS, "atacDDS_6samples.RDS")
atacDDS <- readRDS("atacDDS_6samples.RDS")
atacDDS <- DESeq(atacDDS) 

# get the size factors for normalization to get bigwig files later
size_factors <- colData(atacDDS)
size_factors$scale <- 1/size_factors$sizeFactor
write.csv(size_factors, file=paste0("ATAC_cellline_size_factors.csv"))

saveRDS(atacDDS, file=paste0("atacDDS_ATAC_cellline_allregion_6samples.RDS"))
atacDDS <- readRDS(paste0("atacDDS_ATAC_cellline_allregion_6samples.RDS"))
atac_Rlog <- rlog(atacDDS)
# the above step is time consuming

saveRDS(atac_Rlog, file=paste0("atacDDS_ATAC_cellline_allregion_rlog.RDS"))
atac_Rlog <- readRDS(paste0("atacDDS_ATAC_cellline_allregion_rlog.RDS"))



## step4. Do overall QC and clustering
peak_rlognormalized <- assay(atac_Rlog)
dim(peak_rlognormalized) #99218  *4

pdf(paste0("sample_clustering_DESeq2.pdf"))
plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))
plotPCA(atac_Rlog, intgroup = "pair", ntop = nrow(atac_Rlog))

library("pheatmap")
sample_cor <- cor(peak_rlognormalized, method = "pearson")
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE)
sample_cor <- cor(peak_rlognormalized, method = "pearson")
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE)
dev.off()


# do plotting based on the top 10% variable peaks
top_var_gene <- function(input_matrix, threshold) {
  select <- head(order(rowVars(input_matrix),decreasing=TRUE),threshold)
  topVarGenes <- input_matrix[select, ]
  return(topVarGenes)
}


peak_topvar <- top_var_gene(peak_rlognormalized, as.integer(dim(peak_rlognormalized)[1]/100))

simple_PCA_plot <- function(input_matrix, output_name) {
  PCA_results <- prcomp(t(input_matrix), scale=F)
  percentVar <- round(100*PCA_results$sdev^2/sum(PCA_results$sdev^2),1)
  pdf(paste0(sample_dir, output_name, ".pdf"))
  plot(PCA_results$x[,1], PCA_results$x[,2])
  dev.off()
  return(list(PCA_results, percentVar))
}


peaktopvar_PCA <- simple_PCA_plot(peak_topvar, "peak_top_1_var")

df_out <- data.frame(peaktopvar_PCA[[1]]$x,
                     tumor_type = metaData$Group,
                     sample_type = metaData$pair,
                     name = rownames(metaData))

library(ggrepel)
pca_plot_func <- function(output_name, input_mat=df_out, tumor_type=df_out$tumor_type, sample_type=df_out$sample_type, name_var=df_out$name, pca_percent=peaktopvar_PCA[[2]]) {
  pdf(paste0(output_name, ".pdf"))
  g <- ggplot(input_mat, aes(x=PC1, y=PC2, col=tumor_type, shape=sample_type, label=as.character(name_var))) + 
    geom_point(size=3) +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size = 1)) +
    ggtitle(output_name) +
   # geom_text(aes(label=as.character(name_var)),hjust=0.5,vjust=-2.5, size=2.5) +
    geom_text_repel(
      size=3.5
    ) +
    theme_bw(base_size = 12)+
    xlab(paste0(pca_percent[1], "%Var")) +
    ylab(paste0(pca_percent[2], "%Var"))
  print(g)
  
  dev.off()
}


pca_plot_func("ATAC-seq_peaks_pca_top1", df_out, df_out$tumor_type, df_out$sample_type, df_out$name, peaktopvar_PCA[[2]])

library(pheatmap)
pdf("cor_plot_of_top1varpeaks.pdf")
sample_cor <- cor(peak_topvar, method = "pearson")
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE)
sample_cor <- cor(peak_topvar, method = "spearman")
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE)
dev.off()


##step5. Use DESeq2 for differential peaks #############

#atacDDS <- readRDS(paste0(sample_dir, "results/six_samples/ATACseq_peakcomparison_v2/atacDDS_ATAC_cellline_allregion_6samples.RDS"))

## get results for each independly 
library(DESeq2)
KRPSR_vs_KRPS <- results(atacDDS, 
                              contrast=c("Group", "KRPSR", "KRPS"),
                              independentFiltering = T,
                              alpha = 0.05)
pdf(paste0(sample_dir, "MA_plots.pdf"))
DESeq2::plotMA(KRPSR_vs_KRPS, ylim=c(-4, 5),alpha=0.05, main="MA plot for KRPSR_vs_KRPS")
#DESeq2::plotMA(KRPSR_vs_KRPS)
dev.off()

summary(KRPSR_vs_KRPS)

## Get different peaks, and get an overall heatmap
df_filter <- function(sample1, sample2, atac_res, padj_fil, logfold_fil) {
  #atac_res <- atacDDS
  #sample1 <- "KRPSR"
  #sample2 <- "KRPS"
  #padj_fil <- 0.05
  #logfold_fil <- 1
  res <- results(atac_res, 
                 contrast=c("Group", sample1 , sample2),
                 independentFiltering = T,
                 alpha = padj_fil,
                 format = "GRanges")
  #print(dim(res))
  res <- as.data.frame(res)
  res <- na.omit(res)
  #print(dim(res))
  res_sec <- res[res$padj < padj_fil & (res$log2FoldChange > logfold_fil|res$log2FoldChange < -logfold_fil), ]
  return(res_sec)
}

KRPSR_vs_KRPS_df <- df_filter("KRPSR", "KRPS", atacDDS, 0.05, 1)

KRPSR_vs_KRPS_df_up <- KRPSR_vs_KRPS_df[which(KRPSR_vs_KRPS_df$log2FoldChange>0), ]
KRPSR_vs_KRPS_df_up$peak_ID <- rownames(KRPSR_vs_KRPS_df_up)
KRPSR_vs_KRPS_df_down <- KRPSR_vs_KRPS_df[which(KRPSR_vs_KRPS_df$log2FoldChange<0), ]
KRPSR_vs_KRPS_df_down$peak_ID <- rownames(KRPSR_vs_KRPS_df_down)


write.table(KRPSR_vs_KRPS_df_up[, c(12,1,2,3,4,5)], paste0(sample_dir, "KRPSR_vs_KRPS_df_up.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(KRPSR_vs_KRPS_df_down[, c(12,1,2,3,4,5)], paste0(sample_dir, "KRPSR_vs_KRPS_df_down.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
