sample_dir <- "/Users/fanyingtang/Documents/Teng_project/"

# step1. get the bed file for consensus sequence
# consensus.bed are generated from ATACseq_peakcomparison_v2.R
# cat consensus.bed |cut -f 2,3,4,5,6 |tail -n+2 | bedtools sort -i > consensus_reads.sorted.bed
sample_name <- "consensus_reads.sorted.bed"
consesus_seq <- read.table(sample_name, stringsAsFactors = F, header = F)
vec <- paste0(consesus_seq[,1], "_", consesus_seq[,2],"_", consesus_seq[,3])
consesus_seq <- cbind(consesus_seq[,c(1:3)],vec,consesus_seq[,c(4,5)]) #reorder the consensus seq
consesus_seq$vec <- as.character(vec)
head(consesus_seq)
str(consesus_seq)
write.table(consesus_seq, paste0(sample_dir, "consensus_reads.final.bed"), quote = F, row.names = F, col.names = F, sep = "\t")

# step2. generate the Homer inputs
get_sequences <- function(sample_name_1, sample_name_2, sample_name, sample_dir) {
  #sample_name_1 <- "KRPSR_vs_KRPS_df_up"
  #sample_name_2 <- "KRPSR_vs_KRPS_df_down"
  #sample_name <- "KRPSR_vs_KRPS"
  #sample_dir <- sample_dir
  test_1 <- read.table(paste0(sample_name_1, ".bed"), stringsAsFactors = F, header = F)
  test_1 <- test_1[,-1]
  vec_1 <- paste0(test_1[,1], "_", test_1[,2],"_", test_1[,3])
  test_1 <- cbind(test_1[ ,c(1:3)], vec_1, test_1[ ,c(4,5)])
  test_1$vec_1 <- as.character(vec_1)
  write.table(test_1, paste0(sample_dir, paste0(sample_name_1, ".final.bed")), quote = F, row.names = F, col.names = F, sep = "\t")
  
  test_2 <- read.table(paste0(sample_name_2, ".bed"), stringsAsFactors = F, header = F)
  test_2 <- test_2[,-1]
  vec_2 <- paste0(test_2[,1], "_", test_2[,2],"_", test_2[,3])
  test_2 <- cbind(test_2[,c(1:3)],vec_2,test_2[,c(4,5)])
  test_2$vec_2 <- as.character(vec_2)
  head(test_2)
  write.table(test_2, paste0(sample_dir, paste0(sample_name_2, ".final.bed")), quote = F, row.names = F, col.names = F, sep = "\t")
  
  colnames(consesus_seq) <- c("V1", "V2", "V3", "V4", "V5", "V6")
  colnames(test_1) <- c("V1", "V2", "V3", "V4", "V5", "V6")
  colnames(test_2) <- c("V1", "V2", "V3", "V4", "V5", "V6")
  anti_test <- anti_join(consesus_seq, test_1)
  anti_test <- anti_join(anti_test, test_2)
  write.table(anti_test, paste0(sample_dir, paste0(sample_name, ".bg.bed")), quote = F, row.names = F, col.names = F, sep = "\t")
}
