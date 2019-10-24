setwd("/Users/daylinkuboyama/Desktop/Senior/fall 19/CS 423/Project")

#load file -filename might be different on your laptop (might be best to all use the same filename)
data <- read.delim(file = "R_data.txt", header = TRUE)


#blood - untreated data
BL_UNTR_data <- data.frame(data$BL_UNTR_M38_BL_UNTR_1, data$BL_UNTR_M38_BL_UNTR_2, data$BL_UNTR_BL_UNTR_3_new,data$BL_UNTR_BL_UNTR_4_new, data$BL_UNTR_BL_UNTR_5_new, data$BL_UNTR_BL_UNTR_6_new, data$BL_UNTR_BL_UNTR_7_new)
#blood - treated data
BL_16HR_data <- data.frame(data$BL_16HR_BL_16HR_5_new, data$BL_16HR_M38_BL_16HR_1, data$BL_16HR_M38_BL_16HR_2, data$BL_16HR_M38_BL_16HR_3, data$BL_16HR_M38_BL_16HR_4)
#tumor - untreated data
TIL_UNTR_data <- data.frame(data$TIL_UNTR_M38_TIL_UNTR_1, data$TIL_UNTR_M38_TIL_UNTR_2, data$TIL_UNTR_M38_TIL_UNTR_3, data$TIL_UNTR_M38_TIL_UNTR_4, data$TIL_UNTR_M38_TIL_UNTR_5, data$TIL_UNTR_TIL_UNTR_6_new, data$TIL_UNTR_TIL_UNTR_7_new)
#tumor - treated data
TIL_16HR_data <- data.frame(data$TIL_16HR_M38_TIL_16HR_1, data$TIL_16HR_M38_TIL_16HR_2, data$TIL_16HR_M38_TIL_16HR_3, data$TIL_16HR_M38_TIL_16HR_4, data$TIL_16HR_M38_TIL_16HR_5)


#compiled data set with only numeric column types (i.e. )
all_data <- data.frame(BL_UNTR_data, BL_16HR_data, TIL_UNTR_data, TIL_16HR_data)
#vector with sums of each gene's total counts across all samples
gene_count_sums <- rowSums(all_data) 

#compiled data set with gene expression sum counts attached
data_with_sums <- data.frame(data$Ensembol_id, all_data, gene_count_sums)


#function that extracts gene expression records if sum of gene counts is >= threshold value
#and returns a new dataframe with these extracted records
data_T <- function(thresVal) {
  in_range_data <- data_with_sums[data_with_sums$gene_count_sums >= thresVal,]
  return (in_range_data)
}

