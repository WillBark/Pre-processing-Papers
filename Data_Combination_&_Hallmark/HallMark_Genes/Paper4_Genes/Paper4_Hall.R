#Load Data
paper4 <- read.csv('paper4', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate the absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper4, 1, abs_mean)
row_abs_sd <- apply(paper4, 1, abs_sd)

#Calculate the Coefficient of Variation (CV) for each row
CV_paper4 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper4$CV <- CV_paper4

paper4_sorted <- paper4[order(-paper4$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper4_sorted))
paper4_lowest_5 <- paper4_sorted[(nrow(paper4_sorted) - num_genes + 1):nrow(paper4_sorted), ]

#Mean expression calculated for each row
paper4_mean_expression <- rowMeans(paper4_lowest_5[, -ncol(paper4_lowest_5)])

#Adding column
paper4_lowest_5$mean_expression <- paper4_mean_expression

paper4_sorted_expression <- paper4_lowest_5[order(-paper4_lowest_5$mean_expression), ]

#Calculate number of genes corresponding to top 5%
paper4_top_genes <- ceiling(0.05 * nrow(paper4_sorted_expression))

#Select top 5% with highest mean expression levels
paper4_top_5 <- paper4_sorted_expression[1:paper4_top_genes, ]

paper4_gene_names <- rownames(paper4_top_5)

#Dataframe with these gene names
paper4_gene_names <- data.frame(Gene = paper4_gene_names)

#Save as CSV
write.csv(paper4_top_5, 'paper4_gene_names.csv', row.names = TRUE)
