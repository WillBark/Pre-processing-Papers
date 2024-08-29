#Load Data
paper5 <- read.csv('paper5', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate the absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper5, 1, abs_mean)
row_abs_sd <- apply(paper5, 1, abs_sd)

#Calculate Coefficient of Variation (CV) for each row
CV_paper5 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper5$CV <- CV_paper5

paper5_sorted <- paper5[order(-paper5$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper5_sorted))
paper5_lowest_5 <- paper5_sorted[(nrow(paper5_sorted) - num_genes + 1):nrow(paper5_sorted), ]

#Mean expression calculated for each row
paper5_mean_expression <- rowMeans(paper5_lowest_5[, -ncol(paper5_lowest_5)])

#Adding column
paper5_lowest_5$mean_expression <- paper5_mean_expression

paper5_sorted_expression <- paper5_lowest_5[order(-paper5_lowest_5$mean_expression), ]

#Calculate number of genes corresponding to top 5%
paper5_top_genes <- ceiling(0.05 * nrow(paper5_sorted_expression))

#Select top 5% with highest mean expression
paper5_top_5 <- paper5_sorted_expression[1:paper5_top_genes, ]

paper5_gene_names <- rownames(paper5_top_5)

#Dataframe with these gene names
paper5_gene_names <- data.frame(Gene = paper5_gene_names)

#Save as CSV
write.csv(paper5_top_5, 'paper5_gene_names.csv', row.names = TRUE)
