#Load Data
paper3 <- read.csv('paper3', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate the absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper3, 1, abs_mean)
row_abs_sd <- apply(paper3, 1, abs_sd)

#Calculate Coefficient of Variation (CV) for each row
CV_paper3 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper3$CV <- CV_paper3

paper3_sorted <- paper3[order(-paper3$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper3_sorted))
paper3_lowest_5 <- paper3_sorted[(nrow(paper3_sorted) - num_genes + 1):nrow(paper3_sorted), ]

#Mean expression calculated for each row
paper3_mean_expression <- rowMeans(paper3_lowest_5[, -ncol(paper3_lowest_5)])

#Adding column
paper3_lowest_5$mean_expression <- paper3_mean_expression

paper3_sorted_expression <- paper3_lowest_5[order(-paper3_lowest_5$mean_expression), ]

#Calculate the number of genes corresponding to top 5%
paper3_top_genes <- ceiling(0.05 * nrow(paper3_sorted_expression))

#Select top 5% with highest mean expression levels
paper3_top_5 <- paper3_sorted_expression[1:paper3_top_genes, ]

paper3_gene_names <- rownames(paper3_top_5)

#New dataframe with these gene names
paper3_gene_names <- data.frame(Gene = paper3_gene_names)

#Save as CSV
write.csv(paper3_top_5, 'paper3_gene_names.csv', row.names = TRUE)
