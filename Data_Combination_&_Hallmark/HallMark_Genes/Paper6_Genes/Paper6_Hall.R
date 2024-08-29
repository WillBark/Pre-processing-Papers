#Laod Data
paper6 <- read.csv('paper6', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper6, 1, abs_mean)
row_abs_sd <- apply(paper6, 1, abs_sd)

#Calculate the Coefficient of Variation (CV) for each row
CV_paper1 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper6$CV <- CV_paper6

paper6_sorted <- paper6[order(-paper6$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper6_sorted))
paper6_lowest_5 <- paper6_sorted[(nrow(paper6_sorted) - num_genes + 1):nrow(paper6_sorted), ]

#Mean expression calculated for each row
paper6_mean_expression <- rowMeans(paper6_lowest_5[, -ncol(paper6_lowest_5)])

#Adding column
paper6_lowest_5$mean_expression <- paper6_mean_expression

paper6_sorted_expression <- paper6_lowest_5[order(-paper6_lowest_5$mean_expression), ]

#Calculate number of genes corresponding to top 5%
paper6_top_genes <- ceiling(0.05 * nrow(paper6_sorted_expression))

#Select top 5% with highest mean expression
paper6_top_5 <- paper6_sorted_expression[1:paper6_top_genes, ]

paper6_gene_names <- rownames(paper6_top_5)

#New dataframe with these gene names
paper6_gene_names <- data.frame(Gene = paper6_gene_names)

#Save as CSV
write.csv(paper6_top_5, 'paper6_gene_names.csv', row.names = TRUE)
