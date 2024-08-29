#Load Data
paper9 <- read.csv('paper9', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper9, 1, abs_mean)
row_abs_sd <- apply(paper9, 1, abs_sd)

#Calculate Coefficient of Variation (CV) for each row
CV_paper9 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper9$CV <- CV_paper9

paper9_sorted <- paper9[order(-paper9$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper9_sorted))
paper9_lowest_5 <- paper9_sorted[(nrow(paper9_sorted) - num_genes + 1):nrow(paper9_sorted), ]

#Mean expression calculated for each row
paper9_mean_expression <- rowMeans(paper9_lowest_5[, -ncol(paper9_lowest_5)])

#Adding column
paper9_lowest_5$mean_expression <- paper9_mean_expression

paper9_sorted_expression <- paper9_lowest_5[order(-paper9_lowest_5$mean_expression), ]

#Calculate number of genes corresponding to top 5%
paper9_top_genes <- ceiling(0.05 * nrow(paper9_sorted_expression))

#Select top 5% with highest mean expression
paper9_top_5 <- paper9_sorted_expression[1:paper9_top_genes, ]

paper9_gene_names <- rownames(paper9_top_5)

#Dataframe with these gene names
paper9_gene_names <- data.frame(Gene = paper9_gene_names)

#Save as CSV
write.csv(paper9_top_5, 'paper9_gene_names.csv', row.names = TRUE)
