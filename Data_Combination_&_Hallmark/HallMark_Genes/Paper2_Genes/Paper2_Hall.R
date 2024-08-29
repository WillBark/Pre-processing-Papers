#Load Data
paper2 <- read.csv('paper2', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate the absolute mean and Standard Deviation for each row
row_abs_mean <- apply(paper2, 1, abs_mean)
row_abs_sd <- apply(paper2, 1, abs_sd)

#Calculate the Coefficient of Variation (CV) for each row
CV_paper2 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper2$CV <- CV_paper2

paper2_sorted <- paper2[order(-paper2$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper2_sorted))
paper2_lowest_5 <- paper2_sorted[(nrow(paper2_sorted) - num_genes + 2):nrow(paper2_sorted), ]

#Mean expression calculated for each row
paper2_mean_expression <- rowMeans(paper2_lowest_5[, -ncol(paper2_lowest_5)])

#Adding column
paper2_lowest_5$mean_expression <- paper2_mean_expression

paper2_sorted_expression <- paper2_lowest_5[order(-paper2_lowest_5$mean_expression), ]

#Calculate the number of genes corresponding to top 5%
paper2_top_genes <- ceiling(0.05 * nrow(paper2_sorted_expression))

#Select top 5% with the highest mean expression
paper2_top_5 <- paper2_sorted_expression[1:paper2_top_genes, ]

paper2_gene_names <- rownames(paper2_top_5)

#Create a new dataframe with these gene names
paper2_gene_names <- data.frame(Gene = paper2_gene_names)

#Save as CSV
write.csv(paper2_top_5, 'paper2_gene_names.csv', row.names = TRUE)
