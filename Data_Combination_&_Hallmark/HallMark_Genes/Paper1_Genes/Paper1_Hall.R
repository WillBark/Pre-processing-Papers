#Load Data
paper1 <- read.csv('paper1', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper1, 1, abs_mean)
row_abs_sd <- apply(paper1, 1, abs_sd)

#Calculate Coefficient of Variation (CV) for each row
CV_paper1 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper1$CV <- CV_paper1

paper1_sorted <- paper1[order(-paper1$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper1_sorted))
paper1_lowest_5 <- paper1_sorted[(nrow(paper1_sorted) - num_genes + 1):nrow(paper1_sorted), ] #Selects genes from sorted dataset

#Mean expression calculated for each row
paper1_mean_expression <- rowMeans(paper1_lowest_5[, -ncol(paper1_lowest_5)])

#Adding column
paper1_lowest_5$mean_expression <- paper1_mean_expression

paper1_sorted_expression <- paper1_lowest_5[order(-paper1_lowest_5$mean_expression), ] #Sorting mean expression

#Calculate the number of genes corresponding to top 5%
paper1_top_genes <- ceiling(0.05 * nrow(paper1_sorted_expression))

#Select top 5% with the highest mean expression
paper1_top_5 <- paper1_sorted_expression[1:paper1_top_genes, ]

#Save as a CSV
write.csv(paper1_top_5, 'paper1_gene_names.csv', row.names = TRUE)
