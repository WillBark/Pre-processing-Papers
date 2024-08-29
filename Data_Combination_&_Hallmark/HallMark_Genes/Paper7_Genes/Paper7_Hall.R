#Load Data
paper7 <- read.csv('paper7', row.names = 1)

#Function to calculate absolute Mean
abs_mean <- function(x) {
  abs(mean(x))
}

#Function to calculate absolute Standard Deviation
abs_sd <- function(x) {
  abs(sd(x))
}

#Calculate absolute Mean and Standard Deviation for each row
row_abs_mean <- apply(paper7, 1, abs_mean)
row_abs_sd <- apply(paper7, 1, abs_sd)

#Calculate Coefficient of Variation (CV) for each row
CV_paper7 <- (row_abs_sd / row_abs_mean) * 100

#Add the CV values to the original data frame
paper7$CV <- CV_paper7

paper7_sorted <- paper7[order(-paper7$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper7_sorted))
paper7_lowest_5 <- paper7_sorted[(nrow(paper7_sorted) - num_genes + 1):nrow(paper7_sorted), ]

#Mean expression calculated for each row
paper7_mean_expression <- rowMeans(paper7_lowest_5[, -ncol(paper7_lowest_5)])

#Adding column
paper7_lowest_5$mean_expression <- paper7_mean_expression

paper7_sorted_expression <- paper7_lowest_5[order(-paper7_lowest_5$mean_expression), ]

#Calculate number of genes corresponding to top 5%
paper7_top_genes <- ceiling(0.05 * nrow(paper7_sorted_expression))

#Select top 5% with highest mean expression
paper7_top_5 <- paper7_sorted_expression[1:paper7_top_genes, ]

paper7_gene_names <- rownames(paper7_top_5)

#Dataframe with these gene names
paper7_gene_names <- data.frame(Gene = paper7_gene_names)

#Save as CSV
write.csv(paper7_top_5, 'paper7_gene_names.csv', row.names = TRUE)
