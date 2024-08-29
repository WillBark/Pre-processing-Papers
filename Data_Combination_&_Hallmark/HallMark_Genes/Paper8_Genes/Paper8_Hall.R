#Load Data
paper8 <- read.csv('paper8', row.names = 1)

paper8_abs <- abs(paper8)

#Calculate Mean and Standard Deviation for each row (gene)
mean_values <- rowMeans(paper8_abs)
sd_values <- apply(paper8_abs, 1, sd)

#Calculate coefficient of variance (CV)
CV_paper8 <- sd_values / mean_values

#Add the CV values to the original data frame
paper8$CV <- CV_paper8

paper8_sorted <- paper8[order(-paper8$CV), ]

#Calculates number of genes in lowest 5% of dataset
num_genes <- ceiling(0.05 * nrow(paper8_sorted))
paper8_lowest_5 <- paper8_sorted[(nrow(paper8_sorted) - num_genes + 1):nrow(paper8_sorted), ]

#Mean expression calculated for each row
paper8_mean_expression <- rowMeans(paper8_lowest_5[, -ncol(paper8_lowest_5)])

#Adding column
paper8_lowest_5$mean_expression <- paper8_mean_expression

paper8_sorted_expression <- paper8_lowest_5[order(-paper8_lowest_5$mean_expression), ]

#Calculate number of genes corresponding to top 5%
paper8_top_genes <- ceiling(0.05 * nrow(paper8_sorted_expression))

#Select top 5% with highest mean expression
paper8_top_5 <- paper8_sorted_expression[1:paper8_top_genes, ]

paper8_gene_names <- rownames(paper8_top_5)

#Dataframe with these gene names
paper8_gene_names <- data.frame(Gene = paper8_gene_names)

#Save as CSV
write.csv(paper8_gene_names, 'paper8_gene_names.csv', row.names = FALSE)