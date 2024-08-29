#Load Libraries
library(dplyr)

#Read in all the paper gene names from applying housekeeping procedure
paper1 <- read.csv('paper1_gene_names.csv', row.names = 1)
paper2 <- read.csv('paper2_gene_names.csv', row.names = 1)
paper3 <- read.csv('paper3_gene_names.csv', row.names = 1)
paper4 <- read.csv('paper4_gene_names.csv', row.names = 1)
paper5 <- read.csv('paper5_gene_names.csv', row.names = 1)
paper6 <- read.csv('paper6_gene_names.csv', row.names = 1)
paper7 <- read.csv('paper7_gene_names.csv', row.names = 1)
paper8 <- read.csv('paper8_gene_names.csv', row.names = 1)
paper9 <- read.csv('paper9_gene_names.csv', row.names = 1)

#List of data frames together
papers <- list(paper1, paper2, paper3, paper4, paper5, paper6, paper7, paper8, paper9)

#Find common gene names  across all papers
common_genes <- Reduce(intersect, lapply(papers, rownames))

#Print out common genes
print(common_genes)
