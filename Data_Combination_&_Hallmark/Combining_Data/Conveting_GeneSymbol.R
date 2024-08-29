#Run this code first
#Load Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

paper3 <- read.csv("log2_paper3", row.names = 1)
colnames(paper3) <- sub("^X", "", colnames(paper3))
paper6 <- read.csv("log2_paper6", row.names = 1)
paper9 <- read.csv("log2_paper9", row.names = 1)
colnames(paper9) <- sub("^X", "", colnames(paper9))


#Function to convert gene symbols of paper to ENSEMBL IDs so can be joined with other papers
convert_to_ensembl <- function(data) {
  gene_symbols <- rownames(data)
  conversion <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
  
  #If some GeneSymbols have the same ENSEMBL ID then keep the first occurrence 
  conversion <- conversion %>% 
    distinct(SYMBOL, .keep_all = TRUE) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  #Merge the data frames together through use of the SYMBOL column, so ENSEMBL column is present
  data <- data %>% 
    rownames_to_column(var = "SYMBOL") %>% 
    inner_join(conversion, by = "SYMBOL") %>% 
    select(-SYMBOL) %>% #Remove the SYMBOL column after join
    column_to_rownames(var = "ENSEMBL")
  
  return(data)
}

#Apply the convert_to_ensembl function to each paper
paper3 <- convert_to_ensembl(paper3)
paper6 <- convert_to_ensembl(paper6)
paper9 <- convert_to_ensembl(paper9)

#Save the updated papers with ENSEMBL IDs
write.csv(paper3, 'paper3', row.names = TRUE)
write.csv(paper6, 'paper6', row.names = TRUE)
write.csv(paper9, 'paper9', row.names = TRUE)
