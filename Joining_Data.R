#Run this code last
#Loading Libraries
library(dplyr)
library(purrr)
library(readxl)
library(tibble)

#Removed paper 8 due to outliers

#Loading in the data
paper1 <- read.csv("paper1", row.names = 1)
paper2 <- read.csv("paper2", row.names = 1)
paper3 <- read.csv("paper3", row.names = 1)
colnames(paper3) <- sub("^X", "", colnames(paper3))
paper4 <- read.csv("paper4", row.names = 1)
paper5 <- read.csv("paper5", row.names = 1)
colnames(paper5) <- sub("^X", "", colnames(paper5))
paper6 <- read.csv("paper6", row.names = 1)
paper7 <- read.csv("paper7", row.names = 1)
paper9 <- read.csv("paper9", row.names = 1)
colnames(paper9) <- sub("^X", "", colnames(paper9))

#Converting rownames to column names called 'ENSEMBL'
paper1 <- paper1 %>% rownames_to_column(var = "ENSEMBL")
paper2 <- paper2 %>% rownames_to_column(var = "ENSEMBL")
paper3 <- paper3 %>% rownames_to_column(var = "ENSEMBL")
paper4 <- paper4 %>% rownames_to_column(var = "ENSEMBL")
paper5 <- paper5 %>% rownames_to_column(var = "ENSEMBL")
paper6 <- paper6 %>% rownames_to_column(var = "ENSEMBL")
paper7 <- paper7 %>% rownames_to_column(var = "ENSEMBL")
paper9 <- paper9 %>% rownames_to_column(var = "ENSEMBL")

#Combining papers into a list
papers_list <- list(paper1, paper2, paper3, paper4, paper5, paper6, paper7, paper9)

#Join papers together by merging them on the ENSEMBL column, creating a full join across all papers
joined_papers <- reduce(papers_list, ~full_join(as.data.frame(.x, keep.rownames = "ENSEMBL"), 
                                             as.data.frame(.y, keep.rownames = "ENSEMBL"), by = "ENSEMBL"))

#Get column 'ENSEMBL' back into rownames
joined_papers <- joined_papers %>% column_to_rownames(var = "ENSEMBL")

#Remove any rows containing NA as data analysis will not work otherwise
joined_papers <- na.omit(joined_papers)

#Create a box plot of all papers data
boxplot(joined_papers,
        main = 'Boxplot of Samples',
        ylab = 'Gene Expression')


#Loading metadata for PCA
metadata <- read_excel("Master_for_R.xlsx", sheet = 1)
rownames(metadata) <- metadata$GEO_ID

#prepping metadata and data for PCA

metadata_sel<-metadata[metadata$Original_ID%in%colnames(joined_papers),]
metadata_sel$homoGroup<-metadata_sel$Sample_Type

#check that metadata is in the same order of the data
metadata_sel<-metadata_sel[order(metadata_sel$Original_ID),]#not needed really
joined_papers<-joined_papers[,metadata_sel$Original_ID]
colnames(joined_papers)==metadata_sel$Original_ID

CBF_PCA(data = t(joined_papers),groups = metadata_sel$homoGroup,scale = F,legendName = "group")
