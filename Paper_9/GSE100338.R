#Load Libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(readxl)
library(edgeR)

#Read in Data
data <- read.delim("Data.txt", header = TRUE, sep = "\t")
head(data)
dim(data)

data <- data %>% select(-2)

#Get Metadata
gse <- getGEO(GEO = 'GSE100338', GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

#Rename columns to match metadata
colnames(data) <- c("Gene_symbol", metadata$title)

#Metadata pipeline
metadata_modified  <-  metadata %>%
  select(1, 10, 11, 13) %>% #Selecting required columns
  #Rename for ease of use
  rename(ID = title) %>%
  rename(ethinity = characteristics_ch1) %>%
  rename(age = characteristics_ch1.1) %>%
  rename(tissue = characteristics_ch1.3) %>%
  #Removing prefixss
  mutate(ethinity = gsub("ethinity: ", "", ethinity)) %>%
  mutate(age = gsub("age: ", "", age)) %>%
  mutate(tissue = gsub("tissue subtype: ", "", tissue))

#Reshaping data from wide to long
data_long  <-  data %>%
  gather(key = 'samples', value = 'FPKM', -Gene_symbol)

#Joining data_long and metadata_modified
data_long <- data_long %>%
  left_join(metadata_modified, by = c("samples" = "ID"))

#Explore data

#Make data a matrix
matrix_data <- as.matrix(data[, 2:ncol(data)])
rownames(matrix_data)<-data$Gene_symbol

#Filtering Dirty Data
gene_means <- rowSums(matrix_data) #Mean for each gene

overall_mean <- mean(gene_means) #Overall mean across all genes

threshold <- overall_mean * 0.05 #Threshold as 5% of the overall mean

filtered_matrix <- matrix_data[gene_means > threshold, ]
dim(filtered_matrix)

#Log2 Normalize Data
DGE <- DGEList(counts = filtered_matrix) #DGEList object


DGE <- calcNormFactors(DGE, method = "TMM") #Calculate normalization factors using TMM

log2_normalized_matrix <- cpm(DGE, log = TRUE) #Normalize counts and log2 transform

# Visual comparison
par(mfrow = c(1, 2))
boxplot(filtered_matrix, main = "Before Normalization", las = 2)
boxplot(log2_normalized_matrix,
        main = 'Boxplot of Samples',
        ylab = 'Log2 Gene Number')

#prepping metadata and data for PCA -> run the code CBF_PCA_Paper9.R

metadata_modified_sel<-metadata_modified[metadata_modified$ID%in%colnames(log2_normalized_matrix),]
metadata_modified_sel$homoGroup<-metadata_modified_sel$tissue  

#check that metadata is in the same order than the data
metadata_modified_sel<-metadata_modified_sel[order(metadata_modified_sel$ID),]#not needed really
log2_normalized_matrix<-log2_normalized_matrix[,metadata_modified_sel$ID]
colnames(log2_normalized_matrix)==metadata_modified_sel$ID

CBF_PCA(data = t(log2_normalized_matrix),groups = metadata_modified_sel$homoGroup,scale = F,legendName = "group")

#Saving log2 data as data frame
log2_normalised_data <- as.data.frame(log2_normalized_matrix)
write.csv(log2_normalised_data, 'log2_paper9', row.names = TRUE)
