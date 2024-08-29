#Load Libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(readxl)
library(edgeR)

#Read in Data
data <- read_excel(path = 'GSE224991_normalized_dds.xls')
head(data)
dim(data)

#Get Metadata
gse <- getGEO(GEO = 'GSE224991', GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

#Metadata pipeline
metadata_modified  <-  metadata %>%
  select(1, 2, 10, 11, 12) %>% #Selecting specific columns from metadata
  #Rename for ease of use
  rename(ID = title) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(genotype = characteristics_ch1.1) %>%
  rename(race = characteristics_ch1.2) %>%
  #Remove prefixes
  mutate(tissue = gsub('tissue: ', '', tissue)) %>%
  mutate(genotype = gsub('genotype: ', '', genotype)) %>%
  mutate(race = gsub('race: ', '', race))

#Renames the '...1' column in data to Gene
data <- data %>%
  rename(Gene = `...1`) 

#Reshaping data from wide to long
data_long  <-  data %>%
  gather(key = 'samples', value = 'FPKM', -Gene)

data_long <- data_long %>%
  mutate(samples = str_remove(samples, '^s')) #Removing "s" from samples

#Joining data_long and metadata_modified through samples and ID column
data_long <- data_long %>%
  left_join(metadata_modified, by = c('samples' = 'ID'))

#Explore data

#Make data a matrix
matrix_data <- as.matrix(data[, 2:ncol(data)])
rownames(matrix_data)<-data$Gene

#Removing "s" from data_matrix samples
remove_s <- function(name) {
  gsub('s', '', name)
}

colnames(matrix_data)[1:ncol(matrix_data)] <- sapply(colnames(matrix_data)[1:ncol(matrix_data)], remove_s)

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
boxplot(filtered_matrix, main = 'Before Normalization', las = 2)
boxplot(log2_normalized_matrix,
        main = 'Boxplot of Samples',
        ylab = 'Log2 Gene Number') #Change to ggplot eventually

#prepping metadata and data for PCA -> run the CBF_PCA_Paper3.R

metadata_modified_sel<-metadata_modified[metadata_modified$ID%in%colnames(log2_normalized_matrix),]
metadata_modified_sel$homoGroup<-metadata_modified_sel$tissue  

#check that metadata is in the same order than the data
metadata_modified_sel<-metadata_modified_sel[order(metadata_modified_sel$ID),]
log2_normalized_matrix<-log2_normalized_matrix[,metadata_modified_sel$ID]
stopifnot(colnames(log2_normalized_matrix)==metadata_modified_sel$ID)

CBF_PCA(data = t(log2_normalized_matrix),groups = metadata_modified_sel$homoGroup,scale = F,legendName = "group")

#Saving log2 data as data frame
log2_normalised_data <- as.data.frame(log2_normalized_matrix)
write.csv(log2_normalised_data, 'log2_paper3', row.names = TRUE)
