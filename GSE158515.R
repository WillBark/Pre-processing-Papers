#Clean the sample names for metadata using the cleanTitle first

#Load Libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(edgeR)

#Read in Data
data <- read.table(file = 'GSE158515_data.txt', header=TRUE, sep='\t')
head(data)
dim(data)

data <- data[, -c(2, 3)]#Removes the GeneSymbol and GeneBiotype columns not needed for data processing

#Get Metadata
gse <- getGEO(GEO = 'GSE158515', GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]])) #Extracting the metadata from the first dataset in GSE
head(metadata)

#Metadata pipeline
metadata_modified  <-  metadata %>%
  select(1, 10, 11) %>% #Selecting specific columns
  #Renaming columns with dplyr
  dplyr::rename(sample = title) %>%
  dplyr::rename(cell_type = characteristics_ch1) %>%
  dplyr::rename(group = characteristics_ch1.1) %>%
  #Removing prefix from columns
  mutate(cell_type = gsub('cell type: ', '', cell_type)) %>%
  mutate(group = gsub('group: ', '', group))

#Apply the cleansample code to metadata
metadata_modified$sample <- cleansample(metadata_modified)

#Reshaping data from wide to long
data_long <- data %>%
  gather(key = 'samples', value = 'FPKM', -GeneID)


#Joining data_long and metadata_modified through samples column
data_long <- data_long %>%
  left_join(metadata_modified, by = c('samples' = 'sample'))

#Explore data

#Make data a matrix
matrix_data <- as.matrix(data[,4:ncol(data)])#Converts gene expression into matrix
rownames(matrix_data)<-data$GeneID#Row names to GeneID

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
        main = 'After Normalization',
        ylab = 'Log2 Gene Number')


#Prepping metadata and data for PCA -> Run the CBF_PCA_Paper1.R code

metadata_modified_sel<-metadata_modified[metadata_modified$sample%in%colnames(log2_normalized_matrix),]
metadata_modified_sel$homoGroup<-metadata_modified_sel$group  
gsub(pattern = 'MyoN.SC','myometrium',metadata_modified_sel$homoGroup)  
gsub(pattern = 'MyoF.SC','myo-closeF',metadata_modified_sel$homoGroup)
gsub(pattern = 'F.SC','fibroid',metadata_modified_sel$homoGroup) 

#Check that metadata is in the same order of the data
metadata_modified_sel<-metadata_modified_sel[order(metadata_modified_sel$sample),]
log2_normalized_matrix<-log2_normalized_matrix[,metadata_modified_sel$sample]
colnames(log2_normalized_matrix)==metadata_modified_sel$sample

#PCA analysis
CBF_PCA(data = t(log2_normalized_matrix),groups = metadata_modified_sel$homoGroup,scale = F,legendName = "group")

#Log2 data frame saving
log2_normalised_data <- as.data.frame(log2_normalized_matrix)
write.csv(log2_normalised_data, 'log2_paper1', row.names = TRUE)
