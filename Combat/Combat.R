#Run second
#Loading Libraries
library(sva)
library(Biobase)
library(readxl)

#Read in joined paper data
data <- read.csv('joined_papers.csv', row.names = 1)

#Read in the master metadata sheet
metadata <- read_xlsx('Master_for_R.xlsx', sheet = 1)
rownames(metadata) <- metadata$Original_ID

#Applying combat to the different batches of paper ID
batch <- metadata$`Paper ID`

mod <- model.matrix(~ Sample_Type, data = metadata)

#Correcting Batch Effect
combat_data <- ComBat(dat = as.matrix(data), batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)

#Saving corrected data
write.csv(combat_data, 'combat_data.csv', row.names = TRUE)

#prepping metadata and data for PCA

#BEFORE COMBAT
metadata_modified_sel<-metadata[metadata$Original_ID%in%colnames(data),]
metadata_modified_sel$homoGroup <- metadata_modified_sel$Sample_Type

#check that metadata is in the same order than the data
metadata_modified_sel<-metadata_modified_sel[order(metadata_modified_sel$Original_ID),]#not needed really
data<-data[,metadata_modified_sel$Original_ID]
colnames(data)==metadata_modified_sel$Original_ID


CBF_PCA(data = t(data),groups = metadata_modified_sel$homoGroup,scale = F,legendName = "group")

#AFTER COMBAT
metadata_modified_sel<-metadata[metadata$Original_ID%in%colnames(combat_data),]
metadata_modified_sel$homoGroup <- metadata_modified_sel$Sample_Type

#check that metadata is in the same order than the data
metadata_modified_sel<-metadata_modified_sel[order(metadata_modified_sel$Original_ID),]#not needed really
combat_data<-combat_data[,metadata_modified_sel$Original_ID]
colnames(combat_data)==metadata_modified_sel$Original_ID


CBF_PCA(data = t(combat_data),groups = metadata_modified_sel$homoGroup,scale = F,legendName = "group")
