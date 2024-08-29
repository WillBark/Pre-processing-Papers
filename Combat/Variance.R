#Run last
#Loading Libraries
library(variancePartition)
library(sva)
library(readxl)

#Load data before combat
expression <- read.csv("joined_papers.csv", row.names = 1)

#Load data after combat
combat_data <- read.csv('combat_data.csv', row.names = 1)

#Load metadata
metadata_variance <- read_xlsx('Master_for_R.xlsx', sheet = 1)
rownames(metadata_variance) <- metadata_variance$Original_ID

all(colnames(expression) == rownames(metadata_variance))

#Create formula for model
form <- ~ (1|`Paper ID`) + (1|Sample_Type)

#Ordering the data
combat_data <- combat_data[, order(colnames(combat_data))]
metadata_variance <- metadata_variance[order(rownames(metadata_variance)), ]

#Before Combat
Before_Combat <- fitExtractVarPartModel(expression, form, metadata_variance)
plotVarPart(Before_Combat)
summary(Before_Combat)

#After Combat
After_Combat <- fitExtractVarPartModel(combat_data, form, metadata_variance)
plotVarPart(After_Combat)
summary(After_Combat)
