#Combing the Combat Data and Metadata together

#Load libraries
library(dplyr)
library(tidyr)
library(readxl)

#Loading data
combat_data <- read.csv("combat_data.csv", row.names = 1)
metadata <- read_excel("Master_for_R.xlsx")

#Converting rownames into a column called Gene_ID
combat_data <- combat_data %>%
  tibble::rownames_to_column(var = "Gene_ID")

#Data transformation from wide to long
combat_data_long <- combat_data %>%
  pivot_longer(-Gene_ID, names_to = "Original_ID", values_to = "Expression")

#Joining combat data and metadata together on Original ID
combined_data <- combat_data_long %>%
  inner_join(metadata, by = "Original_ID")

#Saving combined data
write.csv(combined_data, "combined_data.csv", row.names = FALSE)
