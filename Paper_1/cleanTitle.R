#Clean sample code to remove the prefix and [] from the title in metadata
cleansample <- function(metadata_modified) { #metadata_modified is the selection of specific columns from metadata
  finalvec <- rep(0, nrow(metadata_modified)) #finalvec stores cleaned version of samples
  
  for (i in 1:nrow(metadata_modified)) {
    sample <- metadata_modified$sample[i]
    
    if (grepl('\\[.*\\]', sample)) {
      #Remove prefix and []
      sample <- gsub('.*\\[|\\]', '', sample) #removing everything inbetween the [] for formatting
    }
    
    finalvec[i] <- sample
  }
  return(finalvec)
}

