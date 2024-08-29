#Clean sample code to remove the prefix and [] from the title in metadata
cleansample <- function(metadata_modified) {
  finalvec <- rep(0, nrow(metadata_modified))
  
  for (i in 1:nrow(metadata_modified)) {
    sample <- metadata_modified$sample[i]
    
    if (grepl('\\[.*\\]', sample)) {
      #Remove prefix and []
      sample <- gsub('.*\\[|\\]', '', sample)
    }
    
    finalvec[i] <- sample
  }
  return(finalvec)
}

