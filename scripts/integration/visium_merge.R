

sample_ids <- c()

srat_objects <- lapply(sample_ids, function(i) {
  path <- paste0(VISIUM_RDS_PATH, i, "_visium_processed.rds")
  srat <- readRDS(path)
  return(srat)
})

# Merge
# Filter out the QC below 500