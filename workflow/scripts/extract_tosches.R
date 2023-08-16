# Extract raw data from old Seurat object to something Python can handle
library(Seurat)

load(snakemake@input[["robj"]], verbose=TRUE)

# Extract data
counts <- attr(turtle.neurons, "raw.data") # genes x cells
meta.data <- attr(turtle.neurons, "data.info") # cells x features

write.table(meta.data, snakemake@output[["metadata"]])
write.table(counts, snakemake@output[["counts"]])
