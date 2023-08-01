library(Seurat)
library(anndata)
library(dplyr)
DIMS <- 50
FEATURES <- 2000

# Helper functions
path2name <- function(path, ending = ".h5ad" ){
    fname <- path %>% strsplit("/") %>% unlist() %>% tail(n=1)
    name <- fname %>% strsplit(ending) %>% unlist()
    return (name)
}

anndata2seurat <- function(ad, name){
  srt <- CreateSeuratObject(counts = t(ad$X), meta.data = ad$obs)
  srt@meta.data$name <- name
  return (srt)
}

# Load reference dataset
ad <- read_h5ad(snakemake@input[["reference"]])
name <- path2name(snakemake@input[["reference"]])
seurat_list = list()
seurat_list[[paste(name)]] <- anndata2seurat(ad, name)

# Load other datasets
print("Loading other datasets")
for (file in snakemake@input[["others"]]){
    ad <- read_h5ad(file)
    name <- path2name(file)
    seurat_list[[paste(name)]] <- anndata2seurat(ad, name)
}

print(seurat_list)
# Preprocess individual datasets
print("Normalizing individual datasets")
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x %<>% NormalizeData(verbose=FALSE) %>%
    FindVariableFeatures(verbose=FALSE, nfeatures = FEATURES)
})

print("PCA")
# PCA on features that are variable across data sets
features <- SelectIntegrationFeatures(seurat_list, nfeatures = FEATURES)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x %<>% ScaleData(features = features, verbose=FALSE) %>%
    RunPCA(features = features, verbose = FALSE)
})

# Do the integration
print("Find anchors")
print(seurat_list)
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  reduction = snakemake@params[["method"]], reference=c(1),
                                  dims = 1:DIMS, anchor.features=FEATURES)
print("Integrate!")
integrated <- IntegrateData(anchorset = anchors, dims = 1:DIMS)

# Now its business as usual. Cluster at low resolution
integrated %<>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:DIMS) %>%
                FindNeighbors() %>% FindClusters(resolution=0.1)


print("Save")
file_ending = paste0("_integrated_", snakemake@params[["method"]], ".h5ad")
for (file in snakemake@output){
    name <- path2name(file, file_ending)
    print(name)
    srt <- integrated[,integrated@meta.data$name == name]
    # TODO: PCA
    ad <- AnnData(X = GetAssayData(object = srt))
    ad <- t(ad) # obs x features
    ad$obs <- srt@meta.data
    write_h5ad(ad, file)
}
