library(Seurat)
library(anndata)
library(dplyr) 
method <- snakemake@params[["method"]]
dims <- 50
features <- 2000


path2name <- function(path, ending = ".h5ad" ){
    fname <- path %>% strsplit("/") %>% unlist() %>% tail(n=1)
    name <- fname %>% strsplit(ending) %>% unlist()
    return (name)
}
 
ad <- read_h5ad(snakemake@input[["reference"]])
name <- path2name(snakemake@input[["reference"]])
srt <- CreateSeuratObject(counts = t(ad$X), meta.data = ad$obs)
seurat_list = list()
seurat_list[[paste(name)]] <- srt

print("Loading other datasets")
for (file in snakemake@input[["others"]]){
    ad <- read_h5ad(file)
    name <- path2name(file)
    srt <- CreateSeuratObject(counts = t(ad$X), meta.data = ad$obs)
    seurat_list[[paste(name)]] <- srt
}

print(seurat_list)
print("Normalizing individual datasets")
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x %<>% NormalizeData(verbose=FALSE) %>% 
    FindVariableFeatures(verbose=FALSE, nfeatures = 2000)
})

print("PCA")
# PCA on features that are variable across data sets
features <- SelectIntegrationFeatures(seurat_list, nfeatures = features)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x %<>% ScaleData(features = features, verbose=FALSE) %>%
    RunPCA(features = features, verbose = FALSE)
})

# Do the integration based on PCA 
print("Find anchors")
print(seurat_list)
anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                  reduction = method, reference=c(1),
                                  dims = 1:50, anchor.features=2000)
print("Integrate!")
integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

# Now its business as usual. Cluster at low resolution
integrated %<>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:dims) %>%
                FindNeighbors() %>% FindClusters(resolution=0.1)

# print("Saving the reference")
# ad <- AnnData(X = GetAssayData(object = seurat_list[[1]]))
# ad <- t(ad) # obs x features
# ad$obs <- seurat_list[[1]]@meta.data
# write_h5ad(ad, snakemake@output[["reference"]])

print("Save")
for (file in snakemake@output){
    name <- path2name(file, "_integrated.h5ad")
    #paste0("_integrated_", method, ".h5ad"))
    srt <- seurat_list[[paste(name)]]
    ad <- AnnData(X = GetAssayData(object = srt))
    ad <- t(ad) # obs x features
    ad$obs <- srt@meta.data
    write_h5ad(ad, file)
}

