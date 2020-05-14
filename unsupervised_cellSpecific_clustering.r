# function: 
# to do unsupervised clustering within a given cell type. It uses dubstep genes.
# input:
# pbmc.fun is seruat object
# cell.type is the cell type you want to focus on
# res is the resolution for unsupervised clustering
# suffix is any string you want to save your plot names with
# var.toRegress are the variables to regress while doing scTransform on cell-type specific data
# output:
# this function returns a list of annotated cell-type specific seurat object, and list of plots

unsupervised_cellSpecific_clustering <- function(pbmc.fun, cell.type, res, suffix, var.toRegress){
  
  print(paste("inferring subclusters using unsupervised clustering for", cell.type, "at resolutions", res, "..."))
  
  # filter BT cells from raw object
  pbmc.specific = subset(x = pbmc.fun, subset = CELL.TYPE == cell.type)
  print("cell-specific pbmc object prepared")
  print(dim(pbmc.specific))
  
  # normalise
  pbmc.specific = normalise_seurat(pbmc.fun = pbmc.specific, nGenes = nFeaturesGenes, var.toRegress = var.toRegress)
  print("cell-specific normalisation done")
  
  # calculate dubstep feature genes for BT cells 0.01*ncol(seuratObj), k = 10, num.pcs = 15
  dubstep.features = DUBStepR(input.data = pbmc.specific@assays$SCT@data, min.cells = 0.01*nrow(pbmc.specific), k = 5, num.pcs = 15)
  print(paste(length(dubstep.features), "dubstep genes calculated"))

  # apend dubstep variable genes
  pbmc.specific@assays$SCT@var.features = intersect(dubstep.features, rownames(pbmc.specific@assays$SCT@data))
  print("dubstep genes appended to seurat object")
  
  # run PCA and use them for UMAP
  pbmc.specific = RunPCA(object = pbmc.specific, verbose = F, assay = "SCT")
  umap.defaults$input = "data"
  pbmc.specific.umap = umap(d = Embeddings(object = pbmc.specific, reduction = "pca"), config = umap.defaults)
  colnames(pbmc.specific.umap$layout) = paste("UMAP", 1:2, sep = "")
  pbmc.specific[["umap"]] = CreateDimReducObject(embeddings = pbmc.specific.umap$layout, key = "UMAP_", assay = DefaultAssay(pbmc.specific))
  
  # neighbours and clustering
  pbmc.specific = FindNeighbors(object = pbmc.specific, reduction = "pca")
  pbmc.specific = FindClusters(object = pbmc.specific, resolution = res, algorithm = 2)
  print("neighbours found and clustering done")
  print(table(pbmc.specific$seurat_clusters))
  
  # make plots
  p.bt.1 = DimPlot(pbmc.specific, label = T, reduction = "pca", split.by = "COND", pt.size = 1.5, cols = col.arr)
  p.bt.2 = DimPlot(pbmc.specific, label = T, reduction = "pca", group.by = "LIB", pt.size = 1.5, cols = col.arr)
  print("plots prepared")
  
  # print plots
  # ggsave(filename = paste0("plots/umap/", suffix, ".", res, ".subtypes.pdf"), plot = CombinePlots(plots = list(p.bt.1, p.bt.2), ncol = 2), device = "pdf", width = 10, height = 4)
  # print(paste(cell.type,"sub-clusters plotted"))
  
  return(list(pbmc.specific, p.bt.1, p.bt.2))
  
}
