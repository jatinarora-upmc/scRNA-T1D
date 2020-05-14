# function to normalise (SCtransform) data in seurat object and find given number of variable genes
# normalise (sctransform) and get variable genes
# pbmc.fun is seurat object
# nGenes is the number of variable genes you want

normalise_seurat <- function(pbmc.fun, nGenes, var.toRegress){
  # normalise, and regress-out the effect of variable (if given)
  print("normalising...")
  if(is.null(var.toRegress)){
    pbmc.fun = SCTransform(object = pbmc.fun, assay = "RNA", verbose = T, do.scale = F, do.center = F, new.assay.name = "SCT")
  } else {
    pbmc.fun = SCTransform(object = pbmc.fun, assay = "RNA", verbose = T, do.scale = F, do.center = F, new.assay.name = "SCT", vars.to.regress = var.toRegress)
  }
  
  # print summary
  print(paste("finding all variable genes and will keep", nGenes))
  pbmc.fun = FindVariableFeatures(object = pbmc.fun, assay = "SCT", nfeatures = nGenes)
  print(dim(pbmc.fun@assays$SCT@data))
  return(pbmc.fun)
}
