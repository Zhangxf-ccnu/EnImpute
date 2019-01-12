tsne_cluster <- function(x,  pcs = 20, label = NULL, res = 1) {
  
  x.seurat <- CreateSeuratObject(raw.data = x)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- ScaleData(x.seurat)
  x.seurat <- FindVariableGenes(x.seurat, do.plot = FALSE)
  x.seurat <- RunPCA(x.seurat,  pc.genes = x.seurat@var.genes, pcs.compute = pcs, do.print = FALSE)
  
  if(is.null(label)){
    
    x.seurat <- FindClusters(x.seurat, dims.use = 1:pcs, resolution = res,
                             print.output = FALSE, save.SNN = TRUE)
  }else{
    x.seurat@ident =  as.factor(label)
  }
  
  if (is.null(x.seurat@dr$tsne)) {
    x.seurat <- RunTSNE(x.seurat, dims.use = 1:pcs, check_duplicates = FALSE,
                        do.fast = TRUE)
  }
  
  x.seurat
}
