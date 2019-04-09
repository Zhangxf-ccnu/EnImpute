library("Seurat")


#' Run data visualization and cell clustering using the R pacakage Seurat
#'
#' @param x Raw count matrix, where the rows represent genes and the columns denote cells
#' @param pcs Number of principal components
#' @param label Subtype label of cells. If label==NULL, the FindClusters function in the the R pacakage Seurat
#' will be used to generate the subtype label.
#' @param res Value of the resolution paramete
#'
#' @return a Seurat object with the PCA calculation stored in object@dr$pca,
#' with a tSNE embedding in object@dr$tsne@cell.embeddings,
#' and with subtypes stored in object@ident
#'
tsne_cluster <- function(x,  pcs = 20, label = NULL, res = 1) {

  x.seurat <- CreateSeuratObject(raw.data = x)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- ScaleData(x.seurat)
  x.seurat <- FindVariableGenes(x.seurat, do.plot = FALSE)
  x.seurat <- RunPCA(x.seurat,  pc.genes = x.seurat@var.genes, pcs.compute = pcs, do.print = FALSE)


  if (is.null(x.seurat@dr$tsne)) {
    x.seurat <- RunTSNE(x.seurat, dims.use = 1:pcs, check_duplicates = FALSE,
                        do.fast = TRUE)
  }

  
  y.seurat <- x.seurat
  if(is.null(label)){
    
    y.seurat <- FindClusters(y.seurat, dims.use = NULL, resolution = res,
                             print.output = FALSE, save.SNN = TRUE)
    x.seurat@ident  <- y.seurat@ident 
  }else{
    x.seurat@ident =  as.factor(label)
  }
  
  
  x.seurat
}



normalizeData <- function(x, y = x) {
  sf <- colSums(y)/1000000
  return(sweep(x, 2, sf, "/"))
}

###############################################################################
## Correlation with reference plots
get.cor.gene <- function(X, Y) {
  sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ]))
}

get.cor.cell <- function(X, Y) {
  sapply(1:ncol(X), function(i) cor(X[, i], Y[, i]))
}

## Correlation matrix distance
calc_cmd <- function(R1, R2) {
  traceR1R2 <- sum(diag(crossprod(R1, R2)))
  R1.norm <- norm(R1, type = "F")
  R2.norm <- norm(R2, type = "F")
  return(1-traceR1R2/(R1.norm*R2.norm))
}


