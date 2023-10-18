find.nn <- function(q, assay = "RNA", slot="data", signatures=NULL, npca=30,
                    nfeatures=2000, k.param=10, bpp=SerialParam(),
                    min.cells=30, reduction="calculate", genes.blacklist=NULL) {
  
  DefaultAssay(q) <- assay
  ncells <- length(Cells(q))
  ngenes <- nrow(q)
  
  notfound <- signatures[!signatures %in% colnames(q[[]])]
  signatures <- signatures[signatures %in% colnames(q[[]])]
  
  if (length(notfound)>0) {
    message(paste0("Warning: signatures not found: ", notfound))
  }
  
  if(ncells < min.cells){  #Do not do knn-smoothing
    return(q)
  }  
  
  if (reduction=="calculate") {
    if (ngenes < nfeatures) {
      nfeatures <- ngenes
    }
    if (ngenes/2 < npca) {
      npca <- ngenes/2
    }
    
    if (slot=="counts") { 
      q <- NormalizeData(q, verbose = FALSE)
    }
    if (ngenes > 200) {  #only perform this for high-dim data
      q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    } else {
      q@assays[[assay]]@var.features <- rownames(q)
    }
    
    q@assays[[assay]]@var.features <- setdiff(q@assays[[assay]]@var.features, genes.blacklist)
    
    q <- ScaleData(q, verbose=FALSE)
    q <- RunPCA(q, features = q@assays[[assay]]@var.features, npcs=npca, verbose = FALSE, reduction.key = "knnPCA_")
    
    red.use <- 'pca'
  } else {
    red.use <- reduction
  }
  
  #Smooth scores by kNN neighbors
  q <- SmoothKNN(q, signature.names=signatures, reduction=red.use,
                 k=k.param, suffix = NULL, BPPARAM=bpp)
  
  return(q)
  
}
