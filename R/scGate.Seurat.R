scGate.Seurat <- function(data,
                   model,
                   pos.thr=0.2,
                   neg.thr=0.2,
                   assay=NULL,
                   slot="data",
                   ncores=1,
                   seed=123,
                   keep.ranks=FALSE,
                   reduction=c("calculate","pca","umap","harmony","Liors_elephant"),
                   min.cells=30,
                   nfeatures=2000,
                   pca.dim=30,
                   param_decay=0.25,
                   maxRank=1500,
                   output.col.name='is.pure',
                   k.param=30,
                   genes.blacklist="default",
                   multi.asNA = FALSE,
                   additional.signatures=NULL,
                   save.levels=FALSE,
                   verbose=FALSE) {
  
  set.seed(seed)
  
  if (!is.null(assay)) DefaultAssay(data) <- assay
  assay <- DefaultAssay(data)
  
  if (assay == "integrated") { #UCell should not run on integrated assay
    if ('RNA' %in% Assays(data)) {
      assay.ucell <- 'RNA'
    } else if ('SCT' %in% Assays(data)) {
      assay.ucell <- 'SCT'
    } else {
      stop("Cannot find assays with unintegrated data in this Seurat object")
    }
  } else {
    assay.ucell <- assay
  }
  
  reduction <- reduction[1]
  if (is.null(reduction) || tolower(reduction)=="calculate") {
    reduction = "calculate"
  } else {
    if (!reduction %in% Reductions(data)) {
      stop(sprintf("Could not find reduction %s in this object. Set reduction='calculate' to compute a new dimred", reduction))
    }
    pca.dim <- ncol(data@reductions[[reduction]])
  }
    
  #check gene blacklist
  if (!is.null(genes.blacklist)) {
    if (length(genes.blacklist)==1 && genes.blacklist == "default") {  #Default
       genes.blacklist = scGate::genes.blacklist.default
    }  
    if (is.list(genes.blacklist)) {
      genes.blacklist <- unlist(genes.blacklist)
    }
    genes.blacklist <- unique(genes.blacklist)
  }
  
  #With single model, make a list of one element
  if (!inherits(model, "list")) {
    model <- list("Target" = model)
  }
  if (is.null(names(model))) {
    names(model) <- paste(output.col.name, seq_along(model), sep = ".")
  }
  
  if (ncores>1) {
    bpp <- MulticoreParam(workers=ncores)
  } else {
    bpp <- SerialParam()
  }
  
  # compute signature scores using UCell
  if (verbose) {
    message(sprintf("Computing UCell scores for all signatures using %s assay...\n", assay.ucell))
  }
  data <- score.computing.for.scGate(data, model, bpp=bpp, assay=assay.ucell,
                                     slot=slot, maxRank=maxRank, 
                                     keep.ranks=keep.ranks,
                                     add.sign=additional.signatures)
  
  for (m in names(model)) {
    
    col.id <- paste0(output.col.name, "_", m)

    data <- run_scGate_singlemodel(data, model=model[[m]], k.param=k.param,
                           param_decay=param_decay, pca.dim=pca.dim,
                           nfeatures=nfeatures, min.cells=min.cells, bpp=bpp,
                           assay=assay, slot=slot, genes.blacklist=genes.blacklist,
                           pos.thr=pos.thr, neg.thr=neg.thr, verbose=verbose,
                           reduction=reduction, colname=col.id, save.levels=save.levels)
    
    Idents(data) <- col.id
    n_pure <- sum(data[[col.id]]=="Pure")
    frac.to.keep <- n_pure/ncol(data)
    mess <- sprintf("\n### Detected a total of %i pure '%s' cells (%.2f%% of total)",
                    n_pure, m, 100*frac.to.keep)
    message(mess)
  }
 
  #Combine results from multiple model into single cell type annotation 
  data <- combine_scGate_multiclass(data, prefix=paste0(output.col.name,"_"),
                            scGate_classes = names(m), multi.asNA = multi.asNA,
                            min_cells=min.cells, out_column = "scGate_multi")

  #Back-compatibility with previous versions
  if (names(model)[1] == 'Target') {
    cn <- paste0(output.col.name, "_Target")
    data@meta.data[,output.col.name] <- data@meta.data[,cn]
    data@meta.data[,cn] <- NULL
    
    if (save.levels) {
      for (l in unique(model[[1]]$levels)) {
        cn <- paste0(output.col.name, "_Target.",l)
        data@meta.data[,paste0(output.col.name,".",l)] <- data@meta.data[,cn]
        data@meta.data[,cn] <- NULL
      }
    }
  }
  return(data)
}
