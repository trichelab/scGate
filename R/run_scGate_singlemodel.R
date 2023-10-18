## run scGate on a single model (internal function)
run_scGate_singlemodel <- function(data, model, pos.thr=0.2, neg.thr=0.2, assay=NULL, slot="data",
                                   reduction="calculate", nfeatures=2000, pca.dim=30,
                                   param_decay=0.25, min.cells=30, k.param=30, bpp=SerialParam(),
                                   genes.blacklist="default", verbose=FALSE,
                                   colname="is.pure", save.levels=FALSE) {
  
  if (!inherits(model, "data.frame")) {
    stop("Invalid scGate model. Please check the format of your model")
  }
  
  list.model <- table.to.model(model)
  
  q <- data  #local copy to progressively remove cells
  tot.cells <- ncol(q)
  
  ## prepare output object (one pure/impure flag by level)
  output_by_level <- rep("Impure",length(list.model)*tot.cells)
  dim(output_by_level) <- c(tot.cells,length(list.model))
  colnames(output_by_level) <- names(list.model)
  output_by_level <- data.frame(output_by_level)
  rownames(output_by_level) <- rownames(q@meta.data)
  
  for (lev in 1:length(list.model)) {
    
    if (ncol(q) < 2)  break  # if at any level we reach a number of cells below this threshold,
                                # we skip computation, considering 'Impure' by default

    if (verbose) {
      message(sprintf("Running scGate on level %i...", lev))
    }

    pos.names <- sprintf("%s_UCell", names(list.model[[lev]]$positive))
    neg.names <- sprintf("%s_UCell", names(list.model[[lev]]$negative))
    all.names <- c(pos.names, neg.names)
    
    ##Reduce parameter complexity at each iteration
    if (param_decay < 0 | param_decay > 1) {
      stop("Parameter param_decay must be a number between 0 and 1")
    }
    
    k.use <- round((1-param_decay)**(lev-1) * k.param)
    if (reduction=="calculate") {
      pca.use <- round((1-param_decay)**(lev-1) * pca.dim)
      nfeat.use <- round((1-param_decay)**(lev-1) * nfeatures)
    } else {
      pca.use <- pca.dim
    }
    q <- find.nn(q, assay=assay, slot=slot, signatures=all.names,min.cells=min.cells,
                 nfeatures=nfeat.use, reduction=reduction, npca=pca.use, k.param=k.use,
                 bpp=bpp, genes.blacklist=genes.blacklist)
    
    q <- filter_bymean(q, positive=pos.names, negative=neg.names, assay=assay,
                       pos.thr=pos.thr, neg.thr=neg.thr)
    
    Idents(q) <- "is.pure"
    n_rem <- sum(Idents(q)=="Impure")
    frac.to.rem <- n_rem/tot.cells
    mess <- sprintf("scGate: Detected %i non-pure cells at level %i", n_rem, lev)
    if (verbose) { message(mess) }
    
    ## How many cells passed the filter
    pure.cells <- colnames(q)[Idents(q)=="Pure"]
    
    if(length(pure.cells)>0){
      # save layer output
      output_by_level[pure.cells, lev] <- "Pure"
      q <- subset(q, idents="Pure")
    } else {
      break  # in case of all cells became filtered, we do not continue with the next layer
    }
  }
  
  #Add 'pure' labels to metadata
  data <- AddMetaData(data,col.name = colname, metadata = rep("Impure",tot.cells))
  
  if(length(pure.cells)>0){
    data@meta.data[pure.cells, colname] <- "Pure"
  } else {
    message(sprintf("Warning, all cells were removed at level %i. Consider reviewing signatures or model layout...", lev))
  }
  
  data@meta.data[,colname] <- factor(data@meta.data[,colname], levels=c("Pure","Impure"))
  
  # Save output by level
  if (save.levels) {
    for(name.lev in names(list.model)){
      combname <- paste0(colname,".",name.lev)
      data <- AddMetaData(data,col.name = combname, metadata = output_by_level[[name.lev]])
      data@meta.data[,combname] <- factor(data@meta.data[,combname], levels=c("Pure","Impure"))
    }
  }
  return(data)
}
