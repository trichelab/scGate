#' Filter single-cell data by cell type
#'
#' Apply scGate to filter specific cell types in a query dataset.
#' This function dispatches helpers for Seurat or SingleCellExperiment objects. 
#' (Whenever possible, the parameters are kept the same for SCE and Seurat.)
#'
#' @param data object containing query data - filtering will be applied to this
#' @param model A single scGate model, or a list of scGate models. See Details for this format
#' @param pca.dim Number of dimensions for cluster analysis
#' @param assay Seurat assay to use
#' @param slot Data slot in Seurat object
#' @param pos.thr Minimum UCell score value for positive signatures
#' @param neg.thr Maximum UCell score value for negative signatures
#' @param maxRank Maximum number of genes that UCell will rank per cell
#' @param nfeatures Number of variable genes for dimensionality reduction
#' @param k.param Number of nearest neighbors for knn smoothing
#' @param reduction Dimensionality reduction to use for knn smoothing. By default, calculates a new reduction
#'     based on the given \code{assay}; otherwise you may specify a precalculated dimensionality reduction (e.g.
#'     in the case of an integrated dataset after batch-effect correction)
#' @param pca.dim Number of principal components for dimensionality reduction
#' @param param_decay Controls decrease in parameter complexity at each iteration, between 0 and 1.
#'     \code{param_decay == 0} gives no decay, increasingly higher \code{param_decay} gives increasingly stronger decay
#' @param ncores Number of processors for parallel processing
#' @param output.col.name Column name with 'pure/impure' annotation
#' @param min.cells Minimum number of cells to cluster or define cell types
#' @param additional.signatures A list of additional signatures, not included in the model, to be evaluated (e.g. a cycling signature). The scores for this
#'     list of signatures will be returned but not used for filtering.
#' @param save.levels Whether to save in metadata the filtering output for each gating model level
#' @param keep.ranks Store UCell rankings in Seurat object. This will speed up calculations if the same object is applied again with new signatures.
#' @param genes.blacklist Genes blacklisted from variable features. The default loads the list of genes in \code{scGate::genes.blacklist.default};
#'     you may deactivate blacklisting by setting \code{genes.blacklist=NULL}
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param seed Integer seed for random number generator
#' @param verbose Verbose output
#'
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells passed the scGate filter.
#'     The \code{active.ident} is also set to this variable.
#' @details Models for scGate are data frames where each line is a signature for a given filtering level.
#'     A database of models can be downloaded using the function \code{get_scGateDB}.
#'     You may directly use the models from the database, or edit one of these models to generate your own custom gating model.
#'     
#'     Multiple models can also be evaluated at once, by running scGate with a list of models. Gating for each individual model is
#'     returned as metadata, with a consensus annotation stored in \code{scGate_multi} metadata field. This allows using scGate as a
#'     multi-class classifier, where only cells that are "Pure" for a single model are assigned a label, cells that are "Pure" for
#'     more than one gating model are labeled as "Multi", all others cells are annotated as NA.
#'
#' @examples
#' ### Test using a small toy set
#' data(query.seurat)
#' # Define basic gating model for B cells
#' my_scGate_model <- gating_model(name = "Bcell", signature = c("MS4A1")) 
#' query.seurat <- scGate(query.seurat, model = my_scGate_model, reduction="pca")
#' table(query.seurat$is.pure)
#' \donttest{
#' ### Test with larger datasets
#' library(Seurat)
#' testing.datasets <- get_testing_data(version = 'hsa.latest')
#' seurat_object <- testing.datasets[["JerbyArnon"]]
#' # Download pre-defined models
#' models <- get_scGateDB()
#' seurat_object <- scGate(seurat_object, model=models$human$generic$PanBcell)
#' DimPlot(seurat_object)
#' seurat_object_filtered <- subset(seurat_object, subset=is.pure=="Pure")
#'
#' ### Run multiple models at once
#' models <- get_scGateDB()
#' model.list <- list("Bcell" = models$human$generic$Bcell,
#'                    "Tcell" = models$human$generic$Tcell)
#' seurat_object <- scGate(seurat_object, model=model.list)
#' DimPlot(seurat_object, group.by = "scGate_multi")
#' }
#'
#' @seealso \code{\link{load_scGate_model}} \code{\link{get_scGateDB}} \code{\link{plot_tree}}
#'
#' @import Seurat
#' @import ggplot2
#'
#' @importFrom dplyr %>% distinct bind_rows
#' @importFrom UCell AddModuleScore_UCell SmoothKNN
#'
#' @import BiocParallel
#'
#' @export
#'
scGate <- function(data,
                   model,
                   # everything else is defaults 
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

  if (!is.null(assay)) {
    DefaultAssay(data) <- assay
  }
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
