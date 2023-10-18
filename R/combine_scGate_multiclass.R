#' Combine scGate annotations
#'
#' If a single-cell dataset has precomputed results for multiple scGate models, combined them in multi-class annotation
#' 
#' @param obj Seurat object with scGate results for multiple models stored as metadata   
#' @param prefix Prefix in metadata column names for scGate result models
#' @param scGate_classes Vector of scGate model names. If NULL, use all columns that start with "prefix" above.
#' @param min_cells Minimum number of cells for a cell label to be considered
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param out_column The name of the metadata column where to store the multi-class cell labels
#' @return A Seurat object with multi-class annotations based on the combination of multiple models. A new
#'      column (by default "scGate_multi") is added to the metadata of the Seurat object.
#' @import Seurat
#' @examples
#' \donttest{
#' # Define gating models
#' model.B <- gating_model(name = "Bcell", signature = c("MS4A1")) 
#' model.T <- gating_model(name = "Tcell", signature = c("CD2","CD3D","CD3E"))
#' # Apply scGate with these models
#' data(query.seurat)
#' query.seurat <- scGate(query.seurat, model=model.T,
#'     reduction="pca", output.col.name = "is.pure_Tcell")
#' query.seurat <- scGate(query.seurat, model=model.B,
#'     reduction="pca", output.col.name = "is.pure_Bcell")
#' query.seurat <- combine_scGate_multiclass(query.seurat, scGate_class=c("Tcell","Bcell"))      
#' table(query.seurat$scGate_multi)
#' }
#'
#' @export
#'
combine_scGate_multiclass <- function(obj,
                                  prefix="is.pure_",
                                  scGate_classes=NULL,
                                  min_cells=20,
                                  multi.asNA = FALSE,
                                  out_column="scGate_multi"
){
  #Use all columns with given prefix
  if (is.null(scGate_classes)){  
    cols <- grep(prefix, colnames(obj@meta.data), value = TRUE)
    cols <- grep("\\.level\\d+$", cols, invert=TRUE, perl=TRUE, value=TRUE)
  } else {
    cols <- paste0(prefix, scGate_classes)
    cols <- cols[cols %in% colnames(obj@meta.data)]
  }
  if (is.null(cols)) {
    stop("Could not find scGate annotations in this object metadata.")
  }
  
  meta <- obj@meta.data[,cols, drop=FALSE]
  meta[is.na(meta)] <- "Impure"  #Avoid NAs
  obj.logical <- meta=="Pure"
  
  label.sums <- apply(obj.logical,1,sum)
  
  obj.single <- obj.logical[label.sums==1, , drop=FALSE]
  obj.single.labels <- apply(obj.single,1,function(x) names(x)[x])
  #remove prefix
  if (!is.null(prefix)) {
    obj.single.labels <- gsub(prefix, "", obj.single.labels)
  }
  
  #Assign labels to uniquely identified cells
  labs <- rep(NA, ncol(obj))
  names(labs) <- colnames(obj)
  
  labs[names(obj.single.labels)] <- obj.single.labels
 
  #Set to NA classes with too few cells
  tt <- table(labs, useNA = "always")
  labs[labs %in% names(tt)[tt<min_cells]] <- NA
  
  if (multi.asNA) {
    labs[names(label.sums[label.sums>1])] <- NA
  } else { 
    labs[names(label.sums[label.sums>1])] <- "Multi"
  }
  
  obj@meta.data[[out_column]] <- labs
  obj
}
