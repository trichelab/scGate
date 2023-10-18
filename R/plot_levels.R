#' Plot scGate filtering results by level
#'
#' Fast plotting of gating results over each model level.
#' 
#' @param obj Gated Seurat object output of scGate filtering function
#' @param pure.col Color code for pure category 
#' @param impure.col Color code for impure category
#' @return UMAP plots with 'Pure'/'Impure' labels for each level of the scGate model
#' @examples
#' scGate.model.db <- get_scGateDB()
#' model <- scGate.model.db$human$generic$Myeloid
#' # Apply scGate with this model
#' data(query.seurat)
#' query.seurat <- scGate(query.seurat, model=model,
#'     reduction="pca", save.levels=TRUE)
#' library(patchwork)     
#' pll <- plot_levels(query.seurat)
#' wrap_plots(pll)
#' @importFrom Seurat DimPlot
#' @export
#'
plot_levels <- function(obj, pure.col = "green" ,impure.col = "gray"){
  myCols <- grep("^is.pure.", colnames(obj@meta.data),value = TRUE)
  plots <- list()
  for (myCol in myCols){
    plots[[myCol]] <- DimPlot(obj, group.by = myCol, 
                              cols = list(Pure = pure.col,Impure = impure.col)) +
      theme(aspect.ratio = 1)
  }
  return(plots)
}

