#' Model creation and editing
#'
#' Generate an scGate model from scratch or edit an existing one
#' 
#' @param model scGate model to be modified. When is NULL (default) a new model will be initialized.   
#' @param level integer. It refers to the hierarchical level of the model tree in which the signature will be added (level=1 by default)    
#' @param name Arbitrary signature name (i.e. Immune, Tcell, NK etc).   
#' @param signature character vector indicating gene symbols to be included in the signature (e.g. CD3D). If a minus sign is placed to the end of a gene name (e.g. "CD3D-"), this gene will be used as negative in UCell computing. See UCell documentation for details    
#' @param positive Logical indicating if the signature must be used as a positive signature in those model level. Default is TRUE. 
#' @param negative Same as `positive` but negated (negative=TRUE equals to positive=FALSE)
#' @param remove Whether to remove the given signature from the model
#' @return A scGate model that can be used by \code{\link{scGate}} to filter target cell types.
#' @examples
#' # create a simple gating model
#' my_model <- gating_model(level = 1, name = "immune", signature = c("PTPRC"))
#' my_model <- gating_model(model = my_model, level = 1, positive = FALSE,
#'     name = "Epithelial", signature = c("CDH1","FLT1") )
#' # Remove an existing signature
#' dropped_model <- gating_model(model = my_model, remove =TRUE, level = 1, name = "Epithelial")
#' @importFrom stats setNames
#'
#' @export
#'
gating_model <- function(model=NULL, level= 1, name, signature,
                         positive = TRUE, negative = FALSE, remove = FALSE){
  
  template <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("levels","use_as", "name", "signature"))
  
  if(negative){
    positive <- FALSE
  }
  
  if(is.null(model)){
    model <- template 
  }
  
  if(!remove){
    new.signature <- data.frame(levels = paste0("level",level),
                                use_as = ifelse(positive, "positive","negative"),
                                name = name,
                                signature = ifelse(length(signature) >1, paste(signature,collapse = ";") ,signature))
    model <- rbind(model,new.signature)
  }else{
    lev <- paste0("level",level)
    
    model <- model[!((model$levels == lev) & (model$name == name)),]
  }
  return(model)
}
