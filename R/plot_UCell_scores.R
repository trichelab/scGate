#' Plot UCell scores by level
#'
#' Show distribution of UCell scores for each level of a given scGate model
#' 
#' @param obj Gated Seurat object (output of scGate)
#' @param model scGate model used to identify a target population in obj
#' @param pos.thr Threshold for positive signatures used in scGate model (set to NULL to disable)
#' @param neg.thr Threshold for negative signatures used in scGate model (set to NULL to disable) 
#' @param overlay Degree of overlay for ggridges
#' @param ncol Number of columns in output object (passed to wrap_plots)
#' @param combine Whether to combine plots into a single object, or to return a list of plots
#' @return Returns a density plot of UCell scores for the signatures in the scGate model,
#'     for each level of the model  
#' @examples
#' scGate.model.db <- get_scGateDB()
#' model <- scGate.model.db$human$generic$Tcell
#' # Apply scGate with this model
#' data(query.seurat)
#' query.seurat <- scGate(query.seurat, model=model,
#'     reduction="pca", save.levels=TRUE)
#' # View UCell score distribution
#' plot_UCell_scores(query.seurat, model)
#' @return Either a plot combined by patchwork (combine=T) or a list of plots (combine=F)
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @importFrom patchwork wrap_plots
#' 
#' @export
#' 
plot_UCell_scores <- function(obj, model, overlay=5, pos.thr=0.2,
                              neg.thr=0.2, ncol=NULL, combine=TRUE) {
  
  u_cols <- grep('_UCell', colnames(obj@meta.data), value = TRUE)
  
  levs <- unique(model$levels)
  pll <- list()
  
  palette <- c("#00fd0c","#f4340e")
  names(palette) <- c("Positive","Negative")
  
  if (sum(grepl("is.pure.level", colnames(obj@meta.data)))==0) {
     obj$is.pure.level1 <- obj$is.pure
     
     if (length(levs)>1) {
       warning("scGate levels were not stored in this object. Showing results only for top level.")
       levs <- "level1"
     }
  }
  
  for (l in seq_along(levs)) {
    
    lev.name <- levs[l]
    
    sigs <- model[model$levels == lev.name, c("use_as","name")]
    
    col <- sprintf("%s_UCell", sigs$name)
    col <- col[col %in% u_cols]
    
    meta <- obj@meta.data
    if (l>1) {
      meta <- meta[meta[sprintf("is.pure.level%i",l-1)]=="Pure",]
    }
    ncells <- nrow(meta)
    stat <- table(meta[,sprintf("is.pure.level%i",l)])
    
    to.plot <- meta[,col, drop=FALSE]
    colnames(to.plot) <- gsub("_UCell","",colnames(to.plot))
    
    to.plot <- reshape2::melt(to.plot, id=NULL)
    colnames(to.plot) <- c("Signature","Score")
    
    to.plot$Class <- "Positive"
    to.plot$Class[to.plot$Signature %in% sigs[sigs$use_as =="negative","name"]] <- "Negative"
    
    #vertical lines (thresholds)
    to.plot$Thr <- NA
    if (!is.null(pos.thr)) {
       to.plot[to.plot$Class=="Positive","Thr"] <- pos.thr
    }
    if (!is.null(neg.thr)) {
      to.plot[to.plot$Class=="Negative","Thr"] <- neg.thr
    }
    
    #Make ggridges distribution plot
    pll[[l]] <- ggplot(to.plot, aes(x =.data$Score, y =.data$Signature, fill=.data$Class)) + 
      geom_density_ridges(scale = overlay) +
      scale_fill_manual(values = palette) + theme_minimal() +
      theme(axis.title.y=element_blank()) + ggtitle(sprintf("%s - %i/%i pure ",lev.name, stat["Pure"], ncells)) 
    
    #Add threshold lines
    if (!is.null(pos.thr) |  !is.null(neg.thr)) {  
      pll[[l]] <- pll[[l]] + geom_segment(aes(x = .data$Thr, xend = .data$Thr,
                                              y = as.numeric(.data$Signature), 
                                              yend = as.numeric(.data$Signature)+0.9), linetype = "dashed")
    }
  }
  #Return combined plot or list of plots
  if (combine) {
    return(wrap_plots(pll, ncol=ncol))
  } else {
    return(pll)
  }
}
