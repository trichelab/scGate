#' Plot model tree
#'
#' View scGate model as a decision tree (require ggparty package)
#'
#' @param model A scGate model to be visualized
#' @param box.size Box size
#' @param edge.text.size Edge text size
#'
#' @return A plot of the model as a decision tree. At each level, green boxes
#'     indicate the 'positive' (accepted) cell types, red boxed indicate the
#'     'negative' cell types (filtered out). The final Pure population is the
#'     bottom right subset in the tree.
#'
#' \@\import ggparty? 
#'
#' @examples
#' library(ggparty)
#' models <- get_scGateDB()
#' plot_tree(models$human$generic$Tcell)
#'
#' @export
#'
plot_tree <- function(model, box.size = 8, edge.text.size = 4) {
  
  if (!requireNamespace('ggparty', quietly = TRUE)) {  #check whether ggparty is available
    stop("Please install and load package 'ggparty'")
  }
  nlev <- length(unique(model$levels))
  if(nlev <= 1){
    stop("your model must contain at least two levels to be ploted as a tree")
  }
  #restructure data for visualization
  level.list <- list()
  for (i in 1:nlev) {
    level.list[[i]] <- list()
    sub <- subset(model, tolower(model$levels)==paste0("level",i))
    
    level.list[[i]][["positive"]] <- sub[sub$use_as=="positive","name"]
    level.list[[i]][["negative"]] <- sub[sub$use_as=="negative","name"]
  }
  
  #Initialize dataframe for tree
  df <- data.frame(matrix(ncol=nlev+1, nrow=nlev+1, data = 0))
  colnames(df) <- c(paste0("Level_", 1:nlev), "Pure")
  
  for (i in 2:nlev) {
    for (j in 1:(i-1)) {
      df[i,j] <- 1   
    }
  }
  df[nlev+1,] <- 1
  
  ##Construct tree structure
  pn <- list()
  #bottom level
  
  pn[[nlev]] <- partykit::partynode(nlev+1,
                          split = partykit::partysplit(nlev, index=1:2, breaks = 0),
                          kids = list(partykit::partynode(nlev+2),
                                      partykit::partynode(nlev+3))) 
  
  for (i in (nlev-1):1) {
    pn[[i]] <- partykit::partynode(i,
                         split = partykit::partysplit(i, index=1:2, breaks=0),
                         kids = list(partykit::partynode(i+1),
                                     pn[[i+1]]))
  }
  
  #first element in list has complete structure
  py <- partykit::party(pn[[1]], df)
  
  
  sign.annot <- vector(length=2*nlev+1)
  is.pos <- vector(length=2*nlev+1)
  sign.annot[1] <- "Root"
  is.pos[1] <- NA
  
  for (i in 1:nlev) {
    sign.annot[2*i] <- paste0(level.list[[i]]$negative, collapse = "\n")
    sign.annot[2*i+1] <- paste0(level.list[[i]]$positive, collapse = "\n")
    
    is.pos[2*i] <- "Negative"
    is.pos[2*i+1] <- "Positive"
  }
  
  gg <- ggparty::ggparty(py)
  gg$data$info <- sign.annot
  gg$data$p.value <- is.pos
  
  
  gg$data$breaks_label[grep("<=", gg$data$breaks_label)] <- "Negative"
  gg$data$breaks_label[grep(">", gg$data$breaks_label)] <- "Positive"
  
  gg <- gg + ggparty::geom_edge() +
    ggparty::geom_edge_label(size = edge.text.size) +
    ggparty::geom_node_label(ids = "inner",
                    mapping = aes(col = .data$p.value),
                    line_list = list(aes(label= .data$info)),
                    line_gpar = list(list(size = box.size)))  +
    ggparty::geom_node_label(ids = "terminal",
                    mapping = aes(col = .data$p.value),
                    nudge_y=0.01,
                    line_list = list(aes(label= .data$info)),
                    line_gpar = list(list(size = box.size))) +
    scale_color_manual(values=c("#f60a0a", "#00ae60")) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "cm")) 
  
  return(gg)
}
