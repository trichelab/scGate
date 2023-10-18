## Filter by mean
filter_bymean <- function(q, positive, negative, pos.thr=0.1, neg.thr=0.2, assay="RNA") {
  
  DefaultAssay(q) <- assay
  ncells <- ncol(q)
  
  positive <- positive[positive %in% colnames(q[[]])]
  negative <- negative[negative %in% colnames(q[[]])]
  
  cols <- c(positive, negative)
  means <- list()
  
  scores <- q[[]][,cols, drop=FALSE]

  if(length(positive)>1){
    pos <- scores[,positive] |> apply(1,max)
  }else{
    pos <- scores[,positive]
  }
  if(length(negative)>1){
    neg <- scores[,negative] |> apply(1,max)
  }else{
    neg <- scores[,negative]
  }
  ispure <- rep("Impure", ncol(q))
  
  if(length(positive)>0 & length(negative)>0) {
    ispure[pos>pos.thr & neg<neg.thr] <- "Pure"
  } else if (length(positive)>0) {
    ispure[pos>pos.thr] <- "Pure"
  } else if (length(negative)>0) {
    ispure[neg<neg.thr] <- "Pure"
  } else {
    stop("No valid signatures were provided.")
  }
      
  q@meta.data[,"is.pure"] <- ispure
  
  return(q)
}
